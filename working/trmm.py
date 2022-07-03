

def TRMM(A, B, alpha, side, uplo, transpose, diag):

    assert (A!=None)
    assert (B!=None)
    assert side=='L' or side=='R'
    assert uplo=='U' or uplo=='L'
    assert diag=='U' or diag=='N'

    m = len(B)
    n = len(B[0])
    assert m > 0
    assert n > 0

    if transpose:
        A = transpose(A)

    #Quick return
    if m==0 or n==0:
        return B

    #Alpha==0
    if alpha==0:
        for j in range(n):
            for i in range(m):
                B[i][j] = 0
        return B
    
    #Do the operation
    if side=='L':
        if ~transpose:
            if uplo=='U':
                for j in range(n):
                    for k in range(m):
                        if B[k][j]!=0:
                            temp = alpha*B[k][j]
                            for i in range(k-1):
                                B[i][j] = B[i][j] + temp*A[i][k]
                            if diag!='U':
                                temp = temp*A[k][k]
                            B[k][j] = temp
            else:
                for j in range(n):
                    for k in range(m, 1, -1):
                        if B[k][j]!=0:
                            temp = alpha*B[k][j]
                            B[k][j] = temp
                            if diag!='U':
                                B[k][j] = B[k][j]*A[k][k]
                            for i in range(k+1, m):
                                B[i][j] = B[i][j] + temp*A[i][k]
        else: #transpose A
            if uplo=='U':
                for j in range(n):
                    for i in range(m, 1, -1):
                        temp = B[i][j]
                        if diag!='U':
                            temp = temp*A[i][i]
                            for k in range(i-1):
                                temp = temp + A[k][i]*B[k][j]
                            B[i][j] = alpha*temp
            else:
                for j in range(n):
                    for i in range(m):
                        temp = B[i][j]
                        if diag!='U':
                            temp = temp*A[i][i]
                        for k in range(i+1, m):
                            temp = temp + A[k][i]*B[k][j]
                        B[i][j] = alpha*temp
    else:
        if ~transpose:
            if uplo=='U':
                for j in range(n, 1, -1):
                    temp = alpha
                    if diag!='U':
                        temp = temp*A[j][j]
                    for i in range(m):
                        B[i][j] = temp*B[i][j]
                    for k in range(j-1):
                        if A[k][j]!=0:
                            temp = alpha*A[k][j]
                            for i in range(m):
                                B[i][j] = B[i][j] + temp*B[i][k]
            else:
                for j in range(n):
                    temp = alpha
                    if diag!='U':
                        temp = temp*A[j][j]
                    for i in range(m):
                        B[i][j] = temp*B[i][j]
                    for k in range(j+1, n):
                        if A[k][j]!=0:
                            temp = alpha*A[k][j]
                            for i in range(m):
                                B[i][j] = B[i][j]+temp*B[i][k] 
        else:
            if uplo=='U':
                for k in range(n):
                    for j in range(k-1):
                        if A[j][k]!=0:
                            temp = alpha*A[j][k]
                            for i in range(m):
                                B[i][j] = B[i][j] + temp*B[i][k]
                    temp = alpha
                    if diag!='U':
                        temp = temp*A[k][k]
                    if temp!=1:
                        for i in range(m):
                            B[i][k] = temp*B[i][k]
            else:
                #Note, these transposed cases are nearly identical to non transposed, only difference is in the outer loops
                #Seems to be a few types,GEMM/SYMM/HEMM, ryk, and triangular stuff -- lots of general patterns that can be reused with slight loop changes
                for k in range(n, 1, -1):
                    for j in range(k+1, n):
                        if A[j][k]!=0:
                            temp = alpha*A[j][k]
                            for i in range(m):
                                B[i][j] = B[i][j] + temp*B[i][k]
                    temp = alpha
                    if diag!='U':
                        temp = temp*A[k][k]
                    if temp!=1:
                        for i in range(m):
                            B[i][k] = temp*B[i][k]
    return B
