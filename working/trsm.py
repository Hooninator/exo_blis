#It seems like generally, the upper and lower triangular cases are identical save for a few for loop changes

def TRSM(A, B, alpha, uplo, side, transpose, diag):

    assert A
    assert B
    assert side=='L' or side=='R'
    assert uplo=='U' or uplo=='L'
    assert diag=='U' or 'N'
    
    m = len(B)
    n = len(B[0])
    
    assert m > 0
    assert n > 0

    if m==0 or n==0:
        return B
    
    #alpha==0
    if alpha==0:
        for j in range(n):
            for i in range(m):
                B[i][j] = 0
    
    #Operation
    if side=="L":
        if ~transpose:
            if uplo=='U':
                for j in range(n):
                    if alpha!=1: #This reappears
                        for i in range(m):
                            B[i][j] = alpha*B[i][j]
                    for k in range(m, 1, -1):
                        if B[k][j]!=0: #This if statement also reappears elsewhere in the code
                            if diag=='N':
                                B[k][j] = B[k][j]/A[k][k]
                            for i in range(k-1):
                                B[i][j] = B[i][j] - B[k][j]*A[i][k]
            else:
                for j in range(n):
                    if alpha!=1:
                        for i in range(m):
                            B[i][j] = alpha*B[i][j]
                    for k in range(m):
                        if B[k][j]!=0:
                            if diag=='N':
                                B[k][j] = B[k][j]/A[k][k]
                            for i in range(k+1, m): #NOTE: This whole chunk is the same as the non-upper case, just with different loops 
                                B[i][j] = B[i][j] - B[k][j]*A[i][k]
        else:
            if uplo=='U':
                for j in range(n):
                    for i in range(m):
                        temp = alpha*B[i][j]
                        for k in range(i-1):
                            temp -= (A[k][i]*B[k][j])
                        if diag=='N':
                            temp = temp/A[i][i]
                        B[i][j] = temp
            else:
                 for j in range(n):
                    for i in range(m, 1, -1):
                        temp = alpha*B[i][j]
                        for k in range(i+1, m):
                            temp -= (A[k][i]*B[k][j])
                        if diag=='N':
                            temp = temp/A[i][i]
                        B[i][j] = temp
    else:
        if ~transpose:
            if uplo=='U':
                for j in range(n):
                    if alpha!=1: #This reappears
                        for i in range(m):
                            B[i][j] = alpha*B[i][j]
                    for k in range(j-1):
                        if A[k][j]!=0:
                            for i in range(m):
                                B[i][j] = B[i][j] - (A[k][j]*B[i][k])
                    if diag=='N':
                        temp = 1/A[j][j]
                        for i in range(m):
                            B[i][j] = temp*B[i][j]
            else:
                for j in range(n, 1, -1):
                    if alpha!=1: #This reappears
                        for i in range(m):
                            B[i][j] = alpha*B[i][j]
                    for k in range(j+1, n):
                        if A[k][j]!=0:
                            for i in range(m):
                                B[i][j] = B[i][j] - (A[k][j]*B[i][k])
                    if diag=='N':
                        temp = 1/A[j][j]
                        for i in range(m):
                            B[i][j] = temp*B[i][j]
        else:
            if uplo=='U':
                for k in range(n, 1, -1):
                    if diag=='N':
                        temp = 1/A[k][k] #This block reappears too, but in a different spot
                        for i in range(m):
                            B[i][k] = temp*B[i][k]
                    for j in range(k-1):
                        if A[j][k]!=0:
                            temp = A[j][k]
                            for i in range(m):
                                B[i][j] = B[i][j] - (temp*B[i][k])
                    if alpha!=1:
                        for i in range(m):
                            B[i][k] = alpha*B[i][k]
            else:
                for k in range(n):
                    if diag=='N':
                        temp = 1/A[k][k] #This block reappears too, but in a different spot
                        for i in range(m):
                            B[i][k] = temp*B[i][k]
                    for j in range(k+1, n):
                        if A[j][k]!=0:
                            temp = A[j][k]
                            for i in range(m):
                                B[i][j] = B[i][j] - (temp*B[i][k])
                    if alpha!=1:
                        for i in range(m):
                            B[i][k] = alpha*B[i][k]
    return B