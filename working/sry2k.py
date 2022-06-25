

def SYR2K(A, B, C, alpha, beta, uplo, trans):

    assert (A!=None)
    assert (B!=None)
    assert (C!=None)
    assert uplo=='U' or uplo=='L'

    n = len(C)
    k = len(A) if trans else len(A[0])

    assert n > 0
    assert k > 0

    #Quick return
    if n==0 or ((alpha==0 or k==0) and beta==1):
        return C
    
    #Alpha==0
    if alpha==0:
        if uplo=='U':
            if beta==0:
                for j in range(n):
                    for i in range(j):
                        C[i][j] = 0
            else:
                for j in range(n):
                    for i in range(j):
                        C[i][j] = beta*C[i][j]
        else:
            if beta==0:
                for j in range(n):
                    for i in range(j, n):
                        C[i][j] = 0
            else:
                for j in range(n):
                    for i in range(j, n):
                        C[i][j] = beta*C[i][j]
        return C
    
    #Do the operation
    if trans:
        # C = alpha*A*transpose(B) + alpha*B*transpose(A) + C
        if uplo=='U':
            for j in range(n):
                if beta==0:
                    for i in range(j):
                        C[i][j] = 0
                elif beta!=1:
                    for i in range(j):
                        C[i][j] = beta*C[i][j]
                for l in range(k):
                    #Different from sryk
                    if A[j][l]!=0 or B[j][l]!=0:
                        temp1 = alpha*B[j][l]
                        temp2 = alpha*A[j][l]
                        for i in range(j):
                            C[i][j] = C[i][j] + alpha*A[i][l]*temp1 + B[i][l]*temp2
        else:
            for j in range(n):
                if beta==0:
                    for i in range(j, n):
                        C[i][j] = 0
                elif beta!=1:
                    for i in range(j, n):
                        C[i][j] = beta*C[i][j]
                for l in range(k):
                    #Different from sryk
                    if A[j][l]!=0 or B[j][l]!=0:
                        temp1 = alpha*B[j][l]
                        temp2 = alpha*A[j][l]
                        for i in range(j, n):
                            C[i][j] = C[i][j] + alpha*A[i][l]*temp1 + B[i][l]*temp2

    else:
        # C = alpha*transpose(A)*B + alpha*transpose(B)*A + C
        if uplo:
            for j in range(n):
                for i in range(j):
                    temp1=0 #Different from sryk
                    temp2=0
                    for l in range(k):
                        temp1 = temp1 + A[l][i]*B[l][j] 
                        temp2 = temp2 + B[l][i]*A[l][j]
                    if beta==0:
                        C[i][j] = alpha*temp1 + alpha*temp2
                    else:
                        C[i][j] = alpha*temp1 + beta*C[i][j] + alpha*temp2
        else:
            for j in range(n):
                for i in range(j, n):
                    temp1=0 #Different from sryk
                    temp2=0
                    for l in range(k):
                        temp1 = temp1 + A[l][i]*B[l][j] 
                        temp2 = temp2 + B[l][i]*A[l][j]
                    if beta==0:
                        C[i][j] = alpha*temp1 + alpha*temp2
                    else:
                        C[i][j] = alpha*temp1 + beta*C[i][j] + alpha*temp2

    return C