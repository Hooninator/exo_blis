#TODO -- how to support multiple datatypes easier like how BLIS does it? I was confused by the paper

def HEMM(A, B, C, alpha, beta, side, uplo):
    assert (A!=None)
    assert (B!=None)
    assert (C!=None)

    m = len(A)
    n = len(B[0])

    #Check parameters 
    assert m > 0
    assert n > 0
    assert side=='L' or side=='R'
    assert uplo=='U' or uplo=='L'

    #Quick return if certain parameters are 0
    if m==0 or n==0 or (alpha==0 and beta==1):
        return C
    
    #Alpha==0
    if alpha==0:
        if beta==0:
            for j in range(n):
                for i in range(m):
                    C[i][j] = 0
        else:
            for j in range(n):
                for i in range(m):
                    C[i][j]*=beta
        return C
    
    #Do the operation
    if side=='L':
        #C = alpha*A*B + beta*C
        if uplo=='U':
            for j in range(n):
                for i in range(m):
                    temp1 = alpha*B[i][j]
                    temp2 = 0
                    for k in range(i-1):
                        C[k][j]+=temp1*A[k][i]
                        temp2 = temp2+B[k][j]*conj(A[k][i])
                    if beta==0:
                        C[i][j] = temp1*A[i][i].real + alpha*temp2
                    else:
                        C[i][j] = beta*C[i][j] + temp1*A[i][i].real + alpha*temp2
        else:
            for j in range(n):
                for i in range(m, 1, -1):
                    temp1 = alpha*B[i][j]
                    temp2 = 0
                    for k in range(i+1, m):
                        C[k][j] = C[k][j] + temp1*A[k][i]
                        temp2 = temp2 + B[k][j]*conj(A[k][i])
                    if beta==0:
                        C[i][j] = temp1*A[i][i].real + alpha*temp2
                    else:
                        C[i][j] = beta*C[i][j] + temp1*A[i][i].real + alpha*temp2
    else:
        #C = alpha*B*A + beta*C
        for j in range(n):
            temp1 = alpha*A[j][j].real
            if beta==0:
                for i in range(m):
                    C[i][j] = temp1*B[i][j]
            else:
                for i in range(m):
                    C[i][j] = beta*C[i][j] + temp1*B[i][j]
            for k in range(j-1):
                if uplo=='U':
                    temp1 = alpha*A[k][j]
                else:
                    temp1 = alpha*conj(A[j][k])
                for i in range(m):
                    C[i][j] = C[i][j] + temp1*B[i][k]
            for k in range(j+1, n):
                if uplo=='U':
                    temp1 = alpha*conj(A[j][k])
                else:
                    temp1 = alpha*A[k][j]
                for i in range(m):
                    C[i][j] = C[i][j] + temp1*B[i][k]

    return C

def conj(c):
    assert type(c) is complex
    return complex(c.real, (-1*c.imag))