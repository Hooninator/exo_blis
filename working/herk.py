#How to handle complex values in Exo/C? Is it a C struct? 


def HERK(A, C, alpha, beta, uplo, transpose):

    assert A
    assert C
    assert uplo=='U' or uplo=='L'

    n = len(C)
    k = len(A) if transpose else len(A[0])

    assert n > 0
    assert k > 0

    if n==0 or ((alpha==0 or k==0) and beta==1):
        return C
    
    if alpha==0:
        if uplo=='U':
            if beta==0:
                for j in range(n):
                    for i in range(j):
                        C[i][j] = 0
            else:
                for j in range(n):
                    for i in range(j-1):
                        C[i][j] *= beta
                    C[j][j] = beta*C[j][j].real
        else:
            if beta==0:
                for j in range(n):
                    for i in range(j, n):
                        C[i][j] = 0
            else:
                for j in range(n):
                    C[j][j] = beta*C[j][j].real
                    for i in range(j+1, n):
                        C[i][j] *= beta
        return C
    
    if ~transpose:
        #NOTE: for most of these, the upper and lower cases are almost identical, just slight changes in the loops 
        if uplo=='U':
            for j in range(n):
                if beta==0:
                    for i in range(j):
                        C[i][j] = 0 
                elif beta==1:
                    for i in range(j-1):
                        C[i][j] = beta*C[i][j]
                    C[j][j] = beta*C[j][j].real
                else:
                    C[j][j] = C[j][j].real
                for l in range(k):
                    if A[j][l]!=complex(0, 0):
                        temp = alpha*conj(A[j][l])
                        for i in range(j-1):
                            C[i][j] += (temp*A[i][l])
                        C[j][j] = C[j][j].real + (temp*A[i][l]).real
        else:
            for j in range(n):
                if beta==0:
                    for i in range(j, n):
                        C[i][j] = 0
                elif beta==1:
                    C[j][j] = beta*C[j][j].real
                    for i in range(j+1, n):
                        C[i][j]*=beta 
                else:
                    C[j][j] = C[j][j].real
                for l in range(k):
                    if A[j][l]!=complex(0, 0):
                        temp = alpha*conj(A[j][l])
                        C[j][j] = C[j][j].real + (temp*A[j][l]).real
                        for i in range(j+1, n):
                            C[i][j] += (temp*A[i][l])
    else:
        if uplo=='U':
            for j in range(n):
                for i in range(j-1):
                    temp = 0
                    for l in range(k):
                        temp += (conj(A[l][i])*A[l][j])
                    if beta==0:
                        C[i][j] = alpha*temp
                    else:
                        C[i][j] = alpha*temp + beta*C[i][j]
                rtemp = 0
                for l in range(k):
                    rtemp = rtemp + conj(A[l][j])*A[l][j]
        else:
            for j in range(n):
                rtemp = 0
                for l in range(k):
                    rtemp = rtemp + conj(A[l][j])*A[l][j]
                if beta==0:
                    C[j][j] = alpha*rtemp
                else:
                    C[j][j] = alpha*rtemp + beta*C[j][j].real
                for i in range(j+1, n):
                    temp = 0
                    for l in range(k):
                        temp = temp + conj(A[l][i])*A[l][j]
                    if beta==0:
                        C[i][j] = alpha*temp
                    else:
                        C[i][j] = alpha*temp + beta*C[i][j]
    
    return C
    
def conj(c):
    assert type(c) is complex
    return complex(c.real, (-1*c.imag))