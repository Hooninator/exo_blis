

def HER2K(A, B, C, alpha, beta, transpose, uplo):

    assert A!=None
    assert B!=None
    assert C!=None
    assert uplo=='U' or uplo=='L'

    n = len(C)
    k = len(A) if transpose else len(A[0])

    if n==0 or ((alpha==0 or k==0) and beta==1):
        return C
    
    if alpha==0:
        if uplo=='U':
            if beta.real==0:
                for j in range(n):
                    for i in range(j):
                        C[i][j] = 0
            else:
                for j in range(n):
                    for i in range(j-1):
                        C[i][j] *= beta
                    C[j][j] = beta*C[j][j].real
        else:
            if beta.real==0:
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
        if uplo=='U':
            for j in range(n):
                if beta.real==0:
                    for i in range(j):
                        C[i][j] = 0
                elif beta!=1:
                    for i in range(j-1):
                        C[i][j] = beta*C[i][j]
                    C[j][j] = beta*C[j][j].real
                else:
                    C[j][j] = C[j][j].real
                for l in range(k):
                    if A[j][l]!=0 or B[j][l]!=0:
                        temp1 = alpha*conj(B[j][l])
                        temp2 = conj(alpha*A[j][l])
                        for i in range(j-1):
                            C[i][j] = C[i][j] + A[i][l]*temp1 + B[i][l]*temp2
                        C[j][j] = C[j][j].real + (A[j][l]*temp1+B[j][l]*temp2).real
        else:
            for j in range(n):
                if beta.real==0:
                    for i in range(j, n):
                        C[i][j] = 0
                elif beta!=1:
                    #NOTE: this combo is in reversed order in herk
                    for i in range(j+1, n):
                        C[i][j] = beta*C[i][j]
                    C[j][j] = beta*C[j][j].real
                else:
                    C[j][j] = C[j][j].real
                for l in range(k):
                    if A[j][l]!=0 or B[j][l]!=0:
                        temp1 = alpha*conj(B[j][l])
                        temp2 = conj(alpha*A[j][l])
                        for i in range(j+1, n):
                            C[i][j] = C[i][j] + A[i][l]*temp1 + B[i][l]*temp2
    else:
        if uplo=='U':
            for j in range(n):
                for i in range(j):
                    temp1 = temp2 = 0
                    for l in range(k):
                        temp1 += conj(A[l][i])*B[l][j]
                        temp2 += conj(B[l][i])*A[l][j]
                    if i==j:
                        if beta.real==0:
                            C[j][j] = (alpha*temp1 + conj(alpha)*temp2).real
                        else:
                            C[j][j] = beta*C[j][j].real + (alpha*temp1 + conj(alpha)*temp2).real
                    else:
                        if beta.real==0:
                            C[i][j] = alpha*temp1 + conj(alpha)*temp2
                        else:
                            C[i][j] = beta*C[i][j] + alpha*temp1 + conj(alpha)*temp2
        else: #Again, only difference here seems to be in the i loop, think theres a way in exo to like replace part of a procedure with another
            for j in range(n):
                for i in range(j, n):
                    temp1=temp2=0
                    for l in range(k):
                        temp1 += conj(A[l][i])*B[l][j]
                        temp2 += conj(B[l][i])*A[l][j]
                    if i==j:
                        if beta.real==0:
                            C[j][j] = (alpha*temp1 + conj(alpha)*temp2).real
                        else:
                            C[j][j] = beta*C[j][j].real + (alpha*temp1 + conj(alpha)*temp2).real
                    else:
                        if beta.real==0:
                            C[i][j] = alpha*temp1 + conj(alpha)*temp2
                        else:
                            C[i][j] = beta*C[i][j] + alpha*temp1 + conj(alpha)*temp2
    return C
                

def conj(c):
    assert type(c) is complex
    return complex(c.real, (-1*c.imag))