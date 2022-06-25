#NOTE: this is just hemm without complex numbers
#IDEA -- look through all of this code for chunks that can be cast as GEMM

def SYMM(A, B, C, alpha, beta, side, uplo):
    
    assert (A!=None)
    assert (B!=None)
    assert (C!=None)

    m = len(A) #rows of A and C
    n = len(B[0]) #cols of B and C

    #quick return if certain parameters are 0
    if m==0 or n==0 or (alpha==0 and beta==1):
        return C

    #Alpha is zero, so A*B becomes 0
    if alpha==0:
        if beta==0:
            #Zero C
            for j in range(n): 
                for i in range(m):
                    C[i][j] = 0 
        else:
            #Scale C by beta
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
                        temp2 = temp2+B[k][j]*(A[k][i])
                    if beta==0:
                        C[i][j] = temp1*A[i][i] + alpha*temp2
                    else:
                        C[i][j] = beta*C[i][j] + temp1*A[i][i] + alpha*temp2
        else:
            for j in range(n):
                for i in range(m, 1, -1):
                    temp1 = alpha*B[i][j]
                    temp2 = 0
                    for k in range(i+1, m):
                        C[k][j] = C[k][j] + temp1*A[k][i]
                        temp2 = temp2 + B[k][j]*(A[k][i])
                    if beta==0:
                        C[i][j] = temp1*A[i][i]+ alpha*temp2
                    else:
                        C[i][j] = beta*C[i][j] + temp1*A[i][i] + alpha*temp2
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
                    temp1 = alpha*(A[j][k])
                for i in range(m):
                    C[i][j] = C[i][j] + temp1*B[i][k]
            for k in range(j+1, n):
                if uplo=='U':
                    temp1 = alpha*(A[j][k])
                else:
                    temp1 = alpha*A[k][j]
                for i in range(m):
                    C[i][j] = C[i][j] + temp1*B[i][k]

    return C