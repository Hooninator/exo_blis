import numpy as np

# C = alpha*A*transpose(A) + beta*C
# C = alpha*transpose(A)*A + beta*C
def syrk(uplo: bool, trans: bool, alpha: float, A: list, beta: int, C: list):
    
    n = len(C)
    k = len(A[0]) 
    #transposed matrix dimensions
    n_t = len(A[0])
    k_t = len(C) 

    assert n==len(C[0]) # Make sure C is symmetrical 
    assert n >= 0
    assert k >= 0 

    #Quick return if certain parameters are zero
    if n==0 or alpha==0 or (k==0 and beta==1):
        return C 
    
    #Handle alpha=0 cases
    #TODO -- faster to just multiply by beta and remove the if statements?
    if alpha==0:
        #upper triangular case
        if uplo:
            if beta==0:
                for j in range(n):
                    for i in range(j):
                        C[i][j] = 0
            else:
                for j in range(n):
                    for i in range(j):
                        C[i][j] = beta*C[i][j]
        #lower triangular case
        else:
            if beta==0:
                for j in range(n):
                    for i in range(j, n):
                        C[i][j] = 0
            else:
                for j in range(n):
                    for i in range(j, n):
                        C[i][j] = beta*C[i][j]
    
    #Start the operation
    
    # C := alpha*A*transpose(A) + beta*C
    if trans:
        if uplo:
            for j in range(n):
                if beta==0:
                    for i in range(j):
                        C[i][j] = 0
                elif beta!=1:
                    C[i][j] = beta*C[i][j]
                for l in range(k):
                    for i in range(j):
                        C[i][j] = C[i][j] + alpha*A[i][l]*A[j][l]
        else:
            for j in range(n):
                if beta==0:
                    for i in range(j, n):
                        C[i][j] = 0
                elif beta!=1:
                    C[i][j] = beta*C[i][j]
                for l in range(k):
                    for i in range(j, n):
                        C[i][j] = C[i][j] + alpha*A[i][l]*A[j][l]

    # C := alpha*transpose(A)*A + beta*C
    else:
        if uplo:
            for j in range(n):
                for i in range(j):
                    temp=0
                    for l in range(k):
                        temp = temp + A[l][i]*A[l][j]
                    if beta==0:
                        C[i][j] = alpha*temp
                    else:
                        C[i][j] = alpha*temp + beta*C[i][j]
        else:
            for j in range(n):
                for i in range(j, n):
                    temp=0
                    for l in range(k):
                        temp = temp + A[l][i]*A[l][j]
                    if beta==0:
                        C[i][j] = alpha*temp
                    else:
                        C[i][j] = alpha*temp + beta*C[i][j]

    return C
