import numpy as np

def syrk(trans: bool, alpha: float, A: list, beta: int, C: list):
    
    n = len(C)
    k = len(A[0]) if trans else len(A) #k is the number of columns in A or the number of rows if A is transposed

    assert n==len(C[0]) # Make sure C is symmetrical 
    assert n >= 0
    assert k >= 0 

    #Quick return if certain parameters are zero
    if n==0 or alpha==0 or (k==0 and beta==1):
        return C 
    
    #Handle alpha=0 cases
    if alpha==0:
        

    return C



def transpose(M: list) -> list:
    M_t = list(np.zeros((len(M[0]), len(M))))
    for i in range(len(M)):
        for j in range(len(M[i])):
            M_t[j][i] = M[i][j]
    return M_t