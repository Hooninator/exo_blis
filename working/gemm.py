"""

"""

import numpy as np #This is only done for testing purposes 

def gemm(transpose_A: bool, transpose_B: bool, alpha: float, A: list,
        B: list, beta: float, C: list):

        assert (A!=None)
        assert (B!=None)
        assert (C!=None)

        m = len(A) #rows of A and C
        n = len(B[0]) #cols of B and C
        k = len(B) #rows of B, cols of A
        
        
        #quick return if certain parameters are 0
        if m==0 or n==0 or ((k==0 or alpha==0) and beta==1):
            return C

        #Alpha is zero, so A*B becomes 0
        if alpha==0:
            if beta==0:
                #Zero C
                for i in range(m): 
                    for j in range(n):
                        C[i][j] = 0 
            else:
                #Scale C by beta
                for i in range(m):
                    for j in range(n):
                        C[i][j]*=beta
            return C
        
        #Tranpose A or B
        #NOTE: BLAS appears to just handle this by indexing the matrices differently, but BLIS handles it by actually calling a transpose operation
        if transpose_A:
            A = transpose(A)
            #Swap dims
            k = len(A[0])
            m = len(A)
        if transpose_B:
            B = transpose(B)
            #Swap dims
            k = len(B)
            n = len(B[0])
        
        #Do the operation
        for i in range(m):
            for j in range(n):
                for l in range(k):
                    C[i][j] += alpha*A[i][l]*B[l][j]
        return C


def transpose(M: list) -> list:
    M_t = list(np.zeros((len(M[0]), len(M))))
    for i in range(len(M)):
        for j in range(len(M[i])):
            M_t[j][i] = M[i][j]
    return M_t


def test_square():
    n = 32
    m = 32
    k = 32
    A = list(np.random.rand(m, k))
    B = list(np.random.rand(k, n))
    C = list(np.zeros((m, n)))
    C_manual = gemm(False, False, 1, A, B, 1, C)
    C_np = np.dot(A, B)

    for i in range(m):
        for j in range(n):
            if round(C_manual[i][j], 8)!=round(C_np[i][j], 8):
                raise Exception(f"Error: {i, j} should be {C_np[i][j]}\nGot {C_manual[i][j]}")
    print("Square Matrix test passed!")


def test_transpose():
    n = 32
    m = 32
    k = 16

    A = list(np.random.rand(m, k))
    B = list(np.random.rand(k, n))
    C = list(np.zeros((m, n)))
    C_manual = gemm(True, True, 1, A, B, 1, C)
    C_np = np.dot(np.transpose(A), np.transpose(B))

    for i in range(k):
        for j in range(k):
            if round(C_manual[i][j], 8)!=round(C_np[i][j], 8):
                raise Exception(f"Error: {i, j} should be {C_np[i][j]}\nGot {C_manual[i][j]}")
    print("Transposed Matrix test passed!")


def main():
    #Test simple square 32x32 matmul
    test_square()
    #Test transposed 32*16 matmul
    test_transpose()
    

if __name__=="__main__":
    main()