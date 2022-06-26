def T(X):
    rows = len(X)
    cols = len(X[0])

    X_T = []
    for c in range(cols):
        flipCol = []
        for r in range(rows):
            flipCol.append(X[r][c])
        X_T.append(flipCol)

    return X_T

def scale(X, s):
    rows = len(X)
    cols = len(X[0])
    for r in range(rows):
        for c in range(cols):
            X[r][c] *= s
    return X

def add(X, Y):
    rows = len(X)
    cols = len(X[0])
    Z = []
    for r in range(rows):
        row = []
        for c in range(cols):
            row.append(X[r][c] + Y[r][c])
        Z.append(row)  
    return Z

def gemm():
    return None

#params = {UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC}
def SYRK(A, C, a, b):
    '''
    Case 1: C = alpha * A * A^T + beta * C 
    Case 2: C = alpha * A^T * A + beta * C 

    C is n by n 
    A is n by k for case 1
    A is k by n for case 2 

    Params
    UPLO - Specifies upper or lower triangular spot
    TRANS - N = AAT , T or C = ATA 
    N - order of matrix C
    K - # of row or column depending on the transpose
    Alpha -  a scalar
    A - a array of real values
    LDA - first dimension of A
    C - a array of real values
    LDA - first dimension of C 
    '''
    #no numpy
    rows = len(A)
    cols = len(A[0])
    if len(C) == len(A):
        A = A
        #perform A*AT
    elif len(C) == len(A[0]):
        A = T(A)

    newMatrix = []
    for row in A: 
        # for each row in A
        # multiply that by the last to first row of A, and include those values as the first new row
        newrow = []
        for r in range(rows):
            sum = 0
            for c in range(cols):
                sum += row[c] * A[rows-1-r][c]
            newrow.append(sum)
        newMatrix.append(newrow)


    C = add(scale(newMatrix, a), scale(C, b))

    return C


#params = {alpha, lda, ldb, M, N, diag, side, transa, uplo, ALDA, BLDB}
def TRSM(alpha, lda, ldb, M, N, diag, side, transa, uplo, ALDA, BLDB):
    '''
    http://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_ga9893cceb3ffc7ce400eee405970191b3.html
    
    Defintion: 
    X*A = alpha*B
    A*X = alpha*B

    A is an upper or lower triangular matrix
    X overwrites B 

    Functions: 
    lsame, xerbla, max

    Scalars:
    temp, i, info, j, k, nrowa, lside, nountit, upper

    Parameters: 
    one, zero
    '''
    

    return None

def TRSV():
    return None


