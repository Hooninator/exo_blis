import numpy as np #This is only done for testing purposes 

def trsv(upper: bool, lower: bool, transpose: bool,
         diag: bool, diag_unit: bool, A: list, x: list, b: list):
    
    assert A and x and b
    assert not (upper and lower) #Make sure only one case holds true

    m = len(A)
    n = len(A[0])
    if transpose:
        A = transpose(A)

    #Upper triangular matrix
    if upper or (lower and transpose):
        for i in range(m-1, -1, -1):
            temp_b = b[i]
            for j in range(n-1, -1, -1):
                if j==i:
                    x[i] = temp_b/A[i][j]
                    continue
                else:
                    temp_b -= A[i][j]*x[j]
        return x

    #Lower triangular matrix
    elif lower or (upper and transpose):
        for i in range(m):
            temp_b = b[i]
            for j in range(n):
                if j==i:
                    x[i] = temp_b/A[i][j]
                    continue
                else:
                    temp_b -= A[i][j]*x[j]
        return x

    #Diagonal matrix
    elif diag:
        for i in range(n):
            if diag_unit:
                x[i] = b[i]
            else:
                x[i] = b[i]/A[i][i]
        return x


def transpose(M: list) -> list:
    M_t = list(np.zeros((len(M[0]), len(M))))
    for i in range(len(M)):
        for j in range(len(M[i])):
            M_t[j][i] = M[i][j]
    return M_t


def is_correct(a: float, b:float):
    return abs(a-b)<0.0001


def test_upper():
    A = list(np.triu(np.random.rand(32, 32)))
    b = list(np.random.rand(32))
    x = list(np.zeros((32)))
    x_np = np.linalg.solve(A, b)
    x_man = trsv(True, False, False, False, False, A, x, b)
    for i in range(len(x)):
        if not is_correct(x_np[i], x_man[i]):
            raise Exception(f"Error: Element {i} should be {x_np[i]}\nGot{x_man[i]}")
    print("Test upper passed!")


def test_lower():
    A = list(np.tril(np.random.rand(32, 32)))
    b = list(np.random.rand(32))
    x = list(np.zeros((32)))
    x_np = np.linalg.solve(A, b)
    x_man = trsv(False, True, False, False, False, A, x, b)
    for i in range(len(x)):
        if not is_correct(x_np[i], x_man[i]):
            raise Exception(f"Error: Element {i} should be {x_np[i]}\nGot {x_man[i]}")
    print("Test lower passed!")


def test_diag_unit():
    A = list(np.identity(32))
    b = list(np.random.rand(32))
    x = list(np.zeros((32)))
    x_np = np.linalg.solve(A, b)
    x_man = trsv(False, False, False, True, True, A, x, b)
    for i in range(len(x)):
        if not is_correct(x_np[i], x_man[i]):
            raise Exception(f"Error: Element {i} should be {x_np[i]}\nGot {x_man[i]}")
    print("Test diag unit passed!")


def test_diag_nonunit():
    A = list(np.zeros((32, 32)))
    for i in range(32):
            A[i][i] = np.random.rand()
    b = list(np.random.rand(32))
    x = list(np.zeros((32)))
    x_np = np.linalg.solve(A, b)
    x_man = trsv(False, False, False, True, False, A, x, b)
    for i in range(len(x)):
        if not is_correct(x_np[i], x_man[i]):
            raise Exception(f"Error: Element {i} should be {x_np[i]}\nGot {x_man[i]}")
    print("Test diag non-unit passed!")


def main():
    test_upper()
    test_lower()
    test_diag_unit()
    test_diag_nonunit()


if __name__=="__main__":
    main()