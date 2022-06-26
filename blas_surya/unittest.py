from blas import *
from scipy.linalg import blas


X1 = [[1, 2, 3], [1, 2, 3], [1, 2, 3]]
X2 = [[1, 2], [1, 2]]

print("Add:", add(X1, X1))
print("Scale:", scale(X1, 2))
print("SYRK:", SYRK(X1, X1, 2, 2))

print("DSYRK:", blas.dsyrk(2, X1))