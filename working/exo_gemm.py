from __future__ import annotations
from exo import *
from exo.libs.memories import DRAM_STATIC
from exo.platforms.x86 import *
from exo.syntax import *

###Define kernel constants
VEC_W = 16

M_REG_BLK = 6
N_REG_BLK = (4 * VEC_W)

M_L1_FAC = 44
N_L1_FAC = 1

M_L1_BLK = M_REG_BLK * M_L1_FAC
N_L1_BLK = N_REG_BLK * N_L1_FAC
K_L1_BLK = 512

@proc
def SGEMM(M: size, N: size, K: size, A: f32[M, K], B: f32[K, N], C: f32[M, N]):
    assert M >= 1
    assert N >= 1
    assert K >= 1
    assert stride(A, 1) == 1
    assert stride(B, 1) == 1
    assert stride(C, 1) == 1

    for k in par(0, K):
        for i in par(0, M):
            for j in par(0, N):
                C[i, j] += A[i, k]*B[k, j]


###Create the microkernel
gemm_microkernel = (SGEMM
                        .rename('gemm_microkernel')
                        .partial_eval())