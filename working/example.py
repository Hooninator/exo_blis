"""
"""
### Imports
from __future__ import annotations
from exo import proc
from exo.libs.memories import DRAM_STATIC

### Globals
K = 8 # Block size

### Body

"""
SPMV naive implementation
Writes output vector to B
"""
@proc
def SPMV(N: size, 
        A: R[N, N], 
        x: R[N],
        B: R[N]):
        for i in par(0, N):
                for j in par(0, N):
                        B[j] += A[j, i]*x[i]

SPMV_WINDOW = (SPMV.rename("SPMV_WINDOW")
                .set_window('A', True)
                .set_window('x', True)
                .set_window('B', True))


SPMV_TEST = (SPMV_WINDOW.rename("SPMV_TEST")
                .stage_assn('B_reg', 'B[_] += _'))

SPMV_TEST = SPMV_TEST.lift_alloc('B_reg: _') 


#SPMV_WINDOW = SPMV_WINDOW.stage_window('A_cache', 'A[_] #0', DRAM_STATIC)


SPMV_SPLIT = (SPMV.rename("SPMV_SPLIT")
                .split('i', 16, ['io','ii'], tail='cut_and_guard')
                .split('j', 16, ['jo', 'jj'], tail='cut'))

print(SPMV_SPLIT.c_code_str())
