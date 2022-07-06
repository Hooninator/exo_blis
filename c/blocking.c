/* Trying to cast other level 3 BLAS operations in terms of GEMM and special blocking procedures */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define BLOCKSIZE 32
#define REGSIZE 4
#define RAND_HIGH 20
#define RAND_LOW -20


void gemm_pack_A(double *A, double *A_packed, uint32_t M_c, uint32_t K_c, uint32_t M_r) {
    for (int i=0; i<M_c; i+=M_r) {
        for (int j=0; j<K_c; ++j) {
            for (int k=0; k<M_r; ++k) {
                A_packed[k + j*M_r + i*K_c] = A[k + M_c*j + i];
            }
        }
    }
}

void pack_A(double *A, double *A_packed, uint32_t M, uint32_t M_r,  uint32_t M_c, uint32_t K_c) {
	for (int i=0; i<M_c; i+=M_r) {
		for (int j=0; j<K_c; ++j) {
			for (int k=0; k<M_r; ++k) {
				A_packed[k + j*M_r + i*K_c] = A[k + M*j  + i];
			}
		}
	}
}


void gemm_pack_B(double *B, double *B_packed, uint32_t N, uint32_t K_c, uint32_t N_r) {
    for (int i=0; i<N; i+=N_r) {
        for (int j=0; j<K_c; ++j) {
            for (int k=0; k<N_r; ++k) {
                B_packed[k + j*N_r + i*K_c] = B[k*K_c + j + i*N];
            }
        }
    }
}

void pack_B(double *B, double *B_packed, uint32_t N, uint32_t K, uint32_t K_c, uint32_t N_r) {
    for (int i=0; i<N; i+=N_r) {
        for (int j=0; j<K_c; ++j) {
            for (int k=0; k<N_r; ++k) {
                B_packed[k + j*N_r + i*K_c] = B[k*K + j + i*K];
            }
        }
    }
}


void gemm_microkernel(double *A, double *B, double *C, uint32_t K_c) {

    for (int i = 0; i < K_c; ++i) {
		for (int j = 0; j < REGSIZE; ++j) {
			for (int k = 0; k < REGSIZE; ++k) {
				C[k + j*REGSIZE] += A[k + i*REGSIZE]*B[j + i*REGSIZE];
			}
		}
	}

}


void GEBP(double *A, double *B, double *C, uint32_t M_c, uint32_t N_c, uint32_t K_c) {

    double *A_c;
    double *B_c;
    double *C_c;
    
    /* Call the microkernel */
    for (int i=0; i<M_c; i+=REGSIZE) {
        A_c = A + (K_c*i); //TODO, ADD parameters to microkernel that let it jump around in C
        for (int j=0; j<N_c; j+=REGSIZE) {
            B_c = B + (K_c*j);
            for (int k=0; k<K_c; k+=REGSIZE) {
                C_c = C + (M_c*k);
                gemm_microkernel(A_c, B_c, C_c, K_c);
            }
        }
    }
}


void GEPP(double *A, double *B, double *C, uint32_t M, uint32_t N, uint32_t K_c, uint32_t M_c) {
    /* Pack B */
    double B_p[N*K_c];
    gemm_pack_B(B, B_p, N, K_c, REGSIZE);

    for (int i=0; i<M; i+=M_c) {
        double A_p[M_c*K_c];
        gemm_pack_A(A, A_p, M_c, K_c, REGSIZE);
        GEBP(A_p, B_p, C, M_c, N, K_c);
    }
}

/* GEMM Microkernel */
void microkernel(double *A, double *B, double *C, uint32_t M, uint32_t K) {
	for (int i=0; i<REGSIZE; ++i) {
		for (int j=0; j<REGSIZE; ++j) {
			for (int k=0; k<BLOCKSIZE; ++k) {
				//TODO, remember it is packed
				C[k + j*M] = A[k + i*REGSIZE]*B[j + i*REGSIZE]; 
			}
		}
	}
}

/* GEMM */
void GEMM(double *A, double *B, double *C, uint32_t M, int N, uint32_t K) {
	
	/* p loop */
	for (int p=0; p<K; p+=BLOCKSIZE) {
		double *A_panel = A+(p*M);
		double *B_panel = B+(p);
		/* GEPP */
		for (int i=0; i<M; i+=BLOCKSIZE) {
			double *C_panel = C+(i);
			double *B_packed = malloc(sizeof(double)*N*BLOCKSIZE);
			pack_B(B_panel, B_packed, N, K, BLOCKSIZE, REGSIZE);
			double *A_blk = A_panel+(i);
			/* GEBP */
			for (int jr=0; jr<N; jr+=REGSIZE) {
				double *C_strip = C_panel+(jr*M);
				double *B_reg = B_packed+(jr*BLOCKSIZE); //TODO: Change this since it is packed 
				double *A_packed = malloc(sizeof(double)*BLOCKSIZE*BLOCKSIZE);
				pack_A(A_blk, A_packed, M, REGSIZE, BLOCKSIZE, BLOCKSIZE);
				/* Microkernel */
				for (int ir=0; ir<BLOCKSIZE; ir+=REGSIZE) {
					double *C_reg = C_strip+ir;
					double *A_reg = A_packed+ir; //TODO: Change this since it is packed
					microkernel(A_reg, B_reg, C_reg, M, K);	
				}
			}
		}
	}
}


/* Returns a random double between RAND_UPPER and RAND_LOWER, inclusive */
double rand_dbl() {
	double d = (double) rand() / ((double) RAND_MAX+1);
	return (RAND_LOW + d * (RAND_HIGH - RAND_LOW));
}


void naive_matmul(double *A, double *B, double *C, int m, int n, int l) {
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			for (int k = 0; k < l; ++k) {
                printf("%lf + (%lf*%lf)\n", C[i + j*m], A[i + k*m], B[k + j*n]);
				C[i + j*m] += A[i + k*m]*B[k + j*l];
			}
		}
	}
}


/* Initializes a n*n matrix with random doubles from RAND_UPPER to RAND_LOWER 
 * -- the ranges are inclusive.
 */
void init_matrix(int n, int m, double *M) {
	for (int i=0; i < n*m; ++i) {
		M[i] = rand_dbl();
	}
}


void write_matrix(FILE *f, double *M, int m, int n) {
	fprintf(f, "Matrix.\n");
	fprintf(f, "--------\n");
	fprintf(f, "[");
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			fprintf(f, "%lf, ", M[j*m+i]);
		}
		fprintf(f, "\n");
	}
	fprintf(f, "]\n");
	fprintf(f, "--------\n");
}


void test_gebp() {

}


void test_packing() {
    double *A = malloc(sizeof(double)*8*8);
    double *B = malloc(sizeof(double)*8*8);
    init_matrix(8, 8, A);
    init_matrix(8, 8, B);

    double *A_p = malloc(sizeof(double)*8*8);
    double *B_p = malloc(sizeof(double)*8*8);
    gemm_pack_A(A, A_p, 8, 8, REGSIZE);
    gemm_pack_B(B, B_p, 8, 8, REGSIZE);
    char c = 'w';
    FILE *f = fopen("outputA.out", &c);
    write_matrix(f, A, 8, 8);
    write_matrix(f, A_p, 8, 8);
    fclose(f);
    f = fopen("outputB.out", &c);
    write_matrix(f, B, 8, 8);
    write_matrix(f, B_p, 8, 8);
    fclose(f);
    return;

}


void test_microkernel() {
    int n = 16;
    double *A = malloc(sizeof(double)*n*REGSIZE);
    double *B = malloc(sizeof(double)*n*REGSIZE);
    double *A_p = malloc(sizeof(double)*n*REGSIZE);
    double *B_p = malloc(sizeof(double)*n*REGSIZE);
    double *C_n = malloc(sizeof(double)*REGSIZE*REGSIZE);
    double *C_m = malloc(sizeof(double)*REGSIZE*REGSIZE);

    init_matrix(REGSIZE, n, A);
    init_matrix(n, REGSIZE, B);
    gemm_pack_A(A, A_p, REGSIZE, n, REGSIZE);
    gemm_pack_B(B, B_p, REGSIZE, n, REGSIZE);

    naive_matmul(A, B, C_n, REGSIZE, REGSIZE, n);
    printf("==============\n");
   // GEMM_microkernel(A_p, B_p, C_m, n);

    char c = 'w';
    FILE *f = fopen("micro.out", &c);
    write_matrix(f, C_n, REGSIZE, REGSIZE);
    write_matrix(f, C_m, REGSIZE, REGSIZE);
    fclose(f);

    f = fopen("outputA.out", &c);
    write_matrix(f, A, REGSIZE, n);
    write_matrix(f, A_p, REGSIZE, n);
    fclose(f);

    f = fopen("outputB.out", &c);
    write_matrix(f, B, n, REGSIZE);
    write_matrix(f, B_p, n, REGSIZE);
    fclose(f);

    free(A);
    free(B);
    free(A_p);
    free(B_p);
    free(C_n);
    free(C_m);
}


int main(int argc, char **argv) {
    test_microkernel();
}
