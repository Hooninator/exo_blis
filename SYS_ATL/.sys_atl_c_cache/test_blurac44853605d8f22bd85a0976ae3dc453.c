int _ceil_div(int num, int quot) {
  int off = (num>0)? quot-1 : 0;
  return (num+off)/quot;
}


// test_blur( image : R[N,M] @IN, kernel : R[K,K] @IN, res : R[N,M] @OUT )
void test_blur( int N, int M, int K, float* image, float* kernel, float* res) {
for (int x=0; x < N; x++) {
  for (int y=0; y < M; y++) {
    res[(x) * M + (y)] = 0.0;
  }
}
for (int xhi=0; xhi < _ceil_div(N - 4 + 1, 4); xhi++) {
  for (int yhi=0; yhi < _ceil_div(M - 4 + 1, 4); yhi++) {
    for (int i=0; i < K; i++) {
      for (int j=0; j < K; j++) {
        for (int xlo=0; xlo < 4; xlo++) {
          for (int ylo=0; ylo < 4; ylo++) {
            if (4 * xhi + xlo + i < N && 4 * yhi + ylo + j < M) {
              res[(4 * xhi + xlo) * M + (4 * yhi + ylo)] += kernel[(i) * K + (j)] * image[(4 * xhi + xlo + i) * M + (4 * yhi + ylo + j)];
            }
          }
        }
      }
    }
  }
  for (int xlo=0; xlo < 4; xlo++) {
    for (int ylo=0; ylo < M - _ceil_div(M - 4 + 1, 4) * 4; ylo++) {
      for (int i=0; i < K; i++) {
        for (int j=0; j < K; j++) {
          if (4 * xhi + xlo + i < N && ylo + _ceil_div(M - 4 + 1, 4) * 4 + j < M) {
            res[(4 * xhi + xlo) * M + (ylo + _ceil_div(M - 4 + 1, 4) * 4)] += kernel[(i) * K + (j)] * image[(4 * xhi + xlo + i) * M + (ylo + _ceil_div(M - 4 + 1, 4) * 4 + j)];
          }
        }
      }
    }
  }
}
for (int xlo=0; xlo < N - _ceil_div(N - 4 + 1, 4) * 4; xlo++) {
  for (int y=0; y < M; y++) {
    for (int i=0; i < K; i++) {
      for (int j=0; j < K; j++) {
        if (xlo + _ceil_div(N - 4 + 1, 4) * 4 + i < N && y + j < M) {
          res[(xlo + _ceil_div(N - 4 + 1, 4) * 4) * M + (y)] += kernel[(i) * K + (j)] * image[(xlo + _ceil_div(N - 4 + 1, 4) * 4 + i) * M + (y + j)];
        }
      }
    }
  }
}
}