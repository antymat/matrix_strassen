/**
* @file mm_classic.c
* @brief Classic matrix multiplication (mostly for comparison). 
* @author Marcin Wolcendorf
* @version 0.1
* @date 2016-10-28
*/



#include "mm_classic.h"

/**
* @brief Classical matrix multiplication.
*
* @param C - output matrix,
* @param B - multiplier,
* @param A - multiplicand,
* @param N - dimension (like N form NxN).
*/
void mult_classic(data_t *C, data_t *B, data_t *A, const uint32_t N)
{
  uint32_t i,j,k;
  assert(A);
  assert(B);
  assert(C);
  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      C[i*N+j] = 0;
      for(k=0; k<N; k++) {
        C[i*n+j] += A[i*N+k] * B[k*N+j];
      }
    }
  }
}
