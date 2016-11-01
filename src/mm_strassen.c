/**
* @file mm_strassen.c
* @brief Matrix multiplication with Strassen algorithm. 
* @author Marcin Wolcendorf
* @version 0.1
* @date 2016-10-28
*/
#include <assert.h>
#include <stdio.h>
#include <stdint.h>

#define CHECK_POWER_OF_2(_x) (!((_x)&((_x)-1)))
#define LSB_BIT_FILL(x) ({\
  __auto_type _y = (x);     \
  (_y) |= (_y)>>1;          \
  (_y) |= (_y)>>2;          \
  (_y) |= (_y)>>4;          \
  (_y) |= (_y)>>8;          \
  (_y) |= (_y)>>16;         \
  _y;               \
})
#define LOG_BASE_2(x) ({ \
  __auto_type _y = (x); \
  _y = (_y & 0x55555555) + ((_y >> 1) & 0x55555555); \
  _y = (_y & 0x33333333) + ((_y >> 2) & 0x33333333); \
  _y = (_y & 0x0f0f0f0f) + ((_y >> 4) & 0x0f0f0f0f); \
  _y = (_y & 0x00ff00ff) + ((_y >> 8) & 0x00ff00ff); \
  _y = (_y & 0x0000ffff) + ((_y >> 16) & 0x0000ffff); \
  _y; \
})
#define HELPER_ARRAY_CNT 5u

size_t get_extra_array_size(uint32_t N) 
{
  size_t ret = 0;
  //normalize the dimension to be "2^n - 1" (MSBs only 0, LSBs only 1)
  ret = LSB_BIT_FILL(N-1);
  ++ret;                          //make it power of 2 (ret = 2^n)
  assert(CHECK_POWER_OF_2(ret));  //check that 
  ret *= ret;                     //geometric series with base 4, and the power the same as above - so square it (now ret = 4^n)
  --ret;                          //finish the numerator (4^n - 1)
  ret /= (4 - 1);                 //divide by the denominator ((4^n - 1)/(4 - 1))
  //now ret = SUM(i=0,n){4^i}
  return ret * HELPER_ARRAY_CNT; //we need 5 helper arrays. 
}

void mult_strassen_r(data_t *C, data_t *B, data_t *A, uint32_t N_limit, data_t *helper, uint32_t N_norm) 
{
  int32_t *M[HELPER_ARRAY_CNT];//helper arrays, so there is less pointer arithmetic
  uint32_t i;
  N_norm >>= 1;
  for(i=0; i<HELPER_ARRAY_CNT; i++) {
    M[i] = helper + i*N_norm; //there are 7 half-input-dimension helper arrays
  }

}



void matrix_copy(data_t *O, uint32_t O_dim, data_t *N, uint32_t dim) 
{
  uint32_t i, j;
  assert(N); 
  assert(O); 
  for (i=0; i<O_dim, i++) {
    for (j=0; j<O_dim, j++) {
      N[i*dim + j] = O[i*O_dim + j]; 
    }
  }
  return N;
}



/**
* @brief Staring the computations.
*
* This function will normalise the matrices, so they are 2^n x 2^n (as Strassen algorighm expects). Since the task specificaton requeries the matrices to be NxN only, there might be some expanding to be done. 
* This is an inefficiency of this implementation, for the sake of (relative) simplicity.
*
* @param C - output matrix,
* @param B - input matrix 2,
* @param A - input matrix 1,
* @param N - dimension (like N form NxN).
*
* @return No idea yet. 
*/
int32_t mult_strassen(data_t *C, const data_t *B, const data_t *A, const uint32_t N)
{
  data_t *LA = , *LB = B, *LC = C;
  uint32_t dim = N;
  size_t array_offset = 0;
  int32_t *helper;
  assert(A);
  assert(B);
  assert(C);
  assert(N);

  if(!CHECK_POWER_OF_2(dim)) {
    //we have some expanding to do - but only once
    //first - get the dimension to 2^n>N>2^(n-1)
    dim = LSB_BIT_FILL(dim - 1); 
    dim ++; 
    LA = calloc(dim * dim, sizeof(data_t));
    LB = calloc(dim * dim, sizeof(data_t));
    LC = calloc(dim * dim, sizeof(data_t));

    LA = matrix_copy(A, N, LA, dim);
    LB = matrix_copy(A, N, LB, dim);
    LA = matrix_copy(A, N, LC, dim);
  }
  //now we can split the matrices in 4 up to the point they are single data_t element


  //the memory is allocated once only and reused
  helper = calloc(get_extra_array_size(N), sizeof(A[0]));
  assert(helper);
  //do things here. 
  mult_strassen_r(C, B, A, N, helper, LSB_BIT_FILL(N-1)+1);

  free(helper);
  helper=NULL;
  return 0;
}


int main(int argc, char *argv[]) {
  int size = 8;
  printf("Ala ma asa\n");
  printf("Mem needed for array size %d = %d\n", size, get_extra_array_size(size));
  return 0;
}
