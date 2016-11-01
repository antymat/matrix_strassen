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

typedef data_t int32_t; 

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
  return ret * HELPER_ARRAY_CNT; //we need some helper arrays. 
}


void matrix_sub(data_t *result, uint32_t r_rowlen_log, 
    data_t *B, uint32_t B_rowlen_log, 
    data_t *A, uint32_t A_rowlen_log, 
    uint32_t dim_log) 
{
  uint32_t i, j; 
  for (i=0; i < 1UL<<dim_log; i++) {
    for(j=0; j < 1UL<<dim_log; j++) {
      result[i * r_rowlen_log + j] = A[i * A_rowlen_log + j] - B[i * B_rowlen_log + j];
    }
  }
}

void matrix_add(data_t *result, uint32_t r_rowlen_log, 
    data_t *B, uint32_t B_rowlen_log, 
    data_t *A, uint32_t A_rowlen_log, 
    uint32_t dim_log) 
{
  uint32_t i, j; 
  for (i=0; i < 1UL<<dim_log; i++) {
    for(j=0; j < 1UL<<dim_log; j++) {
      result[i * r_rowlen_log + j] = A[i * A_rowlen_log + j] + B[i * B_rowlen_log + j];
    }
  }
}

void mult_strassen_r(data_t *C, uint32_t C_rowlen_log, 
    data_t *B, uint32_t B_rowlen_log, 
    data_t *A, uint32_t A_rowlen_log, 
    uint32_t dim_log, data_t **helper) 
{
  data_t *M[HELPER_ARRAY_CNT];//helper arrays, so there is less pointer arithmetic
  uint32_t i;
  --dim_log;
  if(!dim_log) {
    //the matrices have degenerated to the single data_t elements
    //this is the easy case - just multiply them, and
    //end the recursion
    *C = *A * *B;
    return;
  }
  //this is the hard case

  //ease the pointer operations a bit
  for(i=0; i<HELPER_ARRAY_CNT; i++) {
    M[i] = helper[dim_log] + i*(1UL<<dim_log);
  } 


}



void matrix_copy(data_t *N, uint32_t N_dim, data_t *O, uint32_t O_dim) 
{
  uint32_t i, j;
  uint32_t dim;
  assert(N); 
  assert(O); 
  dim  = O_dim < N_dim ? O_dim : N_dim; 
  for (i=0; i<dim, i++) {
    for (j=0; j<dim, j++) {
      N[i*N_dim + j] = O[i*O_dim + j]; 
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
  uint32_t dim = N; //this will be the 2^n dimension value. 
  uint32_t i; 
  uint32_t dim_log; 
  size_t array_offset = 0;
  data_t **helper;
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

    LA = matrix_copy(LA, dim, A, N);
    LB = matrix_copy(LB, dim, B, N);
    LA = matrix_copy(LC, dim, C, N);
  }
  //now we can split the matrices in 4 up to the point they are single data_t element


  //now we create all the helper matrices we will need in future - they are of the same size, and there are LOG_2(dim) sets of them. 
  dim_log = LOG_BASE_2(dim - 1);
  helper = malloc(dim_log * sizeof(data_t*));
  assert(helper);
  for(i = dim_log; i;) {
    --i;
    helper[i] = malloc(sizeof(data_t) * (1UL<<i) * (1UL << i)i * HELPER_ARRAY_CNT);
    assert(helper[i]);
  }
  //do things here. 
  mult_strassen_r(LC, dim_log, LB, dim_log, LA, dim_log, dim_log, helper);



  //copy the result back to C
  matrix_copy(C, N, LC, dim);

  //free the memory
  for(i = dim_log; i;) {
    free(helper[i]);
  }
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
