/**
* @file mm_strassen.c
* @brief Matrix multiplication with Strassen algorithm.
* @author Marcin Wolcendorf
* @version 0.1
* @date 2016-10-28
*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "mm_strassen.h"


void print_matrix(uint8_t *name, data_t *A, uint32_t row_len, uint32_t dim) 
{
  uint32_t i, j;
  assert(row_len >= dim);
  printf("%s = \n", name);
  for(i=0; i<dim; i++) {
    for(j=0; j<dim; j++) {
      printf("%d\t", (uint32_t)A[i*row_len+j]);
    }
    printf("\n");
  }
}

size_t get_extra_array_size(uint32_t N)
{
  size_t ret = 0;
  //normalize the dimension to be "2^n - 1" (MSBs only 0, LSBs only 1)
  LSB_BIT_FILL(ret, N-1);
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
  uint32_t A_rowlen = 1UL << A_rowlen_log;
  uint32_t B_rowlen = 1UL << B_rowlen_log;
  uint32_t r_rowlen = 1UL << r_rowlen_log;
  for (i=0; i < 1UL<<dim_log; i++) {
    for(j=0; j < 1UL<<dim_log; j++) {
      result[i * r_rowlen + j] = A[i * A_rowlen + j] - B[i * B_rowlen + j];
    }
  }
}

void matrix_add(data_t *result, uint32_t r_rowlen_log,
    data_t *B, uint32_t B_rowlen_log,
    data_t *A, uint32_t A_rowlen_log,
    uint32_t dim_log)
{
  uint32_t i, j;
  uint32_t A_rowlen = 1UL << A_rowlen_log;
  uint32_t B_rowlen = 1UL << B_rowlen_log;
  uint32_t r_rowlen = 1UL << r_rowlen_log;
  for (i=0; i < 1UL<<dim_log; i++) {
    for(j=0; j < 1UL<<dim_log; j++) {
      result[i * r_rowlen + j] = A[i * A_rowlen + j] + B[i * B_rowlen + j];
    }
  }
}

void mult_strassen_r(data_t *C, uint32_t C_rowlen_log,
    data_t *B, uint32_t B_rowlen_log,
    data_t *A, uint32_t A_rowlen_log,
    uint32_t dim_log, data_t **helper)
{

  if(!dim_log) {
    //the matrices have degenerated to the single data_t elements
    //this is the easy case - just multiply them, and
    //end the recursion
    *C = *A * *B;
    return;
  }
  --dim_log; //divide the matrices into 4
  {
    //this is the hard case
    data_t *H[HELPER_ARRAY_CNT];//helper arrays, so there is less pointer arithmetic
    uint32_t dim = 1UL << dim_log;
    uint32_t A_rowlen = 1UL << A_rowlen_log;
    uint32_t B_rowlen = 1UL << B_rowlen_log;
    uint32_t C_rowlen = 1UL << C_rowlen_log;
    uint32_t i;
    //A
    data_t *A11 = MATRIX_PARTITION(A, A_rowlen, 0, 0);
    data_t *A12 = MATRIX_PARTITION(A, A_rowlen, dim, 0);
    data_t *A21 = MATRIX_PARTITION(A, A_rowlen, 0, dim);
    data_t *A22 = MATRIX_PARTITION(A, A_rowlen, dim, dim);
    // B
    data_t *B11 = MATRIX_PARTITION(B, B_rowlen, 0, 0);
    data_t *B12 = MATRIX_PARTITION(B, B_rowlen, dim, 0);
    data_t *B21 = MATRIX_PARTITION(B, B_rowlen, 0, dim);
    data_t *B22 = MATRIX_PARTITION(B, B_rowlen, dim, dim);
    // C
    data_t *C11 = MATRIX_PARTITION(C, C_rowlen, 0, 0);
    data_t *C12 = MATRIX_PARTITION(C, C_rowlen, dim, 0);
    data_t *C21 = MATRIX_PARTITION(C, C_rowlen, 0, dim);
    data_t *C22 = MATRIX_PARTITION(C, C_rowlen, dim, dim);

    data_t *M1, *M2, *M3, *M4, *M5, *M6, *M7;
    uint32_t M1_rowlen_log, M2_rowlen_log, M3_rowlen_log, M4_rowlen_log, M5_rowlen_log, M6_rowlen_log, M7_rowlen_log;


    //ease the pointer operations a bit
    for(i=0; i<HELPER_ARRAY_CNT; i++) {
      H[i] = helper[dim_log] + i*dim*dim;
    }

    M1 = H[0];	  M1_rowlen_log = dim_log;
    M2 = C21;	    M2_rowlen_log = C_rowlen_log;
    M3 = H[1];	  M3_rowlen_log = dim_log;
    M4 = H[2];	  M4_rowlen_log = dim_log;
    M5 = C12;	    M5_rowlen_log = C_rowlen_log;
    M6 = C22;	    M6_rowlen_log = C_rowlen_log;
    M7 = C11;	    M7_rowlen_log = C_rowlen_log;

    {
      //calculate M1
      data_t *T1 = H[3];
      data_t *T2 = H[4];
      matrix_add(T1, dim_log,
          A22, A_rowlen_log,
          A11, A_rowlen_log,
          dim_log);
      matrix_add(T2, dim_log,
          B22, B_rowlen_log,
          B11, B_rowlen_log,
          dim_log);
      mult_strassen_r(M1, M1_rowlen_log,
          T2, dim_log,
          T1, dim_log,
          dim_log, helper);
    }
    {
      //calcualte M2
      data_t *T1 = H[3];
      matrix_add(T1, dim_log,
          A22, A_rowlen_log,
          A21, A_rowlen_log,
          dim_log);
      mult_strassen_r(M2, M2_rowlen_log,
          B11, B_rowlen_log,
          T1, dim_log,
          dim_log, helper);
    }
    {
      //calculate M3
      data_t *T1 = H[3];
      matrix_sub(T1, dim_log,
          B22, B_rowlen_log,
          B12, B_rowlen_log,
          dim_log);
      mult_strassen_r(M3, M3_rowlen_log,
          T1, dim_log,
          A11, A_rowlen_log,
          dim_log, helper);
    }
    {
      //calculate M4
      data_t *T1 = H[3];
      matrix_sub(T1, dim_log,
          B11, B_rowlen_log,
          B21, B_rowlen_log,
          dim_log);
      mult_strassen_r(M4, M4_rowlen_log,
          T1, dim_log,
          A22, A_rowlen_log,
          dim_log, helper);
    }
    {
      //calculate M5
      data_t *T1 = H[3];
      matrix_add(T1, dim_log,
          A12, A_rowlen_log,
          A11, A_rowlen_log,
          dim_log);
      mult_strassen_r(M5, M5_rowlen_log,
          B22, B_rowlen_log,
          T1, dim_log,
          dim_log, helper);
    }
    {
      //calculate M6
      data_t *T1 = H[3];
      data_t *T2 = H[4];
      matrix_sub(T1, dim_log,
          A11, A_rowlen_log,
          A21, A_rowlen_log,
          dim_log);
      matrix_add(T2, dim_log,
          B12, B_rowlen_log,
          B11, B_rowlen_log,
          dim_log);
      mult_strassen_r(M6, M6_rowlen_log,
          T2, dim_log,
          T1, dim_log,
          dim_log, helper);
    }
    {
      //calculate M7
      data_t *T1 = H[3];
      data_t *T2 = H[4];
      matrix_sub(T1, dim_log,
          A22, A_rowlen_log,
          A12, A_rowlen_log,
          dim_log);
      matrix_add(T2, dim_log,
          B22, B_rowlen_log,
          B21, B_rowlen_log,
          dim_log);
      mult_strassen_r(M7, M7_rowlen_log,
          T2, dim_log,
          T1, dim_log,
          dim_log, helper);
    }
    
    //print_matrix("M1",  M1,  1UL << M1_rowlen_log,  1UL << dim_log);
    //print_matrix("M2",  M2,  1UL << M2_rowlen_log,  1UL << dim_log);
    //print_matrix("M3",  M3,  1UL << M3_rowlen_log,  1UL << dim_log);
    //print_matrix("M4",  M4,  1UL << M4_rowlen_log,  1UL << dim_log);
    //print_matrix("M5",  M5,  1UL << M5_rowlen_log,  1UL << dim_log);
    //print_matrix("M6",  M6,  1UL << M6_rowlen_log,  1UL << dim_log);
    //print_matrix("M7",  M7,  1UL << M7_rowlen_log,  1UL << dim_log);
    //printf("--------------------------------------------------------------------------------\n");
    //print_matrix("C11",  C11,  1UL << C_rowlen_log,  1UL << dim_log);
    //print_matrix("C12",  C12,  1UL << C_rowlen_log,  1UL << dim_log);
    //print_matrix("C21",  C21,  1UL << C_rowlen_log,  1UL << dim_log);
    //print_matrix("C22",  C22,  1UL << C_rowlen_log,  1UL << dim_log);
    //printf("--------------------------------------------------------------------------------\n");

    {
      //calculate C11
      matrix_sub(C11, C_rowlen_log,
          M5, M5_rowlen_log,
          C11, C_rowlen_log,
          dim_log);
    //print_matrix("C11",  C11,  1UL << C_rowlen_log,  1UL << dim_log);
    //printf("--------------------------------------------------------------------------------\n");
      matrix_add(C11, C_rowlen_log,
          M4, M4_rowlen_log,
          C11, C_rowlen_log,
          dim_log);
    //print_matrix("C11",  C11,  1UL << C_rowlen_log,  1UL << dim_log);
    //printf("--------------------------------------------------------------------------------\n");
      matrix_add(C11, C_rowlen_log,
          M1, M1_rowlen_log,
          C11, C_rowlen_log,
          dim_log);
    //print_matrix("C11",  C11,  1UL << C_rowlen_log,  1UL << dim_log);
    //printf("--------------------------------------------------------------------------------\n");
    }
    {
      //calculate C12
      matrix_add(C12, C_rowlen_log,
          M3, M3_rowlen_log,
          C12, C_rowlen_log,
          dim_log);
    }
    {
      //calculate C22
      matrix_sub(C22, C_rowlen_log,
          M2, M2_rowlen_log,
          C22, C_rowlen_log,
          dim_log);
      matrix_add(C22, C_rowlen_log,
          M3, M3_rowlen_log,
          C22, C_rowlen_log,
          dim_log);
      matrix_add(C22, C_rowlen_log,
          M1, M1_rowlen_log,
          C22, C_rowlen_log,
          dim_log);
    }
    {
      //calculate C21
      matrix_add(C21, C_rowlen_log,
          M4, M4_rowlen_log,
          C21, C_rowlen_log,
          dim_log);
    }
    //print_matrix("C11",  C11,  1UL << C_rowlen_log,  1UL << dim_log);
    //print_matrix("C12",  C12,  1UL << C_rowlen_log,  1UL << dim_log);
    //print_matrix("C21",  C21,  1UL << C_rowlen_log,  1UL << dim_log);
    //print_matrix("C22",  C22,  1UL << C_rowlen_log,  1UL << dim_log);
    //printf("--------------------------------------------------------------------------------\n");
  }
  return;


}



void matrix_copy(data_t *N, uint32_t N_dim, data_t *O, uint32_t O_dim)
{
  uint32_t i, j;
  uint32_t dim;
  assert(N);
  assert(O);
  dim  = O_dim < N_dim ? O_dim : N_dim;
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      N[i*N_dim + j] = O[i*O_dim + j];
    }
  }
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
int32_t mult_strassen(data_t *C, data_t *B, data_t *A, const uint32_t N)
{
  data_t *LA = A, *LB = B, *LC = C;
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
    //first - get the dimension up to 2^n>N>2^(n-1)
    LSB_BIT_FILL(dim, dim - 1);
    dim ++;
    //need only some of those zeros, but then again - it might still be faster to zero everything.
    LA = calloc(dim * dim, sizeof(data_t));
    LB = calloc(dim * dim, sizeof(data_t));
    LC = calloc(dim * dim, sizeof(data_t));

    matrix_copy(LA, dim, A, N);
    matrix_copy(LB, dim, B, N);
    matrix_copy(LC, dim, C, N);
  }
  //now we can split the matrices in 4 up to the point they are single data_t element


  //now we create all the helper matrices we will need in future - they are of the same size, and there are LOG_2(dim) sets of them.
  LOG_BASE_2(dim_log, dim - 1);
  helper = malloc(dim_log * sizeof(data_t*));
  assert(helper);
  for(i = dim_log; i;) {
    --i;
    helper[i] = malloc(sizeof(data_t) * (1UL<<(i<<1)) * HELPER_ARRAY_CNT);
    assert(helper[i]);
  }
  //do things here.
  mult_strassen_r(LC, dim_log, LB, dim_log, LA, dim_log, dim_log, helper);

  print_matrix("LC", LC, N, N);

  //copy the result back to C
  matrix_copy(C, N, LC, dim);

  //free the memory
  for(i = dim_log; i;) {
    free(helper[--i]);
  }
  free(helper);
  helper=NULL;
  return 0;
}


