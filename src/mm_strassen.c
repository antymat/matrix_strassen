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
#include <stdlib.h>

typedef int32_t data_t;

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


#define MATRIX_PARTITION(M, M_row_len, dx, dy) ((M) + (dx) + (M_row_len) * (dy))

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
    //first - get the dimension to 2^n>N>2^(n-1)
    dim = LSB_BIT_FILL(dim - 1);
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
  dim_log = LOG_BASE_2(dim - 1);
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


int main(int argc, char *argv[]) {
  int size = 8;
  uint32_t A[] = { 1, 2, 3, 4 }, B[] = {1, 0, 0, 1}, C[4] = {0, 0, 0, 0};
  //uint32_t A[] = { 1, 2, 3, 4 }, B[] = {2, 0, 0, 2}, C[4] = {0, 0, 0, 0};
  print_matrix("A", A, 2, 2);
  print_matrix("B", B, 2, 2);
  print_matrix("C", C, 2, 2);

  mult_strassen(C, B, A, 2);
  print_matrix("C", C, 2, 2);
  printf("--------------------------------------------------------------------------------\n");
  printf("--------------------------------------------------------------------------------\n");
  printf("--------------------------------------------------------------------------------\n");

  B[0] = B[3] = 2;
  C[0] = C[1] = C[2] = C[3] = 0;
  print_matrix("A", A, 2, 2);
  print_matrix("B", B, 2, 2);
  print_matrix("C", C, 2, 2);
  mult_strassen(C, B, A, 2);
  print_matrix("C", C, 2, 2);
  {
    uint32_t A[16*16], B[16*16], C[16*16];
    uint32_t i;
    for(i=0; i<16*16; i++) { A[i] = i; B[i] = 0;}
    for(i=0; i<16; i++) B[17*i] = 2; 
    print_matrix("A", A, 16, 16);
    print_matrix("B", B, 16, 16);
    print_matrix("C", C, 16, 16);
    mult_strassen(C, B, A, 16);
    print_matrix("C", C, 16, 16);
  }
  {
    uint32_t A[16], B[16], C[16];
    uint32_t i;
    for(i=0; i<16; i++) { A[i] = i; B[i] = i;}
    print_matrix("A", A, 4, 4);
    print_matrix("B", B, 4, 4);
    print_matrix("C", C, 4, 4);
    mult_strassen(C, B, A, 4);
    print_matrix("C", C, 4, 4);
  }
  {
    uint32_t A[9], B[9], C[9];
    uint32_t i;
    for(i=0; i<9; i++) { A[i] = i; B[i] = i;}
    print_matrix("A", A, 3, 3);
    print_matrix("B", B, 3, 3);
    print_matrix("C", C, 3, 3);
    mult_strassen(C, B, A, 3);
    print_matrix("C", C, 3, 3);
  }




  printf("Ala ma asa\n");
  printf("Mem needed for array size %d = %d\n", size, get_extra_array_size(size));
  return 0;
}

