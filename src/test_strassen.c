/**
* @file test_strassen.c
* @brief Tests for the Strassen algorithm. 
* @author Marcin Wolcendorf
* @version 0.1
* @date 2016-11-02
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "mm_classic.h"
#include "mm_strassen.h"

#define VERBOSE_PRINT_OUTPUT_MATRIX 0x01
#define VERBOSE_PRINT_TIME 0x02
#define VERBOSE_PRINT_COMMENT 0x04
#define VERBOSE_PRINT_INPUT_MATRICES 0x08

/**
* @brief Testing multiplicatoin functions. 
*
* @param name Name of the test
* @param func function to test
* @param C matrix C (C = A x B)
* @param B matrix B
* @param A matrix A
* @param dim matrices' dimension
* @param v verbosity level
*/
void mult_test(char* name, void (*func)(data_t *, data_t *, data_t *, const uint32_t), data_t *C, data_t *B, data_t *A, uint32_t dim, uint32_t v) {
  clock_t before, after;
  if(v&VERBOSE_PRINT_INPUT_MATRICES) {
    print_matrix("A", A, dim, dim);
    print_matrix("B", B, dim, dim);
  }
  if(v&VERBOSE_PRINT_COMMENT) 
    printf("--------------%s----------------------------------------------------------\n", name);
  before = clock();
  func(C, B, A, dim);
  after = clock();
  if(v&VERBOSE_PRINT_TIME)
    printf("C=A*B  (%d x %d) in  %7.2f secs\n", dim, dim, (float)(after - before)/ CLOCKS_PER_SEC);
  if(v&VERBOSE_PRINT_OUTPUT_MATRIX)
    print_matrix("C", C, dim, dim);
  if(v&VERBOSE_PRINT_COMMENT) 
    printf("--------------------------------------------------------------------------------\n");
}




int main(int argc, char *argv[]) 
{


  {
    data_t A[] = { 1, 2, 3, 4 }, B[] = {1, 0, 0, 1}, C[4] = {0, 0, 0, 0};
    //uint32_t A[] = { 1, 2, 3, 4 }, B[] = {2, 0, 0, 2}, C[4] = {0, 0, 0, 0};
    mult_test("Strassen", &mult_strassen, C, B, A, 2, 0x0f);
    mult_test("Classic", &mult_classic, C, B, A, 2, 0x07);
  }
  {
    data_t A[16*16], B[16*16], C[16*16];
    uint32_t i;
    for(i=0; i<16*16; i++) { A[i] = i; B[i] = 0;}
    for(i=0; i<16; i++) B[17*i] = 2; 
    mult_test("Strassen",  &mult_strassen,  C, B, A, 16, 0x0f);
    mult_test("Classic",   &mult_classic,   C, B, A, 16, 0x07);
  }
  {
    data_t A[16], B[16], C[16];
    uint32_t i;
    for(i=0; i<16; i++) { A[i] = i; B[i] = i;}
    mult_test("Strassen",  &mult_strassen,  C,  B,  A,  4,  0x06);
    mult_test("Classic",   &mult_classic,   C,  B,  A,  4,  0x06);
  }
  {
    data_t A[9], B[9], C[9];
    uint32_t i;
    for(i=0; i<9; i++) { A[i] = i; B[i] = i;}
    mult_test("Strassen",  &mult_strassen,  C,  B,  A,  3,  0x0f);
    mult_test("Classic",   &mult_classic,   C,  B,  A,  3,  0x07);
  }
  {
    data_t A[500*500], B[500*500], C[500*500];
    uint32_t i;
    for(i=0; i<500*500; i++) { A[i] = i; B[i] = 0;}
    for(i=0; i<500; i++) {B[501*i] = 1;}
    mult_test("Strassen",  &mult_strassen,  C,  B,  A,  500,  0x06);
    mult_test("Classic",   &mult_classic,   C,  B,  A,  500,  0x06);
  }
  {
    data_t *A, *B, *C;
    uint32_t i;
    A = malloc(sizeof(data_t)*1000*1000);
    B = malloc(sizeof(data_t)*1000*1000);
    C = malloc(sizeof(data_t)*1000*1000);
    assert(A);
    assert(B);
    assert(C);
    for(i=0; i<1000*1000; i++) { A[i] = i; B[i] = 0;}
    for(i=0; i<1000; i++) {B[1001*i] = 1;}
    mult_test("Strassen",  &mult_strassen,  C,  B,  A,  1000,  0x06);
    mult_test("Classic",   &mult_classic,   C,  B,  A,  1000,  0x06);
    free(A); free(B); free(C);
  }
  {
    data_t *A, *B, *C;
    uint32_t i;
    A = malloc(sizeof(data_t)*2049*2049);
    B = malloc(sizeof(data_t)*2049*2049);
    C = malloc(sizeof(data_t)*2049*2049);
    assert(A);
    assert(B);
    assert(C);
    for(i=0; i<2049*2049; i++) { A[i] = i; B[i] = 0;}
    for(i=0; i<2049; i++) {B[2050*i] = 1;}
    mult_test("Strassen",  &mult_strassen,  C,  B,  A,  2049,  0x06);
    mult_test("Classic",   &mult_classic,   C,  B,  A,  2049,  0x06);
    free(A); free(B); free(C);
  }
  {
    data_t *A, *B, *C;
    uint32_t i;
    A = malloc(sizeof(data_t)*4000*4000);
    B = malloc(sizeof(data_t)*4000*4000);
    C = malloc(sizeof(data_t)*4000*4000);
    assert(A);
    assert(B);
    assert(C);
    for(i=0; i<4000*4000; i++) { A[i] = i; B[i] = 0;}
    for(i=0; i<4000; i++) {B[4001*i] = 1;}
    mult_test("Strassen",  &mult_strassen,  C,  B,  A,  4000,  0x06);
    mult_test("Classic",   &mult_classic,   C,  B,  A,  4000,  0x06);
    free(A); free(B); free(C);
  }
  return 0;
}

