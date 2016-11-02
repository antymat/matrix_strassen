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

int main(int argc, char *argv[]) 
{
    clock_t before, after;


  {
    uint32_t A[] = { 1, 2, 3, 4 }, B[] = {1, 0, 0, 1}, C[4] = {0, 0, 0, 0};
    //uint32_t A[] = { 1, 2, 3, 4 }, B[] = {2, 0, 0, 2}, C[4] = {0, 0, 0, 0};
    print_matrix("A", A, 2, 2);
    print_matrix("B", B, 2, 2);
    print_matrix("C", C, 2, 2);

    printf("--------------Strassen----------------------------------------------------------\n");
    mult_strassen(C, B, A, 2);
    print_matrix("C", C, 2, 2);
    printf("--------------Classic-----------------------------------------------------------\n");
    mult_classic(C, B, A, 2);
    print_matrix("C", C, 2, 2);
    printf("--------------------------------------------------------------------------------\n");

    B[0] = B[3] = 2;
    C[0] = C[1] = C[2] = C[3] = 0;
    print_matrix("A", A, 2, 2);
    print_matrix("B", B, 2, 2);
    print_matrix("C", C, 2, 2);
    mult_strassen(C, B, A, 2);
    print_matrix("C", C, 2, 2);
    printf("--------------------------------------------------------------------------------\n");
  }
  {
    uint32_t A[16*16], B[16*16], C[16*16];
    uint32_t i;
    for(i=0; i<16*16; i++) { A[i] = i; B[i] = 0;}
    for(i=0; i<16; i++) B[17*i] = 2; 
    print_matrix("A", A, 16, 16);
    print_matrix("B", B, 16, 16);
    print_matrix("C", C, 16, 16);
    printf("--------------Strassen----------------------------------------------------------\n");
    mult_strassen(C, B, A, 16);
    print_matrix("C", C, 16, 16);
    printf("--------------Classic-----------------------------------------------------------\n");
    mult_classic(C, B, A, 16);
    print_matrix("C", C, 16, 16);
    printf("--------------------------------------------------------------------------------\n");
  }
  {
    uint32_t A[16], B[16], C[16];
    uint32_t i;
    for(i=0; i<16; i++) { A[i] = i; B[i] = i;}
    print_matrix("A", A, 4, 4);
    print_matrix("B", B, 4, 4);
    print_matrix("C", C, 4, 4);
    printf("--------------Strassen----------------------------------------------------------\n");
    mult_strassen(C, B, A, 4);
    print_matrix("C", C, 4, 4);
    printf("--------------Classic-----------------------------------------------------------\n");
    mult_classic(C, B, A, 4);
    print_matrix("C", C, 4, 4);
    printf("--------------------------------------------------------------------------------\n");
  }
  {
    uint32_t A[9], B[9], C[9];
    uint32_t i;
    for(i=0; i<9; i++) { A[i] = i; B[i] = i;}
    print_matrix("A", A, 3, 3);
    print_matrix("B", B, 3, 3);
    print_matrix("C", C, 3, 3);
    printf("--------------Strassen----------------------------------------------------------\n");
    mult_strassen(C, B, A, 3);
    print_matrix("C", C, 3, 3);
    printf("--------------Classic-----------------------------------------------------------\n");
    mult_classic(C, B, A, 3);
    print_matrix("C", C, 3, 3);
    printf("--------------------------------------------------------------------------------\n");
  }
  {
    uint32_t A[500*500], B[500*500], C[500*500];
    uint32_t i;
    for(i=0; i<500*500; i++) { A[i] = i; B[i] = 0;}
    for(i=0; i<500; i++) {B[501*i] = 1;}
    //print_matrix("A", A, 500, 500);
    //print_matrix("B", B, 500, 500);
    //print_matrix("C", C, 500, 500);
    printf("--------------Strassen----------------------------------------------------------\n");
    before = clock();
    mult_strassen(C, B, A, 500);
    after = clock();
    printf("C=A*B in  %7.2f secs\n",(float)(after - before)/ CLOCKS_PER_SEC);
    //print_matrix("C", C, 500, 500);
    printf("--------------Classic-----------------------------------------------------------\n");
    before = clock();
    mult_classic(C, B, A, 500);
    after = clock();
    printf("C=A*B in  %7.2f secs\n",(float)(after - before)/ CLOCKS_PER_SEC);
    //print_matrix("C", C, 500, 500);
    printf("--------------------------------------------------------------------------------\n");
  }
  {
    uint32_t *A, *B, *C;
    uint32_t i;
    A = malloc(sizeof(uint32_t)*1000*1000);
    B = malloc(sizeof(uint32_t)*1000*1000);
    C = malloc(sizeof(uint32_t)*1000*1000);
    assert(A);
    assert(B);
    assert(C);
    for(i=0; i<1000*1000; i++) { A[i] = i; B[i] = 0;}
    for(i=0; i<1000; i++) {B[1001*i] = 1;}
  //  print_matrix("A", A, 1000, 1000);
  //  print_matrix("B", B, 1000, 1000);
  //  print_matrix("C", C, 1000, 1000);
  
    printf("--------------Strassen----------------------------------------------------------\n");
    before = clock();
    mult_strassen(C, B, A, 1000);
    after = clock();
    printf("C=A*B in  %7.2f secs\n",(float)(after - before)/ CLOCKS_PER_SEC);
    //print_matrix("C", C, 1000, 1000);
    printf("--------------Classic-----------------------------------------------------------\n");
    before = clock();
    mult_classic(C, B, A, 1000);
    after = clock();
    printf("C=A*B in  %7.2f secs\n",(float)(after - before)/ CLOCKS_PER_SEC);
    //print_matrix("C", C, 1000, 1000);
    printf("--------------------------------------------------------------------------------\n");
    free(A); free(B); free(C);
  }
  {
    uint32_t *A, *B, *C;
    uint32_t i;
    A = malloc(sizeof(uint32_t)*4000*4000);
    B = malloc(sizeof(uint32_t)*4000*4000);
    C = malloc(sizeof(uint32_t)*4000*4000);
    assert(A);
    assert(B);
    assert(C);
    for(i=0; i<4000*4000; i++) { A[i] = i; B[i] = 0;}
    for(i=0; i<4000; i++) {B[4001*i] = 1;}
  //  print_matrix("A", A, 4000, 4000);
  //  print_matrix("B", B, 4000, 4000);
  //  print_matrix("C", C, 4000, 4000);
  
    printf("--------------Strassen----------------------------------------------------------\n");
    before = clock();
    mult_strassen(C, B, A, 4000);
    after = clock();
    printf("C=A*B in  %7.2f secs\n",(float)(after - before)/ CLOCKS_PER_SEC);
    //print_matrix("C", C, 4000, 4000);
    printf("--------------Classic-----------------------------------------------------------\n");
    before = clock();
    mult_classic(C, B, A, 4000);
    after = clock();
    printf("C=A*B in  %7.2f secs\n",(float)(after - before)/ CLOCKS_PER_SEC);
    //print_matrix("C", C, 4000, 4000);
    printf("--------------------------------------------------------------------------------\n");
    free(A); free(B); free(C);
  }
  return 0;
}

