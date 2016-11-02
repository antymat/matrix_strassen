/**
* @file test_strassen.c
* @brief Tests for the Strassen algorithm. 
* @author Marcin Wolcendorf
* @version 0.1
* @date 2016-11-02
*/

#include <stdio.h>
#include <stdlib.h>
#include "mm_strassen.h"

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




  //printf("Ala ma asa\n");
  //printf("Mem needed for array size %d = %d\n", size, get_extra_array_size(size));
  return 0;
}

