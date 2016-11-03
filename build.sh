#!/bin/sh 

gcc -g -ansi -std=c99 -O2 -msse2 -ftree-vectorize  -Isrc src/mm_strassen.c src/mm_classic.c  src/test_strassen.c -o strassen_test 
