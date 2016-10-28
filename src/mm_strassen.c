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


size_t get_extra_mem_size(uint32_t N) 
{
  size_t ret = 0;
  //normalize the dimension to be "2^n - 1" (MSBs only 0, LSBs only 1)
  N = LSB_BIT_FILL(N-1);
  //find the n from "2^n - 1"
  N = LOG_BASE_2(N) << 1;
  ret  = (1 << N) - 1;
  ret /= (4 - 1);

  return (ret << 3) - ret;
}


int32_t mult_strassen(int32_t *C, int32_t *B, int32_t *A, uint32_t N)
{
  assert(A);
  assert(B);
  assert(C);
  assert(N);

  return 0;
}


int main(int argc, char *argv[]) {
  int size = 8;
  printf("Ala ma asa\n");
  printf("Mem needed for array size %d = %d\n", size, get_extra_mem_size(size));
  return 0;
}
