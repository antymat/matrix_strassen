/**
* @file mm_strassen.h
* @brief Strassen matrix miltiplication API.
* @author Marcin Wolcendorf
* @version 0.1
* @date 2016-11-02
*/
#ifndef __MM_STRASSEN_H__
#define __MM_STRASSEN_H__


#include "mmult_config.h"

#define HELPER_ARRAY_CNT 5u

#define CHECK_POWER_OF_2(_x) (!((_x)&((_x)-1)))

/**
* @brief Set all the bits to the right of the most significant set bit. 
*
* @param r result
* @param x parameter
*
* @return r such that all the bits left of a set MSb are 0, and right of it are 1.
*/
#define LSB_BIT_FILL(r,x) do {\
  __auto_type _y = (x);     \
  (_y) |= (_y)>>1;          \
  (_y) |= (_y)>>2;          \
  (_y) |= (_y)>>4;          \
  (_y) |= (_y)>>8;          \
  (_y) |= (_y)>>16;         \
  (r) = (_y);               \
} while(0)

/**
* @brief r = ceiling(log_2(x))
*
* @param r output
* @param x input
*
* @return r = ceiling(log_2(x))
*/
#define LOG_BASE_2(r,x) do { \
  __auto_type _y = (x); \
  _y = (_y & 0x55555555) + ((_y >> 1) & 0x55555555); \
  _y = (_y & 0x33333333) + ((_y >> 2) & 0x33333333); \
  _y = (_y & 0x0f0f0f0f) + ((_y >> 4) & 0x0f0f0f0f); \
  _y = (_y & 0x00ff00ff) + ((_y >> 8) & 0x00ff00ff); \
  _y = (_y & 0x0000ffff) + ((_y >> 16) & 0x0000ffff); \
  (r) = (_y); \
} while(0)


#define MATRIX_PARTITION(M, M_row_len, dx, dy) ((M) + (dx) + (M_row_len) * (dy))

void print_matrix(uint8_t *name, data_t *A, uint32_t row_len, uint32_t dim);
void mult_strassen(data_t *C, data_t *B, data_t *A, const uint32_t N);


#endif /* __MM_STRASSEN_H__ */

