#include "neural.h"
#include <algorithm>
#include <cmath>
#include <cassert>

namespace nn
{   

  void Matmul(const float *A, const float *B, float *out,  
                      uint32_t m, uint32_t n, uint32_t k)
  {
#ifndef KERNEL_SLICER
    assert(A != out && B != out);
#endif
    std::fill(out, out + m * k, 0.0f);
    for(uint32_t i = 0; i < m; ++i) {
      for(uint32_t j = 0; j < k; ++j) {
        for(uint32_t p = 0; p < n; ++p) {
          out[i * k + j] += A[i * n + p] * B[p * k + j]; // ???
        }
      }
    }
  }

  void Add(const float *A, const float *B, float *out, 
                      uint32_t m, uint32_t n)
  {
    const uint32_t count = m * n; 

    for(uint32_t i = 0; i < count; ++i) {
      out[i] = A[i] + B[i];
    }
  }

  void Sub(const float *A, const float *B, float *out, 
                      uint32_t m, uint32_t n)
  {
    const uint32_t count = m * n; 

    for(uint32_t i = 0; i < count; ++i) {
      out[i] = A[i] - B[i];
    }
  }

  void Neg(const float *A, float *out,
                      uint32_t m, uint32_t n)
  {
    const uint32_t count = m * n; 

    for(uint32_t i = 0; i < count; ++i) {
      out[i] = -A[i];
    }
  }

  void Mul(const float *A, const float *B, float *out,
                      uint32_t m, uint32_t n)
  {
    const uint32_t count = m * n; 

    for(uint32_t i = 0; i < count; ++i) {
      out[i] = A[i] * B[i];
    } 
  }

  void FusedMulAdd(const float *A, const float *B, const float *C, float *out, 
                      uint32_t m, uint32_t n)
  {
    const uint32_t count = m * n; 

    for(uint32_t i = 0; i < count; ++i) {
      out[i] = std::fma(A[i], B[i], C[i]);
    }
  }

  void ReLU(const float *A, float *out,
                      uint32_t m, uint32_t n)
  {
    const uint32_t count = m * n; 

    for(uint32_t i = 0; i < count; ++i) {
      out[i] = std::max(A[i], 0.0f);
    }
  }

  void Transpose(const float *A, float *out,
                      uint32_t m, uint32_t n)
  {
#ifndef KERNEL_SLICER
    assert(A != out);
#endif
    for(uint32_t i = 0; i < m; ++i) {
      for(uint32_t j = 0; j < n; ++j) {
        out[j * m + i] = A[i * n + j];
      }
    }
  }

}