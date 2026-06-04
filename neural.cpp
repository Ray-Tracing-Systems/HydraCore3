#include "neural.h"
#include <algorithm>
#include <cmath>
#include <cassert>

namespace nn
{   

  static inline float _sigmoid(float x)
  {
    if(x >= 0) {
      return 1.0f / (1.0f + std::exp(-x));
    }
    else {
      float e = std::exp(x);
      return e / (1.0f + e);
    }
  }

  void Linear(const float *weights, const float *x, float *out, 
                      uint32_t batch_size, uint32_t in_dim, uint32_t out_dim, bool has_bias)
  {
#ifndef KERNEL_SLICER
    assert(weights != out && x != out);
#endif
    //std::fill(out, out + batch_size * out_dim, 0.0f);

    for(uint32_t i = 0; i < batch_size; ++i) {
      for(uint32_t j = 0; j < out_dim; ++j) {
        out[i * out_dim + j] = 0;
        for(uint32_t p = 0; p < in_dim; ++p) {
          out[i * out_dim + j] += weights[j * in_dim + p] * x[i * in_dim + p]; // ???
        }
        if(has_bias) {
          float bias = weights[in_dim * out_dim + j];
          out[i * out_dim + j] += bias;
        }

      }
    }
  }

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

  void Sigmoid(const float *A, float *out,
                      uint32_t m, uint32_t n)
  {
    for(uint32_t i = 0; i < m; ++i) {
      for(uint32_t j = 0; j < n; ++j) {
        out[i] = _sigmoid(A[i]);
      }
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

  void SiLU(const float *A, float *out,
                      uint32_t m, uint32_t n)
  {
    const uint32_t count = m * n; 

    for(uint32_t i = 0; i < count; ++i) {
      out[i] = A[i] * _sigmoid(A[i]);
    }
  }

}