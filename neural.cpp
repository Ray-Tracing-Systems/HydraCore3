#include "neural.h"
#include <algorithm>
#include <cmath>
#include <cassert>

namespace nn
{   
  static const char HEADER_STRING[] = "hydrann1";


  void WeightsLoader::init()
  {
    char buf[sizeof(HEADER_STRING)];
    file.read(buf, sizeof(HEADER_STRING) - 1);
    buf[sizeof(HEADER_STRING) - 1] = '\0';
    if(std::string(HEADER_STRING) != buf || !file.good()) return;

    file.read(reinterpret_cast<char *>(&layers), sizeof(uint32_t));
    if(!file.good()) {
      layers = 0;
      file.close();
    } 
    else next();
  }

  bool WeightsLoader::load_next(float *weights, float *bias)
  {
    if(layers == 0) return false;

    const uint32_t count = rows * cols;
    if(count) {
      file.read(reinterpret_cast<char *>(weights), count * sizeof(float));
      file.read(reinterpret_cast<char *>(bias), rows * sizeof(float));
    }

    if(!file) {
      layers = 0;
      rows = cols = 0;
      return false;
    }

    layers -= 1;
    next();
    return true;
  }

  void WeightsLoader::next()
  {
    if(layers == 0) {
      rows = cols = 0;
      return;
    }
    

    file.read(reinterpret_cast<char *>(&rows), sizeof(uint32_t));
    file.read(reinterpret_cast<char *>(&cols), sizeof(uint32_t));
    if(!file.good()) {
      layers = 0;
    }
  }



  void Matmul(const float *A, const float *B, float *out,  
                      uint32_t m, uint32_t n, uint32_t k)
  {
#ifndef USE_VULKAN
    assert(A != out && B != out);
#endif
    std::fill_n(out, n * k, 0.0f);
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
#ifndef USE_VULKAN
    assert(A != out);
#endif
    for(uint32_t i = 0; i < m; ++i) {
      for(uint32_t j = 0; j < n; ++j) {
        out[j * m + i] = A[i * n + j];
      }
    }
  }

}