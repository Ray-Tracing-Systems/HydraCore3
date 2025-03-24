#include "neural.h"
#include <algorithm>
#include <cmath>

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

  bool WeightsLoader::load_next(std::vector<float> &weights, std::vector<float> &bias)
  {
    if(layers == 0) return false;

    const uint32_t count = rows * cols;
    weights.resize(count);
    bias.resize(rows);
    if(count) {
      file.read(reinterpret_cast<char *>(weights.data()), count * sizeof(float));
      file.read(reinterpret_cast<char *>(bias.data()), rows * sizeof(float));
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



  void Matmul(const std::vector<float> &A, const std::vector<float> &B, std::vector<float> &out, 
                      uint32_t m, uint32_t n, uint32_t k)
  {
    out.resize(n * k);
    std::fill_n(out.begin(), n * k, 0.0f);
    for(uint32_t i = 0; i < m; ++i) {
      for(uint32_t j = 0; j < k; ++j) {
        for(uint32_t p = 0; p < n; ++p) {
          out[i * k + j] += A[i * n + p] * B[p * k + j]; // ???
        }
      }
    }
  }

  void Add(const std::vector<float> &A, const std::vector<float> &B, std::vector<float> &out, 
                      uint32_t m, uint32_t n)
  {
    if(m == 0 || n == 0) {
      m = A.size();
      n = 1;
    }
    const uint32_t count = m * n; 

    out.resize(count);
    for(uint32_t i = 0; i < count; ++i) {
      out[i] = A[i] + B[i];
    }
  }

  void Sub(const std::vector<float> &A, const std::vector<float> &B, std::vector<float> &out, 
                      uint32_t m, uint32_t n)
  {
    if(m == 0 || n == 0) {
      m = A.size();
      n = 1;
    }
    const uint32_t count = m * n; 

    out.resize(count);
    for(uint32_t i = 0; i < count; ++i) {
      out[i] = A[i] - B[i];
    }
  }

  void Neg(const std::vector<float> &A, std::vector<float> &out,
                      uint32_t m, uint32_t n)
  {
    if(m == 0 || n == 0) {
      m = A.size();
      n = 1;
    }
    const uint32_t count = m * n; 

    out.resize(count);
    for(uint32_t i = 0; i < count; ++i) {
      out[i] = -A[i];
    }
  }

  void Mul(const std::vector<float> &A, const std::vector<float> &B, std::vector<float> &out, 
                      uint32_t m, uint32_t n)
  {
    if(m == 0 || n == 0) {
      m = A.size();
      n = 1;
    }
    const uint32_t count = m * n; 

    out.resize(count);
    for(uint32_t i = 0; i < count; ++i) {
      out[i] = A[i] * B[i];
    } 
  }

  void FusedMulAdd(const std::vector<float> &A, const std::vector<float> &B, const std::vector<float> &C, std::vector<float> &out, 
                      uint32_t m, uint32_t n)
  {
    if(m == 0 || n == 0) {
      m = A.size();
      n = 1;
    }
    const uint32_t count = m * n; 

    out.resize(count);
    for(uint32_t i = 0; i < count; ++i) {
      out[i] = std::fma(A[i], B[i], C[i]);
    } 
  }


}