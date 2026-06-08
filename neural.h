#ifndef NEURAL_H_
#define NEURAL_H_
#include <cstdint>
#include <cmath>

static inline float _sigmoid(float x)
{
  if(x >= 0) {
    return 1.0f / (1.0f + exp(-x));
  }
  else {
    float e = exp(x);
    return e / (1.0f + e);
  }
}

/**
 * A   : m x n
 * B   : n x k
 * res : m x k
 * 
 * assert(res != A && res != B)
 */
static inline void NeuralMatmul(const float *A, const float *B, float *res,  
                    uint32_t m, uint32_t n, uint32_t k)
{
  for(uint32_t i = 0; i < m; ++i) {
    for(uint32_t j = 0; j < k; ++j) {
      res[i * k + j] = 0;
      for(uint32_t p = 0; p < n; ++p) {
        res[i * k + j] += A[i * n + p] * B[p * k + j]; // ???
      }
    }
  }
}


/**
 * A   : m x n
 * B   : m x n
 * res : m x n
 * 
 * (A == res || B == res) is possible
 */
static inline void NeuralAdd(const float *A, const float *B, float *res, 
                    uint32_t m, uint32_t n)
{
  const uint32_t count = m * n; 

  for(uint32_t i = 0; i < count; ++i) {
    res[i] = A[i] + B[i];
  }
}


/**
 * A   : m x n
 * B   : m x n
 * 
 */
static inline void NeuralAdd(float *A, const float *B, uint32_t m, uint32_t n)
{
  const uint32_t count = m * n; 

  for(uint32_t i = 0; i < count; ++i) {
    A[i] += B[i];
  }  
}


/**
 * A   : m x n
 * B   : m x n
 * res : m x n
 * 
 * (A == res || B == res) is possible
 */
static inline void NeuralSub(const float *A, const float *B, float *res, 
                    uint32_t m, uint32_t n)
{
  const uint32_t count = m * n; 

  for(uint32_t i = 0; i < count; ++i) {
    res[i] = A[i] - B[i];
  }
}



 /**
 * weights   : out_dim x in_dim + out_dim (matrix + bias)
 * x   : batch_size x in_dim
 * out : batch_size x out_dim
 * 
 * assert(out != x && out != weights)
 */
static inline void NeuralLinear(const float *weights, const float *x, float *res, 
                    uint32_t batch_size, uint32_t in_dim, uint32_t out_dim, bool has_bias)
{
  for(uint32_t i = 0; i < batch_size; ++i) {
    for(uint32_t j = 0; j < out_dim; ++j) {
      res[i * out_dim + j] = 0;
      for(uint32_t p = 0; p < in_dim; ++p) {
        res[i * out_dim + j] += weights[j * in_dim + p] * x[i * in_dim + p]; // ???
      }
      if(has_bias) {
        float bias = weights[in_dim * out_dim + j];
        res[i * out_dim + j] += bias;
      }

    }
  }
}

/**
 * A   : m x n
 * out : m x n
 * 
 * (A == out) is possible
 */
static inline void NeuralSigmoid(const float *A, float *res,
                    uint32_t m, uint32_t n)
{
  for(uint32_t i = 0; i < m; ++i) {
    for(uint32_t j = 0; j < n; ++j) {
      res[i] = _sigmoid(A[i]);
    }
  }
}

/**
 * A   : m x n
 * out : m x n
 * 
 * (A == out) is possible
 */
static inline void NeuralReLU(const float *A, float *res,
                    uint32_t m, uint32_t n)
{
  const uint32_t count = m * n; 

  for(uint32_t i = 0; i < count; ++i) {
    res[i] = A[i] > 0.0 ? A[i] : 0.0; //max(A[i], 0.0f);
  }
}

/**
 * A   : m x n
 * out : m x n
 * 
 * (A == out) is possible
 */
static inline void NeuralSiLU(const float *A, float *res,
                    uint32_t m, uint32_t n)
{
  const uint32_t count = m * n; 

  for(uint32_t i = 0; i < count; ++i) {
    res[i] = A[i] * _sigmoid(A[i]);
  }
}




#endif