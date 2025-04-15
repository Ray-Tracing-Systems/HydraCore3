#ifndef NEURAL_H_
#define NEURAL_H_
#include <cinttypes>

namespace nn
{

  /**
   * A   : m x n
   * B   : n x k
   * out : m x k
   * 
   * assert(out != A && out != B)
   */
  void Matmul(const float *A, const float *B, float *out, 
                      uint32_t m, uint32_t n, uint32_t k = 1); 


  /**
   * A   : m x n
   * B   : m x n
   * out : m x n
   * 
   * (A == out || B == out) is possible
   */
  void Add(const float *A, const float *B, float *out, 
                      uint32_t m, uint32_t n = 1);

  /**
   * A   : m x n
   * B   : m x n
   * out : m x n
   * 
   * (A == out || B == out) is possible
   */
  void Sub(const float *A, const float *B, float *out, 
                      uint32_t m, uint32_t n = 1);

  /**
   * A   : m x n
   * B   : m x n
   * out : m x n
   * 
   * (A == out || B == out) is possible
   */
  void Mul(const float *A, const float *B, float *out, 
                      uint32_t m, uint32_t n = 1);

  /**
   * A   : m x n
   * B   : m x n
   * C   : m x n
   * out : m x n
   * 
   * (A == out || B == out || C == out) is possible
   */
  void FusedMulAdd(const float *A, const float *B, const float *C, float *out, 
                      uint32_t m, uint32_t n = 1);

  /**
   * A   : m x n
   * out : m x n
   * 
   * (A == out) is possible
   */
  void Neg(const float *A, float *out,
                      uint32_t m, uint32_t n = 1);
  /**
   * A   : m x n
   * out : n x m
   * 
   * assert(out != A)
   */
  void Transpose(const float *A, float *out,
                      uint32_t m, uint32_t n);





   /**
   * weights   : out_dim x in_dim + out_dim (matrix + bias)
   * x   : in_dim
   * out : out_dim
   * 
   * assert(out != A && out != B)
   */
  inline void Linear(const float *weights, const float *x, float *out, 
                      uint32_t in_dim, uint32_t out_dim)
  {
    Matmul(weights, x, out, out_dim, in_dim);
    Add(weights + in_dim * out_dim, out, out, out_dim);
  }


  /**
   * A   : m x n
   * out : m x n
   * 
   * (A == out) is possible
   */
  void ReLU(const float *A, float *out,
                      uint32_t m, uint32_t n = 1);




}



#endif