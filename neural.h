#ifndef NEURAL_H_
#define NEURAL_H_
#include <vector>
#include <string>
#include <cinttypes>
#include <fstream>

namespace nn
{


  /**
   * 
   * "hydrann1"
   * N: u32
   * [
   *  {
   *    r: u32, c: u32
   *    d: (r*c) f32
   *    b: (r) f32
   *  } 
   *  ...
   * ]
   * 
   */
  class WeightsLoader
  {
  public:
    WeightsLoader(const std::string &path)
      : file(path) { init(); }

    uint32_t next_rows() const { return rows; } //returns 0 on error
    uint32_t next_cols() const { return cols; } //returns 0 on error
    bool load_next(float *weights, float *bias);
    bool has_next() const { return layers != 0; }
  private:
    std::ifstream file;
    uint32_t rows = 0;
    uint32_t cols = 0;
    uint32_t layers = 0;

    void init();
    void next();

  };

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
   * out : m x n
   * 
   * (A == out) is possible
   */
  void ReLU(const float *A, float *out,
                      uint32_t m, uint32_t n = 1);


}



#endif