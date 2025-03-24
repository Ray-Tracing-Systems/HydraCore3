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
    bool load_next(std::vector<float> &weights, std::vector<float> &bias);
    bool has_next() const { return layers != 0; }
  private:
    std::ifstream file;
    uint32_t rows = 0;
    uint32_t cols = 0;
    uint32_t layers = 0;

    void init();
    void next();

  };


  void Matmul(const std::vector<float> &A, const std::vector<float> &B, std::vector<float> &out, 
                      uint32_t m, uint32_t n, uint32_t k = 1); // A: m x n, B: n x k, out: m x k

  void Add(const std::vector<float> &A, const std::vector<float> &B, std::vector<float> &out, 
                      uint32_t m = 0, uint32_t n = 0);
  void Sub(const std::vector<float> &A, const std::vector<float> &B, std::vector<float> &out, 
                      uint32_t m = 0, uint32_t n = 0);
  void Mul(const std::vector<float> &A, const std::vector<float> &B, std::vector<float> &out, 
                      uint32_t m = 0, uint32_t n = 0);
  void FusedMulAdd(const std::vector<float> &A, const std::vector<float> &B, const std::vector<float> &C, std::vector<float> &out, 
                      uint32_t m = 0, uint32_t n = 0);

  void Neg(const std::vector<float> &A, std::vector<float> &out,
                      uint32_t m = 0, uint32_t n = 0);

}



#endif