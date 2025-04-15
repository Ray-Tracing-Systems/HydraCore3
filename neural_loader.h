#ifndef NEURAL_LOADER_H_
#define NEURAL_LOADER_H_
#include <fstream>
#include <cinttypes>
#include <string>

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
  
}

#endif