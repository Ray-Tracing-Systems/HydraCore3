#ifndef NEURAL_LOADER_H_
#define NEURAL_LOADER_H_
#include <fstream>
#include <cinttypes>
#include <string>
#include <vector>

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
    WeightsLoader(const std::string &_path)
      : path(_path), file(_path) { init(); }

    uint32_t next_size() const { return layer_size; } //returns 0 on error
    bool load_next(float *weights);
    bool has_next() const { return layers != 0; }
  private:
    enum class Type {
      NBRDF,
      KANBRDF
    };

    std::string path;
    std::ifstream file;
    uint32_t layer_size = 0;
    uint32_t layers = 0;

    size_t next_pos = 0;
    std::vector<uint32_t> shapes;
    Type type;

    void init();
    void next();

  };
  
}

#endif