#include "neural_loader.h"
#include <cerrno>
#include <cstring>


namespace nn
{

  static const char HEADER_STRING[] = "hydrann1";


  void WeightsLoader::init()
  {

    char buf[sizeof(HEADER_STRING)];
    file.read(buf, sizeof(HEADER_STRING) - 1);
    buf[sizeof(HEADER_STRING) - 1] = '\0';
    if(std::string(HEADER_STRING) != buf || !file.good()) {
      return;
    }

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
    
}