#include "neural_loader.h"
#include <cstring>
#include <stdexcept>


namespace nn
{

  static const char HEADER_STRING[] = "hydrann2";
  static constexpr size_t TYPESTR_LEN = 8;
  static const char NBRDF_STRING[] = "nbrdf";
  static const char KANBRDF_STRING[] = "kanbrdf";

  void WeightsLoader::init()
  {

    char buf[sizeof(HEADER_STRING)];
    file.read(buf, sizeof(HEADER_STRING) - 1);
    buf[sizeof(HEADER_STRING) - 1] = '\0';
    if(strcmp(buf, HEADER_STRING) != 0 || !file.good()) {
      throw std::runtime_error("File does not contain nn weigths: " + path);
    }

    char buf1[TYPESTR_LEN + 1];
    file.read(buf1, TYPESTR_LEN);
    buf1[TYPESTR_LEN] = '\0';
    if(!file.good()) {
      throw std::runtime_error("Error occured while reading file: " + path);
    }
    if(strcmp(buf1, NBRDF_STRING) == 0) {
      type = Type::NBRDF;
    }
    else if(strcmp(buf1, KANBRDF_STRING) == 0) {
      type = Type::KANBRDF;
    }
    else {
      throw std::runtime_error("Unknown weigths format in file: " + path);
    }

    file.read(reinterpret_cast<char *>(&layers), sizeof(uint32_t));
    if(!file.good()) {
      throw std::runtime_error("Error occured while reading file: " + path);
    } 
    else {
      uint nshapes = type == Type::KANBRDF ? 4 : 2;

      const size_t pos = file.tellg();
      file.seekg(0, std::ios::end);
      const size_t end = file.tellg();
      if(end - pos < layers * nshapes * sizeof(uint32_t)) {
        throw std::runtime_error("Weigths file has incorrect size: " + path);
      }

      file.seekg(pos, std::ios::beg);

      shapes.resize(layers * nshapes);
      file.read(reinterpret_cast<char *>(shapes.data()), layers * nshapes * sizeof(uint32_t));
      if(!file.good()) {
        throw std::runtime_error("Error occured while reading file: " + path);
      }

      next();
    }
  }

  bool WeightsLoader::load_next(float *weights)
  {
    if(layers == 0) return false;

    if(layer_size > 0) {
      file.read(reinterpret_cast<char *>(weights), layer_size * sizeof(float));
    }

    if(!file) {
      layers = 0;
      layer_size = 0;
      return false;
    }

    layers -= 1;
    next();
    return true;
  }

  void WeightsLoader::next()
  {
    if(layers == 0 || next_pos >= shapes.size()) {
      layer_size = 0;
      return;
    }

    switch(type) {
    case Type::KANBRDF: 
      {
        uint32_t rows1 = shapes[next_pos + 0];
        uint32_t cols1 = shapes[next_pos + 1];
        uint32_t rows2 = shapes[next_pos + 2];
        uint32_t cols2 = shapes[next_pos + 3];
        
        layer_size = rows1 * cols1 + rows2 * (cols2 + 1);

        next_pos += 4;
      } break;
    case Type::NBRDF:
      {
        uint32_t rows = shapes[next_pos + 0];
        uint32_t cols = shapes[next_pos + 1];

        layer_size = rows * (cols + 1);

        next_pos += 2;
      } break;
    }

    const size_t pos1 = file.tellg();
    file.seekg(layer_size * sizeof(float), std::ios::cur);
    const size_t pos2 = file.tellg();

    //Check that file has required weigths
    if(pos2 - pos1 < layer_size * sizeof(float)) {
      throw std::runtime_error("Weigths file has incorrect size: " + path);
    }
    file.seekg(pos1);

  }
    
}