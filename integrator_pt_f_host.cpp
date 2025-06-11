#include "integrator_pt.h"
#include "fourier/fourier.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "utils.h"

#include <chrono>
#include <string>
#include <vector>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

#include <spectral/spec/basic_spectrum.h>
#include <spectral/spec/conversions.h>
#include <iostream>

void Integrator::PathTraceBlockF(uint tid, int channels, float* out_color, uint a_passNum)
{
  std::vector<float> wavelengths;
  for(int i = int(LAMBDA_MIN); i <= int(LAMBDA_MAX); ++i) {
    wavelengths.push_back(i);
  }
  const std::vector<float> phases = fourier::wl_to_phases(wavelengths);
  //float val = 0.0f;
  ConsoleProgressBar progress(tid);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for (int i = 0; i < tid; ++i) {
    FourierSpec spec;
    for (int j = 0; j < a_passNum; ++j) {
      FourierSpec tmp;
      PathTraceF(uint(i), FourierSpec::SIZE, tmp.v);

      spec += tmp;
    }

    //#pragma omp critical
    //for(int i = 0; i < FourierSpec::SIZE; ++i) {
    //  std::cout << spec[i] << " "; 
    //}
    //std::cout << std::endl;

    std::vector<float> stdspec(wavelengths.size());
    fourier::to_std_spectrum(spec, stdspec.data());

    spec::BasicSpectrum spectrum;
    for(int i = 0; i < wavelengths.size(); i += 2) {
      spectrum.set(wavelengths[i], stdspec[i]);
    }

    //val += spectrum.get_or_interpolate(633);

    spec::vec3 rgb = spec::xyz2rgb_unsafe(spectre2xyz0(spectrum)) / 106.856895f;

    const uint XY = m_packedXY[i];
    const uint x  = (XY & 0x0000FFFF);
    const uint y  = (XY & 0xFFFF0000) >> 16;


    //if(i == 0) {
      //for(int i = 0; i < FourierSpec::SIZE; ++i) {
     //   std::cout << spec.v[i] << " ";
     // }
      //std::cout << rgb << std::endl;
    //}

    out_color[(y*m_winWidth+x)*channels + 0] = rgb.x;
    out_color[(y*m_winWidth+x)*channels + 1] = rgb.y;
    out_color[(y*m_winWidth+x)*channels + 2] = rgb.z;

    progress.Update();
    

  }
  //std::cout << val << std::endl;
  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}


