#include "integrator_pt.h"
#include "fourier/fourier.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "utils.h"
#include "spectrum.h"

#include <chrono>
#include <string>
#include <vector>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

#include <iostream>

static float3 fourier_to_rgb(const FourierSpec &spec)
{
  const std::vector<float> &cie_x = Get_CIE_X();
  const std::vector<float> &cie_y = Get_CIE_Y();
  const std::vector<float> &cie_z = Get_CIE_Z();

  std::vector<float> stdspec = fourier::to_std_spectrum(spec);

  float3 xyz;

  for(size_t i = 0; i < cie_x.size(); ++i) {
    xyz.x += cie_x[i] * stdspec[i];
    xyz.y += cie_y[i] * stdspec[i];
    xyz.z += cie_z[i] * stdspec[i];
  }

  return XYZToRGB(xyz) / 106.856895f;
}


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

    float3 rgb = fourier_to_rgb(spec);

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

void Integrator::PathTraceBlockF(uint tid, float* out_spec, uint a_passNum)
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
    for (int j = 0; j < a_passNum; ++j) {
      FourierSpec tmp;
      PathTraceF(uint(i), FourierSpec::SIZE, tmp.v);

      const uint XY = m_packedXY[i];
      const uint x  = (XY & 0x0000FFFF);
      const uint y  = (XY & 0xFFFF0000) >> 16;

      for (int spec_i = 0; spec_i < FourierSpec::SIZE; ++spec_i)
        out_spec[(y*m_winWidth+x)*FourierSpec::SIZE + spec_i] += tmp[spec_i];
    }

    progress.Update();
    

  }
  //std::cout << val << std::endl;
  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}

void Integrator::FSpecToRGB(uint tid, int channels, float* out_color, float* out_spec)
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

    const uint XY = m_packedXY[i];
    const uint x  = (XY & 0x0000FFFF);
    const uint y  = (XY & 0xFFFF0000) >> 16;

    FourierSpec spec;

    for (int spec_i = 0; spec_i < FourierSpec::SIZE; ++spec_i)
      spec[spec_i] = out_spec[(y*m_winWidth+x)*FourierSpec::SIZE + spec_i];

    float3 rgb = fourier_to_rgb(spec);

    out_color[(y*m_winWidth+x)*channels + 0] = rgb.x;
    out_color[(y*m_winWidth+x)*channels + 1] = rgb.y;
    out_color[(y*m_winWidth+x)*channels + 2] = rgb.z;

    progress.Update();
    

  }
  //std::cout << val << std::endl;
  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}

void Integrator::PathTraceNBlock(uint tid, uint channels, float* out_color, uint a_passNum)
{
  ConsoleProgressBar progress(tid);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for (int i = 0; i < tid; ++i) {
    for (int j = 0; j < a_passNum; ++j) {
      PathTraceN(uint(i), channels, out_color);
    }
    progress.Update();
  }
  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}