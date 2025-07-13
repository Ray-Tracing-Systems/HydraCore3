#include "fourier/fspec.h"
#include "include/cglobals.h"
#include "integrator_pt.h"
#include "fourier/fourier.h"

#include "utils.h"
#include "spectrum.h"

#include <LiteMath.h>
#include <chrono>
#include <vector>



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

static float3 fourier_to_rgb_conv(const FourierSpec &spec)
{
  float3 xyz;
  xyz.x = (spec * fourier::CIE_X)[0];
  xyz.y = (spec * fourier::CIE_Y)[0];
  xyz.z = (spec * fourier::CIE_Z)[0];

  return XYZToRGB((LAMBDA_MAX - LAMBDA_MIN) * xyz) / 106.856895f;
}

static inline float conv0(const float a[FourierSpec::SIZE], const float b[FourierSpec::SIZE])
{
  float res = a[0] * b[0];
  for(int i = 0; i < FourierSpec::SIZE; ++i) {
    res += a[i] * b[i];
  }
  return 0.5f * res;
}

static __attribute__((optimize("O3"))) float3 fourier_to_rgb_conv0(const FourierSpec &spec)
{
  float3 xyz;
  xyz.x = conv0(spec.v, fourier::CIE_X.v);
  xyz.y = conv0(spec.v, fourier::CIE_Y.v);
  xyz.z = conv0(spec.v, fourier::CIE_Z.v);

  return XYZToRGB((LAMBDA_MAX - LAMBDA_MIN) * xyz) / 106.856895f;
}

void Integrator::PathTraceBlockF(uint tid, int channels, float* out_color, uint a_passNum, bool use_c0)
{
  std::vector<float> wavelengths;
  for(int i = int(LAMBDA_MIN); i <= int(LAMBDA_MAX); ++i) {
    wavelengths.push_back(float(i));
  }
  const std::vector<float> phases = fourier::wl_to_phases(wavelengths);
  //float val = 0.0f;
  ConsoleProgressBar progress(tid);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for (uint i = 0; i < tid; ++i) {
    FourierSpec spec;
    for (uint j = 0; j < a_passNum; ++j) {
      FourierSpec tmp;
      PathTraceF(uint(i), FourierSpec::SIZE, tmp.v);

      spec += tmp;
    }

    const float3 rgb = use_c0 ? fourier_to_rgb_conv0(spec) : fourier_to_rgb_conv(spec);

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
    wavelengths.push_back(float(i));
  }
  const std::vector<float> phases = fourier::wl_to_phases(wavelengths);
  //float val = 0.0f;
  ConsoleProgressBar progress(tid);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for (uint i = 0; i < tid; ++i) {
    for (uint j = 0; j < a_passNum; ++j) {
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

  //auto start = std::chrono::high_resolution_clock::now();
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

  }
  //std::cout << val << std::endl;
  //shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
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