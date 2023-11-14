#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <vector>
#include <filesystem>
#include "LiteMath.h"
#include "include/cglobals.h"
#ifndef __OPENCL_VERSION__
using namespace LiteMath;
#endif

static constexpr float CIE_Y_integral = 106.856895f;
static constexpr uint32_t nCIESamples = 471;

struct Spectrum
{
  float Sample(float lambda) const;

  // sorted by wavelength
  std::vector<float> wavelengths; 
  std::vector<float> values;
  uint32_t id = 0;
};

inline size_t BinarySearch(const float* array, size_t array_sz, float val) 
{
  int32_t last = (int32_t)array_sz - 2, first = 1;
  while (last > 0) 
  {
    size_t half = (size_t)last >> 1, 
    middle = first + half;
    bool predResult = array[middle] <= val;
    first = predResult ? int32_t(middle + 1) : first;
    last = predResult ? last - int32_t(half + 1) : int32_t(half);
  }
  return (size_t)clamp(int32_t(first - 1), 0, int32_t(array_sz - 2));
}

float4 SampleWavelengths(float u, float a = LAMBDA_MIN, float b = LAMBDA_MAX);

static inline float4 SampleSpectrum(const float* a_spec_wavelengths, const float* a_spec_values, float4 a_wavelengths, uint32_t a_sz)
{
  float4 sample {};
  const uint spectralSamples = uint(sizeof(a_wavelengths.M) / sizeof(a_wavelengths.M[0])); 
  for(uint i = 0; i < spectralSamples; ++i)
  {
    if (a_sz == 0 || a_wavelengths[i] < a_spec_wavelengths[0] || a_wavelengths[i] > a_spec_wavelengths[a_sz - 1])
    {
      sample[i] = 0.0f;
    }
    else
    {
      int32_t last = (int32_t)a_sz - 2, first = 1;
      while (last > 0) 
      {
        size_t half = (size_t)last >> 1, 
        middle = first + half;
        bool predResult = a_spec_wavelengths[middle] <= a_wavelengths[i];
        first = predResult ? int32_t(middle + 1) : first;
        last = predResult ? last - int32_t(half + 1) : int32_t(half);
      }
      size_t o = (size_t)clamp(int32_t(first - 1), 0, int32_t(a_sz - 2));

      float t = (a_wavelengths[i] - a_spec_wavelengths[o]) / (a_spec_wavelengths[o + 1] - a_spec_wavelengths[o]);
      sample[i] =  lerp(a_spec_values[o], a_spec_values[o + 1], t);
    } 
  }
  return sample;
}

static inline float4 SampleCIE(const float4 &lambda, const float* cie, float a = LAMBDA_MIN, float b = LAMBDA_MAX)
{
  float4 res;

  for (uint32_t i = 0; i < SPECTRUM_SAMPLE_SZ; ++i) 
  {
    uint32_t offset = uint32_t(float(std::lround(lambda[i])) - a);
    if (offset < 0 || offset >= nCIESamples)
      res[i] = 0;
    else
      res[i] = cie[offset];
  }
  return res;
}

static inline float SpectrumAverage(float4 spec) 
{
  float sum = spec[0];
  for (uint32_t i = 1; i < SPECTRUM_SAMPLE_SZ; ++i)
    sum += spec[i];
  return sum / SPECTRUM_SAMPLE_SZ;
}

static inline float3 SpectrumToXYZ(float4 spec, const float4 &lambda, float lambda_min, float lambda_max,
                                   const float* a_CIE_X, const float* a_CIE_Y, const float* a_CIE_Z) 
{
  const float pdf = 1.0f / (lambda_max - lambda_min);

  //TODO: fix
  for (uint32_t i = 0; i < SPECTRUM_SAMPLE_SZ; ++i)
    spec[i] = (pdf != 0) ? spec[i] / pdf : 0.0f;

  float4 X = SampleCIE(lambda, a_CIE_X, lambda_min, lambda_max);
  float4 Y = SampleCIE(lambda, a_CIE_Y, lambda_min, lambda_max);
  float4 Z = SampleCIE(lambda, a_CIE_Z, lambda_min, lambda_max);

  for (uint32_t i = 0; i < SPECTRUM_SAMPLE_SZ; ++i)
  {
    X[i] *= spec[i];
    Y[i] *= spec[i];
    Z[i] *= spec[i];
  }

  float x = SpectrumAverage(X) / CIE_Y_integral;
  float y = SpectrumAverage(Y) / CIE_Y_integral;
  float z = SpectrumAverage(Z) / CIE_Y_integral;

  return float3{x ,y, z};
}

inline LiteMath::float3 XYZToRGB(const LiteMath::float3 xyz)
{
  LiteMath::float3 rgb;
  rgb[0] = +3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
  rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
  rgb[2] = +0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];

  return rgb;
}

Spectrum LoadSPDFromFile(const std::filesystem::path &path, uint32_t spec_id);
std::vector<float> Get_CIE_lambda();
std::vector<float> Get_CIE_X();
std::vector<float> Get_CIE_Y();
std::vector<float> Get_CIE_Z();

#endif