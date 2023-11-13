#include <fstream>

#include "spectrum.h"
#include <string>

float Spectrum::Sample(float lambda) const
{
  if (wavelengths.empty() || lambda < wavelengths.front() || lambda > wavelengths.back())
    return 0;
 
  // int o = BinarySearch(wavelengths.size(), [&](int i) { return wavelengths[i] <= lambda; });
  int o = BinarySearch(wavelengths.data(), wavelengths.size(), lambda);

  float t = (lambda - wavelengths[o]) / (wavelengths[o + 1] - wavelengths[o]);
  return lerp(values[o], values[o + 1], t);
}

Spectrum LoadSPDFromFile(const std::filesystem::path &path, uint32_t spec_id)
{
  Spectrum res;

  std::ifstream in(path.string());
  std::string line;
  while(std::getline(in, line))
  {
    if(line.size() > 0 && line[0] == '#')
      continue;

    auto pos = line.find_first_of(' ');
    float lambda = std::stof(line.substr(0, pos));
    float spec   = std::stof(line.substr(pos + 1, line.size() - 1));
    res.wavelengths.push_back(lambda);
    res.values.push_back(spec);
  }
  res.id = spec_id;

  return res;
}

// "stratified" sample wavelengths in [a, b] with random number u
float4 SampleWavelengths(float u, float a, float b) 
{
  // pdf is 1.0f / (b - a)
  float4 res;

  res[0] = lerp(a, b, u);

  float delta = (b - a) / SPECTRUM_SAMPLE_SZ;
  for (uint32_t i = 1; i < SPECTRUM_SAMPLE_SZ; ++i) 
  {
      res[i] = res[i - 1] + delta;
      if (res[i] > b)
        res[i] = a + (res[i] - b);
  }

  return res;
}

float4 SampleCIE(const float4 &lambda, const float* cie, float a = LAMBDA_MIN, float b = LAMBDA_MAX)
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

float SpectrumAverage(float4 spec) 
{
  float sum = spec[0];
  for (uint32_t i = 1; i < SPECTRUM_SAMPLE_SZ; ++i)
    sum += spec[i];
  return sum / SPECTRUM_SAMPLE_SZ;
}

float3 SpectrumToXYZ(float4 spec, const float4 &lambda, float lambda_min, float lambda_max) 
{
  const float pdf = 1.0f / (lambda_max - lambda_min);

  //TODO: fix
  for (uint32_t i = 0; i < SPECTRUM_SAMPLE_SZ; ++i)
    spec[i] = (pdf != 0) ? spec[i] / pdf : 0.0f;

  float4 X = SampleCIE(lambda, CIE_X, lambda_min, lambda_max);
  float4 Y = SampleCIE(lambda, CIE_Y, lambda_min, lambda_max);
  float4 Z = SampleCIE(lambda, CIE_Z, lambda_min, lambda_max);

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