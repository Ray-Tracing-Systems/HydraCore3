#include <fstream>

#include "spectrum.h"


float Spectrum::Sample(float lambda) const
{
  float sample;

  if (wavelengths.empty() || lambda < wavelengths.front() || lambda > wavelengths.back())
    return 0;

  
  int o = BinarySearch(wavelengths.size(), [&](int i) { return wavelengths[i] <= lambda; });

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
LambdaSample SampleWavelengths(float u, float a, float b) 
{
  // pdf is 1.0f / (b - a)
  LambdaSample res;
  const uint32_t sample_sz = sizeof(res.M) / sizeof(res.M[0]);
  
  res[0] = lerp(a, b, u);

  float delta = (b - a) / sample_sz;
  for (int i = 1; i < sample_sz; ++i) 
  {
      res[i] = res[i - 1] + delta;
      if (res[i] > b)
        res[i] = a + (res[i] - b);
  }

  return res;
}

float3 SampleCIE(const LambdaSample &lambda, const float* cie, float a = LAMBDA_MIN, float b = LAMBDA_MAX)
{
  float3 res;
  const uint32_t sample_sz = sizeof(lambda.M) / sizeof(lambda.M[0]);

  for (int i = 0; i < sample_sz; ++i) 
  {
    int offset = std::lround(lambda[i]) - a;
    if (offset < 0 || offset >= nCIESamples)
      res[i] = 0;
    else
      res[i] = cie[offset];
  }
  return res;
}

float SpectrumAverage(float3 spec) 
{
  const uint32_t sample_sz = sizeof(spec.M) / sizeof(spec.M[0]);
  float sum = spec[0];
  for (uint32_t i = 1; i < sample_sz; ++i)
    sum += spec[i];
  return sum / sample_sz;
}

float3 SpectrumToXYZ(float3 spec, const LambdaSample &lambda, float lambda_min, float lambda_max) 
{
  const uint32_t sample_sz = sizeof(lambda.M) / sizeof(lambda.M[0]);
  const float pdf = 1.0f / (lambda_max - lambda_min);

  //TODO: fix
  for (int i = 0; i < sample_sz; ++i)
    spec[i] = (pdf != 0) ? spec[i] / pdf : 0.0f;

  float3 X = SampleCIE(lambda, CIE_X, lambda_min, lambda_max);
  float3 Y = SampleCIE(lambda, CIE_Y, lambda_min, lambda_max);
  float3 Z = SampleCIE(lambda, CIE_Z, lambda_min, lambda_max);

  for (int i = 0; i < sample_sz; ++i)
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