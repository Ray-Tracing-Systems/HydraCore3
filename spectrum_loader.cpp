#include "spectrum_loader.h"
#include <spectral/spec/spectral_util.h>
#include <spectral/spec/conversions.h>
#include <spectral/upsample/sigpoly.h>
#include <unordered_map>
#include <fstream>
#include <iostream>

namespace {
  spec::SigPolyUpsampler upsampler{};
}

std::vector<float> Spectrum::ResampleUniform() const
{
  std::vector<float> res(size_t(LAMBDA_MAX - LAMBDA_MIN));
  for(unsigned c = 0; c < res.size(); c++)
    res[c] = spectrum->get_or_interpolate(LAMBDA_MIN + spec::Float(c));
  return res;
}

float Spectrum::Sample(float lambda) const
{
  return spectrum->get_or_interpolate(lambda);
}

spec::ISpectrum::ptr UpsampleRaw(const spec::vec3 &rgb)
{
  return upsampler.upsample_pixel(spec::Pixel::from_vec3(rgb));
}

spec::ISpectrum::ptr UpsampleAndResample(const spec::vec3 &rgb, float multiplier)
{
  spec::ISpectrum::ptr upsampled = UpsampleRaw(rgb);
  size_t count = size_t(LAMBDA_MAX - LAMBDA_MIN);
  spec::BasicSpectrum *spec = new spec::BasicSpectrum();
  for(unsigned c = 0; c < count; c++) {
    const spec::Float lambda = (LAMBDA_MIN + spec::Float(c));
    spec->set(lambda, upsampled->get_or_interpolate(lambda) * multiplier);
  }
  return spec::ISpectrum::ptr(spec);
}

std::optional<Spectrum> LoadSPDFromFile(const std::filesystem::path &path, uint32_t spec_id)
{
  spec::ISpectrum::ptr sp;
  spec::ISpectrum::csptr illum;
  std::cerr << path.string() << std::endl;

  if(spec::util::load_spectrum(path.string(), sp, illum)) {
    Spectrum res {
      .spectrum = std::move(sp),
      .id = spec_id
    };
    return res;
  }
  else {
    return {};
  }
}

void SpectrumLoader::FileLoader::load(spec::ISpectrum::sptr &ptr) 
{
  spec::ISpectrum::ptr sp;
  spec::ISpectrum::csptr illum;

  spec::util::load_spectrum(path, sp, illum);
  ptr = std::move(sp);
}

void SpectrumLoader::UpsampleLoader::load(spec::ISpectrum::sptr &ptr) 
{
  ptr = UpsampleRaw(rgb);
}

void SpectrumLoader::SimpleLoader::load(spec::ISpectrum::sptr &ptr) 
{
  ptr = std::move(spectrum);
}

std::vector<float> ParseSpectrumStr(const std::string &specStr)
{
  std::vector<float> res;

  std::stringstream ss(specStr);
  float val = 0.0f;
  while (ss >> val) 
  {
    res.push_back(val);
  }

  if(res.size() % 2 != 0)
  {
    std::cerr << "[parseSpectrumStr] : size of spectrum string is not even: " << specStr << std::endl;
  }
  
  return res;
}

void SpectrumLoader::VectorLoader::load(spec::ISpectrum::sptr &ptr)
{
  spec::BasicSpectrum spec;

  if(specVec.size() % 2 != 0)
  {
    std::cerr << "[SpectrumLoader::VectorLoader] : size of spectrum vector is not even: " << specVec.size() << std::endl;
  }

  for(size_t i = 0; i < specVec.size(); i += 2) 
  {
    spec.set(specVec[i], specVec[i + 1]);
  }
  ptr.reset(new spec::BasicSpectrum(spec));
}

std::optional<Spectrum> &SpectrumLoader::load() const
{
  if(loader) {
    Spectrum sp;
    loader->load(sp.spectrum);
    if(sp.spectrum) {
      sp.id = spec_id;
      spectrum = std::move(sp);
    }
    delete loader;
    loader = nullptr;
  }
  return spectrum;
}

namespace {

  std::unordered_map<spec::vec3, uint32_t> spec_cache;

}


uint32_t UpsampleSpectrumFromColor(const float4 &color, std::vector<SpectrumLoader> &loaders, uint32_t &spec_count)
{
  spec::vec3 rgb{color[0], color[1], color[2]};
  auto it = spec_cache.find(rgb);
  if(it != spec_cache.end()) {
    std::cerr << "Reusing spectrum for " << color.x << " " << color.y << " " << color.z << std::endl;
    return it->second;
  }
  else {
    float multiplier = LiteMath::hmax(color);
    uint32_t spec_id = spec_count++;
    std::cerr << "Creating new spectrum of " << color.x << " " << color.y << " " << color.z << std::endl;
    std::cerr << "ID = " << spec_id << std::endl;
    if(multiplier <= 1.0f) {
      loaders.push_back({rgb, spec_id});
      spec_cache[rgb] = spec_id;
    }
    else {
      auto spec = UpsampleAndResample(rgb / multiplier, multiplier);
      loaders.push_back({std::move(spec), spec_id});
      spec_cache[rgb] = spec_id;
    }
    return spec_id;
  }
}

float4 DownsampleSpectrum(const Spectrum &st)
{
  spec::vec3 rgb = spec::xyz2rgb(spec::spectre2xyz(*(st.spectrum)));
  return {rgb.x, rgb.y, rgb.z, 0.0f};
}

