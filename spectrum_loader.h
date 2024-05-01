#ifndef SPECTRUM_LOADER_H
#define SPECTRUM_LOADER_H

#include "include/cglobals.h"
#include <vector>
#include <filesystem>
#include <optional>
#include <variant>
#include "LiteMath.h"
#include <spectrum.h>
#include <spectral/spec/spectrum.h>

#ifndef __OPENCL_VERSION__
using namespace LiteMath;
#endif

void SetEmissLut(const std::string &path);

struct Spectrum
{
  float Sample(float lambda) const;
  //std::vector<float> Resample(int channels, float lambdaOffs = 0.0f);
  std::vector<float> ResampleUniform() const;

  spec::ISpectrum::sptr spectrum;
  uint32_t id = INVALID_SPECTRUM_ID;
};

class SpectrumLoader
{
public:
  SpectrumLoader()
    : loader{}, spec_id{INVALID_SPECTRUM_ID}, spectrum{} {}
  SpectrumLoader(spec::ISpectrum::ptr &&ptr, uint32_t spec_id)
    : loader{new SimpleLoader(std::move(ptr))}, spec_id{spec_id}, spectrum{} {}
  SpectrumLoader(const std::string &str, uint32_t spec_id)
    : loader{new FileLoader(str)}, spec_id{spec_id}, spectrum{} {}
  SpectrumLoader(const spec::vec3 &v, uint32_t spec_id)
    : loader{new UpsampleLoader(v)}, spec_id{spec_id}, spectrum{} {}
  SpectrumLoader(const std::vector<float> &v, uint32_t spec_id)
    : loader{new VectorLoader(v)}, spec_id{spec_id}, spectrum{} {}

  SpectrumLoader(const SpectrumLoader &other) = delete;
  SpectrumLoader &operator=(const SpectrumLoader &other) = delete;

  SpectrumLoader(SpectrumLoader &&other)
    : loader{}, spec_id{}, spectrum{}
  {
    *this = std::move(other);
  }
  
  SpectrumLoader &operator=(SpectrumLoader &&other)
  {
    loader = std::move(other.loader);
    std::swap(spec_id, other.spec_id);
    spectrum = std::move(other.spectrum);
    return *this;
  }

  std::optional<Spectrum> &load() const;

  uint32_t id() const
  {
    return spec_id;
  }

private:
  struct ILoader
  {
    virtual void load(spec::ISpectrum::sptr &ptr) = 0;
    virtual ~ILoader() = default;
  };

  struct FileLoader : public ILoader
  {
    std::string path;
    FileLoader(const std::string &p) : path(p) {}
    void load(spec::ISpectrum::sptr &ptr) override;
  };

  struct UpsampleLoader : public ILoader
  {
    spec::vec3 rgb;
    UpsampleLoader(const spec::vec3 &p) : rgb(p) {}
    void load(spec::ISpectrum::sptr &ptr) override;
  };

  struct SimpleLoader : public ILoader
  {
    spec::ISpectrum::ptr spectrum;
    SimpleLoader(spec::ISpectrum::ptr &&ptr) : spectrum(std::move(ptr)) {}
    void load(spec::ISpectrum::sptr &ptr) override;
  };

  struct VectorLoader : public ILoader
  {
    std::vector<float> specVec;
    VectorLoader(const std::vector<float> &vec) : specVec(std::move(vec)) {}
    void load(spec::ISpectrum::sptr &ptr) override;
  };

  mutable std::unique_ptr<ILoader> loader;
  uint32_t spec_id;
  mutable std::optional<Spectrum> spectrum;
};

extern const Spectrum DUMMY_UNIFORM_SPECTRUM; 

std::optional<Spectrum> LoadSPDFromFile(const std::filesystem::path &path, uint32_t spec_id);

uint32_t UpsampleSpectrumFromColor(const float4 &color, std::vector<SpectrumLoader> &loaders, uint32_t &spec_count);
uint32_t UpsampleEmissionSpectrum(const float4 &color, float power, std::vector<SpectrumLoader> &loaders, uint32_t &spec_count);

float4 DownsampleSpectrum(const Spectrum &st);

spec::ISpectrum::ptr UpsampleRaw(const spec::vec3 &rgb);

spec::ISpectrum::ptr UpsampleEmissionRaw(const spec::vec3 &rgb, float power);

std::vector<float> ParseSpectrumStr(const std::string &specStr);

#endif