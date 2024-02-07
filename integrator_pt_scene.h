#pragma once

#include "integrator_pt.h"
#include "include/crandom.h"

#include "LiteScene/cmesh4.h"
using cmesh4::SimpleMesh;

#include "mi_materials.h"
#include "spectrum.h"
#include "LiteScene/hydraxml.h"

#include <string>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <fstream>
#include <filesystem>

#include "Image2d.h"
using LiteImage::Image2D;
using LiteImage::Sampler;
using LiteImage::ICombinedImageSampler;
using namespace LiteMath;

struct TextureInfo
{
  std::wstring path;   ///< path to file with texture data
  uint32_t     width;  ///< assumed texture width
  uint32_t     height; ///< assumed texture height
  uint32_t     bpp;    ///< assumed texture bytes per pixel, we support 4 (LDR) or 16 (HDR) during loading; Note that HDR texture could be compressed to 8 bytes (half4) on GPU.
};

struct HydraSampler
{
  float4    row0       = float4(1,0,0,0);
  float4    row1       = float4(0,1,0,0);
  float     inputGamma = 2.2f;
  bool      alphaFromRGB = true;
  
  uint32_t  texId = 0;
  Sampler   sampler;

  bool operator==(const HydraSampler& a_rhs) const
  {
    const bool addrAreSame     = (sampler.addressU == a_rhs.sampler.addressU) && (sampler.addressV == a_rhs.sampler.addressV) && (sampler.addressW == a_rhs.sampler.addressW);
    const bool filtersAreSame  = (sampler.filter == a_rhs.sampler.filter);
    const bool hasBorderSam    = (sampler.addressU == Sampler::AddressMode::BORDER || sampler.addressV == Sampler::AddressMode::BORDER || sampler.addressW == Sampler::AddressMode::BORDER);
    const bool sameBorderColor = (length3f(sampler.borderColor - a_rhs.sampler.borderColor) < 1e-5f);
    const bool sameTexId       = (texId == a_rhs.texId);
    return (addrAreSame && filtersAreSame) && (!hasBorderSam || sameBorderColor) && sameTexId;
  }
};

class HydraSamplerHash 
{
public:
  size_t operator()(const HydraSampler& sam) const
  {
    const size_t addressMode1 = size_t(sam.sampler.addressU);
    const size_t addressMode2 = size_t(sam.sampler.addressV) << 4;
    const size_t addressMode3 = size_t(sam.sampler.addressW) << 8;
    const size_t filterMode   = size_t(sam.sampler.filter)   << 12;
    return addressMode1 | addressMode2 | addressMode3 | filterMode | (size_t(sam.texId) << 16);
  }
};

Sampler::AddressMode GetAddrModeFromString(const std::wstring& a_mode);
float4x4 ReadMatrixFromString(const std::string& str);
HydraSampler ReadSamplerFromColorNode(const pugi::xml_node a_colorNodes);

std::shared_ptr<ICombinedImageSampler> MakeWhiteDummy();
std::shared_ptr<ICombinedImageSampler> LoadTextureAndMakeCombined(const TextureInfo& a_texInfo, const Sampler& a_sampler, bool a_disableGamma = false);

struct SpectrumInfo
{
  std::wstring path;   
  uint32_t id;
};

Spectrum LoadSPDFromFile(const std::filesystem::path &path, uint32_t spec_id);

std::optional<Spectrum> LoadSpectrumFromNode(const pugi::xml_node& a_node, const std::vector<SpectrumInfo> &spectraInfo);
uint32_t GetSpectrumIdFromNode(const pugi::xml_node& a_node);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<HydraSampler, uint32_t> LoadTextureFromNode(const pugi::xml_node& node, const std::vector<TextureInfo> &texturesInfo,
                                                      std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                                      std::vector< std::shared_ptr<ICombinedImageSampler> > &textures);

float4 GetColorFromNode(const pugi::xml_node& a_node, bool is_spectral_mode);

Material ConvertGLTFMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                             std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                             std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                             bool is_spectral_mode);

Material ConvertOldHydraMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                                 std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                 std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                                 bool is_spectral_mode);

Material LoadRoughConductorMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                                    std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                    std::vector< std::shared_ptr<ICombinedImageSampler> > &textures);

Material LoadDiffuseMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                             std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                             std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                             bool is_spectral_mode);

Material LoadDielectricMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                                std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                                bool is_spectral_mode);


Material LoadBlendMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                           std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                           std::vector< std::shared_ptr<ICombinedImageSampler> > &textures);

float4 image2D_average(const std::shared_ptr<ICombinedImageSampler> &tex);

Material LoadPlasticMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                             std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache,
                             std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                             std::vector<float> &precomputed_transmittance,
                             bool is_spectral_mode,
                             const std::vector<float> &spectra,
                             const std::vector<uint2> &spec_offsets);

//std::string Integrator::GetFeatureName(uint32_t a_featureId);
//std::vector<uint32_t> Integrator::PreliminarySceneAnalysis(const char* a_scenePath, const char* a_sncDir, SceneInfo* pSceneInfo);
//bool                  Integrator::LoadScene(const char* a_scenePath, const char* a_sncDir);
