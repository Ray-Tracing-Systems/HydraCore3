#include "integrator_pt.h"
#include "include/crandom.h"

#include "LiteScene/cmesh4.h"
using cmesh4::SimpleMesh;

#include "mi_materials.h"
#include "spectrum.h"
#include "include/cmat_film.h"

//#define LAYOUT_STD140 // !!! PLEASE BE CAREFUL WITH THIS !!!
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

Sampler::AddressMode GetAddrModeFromString(const std::wstring& a_mode)
{
  if(a_mode == L"clamp")
    return Sampler::AddressMode::CLAMP;
  else if(a_mode == L"wrap")
    return Sampler::AddressMode::WRAP;
  else if(a_mode == L"mirror")
    return Sampler::AddressMode::MIRROR;
  else if(a_mode == L"border")
    return Sampler::AddressMode::BORDER;
  else if(a_mode == L"mirror_once")
    return Sampler::AddressMode::MIRROR_ONCE;
  else
    return Sampler::AddressMode::WRAP;
}

float4x4 ReadMatrixFromString(const std::string& str)
{
  float4x4 res;
  std::stringstream ss(str);
  for(int i = 0; i < 4; ++i)
  {
    ss >> res.m_col[i].x >> res.m_col[i].y >> res.m_col[i].z >> res.m_col[i].w;
  }

  return res;
}


HydraSampler ReadSamplerFromColorNode(const pugi::xml_node a_colorNodes)
{
  HydraSampler res;
  auto texNode = a_colorNodes.child(L"texture");
  if(texNode == nullptr)
    return res;
  
  res.texId = texNode.attribute(L"id").as_uint();
  
  if(texNode.attribute(L"addressing_mode_u") != nullptr)
  {
    std::wstring addModeU = texNode.attribute(L"addressing_mode_u").as_string();
    res.sampler.addressU  = GetAddrModeFromString(addModeU);
  } 

  if(texNode.attribute(L"addressing_mode_v") != nullptr)
  {
    std::wstring addModeV = texNode.attribute(L"addressing_mode_v").as_string();
    res.sampler.addressV  = GetAddrModeFromString(addModeV);
  }

  if(texNode.attribute(L"addressing_mode_w") == nullptr)
    res.sampler.addressW  = res.sampler.addressV;
  else
  {
    std::wstring addModeW = texNode.attribute(L"addressing_mode_w").as_string();
    res.sampler.addressW  = GetAddrModeFromString(addModeW);
  }

  res.sampler.filter = Sampler::Filter::LINEAR;
  if(texNode.attribute(L"filter") != nullptr)
  {
    std::wstring filterMode = texNode.attribute(L"filter").as_string();
    if(filterMode == L"point" || filterMode == L"nearest")
      res.sampler.filter = Sampler::Filter::NEAREST;
    else if(filterMode == L"cubic" || filterMode == L"bicubic")
      res.sampler.filter = Sampler::Filter::CUBIC;
  }

  if(texNode.attribute(L"input_gamma") != nullptr)
    res.inputGamma = texNode.attribute(L"input_gamma").as_float();

  const std::wstring inputAlphaMode = texNode.attribute(L"input_alpha").as_string();
  if(inputAlphaMode == L"alpha")
    res.alphaFromRGB = false;
  
  // read texture matrix
  //
  std::wstringstream inputStream(texNode.attribute(L"matrix").as_string()); // in HydraXML we store matrices by rows
  for(int i=0;i<4;i++)
    inputStream >> res.row0[i];
  for(int i=0;i<4;i++)
    inputStream >> res.row1[i];
  return res;
}

//bool LoadHDRImageFromFile(const wchar_t* a_fileName, int* pW, int* pH, std::vector<float>& a_data);
//bool LoadLDRImageFromFile(const wchar_t* a_fileName, int* pW, int* pH, std::vector<uint32_t>& a_data);
//bool SaveLDRImageToFile  (const wchar_t* a_fileName, int w, int h, uint32_t* data);

std::shared_ptr<ICombinedImageSampler> MakeWhiteDummy()
{
  constexpr uint32_t WHITE = 0x00FFFFFF;
  std::shared_ptr< Image2D<uint32_t> > pTexture1 = std::make_shared< Image2D<uint32_t> >(1, 1, &WHITE);
  Sampler sampler;
  sampler.filter   = Sampler::Filter::NEAREST; 
  sampler.addressU = Sampler::AddressMode::CLAMP;
  sampler.addressV = Sampler::AddressMode::CLAMP;
  return MakeCombinedTexture2D(pTexture1, sampler);
}

std::shared_ptr<ICombinedImageSampler> LoadTextureAndMakeCombined(const TextureInfo& a_texInfo, const Sampler& a_sampler, bool a_disableGamma = false)
{
  std::shared_ptr<ICombinedImageSampler> pResult = nullptr;
  int wh[2] = {0,0};
  
  #ifdef WIN32
  std::ifstream fin(a_texInfo.path.c_str(), std::ios::binary);
  #else
  std::string   fnameA(a_texInfo.path.begin(), a_texInfo.path.end());
  std::ifstream fin(fnameA.c_str(), std::ios::binary);
  if(!fin.is_open())
    std::cout << "[LoadTextureAndMakeCombined]: can't open '" << fnameA << "'" << std::endl;
  #endif

  fin.read((char*)wh, sizeof(int)*2);
  if(wh[0] == 0 || wh[1] == 0)
  {
    std::cout << "[LoadTextureAndMakeCombined]: can't read texture from file '" << fnameA.c_str() << "'; use white dummy;" << std::endl;
    float4 data[1] = {float4(1.0f, 1.0f, 1.0f, 1.0f)};
    auto pTexture = std::make_shared< Image2D<float4> >(1, 1, data);
    pTexture->setSRGB(false);
    pResult = MakeCombinedTexture2D(pTexture, a_sampler);
  }
  else if(a_texInfo.bpp == 16)
  {
    std::vector<float> data(wh[0]*wh[1]*4);
    fin.read((char*)data.data(), sizeof(float)*4*data.size());
    fin.close();

    auto pTexture = std::make_shared< Image2D<float4> >(wh[0], wh[1], (const float4*)data.data());
    pResult = MakeCombinedTexture2D(pTexture, a_sampler);
  }
  else
  {
    std::vector<uint32_t> data(wh[0]*wh[1]);
    fin.read((char*)data.data(), sizeof(uint32_t)*data.size());
    fin.close();

    //#TODO: if old-version gamma 2.2 is globally enabled, use a trick
    //#TODO: use gamma and invserse sRGB to get same results with sSRB   

    auto pTexture = std::make_shared< Image2D<uint32_t> >(wh[0], wh[1], data.data());
    pTexture->setSRGB(!a_disableGamma);
    pResult = MakeCombinedTexture2D(pTexture, a_sampler);
  }
 
  return pResult;
}

struct SpectrumInfo
{
  std::wstring path;   
  uint32_t id;
};

Spectrum LoadSPDFromFile(const std::filesystem::path &path, uint32_t spec_id);

std::optional<Spectrum> LoadSpectrumFromNode(const pugi::xml_node& a_node, const std::vector<SpectrumInfo> &spectraInfo)
{
  std::optional<Spectrum> spec;
  auto specNode = a_node.child(L"spectrum");
  if(specNode != nullptr)
  {
    uint32_t spec_id = specNode.attribute(L"id").as_uint();
    spec = LoadSPDFromFile(spectraInfo[spec_id].path, spec_id);
  }

  return spec;
}

uint32_t GetSpectrumIdFromNode(const pugi::xml_node& a_node)
{
  uint32_t spec_id = 0xFFFFFFFF;
  auto specNode = a_node.child(L"spectrum");
  if(specNode != nullptr)
  {
    spec_id = specNode.attribute(L"id").as_uint();
  }

  return spec_id;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<HydraSampler, uint32_t> LoadTextureFromNode(const pugi::xml_node& node, const std::vector<TextureInfo> &texturesInfo,
                                                      std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                                      std::vector< std::shared_ptr<ICombinedImageSampler> > &textures)
{
  HydraSampler sampler = ReadSamplerFromColorNode(node);
  auto p = texCache.find(sampler);
  uint32_t texId = 0;
  if(p == texCache.end())
  {
    texCache[sampler] = uint(textures.size());
    auto texNode = node.child(L"texture");
    texId  = texNode.attribute(L"id").as_uint();

    bool disableGamma = false;
    if(texNode.attribute(L"input_gamma").as_int() == 1)
      disableGamma = true;

    textures.push_back(LoadTextureAndMakeCombined(texturesInfo[texId], sampler.sampler, disableGamma));
    p = texCache.find(sampler);
  }

  return {sampler, p->second};
}

float4 GetColorFromNode(const pugi::xml_node& a_node, bool is_spectral_mode)
{
  auto val = hydra_xml::readvalVariant(a_node);
  if(std::holds_alternative<float>(val))
  {
    return float4(std::get<float>(val));
  }
  else if(std::holds_alternative<float3>(val))
  {
    if(is_spectral_mode == true)
    {
      std::cout << "WARNING! Reading float3 color value in spectral mode. Spectral upsampling not implemented yet, results will be incorrect." << std::endl;
    }
    return to_float4(std::get<float3>(val), 0.0f);
  }
  else
  {
    return std::get<float4>(val);
  }
}

Material ConvertOldHydraMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                                 std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                 std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                                 bool is_spectral_mode)
{
  std::wstring name            = materialNode.attribute(L"name").as_string();
  Material mat                 = {};
  mat.mtype                    = MAT_TYPE_GLTF;
  mat.data[GLTF_FLOAT_ALPHA]   = 0.0f;
  mat.colors[GLTF_COLOR_COAT]  = float4(1,1,1,1); 
  mat.colors[GLTF_COLOR_METAL] = float4(0,0,0,0);  
  
  mat.lightId                  = uint(-1);
  
  auto nodeEmiss = materialNode.child(L"emission");

  // read Hydra or GLTF materials
  //
  float4 color(0.0f, 0.0f, 0.0f, 0.0f);

  bool is_emission_color = false;
  if(materialNode.attribute(L"light_id") != nullptr || nodeEmiss != nullptr)
  {
    auto nodeEmissColor = nodeEmiss.child(L"color");
    color               = GetColorFromNode(nodeEmissColor, is_spectral_mode);
    is_emission_color = (materialNode.attribute(L"light_id") != nullptr) || (length(color) > 1e-5f );

    const auto& [emissiveSampler, texID] = LoadTextureFromNode(nodeEmissColor, texturesInfo, texCache, textures);
    
    mat.row0 [0] = emissiveSampler.row0;
    mat.row1 [0] = emissiveSampler.row1;
    mat.texid[0] = texID;
    
    mat.colors[EMISSION_COLOR] = color;
    if(materialNode.attribute(L"light_id") == nullptr)
      mat.lightId = uint(-1);
    else
      mat.lightId = uint(materialNode.attribute(L"light_id").as_int());  // for correct process of "-1"

    auto specId  = GetSpectrumIdFromNode(nodeEmissColor);  
    mat.spdid[0] = specId;
    mat.mtype    = MAT_TYPE_LIGHT_SOURCE;

    auto colorMultNode = nodeEmissColor.child(L"multiplier");
    if(colorMultNode)
    {
      mat.data[EMISSION_MULT] = hydra_xml::readval1f(nodeEmissColor.child(L"multiplier")); 
    }
    else
    {
      mat.data[EMISSION_MULT] = 1.0f;
    }
  }

  auto nodeDiffColor = materialNode.child(L"diffuse").child(L"color");
  if(nodeDiffColor != nullptr)
  {
    color = GetColorFromNode(nodeDiffColor, is_spectral_mode);
    const auto& [diffSampler, texID] = LoadTextureFromNode(nodeDiffColor, texturesInfo, texCache, textures);
    
    mat.row0 [0] = diffSampler.row0;
    mat.row1 [0] = diffSampler.row1;
    mat.texid[0] = texID;
  }

  float4 reflColor     = float4(0, 0, 0, 0);
  float reflGlossiness = 1.0f;
  float fresnelIOR     = 1.5f;
  auto nodeRefl        = materialNode.child(L"reflectivity");
  if(nodeRefl != nullptr)
  {
    reflColor       = GetColorFromNode(nodeRefl.child(L"color"), is_spectral_mode);
    reflGlossiness  = hydra_xml::readval1f(nodeRefl.child(L"glossiness"));  
    fresnelIOR      = hydra_xml::readval1f(nodeRefl.child(L"fresnel_ior"));
  }

  float4 transpColor      = float4(0, 0, 0, 0);
  float  transpGlossiness = 1.0f;
  float  transpIOR        = 1.5f;

  auto nodeTransp = materialNode.child(L"transparency");
  if (nodeTransp != nullptr)
  {
    transpColor      = GetColorFromNode(nodeTransp.child(L"color"), is_spectral_mode);
    transpGlossiness = hydra_xml::readval1f(nodeTransp.child(L"glossiness"));
    transpIOR        = hydra_xml::readval1f(nodeTransp.child(L"ior"));
  }

  const bool hasFresnel  = (nodeRefl.child(L"fresnel").attribute(L"val").as_int() != 0);
  if(!hasFresnel)
    fresnelIOR = 0.0f;
  
  if((length(reflColor) > 1e-5f && length(to_float3(color)) > 1e-5f) || hasFresnel)
  {
    mat.mtype   = MAT_TYPE_GLTF;
    mat.lightId = uint(-1);

    mat.colors[GLTF_COLOR_BASE]  = color;
    mat.colors[GLTF_COLOR_COAT]  = reflColor;

    if(hasFresnel)
    {
      mat.data[GLTF_FLOAT_ALPHA]   = 0.0f;
      mat.colors[GLTF_COLOR_COAT]  = reflColor;
      mat.colors[GLTF_COLOR_METAL] = float4(0,0,0,0); 
      mat.cflags                   = GLTF_COMPONENT_LAMBERT | GLTF_COMPONENT_COAT;
      SetMiPlastic(&mat, fresnelIOR, 1.0f, color, reflColor);
    }
    else
    {
      mat.data[GLTF_FLOAT_ALPHA]   = length(reflColor)/( length(reflColor) + length3f(color) );
      mat.colors[GLTF_COLOR_COAT]  = float4(0,0,0,0); 
      mat.colors[GLTF_COLOR_METAL] = reflColor;   // disable coating for such blend type
      mat.cflags                   = GLTF_COMPONENT_LAMBERT | GLTF_COMPONENT_METAL;
    }
  }
  else if(length(reflColor) > 1e-5f)
  {
    mat.mtype  = MAT_TYPE_GLTF;
    mat.cflags = GLTF_COMPONENT_METAL;
    mat.colors[GLTF_COLOR_METAL]   = reflColor;
    mat.colors[GLTF_COLOR_COAT]    = float4(0,0,0,0); 
    mat.data[GLTF_FLOAT_ALPHA]     = 1.0f;
  }
  else if(length(to_float3(color)) > 1e-5f)
  {
    mat.mtype  = MAT_TYPE_GLTF;
    mat.cflags = GLTF_COMPONENT_LAMBERT;
    mat.colors[GLTF_COLOR_BASE]  = color;
    mat.colors[GLTF_COLOR_COAT]  = float4(0,0,0,0); 
    mat.colors[GLTF_COLOR_METAL] = float4(0,0,0,0);    
    mat.data[GLTF_FLOAT_ALPHA]   = 0.0f;
  }
    
  // Glass
  if (length(transpColor) > 1e-5f)
  {
    mat.mtype                           = MAT_TYPE_GLASS;    
    mat.colors[GLASS_COLOR_REFLECT]     = reflColor;
    mat.colors[GLASS_COLOR_TRANSP]      = transpColor;      
    mat.data[GLASS_FLOAT_GLOSS_REFLECT] = reflGlossiness;
    mat.data[GLASS_FLOAT_GLOSS_TRANSP]  = transpGlossiness;
    mat.data[GLASS_FLOAT_IOR]           = fresnelIOR;
  }

  if(is_emission_color)
  {
    mat.mtype = MAT_TYPE_LIGHT_SOURCE;
  }

  auto nodeDiffRough = materialNode.child(L"diffuse").child(L"roughness");

  if (nodeDiffRough != nullptr)
  {
    mat.data[GLTF_FLOAT_ROUGH_ORENNAYAR] = hydra_xml::readval1f(nodeDiffRough);
    mat.cflags = mat.cflags | GLTF_COMPONENT_ORENNAYAR;
  }

  mat.data[GLTF_FLOAT_GLOSINESS] = reflGlossiness;
  mat.data[GLTF_FLOAT_IOR]       = fresnelIOR;

  return mat;
}

Material LoadRoughConductorMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                                    std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                    std::vector< std::shared_ptr<ICombinedImageSampler> > &textures)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  Material mat = {};
  mat.colors[CONDUCTOR_COLOR]  = float4(1, 1, 1, 1);
  mat.mtype                    = MAT_TYPE_CONDUCTOR;
  mat.lightId                  = uint(-1);

  // auto nodeBSDF = materialNode.child(L"bsdf");

  float alpha_u = 0.0f;
  float alpha_v = 0.0f;

  //auto bsdf_type = nodeBSDF.attribute(L"type").as_string();

  auto nodeAlpha = materialNode.child(L"alpha");
  if(nodeAlpha != nullptr)
  {
    alpha_u = nodeAlpha.attribute(L"val").as_float();
    alpha_v = alpha_u;

    const auto& [sampler, texID] = LoadTextureFromNode(nodeAlpha, texturesInfo, texCache, textures);

    if(texID != 0)
      alpha_u = alpha_v = 1.0f;
    
    mat.row0 [0] = sampler.row0;
    mat.row1 [0] = sampler.row1;
    mat.texid[0] = texID;
  }
  else
  {
    auto nodeAlphaU = materialNode.child(L"alpha_u");
    auto nodeAlphaV = materialNode.child(L"alpha_v");

    alpha_u = nodeAlphaU.attribute(L"val").as_float();
    alpha_v = nodeAlphaV.attribute(L"val").as_float();
  }
  
  auto eta       = materialNode.child(L"eta").attribute(L"val").as_float();
  auto etaSpecId = GetSpectrumIdFromNode(materialNode.child(L"eta"));
  auto k         = materialNode.child(L"k").attribute(L"val").as_float();
  auto kSpecId   = GetSpectrumIdFromNode(materialNode.child(L"k"));
  
  mat.data[CONDUCTOR_ROUGH_U] = alpha_u;
  mat.data[CONDUCTOR_ROUGH_V] = alpha_v; 
  mat.data[CONDUCTOR_ETA]     = eta; 
  mat.data[CONDUCTOR_K]       = k;   
  
  mat.spdid[0] = etaSpecId;
  mat.spdid[1] = kSpecId;

  return mat;
}

static inline void save_to_file(const char* name, float *arr, int x_samples, int y_samples)
{
  std::ofstream precomp_file;
  precomp_file.open(name);
  for (int i = 0; i < x_samples; ++i)
  {
    for (int j = 0; j < y_samples; ++j)
    {
      precomp_file << arr[i * y_samples + j] << " ";
    }
  }
  precomp_file.close();
}

ThinFilmPrecomputed precomputeThinFilm(const float extIOR, const uint* eta_id, const uint* k_id, const std::vector<float> &spec_values, 
        const std::vector<float> &wavelengths, const std::vector<uint2> &spec_offsets, const float* a_thickness, int layers)
{
  ThinFilmPrecomputed res;
  for (int i = 0; i < FILM_LENGTH_RES; ++i)
  {
    float wavelength = (LAMBDA_MAX - LAMBDA_MIN) / (FILM_LENGTH_RES - 1) * i + LAMBDA_MIN;
    std::vector<float> eta, k;
    eta.reserve(layers);
    k.reserve(layers);
    uint2 data;
    uint offset;
    uint size;
    for (int layer = 0; layer < layers; ++layer)
    {
      data  = spec_offsets[eta_id[layer]];
      offset = data.x;
      size   = data.y;
      eta[layer] = SampleSpectrum(wavelengths.data() + offset, spec_values.data() + offset, {wavelength, 0, 0, 0}, size)[0];

      data  = spec_offsets[k_id[layer]];
      offset = data.x;
      size   = data.y;
      k[layer] = SampleSpectrum(wavelengths.data() + offset, spec_values.data() + offset, {wavelength, 0, 0, 0}, size)[0];
    }
    for (int j = 0; j < FILM_ANGLE_RES; ++j)
    {
      float theta = M_PI_2 / float(FILM_ANGLE_RES - 1) * j;
      float ext_eta = 1.f;
      FrReflRefr forward;
      FrReflRefr backward;
      if (ext_eta * sinf(theta) > eta[layers - 1])
      {
        forward = {1.f, 0.f};
      }
      else
      {
        forward = multFrFilmReflRefr(extIOR, cosf(theta), eta.data(), k.data(), a_thickness, layers, wavelength);
      }
      if (eta[layers - 1] * sinf(theta) > ext_eta)
      {
        backward = {1.f, 0.f};
      }
      else
      {
        backward = multFrFilmReflRefr_r(extIOR, cosf(theta), eta.data(), k.data(), a_thickness, layers, wavelength);
      }
 
      res.ext_reflectivity[i * FILM_ANGLE_RES + j] = forward.refl;
      res.ext_transmittivity[i * FILM_ANGLE_RES + j] = forward.refr;
      res.int_reflectivity[i * FILM_ANGLE_RES + j] = backward.refl;
      res.int_transmittivity[i * FILM_ANGLE_RES + j] = backward.refr;  
    }
  }
  //save_to_file("../precomputed_film_refl_ext.txt", res.ext_reflectivity.data(), FILM_LENGTH_RES, FILM_ANGLE_RES);
  //save_to_file("../precomputed_film_refl_int.txt", res.int_reflectivity.data(), FILM_LENGTH_RES, FILM_ANGLE_RES);
  //save_to_file("../precomputed_film_refr_ext.txt", res.ext_transmittivity.data(), FILM_LENGTH_RES, FILM_ANGLE_RES);
  //save_to_file("../precomputed_film_refr_int.txt", res.int_transmittivity.data(), FILM_LENGTH_RES, FILM_ANGLE_RES);
  return res;
}

Material LoadThinFilmMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                              std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                              std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                              std::vector<float> &precomputed_film, std::vector<float> &thickness_vec,
                              std::vector<uint> &spec_id_vec, std::vector<float> &eta_k_vec,
                              const std::vector<float> &spec_values,
                              const std::vector<float> &wavelengths,
                              const std::vector<uint2> &spec_offsets)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  Material mat = {};
  mat.colors[FILM_COLOR]  = float4(1, 1, 1, 0);
  mat.mtype                    = MAT_TYPE_THIN_FILM;
  mat.lightId                  = uint(-1);

  auto nodeBSDF = materialNode.child(L"bsdf");

  float alpha_u = 0.0f;
  float alpha_v = 0.0f;

  auto nodeAlpha = materialNode.child(L"alpha");
  if(nodeAlpha != nullptr)
  {
    alpha_u = nodeAlpha.attribute(L"val").as_float();
    alpha_v = alpha_u;

    const auto& [sampler, texID] = LoadTextureFromNode(nodeAlpha, texturesInfo, texCache, textures);

    if(texID != 0)
      alpha_u = alpha_v = 1.0f;
    
    mat.row0 [0]  = sampler.row0;
    mat.row1 [0]  = sampler.row1;
    mat.data[FILM_TEXID0] = as_float(texID);
  }
  else
  {
    auto nodeAlphaU = materialNode.child(L"alpha_u");
    auto nodeAlphaV = materialNode.child(L"alpha_v");

    alpha_u = nodeAlphaU.attribute(L"val").as_float();
    alpha_v = nodeAlphaV.attribute(L"val").as_float();
  }

  mat.data[FILM_ETA_EXT] = 1.00028f; // air

  auto nodeExtIOR = materialNode.child(L"ext_ior");
  if(nodeExtIOR != nullptr)
  {
    mat.data[FILM_ETA_EXT] = nodeExtIOR.attribute(L"val").as_float();
  }

  mat.data[FILM_THICKNESS_OFFSET] = as_float((uint) thickness_vec.size());
  mat.data[FILM_ETA_SPECID_OFFSET] = as_float((uint) spec_id_vec.size());
  mat.data[FILM_ETA_OFFSET] = as_float((uint) eta_k_vec.size());

  uint layers = 0;
  for (auto layerNode : materialNode.child(L"layers").children())
  {
    layers++;
    if (layerNode.child(L"thickness") != nullptr)
    {
      thickness_vec.push_back(layerNode.child(L"thickness").attribute(L"val").as_float());
    }
    eta_k_vec.push_back(materialNode.child(L"eta").attribute(L"val").as_float());
    spec_id_vec.push_back(GetSpectrumIdFromNode(layerNode.child(L"eta")));
  }
  mat.data[FILM_LAYERS_COUNT] = as_float(layers);
  mat.data[FILM_K_SPECID_OFFSET] = as_float((uint) spec_id_vec.size());
  mat.data[FILM_K_OFFSET] = as_float((uint) eta_k_vec.size());

  for (auto layerNode : materialNode.child(L"layers").children())
  {
    eta_k_vec.push_back(materialNode.child(L"k").attribute(L"val").as_float());
    spec_id_vec.push_back(GetSpectrumIdFromNode(layerNode.child(L"k")));
  }

  uint precompFlag = 0;
  auto nodePrecomp = materialNode.child(L"precompute");
  if(nodePrecomp != nullptr && nodePrecomp.attribute(L"val").as_uint() > 0)
  {
    precompFlag = 1;
    mat.data[FILM_PRECOMP_FLAG] = as_float(precompFlag);
    mat.data[FILM_PRECOMP_ID] = as_float(precomputed_film.size() / (sizeof(ThinFilmPrecomputed) / sizeof(float)));
    auto precomputed = precomputeThinFilm(mat.data[FILM_ETA_EXT], spec_id_vec.data(), spec_id_vec.data() + layers, spec_values, wavelengths, 
          spec_offsets, thickness_vec.data(), layers);
    std::copy(precomputed.ext_reflectivity.begin(), precomputed.ext_reflectivity.end(), std::back_inserter(precomputed_film));
    std::copy(precomputed.ext_transmittivity.begin(), precomputed.ext_transmittivity.end(), std::back_inserter(precomputed_film));
    std::copy(precomputed.int_reflectivity.begin(), precomputed.int_reflectivity.end(), std::back_inserter(precomputed_film));
    std::copy(precomputed.int_transmittivity.begin(), precomputed.int_transmittivity.end(), std::back_inserter(precomputed_film));
  }
  else
  {
    mat.data[FILM_PRECOMP_ID] = 0.f;
  }
  
  mat.data[FILM_ROUGH_U] = alpha_u;
  mat.data[FILM_ROUGH_V] = alpha_v;

  return mat;
}

Material LoadDiffuseMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                             std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                             std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                             bool is_spectral_mode)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  Material mat = {};
  mat.colors[DIFFUSE_COLOR]   = float4(1, 1, 1, 1);

  mat.mtype    = MAT_TYPE_DIFFUSE;
  mat.lightId  = uint(-1);
  mat.texid[0] = 0;
  mat.spdid[0] = uint(-1);
  mat.data[DIFFUSE_ROUGHNESS] = 0.0f;

  static const std::wstring orenNayarNameStr {L"oren-nayar"};

  auto bsdfType = materialNode.child(L"bsdf").attribute(L"type").as_string();
  if(bsdfType == orenNayarNameStr)
  {
    mat.cflags = GLTF_COMPONENT_ORENNAYAR;
  
    auto nodeRoughness = materialNode.child(L"roughness");
    if(nodeRoughness != nullptr)
    {
      auto roughness = hydra_xml::readval1f(nodeRoughness);
      mat.data[DIFFUSE_ROUGHNESS] = roughness;
    }
  }

  auto nodeColor = materialNode.child(L"reflectance");
  if(nodeColor != nullptr)
  {
    mat.colors[DIFFUSE_COLOR] = GetColorFromNode(nodeColor, is_spectral_mode);

    const auto& [sampler, texID] = LoadTextureFromNode(nodeColor, texturesInfo, texCache, textures);
    
    mat.row0 [0]  = sampler.row0;
    mat.row1 [0]  = sampler.row1;
    mat.texid[0]  = texID;

    auto specId  = GetSpectrumIdFromNode(nodeColor);
    mat.spdid[0] = specId;
  }

  return mat;
}


Material LoadDielectricMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                                std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                                bool is_spectral_mode)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  Material mat = {};
  mat.colors[DIELECTRIC_COLOR_REFLECT]  = float4(1, 1, 1, 1);
  mat.colors[DIELECTRIC_COLOR_TRANSMIT] = float4(1, 1, 1, 1);
  mat.mtype                             = MAT_TYPE_DIELECTRIC;  
  mat.lightId                           = uint(-1);
  mat.data[DIELECTRIC_ETA_EXT]          = 1.00028f; // air
  mat.data[DIELECTRIC_ETA_INT]          = 1.5046f;  // bk7 glass
  mat.spdid[0]                          = uint(-1);

  auto nodeIntIOR = materialNode.child(L"int_ior");
  if(nodeIntIOR != nullptr)
  {
    auto specId = GetSpectrumIdFromNode(nodeIntIOR);
    mat.spdid[0] = specId;
    mat.data[DIELECTRIC_ETA_INT] = nodeIntIOR.attribute(L"val").as_float();
  }

  auto nodeExtIOR = materialNode.child(L"ext_ior");
  if(nodeExtIOR != nullptr)
  {
    mat.data[DIELECTRIC_ETA_EXT] = nodeExtIOR.attribute(L"val").as_float();
  }

  auto nodeReflColor = materialNode.child(L"reflectance");
  if(nodeReflColor != nullptr)
  {
    mat.colors[DIELECTRIC_COLOR_REFLECT] = GetColorFromNode(nodeReflColor, is_spectral_mode);
  }

  auto nodeTransColor = materialNode.child(L"transmittance");
  if(nodeTransColor != nullptr)
  {
    mat.colors[DIELECTRIC_COLOR_TRANSMIT] = GetColorFromNode(nodeTransColor, is_spectral_mode);
  }

  return mat;
}


Material LoadBlendMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                           std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                           std::vector< std::shared_ptr<ICombinedImageSampler> > &textures)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  Material mat = {};

  mat.mtype    = MAT_TYPE_BLEND;
  mat.cflags   = uint(-1);
  mat.texid[0] = 0;
  mat.data[BLEND_WEIGHT] = 1.0f;
  
  mat.datai[0] = materialNode.child(L"bsdf_1").attribute(L"id").as_uint();
  mat.datai[1] = materialNode.child(L"bsdf_2").attribute(L"id").as_uint();

  auto nodeWeight = materialNode.child(L"weight");
  if(nodeWeight != nullptr)
  {
    mat.data[BLEND_WEIGHT] = hydra_xml::readval1f(nodeWeight);

    const auto& [sampler, texID] = LoadTextureFromNode(nodeWeight, texturesInfo, texCache, textures);
    
    mat.row0 [0]  = sampler.row0;
    mat.row1 [0]  = sampler.row1;
    mat.texid[0]  = texID;
  }

  return mat;
}

float4 image2D_average(const std::shared_ptr<ICombinedImageSampler> &tex)
{
  float* ptr = (float*)(tex->data());
  float4 res{0.0f};
  size_t tex_sz = tex->width() * tex->height();
  uint32_t channels = tex->bpp() / sizeof(float);
  for(size_t i = 0; i < tex_sz / channels; ++i)
  {  
    for(size_t j = 0; j < channels; ++i)
    {
      res.M[j] += ptr[i * channels + j];
    }
  }

  if(channels == 1)
  {
    res.w = res.x;
    res.z = res.x;
    res.y = res.x;
  }

  res = res / (tex_sz);

  return res;
}

Material LoadPlasticMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                             std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache,
                             std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                             std::vector<float> &precomputed_transmittance,
                             bool is_spectral_mode,
                             const std::vector<float> &spectra,
                             const std::vector<float> &wavelengths,
                             const std::vector<uint2> &spec_offsets)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  Material mat = {};
  
  mat.mtype     = MAT_TYPE_PLASTIC;
  mat.lightId   = uint(-1);
  mat.nonlinear = 0;
  mat.texid[0]  = 0;
  mat.spdid[0]  = uint(-1);

  auto nodeColor = materialNode.child(L"reflectance");
  uint32_t specId = 0xFFFFFFFF;
  if(nodeColor != nullptr)
  {
    mat.colors[PLASTIC_COLOR] = GetColorFromNode(nodeColor, is_spectral_mode);

    const auto& [sampler, texID] = LoadTextureFromNode(nodeColor, texturesInfo, texCache, textures);

    mat.row0 [0]  = sampler.row0;
    mat.row1 [0]  = sampler.row1;
    mat.texid[0]  = texID;

    specId = GetSpectrumIdFromNode(nodeColor);
    mat.spdid[0] = specId;
  }

  float internal_ior = hydra_xml::readval1f(materialNode.child(L"int_ior"), 1.49f);
  float external_ior = hydra_xml::readval1f(materialNode.child(L"ext_ior"), 1.000277);

  mat.data[PLASTIC_IOR_RATIO] = internal_ior / external_ior;

  mat.data[PLASTIC_ROUGHNESS] = hydra_xml::readval1f(materialNode.child(L"alpha"), 0.1f);

  // dirty hack 
  if(mat.data[PLASTIC_ROUGHNESS] == 0.0f)
  {
    mat.data[PLASTIC_ROUGHNESS] = 1e-6f;
  }

  mat.nonlinear = hydra_xml::readval1u(materialNode.child(L"nonlinear"), 0);

  std::vector<float> spectrum;

  if(is_spectral_mode && specId != 0xFFFFFFFF)
  {
    const auto offsets = spec_offsets[specId];
    spectrum.reserve(offsets.y);
    for(size_t i = offsets.x; i < offsets.x + offsets.y; ++i)
    {
      if(wavelengths[i] >= LAMBDA_MIN && wavelengths[i] <= LAMBDA_MAX)
      {
        spectrum.push_back(spectra[i]);
      }
    }

    // std::copy(spectra.begin() + offsets.x, spectra.begin() + offsets.x + offsets.y, std::back_inserter(spectrum));
  }

  float4 diffuse_reflectance = mat.colors[PLASTIC_COLOR];

  if(!is_spectral_mode)
  {
    uint32_t colorTexId = mat.texid[0];
    if(colorTexId > 0 && colorTexId != 0xFFFFFFFF)
      diffuse_reflectance *= image2D_average(textures[colorTexId]);
  }

  auto precomp = mi::fresnel_coat_precompute(mat.data[PLASTIC_ROUGHNESS], internal_ior, external_ior, diffuse_reflectance,
                                            {1.0f, 1.0f, 1.0f, 1.0f}, is_spectral_mode, spectrum);

  mat.data[PLASTIC_PRECOMP_REFLECTANCE] = precomp.internal_reflectance;
  mat.data[PLASTIC_SPEC_SAMPLE_WEIGHT]  = precomp.specular_sampling_weight;

  std::copy(precomp.transmittance.begin(), precomp.transmittance.end(), std::back_inserter(precomputed_transmittance));
  
  mat.datai[0] = (precomputed_transmittance.size() / MI_ROUGH_TRANSMITTANCE_RES) - 1u;

  return mat;
}

std::string Integrator::GetFeatureName(uint32_t a_featureId)
{
  switch(a_featureId)
  {
    case KSPEC_MAT_TYPE_GLTF      : return "GLTF_LITE";
    case KSPEC_MAT_TYPE_GLASS     : return "GLASS";
    case KSPEC_MAT_TYPE_CONDUCTOR : return "CONDUCTOR";
    case KSPEC_MAT_TYPE_THIN_FILM : return "THIN_FILM";
    case KSPEC_MAT_TYPE_DIFFUSE   : return "DIFFUSE";
    case KSPEC_MAT_TYPE_PLASTIC   : return "PLASTIC";
    case KSPEC_SPECTRAL_RENDERING : return "SPECTRAL";
    case KSPEC_MAT_TYPE_BLEND     : return "BLEND";
    case KSPEC_BUMP_MAPPING       : return "BUMP";
    case KSPEC_BLEND_STACK_SIZE   :
    {
      std::stringstream strout;
      strout << "STACK_SIZE = " << m_enabledFeatures[KSPEC_BLEND_STACK_SIZE];
      return strout.str();
    }
    default:
    break;
  };
  return "UNKNOWN";
}

hydra_xml::HydraScene   g_lastScene;
std::string Integrator::g_lastScenePath;
std::string Integrator::g_lastSceneDir;
SceneInfo   Integrator::g_lastSceneInfo;

static const std::wstring hydraOldMatTypeStr       {L"hydra_material"};
static const std::wstring roughConductorMatTypeStr {L"rough_conductor"};
static const std::wstring thinFilmMatTypeStr       {L"thin_film"};
static const std::wstring simpleDiffuseMatTypeStr  {L"diffuse"};
static const std::wstring blendMatTypeStr          {L"blend"};
static const std::wstring plasticMatTypeStr        {L"plastic"};
static const std::wstring dielectricMatTypeStr     {L"dielectric"};

std::vector<uint32_t> Integrator::PreliminarySceneAnalysis(const char* a_scenePath, const char* a_sncDir, SceneInfo* pSceneInfo)
{
  if(pSceneInfo == nullptr)
  {
    std::cout << "[Integrator::PreliminarySceneAnalysis]: nullptr pSceneInfo" << std::endl;
    exit(0);
  }

  g_lastSceneInfo.spectral = pSceneInfo->spectral;

  std::vector<uint32_t> features;
  
  std::string scenePathStr(a_scenePath);
  std::string sceneDirStr(a_sncDir);  
  auto loadRes = g_lastScene.LoadState(scenePathStr, sceneDirStr);
  if(loadRes != 0)
  {
    std::cout << "[Integrator::PreliminarySceneAnalysis]: Load scene xml failed: '" << a_scenePath << "'" << std::endl; 
    exit(0);
  }

  //// initial feature map
  //
  features.resize(TOTAL_FEATURES_NUM);  // disable all features by default
  for(auto& feature : features)         //
    feature = 0;                        //
  features[KSPEC_BLEND_STACK_SIZE]   = 1; // set smallest possible stack size for blends (i.e. blends are disabled!)
  features[KSPEC_SPECTRAL_RENDERING] = (pSceneInfo->spectral == 0) ? 0 : 1;
  
  //// list reauired material features
  //
  for(auto materialNode : g_lastScene.MaterialNodes())
  {
    auto mat_type = materialNode.attribute(L"type").as_string();
    if(mat_type == hydraOldMatTypeStr)
    {
      float4 transpColor = float4(0, 0, 0, 0);
      auto nodeTransp = materialNode.child(L"transparency");
      if (nodeTransp != nullptr)
        transpColor = GetColorFromNode(nodeTransp.child(L"color"), pSceneInfo->spectral);
  
      if(LiteMath::length3f(transpColor) > 1e-5f)
        features[KSPEC_MAT_TYPE_GLASS] = 1;
      else
        features[KSPEC_MAT_TYPE_GLTF] = 1;
    }
    else if(mat_type == roughConductorMatTypeStr)
    {
      features[KSPEC_MAT_TYPE_CONDUCTOR] = 1;
    }
    else if(mat_type == thinFilmMatTypeStr)
    {
      features[KSPEC_MAT_TYPE_THIN_FILM] = 1;
    }
    else if(mat_type == simpleDiffuseMatTypeStr)
    {
      features[KSPEC_MAT_TYPE_DIFFUSE] = 1;
    }
    else if(mat_type == blendMatTypeStr)
    {
      features[KSPEC_MAT_TYPE_BLEND]   = 1;
      features[KSPEC_BLEND_STACK_SIZE] = 4; // set appropriate stack size for blends
    }
    else if(mat_type == plasticMatTypeStr)
    {
      features[KSPEC_MAT_TYPE_PLASTIC] = 1;
    }
    else if(mat_type == dielectricMatTypeStr)
    {
      features[KSPEC_MAT_TYPE_DIELECTRIC] = 1;
    }

    if(materialNode.child(L"displacement") != nullptr)
      features[KSPEC_BUMP_MAPPING] = 1;
  }

  for(auto settings : g_lastScene.Settings())
  {
    g_lastSceneInfo.width  = settings.width;
    g_lastSceneInfo.height = settings.height;
    break; //take first render settings
  }

  g_lastScenePath = scenePathStr;
  g_lastSceneDir  = sceneDirStr;
  g_lastSceneInfo.maxMeshes            = 1024;
  g_lastSceneInfo.maxTotalVertices     = 4'000'000;
  g_lastSceneInfo.maxTotalPrimitives   = 4'000'000;
  g_lastSceneInfo.maxPrimitivesPerMesh = 1'000'000;

  g_lastSceneInfo.memTextures = 0;
  for(auto texNode : g_lastScene.TextureNodes())
  {
    uint32_t width  = texNode.attribute(L"width").as_uint();
    uint32_t height = texNode.attribute(L"height").as_uint();
    size_t byteSize = texNode.attribute(L"bytesize").as_ullong();

    if(width == 0 || height == 0)
    {
      width  = 256;
      height = 256;
      byteSize = 256*256*4;
    }

    //if(byteSize < width*height*4) // what if we have single channel 8 bit texture ...
    //  byteSize = width*height*4;  //

    g_lastSceneInfo.memTextures += uint64_t(byteSize);
  }

  g_lastSceneInfo.memGeom = 0;
  uint32_t meshesNum = 0;
  uint64_t maxTotalVertices = 0;
  uint64_t maxTotalPrimitives = 0;

  for(auto node : g_lastScene.GeomNodes())
  {
    const uint64_t byteSize = std::max<uint64_t>(node.attribute(L"bytesize").as_ullong(), 1024);
    const uint32_t vertNum  = node.attribute(L"vertNum").as_uint();
    const uint32_t trisNum  = node.attribute(L"triNum").as_uint();
    //const uint64_t byteSize = sizeof(float)*8*vertNum + trisNum*4*sizeof(uint32_t) + 1024;
    if(g_lastSceneInfo.maxPrimitivesPerMesh < trisNum)
      g_lastSceneInfo.maxPrimitivesPerMesh = trisNum;
    maxTotalVertices   += uint64_t(vertNum);
    maxTotalPrimitives += uint64_t(trisNum);
    meshesNum          += 1;
    g_lastSceneInfo.memGeom += byteSize;
  }

  g_lastSceneInfo.maxTotalVertices     = maxTotalVertices   + 1024*256;
  g_lastSceneInfo.maxTotalPrimitives   = maxTotalPrimitives + 1024*256;

  g_lastSceneInfo.memGeom     += uint64_t(4*1024*1024); // reserve mem for geom
  g_lastSceneInfo.memTextures += uint64_t(4*1024*1024); // reserve mem for tex

  (*pSceneInfo) = g_lastSceneInfo;
  return features;
}

std::vector<float> CreateSphericalTextureFromIES(const std::string& a_iesData, int* pW, int* pH);

bool Integrator::LoadScene(const char* a_scenePath, const char* a_sncDir)
{ 
  LoadSceneBegin();

  std::string scenePathStr(a_scenePath);
  std::string sceneDirStr(a_sncDir);  
  hydra_xml::HydraScene sceneLocal;
  
  const bool sameSceneAnalyzed = (scenePathStr == g_lastScenePath) && (sceneDirStr == g_lastSceneDir);
  hydra_xml::HydraScene& scene = sameSceneAnalyzed ? g_lastScene : sceneLocal;
  
  if(!sameSceneAnalyzed)
  {
    auto loadRes = scene.LoadState(scenePathStr, sceneDirStr);
    if(loadRes != 0)
    {
      std::cout << "Integrator::LoadScene failed: '" << a_scenePath << "'" << std::endl; 
      exit(0);
    }
  }

  //// init spectral curves
  m_cie_x      = Get_CIE_X();
  m_cie_y      = Get_CIE_Y();
  m_cie_z      = Get_CIE_Z();
  ////
  
  //// init render feature map
  m_actualFeatures.resize(TOTAL_FEATURES_NUM); // disable all features by default
  for(auto& feature : m_actualFeatures)              //
    feature = 0;                                      //
  m_actualFeatures[KSPEC_BLEND_STACK_SIZE] = 1;      // set smallest possible stack size for blends
  m_actualFeatures[KSPEC_SPECTRAL_RENDERING] = (m_spectral_mode == 0) ? 0 : 1;
  //// 

  std::vector<TextureInfo> texturesInfo;
  texturesInfo.resize(0);
  texturesInfo.reserve(100);

  #ifdef WIN32
  size_t endPos = scenePathStr.find_last_of("\\");
  if(endPos == std::string::npos)
    endPos = scenePathStr.find_last_of("/");
  #else
  size_t endPos = scenePathStr.find_last_of('/');
  #endif

  const std::string sceneFolder = (sceneDirStr == "") ? scenePathStr.substr(0, endPos) : sceneDirStr;

  //// (0) load textures info
  //
  for(auto texNode : scene.TextureNodes())
  {
    TextureInfo tex;
    tex.path   = std::wstring(sceneFolder.begin(), sceneFolder.end()) + L"/" + texNode.attribute(L"loc").as_string();
    tex.width  = texNode.attribute(L"width").as_uint();
    tex.height = texNode.attribute(L"height").as_uint();
    if(tex.width != 0 && tex.height != 0)
    {
      const size_t byteSize = texNode.attribute(L"bytesize").as_ullong();
      tex.bpp = uint32_t(byteSize / size_t(tex.width*tex.height));
    }
    texturesInfo.push_back(tex);
  }

  std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> texCache;
  texCache[HydraSampler()] = 0; // zero white texture
  
  m_textures.resize(0);
  m_textures.reserve(256);
  m_textures.push_back(MakeWhiteDummy());

  std::vector<SpectrumInfo> spectraInfo;
  spectraInfo.reserve(100);
  if(m_spectral_mode != 0)
  {  
    for(auto specNode : scene.SpectraNodes())
    {
      auto spec_id   = specNode.attribute(L"id").as_uint();
      auto spec_path = std::filesystem::path(sceneFolder);
      spec_path.append(specNode.attribute(L"loc").as_string());

      auto spec = LoadSPDFromFile(spec_path, spec_id);
      
      uint32_t offset = uint32_t(m_wavelengths.size());
      std::copy(spec.wavelengths.begin(), spec.wavelengths.end(), std::back_inserter(m_wavelengths));
      std::copy(spec.values.begin(), spec.values.end(), std::back_inserter(m_spec_values));
      m_spec_offset_sz.push_back(uint2{offset, static_cast<uint32_t>(spec.wavelengths.size())});
      
      // we expect dense, sorted ids for now
      //assert(m_spectra[spec_id].id == spec_id);
    }

    // if no spectra are loaded add uniform 1.0 spectrum
    if(spectraInfo.empty())
    {
      Spectrum uniform1;
      uniform1.id = 0;
      uniform1.wavelengths = {200.0f, 400.0f, 600.0f, 800.0f};
      uniform1.values = {1.0f, 1.0f, 1.0f, 1.0f};
      
      uint32_t offset = uint32_t(m_wavelengths.size());
      std::copy(uniform1.wavelengths.begin(), uniform1.wavelengths.end(), std::back_inserter(m_wavelengths));
      std::copy(uniform1.values.begin(), uniform1.values.end(), std::back_inserter(m_spec_values));
      m_spec_offset_sz.push_back(uint2{offset, static_cast<uint32_t>(uniform1.wavelengths.size())});
    }
  }

  // (1) load lights
  //
  m_instIdToLightInstId.resize(scene.GetInstancesNum(), -1);

  for(auto lightInst : scene.InstancesLights())
  {
    const std::wstring ltype = lightInst.lightNode.attribute(L"type").as_string();
    const std::wstring shape = lightInst.lightNode.attribute(L"shape").as_string();
    const std::wstring ldist = lightInst.lightNode.attribute(L"distribution").as_string();

    const float sizeX        = lightInst.lightNode.child(L"size").attribute(L"half_width").as_float();
    const float sizeZ        = lightInst.lightNode.child(L"size").attribute(L"half_length").as_float();
    float power              = lightInst.lightNode.child(L"intensity").child(L"multiplier").text().as_float();
    if (power == 0.0f) power = lightInst.lightNode.child(L"intensity").child(L"multiplier").attribute(L"val").as_float();
    if (power == 0.0f) power = 1.0f;

    float4 color = GetColorFromNode(lightInst.lightNode.child(L"intensity").child(L"color"), m_spectral_mode != 0);
    auto matrix  = lightInst.matrix;

    auto lightSpecId = GetSpectrumIdFromNode(lightInst.lightNode.child(L"intensity").child(L"color"));  

    LightSource lightSource{};
    lightSource.specId   = lightSpecId;
    lightSource.mult     = power;
    lightSource.distType = LIGHT_DIST_LAMBERT;
    lightSource.iesId    = uint(-1);
    lightSource.texId    = uint(-1);
    if(ltype == std::wstring(L"sky"))
    {
      m_envColor = color;
    }
    else if(ltype == std::wstring(L"directional"))
    {
      lightSource.pos       = lightInst.matrix * float4(0.0f, 0.0f, 0.0f, 1.0f);
      lightSource.norm      = normalize(lightInst.matrix * float4(0.0f, -1.0f, 0.0f, 0.0f));
      lightSource.intensity = color;
      lightSource.geomType  = LIGHT_GEOM_DIRECT;
    }
    else if(shape == L"rect" || shape == L"disk")
    {
      lightSource.pos       = lightInst.matrix * float4(0.0f, 0.0f, 0.0f, 1.0f);
      lightSource.norm      = normalize(lightInst.matrix * float4(0.0f, -1.0f, 0.0f, 0.0f));
      lightSource.intensity = color;
      lightSource.geomType  = (shape == L"rect") ? LIGHT_GEOM_RECT : LIGHT_GEOM_DISC;

      // extract scale and rotation from transformation matrix
      float3 scale;
      for(int i = 0; i < 3; ++i)
      {
        float4 vec = matrix.col(i);
        scale[i] = length3f(vec);
      }

      lightSource.matrix.set_col(0, matrix.get_col(2)); // why this matrix has swapped (x,z) ?
      lightSource.matrix.set_col(1, matrix.get_col(1)); // why this matrix has swapped (x,z) ?
      lightSource.matrix.set_col(2, matrix.get_col(0)); // why this matrix has swapped (x,z) ?
      lightSource.matrix.set_col(3, float4(0,0,0,1));
      lightSource.size = float2(sizeX, sizeZ);

      if(shape == L"disk")
      {
        lightSource.size.x = lightInst.lightNode.child(L"size").attribute(L"radius").as_float();
        lightSource.pdfA   = 1.0f / (LiteMath::M_PI * lightSource.size.x *lightSource.size.x * scale.x * scale.z);
      }
      else
        lightSource.pdfA   = 1.0f / (4.0f * lightSource.size.x * lightSource.size.y * scale.x * scale.z);
    }
    else if (shape == L"sphere")
    {
      float radius = lightInst.lightNode.child(L"size").attribute(L"radius").as_float();
      float3 scale; 
      for(int i = 0; i < 3; ++i)
      {
        float4 vec = matrix.col(i);
        scale[i] = length3f(vec);
      }

      radius = radius*scale.x; // support for uniform scale, assume scale.x == scale.y == scale.z
      if(std::abs(scale.x - scale.y) > 1e-5f || std::abs(scale.x - scale.z) > 1e-5f)
      {
        std::cout << "[Integrator::LoadScene]: ALERT!" << std::endl;
        std::cout << "[Integrator::LoadScene]: non uniform scale for spherical light instance matrix is not supported: (" << scale.x << ", " << scale.y << ", " << scale.z << ")" << std::endl; 
      }

      lightSource.pos       = lightInst.matrix * float4(0.0f, 0.0f, 0.0f, 1.0f);
      lightSource.norm      = float4(0.0f, -1.0f, 0.0f, 0.0f);
      lightSource.intensity = color;
      lightSource.geomType  = LIGHT_GEOM_SPHERE;

      lightSource.matrix    = float4x4{};
      lightSource.size      = float2(radius, radius);
      lightSource.pdfA      = 1.0f / (4.0f*LiteMath::M_PI*radius*radius);
    }
    else if (shape == L"point")
    {
      lightSource.pos       = lightInst.matrix * float4(0.0f, 0.0f, 0.0f, 1.0f);
      lightSource.norm      = normalize(lightInst.matrix * float4(0.0f, -1.0f, 0.0f, 0.0f));
      lightSource.intensity = color;
      lightSource.geomType  = LIGHT_GEOM_POINT;
      lightSource.distType  = (ldist == L"uniform" || ldist == L"omni" || ldist == L"ies") ? LIGHT_DIST_OMNI : LIGHT_DIST_LAMBERT;
      lightSource.pdfA      = 1.0f;
      lightSource.size      = float2(0,0);
      lightSource.matrix    = float4x4{};
    }

    auto iesNode = lightInst.lightNode.child(L"ies");
    if(iesNode != nullptr)
    {
      const std::wstring iesFileW = std::wstring(sceneFolder.begin(), sceneFolder.end()) + L"/" + iesNode.attribute(L"loc").as_string();\
      const std::string  iesFileA = hydra_xml::ws2s(iesFileW);
      
      int w,h;
      std::vector<float> sphericalTexture = CreateSphericalTextureFromIES(iesFileA.c_str(), &w, &h);
      
      // normalize ies texture
      //
      float maxVal = 0.0f;
      for (auto i = 0; i < sphericalTexture.size(); i++)
        maxVal = std::max(maxVal, sphericalTexture[i]);
  
      if(maxVal == 0.0f)
      {
        std::cerr << "[ERROR]: broken IES file (maxVal = 0.0): " << iesFileA.c_str() << std::endl;
        maxVal = 1.0f;
      }
  
      float invMax = 1.0f / maxVal;
      for (auto i = 0; i < sphericalTexture.size(); i++)
      {
        float val = invMax*sphericalTexture[i];
        sphericalTexture[i] = val;
      }
      ////

      auto pTexture = std::make_shared< Image2D<float> >(w, h, sphericalTexture.data());
      pTexture->setSRGB(false);
      
      Sampler sampler;
      sampler.filter   = Sampler::Filter::LINEAR; 
      sampler.addressU = Sampler::AddressMode::WRAP;
      sampler.addressV = Sampler::AddressMode::WRAP;

      m_textures.push_back(MakeCombinedTexture2D(pTexture, sampler));
      lightSource.iesId = uint(m_textures.size()-1);
      //std::cout << "lightSource.iesId = " << lightSource.iesId << std::endl;
    }

    m_lights.push_back(lightSource);
  }

  //// (2) load materials
  //
  m_materials.resize(0);
  m_materials.reserve(100);

  for(auto materialNode : scene.MaterialNodes())
  {
    Material mat = {};
    auto mat_type = materialNode.attribute(L"type").as_string();
    
    mat.data[EMISSION_MULT] = 1.0f;

    if(mat_type == hydraOldMatTypeStr)
    {
      mat = ConvertOldHydraMaterial(materialNode, texturesInfo, texCache, m_textures, m_spectral_mode);
      if(mat.mtype == MAT_TYPE_GLASS)
        m_actualFeatures[KSPEC_MAT_TYPE_GLASS] = 1;
      else
        m_actualFeatures[KSPEC_MAT_TYPE_GLTF] = 1;
    }
    else if(mat_type == roughConductorMatTypeStr)
    {
      mat = LoadRoughConductorMaterial(materialNode, texturesInfo, texCache, m_textures);
      m_actualFeatures[KSPEC_MAT_TYPE_CONDUCTOR] = 1;
    }
    else if (mat_type == thinFilmMatTypeStr)
    {
      mat = LoadThinFilmMaterial(materialNode, texturesInfo, texCache, m_textures, m_precomp_thin_films, m_films_thickness_vec, m_films_spec_id_vec, m_films_eta_k_vec,
                                 m_spec_values, m_wavelengths, m_spec_offset_sz);
      m_actualFeatures[KSPEC_MAT_TYPE_THIN_FILM] = 1;
    }
    else if(mat_type == simpleDiffuseMatTypeStr)
    {
      mat = LoadDiffuseMaterial(materialNode, texturesInfo, texCache, m_textures, m_spectral_mode);
      m_actualFeatures[KSPEC_MAT_TYPE_DIFFUSE] = 1;
    }
    else if(mat_type == blendMatTypeStr)
    {
      mat = LoadBlendMaterial(materialNode, texturesInfo, texCache, m_textures);
      m_actualFeatures[KSPEC_MAT_TYPE_BLEND]   = 1;
      m_actualFeatures[KSPEC_BLEND_STACK_SIZE] = 4; // set appropriate stack size for blends
    }
    else if(mat_type == plasticMatTypeStr)
    {
      mat = LoadPlasticMaterial(materialNode, texturesInfo, texCache, m_textures, m_precomp_coat_transmittance, m_spectral_mode,
                                m_spec_values, m_wavelengths, m_spec_offset_sz);
      m_actualFeatures[KSPEC_MAT_TYPE_PLASTIC] = 1;
    }
    else if(mat_type == dielectricMatTypeStr)
    {
      mat = LoadDielectricMaterial(materialNode, texturesInfo, texCache, m_textures, m_spectral_mode);
      m_actualFeatures[KSPEC_MAT_TYPE_DIELECTRIC] = 1;
    }

    if(materialNode.attribute(L"light_id") != nullptr)
    {
      int lightId = materialNode.attribute(L"light_id").as_int();
      if(lightId >= 0 && lightId < static_cast<int>(m_lights.size()))
      {
        auto tmp = mat.colors[EMISSION_COLOR] != m_lights[lightId].intensity;
        if(tmp.x == 0xFFFFFFFF && tmp.y == 0xFFFFFFFF && tmp.z == 0xFFFFFFFF && tmp.w == 0xFFFFFFFF)
          std::cout << "Color in material for light geom and color in light intensity node are different! " 
                    << "Using values from light intensity node. lightId = " << lightId << std::endl;

        mat.colors[EMISSION_COLOR] = m_lights[lightId].intensity;

        if(mat.data[EMISSION_MULT] != m_lights[lightId].mult)
          std::cout << "Color multiplier in material for light geom and in light intensity node are different! " 
                    << "Using values from light intensity node. lightId = " << lightId << std::endl;

        mat.data[EMISSION_MULT] = m_lights[lightId].mult;

        if(mat.spdid[0] != m_lights[lightId].specId)
          std::cout << "Spectrum in material for light geom and in light intensity node are different! " 
                    << "Using values from light intensity node. lightId = " << lightId << std::endl;
        
        mat.spdid[0] = m_lights[lightId].specId;
      }
    }

    // setup normal map
    //
    mat.texid[1] = 0xFFFFFFFF;

    if(materialNode.child(L"displacement") != nullptr)
    {
      auto dispNode = materialNode.child(L"displacement");
      if(dispNode.attribute(L"type").as_string() != std::wstring(L"normal_bump"))
      {
        std::string bumpType = hydra_xml::ws2s(dispNode.attribute(L"type").as_string());
        std::cout << "[Integrator::LoadScene]: bump type '" << bumpType.c_str() << "' is not supported! only 'normal_bump' is allowed."  << std::endl;
      }
      else
      {
        auto normalNode  = dispNode.child(L"normal_map");
        auto invertNode  = normalNode.child(L"invert");   // todo read swap flag also
        auto textureNode = normalNode.child(L"texture");

        const auto& [sampler, texID] = LoadTextureFromNode(normalNode, texturesInfo, texCache, m_textures);

        mat.row0 [1] = sampler.row0;
        mat.row1 [1] = sampler.row1;
        mat.texid[1] = texID;

        const bool invertX = (invertNode.attribute(L"x").as_int() == 1);
        const bool invertY = (invertNode.attribute(L"y").as_int() == 1);
        const bool swapXY  = (invertNode.attribute(L"swap_xy").as_int() == 1);

        if(invertX)
          mat.cflags |= uint(FLAG_NMAP_INVERT_X);
        if(invertY)
          mat.cflags |= uint(FLAG_NMAP_INVERT_Y);
        if(swapXY)
          mat.cflags |= uint(FLAG_NMAP_SWAP_XY);
        
        m_actualFeatures[KSPEC_BUMP_MAPPING] = 1; // enable bump mapping feature
      }
    }

    m_materials.push_back(mat);
  }

  // load first camera and update matrix
  //
  for(auto cam : scene.Cameras())
  {
    float aspect   = float(m_winWidth) / float(m_winHeight);
    auto proj      = perspectiveMatrix(cam.fov, aspect, cam.nearPlane, cam.farPlane);
    auto worldView = lookAt(float3(cam.pos), float3(cam.lookAt), float3(cam.up));
      
    m_exposureMult = cam.exposureMult;
    m_proj         = proj;
    m_worldView    = worldView;
    m_projInv      = inverse4x4(proj);
    m_worldViewInv = inverse4x4(worldView);

    auto sensorNode = cam.node.child(L"sensor");
    if(sensorNode != nullptr)
    {
      auto responceNode = sensorNode.child(L"response");
      if(responceNode != nullptr)
      {
        std::wstring responceType = responceNode.attribute(L"type").as_string();
        if(responceType == L"xyz" || responceType == L"XYZ")
          m_camResponseType = CAM_RESPONCE_XYZ;
        else
          m_camResponseType = CAM_RESPONCE_RGB;
        
        int id = 0;
        for(auto spec : responceNode.children(L"spectrum")) {
          m_camResponseSpectrumId[id] = spec.attribute(L"id").as_int();
          id++;
          if(id >= 3)
            break;
        }

        m_camRespoceRGB   = GetColorFromNode(responceNode.child(L"color"), false);
        m_camRespoceRGB.w = 1.0f;
      }
    }

    break; // take first cam
  }

  
  //// (2) load meshes
  //
  m_matIdOffsets.reserve(1024);
  m_vertOffset.reserve(1024);
  m_matIdByPrimId.reserve(128000);
  m_triIndices.reserve(128000*3);

  m_vNorm4f.resize(0);
  m_vTang4f.resize(0);
  //m_vTexc2f.resize(0);

  m_pAccelStruct->ClearGeom();
  for(auto meshPath : scene.MeshFiles())
  {
    #ifdef _DEBUG
    std::cout << "[LoadScene]: mesh = " << meshPath.c_str() << std::endl;
    #endif
    auto currMesh = cmesh4::LoadMeshFromVSGF(meshPath.c_str());
    auto geomId   = m_pAccelStruct->AddGeom_Triangles3f((const float*)currMesh.vPos4f.data(), currMesh.vPos4f.size(), currMesh.indices.data(), currMesh.indices.size(), BUILD_HIGH, sizeof(float)*4);

    (void)geomId; // silence unused var. warning

    m_matIdOffsets.push_back(static_cast<unsigned int>(m_matIdByPrimId.size()));
    m_vertOffset.push_back(static_cast<unsigned int>(m_vNorm4f.size()));
    const size_t lastVertex = m_vNorm4f.size();

    m_matIdByPrimId.insert(m_matIdByPrimId.end(), currMesh.matIndices.begin(), currMesh.matIndices.end() );
    m_triIndices.insert(m_triIndices.end(), currMesh.indices.begin(), currMesh.indices.end());

    m_vNorm4f.insert(m_vNorm4f.end(), currMesh.vNorm4f.begin(), currMesh.vNorm4f.end());
    m_vTang4f.insert(m_vTang4f.end(), currMesh.vTang4f.begin(), currMesh.vTang4f.end());

    for(size_t i = 0; i<currMesh.VerticesNum(); i++) {          // pack texture coords
      m_vNorm4f[lastVertex + i].w = currMesh.vTexCoord2f[i].x;
      m_vTang4f[lastVertex + i].w = currMesh.vTexCoord2f[i].y;
    }
    //m_vTexc2f.insert(m_vTexc2f.end(), currMesh.vTexCoord2f.begin(), currMesh.vTexCoord2f.end()); // #TODO: store quantized texture coordinates
  }

  //// (3) make instances of created meshes
  //
  m_normMatrices.clear(); m_normMatrices.reserve(1000);
  m_remapInst.clear();    m_remapInst.reserve(1000);
  m_pAccelStruct->ClearScene();
  uint32_t realInstId = 0;
  for(auto inst : scene.InstancesGeom())
  {
    if(inst.instId != realInstId)
    {
      #ifdef _DEBUG
      std::cout << "[Integrator::LoadScene]: WARNING, bad instance id: written in xml: inst.instId is '" <<  inst.instId << "', realInstId by node order is '" << realInstId << "'" << std::endl;
      std::cout << "[Integrator::LoadScene]: -->      instances must be written in a sequential order, perform 'inst.instId = realInstId'" << std::endl;
      #endif
      inst.instId = realInstId;
    }
    m_pAccelStruct->AddInstance(inst.geomId, inst.matrix);
    m_normMatrices.push_back(transpose(inverse4x4(inst.matrix)));
    m_remapInst.push_back(inst.rmapId);

    m_instIdToLightInstId[inst.instId] = inst.lightInstId;
    realInstId++;
  }

  m_pAccelStruct->CommitScene(); // to enable more anync may call CommitScene later, but need acync API: CommitSceneStart() ... CommitSceneFinish()
  
  // (4) load remap lists and put all of the to the flat data structure
  // 
  m_allRemapLists.clear();
  m_allRemapListsOffsets.clear();
  m_allRemapLists.reserve(m_normMatrices.size()*10);     // approx size for all reamp lists based on number of instances; may not do this reserve in fact, or make it more precise
  m_allRemapListsOffsets.reserve(m_normMatrices.size()); // approx size for all reamp lists ... 
  for(auto remapList : scene.RemapLists())
  {
    m_allRemapListsOffsets.push_back(static_cast<int>(m_allRemapLists.size()));
    m_allRemapLists.insert(m_allRemapLists.end(), remapList.begin(), remapList.end());
  }
  m_allRemapListsOffsets.push_back(static_cast<int>(m_allRemapLists.size())); // put size of the list remap list
  
  // (5) load render settings
  //
  for(const auto& sett : scene.Settings())
  {
    m_traceDepth = sett.depth;
    m_spp        = sett.spp;

    if(m_traceDepth == 0)
      m_traceDepth = 6;
    
    if(m_spp == 0)
      m_spp = 1;

    break; // take first render settings
  }

  // (6) print enabled features in scene
  //
  for(size_t i=0; i<m_enabledFeatures.size();i++)
  {
    if(m_actualFeatures[i] != m_enabledFeatures[i])
    {
      std::string featureName = GetFeatureName(uint32_t(i));
      std::cout << "[Integrator::LoadScene]: feature '" << featureName.c_str() << "' has different values in 'enabled' and 'actual' features array" << std::endl;
      std::cout << "[Integrator::LoadScene]: enabled = " << m_enabledFeatures[i] << ", actual = " << m_actualFeatures[i] << std::endl;
    }
  }

  std::cout << "features = {";
  bool firstFeature = true;
  for(size_t i=0; i<m_enabledFeatures.size();i++)
  {
    if(m_enabledFeatures[i] != 0)
    {
      std::string featureName = GetFeatureName(uint32_t(i));
      if(!firstFeature)
        std::cout << ",";
      std::cout << featureName.c_str();
      firstFeature = false;
    }
  }
  std::cout << "};" << std::endl;

  LoadSceneEnd();
  return true;
}
