#include "integrator_pt.h"
#include "include/crandom.h"

#include "LiteScene/cmesh4.h"
using cmesh4::SimpleMesh;

#include "mi_materials.h"
#include "spectrum.h"

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

std::shared_ptr<ICombinedImageSampler> LoadTextureAndMakeCombined(const TextureInfo& a_texInfo, const Sampler& a_sampler)
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
  if(a_texInfo.bpp == 16)
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
    pTexture->setSRGB(true);
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
    texId  = node.child(L"texture").attribute(L"id").as_uint();
    textures.push_back(LoadTextureAndMakeCombined(texturesInfo[texId], sampler.sampler));
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
  mat.data[UINT_MTYPE]         = as_float(MAT_TYPE_GLTF);
  mat.data[GLTF_FLOAT_ALPHA]   = 0.0f;
  mat.colors[GLTF_COLOR_COAT]  = float4(1,1,1,1); 
  mat.colors[GLTF_COLOR_METAL] = float4(0,0,0,0);  
  mat.data[UINT_LIGHTID]       = as_float(uint(-1));
  
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
    
    mat.row0 [0]  = emissiveSampler.row0;
    mat.row1 [0]  = emissiveSampler.row1;
    mat.data[EMISSION_TEXID0] = as_float(texID);
    
    mat.colors[EMISSION_COLOR] = color;
    if(materialNode.attribute(L"light_id") == nullptr)
      mat.data[UINT_LIGHTID] = as_float(uint(-1));
    else
      mat.data[UINT_LIGHTID] = as_float(uint(materialNode.attribute(L"light_id").as_int())); // for correct process of "-1"

    auto specId = GetSpectrumIdFromNode(nodeEmissColor);  
    mat.data[EMISSION_SPECID0] = as_float(specId);
    mat.data[UINT_MTYPE] = as_float(MAT_TYPE_LIGHT_SOURCE);

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
    
    mat.row0 [0]  = diffSampler.row0;
    mat.row1 [0]  = diffSampler.row1;
    mat.data[GLTF_UINT_TEXID0] = as_float(texID);
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
    mat.data[UINT_MTYPE]         = as_float(MAT_TYPE_GLTF);
    mat.colors[GLTF_COLOR_BASE]  = color;
    mat.colors[GLTF_COLOR_COAT]  = reflColor;
    mat.data[UINT_LIGHTID]       = as_float(uint(-1));

    if(hasFresnel)
    {
      mat.data[GLTF_FLOAT_ALPHA]   = 0.0f;
      mat.colors[GLTF_COLOR_COAT]  = reflColor;
      mat.colors[GLTF_COLOR_METAL] = float4(0,0,0,0); 
      mat.data[UINT_CFLAGS]        = as_float(GLTF_COMPONENT_LAMBERT | GLTF_COMPONENT_COAT);
      SetMiPlastic(&mat, fresnelIOR, 1.0f, color, reflColor);
    }
    else
    {
      mat.data[GLTF_FLOAT_ALPHA]   = length(reflColor)/( length(reflColor) + length3f(color) );
      mat.colors[GLTF_COLOR_COAT]  = float4(0,0,0,0); 
      mat.colors[GLTF_COLOR_METAL] = reflColor;   // disable coating for such blend type
      mat.data[UINT_CFLAGS]        = as_float(GLTF_COMPONENT_LAMBERT | GLTF_COMPONENT_METAL);
    }
  }
  else if(length(reflColor) > 1e-5f)
  {
    mat.data[UINT_MTYPE]           = as_float(MAT_TYPE_GLTF);
    mat.data[UINT_CFLAGS]          = as_float(GLTF_COMPONENT_METAL);
    mat.colors[GLTF_COLOR_METAL]   = reflColor;
    mat.colors[GLTF_COLOR_COAT]    = float4(0,0,0,0); 
    mat.data[GLTF_FLOAT_ALPHA]     = 1.0f;
  }
  else if(length(to_float3(color)) > 1e-5f)
  {
    mat.data[UINT_MTYPE]         = as_float(MAT_TYPE_GLTF);
    mat.data[UINT_CFLAGS]        = as_float(GLTF_COMPONENT_LAMBERT);
    mat.colors[GLTF_COLOR_BASE]  = color;
    mat.colors[GLTF_COLOR_COAT]  = float4(0,0,0,0); 
    mat.colors[GLTF_COLOR_METAL] = float4(0,0,0,0);    
    mat.data[GLTF_FLOAT_ALPHA]   = 0.0f;
  }
    
  // Glass
  if (length(transpColor) > 1e-5f)
  {
    mat.data[UINT_MTYPE]                = as_float(MAT_TYPE_GLASS);      
    mat.colors[GLASS_COLOR_REFLECT]     = reflColor;
    mat.colors[GLASS_COLOR_TRANSP]      = transpColor;      
    mat.data[GLASS_FLOAT_GLOSS_REFLECT] = reflGlossiness;
    mat.data[GLASS_FLOAT_GLOSS_TRANSP]  = transpGlossiness;
    mat.data[GLASS_FLOAT_IOR]           = fresnelIOR;
  }

  if(is_emission_color)
    mat.data[UINT_MTYPE] = as_float(MAT_TYPE_LIGHT_SOURCE);

  auto nodeDiffRough = materialNode.child(L"diffuse").child(L"roughness");

  if (nodeDiffRough != nullptr)
  {
    mat.data[GLTF_FLOAT_ROUGH_ORENNAYAR] = hydra_xml::readval1f(nodeDiffRough);
    mat.data[UINT_CFLAGS] = as_float(as_uint(mat.data[UINT_CFLAGS]) | GLTF_COMPONENT_ORENNAYAR);
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
  mat.data[UINT_MTYPE]         = as_float(MAT_TYPE_CONDUCTOR);  
  mat.data[UINT_LIGHTID]       = as_float(uint(-1));

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
    
    mat.row0 [0]  = sampler.row0;
    mat.row1 [0]  = sampler.row1;
    mat.data[CONDUCTOR_TEXID0] = as_float(texID);
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

  mat.data[CONDUCTOR_ETA_SPECID] = as_float(etaSpecId);
  mat.data[CONDUCTOR_K_SPECID]   = as_float(kSpecId);

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
  mat.data[UINT_MTYPE]        = as_float(MAT_TYPE_DIFFUSE);  
  mat.data[UINT_LIGHTID]      = as_float(uint(-1));
  mat.data[DIFFUSE_ROUGHNESS] = 0.0f;
  mat.data[DIFFUSE_TEXID0]    = as_float(0);
  mat.data[DIFFUSE_SPECID]    = as_float(uint(-1));

  static const std::wstring orenNayarNameStr {L"oren-nayar"};

  auto bsdfType = materialNode.child(L"bsdf").attribute(L"type").as_string();
  if(bsdfType == orenNayarNameStr)
  {
    mat.data[UINT_CFLAGS] = as_float(GLTF_COMPONENT_ORENNAYAR);
  
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
    mat.data[DIFFUSE_TEXID0] = as_float(texID);

    auto specId = GetSpectrumIdFromNode(nodeColor);
    mat.data[DIFFUSE_SPECID] = as_float(specId);
  }

  return mat;
}


Material LoadBlendMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                           std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                           std::vector< std::shared_ptr<ICombinedImageSampler> > &textures)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  Material mat = {};
  mat.data[UINT_MTYPE]        = as_float(MAT_TYPE_BLEND);  
  mat.data[UINT_LIGHTID]      = as_float(uint(-1));
  mat.data[BLEND_WEIGHT]      = 1.0f;
  mat.data[BLEND_TEXID0]      = as_float(0);

  mat.data[BLEND_MAT_ID_1]    = as_float(materialNode.child(L"bsdf_1").attribute(L"id").as_uint());
  mat.data[BLEND_MAT_ID_2]    = as_float(materialNode.child(L"bsdf_2").attribute(L"id").as_uint());

  auto nodeWeight = materialNode.child(L"weight");
  if(nodeWeight != nullptr)
  {
    mat.data[BLEND_WEIGHT] = (hydra_xml::readval1f(nodeWeight), 1.0f);

    const auto& [sampler, texID] = LoadTextureFromNode(nodeWeight, texturesInfo, texCache, textures);
    
    mat.row0 [0]  = sampler.row0;
    mat.row1 [0]  = sampler.row1;
    mat.data[BLEND_TEXID0] = as_float(texID);
  }

  return mat;
}

std::string Integrator::GetFeatureName(uint32_t a_featureId)
{
  switch(a_featureId)
  {
    case KSPEC_MAT_TYPE_GLTF      : return "GLTF_LITE";
    case KSPEC_MAT_TYPE_GLASS     : return "GLASS";
    case KSPEC_MAT_TYPE_CONDUCTOR : return "CONDUCTOR";
    case KSPEC_MAT_TYPE_DIFFUSE   : return "DIFFUSE";
    case KSPEC_SOME_FEATURE_DUMMY : return "DUMMY";
    case KSPEC_SPECTRAL_RENDERING : return "SPECTRAL";
    case KSPEC_MAT_TYPE_BLEND     : return "BLEND";
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
static const std::wstring hydraOldMatTypeStr       {L"hydra_material"};
static const std::wstring roughConductorMatTypeStr {L"rough_conductor"};
static const std::wstring simpleDiffuseMatTypeStr  {L"diffuse"};
static const std::wstring blendMatTypeStr          {L"blend"};

std::vector<uint32_t> Integrator::PreliminarySceneAnalysis(const char* a_scenePath, const char* a_sncDir,
                                                           int& width, int& height, int& spectral_mode)
{
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
  features.resize(TOTAL_FEATURES_NUM); // disable all features by default
  for(auto& feature : features)              //
    feature = 0;                             //
  features[KSPEC_BLEND_STACK_SIZE] = 1;      // set smallest possible stack size for blends (i.e. blends are disabled!)
  features[KSPEC_SPECTRAL_RENDERING] = (spectral_mode == 0) ? 0 : 1;
  
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
        transpColor = GetColorFromNode(nodeTransp.child(L"color"), spectral_mode);
  
      if(LiteMath::length3f(transpColor) > 1e-5f)
        features[KSPEC_MAT_TYPE_GLASS] = 1;
      else
        features[KSPEC_MAT_TYPE_GLTF] = 1;
    }
    else if(mat_type == roughConductorMatTypeStr)
    {
      features[KSPEC_MAT_TYPE_CONDUCTOR] = 1;
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
  }

  for(auto settings : g_lastScene.Settings())
  {
    width  = settings.width;
    height = settings.height;
    break; //take ferst render settings
  }

  g_lastScenePath = scenePathStr;
  g_lastSceneDir  = sceneDirStr;

  return features;
}

bool Integrator::LoadScene(const char* a_scenePath, const char* a_sncDir)
{ 
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
      
      uint32_t offset = m_wavelengths.size();
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
      
      uint32_t offset = m_wavelengths.size();
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
      lightSource.distType  = (ldist == L"uniform" || ldist == L"omni") ? LIGHT_DIST_OMNI : LIGHT_DIST_LAMBERT;
      lightSource.pdfA      = 1.0f;
      lightSource.size      = float2(0,0);
      lightSource.matrix    = float4x4{};
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
      mat = ConvertOldHydraMaterial(materialNode, texturesInfo, texCache, m_textures, m_spectral_mode != 0);
      if(as_uint(mat.data[UINT_MTYPE]) == MAT_TYPE_GLASS)
        m_actualFeatures[KSPEC_MAT_TYPE_GLASS] = 1;
      else
        m_actualFeatures[KSPEC_MAT_TYPE_GLTF] = 1;
    }
    else if(mat_type == roughConductorMatTypeStr)
    {
      mat = LoadRoughConductorMaterial(materialNode, texturesInfo, texCache, m_textures);
      m_actualFeatures[KSPEC_MAT_TYPE_CONDUCTOR] = 1;
    }
    else if(mat_type == simpleDiffuseMatTypeStr)
    {
      mat = LoadDiffuseMaterial(materialNode, texturesInfo, texCache, m_textures, m_spectral_mode != 0);
      m_actualFeatures[KSPEC_MAT_TYPE_DIFFUSE] = 1;
    }
    else if(mat_type == blendMatTypeStr)
    {
      mat = LoadBlendMaterial(materialNode, texturesInfo, texCache, m_textures);
      m_actualFeatures[KSPEC_MAT_TYPE_BLEND]   = 1;
      m_actualFeatures[KSPEC_BLEND_STACK_SIZE] = 4; // set appropriate stack size for blends
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

        if(as_uint(mat.data[EMISSION_SPECID0]) != m_lights[lightId].specId)
          std::cout << "Spectrum in material for light geom and in light intensity node are different! " 
                    << "Using values from light intensity node. lightId = " << lightId << std::endl;

        mat.data[EMISSION_SPECID0] = as_float(m_lights[lightId].specId);
      }
    }
    
    // setup normal map
    //
    mat.data[UINT_NMAP_ID] = as_float(0xFFFFFFFF);
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
        auto invertNode  = normalNode.child(L"invert");
        auto textureNode = normalNode.child(L"texture");

        const bool invertX = (invertNode.attribute(L"x").as_int() == 1);
        const bool invertY = (invertNode.attribute(L"y").as_int() == 1);
        
        const auto& [sampler, texID] = LoadTextureFromNode(normalNode, texturesInfo, texCache, m_textures); // todo: gamma?
        
        mat.row0[1] = sampler.row0;
        mat.row1[1] = sampler.row1;
        mat.data[UINT_NMAP_ID] = as_float(texID);
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
    std::cout << "[LoadScene]: mesh = " << meshPath.c_str() << std::endl;
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
      std::cout << "[Integrator::LoadScene]: WARNING, bad instance id: written in xml: inst.instId is '" <<  inst.instId << "', realInstId by node order is '" << realInstId << "'" << std::endl;
      std::cout << "[Integrator::LoadScene]: -->      instances must be written in a sequential order, perform 'inst.instId = realInstId'" << std::endl;
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
  
  // we should not leave empty vectors for data which are used on GPU, kernel slicer does not handle this yet
  //
  if(m_spec_offset_sz.capacity() == 0)
    m_spec_offset_sz.reserve(16);
  if(m_spec_values.capacity() == 0)
    m_spec_values.reserve(16);
  if(m_wavelengths.capacity() == 0)
    m_wavelengths.reserve(16);
  return true;
}
