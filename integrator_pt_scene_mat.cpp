#include "integrator_pt_scene.h"

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
    ss >> res.m_col[i].x >> res.m_col[i].y >> res.m_col[i].z >> res.m_col[i].w;
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

void LoadSpectralTextures(const uint32_t specId,
                          const std::vector<TextureInfo> &texturesInfo,
                          std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                          std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                          std::vector<uint2> &spec_tex_ids_wavelengths,
                          const std::vector<uint2> &spec_tex_offset_sz, 
                          std::set<uint32_t> &loadedSpectralTextures)
{
  auto textures_sz = spec_tex_offset_sz[specId].y;
  auto offset = spec_tex_offset_sz[specId].x;
  if(textures_sz > 0 && loadedSpectralTextures.count(specId) == 0)
  {
    for(uint32_t i = 0; i < textures_sz; ++i)
    {
      uint32_t xml_tex_id = spec_tex_ids_wavelengths[offset + i].x;

      //TODO: put sampler somewhere in XML
      HydraSampler sampler;
      sampler.inputGamma = 1.0f;
      sampler.texId = xml_tex_id;

      const auto& [sampler_out, loaded_tex_id] = LoadTextureById(xml_tex_id, texturesInfo, sampler, texCache, textures);
      spec_tex_ids_wavelengths[offset + i].x = loaded_tex_id;     
    }
    loadedSpectralTextures.insert(specId);
  }
}


Material ConvertGLTFMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                             std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                             std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                             bool is_spectral_mode)
{
  std::wstring name              = materialNode.attribute(L"name").as_string();
  Material mat                   = {};
  mat.mtype                      = MAT_TYPE_GLTF;
  mat.cflags                     = GLTF_COMPONENT_LAMBERT | GLTF_COMPONENT_COAT;
  mat.data[GLTF_FLOAT_ALPHA]     = 0.0f;
  mat.data[GLTF_FLOAT_REFL_COAT] = 1.0f;
  mat.colors[GLTF_COLOR_COAT]    = float4(1,1,1,1); 
  mat.colors[GLTF_COLOR_METAL]   = float4(1,1,1,1); 

  for(int i=0;i<4;i++) {
    mat.row0 [i] = float4(1,0,0,0);
    mat.row1 [i] = float4(0,1,0,0);
    mat.texid[i] = 0;
  }

  float fresnelIOR     = 1.5f;
  float reflGlossiness = 1.0f;
  float metalness      = 0.0f;
  float4 baseColor(1,1,1,1);
  if(materialNode != nullptr)
  {
    if(materialNode.child(L"color") != nullptr) {
      baseColor = GetColorFromNode(materialNode.child(L"color"), is_spectral_mode);
      if(materialNode.child(L"color").child(L"texture") != nullptr) {
        const auto& [sampler, texID] = LoadTextureFromNode(materialNode.child(L"color"), texturesInfo, texCache, textures);
        mat.row0 [0] = sampler.row0;
        mat.row1 [0] = sampler.row1;
        mat.texid[0] = texID;
      }
    }

    if(materialNode.child(L"glossiness") != nullptr)
    {
      reflGlossiness = hydra_xml::readval1f(materialNode.child(L"glossiness"));  
      if(materialNode.child(L"glossiness").child(L"texture") != nullptr) {
        const auto& [sampler, texID] = LoadTextureFromNode(materialNode.child(L"glossiness"), texturesInfo, texCache, textures);
        mat.row0 [2] = sampler.row0;
        mat.row1 [2] = sampler.row1;
        mat.texid[2] = texID;
        mat.cflags |= FLAG_FOUR_TEXTURES;
      }
    }
    else if (materialNode.child(L"roughness") != nullptr)
    {
      reflGlossiness = hydra_xml::readval1f(materialNode.child(L"roughness")); 
      mat.cflags |= FLAG_INVERT_GLOSINESS; 
      if(materialNode.child(L"roughness").child(L"texture") != nullptr) {
        const auto& [sampler, texID] = LoadTextureFromNode(materialNode.child(L"roughness"), texturesInfo, texCache, textures);
        mat.row0 [2] = sampler.row0;
        mat.row1 [2] = sampler.row1;
        mat.texid[2] = texID;
        mat.cflags |= FLAG_FOUR_TEXTURES;
      }
    }
    
    if(materialNode.child(L"metalness") != nullptr)  
    {
      metalness = hydra_xml::readval1f(materialNode.child(L"metalness"));
      if(materialNode.child(L"metalness").child(L"texture") != nullptr) {
        const auto& [sampler, texID] = LoadTextureFromNode(materialNode.child(L"metalness"), texturesInfo, texCache, textures);
        mat.row0 [3] = sampler.row0;
        mat.row1 [3] = sampler.row1;
        mat.texid[3] = texID;
        mat.cflags |= FLAG_FOUR_TEXTURES;
      }
    }

    if(materialNode.child(L"fresnel_ior") != nullptr)  
      fresnelIOR     = hydra_xml::readval1f(materialNode.child(L"fresnel_ior"));
    if(materialNode.child(L"coat") != nullptr)  
      mat.data[GLTF_FLOAT_REFL_COAT] = hydra_xml::readval1f(materialNode.child(L"coat"));

    if(materialNode.child(L"glossiness_metalness_coat") != nullptr)
    {
      const float val = hydra_xml::readval1f(materialNode.child(L"glossiness_metalness_coat"));  
      metalness       = val;
      reflGlossiness  = val;
      mat.data[GLTF_FLOAT_REFL_COAT] = val;
      if(materialNode.child(L"glossiness_metalness_coat").child(L"texture") != nullptr) {
        const auto& [sampler, texID] = LoadTextureFromNode(materialNode.child(L"glossiness_metalness_coat"), texturesInfo, texCache, textures);
        mat.row0 [2] = sampler.row0;
        mat.row1 [2] = sampler.row1;
        mat.texid[2] = texID;
        mat.cflags |= (FLAG_FOUR_TEXTURES | FLAG_PACK_FOUR_PARAMS_IN_TEXTURE);
      }
    }
  }

  mat.colors[GLTF_COLOR_BASE]  = baseColor; 
  mat.colors[GLTF_COLOR_METAL] = float4(1.0f); 
  mat.colors[GLTF_COLOR_COAT]  = float4(1.0f); 
  mat.data  [GLTF_FLOAT_ALPHA] = metalness;
  mat.data  [GLTF_FLOAT_GLOSINESS] = reflGlossiness;
  mat.data  [GLTF_FLOAT_IOR]       = fresnelIOR;
  SetMiPlastic(&mat, fresnelIOR, 1.0f, baseColor, float4(1,1,1,1));

  return mat;
}

Material ConvertOldHydraMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                                 std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                 std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                                 bool is_spectral_mode)
{
  std::wstring name              = materialNode.attribute(L"name").as_string();
  Material mat                   = {};
  mat.mtype                      = MAT_TYPE_GLTF;
  mat.data[GLTF_FLOAT_ALPHA]     = 0.0f;
  mat.data[GLTF_FLOAT_REFL_COAT] = 1.0f;
  mat.colors[GLTF_COLOR_COAT]    = float4(1,1,1,1); 
  mat.colors[GLTF_COLOR_METAL]   = float4(0,0,0,0);  
  mat.lightId                    = uint(-1);
  
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
      mat.data[GLTF_FLOAT_ALPHA]     = 0.0f;
      mat.data[GLTF_FLOAT_REFL_COAT] = 1.0f;
      mat.colors[GLTF_COLOR_COAT]  = reflColor;
      mat.colors[GLTF_COLOR_METAL] = float4(0,0,0,0); 
      mat.cflags                   = GLTF_COMPONENT_LAMBERT | GLTF_COMPONENT_COAT;
      SetMiPlastic(&mat, fresnelIOR, 1.0f, color, reflColor);
    }
    else
    {
      mat.data[GLTF_FLOAT_ALPHA]     = length(reflColor)/( length(reflColor) + length3f(color) );
      mat.data[GLTF_FLOAT_REFL_COAT] = 0.0f;
      mat.colors[GLTF_COLOR_COAT]  = float4(0,0,0,0); 
      mat.colors[GLTF_COLOR_METAL] = reflColor;   // disable coating for such blend type
      mat.cflags                   = GLTF_COMPONENT_LAMBERT | GLTF_COMPONENT_METAL;
    }
  }
  else if(length(reflColor) > 1e-5f)
  {
    mat.mtype  = MAT_TYPE_GLTF;
    mat.cflags = GLTF_COMPONENT_METAL;
    mat.colors[GLTF_COLOR_BASE]    = reflColor;
    mat.colors[GLTF_COLOR_METAL]   = float4(1.0f);
    mat.colors[GLTF_COLOR_COAT]    = float4(0.0f); 
    mat.data[GLTF_FLOAT_ALPHA]     = 1.0f;
  }
  else if(length(to_float3(color)) > 1e-5f)
  {
    mat.mtype  = MAT_TYPE_GLTF;
    mat.cflags = GLTF_COMPONENT_LAMBERT;
    mat.colors[GLTF_COLOR_BASE]    = color;
    mat.colors[GLTF_COLOR_COAT]    = float4(0.0f); 
    mat.colors[GLTF_COLOR_METAL]   = float4(0.0f);    
    mat.data[GLTF_FLOAT_ALPHA]     = 0.0f;
    mat.data[GLTF_FLOAT_REFL_COAT] = 0.0f;
  }
    
  // Glass
  if (length(transpColor) > 1e-5f)
  {
    mat.mtype                           = MAT_TYPE_GLASS;
    mat.colors[GLTF_COLOR_BASE]         = reflColor;   
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


Material LoadDiffuseMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                             std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                             std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                             std::vector<uint2> &spec_tex_ids_wavelengths,
                             const std::vector<uint2> &spec_tex_offset_sz, std::set<uint32_t> &loadedSpectralTextures,
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

    if(is_spectral_mode && specId != 0xFFFFFFFF)
    {
      LoadSpectralTextures(specId, texturesInfo, texCache, textures, spec_tex_ids_wavelengths, spec_tex_offset_sz, 
                           loadedSpectralTextures);
    }
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
                             const std::vector<uint2> &spec_offsets, std::vector<uint2> &spec_tex_ids_wavelengths,
                             const std::vector<uint2> &spec_tex_offset_sz, std::set<uint32_t> &loadedSpectralTextures)
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

    if(is_spectral_mode && specId != 0xFFFFFFFF)
    {
      LoadSpectralTextures(specId, texturesInfo, texCache, textures, spec_tex_ids_wavelengths, spec_tex_offset_sz, 
                           loadedSpectralTextures);
    }
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
    if(offsets.x != 0xFFFFFFFF)
    {
      spectrum.reserve(offsets.y);
      
      std::copy(spectra.begin() + offsets.x, spectra.begin() + offsets.x + offsets.y, std::back_inserter(spectrum));
    }
  }

  float4 diffuse_reflectance = mat.colors[PLASTIC_COLOR];

  // if(!is_spectral_mode)
  // {
  //   uint32_t colorTexId = mat.texid[0];
  //   if(colorTexId > 0 && colorTexId != 0xFFFFFFFF)
  //     diffuse_reflectance *= image2D_average(textures[colorTexId]);
  // }

  auto precomp = mi::fresnel_coat_precompute(mat.data[PLASTIC_ROUGHNESS], internal_ior, external_ior, diffuse_reflectance,
                                            {1.0f, 1.0f, 1.0f, 1.0f}, is_spectral_mode, spectrum);

  mat.data[PLASTIC_PRECOMP_REFLECTANCE] = precomp.internal_reflectance;
  mat.data[PLASTIC_SPEC_SAMPLE_WEIGHT]  = precomp.specular_sampling_weight;

  std::copy(precomp.transmittance.begin(), precomp.transmittance.end(), std::back_inserter(precomputed_transmittance));
  
  mat.datai[0] = (precomputed_transmittance.size() / MI_ROUGH_TRANSMITTANCE_RES) - 1u;

  return mat;
}
