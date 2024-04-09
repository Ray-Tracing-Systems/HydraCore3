#include "integrator_pt_scene.h"
#include <unordered_map>
#include <spectral/spec/basic_spectrum.h>
#include <spectral/spec/conversions.h>
#include "spectrum_loader.h"
#include <algorithm>

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

HydraSampler ReadSamplerFromColorNode(const pugi::xml_node a_colorNodes, bool from_spectrum)
{
  HydraSampler res;
  pugi::xml_node texNode;

  if(from_spectrum)
    texNode = a_colorNodes.child(L"spectrum");
  else
    texNode = a_colorNodes.child(L"texture");
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

std::optional<Spectrum> LoadSpectrumFromNode(const pugi::xml_node& a_node, const std::vector<SpectrumLoader> &spectraInfo)
{
  auto specNode = a_node.child(L"spectrum");
  if(specNode)
  {
    uint32_t spec_id = specNode.attribute(L"id").as_uint();
    if(spec_id >= spectraInfo.size()) {
      std::cerr << "Unknown spectum error" << std::endl;
      return {};
    } 

    const auto &ref = spectraInfo[spec_id];
    return ref.load();
  }

  return {};
}

uint32_t GetSpectrumIdFromNode(const pugi::xml_node& a_node)
{
  uint32_t spec_id = INVALID_SPECTRUM_ID;
  auto specNode = a_node.child(L"spectrum");
  if(specNode != nullptr)
  {
    spec_id = specNode.attribute(L"id").as_uint();
  }

  return spec_id;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::optional<float4> GetColorFromNode(const pugi::xml_node& a_node, ResourceContext &resources)
{
  auto val = hydra_xml::readvalVariant(a_node);

  auto node = a_node.child(L"spectrum");
  if(node) 
  {
    const auto sp_struct = LoadSpectrumFromNode(a_node, resources.spectraInfo);
    if(!sp_struct) 
    {
      std::cerr << "Spectral value not found (LoadSpectrumFromNode)" << std::endl;
      return {};
    }
    return DownsampleSpectrum(*sp_struct);
  }

  if(std::holds_alternative<float>(val))
  {
    return float4(std::get<float>(val));
  }
  else if(std::holds_alternative<float3>(val))
  {
    float3 rgb = std::get<float3>(val);
    return to_float4(rgb, 0.0f);
  }
  return {};
}

ColorHolder GetVariableColorFromNode(const pugi::xml_node& a_node, ResourceContext &resources, bool is_spectral_mode)
{ 
  if(is_spectral_mode) {
    uint32_t spec_id = GetSpectrumIdFromNode(a_node);
    if(spec_id != INVALID_SPECTRUM_ID) {
      return {spec_id};
    }
    else if(!a_node.child(L"texture")) {//Register new spectrum
      auto coloropt = GetColorFromNode(a_node, resources);
      if(coloropt) {
        return UpsampleSpectrumFromColor(*coloropt, resources.spectraInfo, resources.loadedSpectrumCount);
      }
    }
    return {};
  }
  else {
    auto coloropt = GetColorFromNode(a_node, resources);
    if(coloropt) {
      return {*coloropt};
    }
    return {};
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void LoadSpectralTextures(const uint32_t specId, HydraSampler& a_sampler,
                          const ResourceContext &resources,
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
      // HydraSampler sampler;
      // sampler.inputGamma = 1.0f;
      a_sampler.texId = xml_tex_id;
      // sampler.sampler.addressU = Sampler::AddressMode::CLAMP;
      // sampler.sampler.addressV = Sampler::AddressMode::CLAMP;

      const auto& [sampler_out, loaded_tex_id] = LoadTextureById(xml_tex_id, resources, a_sampler, texCache, textures);
      spec_tex_ids_wavelengths[offset + i].x = loaded_tex_id;     
    }
    loadedSpectralTextures.insert(specId);
  }
}


Material ConvertGLTFMaterial(const pugi::xml_node& materialNode, ResourceContext &resources,
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
      auto baseColorOpt = GetColorFromNode(materialNode.child(L"color"), resources);
      if(baseColorOpt) {
        baseColor = *baseColorOpt;
      }
      if(materialNode.child(L"color").child(L"texture") != nullptr) {
        const auto& [sampler, texID] = LoadTextureFromNode(materialNode.child(L"color"), resources, texCache, textures);
        mat.row0 [0] = sampler.row0;
        mat.row1 [0] = sampler.row1;
        mat.texid[0] = texID;
      }
    }

    if(materialNode.child(L"glossiness") != nullptr)
    {
      reflGlossiness = hydra_xml::readval1f(materialNode.child(L"glossiness"));  
      if(materialNode.child(L"glossiness").child(L"texture") != nullptr) {
        const auto& [sampler, texID] = LoadTextureFromNode(materialNode.child(L"glossiness"), resources, texCache, textures);
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
        const auto& [sampler, texID] = LoadTextureFromNode(materialNode.child(L"roughness"), resources, texCache, textures);
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
        const auto& [sampler, texID] = LoadTextureFromNode(materialNode.child(L"metalness"), resources, texCache, textures);
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
        const auto& [sampler, texID] = LoadTextureFromNode(materialNode.child(L"glossiness_metalness_coat"), resources, texCache, textures);
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

Material ConvertOldHydraMaterial(const pugi::xml_node& materialNode, ResourceContext &resources,
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
  mat.spdid[0]                   = INVALID_SPECTRUM_ID;
  auto nodeEmiss = materialNode.child(L"emission");

  // read Hydra or GLTF materials
  //
  float4 color(0.0f, 0.0f, 0.0f, 0.0f);

  bool is_emission_color = false;
  if(materialNode.attribute(L"light_id") != nullptr || nodeEmiss != nullptr)
  {
    auto nodeEmissColor = nodeEmiss.child(L"color");

    /*auto var_color = GetVariableColorFromNode(nodeEmissColor, resources, is_spectral_mode);
    if(var_color.isRGB()) {
        color = var_color.getRGB();
    }
    else if(var_color) {
      mat.spdid[0] = var_color.getSpectrumId();
      std::cerr << "Mat specid = " << mat.spdid[0] << std::endl; 
    }
    */

    color = GetColorFromNode(nodeEmissColor, resources).value_or(float4(1.0f));
    mat.spdid[0] = GetSpectrumIdFromNode(nodeEmissColor);
    is_emission_color = (materialNode.attribute(L"light_id") != nullptr) || (length(color) > 1e-5f );

    const auto& [emissiveSampler, texID] = LoadTextureFromNode(nodeEmissColor, resources, texCache, textures);
    
    mat.row0 [0] = emissiveSampler.row0;
    mat.row1 [0] = emissiveSampler.row1;
    mat.texid[0] = texID;
    
    mat.colors[EMISSION_COLOR] = color;
    if(materialNode.attribute(L"light_id") == nullptr)
      mat.lightId = uint(-1);
    else
      mat.lightId = uint(materialNode.attribute(L"light_id").as_int());  // for correct process of "-1"

    //auto specId  = GetSpectrumIdFromNode(nodeEmissColor);  
   // mat.spdid[0] = specId;
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
    color = GetColorFromNode(nodeDiffColor, resources).value_or(float4(1.0f));
    const auto& [diffSampler, texID] = LoadTextureFromNode(nodeDiffColor, resources, texCache, textures);
    
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
    reflColor       = GetColorFromNode(nodeRefl.child(L"color"), resources).value_or(float4(0.0f));
    reflGlossiness  = hydra_xml::readval1f(nodeRefl.child(L"glossiness"));  
    fresnelIOR      = hydra_xml::readval1f(nodeRefl.child(L"fresnel_ior"));
  }

  float4 transpColor      = float4(0, 0, 0, 0);
  float  transpGlossiness = 1.0f;
  float  transpIOR        = 1.5f;

  auto nodeTransp = materialNode.child(L"transparency");
  if (nodeTransp != nullptr)
  {
    transpColor      = GetColorFromNode(nodeTransp.child(L"color"), resources).value_or(float4(0.0f));
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

Material LoadRoughConductorMaterial(const pugi::xml_node& materialNode, ResourceContext &resources,
                                    std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                    std::vector< std::shared_ptr<ICombinedImageSampler> > &textures, bool is_spectral_mode)
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

    const auto& [sampler, texID] = LoadTextureFromNode(nodeAlpha, resources, texCache, textures);

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

  float eta       = materialNode.child(L"eta").attribute(L"val").as_float();
  uint32_t etaSpecId = GetSpectrumIdFromNode(materialNode.child(L"eta"));
  float k         = materialNode.child(L"k").attribute(L"val").as_float();
  uint32_t kSpecId   = GetSpectrumIdFromNode(materialNode.child(L"k"));

  auto nodeColor = materialNode.child(L"reflectance");
  if(nodeColor)
  {
    if(!is_spectral_mode) {
      mat.colors[CONDUCTOR_COLOR] = GetColorFromNode(nodeColor, resources).value_or(float4(1.0f));
    }
    else if(etaSpecId == INVALID_SPECTRUM_ID && kSpecId == INVALID_SPECTRUM_ID) {
      auto coloropt = GetColorFromNode(nodeColor, resources);
      if(coloropt) {

        /* PREMP
        float etaMp = eta > k ? 1.0f : eta / k;
        float kMp = eta > k ? k / eta : 1.0f;


        spec::ISpectrum::ptr etaupsampled = UpsampleRaw(spec::vec3(coloropt->x, coloropt->y, coloropt->z) * etaMp);
        spec::ISpectrum::ptr kupsampled = UpsampleRaw(spec::vec3(coloropt->x, coloropt->y, coloropt->z) * kMp);

        size_t count = size_t(LAMBDA_MAX - LAMBDA_MIN);
        spec::BasicSpectrum *etaspec = new spec::BasicSpectrum();
        spec::BasicSpectrum *kspec = new spec::BasicSpectrum();
        for(unsigned c = 0; c < count; c++) {
          spec::Float lambda = (LAMBDA_MIN + spec::Float(c));
          auto p = spec::color2ior(etaupsampled->get_or_interpolate(lambda), kupsampled->get_or_interpolate(lambda));
          etaspec->set(lambda, p.first);
          kspec->set(lambda, p.second);
        }
        etaSpecId = uint32_t(resources.spectraInfo.size());
        resources.spectraInfo.push_back({spec::ISpectrum::ptr(etaspec), etaSpecId});
        kSpecId = uint32_t(resources.spectraInfo.size());
        resources.spectraInfo.push_back({spec::ISpectrum::ptr(kspec), kSpecId});
        */
        
        spec::ISpectrum::ptr upsampled = UpsampleRaw({coloropt->x, coloropt->y, coloropt->z});
        size_t count = size_t(LAMBDA_MAX - LAMBDA_MIN);
        spec::BasicSpectrum *etaspec = new spec::BasicSpectrum();
        spec::BasicSpectrum *kspec = new spec::BasicSpectrum();
        for(unsigned c = 0; c < count; c++) {
          spec::Float lambda = (LAMBDA_MIN + spec::Float(c));
          auto p = spec::color2ior(upsampled->get_or_interpolate(lambda));
          etaspec->set(lambda, p.first);
          kspec->set(lambda, p.second);
        }
        etaSpecId = uint32_t(resources.spectraInfo.size());
        resources.spectraInfo.push_back({spec::ISpectrum::ptr(etaspec), etaSpecId});
        kSpecId = uint32_t(resources.spectraInfo.size());
        resources.spectraInfo.push_back({spec::ISpectrum::ptr(kspec), kSpecId});

        /* POSTMP
        size_t count = size_t(LAMBDA_MAX - LAMBDA_MIN);
        spec::BasicSpectrum *etaspec = new spec::BasicSpectrum();
        spec::BasicSpectrum *kspec = new spec::BasicSpectrum();
        for(unsigned c = 0; c < count; c++) {
          spec::Float lambda = (LAMBDA_MIN + spec::Float(c));
          auto p = spec::color2ior(upsampled->get_or_interpolate(lambda) * eta, upsampled->get_or_interpolate(lambda) * k);
          etaspec->set(lambda, p.first);
          kspec->set(lambda, p.second);
        }
        etaSpecId = uint32_t(resources.spectraInfo.size());
        resources.spectraInfo.push_back({spec::ISpectrum::ptr(etaspec), etaSpecId});
        kSpecId = uint32_t(resources.spectraInfo.size());
        resources.spectraInfo.push_back({spec::ISpectrum::ptr(kspec), kSpecId});
        */
      }
    }
  }

  mat.data[CONDUCTOR_ROUGH_U] = alpha_u;
  mat.data[CONDUCTOR_ROUGH_V] = alpha_v; 
  mat.data[CONDUCTOR_ETA]     = eta; 
  mat.data[CONDUCTOR_K]       = k;   
  
  mat.spdid[0] = etaSpecId;
  mat.spdid[1] = kSpecId;

  return mat;
}


Material LoadDiffuseMaterial(const pugi::xml_node& materialNode, ResourceContext &resources,
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
  mat.spdid[0] = INVALID_SPECTRUM_ID;
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
    ColorHolder var_color = GetVariableColorFromNode(nodeColor, resources, is_spectral_mode);
    if(var_color.isRGB()) {
      mat.colors[DIFFUSE_COLOR] = var_color.getRGB();
    }
    else if(var_color){
      mat.spdid[0] = var_color.getSpectrumId();
    }

    const auto& [sampler, texID] = LoadTextureFromNode(nodeColor, resources, texCache, textures);
    
    mat.row0 [0]  = sampler.row0;
    mat.row1 [0]  = sampler.row1;
    mat.texid[0]  = texID;

    if(is_spectral_mode && resources.specTexIds.find(mat.spdid[0]) != resources.specTexIds.end())
    {
      auto sampler = ReadSamplerFromColorNode(nodeColor, true);
      LoadSpectralTextures(mat.spdid[0], sampler, resources, texCache, textures, spec_tex_ids_wavelengths, spec_tex_offset_sz, 
                           loadedSpectralTextures);
    }
  }

  return mat;
}


Material LoadDielectricMaterial(const pugi::xml_node& materialNode, ResourceContext &resources,
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
    mat.colors[DIELECTRIC_COLOR_REFLECT] = GetColorFromNode(nodeReflColor, resources).value_or(float4(1.0f));
  }

  auto nodeTransColor = materialNode.child(L"transmittance");
  if(nodeTransColor != nullptr)
  {
    mat.colors[DIELECTRIC_COLOR_TRANSMIT] = GetColorFromNode(nodeTransColor, resources).value_or(float4(1.0f));
  }

  return mat;
}


Material LoadBlendMaterial(const pugi::xml_node& materialNode, ResourceContext &resources,
                           std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                           std::vector< std::shared_ptr<ICombinedImageSampler> > &textures)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  Material mat = {};

  mat.mtype    = MAT_TYPE_BLEND;
  mat.cflags   = 0;
  mat.texid[0] = 0;
  mat.data[BLEND_WEIGHT] = 1.0f;
  
  mat.datai[0] = materialNode.child(L"bsdf_1").attribute(L"id").as_uint();
  mat.datai[1] = materialNode.child(L"bsdf_2").attribute(L"id").as_uint();

  auto nodeWeight = materialNode.child(L"weight");
  if(nodeWeight != nullptr)
  {
    mat.data[BLEND_WEIGHT] = hydra_xml::readval1f(nodeWeight);

    const auto& [sampler, texID] = LoadTextureFromNode(nodeWeight, resources, texCache, textures);
    
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
    for(size_t j = 0; j < channels; ++j)
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

Material LoadPlasticMaterial(const pugi::xml_node& materialNode, ResourceContext &resources,
                             std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache,
                             std::vector< std::shared_ptr<ICombinedImageSampler> > &textures,
                             std::vector<float> &precomputed_transmittance,
                             bool is_spectral_mode,
                             std::vector<uint2> &spec_tex_ids_wavelengths,
                             const std::vector<uint2> &spec_tex_offset_sz, std::set<uint32_t> &loadedSpectralTextures)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  Material mat = {};
  
  mat.mtype     = MAT_TYPE_PLASTIC;
  mat.lightId   = uint(-1);
  mat.nonlinear = 0;
  mat.texid[0]  = 0;
  //std::fill(mat.spdid, mat.spdid + 4, INVALID_SPECTRUM_ID);
  mat.spdid[0]  = INVALID_SPECTRUM_ID;

  auto nodeColor = materialNode.child(L"reflectance");
  uint32_t specId = INVALID_SPECTRUM_ID;
  if(nodeColor)
  {
    ColorHolder var_color = GetVariableColorFromNode(nodeColor, resources, is_spectral_mode);

    if(var_color.isRGB()) {
      mat.colors[PLASTIC_COLOR] = var_color.getRGB();
    }
    else if(var_color) {
      specId = var_color.getSpectrumId();
    }
    const auto& [sampler, texID] = LoadTextureFromNode(nodeColor, resources, texCache, textures);

    mat.row0 [0]  = sampler.row0;
    mat.row1 [0]  = sampler.row1;
    mat.texid[0]  = texID;

  //  specId = GetSpectrumIdFromNode(nodeColor);
    mat.spdid[0] = specId;

    if(is_spectral_mode && resources.specTexIds.find(specId) != resources.specTexIds.end())
    {
      auto sampler = ReadSamplerFromColorNode(nodeColor, true);
      LoadSpectralTextures(specId, sampler, resources, texCache, textures, spec_tex_ids_wavelengths, spec_tex_offset_sz, 
                           loadedSpectralTextures);
    }

  }

  float internal_ior = hydra_xml::readval1f(materialNode.child(L"int_ior"), 1.49f);
  float external_ior = hydra_xml::readval1f(materialNode.child(L"ext_ior"), 1.000277f);

  mat.data[PLASTIC_IOR_RATIO] = internal_ior / external_ior;

  mat.data[PLASTIC_ROUGHNESS] = hydra_xml::readval1f(materialNode.child(L"alpha"), 0.1f);

  // dirty hack 
  if(mat.data[PLASTIC_ROUGHNESS] == 0.0f)
  {
    mat.data[PLASTIC_ROUGHNESS] = 1e-6f;
  }

  mat.nonlinear = hydra_xml::readval1u(materialNode.child(L"nonlinear"), 0);

  std::vector<float> spectrum;

  if(is_spectral_mode && specId != INVALID_SPECTRUM_ID)
  {
    const auto &spec = resources.spectraInfo[specId].load();
    if(spec) //Might need to check
    {
      spectrum = spec->ResampleUniform();
    }
    else {
      std::cerr << "Using dummy spectrum" << std::endl;
      Spectrum uniform1;
      uniform1.id = 0;

      spec::BasicSpectrum *sp = new spec::BasicSpectrum{
        {200.0f, 1.0f},
        {400.0f, 1.0f},
        {600.0f, 1.0f},
        {800.0f, 1.0f}
      };

      uniform1.spectrum.reset(sp);
      auto specValsUniform = uniform1.ResampleUniform();
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
