#include "integrator_pt.h"
#include "include/crandom.h"

#include "LiteScene/cmesh4.h"
using cmesh4::SimpleMesh;

#include "mi_materials.h"

//#define LAYOUT_STD140 // !!! PLEASE BE CAREFUL WITH THIS !!!
#include "LiteScene/hydraxml.h"

#include <string>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

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



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

GLTFMaterial ConvertOldHydraMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                                     std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                     std::vector< std::shared_ptr<ICombinedImageSampler> > &textures)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  GLTFMaterial mat = {};
  mat.mtype        = MAT_TYPE_GLTF;
  mat.alpha        = 0.0f;
  mat.coatColor    = float4(1,1,1,0); 
  mat.metalColor   = float4(0,0,0,0);  
  mat.lightId      = uint(-1);
  
  auto nodeEmiss = materialNode.child(L"emission");

  // read Hydra or GLTF materials
  //
  float4 color(0.0f, 0.0f, 0.0f, 0.0f);

  if(materialNode.attribute(L"light_id") != nullptr || nodeEmiss != nullptr)
  {
    auto nodeEmissColor = nodeEmiss.child(L"color");
    color   = to_float4(hydra_xml::readval3f(nodeEmissColor), 1.0f);

    HydraSampler emissiveSampler = ReadSamplerFromColorNode(nodeEmissColor);
    auto p = texCache.find(emissiveSampler);
    if(p == texCache.end())
    {
      texCache[emissiveSampler] = uint(textures.size());
      const uint32_t texId = nodeEmissColor.child(L"texture").attribute(L"id").as_uint();
      textures.push_back(LoadTextureAndMakeCombined(texturesInfo[texId], emissiveSampler.sampler));
      p = texCache.find(emissiveSampler);
    }
    
    mat.row0 [0]  = emissiveSampler.row0;
    mat.row1 [0]  = emissiveSampler.row1;
    mat.texId[0]  = p->second;
    
    mat.baseColor = color;
    if(materialNode.attribute(L"light_id") == nullptr)
      mat.lightId = uint(-1);
    else
      mat.lightId = uint(materialNode.attribute(L"light_id").as_int()); // for correct process of "-1"
    
    mat.mtype = MAT_TYPE_LIGHT_SOURCE;
  }

  auto nodeDiffColor = materialNode.child(L"diffuse").child(L"color");
  if(nodeDiffColor != nullptr)
  {
    color = to_float4(hydra_xml::readval3f(nodeDiffColor), 0.0f);
    HydraSampler diffSampler = ReadSamplerFromColorNode(nodeDiffColor);
    auto p = texCache.find(diffSampler);
    if(p == texCache.end())
    {
      texCache[diffSampler] = uint(textures.size());
      const uint32_t texId  = nodeDiffColor.child(L"texture").attribute(L"id").as_uint();
      textures.push_back(LoadTextureAndMakeCombined(texturesInfo[texId], diffSampler.sampler));
      p = texCache.find(diffSampler);
    }
    
    mat.row0 [0]  = diffSampler.row0;
    mat.row1 [0]  = diffSampler.row1;
    mat.texId[0]  = p->second;
  }

  float3 reflColor = float3(0,0,0);
  float glosiness  = 1.0f;
  float fresnelIOR = 1.5f;
  auto nodeRefl = materialNode.child(L"reflectivity");
  if(nodeRefl != nullptr)
  {
    reflColor  = hydra_xml::readval3f(nodeRefl.child(L"color"));
    glosiness  = hydra_xml::readval1f(nodeRefl.child(L"glossiness"));  
    fresnelIOR = hydra_xml::readval1f(nodeRefl.child(L"fresnel_ior"));
  }

  const bool hasFresnel  = (nodeRefl.child(L"fresnel").attribute(L"val").as_int() != 0);
  if(!hasFresnel)
    fresnelIOR = 0.0f;
  
  if((length(reflColor) > 1e-5f && length(to_float3(color)) > 1e-5f) || hasFresnel)
  {
    mat.mtype      = MAT_TYPE_GLTF;
    mat.baseColor  = color;
    mat.coatColor  = to_float4(reflColor, 0.0f);
    mat.lightId    = uint(-1);

    if(hasFresnel)
    {
      mat.alpha      = 0.0f;
      mat.coatColor  = to_float4(reflColor, 0.0f);
      mat.metalColor = float4(0,0,0,0); 
      mat.cflags     = GLTF_COMPONENT_LAMBERT | GLTF_COMPONENT_COAT;
    }
    else
    {
      mat.alpha      = length(reflColor)/( length(reflColor) + length3f(color) );
      mat.coatColor  = float4(0,0,0,0); 
      mat.metalColor = to_float4(reflColor, 0.0f);                                           // disable coating for such blend type
      mat.cflags     = GLTF_COMPONENT_LAMBERT | GLTF_COMPONENT_METAL;
    }

    SetMiPlastic(&mat, fresnelIOR, 1.0f, to_float3(color), reflColor);
  }
  else if(length(reflColor) > 1e-5f)
  {
    mat.mtype      = MAT_TYPE_GLTF;
    mat.cflags     = GLTF_COMPONENT_METAL;
    mat.metalColor = to_float4(reflColor, 0.0f);
    mat.coatColor  = float4(0,0,0,0); 
    mat.alpha      = 1.0f;
  }
  else if(length(to_float3(color)) > 1e-5f)
  {
    mat.mtype      = MAT_TYPE_GLTF;
    mat.cflags     = GLTF_COMPONENT_LAMBERT;
    mat.baseColor  = color;
    mat.coatColor  = float4(0,0,0,0); 
    mat.metalColor = float4(0,0,0,0);    
    mat.alpha      = 0.0f;
  }
  
  if(color[3] > 1e-5f)
    mat.mtype    = MAT_TYPE_LIGHT_SOURCE;

  mat.glosiness  = glosiness;
  mat.ior        = fresnelIOR;
  
  return mat;
}

GLTFMaterial LoadRoughConductorMaterial(const pugi::xml_node& materialNode, const std::vector<TextureInfo> &texturesInfo,
                                        std::unordered_map<HydraSampler, uint32_t, HydraSamplerHash> &texCache, 
                                        std::vector< std::shared_ptr<ICombinedImageSampler> > &textures)
{
  std::wstring name = materialNode.attribute(L"name").as_string();
  GLTFMaterial mat = {};
  mat.mtype        = MAT_TYPE_CONDUCTOR;
  mat.alpha        = 0.0f;
  mat.coatColor    = float4(1,1,1,0); 
  mat.metalColor   = float4(0,0,0,0);  
  mat.lightId      = uint(-1);

  auto nodeBSDF = materialNode.child(L"bsdf");

  auto bsdf_type = nodeBSDF.attribute(L"type").as_string();
  auto alpha_u   = materialNode.child(L"alpha_u").attribute(L"val").as_float();
  auto alpha_v   = materialNode.child(L"alpha_v").attribute(L"val").as_float();
  auto eta       = materialNode.child(L"eta").attribute(L"val").as_float();
  auto k         = materialNode.child(L"k").attribute(L"val").as_float();

  mat.metalColor = {alpha_u, alpha_v, eta, k};

  return mat;
}

bool Integrator::LoadScene(const char* a_scehePath, const char* a_sncDir)
{ 
  std::string scenePathStr(a_scehePath);
  std::string sceneDirStr(a_sncDir);  
  hydra_xml::HydraScene scene;
  
  auto loadRes = scene.LoadState(scenePathStr, sceneDirStr);
  if(loadRes != 0)
  {
    std::cout << "Integrator::LoadScene failed: '" << a_scehePath << "'" << std::endl; 
    exit(0);
  }

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

  //// (1) load materials
  //
  m_materials.resize(0);
  m_materials.reserve(100);

  const std::wstring hydraOldMatTypeStr {L"hydra_material"};
  const std::wstring roughConductorMatTypeStr {L"rough_conductor"};

  for(auto materialNode : scene.MaterialNodes())
  {
    GLTFMaterial mat = {};
    auto mat_type = materialNode.attribute(L"type").as_string();
    if(mat_type == hydraOldMatTypeStr)
    {
      mat = ConvertOldHydraMaterial(materialNode, texturesInfo, texCache, m_textures);
    }
    else if(mat_type == roughConductorMatTypeStr)
    {
      mat = LoadRoughConductorMaterial(materialNode, texturesInfo, texCache, m_textures);
    }
    m_materials.push_back(mat);
  }

  // load first camera and update matrix
  //
  for(auto cam : scene.Cameras())
  {
    float aspect   = 1.0f;
    auto proj      = perspectiveMatrix(cam.fov, aspect, cam.nearPlane, cam.farPlane);
    auto worldView = lookAt(float3(cam.pos), float3(cam.lookAt), float3(cam.up));
      
    m_projInv      = inverse4x4(proj);
    m_worldViewInv = inverse4x4(worldView);
    break; // take first cam
  }

  // load lights
  //
  m_instIdToLightInstId.resize(scene.GetInstancesNum(), -1);

  for(auto lightInst : scene.InstancesLights())
  {
    const std::wstring ltype = lightInst.lightNode.attribute(L"type").as_string();
    const std::wstring shape = lightInst.lightNode.attribute(L"shape").as_string();
    const float sizeX = lightInst.lightNode.child(L"size").attribute(L"half_width").as_float();
    const float sizeZ = lightInst.lightNode.child(L"size").attribute(L"half_length").as_float();
    float power = lightInst.lightNode.child(L"intensity").child(L"multiplier").text().as_float();
    if(power == 0.0f)
      power = lightInst.lightNode.child(L"intensity").child(L"multiplier").attribute(L"val").as_float();

    float3 color = hydra_xml::readval3f(lightInst.lightNode.child(L"intensity").child(L"color"));
    auto matrix  = lightInst.matrix;

    if(ltype == std::wstring(L"sky"))
    {
      m_envColor = to_float4(color * power, 1.0f); // set pdf to 1.0f
    }
    else if(ltype == std::wstring(L"directional"))
    {
      LightSource lightSource{};
      lightSource.pos       = lightInst.matrix * float4(0.0f, 0.0f, 0.0f, 1.0f);
      lightSource.norm      = normalize(lightInst.matrix * float4(0.0f, -1.0f, 0.0f, 0.0f));
      lightSource.intensity = to_float4(color*power,0);
      lightSource.geomType  = LIGHT_GEOM_DIRECT;
      m_lights.push_back(lightSource);
    }
    else if(shape == L"rect" || shape == L"disk")
    {
      LightSource lightSource{};

      lightSource.pos       = lightInst.matrix * float4(0.0f, 0.0f, 0.0f, 1.0f);
      lightSource.norm      = normalize(lightInst.matrix * float4(0.0f, -1.0f, 0.0f, 0.0f));
      lightSource.intensity = to_float4(color*power,0);
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

      m_lights.push_back(lightSource);
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

      LightSource lightSource{};
      {
        lightSource.pos       = lightInst.matrix * float4(0.0f, 0.0f, 0.0f, 1.0f);
        lightSource.norm      = float4(0.0f, -1.0f, 0.0f, 0.0f);
        lightSource.intensity = to_float4(color*power,0);
        lightSource.geomType  = LIGHT_GEOM_SPHERE;
  
        lightSource.matrix    = float4x4{};
        lightSource.size      = float2(radius, radius);
        lightSource.pdfA      = 1.0f / (4.0f*LiteMath::M_PI*radius*radius);
      }
      m_lights.push_back(lightSource);
    }
  }

  //// (2) load meshes
  //
  m_matIdOffsets.reserve(1024);
  m_vertOffset.reserve(1024);
  m_matIdByPrimId.reserve(128000);
  m_triIndices.reserve(128000*3);

  m_vNorm4f.resize(0);
  m_vTexc2f.resize(0);

  m_pAccelStruct->ClearGeom();
  for(auto meshPath : scene.MeshFiles())
  {
    std::cout << "[LoadScene]: mesh = " << meshPath.c_str() << std::endl;
    auto currMesh = cmesh4::LoadMeshFromVSGF(meshPath.c_str());
    auto geomId   = m_pAccelStruct->AddGeom_Triangles3f((const float*)currMesh.vPos4f.data(), currMesh.vPos4f.size(), currMesh.indices.data(), currMesh.indices.size(), BUILD_HIGH, sizeof(float)*4);

    (void)geomId; // silence unused var. warning

    m_matIdOffsets.push_back(static_cast<unsigned int>(m_matIdByPrimId.size()));
    m_vertOffset.push_back(static_cast<unsigned int>(m_vNorm4f.size()));

    m_matIdByPrimId.insert(m_matIdByPrimId.end(), currMesh.matIndices.begin(), currMesh.matIndices.end() );
    m_triIndices.insert(m_triIndices.end(), currMesh.indices.begin(), currMesh.indices.end());

    //for(size_t i=0;i<currMesh.vPos4f.size();i++) // pack texture coords to fourth components
    //{
    //  currMesh.vPos4f [i].w = currMesh.vTexCoord2f[i].x;
    //  currMesh.vNorm4f[i].w = currMesh.vTexCoord2f[i].y;
    //}

    m_vNorm4f.insert(m_vNorm4f.end(), currMesh.vNorm4f.begin(),     currMesh.vNorm4f.end());     // #TODO: store compressed normal and tangent together
    m_vTexc2f.insert(m_vTexc2f.end(), currMesh.vTexCoord2f.begin(), currMesh.vTexCoord2f.end()); // #TODO: store quantized texture coordinates
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
  for(auto sett : scene.Settings())
  {
    m_traceDepth = sett.depth;
    m_spp        = sett.spp;

    if(m_traceDepth == 0)
      m_traceDepth = 6;
    
    if(m_spp == 0)
      m_spp = 1;

    break; // take first render settings
  }

  return true;
}