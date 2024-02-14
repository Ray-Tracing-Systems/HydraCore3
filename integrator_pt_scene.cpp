#include "integrator_pt_scene.h"

std::string Integrator::GetFeatureName(uint32_t a_featureId)
{
  switch(a_featureId)
  {
    case KSPEC_MAT_TYPE_GLTF      : return "GLTF_LITE";
    case KSPEC_MAT_TYPE_GLASS     : return "GLASS";
    case KSPEC_MAT_TYPE_CONDUCTOR : return "CONDUCTOR";
    case KSPEC_MAT_TYPE_DIFFUSE   : return "DIFFUSE";
    case KSPEC_MAT_TYPE_PLASTIC   : return "PLASTIC";
    case KSPEC_SPECTRAL_RENDERING : return "SPECTRAL";
    case KSPEC_MAT_TYPE_BLEND     : return "BLEND";
    case KSPEC_BUMP_MAPPING       : return "BUMP";
    case KSPEC_MAT_FOUR_TEXTURES  : return "4TEX";
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
static const std::wstring hydraGLTFTypeStr         {L"gltf"};
static const std::wstring roughConductorMatTypeStr {L"rough_conductor"};
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
    else if(mat_type == hydraGLTFTypeStr)
    {
      features[KSPEC_MAT_TYPE_GLTF] = 1;

      auto gloss = materialNode.child(L"glossiness");
      auto rough = materialNode.child(L"roughness");
      auto metal = materialNode.child(L"metalness");
      auto alltm = materialNode.child(L"glossiness_metalness_coat");
      auto nodes = {gloss, rough, metal, alltm};
      for(auto node : nodes)
        if(node.child(L"texture") != nullptr)
          features[KSPEC_MAT_FOUR_TEXTURES] = 1;
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

  for(auto lightInst : g_lastScene.InstancesLights())
  {
    if(lightInst.lightNode.child(L"ies") != nullptr)
      features[KSPEC_LIGHT_IES] = 1;
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
//bool SaveImage4fToEXR(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, bool a_invertY = false);

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

    if (texNode.attribute(L"loc").empty())
      tex.path = std::wstring(sceneFolder.begin(), sceneFolder.end()) + L"/" + texNode.attribute(L"path").as_string();
    else
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
      auto specValsUniform = spec.ResampleUniform();
      
      uint32_t offset = uint32_t(m_spec_values.size());
      std::copy(specValsUniform.begin(),    specValsUniform.end(),    std::back_inserter(m_spec_values));
      m_spec_offset_sz.push_back(uint2{offset, uint32_t(specValsUniform.size())});

    }

    // if no spectra are loaded add uniform 1.0 spectrum
    if(spectraInfo.empty())
    {
      Spectrum uniform1;
      uniform1.id = 0;
      uniform1.wavelengths = {200.0f, 400.0f, 600.0f, 800.0f};
      uniform1.values = {1.0f, 1.0f, 1.0f, 1.0f};
      auto specValsUniform = uniform1.ResampleUniform();
      
      uint32_t offset = uint32_t(m_spec_values.size());
      std::copy(specValsUniform.begin(),    specValsUniform.end(),    std::back_inserter(m_spec_values));
      m_spec_offset_sz.push_back(uint2{offset, uint32_t(specValsUniform.size())});
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
    lightSource.flags    = 0;
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

      lightSource.matrix = matrix;
      lightSource.matrix.set_col(3, float4(0,0,0,1));
      lightSource.size = float2(sizeZ, sizeX);       ///<! Please note tha we HAVE MISTAKEN with ZX order in Hydra2 implementation

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
      sampler.addressU = Sampler::AddressMode::CLAMP;
      sampler.addressV = Sampler::AddressMode::CLAMP;

      m_textures.push_back(MakeCombinedTexture2D(pTexture, sampler));
      lightSource.iesId = uint(m_textures.size()-1);
      
      auto matrixAttrib = iesNode.attribute(L"matrix");
      if(matrixAttrib != nullptr)
      {
        float4x4 mrot           = LiteMath::rotate4x4Y(DEG_TO_RAD*90.0f);
        float4x4 matrixFromNode = hydra_xml::float4x4FromString(matrixAttrib.as_string());
        float4x4 instMatrix     = matrix;
        instMatrix.set_col(3, float4(0,0,0,1));
        lightSource.iesMatrix = mrot*transpose(transpose(matrixFromNode)*instMatrix);
        lightSource.iesMatrix.set_col(3, float4(0,0,0,1));
      }
      
      int pointArea = iesNode.attribute(L"point_area").as_int();
      if(pointArea != 0)
        lightSource.flags |= LIGHT_FLAG_POINT_AREA;

      m_actualFeatures[KSPEC_LIGHT_IES] = 1;
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
    else if(mat_type == hydraGLTFTypeStr)
    {
      mat = ConvertGLTFMaterial(materialNode, texturesInfo, texCache, m_textures, m_spectral_mode);
      m_actualFeatures[KSPEC_MAT_TYPE_GLTF] = 1;
    }
    else if(mat_type == roughConductorMatTypeStr)
    {
      mat = LoadRoughConductorMaterial(materialNode, texturesInfo, texCache, m_textures);
      m_actualFeatures[KSPEC_MAT_TYPE_CONDUCTOR] = 1;
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
                                m_spec_values, m_spec_offset_sz);
      m_actualFeatures[KSPEC_MAT_TYPE_PLASTIC] = 1;
    }
    else if(mat_type == dielectricMatTypeStr)
    {
      mat = LoadDielectricMaterial(materialNode, texturesInfo, texCache, m_textures, m_spectral_mode);
      m_actualFeatures[KSPEC_MAT_TYPE_DIELECTRIC] = 1;
    }

    if((mat.cflags & FLAG_FOUR_TEXTURES) != 0 )
      m_actualFeatures[KSPEC_MAT_FOUR_TEXTURES] = 1;

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
