#include <iostream>
#include <fstream>
#include <filesystem>

#include "integrator_pt.h"
#include "ArgParser.h"

bool SaveImage4fToEXR(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, bool a_invertY = false);
bool SaveImage4fToBMP(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, float a_gamma = 2.2f);
float4x4 ReadMatrixFromString(const std::string& str);

#ifdef USE_VULKAN
#include "vk_context.h"
#include "integrator_pt1_generated.h" // advanced way
//std::shared_ptr<Integrator> CreateIntegrator_Generated(int a_maxThreads, int a_spectral_mode, std::vector<uint32_t> a_features, vk_utils::VulkanContext a_ctx, size_t a_maxThreadsGenerated); // simple way
#endif

int main(int argc, const char** argv)
{
  #ifndef NDEBUG
  bool enableValidationLayers = true;
  #else
  bool enableValidationLayers = false;
  #endif

  int WIN_WIDTH       = 1024;
  int WIN_HEIGHT      = 1024;
  int PASS_NUMBER     = 1024;
  int NAIVE_PT_REPEAT = 1; // make more samples for naivept which is quite useful for testing cases to get less noise for 

  std::string scenePath      = "../resources/HydraCore/hydra_app/tests/test_42/statex_00001.xml"; 
  std::string sceneDir       = "";          // alternative path of scene library root folder (by default it is the folder where scene xml is located)
  std::string imageOut       = "z_out.bmp";
  std::string integratorType = "mispt";
  float gamma                = 2.4f; // out gamma, special value, see save image functions

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  std::shared_ptr<Integrator> pImpl = nullptr;
  ArgParser args(argc, argv);
  
  if(args.hasOption("-in"))
    scenePath = args.getOptionValue<std::string>("-in");

  if(args.hasOption("-out"))
    imageOut = args.getOptionValue<std::string>("-out");

  std::filesystem::path out_path {imageOut};
  auto dir = out_path.parent_path();
  if(!dir.empty() && !std::filesystem::exists(dir))
    std::filesystem::create_directories(dir);
 
  if(args.hasOption("-scn_dir"))
    sceneDir = args.getOptionValue<std::string>("-scn_dir");

  const bool saveHDR = imageOut.find(".exr") != std::string::npos;
  const std::string imageOutClean = imageOut.substr(0, imageOut.find_last_of("."));

  if(args.hasOption("-integrator"))
    integratorType = args.getOptionValue<std::string>("-integrator");

  if(args.hasOption("-spp-naive-mul"))
    NAIVE_PT_REPEAT = std::max(args.getOptionValue<int>("-spp-naive-mul"),1);
  
  if(args.hasOption("-gamma")) {
    std::string gammaText = args.getOptionValue<std::string>("-gamma");
    if(gammaText == "srgb" || gammaText == "sSRGB")
      gamma = 2.4f;
    else
      gamma = args.getOptionValue<float>("-gamma");
  }
  
  int spectral_mode = args.hasOption("--spectral") ? 1 : 0;

  float4x4 look_at;
  auto override_camera_pos = args.hasOption("-look_at");
  if(override_camera_pos)
  {
    auto str = args.getOptionValue<std::string>("-look_at");
    std::cout << str << std::endl;
    look_at = ReadMatrixFromString(str);
  }
  
  const bool enableNaivePT  = (integratorType == "naivept" || integratorType == "all");
  const bool enableShadowPT = (integratorType == "shadowpt" || integratorType == "all");
  const bool enableMISPT    = (integratorType == "mispt" || integratorType == "all");
  const bool enableRT       = (integratorType == "raytracing" || integratorType == "rt" || integratorType == "whitted_rt");

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  std::cout << "[main]: loading xml ... " << scenePath.c_str() << std::endl;

  auto features = Integrator::PreliminarySceneAnalysis(scenePath.c_str(), sceneDir.c_str(), 
                                                       WIN_WIDTH, WIN_HEIGHT, spectral_mode);

  //// override parameters which are explicitly defined in command line
  //
  if(args.hasOption("-width"))
    WIN_WIDTH = args.getOptionValue<int>("-width");
  if(args.hasOption("-height"))
    WIN_HEIGHT = args.getOptionValue<int>("-height");
  if(args.hasOption("--spectral"))
    spectral_mode = 1;
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  std::vector<float4> realColor(WIN_WIDTH*WIN_HEIGHT);

  bool onGPU = args.hasOption("--gpu");
  #ifdef USE_VULKAN
  if(onGPU)
  { 
    unsigned int a_preferredDeviceId = args.getOptionValue<int>("-gpu_id", 0);

    // simple way
    //
    //auto ctx = vk_utils::globalContextGet(enableValidationLayers, a_preferredDeviceId);
    //pImpl = CreateIntegrator_Generated(WIN_WIDTH*WIN_HEIGHT, spectral_mode, features, ctx, WIN_WIDTH*WIN_HEIGHT);
    
    // advanced way, init device with features which is required by generated class
    //
    std::vector<const char*> requiredExtensions;
    auto deviceFeatures = Integrator_Generated::ListRequiredDeviceFeatures(requiredExtensions);                                          
    auto ctx            = vk_utils::globalContextInit(requiredExtensions, enableValidationLayers, a_preferredDeviceId, &deviceFeatures); 
     
    // advanced way, you can disable some pipelines creation which you don't actually need;
    // this will make application start-up faster
    //
    Integrator_Generated::EnabledPipelines().enableRayTraceMega               = enableRT;
    Integrator_Generated::EnabledPipelines().enableCastSingleRayMega          = false; // not used, for testing only
    Integrator_Generated::EnabledPipelines().enablePackXYMega                 = true;  // always true for this main.cpp;
    Integrator_Generated::EnabledPipelines().enablePathTraceFromInputRaysMega = false; // always false in this main.cpp; see cam_plugin main
    Integrator_Generated::EnabledPipelines().enablePathTraceMega              = enableShadowPT || enableMISPT;
    Integrator_Generated::EnabledPipelines().enableNaivePathTraceMega         = enableNaivePT;

    // advanced way
    //
    auto pObj = std::make_shared<Integrator_Generated>(WIN_WIDTH*WIN_HEIGHT, spectral_mode, features); 
    pObj->SetVulkanContext(ctx);
    pObj->InitVulkanObjects(ctx.device, ctx.physicalDevice, WIN_WIDTH*WIN_HEIGHT); 
    pImpl = pObj;
  }
  else
  #endif
    pImpl = std::make_shared<Integrator>(WIN_WIDTH*WIN_HEIGHT, spectral_mode, features);
  
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  pImpl->SetViewport(0,0,WIN_WIDTH,WIN_HEIGHT);
  std::cout << "[main]: Loading scene ... " << scenePath.c_str() << std::endl;
  pImpl->LoadScene(scenePath.c_str(), sceneDir.c_str());

  if(override_camera_pos)
  {
    pImpl->SetWorldView(look_at);
  }

  pImpl->CommitDeviceData();

  PASS_NUMBER = pImpl->GetSPP();                     // read target spp from scene
  if(args.hasOption("-spp"))                         // override it if spp is specified via command line
    PASS_NUMBER = args.getOptionValue<int>("-spp");

  // remember (x,y) coords for each thread to make our threading 1D
  //
  std::cout << "[main]: PackXYBlock() ... " << std::endl; 
  pImpl->PackXYBlock(WIN_WIDTH, WIN_HEIGHT, 1);

  float timings[4] = {0,0,0,0};
  
  // now test path tracing
  //
  if(enableNaivePT)
  {
    std::cout << "[main]: NaivePathTraceBlock() ... " << std::endl;
    std::fill(realColor.begin(), realColor.end(), LiteMath::float4{});

    pImpl->SetIntegratorType(Integrator::INTEGRATOR_STUPID_PT);
    pImpl->UpdateMembersPlainData();
    pImpl->NaivePathTraceBlock(WIN_WIDTH*WIN_HEIGHT, realColor.data(), PASS_NUMBER*NAIVE_PT_REPEAT);
    
    std::cout << std::endl;
    pImpl->GetExecutionTime("NaivePathTraceBlock", timings);
    std::cout << "NaivePathTraceBlock(exec)  = " << timings[0]              << " ms " << std::endl;
    std::cout << "NaivePathTraceBlock(copy)  = " << timings[1] + timings[2] << " ms " << std::endl;
    std::cout << "NaivePathTraceBlock(ovrh)  = " << timings[3]              << " ms " << std::endl;
    std::cout << std::endl;

    const float normConst = 1.0f/float(PASS_NUMBER*NAIVE_PT_REPEAT);

    if(saveHDR)
    {
      const std::string outName = (integratorType == "naivept") ? imageOut : imageOutClean + "_naivept.exr"; 
      SaveImage4fToEXR((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, outName.c_str(), normConst, true);
    }
    else
    {
      const std::string outName = (integratorType == "naivept") ? imageOut : imageOutClean + "_naivept.bmp"; 
      SaveImage4fToBMP((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, outName.c_str(), normConst, gamma);
    }
  }
  
  const float normConst = 1.0f/float(PASS_NUMBER);
  if(enableShadowPT)
  {
    std::cout << "[main]: PathTraceBlock(Shadow-PT) ... " << std::endl;
    
    std::fill(realColor.begin(), realColor.end(), LiteMath::float4{});

    pImpl->SetIntegratorType(Integrator::INTEGRATOR_SHADOW_PT);
    pImpl->UpdateMembersPlainData();
    pImpl->PathTraceBlock(WIN_WIDTH*WIN_HEIGHT, realColor.data(), PASS_NUMBER);
    
    if(saveHDR) 
    {
      const std::string outName = (integratorType == "shadowpt") ? imageOut : imageOutClean + "_shadowpt.exr"; 
      SaveImage4fToEXR((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, outName.c_str(), normConst, true);
    }
    else
    {
      const std::string outName = (integratorType == "shadowpt") ? imageOut : imageOutClean + "_shadowpt.bmp"; 
      SaveImage4fToBMP((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, outName.c_str(), normConst, gamma);
    }
  }

  if(enableMISPT)
  {
    std::cout << "[main]: PathTraceBlock(MIS-PT) ... " << std::endl;
    
    std::fill(realColor.begin(), realColor.end(), LiteMath::float4{});

    pImpl->SetIntegratorType(Integrator::INTEGRATOR_MIS_PT);
    pImpl->UpdateMembersPlainData();
    pImpl->PathTraceBlock(WIN_WIDTH*WIN_HEIGHT, realColor.data(), PASS_NUMBER);
    
    pImpl->GetExecutionTime("PathTraceBlock", timings);
    std::cout << "PathTraceBlock(exec) = " << timings[0]              << " ms " << std::endl;
    std::cout << "PathTraceBlock(copy) = " << timings[1] + timings[2] << " ms " << std::endl;
    std::cout << "PathTraceBlock(ovrh) = " << timings[3]              << " ms " << std::endl;

    if(saveHDR) 
    {
      const std::string outName = (integratorType == "mispt") ? imageOut : imageOutClean + "_mispt.exr";
      SaveImage4fToEXR((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, outName.c_str(), normConst, true);
    }
    else
    {  
      const std::string outName = (integratorType == "mispt") ? imageOut : imageOutClean + "_mispt.bmp"; 
      SaveImage4fToBMP((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, outName.c_str(), normConst, gamma);
    }
  }
  
  if(enableRT)
  {
    PASS_NUMBER = 1;               // must be always one for RT currently
    const float normConstRT = 1.0f;  // must be always one for RT currently
    std::cout << "[main]: RayBlock ... " << std::endl;

    std::fill(realColor.begin(), realColor.end(), LiteMath::float4{});
   
    pImpl->UpdateMembersPlainData();
    pImpl->RayTraceBlock(WIN_WIDTH*WIN_HEIGHT, realColor.data(), PASS_NUMBER);

    pImpl->GetExecutionTime("RayTraceBlock", timings);
    std::cout << "RayTraceBlock(exec) = " << timings[0]              << " ms " << std::endl;
    std::cout << "RayTraceBlock(copy) = " << timings[1] + timings[2] << " ms " << std::endl;
    std::cout << "RayTraceBlock(ovrh) = " << timings[3]              << " ms " << std::endl;

    if(saveHDR)
    {
      const std::string outName = (integratorType == "raytracing") ? imageOut : imageOutClean + "_rt.exr";
      SaveImage4fToEXR((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, outName.c_str(), normConstRT, true);
    }
    else
    {
      const std::string outName = (integratorType == "raytracing") ? imageOut : imageOutClean + "_rt.bmp";
      SaveImage4fToBMP((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, outName.c_str(), normConstRT, gamma);
    }
  }

  return 0;
}
