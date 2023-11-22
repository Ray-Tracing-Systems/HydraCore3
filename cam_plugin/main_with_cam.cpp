#include <iostream>
#include <fstream>
#include <filesystem>

#include "integrator_pt.h"
#include "ArgParser.h"
#include "CamPluginAPI.h"
#include "CamPinHole.h"
#include "CamTableLens.h"

bool SaveImage4fToEXR(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, bool a_invertY = false);
bool SaveImage4fToBMP(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, float a_gamma = 2.2f);

#ifdef USE_VULKAN
#include "vk_context.h"
#include "integrator_pt1_generated.h"              // advanced way of woking with hydra
#include "cam_plugin/CamPinHole_pinhole_gpu.h"     // same way for camera plugins
#include "cam_plugin/CamTableLens_tablelens_gpu.h" // same way for camera plugins
#endif

int main(int argc, const char** argv)
{
  #ifndef NDEBUG
  bool enableValidationLayers = true;
  #else
  bool enableValidationLayers = false;
  #endif

  int WIN_WIDTH  = 1024;
  int WIN_HEIGHT = 1024;
  int SPP_TOTAL  = 1024;

  std::string scenePath      = "../resources/HydraCore/hydra_app/tests/test_42/statex_00001.xml"; 
  std::string sceneDir       = "";          // alternative path of scene library root folder (by default it is the folder where scene xml is located)
  std::string imageOut       = "z_out.bmp";
  std::string integratorType = "mispt";
  float gamma                = 2.4f; // out gamma, special value, see save image functions

  std::shared_ptr<Integrator>  pRender  = nullptr;
  std::shared_ptr<ICamRaysAPI> pCamImpl = nullptr;

  ArgParser args(argc, argv);
  
  if(args.hasOption("-in"))
    scenePath = args.getOptionValue<std::string>("-in");

  if(args.hasOption("-out"))
    imageOut = args.getOptionValue<std::string>("-out");

  if(args.hasOption("-scn_dir"))
    sceneDir = args.getOptionValue<std::string>("-scn_dir");

  int camType = 1;

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  int spectral_mode = args.hasOption("--spectral") ? 1 : 0;
 
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

  std::filesystem::path out_path {imageOut};
  auto dir = out_path.parent_path();
  if(!dir.empty() && !std::filesystem::exists(dir))
    std::filesystem::create_directories(dir);

  const bool saveHDR = imageOut.find(".exr") != std::string::npos;
  const std::string imageOutClean = imageOut.substr(0, imageOut.find_last_of("."));
  
  if(args.hasOption("-gamma")) {
    std::string gammaText = args.getOptionValue<std::string>("-gamma");
    if(gammaText == "srgb" || gammaText == "sSRGB")
      gamma = 2.4f;
    else
      gamma = args.getOptionValue<float>("-gamma");
  }
  
  if(args.hasOption("-width"))
    WIN_WIDTH = args.getOptionValue<int>("-width");
  if(args.hasOption("-height"))
    WIN_HEIGHT = args.getOptionValue<int>("-height");
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  
  const int MEGA_TILE_SIZE = 512*512;           ///<! tile size

  std::vector<RayPart1> rayPos(MEGA_TILE_SIZE); ///<! per tile data, input 
  std::vector<RayPart2> rayDir(MEGA_TILE_SIZE); ///<! per tile data, input
  std::vector<float4>   rayCol(MEGA_TILE_SIZE); ///<! per tile data, output 

  std::vector<float4>   realColor(WIN_WIDTH*WIN_HEIGHT);             ///<! frame buffer (TODO: monochrome FB need also)
  std::fill(realColor.begin(), realColor.end(), LiteMath::float4{}); // clear frame buffer

  bool onGPU = args.hasOption("--gpu");
  #ifdef USE_VULKAN
  if(onGPU)
  {
    // (1) advanced way, you may disable unused features in shader code via spec constants.
    //     To do this, you have to know what materials, lights and e.t.c. is actualle presented in scene 
    //
    std::cout << "[main]: loading xml ... " << scenePath.c_str() << std::endl;
    auto hydraFeatures = Integrator::PreliminarySceneAnalysis(scenePath.c_str(), sceneDir.c_str(), WIN_WIDTH, WIN_HEIGHT, spectral_mode);
    
    // (2) init device with apropriate features for both hydra and camera plugin
    //
    unsigned int preferredDeviceId = args.getOptionValue<int>("-gpu_id", 0);
    std::vector<const char*> requiredExtensions;
    
    auto devFeaturesCam = (camType == 0) ? CamPinHole_PINHOLE_GPU::ListRequiredDeviceFeatures(requiredExtensions) :
                                           CamTableLens_TABLELENS_GPU::ListRequiredDeviceFeatures(requiredExtensions);
                                           
    auto devFeaturesHydra = Integrator_Generated::ListRequiredDeviceFeatures(requiredExtensions); 
    
    // TBD: you actually need to carefully join all required device features structures and Vulkan lists 
    //
    if(devFeaturesCam.features.shaderFloat64 == VK_TRUE) // in this example we know that hydra3 don't use double precition  
      devFeaturesHydra.features.shaderFloat64 = VK_TRUE; // while cam plugin probably uses it ... 
    
    auto ctx = vk_utils::globalContextInit(requiredExtensions, enableValidationLayers, preferredDeviceId, &devFeaturesHydra); 

    // (3) Explicitly disable all pipelines which you don't need.
    //     This will make application start-up faster.
    //
    Integrator_Generated::EnabledPipelines().enableRayTraceMega               = false;
    Integrator_Generated::EnabledPipelines().enableCastSingleRayMega          = false; 
    Integrator_Generated::EnabledPipelines().enablePackXYMega                 = false; 
    Integrator_Generated::EnabledPipelines().enablePathTraceFromInputRaysMega = true;  // you need only this pipeline!
    Integrator_Generated::EnabledPipelines().enablePathTraceMega              = false;
    Integrator_Generated::EnabledPipelines().enableNaivePathTraceMega         = false;

    // advanced way, init renderer
    //
    {
      auto pObj = std::make_shared<Integrator_Generated>(MEGA_TILE_SIZE, spectral_mode, hydraFeatures); 
      pObj->SetVulkanContext(ctx);
      pObj->InitVulkanObjects(ctx.device, ctx.physicalDevice, MEGA_TILE_SIZE); 
      pRender = pObj;
    }

    // init appropriate camera plugin and put it to 'pCamImpl'
    //
    if(camType == 0)
    {
      auto pObj = std::make_shared<CamPinHole_PINHOLE_GPU>(); 
      pObj->SetVulkanContext(ctx);
      pObj->InitVulkanObjects(ctx.device, ctx.physicalDevice, MEGA_TILE_SIZE); 
      pCamImpl = pObj;
    }
    else if(camType == 1)
    {
      auto pObj = std::make_shared<CamTableLens_TABLELENS_GPU>(); 
      pObj->SetVulkanContext(ctx);
      pObj->InitVulkanObjects(ctx.device, ctx.physicalDevice, MEGA_TILE_SIZE); 
      pCamImpl = pObj;
    }
  }
  else
  #endif
  {
    pRender = std::make_shared<Integrator>(MEGA_TILE_SIZE, spectral_mode);

    if(camType == 0)
      pCamImpl = std::make_shared<CamPinHole>(); // (WIN_WIDTH*WIN_HEIGHT);
    else if(camType == 1)
      pCamImpl = std::make_shared<CamTableLens>(); // (WIN_WIDTH*WIN_HEIGHT);
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  
  std::cout << "[main_with_cam]: Loading scene ... " << scenePath.c_str() << std::endl;

  pCamImpl->SetParameters(WIN_WIDTH, WIN_HEIGHT, {45.0f, 1.0f, 0.01f, 100.0f, spectral_mode});
  pCamImpl->SetBatchSize(MEGA_TILE_SIZE);

  pRender->LoadScene(scenePath.c_str(), sceneDir.c_str());
  pRender->SetIntegratorType(Integrator::INTEGRATOR_MIS_PT);

  pRender->CommitDeviceData();
  pCamImpl->CommitDeviceData();

  SPP_TOTAL = pRender->GetSPP();                     // read target spp from scene
  if(args.hasOption("-spp"))                         // override it if spp is specified via command line
    SPP_TOTAL = args.getOptionValue<int>("-spp");

  int SAMPLES_PER_RAY = 1;                          
  int CAM_PASSES_NUM  = SPP_TOTAL/SAMPLES_PER_RAY; 
  if(SPP_TOTAL == 1) {
    SAMPLES_PER_RAY = 1;
    CAM_PASSES_NUM = 1;
  }

  std::cout << "[main_with_cam]: spp     = " << SPP_TOTAL << std::endl;
  std::cout << "[main_with_cam]: passNum = " << CAM_PASSES_NUM << std::endl;

  float timings   [4] = {0,0,0,0};
  float timingSum[4] = {0,0,0,0};
  const float normConst = 1.0f/float(SPP_TOTAL);

  // do rendering
  //
  for(int passId = 0; passId < CAM_PASSES_NUM; passId++)
  {
    std::cout << "rendering, pass " << passId << " / " << CAM_PASSES_NUM  << "\r"; 
    std::cout.flush();
    
    const int passNum = (WIN_WIDTH*WIN_HEIGHT/MEGA_TILE_SIZE);
    for(int subPassId = 0; subPassId < passNum; subPassId++) 
    {
      std::fill(rayCol.begin(), rayCol.end(), LiteMath::float4{}); // clear temp color buffer, gpu ver should do this automaticly, please check(!!!)

      pCamImpl->MakeRaysBlock(rayPos.data(), rayDir.data(), MEGA_TILE_SIZE, subPassId);
      pRender->PathTraceFromInputRaysBlock(MEGA_TILE_SIZE, rayPos.data(), rayDir.data(), rayCol.data(), SAMPLES_PER_RAY);
      pCamImpl->AddSamplesContributionBlock((float*)realColor.data(), (const float*)rayCol.data(), MEGA_TILE_SIZE, WIN_WIDTH, WIN_HEIGHT, subPassId);
      
      pRender->GetExecutionTime("PathTraceFromInputRaysBlock", timings);
      for(int i=0;i<4;i++)
        timingSum[i] += timings[i];
    }
  }

  std::cout << std::endl << std::endl;
  std::cout << "PathTraceFromInputRays(exec, total) = " << timingSum[0]                << " ms " << std::endl;
  std::cout << "PathTraceFromInputRays(copy, total) = " << timingSum[1] + timingSum[2] << " ms " << std::endl;
  std::cout << "PathTraceFromInputRays(ovrh, total) = " << timingSum[3]                << " ms " << std::endl;

  if(saveHDR) 
    SaveImage4fToEXR((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, imageOut.c_str(), normConst, true);
  else
    SaveImage4fToBMP((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, imageOut.c_str(), normConst, gamma);
  
  return 0;
}
