#include <iostream>
#include <fstream>
#include <filesystem>

#include "integrator_pt.h"
#include "ArgParser.h"
#include "CamPluginAPI.h"
#include "CamPinHole.h"

bool SaveImage4fToEXR(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, bool a_invertY = false);
bool SaveImage4fToBMP(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, float a_gamma = 2.2f);

#ifdef USE_VULKAN
#include "vk_context.h"
std::shared_ptr<Integrator> CreateIntegrator_Generated(int a_maxThreads, vk_utils::VulkanContext a_ctx, size_t a_maxThreadsGenerated);
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

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  int spectral_mode = args.hasOption("--spectral") ? 1 : 0;
  std::cout << "[main_with_cam]: loading xml ... " << scenePath.c_str() << std::endl;

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
  
  std::vector<float4>     rayPos(WIN_WIDTH*WIN_HEIGHT); ///<! per tile data, input 
  std::vector<float4>     rayDir(WIN_WIDTH*WIN_HEIGHT); ///<! per tile data, input
  std::vector<AuxRayData> auxDat(WIN_WIDTH*WIN_HEIGHT); ///<! per tile data, input
  std::vector<float4>     rayCol(WIN_WIDTH*WIN_HEIGHT); ///<! per tile data, output 
  
  std::vector<float4> realColor(WIN_WIDTH*WIN_HEIGHT);  ///<! frame buffer (TODO: spectral FB?)
  std::fill(realColor.begin(), realColor.end(), LiteMath::float4{}); // clear frame buffer

  bool onGPU = args.hasOption("--gpu");
  #ifdef USE_VULKAN
  if(onGPU)
  {
    unsigned int a_preferredDeviceId = args.getOptionValue<int>("-gpu_id", 0);
    auto ctx = vk_utils::globalContextGet(enableValidationLayers, a_preferredDeviceId);
    pRender = CreateIntegrator_Generated(WIN_WIDTH*WIN_HEIGHT, ctx, WIN_WIDTH*WIN_HEIGHT);
  }
  else
  #endif
  {
    pRender  = std::make_shared<Integrator>(WIN_WIDTH*WIN_HEIGHT, spectral_mode);
    pCamImpl = std::make_shared<CamPinHole>(); // (WIN_WIDTH*WIN_HEIGHT);
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  
  std::cout << "[main_with_cam]: Loading scene ... " << scenePath.c_str() << std::endl;

  pCamImpl->SetParameters(WIN_WIDTH, WIN_HEIGHT, {45.0f, 1.0f, 0.01f, 100.0f, spectral_mode});
  pCamImpl->SetBatchSize(WIN_WIDTH*WIN_HEIGHT);

  pRender->LoadScene(scenePath.c_str(), sceneDir.c_str());
  pRender->SetIntegratorType(Integrator::INTEGRATOR_MIS_PT);
  pRender->CommitDeviceData();

  SPP_TOTAL = pRender->GetSPP();                     // read target spp from scene
  if(args.hasOption("-spp"))                         // override it if spp is specified via command line
    SPP_TOTAL = args.getOptionValue<int>("-spp");

  int SAMPLES_PER_RAY = 4;
  int CAM_PASSES_NUM  = SPP_TOTAL/SAMPLES_PER_RAY;
  if(SPP_TOTAL == 1) {
    SAMPLES_PER_RAY = 1;
    CAM_PASSES_NUM = 1;
  }

  std::cout << "[main_with_cam]: spp     = " << SPP_TOTAL << std::endl;
  std::cout << "[main_with_cam]: passNum = " << CAM_PASSES_NUM << std::endl;

  float timings   [4] = {0,0,0,0};
  float timingsAvg[4] = {0,0,0,0};
  const float normConst = 1.0f/float(SPP_TOTAL);

  // do rendering
  //
  for(int passId = 0; passId < CAM_PASSES_NUM; passId++)
  {
    std::cout << "rendering, pass " << passId << " / " << CAM_PASSES_NUM  << "\r"; 
    std::cout.flush();
    std::fill(rayCol.begin(), rayCol.end(), LiteMath::float4{}); // clear temp color buffer, gpu ver should do this automaticly, please check(!!!)
    
    pCamImpl->MakeRaysBlock((float*)rayPos.data(), (float*)rayDir.data(), auxDat.data(), WIN_WIDTH*WIN_HEIGHT, passId);
    pRender->PathTraceFromInputRaysBlock(WIN_WIDTH*WIN_HEIGHT, rayPos.data(), rayDir.data(), (const LiteMath::ushort4*)auxDat.data(), rayCol.data(), SAMPLES_PER_RAY);
    pCamImpl->AddSamplesContributionBlock((float*)realColor.data(), (const float*)rayCol.data(), auxDat.data(), WIN_WIDTH*WIN_HEIGHT, WIN_WIDTH, WIN_HEIGHT, passId);
    
    pRender->GetExecutionTime("PathTraceFromInputRaysBlock", timings);
    for(int i=0;i<4;i++)
      timingsAvg[i] += timings[i];
  }

  std::cout << std::endl << std::endl;
  std::cout << "PathTraceFromInputRays(exec, total) = " << timingsAvg[0]                 << " ms " << std::endl;
  std::cout << "PathTraceFromInputRays(copy, total) = " << timingsAvg[1] + timingsAvg[2] << " ms " << std::endl;
  std::cout << "PathTraceFromInputRays(ovrh, total) = " << timingsAvg[3]                 << " ms " << std::endl;

  if(saveHDR) 
    SaveImage4fToEXR((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, imageOut.c_str(), normConst, true);
  else
    SaveImage4fToBMP((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, imageOut.c_str(), normConst, gamma);
  
  return 0;
}
