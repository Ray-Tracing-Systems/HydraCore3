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

  int WIN_WIDTH       = 1024;
  int WIN_HEIGHT      = 1024;
  int PASS_NUMBER     = 1024;

  std::string scenePath      = "../resources/HydraCore/hydra_app/tests/test_42/statex_00001.xml"; 
  std::string sceneDir       = "";          // alternative path of scene library root folder (by default it is the folder where scene xml is located)
  std::string imageOut       = "z_out.bmp";
  std::string integratorType = "mispt";
  float gamma                = 2.4f; // out gamma, special value, see save image functions

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  std::shared_ptr<Integrator>  pImpl    = nullptr;
  std::shared_ptr<ICamRaysAPI> pCamImpl = nullptr;

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
  
  std::vector<float4> rayPos(WIN_WIDTH*WIN_HEIGHT); ///<! per tile data 
  std::vector<float4> rayDir(WIN_WIDTH*WIN_HEIGHT); ///<! per tile data
  std::vector<float4> rayCol(WIN_WIDTH*WIN_HEIGHT); ///<! per tile data
  
  std::vector<float4> realColor(WIN_WIDTH*WIN_HEIGHT); ///<! frame buffer (TODO: spectral FB?)

  bool onGPU = args.hasOption("--gpu");
  #ifdef USE_VULKAN
  if(onGPU)
  {
    unsigned int a_preferredDeviceId = args.getOptionValue<int>("-gpu_id", 0);
    auto ctx = vk_utils::globalContextGet(enableValidationLayers, a_preferredDeviceId);
    pImpl = CreateIntegrator_Generated(WIN_WIDTH*WIN_HEIGHT, ctx, WIN_WIDTH*WIN_HEIGHT);
  }
  else
  #endif
  {
    pImpl    = std::make_shared<Integrator>(WIN_WIDTH*WIN_HEIGHT);
    pCamImpl = std::make_shared<CamPinHole>(); // (WIN_WIDTH*WIN_HEIGHT);
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  pImpl->SetViewport(0,0,WIN_WIDTH,WIN_HEIGHT);                   /////////////////////////////// TODO: remove it later whet cam API is ready (!!!)
  pCamImpl->SetParameters(WIN_WIDTH, WIN_HEIGHT, {45.0f, 1.0f, 0.01f, 100.0f});
  
  std::cout << "[main_with_cam]: Loading scene ... " << scenePath.c_str() << std::endl;
  pImpl->LoadScene(scenePath.c_str(), sceneDir.c_str());
  pImpl->CommitDeviceData();

  PASS_NUMBER = pImpl->GetSPP();                     // read target spp from scene
  if(args.hasOption("-spp"))                         // override it if spp is specified via command line
    PASS_NUMBER = args.getOptionValue<int>("-spp");

  std::cout << "[main_with_cam]: spp = " << PASS_NUMBER << std::endl;

  // remember (x,y) coords for each thread to make our threading 1D
  //
  std::cout << "[main_with_cam]: PackXYBlock() ... " << std::endl; /////////////////////////////// TODO: remove it later whet cam API is ready (!!!)
  pImpl->PackXYBlock(WIN_WIDTH, WIN_HEIGHT, 1);

  float timings[4] = {0,0,0,0};
  const float normConst = 1.0f/float(PASS_NUMBER);

  //if(integratorType == "mispt" || integratorType == "all")
  // do rendering
  {
    std::cout << "[main_with_cam]: PathTraceBlock(MIS-PT) ... " << std::endl;
    
    std::fill(realColor.begin(), realColor.end(), LiteMath::float4{});

    pImpl->SetIntegratorType(Integrator::INTEGRATOR_MIS_PT);
    pImpl->UpdateMembersPlainData();
    
    int passId = 0;
    pCamImpl->MakeRaysBlock((float*)rayPos.data(), (float*)rayDir.data(), WIN_WIDTH*WIN_HEIGHT, passId);
    pImpl->PathTraceFromInputRaysBlock(WIN_WIDTH*WIN_HEIGHT, rayPos.data(), rayDir.data(), rayCol.data(), PASS_NUMBER);
    pCamImpl->AddSamplesContributionBlock((float*)realColor.data(), (const float*)rayCol.data(), WIN_WIDTH*WIN_HEIGHT, WIN_WIDTH, WIN_HEIGHT, passId);
    //pImpl->PathTraceBlock(WIN_WIDTH*WIN_HEIGHT, realColor.data(), PASS_NUMBER);
    
    pImpl->GetExecutionTime("PathTraceBlock", timings);
    std::cout << "PathTraceBlock(exec) = " << timings[0]              << " ms " << std::endl;
    std::cout << "PathTraceBlock(copy) = " << timings[1] + timings[2] << " ms " << std::endl;
    std::cout << "PathTraceBlock(ovrh) = " << timings[3]              << " ms " << std::endl;

    if(saveHDR) 
      SaveImage4fToEXR((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, imageOut.c_str(), normConst, true);
    else
      SaveImage4fToBMP((const float*)realColor.data(), WIN_WIDTH, WIN_HEIGHT, imageOut.c_str(), normConst, gamma);
  }
  
  return 0;
}
