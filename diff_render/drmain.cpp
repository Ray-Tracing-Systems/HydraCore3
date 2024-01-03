#include <iostream>
#include <fstream>
#include <filesystem>

#include "integrator_pt.h"
#include "ArgParser.h"
#include "mi_materials.h"

void SaveFrameBufferToEXR(float* data, int width, int height, int channels, const char* outfilename, float a_normConst = 1.0f);
bool SaveImage4fToBMP(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, float a_gamma = 2.2f);

float4x4 ReadMatrixFromString(const std::string& str);

int main(int argc, const char** argv)
{
  int FB_WIDTH        = 1024;
  int FB_HEIGHT       = 1024;
  int FB_CHANNELS     = 4;

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

  const bool saveHDR = imageOut.find(".exr") != std::string::npos || 
                       imageOut.find(".image1f") != std::string::npos || 
                       imageOut.find(".image4f") != std::string::npos || 
                       imageOut.find(".image3d1f") != std::string::npos;
                       
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
  std::cout << "[drmain]: loading xml ... " << scenePath.c_str() << std::endl;
  
  SceneInfo sceneInfo = {};
  sceneInfo.spectral  = spectral_mode;
  auto features = Integrator::PreliminarySceneAnalysis(scenePath.c_str(), sceneDir.c_str(), &sceneInfo); 
  FB_WIDTH      = sceneInfo.width;
  FB_HEIGHT     = sceneInfo.height;
  spectral_mode = sceneInfo.spectral;

  //// override parameters which are explicitly defined in command line
  //
  if(args.hasOption("-width"))
    FB_WIDTH = args.getOptionValue<int>("-width");
  if(args.hasOption("-height"))
    FB_HEIGHT = args.getOptionValue<int>("-height");
  if(args.hasOption("-channels"))
    FB_CHANNELS = args.getOptionValue<int>("-channels");
  if(args.hasOption("--spectral"))
    spectral_mode = 1;
  
  if(FB_CHANNELS == 2 || FB_CHANNELS == 3) // we don't support these values currently
    FB_CHANNELS = 4;
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  std::vector<float> realColor(FB_WIDTH*FB_HEIGHT*FB_CHANNELS);
  
  pImpl = std::make_shared<Integrator>(FB_WIDTH*FB_HEIGHT, spectral_mode, features);
  
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  pImpl->SetViewport(0,0,FB_WIDTH,FB_HEIGHT);
  std::cout << "[drmain]: Loading scene ... " << scenePath.c_str() << std::endl;
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
  std::cout << "[drmain]: PackXYBlock() ... " << std::endl; 
  pImpl->PackXYBlock(FB_WIDTH, FB_HEIGHT, 1);

  float timings[4] = {0,0,0,0};
  const float normConst = 1.0f/float(PASS_NUMBER);

  // now test path tracing
  //
  {
    std::cout << "[drmain]: PathTraceBlock(MIS-PT) ... " << std::endl;
    
    std::fill(realColor.begin(), realColor.end(), 0.0f);

    pImpl->SetIntegratorType(Integrator::INTEGRATOR_MIS_PT);
    pImpl->UpdateMembersPlainData();
    pImpl->PathTraceBlock(FB_WIDTH*FB_HEIGHT, FB_CHANNELS, realColor.data(), PASS_NUMBER);
    
    pImpl->GetExecutionTime("PathTraceBlock", timings);
    std::cout << "PathTraceBlock(exec) = " << timings[0]              << " ms " << std::endl;
    std::cout << "PathTraceBlock(copy) = " << timings[1] + timings[2] << " ms " << std::endl;
    std::cout << "PathTraceBlock(ovrh) = " << timings[3]              << " ms " << std::endl;

    if(saveHDR) 
    {
      const std::string outName = (integratorType == "mispt") ? imageOut : imageOutClean + "_mispt.exr";
      SaveFrameBufferToEXR(realColor.data(), FB_WIDTH, FB_HEIGHT, FB_CHANNELS, outName.c_str(), normConst);
    }
    else
    {  
      const std::string outName = (integratorType == "mispt") ? imageOut : imageOutClean + "_mispt.bmp"; 
      SaveImage4fToBMP(realColor.data(), FB_WIDTH, FB_HEIGHT, outName.c_str(), normConst, gamma);
    }
  }

  return 0;
}
