#include <iostream>
#include <fstream>
#include <filesystem>
#include <memory>
#include <iomanip>

#include "integrator_pt.h"
#include "diff_render/integrator_dr.h"
#include "ArgParser.h"

#include "adam.h"

void SaveFrameBufferToEXR(float* data, int width, int height, int channels, const char* outfilename, float a_normConst = 1.0f);
bool SaveImage4fToBMP(const float* rgb, int width, int height, const char* outfilename, float a_normConst = 1.0f, float a_gamma = 2.2f);
std::vector<float> LoadImage4fFromEXR(const char* infilename, int* pW, int* pH);

float4x4 ReadMatrixFromString(const std::string& str);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// // https://www.shadertoy.com/view/WlG3zG
// inline float4 exp2m1(float4 v) { return float4(std::exp2(v.x), std::exp2(v.y), std::exp2(v.z), std::exp2(v.w)) - float4(1.0f); }
// inline float4 pow_22(float4 x) { return (exp2m1(0.718151f*x)-0.503456f*x)*7.07342f; }
// //inline float4 pow_22(float4 x) { x*x*(float4(0.75f) + 0.25f*x); }
                                          // not so fast as pow_22, but more correct implementation
static inline float sRGBToLinear(float s) // https://entropymine.com/imageworsener/srgbformula/
{
  if(s <= 0.0404482362771082f)
    return s*0.077399381f;
  else 
    return std::pow((s+0.055f)*0.947867299f, 2.4f);
}

static inline float4 read_array_uchar4(const uchar4* a_data, int offset)
{
  const float mult = 0.003921568f; // (1.0f/255.0f);
  const uchar4 c0  = a_data[offset];
  return mult*float4((float)c0.x, (float)c0.y, (float)c0.z, (float)c0.w);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  auto pImpl = std::make_shared<IntegratorDR>(FB_WIDTH*FB_HEIGHT, spectral_mode, features);
  
  std::string refImgpath = "z_ref_rt.exr";
  std::cout << "[drmain]: Loading reference image ... " << refImgpath.c_str() << std::endl;

  int refW = 0, refH = 0;
  std::vector<float> refColor = LoadImage4fFromEXR(refImgpath.c_str(), &refW, &refH);
  
  if(refW  == 0 || refH == 0)
  {
    std::cout << "[drmain]: can't load reference image '" << refImgpath.c_str() << "'" << std::endl;
    exit(0);
  }
  else if(refW != FB_WIDTH || refH != FB_HEIGHT)
  {
    std::cout << "[drmain]: bad resolution of reference image, must   be " << FB_WIDTH << ", " << FB_HEIGHT << std::endl;
    std::cout << "[drmain]: bad resolution of reference image, actual is " << refW << ", " << refH << std::endl;
    exit(0);
  }

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

  pImpl->SetIntegratorType(Integrator::INTEGRATOR_MIS_PT);
  pImpl->UpdateMembersPlainData();

  float timings[4] = {0,0,0,0};
  const float normConst = 1.0f/float(PASS_NUMBER);

  std::vector<float> imgData(256*256*4);
  std::vector<float> imgGrad(256*256*4);

  std::fill(imgData.begin(), imgData.end(), 1.0f);
  std::fill(imgGrad.begin(), imgGrad.end(), 0.0f);
  
  if(false)
  {
    std::vector<uchar4> img;
    unsigned wh[2] = { 0,0};
    std::ifstream fin("/home/frol/PROG/HydraRepos/HydraAPI-tests/tests/test_35/data/chunk_00001.image4ub", std::ios::binary);
    if(!fin.is_open())
    {
      std::cout << "[LoadImage<uint>]: can't open file '" << "/home/frol/PROG/HydraRepos/HydraAPI-tests/tests/test_35/data/chunk_00001.image4ub" << "' " << std::endl;
      exit(0);
    }
    fin.read((char*)wh, sizeof(unsigned)* 2);
    img.resize(wh[0]*wh[1]);
    fin.read((char*)img.data(), size_t(wh[0]*wh[1])*sizeof(uint32_t));
    fin.close();
    
    for(size_t i=0;i<imgData.size()/4;i++)
    {
      float4 color = read_array_uchar4(img.data(), i);
      imgData[i*4+0] = sRGBToLinear(color.x);
      imgData[i*4+1] = sRGBToLinear(color.y);
      imgData[i*4+2] = sRGBToLinear(color.z);
      imgData[i*4+3] = sRGBToLinear(color.w);
    }
  }

  std::shared_ptr< IGradientOptimizer<float> > pOpt = std::make_shared< AdamOptimizer2<float> >(imgGrad.size());

  // now run opt loop
  //
  for(int iter = 0; iter < 50; iter++) 
  {
    std::cout << "[drmain]: Render(" << std::setfill('0') << std::setw(2) << iter << ").., ";
    std::cout.flush();
    
    std::fill(realColor.begin(), realColor.end(), 0.0f);

    float loss = pImpl->PathTraceDR(FB_WIDTH*FB_HEIGHT, FB_CHANNELS, realColor.data(), PASS_NUMBER,
                                    refColor.data(), imgData.data(), imgGrad.data(), imgGrad.size());
    
    std::cout << "loss = " << loss << std::endl;
    std::cout.flush();
    
    //pImpl->GetExecutionTime("PathTraceBlock", timings);
    //std::cout << "PathTraceBlock(exec) = " << timings[0] << " ms " << std::endl;
    
    pOpt->step(imgData.data(), imgGrad.data(), iter);
    
    std::stringstream strOut;
    strOut << imageOutClean << std::setfill('0') << std::setw(2) << iter << ".bmp";
    auto outName = strOut.str();
    SaveImage4fToBMP(realColor.data(), FB_WIDTH, FB_HEIGHT, outName.c_str(), normConst, 2.4f);
  }

  return 0;
}
