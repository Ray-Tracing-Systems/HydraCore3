/*

Minimal `main` specifically for langevin metropolis integrator
without endless `qmcIsEnabled`, `mltIsEnabled` `enableMISPT` flags.

And without almost all of the features too.

*/

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "ArgParser.h"
#include "imageutils.h"
#include "integrator_pt.h"
#include "mi_materials.h"

std::shared_ptr<Integrator> CreateIntegratorKMLT(int a_maxThreads = 1, int a_spectral_mode = 0, std::vector<uint32_t> a_features = {});
std::shared_ptr<Integrator> CreateIntegratorLangevin(int a_maxThreads = 1, int a_spectral_mode = 0, std::vector<uint32_t> a_features = {});

int main(int argc, const char **argv) {
    int FB_WIDTH = 1024;
    int FB_HEIGHT = 1024;
    int FB_CHANNELS = 4;
    int PASS_NUMBER = 1024;
    std::string scenePath = "";
    std::string sceneDir = ""; // alternative path of scene library root folder (by default it is the folder where scene xml is located)
    std::string imageOut = "z_out.bmp";
    std::string integratorType = "mispt";
    std::string fbLayer = "split_direct_indirect";
    float gamma = 2.4f; // out gamma, special value, see save image functions.

    ///////////////////////////////////////////////////////////////////////////////////////

    ArgParser args(argc, argv);

    if (args.hasOption("-in"))
        scenePath = args.getOptionValue<std::string>("-in");
    if (args.hasOption("-out"))
        imageOut = args.getOptionValue<std::string>("-out");

    std::filesystem::path out_path{imageOut};
    auto dir = out_path.parent_path();
    if (!dir.empty() && !std::filesystem::exists(dir))
        std::filesystem::create_directories(dir);

    if (args.hasOption("-scn_dir"))
        sceneDir = args.getOptionValue<std::string>("-scn_dir");

    const std::string imageOutClean = imageOut.substr(0, imageOut.find_last_of("."));  // "image.png" ==> "image"
    const std::string imageOutFiExt = imageOut.substr(imageOut.find_last_of(".") + 1); // "image.png" ==> ".png"

    ///////////////////////////////////////////////////////////////////////////////////////

    std::cout << "[main]: loading xml ... " << scenePath.c_str() << std::endl;

    SceneInfo sceneInfo = {};
    sceneInfo.spectral = 0;
    auto features = Integrator::PreliminarySceneAnalysis(scenePath.c_str(), sceneDir.c_str(), &sceneInfo);
    FB_WIDTH = sceneInfo.width;
    FB_HEIGHT = sceneInfo.height;
    int spectral_mode = sceneInfo.spectral;

    if (args.hasOption("-width"))
        FB_WIDTH = args.getOptionValue<int>("-width");
    if (args.hasOption("-height"))
        FB_HEIGHT = args.getOptionValue<int>("-height");

    ///////////////////////////////////////////////////////////////////////////////////////

    std::vector<float> realColor(FB_WIDTH * FB_HEIGHT * FB_CHANNELS);

    std::shared_ptr<Integrator> pImpl;
    pImpl = CreateIntegratorLangevin(FB_WIDTH * FB_HEIGHT, spectral_mode, features);
    pImpl->SetViewport(0, 0, FB_WIDTH, FB_HEIGHT);
    std::cout << "[main]: Loading scene ... " << scenePath.c_str() << std::endl;
    pImpl->LoadScene(scenePath.c_str(), sceneDir.c_str());

    PASS_NUMBER = pImpl->GetSPP();
    if (args.hasOption("-spp"))
        PASS_NUMBER = args.getOptionValue<int>("-spp");

    // remember (x,y) coords for each thread to make our threading 1D
    std::cout << "[main]: PackXYBlock() ... " << std::endl;
    pImpl->PackXYBlock(FB_WIDTH, FB_HEIGHT, 1);

    const float normConst = 1.0f / float(PASS_NUMBER);

    pImpl->m_traceDepth = 10;
    pImpl->SetFrameBufferLayer(Integrator::FB_COLOR);
    std::cout << "[main]: PathTraceBlock(MIS-PT) ... " << std::endl;
    std::fill(realColor.begin(), realColor.end(), 0.0f);
    pImpl->SetIntegratorType(Integrator::INTEGRATOR_MIS_PT);
    pImpl->UpdateMembersPlainData();
    pImpl->PathTraceBlock(FB_WIDTH * FB_HEIGHT, FB_CHANNELS, realColor.data(), PASS_NUMBER);

    std::cout << "[main]: save final image to " << imageOut.c_str() << std::endl;
    SaveImage4fToBMP(realColor.data(), FB_WIDTH, FB_HEIGHT, FB_CHANNELS, imageOut.c_str(), normConst, gamma);

    return 0;
}
