#include "integrator_pt.h"
#include <memory>
#include <omp.h>

#include "utils.h" // for progress bar

class IntegratorQMC : public Integrator
{
public:

  IntegratorQMC(int a_maxThreads, int a_spectral_mode, std::vector<uint32_t> a_features) : Integrator(a_maxThreads, a_spectral_mode, a_features)
  {
    int maxThreadNum = 0;
    #pragma omp parallel
    {
      #pragma omp single
      maxThreadNum = omp_get_max_threads();
    }

    std::cout << "[IntegratorQMC]: omp_get_num_threads() = " << maxThreadNum << std::endl;
  }

  float  GetRandomNumbersSpec(RandomGen* a_gen) override;
  float4 GetRandomNumbersLens(RandomGen* a_gen) override;
  float4 GetRandomNumbersMats(RandomGen* a_gen, int a_bounce) override;
  float4 GetRandomNumbersLgts(RandomGen* a_gen, int a_bounce) override;
  float  GetRandomNumbersMatB(RandomGen* a_gen, int a_bounce, int a_layer) override;

  void   PathTraceBlock(uint pixelsNum, uint channels, float* out_color, uint a_passNum) override;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float  IntegratorQMC::GetRandomNumbersSpec(RandomGen* a_gen) 
{ 
  return rndFloat1_Pseudo(a_gen); 
}

float4 IntegratorQMC::GetRandomNumbersLens(RandomGen* a_gen) 
{ 
  return rndFloat4_Pseudo(a_gen); 
}

float4 IntegratorQMC::GetRandomNumbersMats(RandomGen* a_gen, int a_bounce) 
{ 
  return rndFloat4_Pseudo(a_gen); 
}

float4 IntegratorQMC::GetRandomNumbersLgts(RandomGen* a_gen, int a_bounce)
{
  const float4 rands = rndFloat4_Pseudo(a_gen); // don't use single rndFloat4 (!!!)
  const float  rndId = rndFloat1_Pseudo(a_gen); // don't use single rndFloat4 (!!!)
  return float4(rands.x, rands.y, rands.z, rndId);
}

float IntegratorQMC::GetRandomNumbersMatB(RandomGen* a_gen, int a_bounce, int a_layer) 
{ 
  return rndFloat1_Pseudo(a_gen); 
}

void IntegratorQMC::PathTraceBlock(uint pixelsNum, uint channels, float* out_color, uint a_passNum)
{
  ConsoleProgressBar progress(pixelsNum);
  progress.Start();
  auto start = std::chrono::high_resolution_clock::now();
  #ifndef _DEBUG
  #pragma omp parallel for default(shared)
  #endif
  for (int i = 0; i < pixelsNum; ++i) {
    for (int j = 0; j < a_passNum; ++j) {
      PathTrace(uint(i), channels, out_color);
    }
    progress.Update();
  }
  progress.Done();
  shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count())/1000.f;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Integrator> CreateIntegratorQMC(int a_maxThreads = 1, int a_spectral_mode = 0, std::vector<uint32_t> a_features = {}) 
{ 
  return std::make_shared<IntegratorQMC>(a_maxThreads, a_spectral_mode, a_features); 
}