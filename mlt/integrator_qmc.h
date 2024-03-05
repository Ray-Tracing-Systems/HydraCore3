#pragma once 

#include "integrator_pt.h"
#include "rnd_qmc.h"

#include <memory>
#include <limits>

class IntegratorQMC : public Integrator
{
public:

  IntegratorQMC(int a_maxThreads, int a_spectral_mode, std::vector<uint32_t> a_features) : Integrator(a_maxThreads, a_spectral_mode, a_features)
  {
    qmc::init(m_qmcTable);
  }

  float  GetRandomNumbersSpec(uint tid, RandomGen* a_gen) override;
  float4 GetRandomNumbersLens(uint tid, RandomGen* a_gen) override;
  float4 GetRandomNumbersMats(uint tid, RandomGen* a_gen, int a_bounce) override;
  float4 GetRandomNumbersLgts(uint tid, RandomGen* a_gen, int a_bounce) override;
  float  GetRandomNumbersMatB(uint tid, RandomGen* a_gen, int a_bounce, int a_layer) override;

  void   kernel_InitEyeRay2(uint tid, const uint* packedXY, 
                            float4* rayPosAndNear, float4* rayDirAndFar, float4* wavelengths, 
                            float4* accumColor,    float4* accumuThoroughput,
                            RandomGen* gen, uint* rayFlags, MisData* misData) override;

  void   kernel_ContributeToImage(uint tid, const uint* rayFlags, uint channels, const float4* a_accumColor, const RandomGen* gen,
                                  const uint* in_pakedXY, const float4* wavelengths, float* out_color) override;

  void   PathTraceBlock(uint pixelsNum, uint channels, float* out_color, uint a_passNum) override;
  unsigned int m_qmcTable[qmc::QRNG_DIMENSIONS][qmc::QRNG_RESOLUTION];
};
