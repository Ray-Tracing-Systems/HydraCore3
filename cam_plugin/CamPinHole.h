#pragma once
#include "CamPluginAPI.h"

#include "LiteMath.h"
#include "../include/cglobals.h"
#include "../include/crandom.h"
#include "../spectrum.h"

#include <vector>

using LiteMath::float4x4;
using LiteMath::float4;
using LiteMath::float3;

class CamPinHole : public ICamRaysAPI
{
public:
  CamPinHole();
  virtual ~CamPinHole();

  void SetParameters(int a_width, int a_height, const CamParameters& a_params) override;
  void SetBatchSize(int a_tileSize) override { Init(a_tileSize); };
  
  void MakeRaysBlock(RayPart1* out_rayPosAndNear4f, RayPart2* out_rayDirAndFar4f, size_t in_blockSize, int passId)  override;
  void AddSamplesContributionBlock(float* out_color4f, const float* colors4f, size_t in_blockSize, 
                                   uint32_t a_width, uint32_t a_height, int passId);

protected:

  void kernel1D_MakeEyeRay   (int in_blockSize, RayPart1* out_rayPosAndNear4f, RayPart2* out_rayDirAndFar4f);
  void kernel1D_ContribSample(int in_blockSize, const float4* in_color, float4* out_color);

  float4x4 m_proj;
  float4x4 m_projInv;

  int m_width;
  int m_height;
  int m_spectral_mode;

  std::vector<RandomGen>  m_randomGens;
  std::vector<uint2>      m_storedWaves;
  void Init(int a_maxThreads);

  std::vector<float> m_cie_lambda;
  std::vector<float> m_cie_x;
  std::vector<float> m_cie_y;
  std::vector<float> m_cie_z;
};
