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

struct LensElementInterface 
{
  float curvatureRadius;
  float thickness;
  float eta;
  float apertureRadius;
};

class CamTableLens : public ICamRaysAPI
{
public:
  CamTableLens();
  virtual ~CamTableLens();

  void SetParameters(int a_width, int a_height, const CamParameters& a_params) override;
  void SetBatchSize(int a_tileSize) override { Init(a_tileSize); };
  
  void MakeRaysBlock(RayPart1* out_rayPosAndNear4f [[size("in_blockSize")]], 
                     RayPart2* out_rayDirAndFar4f  [[size("in_blockSize")]], uint32_t in_blockSize, int subPassId)  override;
                     
  void AddSamplesContributionBlock(float* out_color4f [[size("a_width*a_height*4")]], const float* colors4f [[size("in_blockSize*4")]], uint32_t in_blockSize, 
                                   uint32_t a_width, uint32_t a_height, int subPassId) override;

  virtual void CommitDeviceData(){}
  virtual void GetExecutionTime(const char* a_funcName, float a_out[4]){}

protected:

  void kernel1D_MakeEyeRay   (int in_blockSize, RayPart1* out_rayPosAndNear4f, RayPart2* out_rayDirAndFar4f, int subPassId);
  void kernel1D_ContribSample(int in_blockSize, const float4* in_color, float4* out_color, int subPassId);

  bool IntersectSphericalElement(float radius, float zCenter, float3 rayPos, float3 rayDir, 
                                 float *t, float3 *n) const;

  bool TraceLensesFromFilm(float3 inRayPos, float3 inRayDir, 
                           float3* outRayPos, float3* outRayDir) const;

  inline float LensRearZ()      const { return lines[0].thickness; }
  inline float LensRearRadius() const { return lines[0].apertureRadius; }         

  float4x4 m_proj;
  float4x4 m_projInv;

  int m_width;
  int m_height;
  int m_spectral_mode;

  std::vector<RandomGen>   m_randomGens;
  std::vector<uint2>       m_storedWaves;
  std::vector<float>       m_storedCos4;
  void Init(int a_maxThreads);

  std::vector<float> m_cie_x;
  std::vector<float> m_cie_y;
  std::vector<float> m_cie_z;

  std::vector<LensElementInterface> lines;
  float2 m_physSize;
  float m_diagonal;
  float m_aspect;
};
