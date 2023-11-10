#pragma once
#include "CamPluginAPI.h"

class CamPinHole : public ICamRaysAPI
{
  CamPinHole();
  virtual ~CamPinHole();

  void SetParameters(int a_width, int a_height, const float a_projInvMatrix[16], const CamParameters& a_params) override;
  void MakeRaysBlock(float* out_rayPosAndNear4f, float* out_rayDirAndFar4f, size_t in_blockSize, int passId)    override;
  void AddSamplesContributionBlock(float* out_color4f, const float* colors4f, size_t in_blockSize, 
                                   uint32_t a_width, uint32_t a_height, int passId);
};
