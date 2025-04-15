#ifndef INCLUDE__CMAT_NEURAL_SPEC_H_
#define INCLUDE__CMAT_NEURAL_SPEC_H_
#include "../neural.h" 
#include "cmaterial.h"

#include <iostream>
static constexpr uint NSPEC_INPUT_DIM = 7;
static constexpr uint NSPEC_HIDDEN_DIM = 24;

static constexpr uint NSPEC_OFFSET_BIAS0 = NSPEC_INPUT_DIM * NSPEC_HIDDEN_DIM;
static constexpr uint NSPEC_OFFSET1      = NSPEC_OFFSET_BIAS0 + NSPEC_HIDDEN_DIM;
static constexpr uint NSPEC_SIZE_MAT     = NSPEC_HIDDEN_DIM * NSPEC_HIDDEN_DIM;
static constexpr uint NSPEC_SIZE_LAYER   = NSPEC_SIZE_MAT + NSPEC_HIDDEN_DIM;

static constexpr uint NSPEC_SIZE_MAT6    = NSPEC_HIDDEN_DIM * 6;

static inline void neuralSpecSmoothSampleAndEval(const Material* a_materials, const float *weights,
                                    float3 v, float3 n, float4 wavelengths, float *tex, BsdfSample* pRes)
{
  const float3 l = reflect((-1.0f)*v, n);
  const float cosThetaOut = dot(l, n);
  float3 dir              = l;
  float  pdf              = 1.0f;

  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);
  const float3 wo = LiteMath::normalize(float3(dot(l, s), dot(l, t), dot(l, n)));
  const float3 wi = LiteMath::normalize(float3(dot(v, s), dot(v, t), dot(v, n)));
  const float3 wm = normalize(wo + wi);

  if (wi.z * wo.z < 0.0f)
  {
    return;
  }

  float4 wl = LiteMath::clamp((wavelengths - 360.0f) / (830.0f - 360.0f), 0.0f, 1.0f);


  float x[NSPEC_HIDDEN_DIM];
  float buf[NSPEC_HIDDEN_DIM];
  
  float4 val;
  for(uint32_t i = 0; i < SPECTRUM_SAMPLE_SZ; ++i)
  {
    x[0] = wo.x;
    x[1] = wo.y;
    x[2] = wo.z;
    x[3] = wi.x;
    x[4] = wi.y;
    x[5] = wi.z;
    x[6] = wl[i];    

    //Layer0
    nn::Linear(weights, x, buf, NSPEC_INPUT_DIM, NSPEC_HIDDEN_DIM);
    nn::SiLU(buf, buf, NSPEC_HIDDEN_DIM);

    uint32_t current_offset = NSPEC_OFFSET1;

    //Layer1
    nn::Linear(weights + current_offset,
                buf, x, NSPEC_HIDDEN_DIM, NSPEC_HIDDEN_DIM);
    nn::SiLU(x, x, NSPEC_HIDDEN_DIM);
    current_offset += NSPEC_SIZE_LAYER;

    //Layer2
    nn::Linear(weights + current_offset,
                x, buf, NSPEC_HIDDEN_DIM, NSPEC_HIDDEN_DIM);
    nn::SiLU(buf, buf, NSPEC_HIDDEN_DIM);
    current_offset += NSPEC_SIZE_LAYER;

    //Layer3
    nn::Linear(weights + current_offset,
                buf, x, NSPEC_HIDDEN_DIM, NSPEC_HIDDEN_DIM);
    nn::SiLU(x, x, NSPEC_HIDDEN_DIM);
    current_offset += NSPEC_SIZE_LAYER;

    //Layer4
    nn::Linear(weights + current_offset,
                x, buf, NSPEC_HIDDEN_DIM, 6);
    nn::SiLU(buf, buf, 6);
    current_offset += NSPEC_SIZE_MAT6 + 6;

    //Layer5
    nn::Linear(weights + current_offset,
                buf, x, 6, 1);
    nn::ReLU(buf, buf, 1);

    val[i] = x[0];
    //std::cout << val[i] << " ";
    val[i] = (cosThetaOut <= 1e-6f) ? 0.0f : (val[i] / std::max(cosThetaOut, 1e-6f));  
  }
  //std::cout << std::endl;
  pRes->val = val; 
  pRes->dir = dir;
  pRes->pdf = pdf;
  pRes->flags = RAY_EVENT_S;
}

#endif