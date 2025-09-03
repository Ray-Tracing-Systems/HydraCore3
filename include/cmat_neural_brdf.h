#ifndef INCLUDE__CMAT_NEURAL_BRDF_H_
#define INCLUDE__CMAT_NEURAL_BRDF_H_
#include "../neural.h"
#include "cmaterial.h"
#include "include/cglobals.h"

#include <iostream>
static constexpr uint NBRDF_LATENT_CHANNELS = 0;
static constexpr uint NBRDF_INPUT_DIM = NBRDF_LATENT_CHANNELS + 6;
static constexpr uint NBRDF_HIDDEN_DIM = 21;

static constexpr uint NBRDF_OFFSET_BIAS0 = NBRDF_INPUT_DIM * NBRDF_HIDDEN_DIM;
static constexpr uint NBRDF_OFFSET1      = NBRDF_OFFSET_BIAS0 + NBRDF_HIDDEN_DIM;
static constexpr uint NBRDF_SIZE_MAT     = NBRDF_HIDDEN_DIM * NBRDF_HIDDEN_DIM;
static constexpr uint NBRDF_SIZE_LAYER   = NBRDF_SIZE_MAT + NBRDF_HIDDEN_DIM;

static constexpr uint NBRDF_SIZE_MAT6    = NBRDF_HIDDEN_DIM * 6;

static bool isUsed = false;

static constexpr float ALPHA = 0.4;

static inline void neuralBrdfEval(const Material* a_materials, const float *weights,
                                    float3 l, float3 v, float3 n, float *tex, BsdfEval *pRes)
{
  const float2 alpha = float2(ALPHA, ALPHA);
  float3 nx, ny, nz = n;
  CoordinateSystemV2(nz, &nx, &ny);
  const float cosThetaOut = dot(l, n); 
  // v = (-1.0f) * v;
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  const float3 wi = float3(dot(l, nx), dot(l, ny), dot(l, nz));

  if(wo.z * wi.z < 0.0f)
    return;

  float3 wm = wo + wi;
  if (dot(wm, wm) == 0)
      return;

  wm = normalize(wm);
  if (wi.z * wo.z < 0.0f)
  {
    return;
  }

  float3 a, b;
  RusinkiewiczTransform(wi, wo, &a, &b);

  float x[NBRDF_HIDDEN_DIM];
  float buf[NBRDF_HIDDEN_DIM];
  x[0] = a.x;
  x[1] = a.y;
  x[2] = a.z;
  x[3] = b.x;
  x[4] = b.y;
  x[5] = b.z;
  for(uint i = 0; i < NBRDF_LATENT_CHANNELS; ++i) {
      x[6 + i] = tex[i];
  }



  //Layer0
  nn::Linear(weights, x, buf, NBRDF_INPUT_DIM, NBRDF_HIDDEN_DIM);
  nn::ReLU(buf, buf, NBRDF_HIDDEN_DIM);

  uint32_t current_offset = NBRDF_OFFSET1;

  //Layer1
  nn::Linear(weights + current_offset,
              buf, x, NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
  nn::ReLU(x, x, NBRDF_HIDDEN_DIM);
  current_offset += NBRDF_SIZE_LAYER;

  //Layer2
  nn::Linear(weights + current_offset,
              x, buf, NBRDF_HIDDEN_DIM, 3);

  nn::Exp(buf, buf, 3);
  nn::Add(buf, -1.0f, buf, 3);
  nn::ReLU(buf, buf, 3);

  wm        = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
  pRes->pdf = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  pRes->val = float4(buf[0], buf[1], buf[2], 1.0f);
  pRes->val = (cosThetaOut <= 1e-6f) ? float4(0.0f) : (pRes->val * std::max(cosThetaOut, 1e-6f));  
}


float ggxDistribution(float3 n, float3 h, float roughness) 
{
    float alpha = roughness * roughness;
    float alpha2 = alpha * alpha;
    float cosThetaH = std::max(LiteMath::dot(n, h), 0.0f);
    float denom = (cosThetaH * cosThetaH) * (alpha2 - 1.0f) + 1.0f;
    return alpha2 / (LiteMath::M_PI * denom * denom + EPSILON);
}

float schlickSmithG1(float3 n, float3 v, float roughness) 
{
    float k = (roughness + 1.0f) * (roughness + 1.0f) / 8.0f;
    float nv = std::max(LiteMath::dot(n, v), 0.0f);
    return nv / (nv * (1.0f - k) + k + EPSILON);
}


static inline void neuralBrdfSampleAndEval(const Material* a_materials, const float *weights, float4 rands, 
                                            float3 v, float3 n, float *tex, BsdfSample* pRes)
{
  const float2 alpha = float2(ALPHA, ALPHA);
  const float3 l = reflect((-1.0f)*v, n);
  const float cosThetaOut = dot(l, n); 

  float3 nx, ny, nz = n;
  CoordinateSystemV2(nz, &nx, &ny);
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  if(wo.z == 0)
    return;

  if(wo.z == 0)
    return;

  float3 wm = trSample(wo, float2(rands.x, rands.y), alpha);
  float3 wi = reflect((-1.0f) * wo, wm);

  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);

  if (wi.z * wo.z < 0.0f)
  {
    return;
  }

  float3 a, b;
  RusinkiewiczTransform(wi, wo, &a, &b);

  float x[NBRDF_HIDDEN_DIM];
  float buf[NBRDF_HIDDEN_DIM];
  x[0] = a.x;
  x[1] = a.y;
  x[2] = a.z;
  x[3] = b.x;
  x[4] = b.y;
  x[5] = b.z;
  for(uint i = 0; i < NBRDF_LATENT_CHANNELS; ++i) {
      x[6 + i] = tex[i];
  }

  //Layer0
  nn::Linear(weights, x, buf, NBRDF_INPUT_DIM, NBRDF_HIDDEN_DIM);
  nn::ReLU(buf, buf, NBRDF_HIDDEN_DIM);

  uint32_t current_offset = NBRDF_OFFSET1;

  //Layer1
  nn::Linear(weights + current_offset,
              buf, x, NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
  nn::ReLU(x, x, NBRDF_HIDDEN_DIM);
  current_offset += NBRDF_SIZE_LAYER;

  //Layer2
  nn::Linear(weights + current_offset,
              x, buf, NBRDF_HIDDEN_DIM, 3);

  nn::Exp(buf, buf, 3);
  nn::Add(buf, -1.0f, buf, 3);
  nn::ReLU(buf, buf, 3);

  pRes->val = float4(buf[0], buf[1], buf[2], 1.0f);    
  //std::cout << a.x << " " << a.y << " " << a.z << std::endl;
  //std::cout << b.z << " " << b.z << " " << b.z << std::endl;
  //std::cout << b.z << " " << pRes->val.x << std::endl;
  //pRes->val = (cosThetaOut <= 1e-6f) ? float4(0.0f) : (pRes->val / std::max(cosThetaOut, 1e-6f));  
  
  pRes->dir   = normalize(wi.x * nx + wi.y * ny + wi.z * nz);
  pRes->pdf   = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}



#endif