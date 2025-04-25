#ifndef INCLUDE__CMAT_NEURAL_BRDF_H_
#define INCLUDE__CMAT_NEURAL_BRDF_H_
#include "../neural.h" 

#include <iostream>
static constexpr uint NBRDF_LATENT_CHANNELS = 0;
static constexpr uint NBRDF_INPUT_DIM = NBRDF_LATENT_CHANNELS + 6;
static constexpr uint NBRDF_HIDDEN_DIM = 64;

static constexpr uint NBRDF_OFFSET_BIAS0 = NBRDF_INPUT_DIM * NBRDF_HIDDEN_DIM;
static constexpr uint NBRDF_OFFSET1      = NBRDF_OFFSET_BIAS0 + NBRDF_HIDDEN_DIM;
static constexpr uint NBRDF_SIZE_MAT     = NBRDF_HIDDEN_DIM * NBRDF_HIDDEN_DIM;
static constexpr uint NBRDF_SIZE_LAYER   = NBRDF_SIZE_MAT + NBRDF_HIDDEN_DIM;

static constexpr uint NBRDF_SIZE_MAT6    = NBRDF_HIDDEN_DIM * 6;

static bool isUsed = false;

static inline void neuralBrdfEval(const Material* a_materials, const float *weights,
                                    float3 l, float3 v, float3 n, float *tex, BsdfEval *pRes)
{
  const float cosThetaOut = dot(l, n);
  l = LiteMath::normalize(l);
  v = LiteMath::normalize(v);
  n = LiteMath::normalize(n);
  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);
  const float3 wo = LiteMath::normalize(float3(dot(l, s), dot(l, t), dot(l, n)));
  const float3 wi = LiteMath::normalize(float3(dot(v, s), dot(v, t), dot(v, n)));
  const float3 wm = normalize(wo + wi);



  if (wi.z * wo.z < 0.0f)
  {
    return;
  }

  float x[NBRDF_HIDDEN_DIM];
  float buf[NBRDF_HIDDEN_DIM];
  x[0] = wo.x;
  x[1] = wo.y;
  x[2] = wo.z;
  x[3] = wi.x;
  x[4] = wi.y;
  x[5] = wi.z;
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
              x, buf, NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
  nn::ReLU(buf, buf, NBRDF_HIDDEN_DIM);
  current_offset += NBRDF_SIZE_LAYER;

  //Layer3
  nn::Linear(weights + current_offset,
              buf, x, NBRDF_HIDDEN_DIM, 6);
  nn::ReLU(x, x, 6);
  current_offset += NBRDF_SIZE_MAT6 + 6;

  //Layer4
  nn::Linear(weights + current_offset,
              x, buf, 6, 3);


  pRes->pdf = lambertEvalPDF(l, v, n);
  pRes->val = float4(buf[2], buf[1], buf[0], 1.0f);

}


static inline void neuralBrdfSampleAndEval(const Material* a_materials, const float *weights, float4 rands, 
                                            float3 vec, float3 n, float *tex, BsdfSample* pRes)
{
  const uint   cflags     = a_materials[0].cflags;
  const float3 lambertDir = lambertSample(float2(rands.x, rands.y), vec, n);
  const float  lambertPdf = lambertEvalPDF(lambertDir, vec, n);
  BsdfEval tRes;
  neuralBrdfEval(a_materials, weights, lambertDir, vec, n, tex, &tRes);

  pRes->dir   = lambertDir;
  pRes->val   = tRes.val;
  pRes->pdf   = lambertPdf;
  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}



#endif