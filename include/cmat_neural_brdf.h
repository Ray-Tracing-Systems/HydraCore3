#ifndef INCLUDE__CMAT_NEURAL_BRDF_H_
#define INCLUDE__CMAT_NEURAL_BRDF_H_

#include "../neural.h"
#include "cglobals.h"
#include "cmaterial.h"

#include <iostream>

static constexpr uint NBRDF_INPUT_DIM = 6;
static constexpr uint NBRDF_HIDDEN_DIM = 64;

static constexpr uint NBRDF_BATCH_SIZE = 1;
static constexpr size_t NBRDF_WEIGTH_OFFSETS[] = {0, 448, 4608, 8768, 12928};

static constexpr uint NBRDF_SIZE_MAT6 = NBRDF_HIDDEN_DIM * 6;


/**
 * x is of size max(NBRDF_INPUT_DIM, out_dim))
 */
static inline void evalNeuralNetwork(const float *weights, float *x, uint out_dim)
{
  float buf0[NBRDF_HIDDEN_DIM];
  float buf1[NBRDF_HIDDEN_DIM];

  //Layer0
  nn::Linear(weights + NBRDF_WEIGTH_OFFSETS[0], x, buf1,
             NBRDF_BATCH_SIZE, NBRDF_INPUT_DIM, NBRDF_HIDDEN_DIM);
  nn::ReLU(buf1, buf1, NBRDF_BATCH_SIZE * NBRDF_HIDDEN_DIM);

  //Layer1
  nn::Linear(weights + NBRDF_WEIGTH_OFFSETS[1], buf1, buf0,
             NBRDF_BATCH_SIZE, NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
  nn::ReLU(buf0, buf0, NBRDF_BATCH_SIZE * NBRDF_HIDDEN_DIM);

  //Layer2
  nn::Linear(weights + NBRDF_WEIGTH_OFFSETS[2], buf0, buf1,
             NBRDF_BATCH_SIZE, NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
  nn::ReLU(buf1, buf1, NBRDF_BATCH_SIZE * NBRDF_HIDDEN_DIM);

  //Layer3
  nn::Linear(weights + NBRDF_WEIGTH_OFFSETS[3], buf1, buf0,
             NBRDF_BATCH_SIZE, NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
  nn::ReLU(buf0, buf0, NBRDF_BATCH_SIZE * NBRDF_HIDDEN_DIM);

  //Layer4
  nn::Linear(weights + NBRDF_WEIGTH_OFFSETS[4], buf0, x, 
             NBRDF_BATCH_SIZE, NBRDF_HIDDEN_DIM, out_dim);
}

static inline void neuralBrdfEval(const Material* a_materials, const float *weights, float4 wavelengths,
                                    float3 l, float3 v, float3 n, BsdfEval *pRes, int spectral_mode)
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

  uint out_dim = spectral_mode == 0 ? 3 : 32;

  float x[NBRDF_HIDDEN_DIM];
  x[0] = wo.x;
  x[1] = wo.y;
  x[2] = wo.z;
  x[3] = wi.x;
  x[4] = wi.y;
  x[5] = wi.z;
  evalNeuralNetwork(weights, x, out_dim)

  

  pRes->pdf = lambertEvalPDF(l, v, n); //TODO
  pRes->val = float4(buf[2], buf[1], buf[0], 1.0f);

}


static inline void neuralBrdfSampleAndEval(const Material* a_materials, const float *weights, float4 wavelengths, float4 rands, 
                                            float3 vec, float3 n, BsdfSample* pRes)
{
  const uint   cflags     = a_materials[0].cflags;
  const float3 lambertDir = lambertSample(float2(rands.x, rands.y), vec, n);
  const float  lambertPdf = lambertEvalPDF(lambertDir, vec, n);
  BsdfEval tRes;
  neuralBrdfEval(a_materials, weights, wavelengths, lambertDir, vec, n, &tRes);

  pRes->dir   = lambertDir;
  pRes->val   = tRes.val;
  pRes->pdf   = lambertPdf; //TODO
  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}



#endif