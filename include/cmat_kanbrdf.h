#ifndef INCLUDE__CMAT_KANBRDF_H_
#define INCLUDE__CMAT_KANBRDF_H_

#include "../neural.h"
#include "cglobals.h"
#include "cmaterial.h"
#include "../spectrum.h"

#include <cmath>
#include <iostream>
#include <string>


static constexpr float KANBRDF_GRID_MIN = -1.0;
static constexpr float KANBRDF_GRID_MAX = 1.0;
static constexpr uint KANBRDF_GRID_SIZE = 6;

static constexpr size_t KANBRDF_MAX_SIZE = 6;

static constexpr uint KANBRDF_LAYER_COUNT = 3;
static constexpr size_t KANBRDF_LAYER_SIZES[KANBRDF_LAYER_COUNT + 1] = {6, 5, 5, 3};
static constexpr size_t KANBRDF_WEIGTH_OFFSETS[KANBRDF_LAYER_COUNT] = {
    ((KANBRDF_GRID_SIZE + 1) * KANBRDF_LAYER_SIZES[0] + 1) * KANBRDF_LAYER_SIZES[1],
    ((KANBRDF_GRID_SIZE + 1) * KANBRDF_LAYER_SIZES[1] + 1) * KANBRDF_LAYER_SIZES[2],
    ((KANBRDF_GRID_SIZE + 1) * KANBRDF_LAYER_SIZES[2] + 1) * KANBRDF_LAYER_SIZES[3]
};


/**
 * x is of size max(KANBRDF_INPUT_DIM, out_dim))
 */
static inline void evalKanLayer(const float *weights, const float *x, float *y, uint in_dim, uint out_dim)
{
  float grid_step = (KANBRDF_GRID_MAX - KANBRDF_GRID_MIN) / (KANBRDF_GRID_SIZE - 1);

  float spline_basis[KANBRDF_GRID_SIZE * KANBRDF_MAX_SIZE];
  for(uint i = 0; i < in_dim; ++i) {
    for(uint j = 0; j < KANBRDF_GRID_SIZE; ++j) {
      float grid_val = KANBRDF_GRID_MIN + grid_step * j;
      float t = (x[i] - grid_val) / grid_step;

      spline_basis[i * KANBRDF_GRID_SIZE + j] = expf(-t * t);
    }
  }
  nn::Matmul(spline_basis, weights, y, 1, KANBRDF_GRID_SIZE * in_dim, out_dim);


  float activated[KANBRDF_MAX_SIZE];
  float base[KANBRDF_MAX_SIZE];
  nn::SiLU(x, activated, in_dim);
  nn::Linear(weights + KANBRDF_GRID_SIZE * in_dim, activated, base, 1, in_dim, out_dim);
  nn::Add(y, base, y, out_dim);
}

static inline float kanbrdf_output_transform(float x) {
  float val = expf(x) - 1;
  return val < 0 ? 0 : val;
}

static inline void kanBrdfEval(const Material* a_materials, const float *weights,
                                    float3 l, float3 v, float3 n, BsdfEval *pRes, int spectral_mode)
{

  //const float cosThetaOut = dot(l, n);
  l = LiteMath::normalize(l);
  v = LiteMath::normalize(v);
  n = LiteMath::normalize(n);
  float3 s, t = n;
  CoordinateSystemV2(n, &s, &t);
  const float3 wo = LiteMath::normalize(float3(dot(l, s), dot(l, t), dot(l, n)));
  const float3 wi = LiteMath::normalize(float3(dot(v, s), dot(v, t), dot(v, n)));
  //const float3 wm = normalize(wo + wi);

  if (wi.z * wo.z < 0.0f)
  {
    return;
  }

  float3 half, diff;
  RusinkiewiczTransform(l, v, &half, &diff);

  float x[6];
  x[0] = half.x;
  x[1] = half.y;
  x[2] = half.z;
  x[3] = diff.x;
  x[4] = diff.y;
  x[5] = diff.z;

  float buf0[KANBRDF_MAX_SIZE];
  float buf1[KANBRDF_MAX_SIZE];

  
  evalKanLayer(weights + KANBRDF_WEIGTH_OFFSETS[0], x, buf0, KANBRDF_LAYER_SIZES[0], KANBRDF_LAYER_SIZES[1]);
  evalKanLayer(weights + KANBRDF_WEIGTH_OFFSETS[1], buf0, buf1, KANBRDF_LAYER_SIZES[1], KANBRDF_LAYER_SIZES[2]);
  evalKanLayer(weights + KANBRDF_WEIGTH_OFFSETS[1], buf1, buf0, KANBRDF_LAYER_SIZES[1], KANBRDF_LAYER_SIZES[2]);

  pRes->val = float4(kanbrdf_output_transform(buf0[0]), kanbrdf_output_transform(buf0[1]), kanbrdf_output_transform(buf0[2]), 1.0f);
  pRes->pdf = lambertEvalPDF(l, v, n); //TODO
}


static inline void kanBrdfSampleAndEval(const Material* a_materials, const float *weights, float4 rands, 
                                            float3 vec, float3 n, BsdfSample* pRes, int spectral_mode)
{
  const float3 lambertDir = MapSampleToCosineDistribution(rands.x, rands.y, vec, n, 1.0f);//lambertSample(float2(rands.x, rands.y), vec, n);
  const float  lambertPdf = lambertEvalPDF(lambertDir, vec, n);
  BsdfEval tRes;
  kanBrdfEval(a_materials, weights, lambertDir, vec, n, &tRes, spectral_mode);

  pRes->dir   = lambertDir;
  pRes->val   = tRes.val;
  pRes->pdf   = lambertPdf; //TODO
  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}


#endif