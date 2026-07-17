#include "include/cglobals.h"
#include "integrator_pt.h"

#include "include/cmaterial.h"
#include <cstdint>

using namespace LiteMath;


static inline void NeuralLinear1(const float *weights, const float x[KANBRDF_MAX_SIZE * KANBRDF_GRID_SIZE], float res[KANBRDF_MAX_SIZE], 
                    uint32_t batch_size, uint32_t in_dim, uint32_t out_dim)
{
  for(uint32_t i = 0; i < batch_size; ++i) {
    for(uint32_t j = 0; j < out_dim; ++j) {
      res[i * out_dim + j] = 0;
      for(uint32_t p = 0; p < in_dim; ++p) {
        res[i * out_dim + j] += weights[j * in_dim + p] * x[i * in_dim + p]; // ???
      }
      // if(has_bias) {
      //  float bias = weights[weights_offset + in_dim * out_dim + j];
      //  res[i * out_dim + j] += bias;
      // }
    }
  }
}

static inline void NeuralLinear2(const float *weights, const float x[KANBRDF_MAX_SIZE], float res[KANBRDF_MAX_SIZE], 
                    uint32_t batch_size, uint32_t in_dim, uint32_t out_dim)
{
  for(uint32_t i = 0; i < batch_size; ++i) {
    for(uint32_t j = 0; j < out_dim; ++j) {
      res[i * out_dim + j] = 0;
      for(uint32_t p = 0; p < in_dim; ++p) {
        res[i * out_dim + j] += weights[j * in_dim + p] * x[i * in_dim + p]; // ???
      }
     // if(has_bias) {
        float bias = weights[in_dim * out_dim + j];
        res[i * out_dim + j] += bias;
     // }

    }
  }
}



static inline void NeuralAdd(float A[KANBRDF_MAX_SIZE], const float B[KANBRDF_MAX_SIZE], 
                    uint32_t m, uint32_t n)
{
  const uint32_t count = m * n; 

  for(uint32_t i = 0; i < count; ++i) {
    A[i] += B[i];
  }
}

static inline float _sigmoid(float x)
{
  if(x >= 0) {
    return 1.0f / (1.0f + std::exp(-x));
  }
  else {
    float e = std::exp(x);
    return e / (1.0f + e);
  }
}

/**
 * A   : m x n
 * out : m x n
 * 
 * (A == out) is possible
 */
static inline void NeuralSiLU(const float A[KANBRDF_MAX_SIZE], float res[KANBRDF_MAX_SIZE],
                    uint32_t m, uint32_t n)
{
  const uint32_t count = m * n; 

  for(uint32_t i = 0; i < count; ++i) {
    float v = A[i];
    res[i] = v * _sigmoid(v);
  }
}


static inline float kanbrdf_output_transform(float x) {
  float val = std::exp(x) - 1.0f;

  return max(0.0f, val);
}


//------------------------KANBRDF------------------------

void Integrator::EvalKANLayer(uint weights_offset, const float x[KANBRDF_MAX_SIZE], float y[KANBRDF_MAX_SIZE], uint in_dim, uint out_dim)
{
  const float grid_step = (KANBRDF_GRID_MAX - KANBRDF_GRID_MIN) / (KANBRDF_GRID_SIZE - 1);
  const float denom_inv = 1.0f / grid_step;

  float spline_basis[KANBRDF_GRID_SIZE * KANBRDF_MAX_SIZE];
  for(uint i = 0; i < in_dim; ++i) {
    for(uint j = 0; j < KANBRDF_GRID_SIZE; ++j) {
      float grid_val = KANBRDF_GRID_MIN + grid_step * float(j);
      float t = (x[i] - grid_val) * denom_inv;

      spline_basis[i * KANBRDF_GRID_SIZE + j] = std::exp(-t * t);
    }
  }

  const uint offset1 = weights_offset;
  NeuralLinear1(m_neural_weights.data() + offset1, spline_basis, y, 1, KANBRDF_GRID_SIZE * in_dim, out_dim);


  float activated[KANBRDF_MAX_SIZE];
  float base[KANBRDF_MAX_SIZE];
  NeuralSiLU(x, activated, in_dim, 1);

  const uint offset2 = weights_offset + KANBRDF_GRID_SIZE * in_dim * out_dim;
  NeuralLinear2(m_neural_weights.data() + offset2, activated, base, 1, in_dim, out_dim);

  NeuralAdd(y, base, out_dim, 1);
}


float4 Integrator::KanBrdfEvalInternal(uint weights_offset, float3 wo, float3 wi, int spectral_mode)
{
  float3 hf, diff;
  RusinkiewiczTransform(wi, wo, &hf, &diff);

  float x[KANBRDF_MAX_SIZE];
  x[0] = hf.x;
  x[1] = hf.y;
  x[2] = hf.z;
  x[3] = diff.x;
  x[4] = diff.y;
  x[5] = diff.z;

  float buf0[KANBRDF_MAX_SIZE];
  float buf1[KANBRDF_MAX_SIZE];

  uint32_t offset = weights_offset + KANBRDF_WEIGTH_OFFSETS[0];
  EvalKANLayer(offset, x, buf0, KANBRDF_LAYER_SIZES[0], KANBRDF_LAYER_SIZES[1]);

  offset += KANBRDF_WEIGTH_OFFSETS[1];
  EvalKANLayer(offset, buf0, buf1, KANBRDF_LAYER_SIZES[1], KANBRDF_LAYER_SIZES[2]);

  offset += KANBRDF_WEIGTH_OFFSETS[2];
  EvalKANLayer(offset, buf1, buf0, KANBRDF_LAYER_SIZES[2], KANBRDF_LAYER_SIZES[3]);

  return float4(kanbrdf_output_transform(buf0[0]), kanbrdf_output_transform(buf0[1]), kanbrdf_output_transform(buf0[2]), 1.0f);
}

void Integrator::KanBrdfEval(uint32_t a_matId, uint weights_offset, float3 l, float3 v, float3 n, BsdfEval *pRes, int spectral_mode)
{
  const float alpha0 = m_materials[a_matId].data[KANBRDF_ALPHA];
  const float2 alpha = float2(alpha0, alpha0);
  float3 nx, ny, nz = n;
  CoordinateSystemV2(nz, &nx, &ny);

  // v = (-1.0f) * v;
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  const float3 wi = float3(dot(l, nx), dot(l, ny), dot(l, nz));

  if(wo.z * wi.z < 0.0f)
    return;

  float3 wm = wo + wi;
  if (dot(wm, wm) == 0)
      return;

  wm = normalize(wm);
  pRes->val = KanBrdfEvalInternal(weights_offset, wo, wi, spectral_mode);
  wm        = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
  pRes->pdf = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
}


void Integrator::KanBrdfSampleAndEval(uint32_t a_matId, uint weights_offset, float4 rands, 
                                            float3 v, float3 n, BsdfSample* pRes, int spectral_mode)
{
  const float alpha0 = m_materials[a_matId].data[KANBRDF_ALPHA];
  const float2 alpha = float2(alpha0, alpha0);

  float3 nx, ny, nz = n;
  CoordinateSystemV2(nz, &nx, &ny);
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  if(wo.z == 0) return;

  float3 wm = trSample(wo, float2(rands.x, rands.y), alpha);
  float3 wi = reflect((-1.0f) * wo, wm);

  if(wo.z * wi.z < 0) return;// not in the same hemisphere


  pRes->val   = KanBrdfEvalInternal(weights_offset, wo, wi, spectral_mode);
  pRes->pdf   = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  pRes->dir   = normalize(wi.x * nx + wi.y * ny + wi.z * nz);

  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}

//------------------------MEASURED------------------------

static inline float4 RvectorsToRangles(float3 half, float3 diff)
{
  float4 res;
  res.x = atan2(sqrt(half.x * half.x + half.y * half.y), half.z); //theta_h
  res.y = atan2(sqrt(diff.x * diff.x + diff.y * diff.y), diff.z); //theta_d

  res.z = atan2(half.y, half.x); //phi_h
  res.w = atan2(diff.y, diff.x); //phi_d

  return res;
}


static inline void GetMeasuredInterpParams(uint4 dim, float3 wo, float3 wi, uint4 *idx0, uint4 *idx1, float4 *coeff)
{

  const float4 MAX_ANGLES = float4(0.5f * M_PI, 0.5f * M_PI, M_PI, M_PI);


  float4 maxIdx = float4(dim - 1);

  float3 half, diff;
  RusinkiewiczTransform(wi, wo, &half, &diff);
  float4 angles = RvectorsToRangles(half, diff);
  if(angles.z < 0) angles.z += M_PI;
  if(angles.w < 0) angles.w += M_PI;


  float4 idxF = clamp(angles / MAX_ANGLES, 0.0f, 1.0f);
  idxF.x = sqrt(idxF.x);
  idxF *= maxIdx;

  uint4 idxI0 = uint4(idxF);
  uint4 idxI1 = min(idxI0 + 1, dim - 1);

  float x00 = sqr(float(idxI0.x) / maxIdx.x) * MAX_ANGLES.x;
  float x01 = sqr(float(idxI1.x) / maxIdx.x) * MAX_ANGLES.x;
  float dx0_safe = x01 - x00;
  if(dx0_safe == 0.0f) dx0_safe = 1.0f;

  *coeff = idxF - float4(idxI0);
  coeff->x = (angles.x - x00) / dx0_safe;

  *idx0 = idxI0; 
  *idx1 = idxI1;
}

static inline uint32_t CalcOffset4d(uint i0, uint i1, uint i2, uint i3, uint4 dim, uint n_channels)
{
  return (((i0 * dim.y + i1) * dim.z + i2) * dim.w + i3) * n_channels;
}

float4 Integrator::MeasuredEvalInternal(uint32_t entry_id, float3 wo, float3 wi, int spectral_mode)
{
  if(entry_id == uint32_t(-1)) {
    return float4(1.0f, 1.0f, 1.0f, 1.0f);
  }


  MeasuredBrdfEntry entry = m_measured_brdf_data[entry_id];

  uint4 idx0, idx1;
  float4 coeff0;
  bool aniso = entry.dim.z > 1;
  GetMeasuredInterpParams(entry.dim, wo, wi, &idx0, &idx1, &coeff0);
  float4 coeff1 = 1 - coeff0;

  uint offset = uint(entry.offset);



  if(spectral_mode == 0) {
    uint nch = 3;

    float3 res = float3(0, 0, 0);
    uint i;

    i = offset + CalcOffset4d(idx0.x, idx0.y, idx0.z, idx0.w, entry.dim, nch); //0000
    res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
         * coeff1.x * coeff1.y * coeff1.z * coeff1.w;

    i = offset + CalcOffset4d(idx0.x, idx0.y, idx0.z, idx1.w, entry.dim, nch); //0001
    res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
         * coeff1.x * coeff1.y * coeff1.z * coeff0.w;

    i = offset + CalcOffset4d(idx0.x, idx1.y, idx0.z, idx0.w, entry.dim, nch); //0100
    res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
         * coeff1.x * coeff0.y * coeff1.z * coeff1.w;

    i = offset + CalcOffset4d(idx0.x, idx1.y, idx0.z, idx1.w, entry.dim, nch); //0101
    res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
         * coeff1.x * coeff0.y * coeff1.z * coeff0.w;

    i = offset + CalcOffset4d(idx1.x, idx0.y, idx0.z, idx0.w, entry.dim, nch); //1000
    res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
         * coeff0.x * coeff1.y * coeff1.z * coeff1.w;

    i = offset + CalcOffset4d(idx1.x, idx0.y, idx0.z, idx1.w, entry.dim, nch); //1001
    res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
         * coeff0.x * coeff1.y * coeff1.z * coeff0.w;

    i = offset + CalcOffset4d(idx1.x, idx1.y, idx0.z, idx0.w, entry.dim, nch); //1100
    res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
         * coeff0.x * coeff0.y * coeff1.z * coeff1.w;

    i = offset + CalcOffset4d(idx1.x, idx1.y, idx0.z, idx1.w, entry.dim, nch); //1101
    res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
         * coeff0.x * coeff0.y * coeff1.z * coeff0.w;


    if(aniso) {

      i = offset + CalcOffset4d(idx0.x, idx0.y, idx1.z, idx0.w, entry.dim, nch); //0010
      res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
           * coeff1.x * coeff1.y * coeff0.z * coeff1.w;

      i = offset + CalcOffset4d(idx0.x, idx0.y, idx1.z, idx1.w, entry.dim, nch); //0011
      res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
           * coeff1.x * coeff1.y * coeff0.z * coeff0.w;

      i = offset + CalcOffset4d(idx0.x, idx1.y, idx1.z, idx0.w, entry.dim, nch); //0110
      res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
           * coeff1.x * coeff0.y * coeff0.z * coeff1.w;

      i = offset + CalcOffset4d(idx0.x, idx1.y, idx1.z, idx1.w, entry.dim, nch); //0111
      res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
           * coeff1.x * coeff0.y * coeff0.z * coeff0.w;

      i = offset + CalcOffset4d(idx1.x, idx0.y, idx1.z, idx0.w, entry.dim, nch); //1010
      res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
           * coeff0.x * coeff1.y * coeff0.z * coeff1.w;

      i = offset + CalcOffset4d(idx1.x, idx0.y, idx1.z, idx1.w, entry.dim, nch); //1011
      res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
           * coeff0.x * coeff1.y * coeff0.z * coeff0.w;

      i = offset + CalcOffset4d(idx1.x, idx1.y, idx1.z, idx0.w, entry.dim, nch); //1110
      res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
           * coeff0.x * coeff0.y * coeff0.z * coeff1.w;

      i = offset + CalcOffset4d(idx1.x, idx1.y, idx1.z, idx1.w, entry.dim, nch); //1111
      res += float3(m_measured_brdfs[i], m_measured_brdfs[i + 1], m_measured_brdfs[i + 2])
           * coeff0.x * coeff0.y * coeff0.z * coeff0.w;
    }

    return float4(res.x, res.y, res.z, 1.0f);
  }
  


  return float4(1.0f, 1.0f, 1.0f, 1.0f);
}

void Integrator::MeasuredEval(uint32_t a_matId, float3 l, float3 v, float3 n, BsdfEval *pRes, int spectral_mode)
{
  const uint32_t entry_id = m_materials[a_matId].datai[MEASURED_DATAIDX];

  const float alpha0 = m_materials[a_matId].data[MEASURED_ALPHA];
  const float2 alpha = float2(alpha0, alpha0);
  float3 nx, ny, nz = n;
  CoordinateSystemV2(nz, &nx, &ny);

  // v = (-1.0f) * v;
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  const float3 wi = float3(dot(l, nx), dot(l, ny), dot(l, nz));

  if(wo.z * wi.z < 0.0f)
    return;

  float3 wm = wo + wi;
  if (dot(wm, wm) == 0)
      return;

  wm = normalize(wm);
  pRes->val = MeasuredEvalInternal(entry_id, wo, wi, spectral_mode);
  wm        = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
  pRes->pdf = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
}


void Integrator::MeasuredSampleAndEval(uint32_t a_matId, float4 rands, 
                                            float3 v, float3 n, BsdfSample* pRes, int spectral_mode)
{
  const uint32_t entry_id = m_materials[a_matId].datai[MEASURED_DATAIDX];

  const float alpha0 = m_materials[a_matId].data[MEASURED_ALPHA];
  const float2 alpha = float2(alpha0, alpha0);

  float3 nx, ny, nz = n;
  CoordinateSystemV2(nz, &nx, &ny);
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  if(wo.z == 0) return;

  float3 wm = trSample(wo, float2(rands.x, rands.y), alpha);
  float3 wi = reflect((-1.0f) * wo, wm);

  if(wo.z * wi.z < 0) return;// not in the same hemisphere


  pRes->val   = MeasuredEvalInternal(entry_id, wo, wi, spectral_mode);
  pRes->pdf   = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  pRes->dir   = normalize(wi.x * nx + wi.y * ny + wi.z * nz);

  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}