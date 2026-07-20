#include "include/cglobals.h"
#include "integrator_pt.h"

#include "include/cmaterial.h"
#include <cstdint>


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


//-----------------------NEURALBRDF----------------------
static constexpr uint32_t NBRDF_WEIGTH_OFFSETS[] = {0, 448, 4608, 8768, 12928};


static inline float4 invLogMapping(float4 x, float p_ref)
{
  float4 eX = float4(exp(x.x), exp(x.y), exp(x.z), exp(x.w));
  return eX * (p_ref + NBRDF_INVMAP_EPS) - NBRDF_INVMAP_EPS; 
}


static inline void NeuralBrdfLinear(const float *weights, const float x[NBRDF_MAX_SIZE], float res[NBRDF_MAX_SIZE], 
                                    uint32_t in_dim, uint32_t out_dim)
{
  for(uint32_t j = 0; j < out_dim; ++j) {
    res[j] = 0;
    for(uint32_t p = 0; p < in_dim; ++p) {
      res[j] += weights[j * in_dim + p] * x[p];
    }
    float bias = weights[in_dim * out_dim + j];
    res[j] += bias;
  }
}

static inline float4 NeuralBrdfCroppedLinear(const float *weights, const float x[NBRDF_MAX_SIZE], 
                                                   uint4 idx, uint32_t in_dim, uint32_t uncropped_out_dim)
{
  float4 res = float4(0, 0, 0, 0);
  for(uint32_t j = 0; j < 4; ++j) {
    uint cidx = idx[j];

    for(uint32_t p = 0; p < in_dim; ++p) {
      res[j] += weights[cidx * in_dim + p] * x[p];
    }
    float bias = weights[in_dim * uncropped_out_dim + cidx];
    res[j] += bias;
  }
  return res;
}

static inline void NeuralBrdfReLU(float x[NBRDF_MAX_SIZE], uint32_t size)
{
  for(uint32_t i = 0; i < size; ++i) {
    x[i] = max(x[i], 0.0f);
  }
}

inline uint NeuralBrdfBinarySearch(const float array[NBRDF_SPECTRUM_SIZE], float val) 
{
  int last  = int(NBRDF_SPECTRUM_SIZE) - 2;
  int first = 1;
  while (last > 0) 
  {
    uint half = uint(last) >> 1; 
    int middle = first + int(half);
    bool predResult = array[middle] <= val;
    first = predResult ? int(middle + 1) : first;
    last = predResult ? last - int(half + 1) : int(half);
  }
  return uint32_t(clamp(int(first - 1), 0, int(NBRDF_SPECTRUM_SIZE - 2)));
}


float4 Integrator::NeuralBrdfEvalInternal(uint weights_offset, uint median_entry_id, float4 wavelengths, float3 wi, float3 wo, int spectral_mode)
{
  uint out_dim = spectral_mode == 0 ? 3 : NBRDF_SPECTRUM_SIZE;

  float3 half, diff;
  RusinkiewiczTransform(wi, wo, &half, &diff);

  float buf0[NBRDF_MAX_SIZE];
  buf0[0] = half.x;
  buf0[1] = half.y;
  buf0[2] = half.z;
  buf0[3] = diff.x;
  buf0[4] = diff.y;
  buf0[5] = diff.z;

  float buf1[NBRDF_MAX_SIZE];


  uint32_t offset;

  //Layer0
  offset = weights_offset + NBRDF_WEIGTH_OFFSETS[0];
  NeuralBrdfLinear(m_neural_weights.data() + offset, buf0, buf1,
                   NBRDF_INPUT_DIM, NBRDF_HIDDEN_DIM);
  NeuralBrdfReLU(buf1, NBRDF_HIDDEN_DIM);

  //Layer1
  offset = weights_offset + NBRDF_WEIGTH_OFFSETS[1];
  NeuralBrdfLinear(m_neural_weights.data() + offset, buf1, buf0,
                   NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
  NeuralBrdfReLU(buf0, NBRDF_HIDDEN_DIM);

  //Layer2
  offset = weights_offset + NBRDF_WEIGTH_OFFSETS[2];
  NeuralBrdfLinear(m_neural_weights.data() + offset, buf0, buf1,
                   NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
  NeuralBrdfReLU(buf1, NBRDF_HIDDEN_DIM);

  //Layer3
  offset = weights_offset + NBRDF_WEIGTH_OFFSETS[3];
  NeuralBrdfLinear(m_neural_weights.data() + offset, buf1, buf0,
                   NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
  NeuralBrdfReLU(buf0, NBRDF_HIDDEN_DIM);


  float4 t;
  uint4 lerpIdx;
  for(int i = 0; i < 4; ++i) {
    float lambda = wavelengths[i];
    uint idx = NeuralBrdfBinarySearch(NBRDF_SPECTRAL_WAVELENGTHS, lambda);
    t[i] = (lambda - NBRDF_SPECTRAL_WAVELENGTHS[idx]) / (NBRDF_SPECTRAL_WAVELENGTHS[idx + 1] - NBRDF_SPECTRAL_WAVELENGTHS[idx]);

    lerpIdx[i] = idx; 
  }


  //Layer4
  offset = weights_offset + NBRDF_WEIGTH_OFFSETS[4];
  float4 y0 = NeuralBrdfCroppedLinear(m_neural_weights.data() + offset, buf0, 
                                       lerpIdx, NBRDF_HIDDEN_DIM, out_dim);
  y0 = clamp(y0, -20.0f, 20.0f);

  float4 y1 = NeuralBrdfCroppedLinear(m_neural_weights.data() + offset, buf0, 
                                       lerpIdx + 1, NBRDF_HIDDEN_DIM, out_dim);
  y1 = clamp(y1, -20.0f, 20.0f);


  float p_ref = MeasuredInterpIso1D(median_entry_id, wi, wo);

 // std::cout << "pref=" + std::to_string(p_ref) + " " << std::endl;

  y0 = invLogMapping(y0, p_ref);
  y1 = invLogMapping(y1, p_ref);

  float4 res = y0 + t * (y1 - y0);

 // std::cout << std::to_string(res[0]) + " " + std::to_string(res[1]) + " " + std::to_string(res[2]) + " " + std::to_string(res[3]) << std::endl; 

  return res;
}

void Integrator::NeuralBrdfEval(uint32_t matId, float4 wavelengths,
                                float3 l, float3 v, float3 n, BsdfEval *pRes, int spectral_mode)
{
  const float alpha0 = m_materials[matId].data[NBRDF_ALPHA];
  const float2 alpha = float2(alpha0, alpha0);
  const uint32_t median_idx = m_materials[matId].datai[NBRDF_MEDIANIDX];
  const uint32_t median_entry_id = m_materials[median_idx].datai[MEASURED_DATAIDX];
  const uint32_t weights_offset = m_neural_weights_offsets[matId];

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
  pRes->val = NeuralBrdfEvalInternal(weights_offset, median_entry_id, wavelengths, wo, wi, spectral_mode);
  wm        = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
  pRes->pdf = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
}


void Integrator::NeuralBrdfSampleAndEval(uint32_t matId, float4 rands, float4 wavelengths,
                                         float3 v, float3 n, BsdfSample* pRes, int spectral_mode)
{
  const float alpha0 = m_materials[matId].data[NBRDF_ALPHA];
  const float2 alpha = float2(alpha0, alpha0);
  const uint32_t median_idx = m_materials[matId].datai[NBRDF_MEDIANIDX];
  const uint32_t median_entry_id = m_materials[median_idx].datai[MEASURED_DATAIDX];
  const uint32_t weights_offset = m_neural_weights_offsets[matId];

  float3 nx, ny, nz = n;
  CoordinateSystemV2(nz, &nx, &ny);
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  if(wo.z == 0) return;

  float3 wm = trSample(wo, float2(rands.x, rands.y), alpha);
  float3 wi = reflect((-1.0f) * wo, wm);

  if(wo.z * wi.z < 0) return;// not in the same hemisphere


  pRes->val   = NeuralBrdfEvalInternal(weights_offset, median_entry_id, wavelengths, wo, wi, spectral_mode);
  pRes->pdf   = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  pRes->dir   = normalize(wi.x * nx + wi.y * ny + wi.z * nz);

  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}

//------------------------KANBRDF------------------------

static inline float kanbrdf_output_transform(float x) {
  float val = std::exp(x) - 1.0f;

  return max(0.0f, val);
}

static inline void KanBrdfNeuralLinear1(const float *weights, const float x[KANBRDF_MAX_SIZE * KANBRDF_GRID_SIZE], float res[KANBRDF_MAX_SIZE], 
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

static inline void KanBrdfNeuralLinear2(const float *weights, const float x[KANBRDF_MAX_SIZE], float res[KANBRDF_MAX_SIZE], 
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

static inline void KanBrdfNeuralAdd(float A[KANBRDF_MAX_SIZE], const float B[KANBRDF_MAX_SIZE], 
                    uint32_t m, uint32_t n)
{
  const uint32_t count = m * n; 

  for(uint32_t i = 0; i < count; ++i) {
    A[i] += B[i];
  }
}

static inline void KanBrdfNeuralSiLU(const float A[KANBRDF_MAX_SIZE], float res[KANBRDF_MAX_SIZE],
                    uint32_t m, uint32_t n)
{
  const uint32_t count = m * n; 

  for(uint32_t i = 0; i < count; ++i) {
    float v = A[i];
    res[i] = v * _sigmoid(v);
  }
}

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
  KanBrdfNeuralLinear1(m_neural_weights.data() + offset1, spline_basis, y, 1, KANBRDF_GRID_SIZE * in_dim, out_dim);


  float activated[KANBRDF_MAX_SIZE];
  float base[KANBRDF_MAX_SIZE];
  KanBrdfNeuralSiLU(x, activated, in_dim, 1);

  const uint offset2 = weights_offset + KANBRDF_GRID_SIZE * in_dim * out_dim;
  KanBrdfNeuralLinear2(m_neural_weights.data() + offset2, activated, base, 1, in_dim, out_dim);

  KanBrdfNeuralAdd(y, base, out_dim, 1);
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

void Integrator::KanBrdfEval(uint32_t matId, uint weights_offset, float3 l, float3 v, float3 n, BsdfEval *pRes, int spectral_mode)
{
  const float alpha0 = m_materials[matId].data[KANBRDF_ALPHA];
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


void Integrator::KanBrdfSampleAndEval(uint32_t matId, uint weights_offset, float4 rands, 
                                            float3 v, float3 n, BsdfSample* pRes, int spectral_mode)
{
  const float alpha0 = m_materials[matId].data[KANBRDF_ALPHA];
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

void Integrator::MeasuredEval(uint32_t matId, float3 l, float3 v, float3 n, BsdfEval *pRes, int spectral_mode)
{
  const uint32_t entry_id = m_materials[matId].datai[MEASURED_DATAIDX];

  const float alpha0 = m_materials[matId].data[MEASURED_ALPHA];
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

  if(spectral_mode == 0) {
    float3 t = MeasuredInterpRGB(entry_id, wi, wo);
    pRes->val = float4(t.x, t.y, t.z, 1.0);
  }
  else {
    pRes->val = float4(1.0f, 1.0f, 1.0f, 1.0f); //TODO
  }

  wm = normalize(wm);
  wm        = FaceForward(wm, float3(0.0f, 0.0f, 1.0f));
  pRes->pdf = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
}


void Integrator::MeasuredSampleAndEval(uint32_t matId, float4 rands, 
                                            float3 v, float3 n, BsdfSample* pRes, int spectral_mode)
{
  const uint32_t entry_id = m_materials[matId].datai[MEASURED_DATAIDX];

  const float alpha0 = m_materials[matId].data[MEASURED_ALPHA];
  const float2 alpha = float2(alpha0, alpha0);

  float3 nx, ny, nz = n;
  CoordinateSystemV2(nz, &nx, &ny);
  const float3 wo = float3(dot(v, nx), dot(v, ny), dot(v, nz));
  if(wo.z == 0) return;

  float3 wm = trSample(wo, float2(rands.x, rands.y), alpha);
  float3 wi = reflect((-1.0f) * wo, wm);

  if(wo.z * wi.z < 0) return;// not in the same hemisphere


  if(spectral_mode == 0) {
    float3 t = MeasuredInterpRGB(entry_id, wi, wo);
    pRes->val = float4(t.x, t.y, t.z, 1.0);
  }
  else {
    pRes->val = float4(1.0f, 1.0f, 1.0f, 1.0f); //TODO
  }

  pRes->pdf   = trPDF(wo, wm, alpha) / (4.0f * std::abs(dot(wo, wm)));
  pRes->dir   = normalize(wi.x * nx + wi.y * ny + wi.z * nz);

  pRes->flags = RAY_FLAG_HAS_NON_SPEC;
}