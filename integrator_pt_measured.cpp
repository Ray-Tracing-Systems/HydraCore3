#include "include/cglobals.h"
#include "integrator_pt.h"

#include "include/cmaterial.h"
#include <cstdint>

using namespace LiteMath;

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

float3 Integrator::MeasuredInterpRGB(uint32_t entry_id, float3 wo, float3 wi)
{
  if(entry_id == uint32_t(-1)) {
    return float3(1.0f, 1.0f, 1.0f);
  }

  MeasuredBrdfEntry entry = m_measured_brdf_data[entry_id];

  if(entry.nchannels != 3) {
    return float3(1.0f, 1.0f, 1.0f);
  }

  uint4 idx0, idx1;
  float4 coeff0;
  bool aniso = entry.dim.z > 1;
  GetMeasuredInterpParams(entry.dim, wo, wi, &idx0, &idx1, &coeff0);
  float4 coeff1 = 1 - coeff0;

  uint offset = uint(entry.offset);

  const uint nch = 3;

  float3 res = float3(0.0f, 0.0f, 0.0f);
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

  return float3(res.x, res.y, res.z);
}


float Integrator::MeasuredInterpIso1D(uint32_t entry_id, float3 wo, float3 wi)
{
  if(entry_id == uint32_t(-1)) {
    return 0.0f;
  }
  MeasuredBrdfEntry entry = m_measured_brdf_data[entry_id];
  if(entry.nchannels != 1) {
    return 0.0f;
  }

  uint4 idx0, idx1;
  float4 coeff0;
  GetMeasuredInterpParams(entry.dim, wo, wi, &idx0, &idx1, &coeff0);
  float4 coeff1 = 1 - coeff0;

  uint offset = uint(entry.offset);

  float res = 0.0f;
  uint i;

  i = offset + CalcOffset4d(idx0.x, idx0.y, idx0.z, idx0.w, entry.dim, 1); //0000
  res += m_measured_brdfs[i] * coeff1.x * coeff1.y * coeff1.z * coeff1.w;

  i = offset + CalcOffset4d(idx0.x, idx0.y, idx0.z, idx1.w, entry.dim, 1); //0001
  res += m_measured_brdfs[i] * coeff1.x * coeff1.y * coeff1.z * coeff0.w;

  i = offset + CalcOffset4d(idx0.x, idx1.y, idx0.z, idx0.w, entry.dim, 1); //0100
  res += m_measured_brdfs[i] * coeff1.x * coeff0.y * coeff1.z * coeff1.w;

  i = offset + CalcOffset4d(idx0.x, idx1.y, idx0.z, idx1.w, entry.dim, 1); //0101
  res += m_measured_brdfs[i] * coeff1.x * coeff0.y * coeff1.z * coeff0.w;

  i = offset + CalcOffset4d(idx1.x, idx0.y, idx0.z, idx0.w, entry.dim, 1); //1000
  res += m_measured_brdfs[i] * coeff0.x * coeff1.y * coeff1.z * coeff1.w;

  i = offset + CalcOffset4d(idx1.x, idx0.y, idx0.z, idx1.w, entry.dim, 1); //1001
  res += m_measured_brdfs[i] * coeff0.x * coeff1.y * coeff1.z * coeff0.w;

  i = offset + CalcOffset4d(idx1.x, idx1.y, idx0.z, idx0.w, entry.dim, 1); //1100
  res += m_measured_brdfs[i] * coeff0.x * coeff0.y * coeff1.z * coeff1.w;

  i = offset + CalcOffset4d(idx1.x, idx1.y, idx0.z, idx1.w, entry.dim, 1); //1101
  res += m_measured_brdfs[i] * coeff0.x * coeff0.y * coeff1.z * coeff0.w;

  return res;
}