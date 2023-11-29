#ifndef BASIC_PROJ_LOGIC_H
#define BASIC_PROJ_LOGIC_H

#include "LiteMath.h"
#ifndef __OPENCL_VERSION__
using namespace LiteMath;
#endif

static constexpr uint RAY_FLAG_IS_DEAD        = 0x80000000;
static constexpr uint RAY_FLAG_OUT_OF_SCENE   = 0x40000000;
static constexpr uint RAY_FLAG_HIT_LIGHT      = 0x20000000;
static constexpr uint RAY_FLAG_HAS_NON_SPEC   = 0x10000000; // at least one bounce was non specular
static constexpr uint RAY_FLAG_HAS_INV_NORMAL = 0x08000000;
//static constexpr uint RAY_FLAG_DUMMY        = 0x04000000;
//static constexpr uint RAY_FLAG_DUMMY        = 0x02000000;
//static constexpr uint RAY_FLAG_DUMMY        = 0x01000000;


static constexpr float LAMBDA_MIN = 360.0f;
static constexpr float LAMBDA_MAX = 830.0f;

using float4 = float4;
static constexpr uint32_t SPECTRUM_SAMPLE_SZ = 4; //sizeof(float4) / sizeof(float); // srry, sizeof() evaluation not yet supported ... 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct Lite_HitT
{
  float t;
  int   primId; 
  int   instId;
  int   geomId;
} Lite_Hit;

typedef struct SurfaceHitT
{
  float3 pos;
  float3 norm;
  float2 uv;
}SurfaceHit;


static inline float3 EyeRayDirNormalized(float x, float y, float4x4 a_mViewProjInv)
{
  float4 pos = float4(2.0f*x - 1.0f, 2.0f*y - 1.0f, 0.0f, 1.0f );
  pos = a_mViewProjInv * pos;
  pos /= pos.w;
  return normalize(to_float3(pos));
}

static inline uint RealColorToUint32_f3(float3 real_color)
{
  float  r = real_color.x*255.0f;
  float  g = real_color.y*255.0f;
  float  b = real_color.z*255.0f;
  unsigned int red = (unsigned int)r, green = (unsigned int)g, blue = (unsigned int)b;
  return red | (green << 8) | (blue << 16) | 0xFF000000;
}

static inline uint RealColorToUint32(float4 real_color)
{
  float  r = real_color.x*255.0f;
  float  g = real_color.y*255.0f;
  float  b = real_color.z*255.0f;
  float  a = real_color.w*255.0f;

  unsigned int red   = (unsigned int)r;
  unsigned int green = (unsigned int)g;
  unsigned int blue  = (unsigned int)b;
  unsigned int alpha = (unsigned int)a;

  return red | (green << 8) | (blue << 16) | (alpha << 24);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline void CoordinateSystem(float3 v1, float3* v2, float3* v3)
{
  float invLen = 1.0f;

  if (std::abs(v1.x) > std::abs(v1.y))
  {
    invLen = 1.0f / std::sqrt(v1.x*v1.x + v1.z*v1.z);
    (*v2)  = float3((-1.0f) * v1.z * invLen, 0.0f, v1.x * invLen);
  }
  else
  {
    invLen = 1.0f / sqrt(v1.y * v1.y + v1.z * v1.z);
    (*v2)  = float3(0.0f, v1.z * invLen, (-1.0f) * v1.y * invLen);
  }

  (*v3) = cross(v1, (*v2));
}

//constexpr float M_PI     = 3.14159265358979323846f;
//constexpr float M_TWOPI  = 6.28318530717958647692f;
//constexpr float INV_PI   = 0.31830988618379067154f

constexpr float GEPSILON = 1e-5f ;
constexpr float DEPSILON = 1e-20f;

//enum THREAD_FLAGS { THREAD_IS_DEAD = 2147483648};

static inline float3 MapSampleToCosineDistribution(float r1, float r2, float3 direction, float3 hit_norm, float power)
{
  if(power >= 1e6f)
    return direction;

  const float sin_phi = std::sin(M_TWOPI * r1);
  const float cos_phi = std::cos(M_TWOPI * r1);

  //sincos(2.0f*r1*3.141592654f, &sin_phi, &cos_phi);

  const float cos_theta = std::pow(1.0f - r2, 1.0f / (power + 1.0f));
  const float sin_theta = std::sqrt(1.0f - cos_theta*cos_theta);

  float3 deviation;
  deviation.x = sin_theta*cos_phi;
  deviation.y = sin_theta*sin_phi;
  deviation.z = cos_theta;

  float3 ny = direction, nx, nz;
  CoordinateSystem(ny, &nx, &nz);

  {
    float3 temp = ny;
    ny = nz;
    nz = temp;
  }

  float3 res = nx*deviation.x + ny*deviation.y + nz*deviation.z;

  float invSign = dot(direction, hit_norm) > 0.0f ? 1.0f : -1.0f;

  if (invSign*dot(res, hit_norm) < 0.0f) // reflected ray is below surface #CHECK_THIS
  {
    res = (-1.0f)*nx*deviation.x + ny*deviation.y - nz*deviation.z;
    //belowSurface = true;
  }

  return res;
}

/**
\brief  transform float2 sample in rect [-1,1]x[-1,1] to disc centered at (0,0) with radius == 1. 
\param  xy - input sample in rect [-1,1]x[-1,1]
\return position in disc
*/
static inline float2 MapSamplesToDisc(float2 xy)
{
  float x = xy.x;
  float y = xy.y;

  float r = 0;
  float phi = 0;

  float2 res = xy;

  if (x>y && x>-y)
  {
    r = x;
    phi = 0.25f*3.141592654f*(y / x);
  }

  if (x < y && x > -y)
  {
    r = y;
    phi = 0.25f*3.141592654f*(2.0f - x / y);
  }

  if (x < y && x < -y)
  {
    r = -x;
    phi = 0.25f*3.141592654f*(4.0f + y / x);
  }

  if (x >y && x<-y)
  {
    r = -y;
    phi = 0.25f*3.141592654f*(6 - x / y);
  }

  //float sin_phi, cos_phi;
  //sincosf(phi, &sin_phi, &cos_phi);
  float sin_phi = std::sin(phi);
  float cos_phi = std::cos(phi);

  res.x = r*sin_phi;
  res.y = r*cos_phi;

  return res;
}

static inline float epsilonOfPos(float3 hitPos) { return std::max(std::max(std::abs(hitPos.x), std::max(std::abs(hitPos.y), std::abs(hitPos.z))), 2.0f*GEPSILON)*GEPSILON; }

/**
\brief offset reflected ray position by epsilon;
\param  a_hitPos      - world space position on surface
\param  a_surfaceNorm - surface normal at a_hitPos
\param  a_sampleDir   - ray direction in which we are going to trace reflected ray
\return offseted ray position
*/
static inline float3 OffsRayPos(const float3 a_hitPos, const float3 a_surfaceNorm, const float3 a_sampleDir)
{
  const float signOfNormal2 = dot(a_sampleDir, a_surfaceNorm) < 0.0f ? -1.0f : 1.0f;
  const float offsetEps     = epsilonOfPos(a_hitPos);
  return a_hitPos + signOfNormal2*offsetEps*a_surfaceNorm;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


static inline void transform_ray3f(float4x4 a_mWorldViewInv, float3* ray_pos, float3* ray_dir) 
{
  float3 pos  = mul4x3(a_mWorldViewInv, (*ray_pos));
  float3 pos2 = mul4x3(a_mWorldViewInv, ((*ray_pos) + 100.0f*(*ray_dir)));

  float3 diff = pos2 - pos;

  (*ray_pos)  = pos;
  (*ray_dir)  = normalize(diff);
}

static inline float PdfAtoW(const float aPdfA, const float aDist, const float aCosThere)
{
  return (aPdfA*aDist*aDist) / std::max(aCosThere, 1e-30f);
}

static inline float PdfWtoA(const float aPdfW, const float aDist, const float aCosThere)
{
  return aPdfW * std::abs(aCosThere) / std::max(aDist*aDist, 1e-30f);
}

static inline float maxcomp(float3 v) { return std::max(v.x, std::max(v.y, v.z)); }

static inline float misHeuristicPower1(float p) { return std::isfinite(p) ? std::abs(p) : 0.0f; }
static inline float misWeightHeuristic(float a, float b)
{
  const float w = misHeuristicPower1(a) / std::max(misHeuristicPower1(a) + misHeuristicPower1(b), 1e-30f);
  return std::isfinite(w) ? w : 0.0f;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
\brief This structure is used as transit to pass MIS-weights-important-data from previouce bounce to current (or from current to next).

*/
struct MisData
{
  float matSamplePdf; ///< previous angle pdf (pdfW) that were used for sampling material. if < 0, then material sample was pure specular 
  float cosTheta;     ///< previous dot(matSam.dir, hit.norm)
  float ior;          ///< previous ior
  float dummy;        ///< dummy for 4 float
};

static inline bool isSpecular(const MisData* data) { return (data->matSamplePdf < 0.0f); }

static inline MisData makeInitialMisData()
{
  MisData data;
  data.matSamplePdf = 1.0f;
  data.ior          = 1.0f; // start from air
  return data;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


static inline float2 mulRows2x4(const float4 row0, const float4 row1, float2 v)
{
  float2 res;
  res.x = row0.x*v.x + row0.y*v.y + row0.w;
  res.y = row1.x*v.x + row1.y*v.y + row1.w;
  return res;
}

static inline uint  packXY1616(uint x, uint y) { return (y << 16u) | (x & 0x0000FFFF); }
static inline uint2 unpackXY1616(uint packedIndex) 
{
  uint2 res; 
  res.x = (packedIndex & 0x0000FFFF);         
  res.y = (packedIndex & 0xFFFF0000) >> 16;   
  return res;
}


#endif
