/////////////////////////////////////////////////////////////////////
/////////////  Required  Shader Features ////////////////////////////
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
/////////////////// include files ///////////////////////////////////
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
/////////////////// declarations in class ///////////////////////////
/////////////////////////////////////////////////////////////////////
#ifndef uint32_t
#define uint32_t uint
#endif
#define FLT_MAX 1e37f
#define FLT_MIN -1e37f
#define FLT_EPSILON 1e-6f
#define DEG_TO_RAD  0.017453293f
#define unmasked
#define half  float16_t
#define half2 f16vec2
#define half3 f16vec3
#define half4 f16vec4
bool  isfinite(float x)            { return !isinf(x); }
float copysign(float mag, float s) { return abs(mag)*sign(s); }
struct CamParameters    ///<! add any parameter you like to this structure
{
  float fov;
  float aspect;
  float nearPlane;
  float farPlane;
  int   spectralMode;
};
struct RayPosAndW 
{
  float origin[3]; ///<! ray origin, x,y,z
  float wave;      ///<! wavelength
};
struct RayDirAndT 
{
  float direction[3]; ///<! normalized ray direction, x,y,z
  float time;         ///<! time in ... 
};
const uint RAY_FLAG_IS_DEAD = 0x80000000;
const uint RAY_FLAG_OUT_OF_SCENE = 0x40000000;
const uint RAY_FLAG_HIT_LIGHT = 0x20000000;
const uint RAY_FLAG_HAS_NON_SPEC = 0x10000000;
const uint RAY_FLAG_HAS_INV_NORMAL = 0x08000000;
const float LAMBDA_MIN = 360.0f;
const float LAMBDA_MAX = 830.0f;
const uint SPECTRUM_SAMPLE_SZ = 4;
struct Lite_HitT
{
  float t;
  int   primId; 
  int   instId;
  int   geomId;
};
#define Lite_Hit Lite_HitT
struct SurfaceHitT
{
  vec3 pos;
  vec3 norm;
  vec2 uv;
};
#define SurfaceHit SurfaceHitT
const float GEPSILON = 1e-5f;
const float DEPSILON = 1e-20f;
struct MisData
{
  float matSamplePdf; ///< previous angle pdf (pdfW) that were used for sampling material. if < 0, then material sample was pure specular 
  float cosTheta;     ///< previous dot(matSam.dir, hit.norm)
  float ior;          ///< previous ior
  float dummy;        ///< dummy for 4 float
};
struct RandomGenT
{
  uvec2 state;

};
#define RandomGen RandomGenT
struct LensElementInterface 
{
  float curvatureRadius;
  float thickness;
  float eta;
  float apertureRadius;
};
const float CAM_LAMBDA_MIN = 360.0f;
const float CAM_LAMBDA_MAX = 830.0f;

#ifndef SKIP_UBO_INCLUDE
#include "include/CamTableLens_tablelens_gpu_ubo.h"
#endif

/////////////////////////////////////////////////////////////////////
/////////////////// local functions /////////////////////////////////
/////////////////////////////////////////////////////////////////////

mat4 translate4x4(vec3 delta)
{
  return mat4(vec4(1.0, 0.0, 0.0, 0.0),
              vec4(0.0, 1.0, 0.0, 0.0),
              vec4(0.0, 0.0, 1.0, 0.0),
              vec4(delta, 1.0));
}

mat4 rotate4x4X(float phi)
{
  return mat4(vec4(1.0f, 0.0f,  0.0f,           0.0f),
              vec4(0.0f, +cos(phi),  +sin(phi), 0.0f),
              vec4(0.0f, -sin(phi),  +cos(phi), 0.0f),
              vec4(0.0f, 0.0f,       0.0f,      1.0f));
}

mat4 rotate4x4Y(float phi)
{
  return mat4(vec4(+cos(phi), 0.0f, -sin(phi), 0.0f),
              vec4(0.0f,      1.0f, 0.0f,      0.0f),
              vec4(+sin(phi), 0.0f, +cos(phi), 0.0f),
              vec4(0.0f,      0.0f, 0.0f,      1.0f));
}

mat4 rotate4x4Z(float phi)
{
  return mat4(vec4(+cos(phi), sin(phi), 0.0f, 0.0f),
              vec4(-sin(phi), cos(phi), 0.0f, 0.0f),
              vec4(0.0f,      0.0f,     1.0f, 0.0f),
              vec4(0.0f,      0.0f,     0.0f, 1.0f));
}

mat4 inverse4x4(mat4 m) { return inverse(m); }
vec3 mul4x3(mat4 m, vec3 v) { return (m*vec4(v, 1.0f)).xyz; }
vec3 mul3x3(mat4 m, vec3 v) { return (m*vec4(v, 0.0f)).xyz; }

bool Quadratic(float A, float B, float C, inout float t0, inout float t1) {
  // Find quadratic discriminant
  double discrim = double(B) * double(B) - 4. * double(A) * double(C);
  if (discrim < 0.) 
    return false;
  double rootDiscrim = sqrt(discrim);
  float floatRootDiscrim   = float(rootDiscrim);
  // Compute quadratic _t_ values
  float q;
  if (float(B) < 0)
      q = -.5 * (B - floatRootDiscrim);
  else
      q = -.5 * (B + floatRootDiscrim);
  t0 = q / A;
  t1 = C / q;
  if (float(t0) > float(t1)) 
  {
    // std::swap(*t0, *t1);
    float temp = t0;
    t0 = t1;
    t1 = temp;
  }
  return true;
}

vec3 faceforward(const vec3 n, const vec3 v) { return (dot(n, v) < 0.f) ? (-1.0f)*n : n; }

uint NextState(inout RandomGen gen) {
  const uint x = (gen.state).x * 17 + (gen.state).y * 13123;
  (gen.state).x = (x << 13) ^ x;
  (gen.state).y ^= (x << 7);
  return x;
}

bool Refract(const vec3 wi, const vec3 n, float eta, inout vec3 wt) {
  // Compute $\cos \theta_\roman{t}$ using Snell's law
  float cosThetaI  = dot(n, wi);
  float sin2ThetaI = max(float(0), float(1.0f - cosThetaI * cosThetaI));
  float sin2ThetaT = eta * eta * sin2ThetaI;
  // Handle total internal reflection for transmission
  if (sin2ThetaT >= 1) return false;
  float cosThetaT = sqrt(1 - sin2ThetaT);
  wt = eta * (-1.0f)*wi + (eta * cosThetaI - cosThetaT) * n;
  return true;
}

float SpectrumAverage(vec4 spec) {
  float sum = spec[0];
  for (uint i = 1; i < SPECTRUM_SAMPLE_SZ; ++i)
    sum += spec[int(i)];
  return sum / float(SPECTRUM_SAMPLE_SZ);
}

vec3 XYZToRGB(vec3 xyz) {
  vec3 rgb;
  rgb[0] = +3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
  rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
  rgb[2] = +0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];

  return rgb;
}

vec4 SampleWavelengths(float u, float a, float b) {
  // pdf is 1.0f / (b - a)
  vec4 res;

  res[0] = mix(a, b, u);

  float delta = (b - a) / float(SPECTRUM_SAMPLE_SZ);
  for (uint i = 1; i < SPECTRUM_SAMPLE_SZ; ++i) 
  {
      res[int(i)] = res[i - 1] + delta;
      if (res[int(i)] > b)
        res[int(i)] = a + (res[int(i)] - b);
  }

  return res;
}

vec2 MapSamplesToDisc(vec2 xy) {
  float x = xy.x;
  float y = xy.y;

  float r = 0;
  float phi = 0;

  vec2 res = xy;

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
  float sin_phi = sin(phi);
  float cos_phi = cos(phi);

  res.x = r*sin_phi;
  res.y = r*cos_phi;

  return res;
}

vec4 rndFloat4_Pseudo(inout RandomGen gen) {
  uint x = NextState(gen);

  const uint x1 = (x * (x * x * 15731 + 74323) + 871483);
  const uint y1 = (x * (x * x * 13734 + 37828) + 234234);
  const uint z1 = (x * (x * x * 11687 + 26461) + 137589);
  const uint w1 = (x * (x * x * 15707 + 789221) + 1376312589);

  const float scale = (1.0f / 4294967296.0f);

  return vec4(float((x1)), float((y1)), float((z1)), float((w1)))*scale;
}

#define KGEN_FLAG_RETURN            1
#define KGEN_FLAG_BREAK             2
#define KGEN_FLAG_DONT_SET_EXIT     4
#define KGEN_FLAG_SET_EXIT_NEGATIVE 8
#define KGEN_REDUCTION_LAST_STEP    16
#define SPECTRUM_H 
#define RTC_RANDOM 
#define BASIC_PROJ_LOGIC_H 
#define CFLOAT_GUARDIAN 
#define MAXFLOAT FLT_MAX

