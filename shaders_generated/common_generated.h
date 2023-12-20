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

struct complex
{
  float re, im;
};

complex make_complex(float re, float im) { 
  complex res;
  res.re = re;
  res.im = im;
  return res;
}

complex to_complex(float re)              { return make_complex(re, 0.0f);}
complex complex_add(complex a, complex b) { return make_complex(a.re + b.re, a.im + b.im); }
complex complex_sub(complex a, complex b) { return make_complex(a.re - b.re, a.im - b.im); }
complex complex_mul(complex a, complex b) { return make_complex(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
complex complex_div(complex a, complex b) {
  const float scale = 1 / (b.re * b.re + b.im * b.im);
  return make_complex(scale * (a.re * b.re + a.im * b.im), scale * (a.im * b.re - a.re * b.im));
}

complex real_add_complex(float value, complex z) { return complex_add(to_complex(value),z); }
complex real_sub_complex(float value, complex z) { return complex_sub(to_complex(value),z); }
complex real_mul_complex(float value, complex z) { return complex_mul(to_complex(value),z); }
complex real_div_complex(float value, complex z) { return complex_div(to_complex(value),z); }

complex complex_add_real(complex z, float value) { return complex_add(z, to_complex(value)); }
complex complex_sub_real(complex z, float value) { return complex_sub(z, to_complex(value)); }
complex complex_mul_real(complex z, float value) { return complex_mul(z, to_complex(value)); }
complex complex_div_real(complex z, float value) { return complex_div(z, to_complex(value)); }

float real(complex z) { return z.re;}
float imag(complex z) { return z.im; }
float complex_norm(complex z) { return z.re * z.re + z.im * z.im; }
float complex_abs(complex z) { return sqrt(complex_norm(z)); }
complex complex_sqrt(complex z) 
{
  float n = complex_abs(z);
  float t1 = sqrt(0.5f * (n + abs(z.re)));
  float t2 = 0.5f * z.im / t1;
  if (n == 0.0f)
    return to_complex(0.0f);
  if (z.re >= 0.0f)
    return make_complex(t1, t2);
  else
    return make_complex(abs(t2), copysign(t1, z.im));
}

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
  vec3 tang;
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
const uint LIGHT_GEOM_RECT = 1;
const uint LIGHT_GEOM_DISC = 2;
const uint LIGHT_GEOM_SPHERE = 3;
const uint LIGHT_GEOM_DIRECT = 4;
const uint LIGHT_GEOM_POINT = 5;
const uint LIGHT_DIST_LAMBERT = 0;
const uint LIGHT_DIST_OMNI = 1;
struct LightSource
{
  mat4 matrix;    ///<! translation in matrix is always (0,0,0,1)
  vec4   pos;       ///<! translation aclually stored here
  vec4   intensity; ///<! brightress, i.e. screen value if light is visable directly
  vec4   norm;
  
  vec2   size;
  float    pdfA;
  uint     geomType;  ///<! LIGHT_GEOM_RECT, LIGHT_GEOM_DISC, LIGHT_GEOM_SPHERE, ...
  
  uint     distType;  ///<! LIGHT_DIST_LAMBERT, LIGHT_DIST_OMNI, ...
  uint     dummy1;
  uint     dummy2;
  uint     dummy3;

  uint     specId;
  uint     texId;
  uint     iesId;
  float    mult;
};
struct LightSample
{
  vec3 pos;
  vec3 norm;
  bool   isOmni;
};
struct BsdfSample
{
  vec4 val;
  vec3 dir;
  float  pdf; 
  uint   flags;
  float  ior;
};
struct BsdfEval
{
  vec4 val;
  float  pdf; 
};
const uint GLTF_COMPONENT_LAMBERT = 1;
const uint GLTF_COMPONENT_COAT = 2;
const uint GLTF_COMPONENT_METAL = 4;
const uint GLTF_METAL_PERF_MIRROR = 8;
const uint GLTF_COMPONENT_ORENNAYAR = 16;
const uint FLAG_NMAP_INVERT_X = 32;
const uint FLAG_NMAP_INVERT_Y = 64;
const uint FLAG_NMAP_SWAP_XY = 128;
const uint MAT_TYPE_GLTF = 1;
const uint MAT_TYPE_GLASS = 2;
const uint MAT_TYPE_CONDUCTOR = 3;
const uint MAT_TYPE_DIFFUSE = 4;
const uint MAT_TYPE_BLEND = 5;
const uint MAT_TYPE_LIGHT_SOURCE = 0xEFFFFFFF;
const uint RAY_EVENT_S = 1;
const uint RAY_EVENT_D = 2;
const uint RAY_EVENT_G = 4;
const uint RAY_EVENT_T = 8;
const uint RAY_EVENT_V = 16;
const uint RAY_EVENT_TOUT = 32;
const uint RAY_EVENT_TNINGLASS = 64;
const uint UINT_MTYPE = 0;
const uint UINT_CFLAGS = 1;
const uint UINT_LIGHTID = 2;
const uint UINT_NMAP_ID = 3;
const uint UINT_MAIN_LAST_IND = 4;
const uint GLTF_COLOR_BASE = 0;
const uint GLTF_COLOR_COAT = 1;
const uint GLTF_COLOR_METAL = 2;
const uint GLTF_COLOR_LAST_IND = GLTF_COLOR_METAL;
const uint GLTF_FLOAT_MI_FDR_INT = UINT_MAIN_LAST_IND + 0;
const uint GLTF_FLOAT_MI_FDR_EXT = UINT_MAIN_LAST_IND + 1;
const uint GLTF_FLOAT_MI_SSW = UINT_MAIN_LAST_IND + 2;
const uint GLTF_FLOAT_ALPHA = UINT_MAIN_LAST_IND + 3;
const uint GLTF_FLOAT_GLOSINESS = UINT_MAIN_LAST_IND + 4;
const uint GLTF_FLOAT_IOR = UINT_MAIN_LAST_IND + 5;
const uint GLTF_FLOAT_ROUGH_ORENNAYAR = UINT_MAIN_LAST_IND + 6;
const uint GLTF_UINT_TEXID0 = UINT_MAIN_LAST_IND + 7;
const uint GLTF_CUSTOM_LAST_IND = GLTF_UINT_TEXID0;
const uint GLASS_COLOR_REFLECT = 0;
const uint GLASS_COLOR_TRANSP = 1;
const uint GLASS_COLOR_LAST_IND = GLASS_COLOR_TRANSP;
const uint GLASS_FLOAT_GLOSS_REFLECT = UINT_MAIN_LAST_IND + 0;
const uint GLASS_FLOAT_GLOSS_TRANSP = UINT_MAIN_LAST_IND + 1;
const uint GLASS_FLOAT_IOR = UINT_MAIN_LAST_IND + 2;
const uint GLASS_CUSTOM_LAST_IND = GLASS_FLOAT_IOR;
const uint EMISSION_COLOR = 0;
const uint EMISSION_COLOR_LAST_IND = EMISSION_COLOR;
const uint EMISSION_MULT = UINT_MAIN_LAST_IND + 0;
const uint EMISSION_TEXID0 = UINT_MAIN_LAST_IND + 1;
const uint EMISSION_SPECID0 = UINT_MAIN_LAST_IND + 2;
const uint EMISSION_CUSTOM_LAST_IND = EMISSION_SPECID0;
const uint CONDUCTOR_COLOR = 0;
const uint CONDUCTOR_COLOR_LAST_IND = CONDUCTOR_COLOR;
const uint CONDUCTOR_ROUGH_U = UINT_MAIN_LAST_IND + 0;
const uint CONDUCTOR_ROUGH_V = UINT_MAIN_LAST_IND + 1;
const uint CONDUCTOR_ETA = UINT_MAIN_LAST_IND + 2;
const uint CONDUCTOR_K = UINT_MAIN_LAST_IND + 3;
const uint CONDUCTOR_TEXID0 = UINT_MAIN_LAST_IND + 4;
const uint CONDUCTOR_ETA_SPECID = UINT_MAIN_LAST_IND + 5;
const uint CONDUCTOR_K_SPECID = UINT_MAIN_LAST_IND + 6;
const uint CONDUCTOR_CUSTOM_LAST_IND = CONDUCTOR_K_SPECID;
const uint DIFFUSE_COLOR = 0;
const uint DIFFUSE_COLOR_LAST_IND = DIFFUSE_COLOR;
const uint DIFFUSE_ROUGHNESS = UINT_MAIN_LAST_IND + 0;
const uint DIFFUSE_TEXID0 = UINT_MAIN_LAST_IND + 1;
const uint DIFFUSE_SPECID = UINT_MAIN_LAST_IND + 2;
const uint DIFFUSE_CUSTOM_LAST_IND = DIFFUSE_SPECID;
const uint BLEND_COLOR_LAST_IND = 0;
const uint BLEND_WEIGHT = UINT_MAIN_LAST_IND + 0;
const uint BLEND_MAT_ID_1 = UINT_MAIN_LAST_IND + 1;
const uint BLEND_MAT_ID_2 = UINT_MAIN_LAST_IND + 2;
const uint BLEND_TEXID0 = UINT_MAIN_LAST_IND + 3;
const uint BLEND_CUSTOM_LAST_IND = BLEND_TEXID0;
const uint COLOR_DATA_SIZE = 3;
const uint CUSTOM_DATA_SIZE = 12;
struct Material
{
  vec4 colors[COLOR_DATA_SIZE]; ///< colors data

  vec4 row0[2];     ///< texture matrix
  vec4 row1[2];     ///< texture matrix
      
  float data[CUSTOM_DATA_SIZE]; ///< float, uint and custom data. Read uint: uint x = as_uint(data[INDEX]), write: data[INDEX] = as_float(x)
};
struct MatIdWeight
{
  uint id;
  float weight;
};
struct MatIdWeightPair
{
  MatIdWeight first;
  MatIdWeight second;
};
const uint BUILD_LOW = 0;
const uint BUILD_MEDIUM = 1;
const uint BUILD_HIGH = 2;
const uint BUILD_REFIT = 3;
struct CRT_Hit 
{
  float    t;         ///< intersection distance from ray origin to object
  uint primId; 
  uint instId;
  uint geomId;    ///< use 4 most significant bits for geometry type; thay are zero for triangles 
  float    coords[4]; ///< custom intersection data; for triangles coords[0] and coords[1] stores baricentric coords (u,v)
};
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
struct RefractResultT
{
  vec3 rayDir;
  bool   success;
  float  eta;

};
#define RefractResult RefractResultT
const uint INTEGRATOR_STUPID_PT = 0;
const uint INTEGRATOR_SHADOW_PT = 1;
const uint INTEGRATOR_MIS_PT = 2;
const uint TOTAL_FEATURES_NUM = 10;

#ifndef SKIP_UBO_INCLUDE
#include "include/Integrator_generated_ubo.h"
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

mat3 make_float3x3(vec3 a, vec3 b, vec3 c) { // different way than mat3(a,b,c)
  return mat3(a.x, b.x, c.x,
              a.y, b.y, c.y,
              a.z, b.z, c.z);
}

float Cos2Theta(vec3 w) {
  return w.z * w.z;
}

float Sin2Theta(vec3 w) {
  return max(0.0f, 1.0f - Cos2Theta(w));
}

float SinTheta(vec3 w) {
  return sqrt(Sin2Theta(w));
}

float CosPhi(vec3 w) {
  float sinTheta = SinTheta(w);
  return (sinTheta == 0) ? 1 : clamp(w.x / sinTheta, -1.0f, 1.0f);
}

float SinPhi(vec3 w) {
  float sinTheta = SinTheta(w);
  return (sinTheta == 0) ? 0 : clamp(w.y / sinTheta, -1.0f, 1.0f);
}

float Tan2Theta(vec3 w) {
  return Sin2Theta(w) / Cos2Theta(w);
}

float trLambda(vec3 w, vec2 alpha) {
  float tan2Theta = Tan2Theta(w);
  if (isinf(tan2Theta))
    return 0;
  float alpha2 = (CosPhi(w) * alpha.x) * (CosPhi(w) * alpha.x) + (SinPhi(w) * alpha.y) * (SinPhi(w) * alpha.y);
  return (sqrt(1.0f + alpha2 * tan2Theta) - 1.0f) / 2.0f;
}

float trD(vec3 wm, vec2 alpha) {
  float tan2Theta = Tan2Theta(wm);
  if (isinf(tan2Theta))
      return 0;
  float cos4Theta = Cos2Theta(wm) * Cos2Theta(wm);
  if (cos4Theta < 1e-16f)
      return 0;
  float e = tan2Theta * ((CosPhi(wm) / alpha.x) * (CosPhi(wm) / alpha.x) + (SinPhi(wm) / alpha.y) * (SinPhi(wm) / alpha.y));
  return 1.0f / (M_PI * alpha.x * alpha.y * cos4Theta * (1 + e) * (1 + e));
}

float AbsCosTheta(vec3 w) {
  return abs(w.z);
}

float trG1(vec3 w, vec2 alpha) { 
  return 1.0f / (1.0f + trLambda(w, alpha)); 
}

void CoordinateSystem(vec3 v1, inout vec3 v2, inout vec3 v3) {
  float invLen = 1.0f;

  if (abs(v1.x) > abs(v1.y))
  {
    invLen = 1.0f / sqrt(v1.x*v1.x + v1.z*v1.z);
    (v2)  = vec3((-1.0f) * v1.z * invLen,0.0f,v1.x * invLen);
  }
  else
  {
    invLen = 1.0f / sqrt(v1.y * v1.y + v1.z * v1.z);
    (v2)  = vec3(0.0f,v1.z * invLen,(-1.0f) * v1.y * invLen);
  }

  (v3) = cross(v1, (v2));
}

uint NextState(inout RandomGen gen) {
  const uint x = (gen.state).x * 17 + (gen.state).y * 13123;
  (gen.state).x = (x << 13) ^ x;
  (gen.state).y ^= (x << 7);
  return x;
}

float sinPhiPBRT(const vec3 w, const float sintheta) {
  if (sintheta == 0.0f)
    return 0.0f;
  else
    return clamp(w.y / sintheta, -1.0f, 1.0f);
}

float trD(vec3 w, vec3 wm, vec2 alpha) {
  return trG1(w, alpha) / AbsCosTheta(w) * trD(wm, alpha) * abs(dot(w, wm));
}

float FrComplexConductor(float cosThetaI, complex eta) {
  float sinThetaI = 1.0f - cosThetaI * cosThetaI;
  complex sinThetaT = real_div_complex(sinThetaI,(complex_mul(eta,eta)));
  complex cosThetaT = complex_sqrt(real_sub_complex(1.0f,sinThetaT));

  complex r_parl = complex_div((complex_sub(complex_mul(eta,to_complex(cosThetaI)),cosThetaT)),(complex_add(complex_mul(eta,to_complex(cosThetaI)),cosThetaT)));
  complex r_perp = complex_div((real_sub_complex(cosThetaI,complex_mul(eta,cosThetaT))),(real_add_complex(cosThetaI,complex_mul(eta,cosThetaT))));
  return (complex_norm(r_parl) + complex_norm(r_perp)) / 2.0f;
}

float trG(vec3 wo, vec3 wi, vec2 alpha) { 
  return 1.0f / (1.0f + trLambda(wo, alpha) + trLambda(wi, alpha)); 
}

vec2 SampleUniformDiskPolar(vec2 u) {
  float r = sqrt(u[0]);
  float theta = M_TWOPI * u[1];
  return vec2(r * cos(theta),r * sin(theta));
}

vec3 reflect2(const vec3 dir, const vec3 n) {  
  return normalize(dir - 2.0f * dot(dir, n) * n);  // dir - vector from light
}

float cosPhiPBRT(const vec3 w, const float sintheta) {
  if (sintheta == 0.0f)
    return 1.0f;
  else
    return clamp(w.x / sintheta, -1.0f, 1.0f);
}

float fresnelSlick(const float VdotH) {
  const float tmp = 1.0f - abs(VdotH);
  return (tmp*tmp)*(tmp*tmp)*tmp;
}

float GGX_GeomShadMask(const float cosThetaN, const float alpha) {
  // Height - Correlated G.
  //const float tanNV      = sqrt(1.0f - cosThetaN * cosThetaN) / cosThetaN;
  //const float a          = 1.0f / (alpha * tanNV);
  //const float lambda     = (-1.0f + sqrt(1.0f + 1.0f / (a*a))) / 2.0f;
  //const float G          = 1.0f / (1.0f + lambda);

  // Optimized and equal to the commented-out formulas on top.
  const float cosTheta_sqr = clamp(cosThetaN * cosThetaN, 0.0f, 1.0f);
  const float tan2         = (1.0f - cosTheta_sqr) / max(cosTheta_sqr, 1e-6f);
  const float GP           = 2.0f / (1.0f + sqrt(1.0f + alpha * alpha * tan2));
  return GP;
}

float GGX_Distribution(const float cosThetaNH, const float alpha) {
  const float alpha2 = alpha * alpha;
  const float NH_sqr = clamp(cosThetaNH * cosThetaNH, 0.0f, 1.0f);
  const float den    = NH_sqr * alpha2 + (1.0f - NH_sqr);
  return alpha2 / max(float((M_PI)) * den * den, 1e-6f);
}

vec3 MapSampleToCosineDistribution(float r1, float r2, vec3 direction, vec3 hit_norm, float power) {
  if(power >= 1e6f)
    return direction;

  const float sin_phi = sin(M_TWOPI * r1);
  const float cos_phi = cos(M_TWOPI * r1);

  //sincos(2.0f*r1*3.141592654f, &sin_phi, &cos_phi);

  const float cos_theta = pow(1.0f - r2, 1.0f / (power + 1.0f));
  const float sin_theta = sqrt(1.0f - cos_theta*cos_theta);

  vec3 deviation;
  deviation.x = sin_theta*cos_phi;
  deviation.y = sin_theta*sin_phi;
  deviation.z = cos_theta;

  vec3 ny = direction,  nx,  nz;
  CoordinateSystem(ny, nx, nz);

  {
    vec3 temp = ny;
    ny = nz;
    nz = temp;
  }

  vec3 res = nx*deviation.x + ny*deviation.y + nz*deviation.z;

  float invSign = dot(direction, hit_norm) > 0.0f ? 1.0f : -1.0f;

  if (invSign*dot(res, hit_norm) < 0.0f) // reflected ray is below surface #CHECK_THIS
  {
    res = (-1.0f)*nx*deviation.x + ny*deviation.y - nz*deviation.z;
    //belowSurface = true;
  }

  return res;
}

vec3 SphericalDirectionPBRT(const float sintheta, const float costheta, const float phi) { 
  return vec3(sintheta * cos(phi),sintheta * sin(phi),costheta); 
}

vec3 ggxSample(const vec2 rands, const vec3 v, const vec3 n, const float roughness) {
  const float roughSqr = roughness * roughness;
    
  vec3 nx,  ny, nz = n;
  CoordinateSystem(nz, nx, ny);
    
  const vec3 wo = vec3(dot(v, nx),dot(v, ny),dot(v, nz));
  const float phi       = rands.x * M_TWOPI;
  const float cosTheta  = clamp(sqrt((1.0f - rands.y) / (1.0f + roughSqr * roughSqr * rands.y - rands.y)), 0.0f, 1.0f);
  const float sinTheta  = sqrt(1.0f - cosTheta * cosTheta);
  const vec3 wh = SphericalDirectionPBRT(sinTheta, cosTheta, phi);
    
  const vec3 wi = 2.0f * dot(wo, wh) * wh - wo;      // Compute incident direction by reflecting about wm  
  return normalize(wi.x * nx + wi.y * ny + wi.z * nz); // back to normal coordinate system
}

float ggxEvalBSDF(const vec3 l, const vec3 v, const vec3 n, const float roughness) {
  if(abs(dot(l, n)) < 1e-5f)
    return 0.0f; 
 
  const float dotNV = dot(n, v);  
  const float dotNL = dot(n, l);
  if (dotNV < 1e-6f || dotNL < 1e-6f)
    return 0.0f; 

  const float  roughSqr = roughness * roughness;
  const vec3 h = normalize(v + l); // half vector.
  const float dotNH = dot(n, h);
  const float D     = GGX_Distribution(dotNH, roughSqr);
  const float G     = GGX_GeomShadMask(dotNV, roughSqr)*GGX_GeomShadMask(dotNL, roughSqr);      

  return (D * G / max(4.0f * dotNV * dotNL, 1e-6f));  // Pass single-scattering
}

float conductorRoughEvalInternal(vec3 wo, vec3 wi, vec3 wm, vec2 alpha, complex ior) {
  if(wo.z * wi.z < 0) // not in the same hemisphere
  {
    return 0.0f;
  }

  float cosTheta_o = AbsCosTheta(wo);
  float cosTheta_i = AbsCosTheta(wi);
  if (cosTheta_i == 0 || cosTheta_o == 0)
    return 0.0f;

  float F = FrComplexConductor(abs(dot(wo, wm)), ior);
  float val = trD(wm, alpha) * F * trG(wo, wi, alpha) / (4.0f * cosTheta_i * cosTheta_o);

  return val;
}

float epsilonOfPos(vec3 hitPos) { return max(max(abs(hitPos.x), max(abs(hitPos.y), abs(hitPos.z))), 2.0f*GEPSILON)*GEPSILON; }

float fresnel2(vec3 v, vec3 n, float ior) {
  // Calculating the angle of incidence of light
  const float cosi = dot(v, n);
  // We calculate the angle of refraction of light according to the Snellius law
  const float sint = sqrt(1.0f - cosi * cosi) / ior;
  // Check if there is a complete internal reflection
  if (sint > 1.0f) 
  {
    // If yes, then we return the reflection coefficient equal to 1
    return 1.0f;
  }
  else 
  {
    // Otherwise we calculate the angle of refraction of light
    const float cost = sqrt(1.0f - sint * sint);
    // We calculate the reflection coefficients for parallel and perpendicular polarization using Fresnel formulas
    const float Rp   = (ior * cosi - cost) / (ior * cosi + cost);
    const float Rs   = (cosi - ior * cost) / (cosi + ior * cost);
    // We return the average value of these coefficients
    return (Rp * Rp + Rs * Rs) * 0.5f;
  }
}

float rndFloat1_Pseudo(inout RandomGen gen) {
  const uint x = NextState(gen);
  const uint tmp = (x * (x * x * 15731 + 74323) + 871483);
  const float scale      = (1.0f / 4294967296.0f);
  return (float((tmp)))*scale;
}

vec4 hydraFresnelCond(vec4 f0, float VdotH, float ior, float roughness) {  
  if(ior == 0.0f) // fresnel reflactance is disabled
    return f0;

  return f0 + (vec4(1.0f) - f0) * fresnelSlick(VdotH); // return bsdf * (f0 + (1 - f0) * (1 - abs(VdotH))^5)
}

float trPDF(vec3 w, vec3 wm, vec2 alpha) { 
  return trD(w, wm, alpha); 
}

vec3 FaceForward(vec3 v, vec3 n2) {
    return (dot(v, n2) < 0.f) ? (-1.0f) * v : v;
}

MatIdWeightPair make_weight_pair(MatIdWeight a, MatIdWeight b) {
  MatIdWeightPair res;
  res.first  = a;
  res.second = b;
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

vec3 refract2(const vec3 dir, const vec3 n, const float relativeIor) {  
  const float cosi = dot(dir, n);        // dir - vector from light. The normal should always look at the light vector.
  const float eta  = 1.0f / relativeIor; // Since the incoming vector and the normal are directed in the same direction.
  const float k    = 1.0f - eta * eta * (1.0f - cosi * cosi);
  if (k < 0)       
    return reflect2(dir, n); // full internal reflection 
  else         
    return normalize(eta * dir - (eta * cosi + sqrt(k)) * n); // the refracted vector    
}

vec3 NormalMapTransform(const uint materialFlags, vec3 normalFromTex) {
  vec3 normalTS = vec3(2.0f * normalFromTex.x - 1.0f, 2.0f * normalFromTex.y - 1.0f, normalFromTex.z);

  if((materialFlags & FLAG_NMAP_INVERT_X) != 0)
    normalTS.x *= (-1.0f);

  if((materialFlags & FLAG_NMAP_INVERT_Y) != 0)
    normalTS.y *= (-1.0f);

  if((materialFlags & FLAG_NMAP_SWAP_XY) != 0)
  {
    float temp = normalTS.x;
    normalTS.x = normalTS.y;
    normalTS.y = temp;
  }

  return normalTS; // normalize(normalTS); // do we nedd this normalize here?
}

vec3 trSample(vec3 wo, vec2 rands, vec2 alpha) {
  // Transform _w_ to hemispherical configuration
  vec3 wh = normalize(vec3(alpha.x * wo.x,alpha.y * wo.y,wo.z));
  if (wh.z < 0)
  {
    wh = (-1.0f) * wh;
  }

  // Find orthonormal basis for visible normal sampling
  vec3 T1 = (wh.z < 0.99999f) ? normalize(cross(vec3(0,0,1), wh)) : vec3(1,0,0);
  vec3 T2 = cross(wh, T1);

  // Generate uniformly distributed points on the unit disk
  vec2 p = SampleUniformDiskPolar(rands);

  // Warp hemispherical projection for visible normal sampling
  float h = sqrt(1 - p.x * p.x);
  p.y = mix(h, p.y, (1 + wh.z) / 2);

  // Reproject to hemisphere and transform normal to ellipsoid configuration
  float pz = sqrt(max(0.0f, 1.0f - dot(p, p)));
  vec3 nh = p.x * T1 + p.y * T2 + pz * wh;
  return normalize(vec3(alpha.x * nh.x,alpha.y * nh.y,max(1e-6f, nh.z)));
}

float ggxEvalPDF(const vec3 l, const vec3 v, const vec3 n, const float roughness) { 
  const float dotNV = dot(n, v);
  const float dotNL = dot(n, l);
  if (dotNV < 1e-6f || dotNL < 1e-6f)
    return 1.0f;

  const float  roughSqr  = roughness * roughness;
    
  const vec3 h = normalize(v + l); // half vector.
  const float dotNH = dot(n, h);
  const float dotHV = dot(h, v);
  const float D     = GGX_Distribution(dotNH, roughSqr);
  return  D * dotNH / (4.0f * max(dotHV, 1e-6f));
}

float lambertEvalPDF(vec3 l, vec3 v, vec3 n) { 
  return abs(dot(l, n)) * INV_PI;
}

float FrDielectricPBRT(float cosThetaI, float etaI, float etaT) {
  cosThetaI = clamp(cosThetaI, -1.0f, 1.0f);
  // Potentially swap indices of refraction
  bool entering = cosThetaI > 0.0f;
  if (!entering) 
  {
    const float tmp = etaI;
    etaI = etaT;
    etaT = tmp;
    cosThetaI = abs(cosThetaI);
  }

  // Compute _cosThetaT_ using Snell's law
  float sinThetaI = sqrt(max(0.0f, 1.0f - cosThetaI * cosThetaI));
  float sinThetaT = etaI / etaT * sinThetaI;

  // Handle total internal reflection
  if (sinThetaT >= 1.0f) 
    return 1.0f;

  const float cosThetaT = sqrt(max(0.0f, 1.0f - sinThetaT * sinThetaT));
  const float Rparl     = ((etaT * cosThetaI) - (etaI * cosThetaT)) / ((etaT * cosThetaI) + (etaI * cosThetaT));
  const float Rperp     = ((etaI * cosThetaI) - (etaT * cosThetaT)) / ((etaI * cosThetaI) + (etaT * cosThetaT));
  return 0.5f*(Rparl * Rparl + Rperp * Rperp);
}

float orennayarFunc(const vec3 a_l, const vec3 a_v, const vec3 a_n, const float a_roughness) {
  const float cosTheta_wi = dot(a_l, a_n);
  const float cosTheta_wo = dot(a_v, a_n);

  const float sinTheta_wi = sqrt(max(0.0f, 1.0f - cosTheta_wi * cosTheta_wi));
  const float sinTheta_wo = sqrt(max(0.0f, 1.0f - cosTheta_wo * cosTheta_wo));

  const float sigma  = a_roughness * M_PI * 0.5f; //Radians(sig)
  const float sigma2 = sigma * sigma;
  const float A      = 1.0f - (sigma2 / (2.0f * (sigma2 + 0.33f)));
  const float B      = 0.45f * sigma2 / (sigma2 + 0.09f);

  ///////////////////////////////////////////////////////////////////////////// to PBRT coordinate system
  // wo = a_v = -ray_dir
  // wi = a_l = newDir
  //
  vec3 nx,  ny, nz = a_n;
  CoordinateSystem(nz, nx, ny);

  ///////////////////////////////////////////////////////////////////////////// to PBRT coordinate system

  // Compute cosine term of Oren-Nayar model
  float maxcos = 0.0f;

  if (sinTheta_wi > 1e-4f && sinTheta_wo > 1e-4f)
  {
    const vec3 wo = vec3(-dot(a_v, nx),-dot(a_v, ny),-dot(a_v, nz));
    const vec3 wi = vec3(-dot(a_l, nx),-dot(a_l, ny),-dot(a_l, nz));
    const float sinphii = sinPhiPBRT(wi, sinTheta_wi);
    const float cosphii = cosPhiPBRT(wi, sinTheta_wi);
    const float sinphio = sinPhiPBRT(wo, sinTheta_wo);
    const float cosphio = cosPhiPBRT(wo, sinTheta_wo);
    const float dcos    = cosphii * cosphio + sinphii * sinphio;
    maxcos              = max(0.0f, dcos);
  }

  // Compute sine and tangent terms of Oren-Nayar model
  float sinalpha = 0.0f, tanbeta = 0.0f;

  if (abs(cosTheta_wi) > abs(cosTheta_wo))
  {
    sinalpha = sinTheta_wo;
    tanbeta  = sinTheta_wi / max(abs(cosTheta_wi), DEPSILON);
  }
  else
  {
    sinalpha = sinTheta_wi;
    tanbeta  = sinTheta_wo / max(abs(cosTheta_wo), DEPSILON);
  }

  return (A + B * maxcos * sinalpha * tanbeta);
}

float lambertEvalBSDF(vec3 l, vec3 v, vec3 n) {
  return INV_PI;
}

vec3 lambertSample(const vec2 rands, const vec3 v, const vec3 n) {
  return MapSampleToCosineDistribution(rands.x, rands.y, n, n, 1.0f);
}

vec2 mulRows2x4(const vec4 row0, const vec4 row1, vec2 v) {
  vec2 res;
  res.x = row0.x*v.x + row0.y*v.y + row0.w;
  res.y = row1.x*v.x + row1.y*v.y + row1.w;
  return res;
}

float PdfAtoW(const float aPdfA, const float aDist, const float aCosThere) {
  return (aPdfA*aDist*aDist) / max(aCosThere, 1e-30f);
}

float misHeuristicPower1(float p) { return isfinite(p) ? abs(p) : 0.0f; }

vec4 rndFloat4_Pseudo(inout RandomGen gen) {
  uint x = NextState(gen);

  const uint x1 = (x * (x * x * 15731 + 74323) + 871483);
  const uint y1 = (x * (x * x * 13734 + 37828) + 234234);
  const uint z1 = (x * (x * x * 11687 + 26461) + 137589);
  const uint w1 = (x * (x * x * 15707 + 789221) + 1376312589);

  const float scale = (1.0f / 4294967296.0f);

  return vec4(float((x1)), float((y1)), float((z1)), float((w1)))*scale;
}

float SpectrumAverage(vec4 spec) {
  float sum = spec[0];
  for (uint i = 1; i < SPECTRUM_SAMPLE_SZ; ++i)
    sum += spec[int(i)];
  return sum / float(SPECTRUM_SAMPLE_SZ);
}

bool trEffectivelySmooth(vec2 alpha) { 
  return max(alpha.x, alpha.y) < 1e-3f; 
}

MatIdWeight make_id_weight(uint a, float b) {
  MatIdWeight res;
  res.id  = a;
  res.weight = b;
  return res;
}

vec3 EyeRayDirNormalized(float x, float y, mat4 a_mViewProjInv) {
  vec4 pos = vec4(2.0f*x - 1.0f,2.0f*y - 1.0f,0.0f,1.0f);
  pos = a_mViewProjInv * pos;
  pos /= pos.w;
  return normalize(pos.xyz);
}

vec3 XYZToRGB(vec3 xyz) {
  vec3 rgb;
  rgb[0] = +3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
  rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
  rgb[2] = +0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];

  return rgb;
}

float misWeightHeuristic(float a, float b) {
  const float w = misHeuristicPower1(a) / max(misHeuristicPower1(a) + misHeuristicPower1(b), 1e-30f);
  return isfinite(w) ? w : 0.0f;
}

void transform_ray3f(mat4 a_mWorldViewInv, inout vec3 ray_pos, inout vec3 ray_dir) {
  vec3 pos = mul4x3(a_mWorldViewInv, (ray_pos));
  vec3 pos2 = mul4x3(a_mWorldViewInv, ((ray_pos) + 100.0f*(ray_dir)));

  vec3 diff = pos2 - pos;

  (ray_pos)  = pos;
  (ray_dir)  = normalize(diff);
}

vec3 OffsRayPos(const vec3 a_hitPos, const vec3 a_surfaceNorm, const vec3 a_sampleDir) {
  const float signOfNormal2 = dot(a_sampleDir, a_surfaceNorm) < 0.0f ? -1.0f : 1.0f;
  const float offsetEps     = epsilonOfPos(a_hitPos);
  return a_hitPos + signOfNormal2*offsetEps*a_surfaceNorm;
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

vec2 rndFloat2_Pseudo(inout RandomGen gen) {
  uint x = NextState(gen); 

  const uint x1 = (x * (x * x * 15731 + 74323) + 871483);
  const uint y1 = (x * (x * x * 13734 + 37828) + 234234);

  const float scale     = (1.0f / 4294967296.0f);

  return vec2(float((x1)), float((y1)))*scale;
}

MisData makeInitialMisData() {
  MisData data;
  data.matSamplePdf = 1.0f;
  data.ior          = 1.0f; // start from air
  return data;
}

float maxcomp(vec3 v) { return max(v.x, max(v.y, v.z)); }

uint RealColorToUint32_f3(vec3 real_color) {
  float  r = real_color.x*255.0f;
  float  g = real_color.y*255.0f;
  float  b = real_color.z*255.0f;
  uint red = uint(r), green = uint(g), blue = uint(b);
  return red | (green << 8) | (blue << 16) | 0xFF000000;
}

uint fakeOffset(uint x, uint y, uint pitch) { return y*pitch + x; }  // RTV pattern, for 2D threading

#define KGEN_FLAG_RETURN            1
#define KGEN_FLAG_BREAK             2
#define KGEN_FLAG_DONT_SET_EXIT     4
#define KGEN_FLAG_SET_EXIT_NEGATIVE 8
#define KGEN_REDUCTION_LAST_STEP    16
#define SPECTRUM_H 
#define BASIC_PROJ_LOGIC_H 
#define TEST_CLASS_H 
#define IMAGE2D_H 
#define RTC_RANDOM 
#define CFLOAT_GUARDIAN 
#define RTC_MATERIAL 
#define MAXFLOAT FLT_MAX

