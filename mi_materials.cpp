#include "mi_materials.h"

#include <array>

namespace dr
{
  static inline float fmadd(float a, float b, float c) { return a*b+c; }

  template <size_t n>
  static inline float horner(float x, const std::array<float,n>& coeff) 
  {
    float accum = coeff[n - 1];
    for (size_t i = 1; i < n; ++i)
      accum = fmadd(x, accum, coeff[n - 1 - i]);
    return accum;
  }
};

namespace mi
{
  static float fresnel_diffuse_reflectance(float eta) 
  {
    /* Fast mode: the following code approximates the diffuse Frensel reflectance
       for the eta<1 and eta>1 cases. An evaluation of the accuracy led to the
       following scheme, which cherry-picks fits from two papers where they are
       best. */
    float inv_eta = 1.0f/eta;

    /* Fit by Egan and Hilgeman (1973). Works reasonably well for
       "normal" IOR values (<2).
       Max rel. error in 1.0 - 1.5 : 0.1%
       Max rel. error in 1.5 - 2   : 0.6%
       Max rel. error in 2.0 - 5   : 9.5%
    */
    float approx_1 = dr::fmadd(0.0636f, inv_eta, dr::fmadd(eta, dr::fmadd(eta, -1.4399f, 0.7099f), 0.6681f));

    /* Fit by d'Eon and Irving (2011)

       Maintains a good accuracy even for unrealistic IOR values.

       Max rel. error in 1.0 - 2.0   : 0.1%
       Max rel. error in 2.0 - 10.0  : 0.2%  */
    float approx_2 = dr::horner(inv_eta, std::array<float,6>{0.919317f, -3.4793f, 6.75335f, -7.80989f, 4.98554f, -1.36881f});

    return eta < 1.f ? approx_1 : approx_2;
  }
};

void SetMiPlastic(GLTFMaterial* material, float int_ior, float ext_ior, float3 diffuse_reflectance, float3 specular_reflectance)
{
  material->brdfType  = BRDF_TYPE_GLTF;
  material->lightId   = uint(-1); 
  material->baseColor = to_float4(diffuse_reflectance,0);
  material->coatColor = to_float4(specular_reflectance,0);
  
  const float m_eta = int_ior / ext_ior;
  material->ior                = m_eta;
  material->data[MI_ETA]       = m_eta;
  material->data[MI_INV_ETA_2] = 1.f / (m_eta * m_eta);
  
  material->data[MI_FDR_INT] = mi::fresnel_diffuse_reflectance(1.f / m_eta);
  material->data[MI_FDR_EXT] = mi::fresnel_diffuse_reflectance(m_eta);
 
  const float d_mean = 0.3333333f*(diffuse_reflectance.x  + diffuse_reflectance.y  + diffuse_reflectance.z); 
  const float s_mean = 0.3333333f*(specular_reflectance.x + specular_reflectance.y + specular_reflectance.z); 

  material->data[MI_SSW] = s_mean / (d_mean + s_mean);
}
