/*

  Some of the code in this file was taken directly from Mitsuba3 render.
  https://github.com/mitsuba-renderer/mitsuba3

  Mitsuba3 lincense is as follows:

  Copyright (c) 2017 Wenzel Jakob <wenzel.jakob@epfl.ch>, All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  You are under no obligation whatsoever to provide any bug fixes, patches, or
  upgrades to the features, functionality or performance of the source code
  ("Enhancements") to anyone; however, if you choose to make your Enhancements
  available either publicly, or directly to the author of this software, without
  imposing a separate written license agreement for such Enhancements, then you
  hereby grant the following license: a non-exclusive, royalty-free perpetual
  license to install, use, modify, prepare derivative works, incorporate into
  other computer software, distribute, and sublicense such enhancements or
  derivative works thereof, in binary and source code form.

*/

#include "mi_materials.h"
#include "include/cmaterial.h"
#include <array>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <tuple>

namespace dr
{
  static inline float fmadd(float a, float b, float c) { return a*b+c; }
  static inline float fmsub(float a, float b, float c) { return a*b-c; }

  template <size_t n>
  static inline float horner(float x, const std::array<float,n>& coeff) 
  {
    float accum = coeff[n - 1];
    for (size_t i = 1; i < n; ++i)
      accum = fmadd(x, accum, coeff[n - 1 - i]);
    return accum;
  }

  template <size_t n>
  static inline std::array<float, n> linspace(float start, float end, size_t points)
  {
    float delta = (end - start) / (points - 1);

    std::array<float, n> res;
    for(size_t i = 0; i < points; ++i)
    {
      res[i] = fmadd(i, delta, start);
    }
    
    return res;
  }

  template <typename T>
  std::vector<T> arange(T start, T stop, T step)
  {
    size_t size = size_t((stop - start + step - (step > 0 ? 1 : -1)) / step);

    std::vector<T> res(size);
    for(size_t i = 0; i < size; ++i)
    {
      res[i] = fmadd(i, step, start);
    }
    
    return res;
  }
};

namespace mi
{
  static constexpr uint32_t MI_ROUGH_TRANSMITTANCE_RES = 64;
  static constexpr double M_PI_64 = 3.14159265358979323846;
  static constexpr double EPSILON_64 = 0x1p-53;

  float fresnel_diffuse_reflectance(float eta) 
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


  /// Evaluate the l-th Legendre polynomial and its derivative using recurrence
  template <typename Value>
  std::pair<Value, Value> legendre_pd(int l, Value x) 
  {
    using Scalar = Value;

    assert(l >= 0);
    Value l_cur = Value(0), d_cur = Value(0);

    if (l > 1) 
    {
      Value l_p_pred = Value(1), l_pred = x,
            d_p_pred = Value(0), d_pred = Value(1);
      Value k0 = Value(3), k1 = Value(2), k2 = Value(1);

      for (int ki = 2; ki <= l; ++ki) 
      {
        l_cur = (k0 * x * l_pred - k2 * l_p_pred) / k1;
        d_cur = d_p_pred + k0 * l_pred;
        l_p_pred = l_pred; l_pred = l_cur;
        d_p_pred = d_pred; d_pred = d_cur;
        k2 = k1; k0 += Value(2); k1 += Value(1);
      }
    } 
    else 
    {
      if (l == 0) 
      {
        l_cur = Value(1); d_cur = Value(0);
      } 
      else 
      {
        l_cur = x; d_cur = Value(1);
      }
    }

    return { l_cur, d_cur };
  }

  std::pair<std::vector<float>, std::vector<float>> gauss_legendre(int n) 
  {
    if (n < 1)
      throw std::runtime_error("gauss_legendre(): n must be >= 1");

    std::vector<float> nodes(n), weights(n);

    n--;

    if (n == 0) 
    {
      nodes[0] = 0.0f;
      weights[0] = 2.0f;
    } 
    else if (n == 1) 
    {
      nodes[0] = -std::sqrt(1.0f / 3.0f);
      nodes[1] = -nodes[0];
      weights[0] = weights[1] = 1.0f;
    }

    int m = (n + 1) / 2;
    for (int i = 0; i < m; ++i) 
    {
      // Initial guess for this root using that of a Chebyshev polynomial
      double x = -std::cos((double) (2*i + 1) / (double) (2*n + 2) * M_PI_64);
      int it = 0;

      while (true) 
      {
        if (++it > 20)
          throw std::runtime_error("gauss_lobatto(" + std::to_string(n) + "): did not converge after 20 iterations!");

        // Search for the interior roots of P_{n+1}(x) using Newton's method.
        std::pair<double, double> L = legendre_pd(n+1, x);
        double step = L.first / L.second;
        x -= step;

        if (std::abs(step) <= 4 * std::abs(x) * EPSILON_64)
            break;
      }

      std::pair<double, double> L = legendre_pd(n+1, x);
      weights[i] = weights[n - i] =
          (float)(2 / ((1 - x * x) * (L.second * L.second)));
      nodes[i] = (float) x;
      nodes[n - i] = (float) -x;
      assert(i == 0 || (float) x > nodes[i-1]);
    }

    if ((n % 2) == 0) 
    {
        std::pair<double, double> L = legendre_pd(n+1, 0.0);
        weights[n / 2] = (float) (2.0 / (L.second * L.second));
        nodes[n / 2] = 0.0f;
    }

    return {nodes, weights};
  }


  float sin_theta_2(const float3 &v) 
  { 
    return dr::fmadd(v.x, v.x, v.y * v.y); 
  }

  /** \brief Give a unit direction, this function returns the sine and cosine
     * of the azimuth in a reference spherical coordinate system 
     */
  std::pair<float, float> sincos_phi(const float3 &v) 
  {
    float sin_theta2 = sin_theta_2(v);
    float inv_sin_theta = 1.f / std::sqrt(sin_theta_2(v));

    float2 result = {v.x * inv_sin_theta, v.y * inv_sin_theta};

    result = std::abs(sin_theta2) <= 4.f * EPSILON ? float2(1.f, 0.f) : clamp(result, -1.f, 1.f);

    return { result.y, result.x };
  }

  /// Low-distortion concentric square to disk mapping by Peter Shirley
  float2 square_to_uniform_disk_concentric(const float2 &sample) 
  {
    float x = dr::fmsub(2.f, sample.x, 1.f),
          y = dr::fmsub(2.f, sample.y, 1.f);

    float phi, r;
    if (x == 0 && y == 0) 
    {
      r = phi = 0;
    } 
    else if (x * x > y * y) 
    {
      r = x;
      phi = (M_PI / 4.f) * (y / x);
    } 
    else 
    {
      r = y;
      phi = (M_PI / 2.f) - (x / y) * (M_PI / 4.f);
    }
    return {r * std::cos(phi), r * std::sin(phi)};
  }

  /**
     * \brief Smith's shadowing-masking function for a single direction
     *
     * \param v
     *     An arbitrary direction
     * \param m
     *     The microfacet normal
  */
  float smith_g1(const float3 &v, const float3 &m, float2 alpha) 
  {
    float xy_alpha_2 = alpha.x * v.x * alpha.x * v.x + alpha.y * v.y * alpha.y * v.y,
          tan_theta_alpha_2 = xy_alpha_2 / (v.z * v.z),
          result;


    result = 2.f / (1.f + std::sqrt(1.f + tan_theta_alpha_2));

    // Perpendicular incidence -- no shadowing/masking
    if(xy_alpha_2 == 0.f)
      result = 1.f;

    /* Ensure consistent orientation (can't see the back
        of the microfacet from the front and vice versa) */

    if(v.z * dot(v, m) <= 0.f)
      result = 0.f;

    return result;
  }

  float2 sample_visible_11(float cos_theta_i, float2 sample)
  {
    float2 p = square_to_uniform_disk_concentric(sample);

    float s = 0.5f * (1.f + cos_theta_i);
    p.y = lerp(std::sqrt(1.f - p.x * p.x), p.y, s);

    // Project onto chosen side of the hemisphere
    float x = p.x, y = p.y,
          z = std::sqrt(1.f - dot(p, p));

    // Convert to slope
    float sin_theta_i = std::sqrt(1.f - cos_theta_i * cos_theta_i);
    float norm = 1.f / (dr::fmadd(sin_theta_i, y, cos_theta_i * z));
    return float2(dr::fmsub(cos_theta_i, y, sin_theta_i * z), x) * norm;
  }

  float3 sample_visible_normal(float3 wi, float2 rands, float2 alpha)
  {
    float sin_phi, cos_phi, cos_theta;

    // Step 1: stretch wi
    float3 wi_p = normalize(float3(alpha.x * wi.x, alpha.y * wi.y, wi.z));

    std::tie(sin_phi, cos_phi) = sincos_phi(wi_p);
    cos_theta = wi_p.z;

    // Step 2: simulate P22_{wi}(slope.x, slope.y, 1, 1)
    float2 slope = sample_visible_11(cos_theta, rands);

    // Step 3: rotate & unstretch
    slope = float2(
        dr::fmsub(cos_phi, slope.x, sin_phi * slope.y) * alpha.x,
        dr::fmadd(sin_phi, slope.x, cos_phi * slope.y) * alpha.y);

    // Step 4: compute normal
    float3 m = normalize(float3(-slope.x, -slope.y, 1));


    return m;
  }

  float3 refract(const float3 &wi, const float3 &m, float cos_theta_t, float eta_ti) 
  {
    return m * dr::fmadd(dot(wi, m), eta_ti, cos_theta_t) - wi * eta_ti;
  }

  float3 reflect(const float3 &wi, const float3 &m) 
  {
    return m * 2.f * dot(wi, m) - wi;
  }

  template <size_t n, size_t M>
  std::array<float, n> eval_transmittance(float alpha, float eta, const std::array<float, M>& mu, const std::array<float, M>& one_minus_mu_sqr,
                                          const std::array<float, M>& zeros)
  {
    size_t grid_sz = eta > 1.0f ? 32 : 128;

    auto [nodes, weights] = gauss_legendre(grid_sz); 

    std::array<float, n> result;
    for(size_t i = 0; i < M; ++i) 
    {
      float3 wi = {one_minus_mu_sqr[i], zeros[i], mu[i]};
      result[i] = 0.f;
      for(size_t j = 0; j < grid_sz * grid_sz; ++j)
      {
        size_t idx_x = j % grid_sz;
        size_t idx_y = j / grid_sz;
        float2 node   = {nodes[idx_x], nodes[idx_y]};
        float2 weight = {weights[idx_x], weights[idx_y]};

        node.x = node.x * 0.5f + 0.5f;
        node.y = node.y * 0.5f + 0.5f;

        // float3 normal_pbrt = trSample(wi, node, {alpha, alpha});
        float3 normal = sample_visible_normal(wi, node, {alpha, alpha});
        auto fres = FrDielectricDetailed(dot(wi, normal), eta);
        auto f = fres.x;
        auto cos_theta_t = fres.y;
        // auto eta_it = fres.z;
        auto eta_ti = fres.w;

        float3 wo   = mi::refract(wi, normal, cos_theta_t, eta_ti);
        // float smith = trD(wo, normal, {alpha, alpha}) * (1.f - f);
        float smith = smith_g1(wo, normal, {alpha, alpha}) * (1.f - f);;
        if (wo.z * wi.z >= 0.f)
          smith = 0.f;

        result[i] += smith * weight.x * weight.y * 0.25f;
        // std::cout << "i = " << i << " , j = " << j << " :" <<  result[i] << " : " 
        //           << "smith = " << smith << " , weight.x = " << weight.x << ", weight.y = " << weight.y<< std::endl;
      }
      
    }
    return result;
  }
  
  template <size_t n, size_t M>
  std::array<float, n> eval_reflectance(float alpha, float eta, const std::array<float, M>& mu, const std::array<float, M>& one_minus_mu_sqr,
                                        const std::array<float, M>& zeros)
  {
    size_t grid_sz = eta > 1.0f ? 32 : 128;

    auto [nodes, weights] = gauss_legendre(grid_sz); 

    std::array<float, n> result;
    for(size_t i = 0; i < M; ++i) 
    {
      float3 wi = {one_minus_mu_sqr[i], zeros[i], mu[i]};
      result[i] = 0.f;
      for(size_t j = 0; j < grid_sz * grid_sz; ++j)
      {
        size_t idx_x = j % grid_sz;
        size_t idx_y = j / grid_sz;
        float2 node   = {nodes[idx_x], nodes[idx_y]};
        float2 weight = {weights[idx_x], weights[idx_y]};

        node.x = node.x * 0.5f + 0.5f;
        node.y = node.y * 0.5f + 0.5f;

        float3 normal = sample_visible_normal(wi, node, {alpha, alpha});
        float3 wo     = mi::reflect(wi, normal);
        auto fres     = FrDielectricDetailed(dot(wi, normal), eta);
        auto f = fres.x;

        float smith = smith_g1(wo, normal, {alpha, alpha}) * f;
        if (wo.z <= 0.f)
          smith = 0.f;
        if (wi.z <= 0.f)
          smith = 0.f;

        result[i] += smith * weight.x * weight.y * 0.25f;
        // std::cout << "i = " << i << " , j = " << j << " :" <<  result[i] << " : " 
        //           << "smith = " << smith << " , weight = " << weight.x * weight.y
        //           << " , fresnel = " << f << std::endl;
      }
      std::cout << result[i] << std::endl;
    }
    return result;
  }

  void fresnel_coat_precompute(float alpha, float int_ior, float ext_ior, float4 diffuse_reflectance, float4 specular_reflectance,
                               bool is_spectral)
  {
    float eta = int_ior / ext_ior;
    float inv_eta_2 = 1.f / (eta * eta);

    uint32_t sz = is_spectral ? 4 : 3;

    float d_mean = 0.f, s_mean = 0.f;
    
    for(uint32_t i = 0; i < sz; ++i)
    {
      d_mean += diffuse_reflectance.M[i];
      s_mean += specular_reflectance.M[i];
    }
    d_mean /= sz;
    s_mean /= sz;

    float spec_sampling_weight = s_mean / (d_mean + s_mean);

    std::array<float, MI_ROUGH_TRANSMITTANCE_RES> zeros;
    zeros.fill(0.0f);

    std::array<float, MI_ROUGH_TRANSMITTANCE_RES> mu = dr::linspace<MI_ROUGH_TRANSMITTANCE_RES>(0.0f, 1.0f, MI_ROUGH_TRANSMITTANCE_RES);
    std::array<float, MI_ROUGH_TRANSMITTANCE_RES> one_minus_mu_sqr;
    for(size_t i = 0; i < MI_ROUGH_TRANSMITTANCE_RES; ++i)
    {
      mu[i] = std::max(mu[i], 1e-6f);
      one_minus_mu_sqr[i] = std::sqrt(1.0f - mu[i] * mu[i]);
    }

    auto transmittance = eval_transmittance<MI_ROUGH_TRANSMITTANCE_RES>(alpha, eta, mu, one_minus_mu_sqr, zeros);

    auto reflectance = eval_reflectance<MI_ROUGH_TRANSMITTANCE_RES>(alpha, 1.f / eta, mu, one_minus_mu_sqr, zeros);
    float internal_reflectance = 0.0f;
    for(auto i = 0; i < reflectance.size(); ++i)
    {
      internal_reflectance += reflectance[i] * mu[i];
    }
    internal_reflectance = (internal_reflectance / reflectance.size()) * 2.f;

    std::cout << "internal_transmittance = ";
    for(auto i = 0; i < transmittance.size(); ++i)
    {
      std::cout << transmittance[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "internal_reflectance = " << internal_reflectance << std::endl;
  }

};

void SetMiPlastic(Material* material, float int_ior, float ext_ior, float4 diffuse_reflectance, float4 specular_reflectance)
{ 
  material->colors[GLTF_COLOR_BASE] = diffuse_reflectance;
  material->colors[GLTF_COLOR_COAT] = specular_reflectance;
  
  const float m_eta = int_ior / ext_ior;
  material->data[GLTF_FLOAT_IOR]        = m_eta;  
  material->data[GLTF_FLOAT_MI_FDR_INT] = mi::fresnel_diffuse_reflectance(1.f / m_eta);
  material->data[GLTF_FLOAT_MI_FDR_EXT] = mi::fresnel_diffuse_reflectance(m_eta);
 
  const float d_mean = 0.3333333f * (diffuse_reflectance.x  + diffuse_reflectance.y  + diffuse_reflectance.z); 
  const float s_mean = 0.3333333f * (specular_reflectance.x + specular_reflectance.y + specular_reflectance.z); 

  material->data[GLTF_FLOAT_MI_SSW] = s_mean / (d_mean + s_mean);
}
