#include "integrator_pt.h"
#include "include/cglobals.h"
#include "LiteMath.h"
#include <iostream>
#include <spectral/internal/common/format.h>

float fmaf(float a, float b, float c) 
{
    return a * b + c;
}

int safe_int(int i, int size)
{
    return i >= size ? size - 1 : i;
}

float smoothstep2(float x) {
    return LiteMath::smoothstep(0.0f, 1.0f, LiteMath::smoothstep(0.0f, 1.0f, x));
}

float aid_to_alpha(int alpha_id, int size)
{
    return smoothstep2(float(alpha_id) / float(size - 1));
}

float4 sigmoid(float4 x)
{
    return 0.5f * x / sqrt(x * x + 1) + 0.5f;
}

float4 sigmoid_polynomial(float4 x, float3 coef)
{
    return sigmoid((coef[0] * x + coef[1]) * x + coef[2]);
}

float inv_smoothstep(float x) 
{
    return 0.5f - sinf(asin(fmaf(-2.0f, x, 1.0f)) / 3.0f);
}

void Integrator::Upsample(const float4 *in_color, const float4 *in_wavelenghts, float4 *out_spectrum)
{   
    constexpr uint _size = UPSAMPLER_LUT_SIZE;
    constexpr uint _step = UPSAMPLER_LUT_STEP;

    float4 color = in_color[0];
    uint amax = color[0] >= color[1] ? 0 : 1;
    amax = color[amax] >= color[2] ? amax : 2;
    const uint offset = amax * _size * _size * _size;

    int a = 0, b = 0;

    const float alphaf = color[amax];

    if(alphaf > EPSILON_32) {
        a = int(color[(amax + 1) % 3] / alphaf * 255.0f);
        b = int(color[(amax + 2) % 3] / alphaf * 255.0f);
    }

    const int a1_id = a / _step;
    const int a2_id = safe_int(a1_id + 1, _size);
    const int b1_id = b / _step;
    const int b2_id = safe_int(b1_id + 1, _size);

    const int a1 = safe_int(a1_id * _step, 256);
    const int a2 = a == 255 ? 256 : safe_int(a2_id * _step, 256);
    const int b1 = safe_int(b1_id * _step, 256);
    const int b2 = b == 255 ? 256 : safe_int(b2_id * _step, 256);

    const unsigned alpha1_id = int((_size - 1) * inv_smoothstep(inv_smoothstep(alphaf))); 
    const unsigned alpha2_id = alpha1_id == _size - 1 ? _size - 1 : alpha1_id + 1;
    const float alphaf1 = aid_to_alpha(alpha1_id, _size);
    const float alphaf2 = alpha1_id == alpha2_id ? 2.0f : aid_to_alpha(alpha2_id, _size);

    const float daf = float(a2 - a1) / 255.0f;
    const float daf1 = float(a - a1) / 255.0f;
    const float daf2 = float(a2 - a) / 255.0f;

    const float dbf = float(b2 - b1) / 255.0f;
    const float dbf1 = float(b - b1) / 255.0f;
    const float dbf2 = float(b2 - b) / 255.0f; 

    const float dalphaf = alphaf2 - alphaf1;
    const float dalphaf1 = alphaf - alphaf1;
    const float dalphaf2 = alphaf2 - alphaf;

    const float div = 1.0f / (daf * dbf * dalphaf);
    
    float3 res = m_spec_lut[offset + ((alpha1_id * _size) + a1_id) * _size + b1_id] * daf2 * dbf2 * dalphaf2 * div
               + m_spec_lut[offset + ((alpha2_id * _size) + a1_id) * _size + b1_id] * daf2 * dbf2 * dalphaf1 * div
               + m_spec_lut[offset + ((alpha1_id * _size) + a1_id) * _size + b2_id] * daf2 * dbf1 * dalphaf2 * div
               + m_spec_lut[offset + ((alpha2_id * _size) + a1_id) * _size + b2_id] * daf2 * dbf1 * dalphaf1 * div
               + m_spec_lut[offset + ((alpha1_id * _size) + a2_id) * _size + b1_id] * daf1 * dbf2 * dalphaf2 * div
               + m_spec_lut[offset + ((alpha2_id * _size) + a2_id) * _size + b1_id] * daf1 * dbf2 * dalphaf1 * div
               + m_spec_lut[offset + ((alpha1_id * _size) + a2_id) * _size + b2_id] * daf1 * dbf1 * dalphaf2 * div
               + m_spec_lut[offset + ((alpha2_id * _size) + a2_id) * _size + b2_id] * daf1 * dbf1 * dalphaf1 * div;


    std::cerr << spec::format("[%f %f %f]", res.x, res.y, res.z) << std::endl; 
    out_spectrum[0] = sigmoid_polynomial(in_wavelenghts[0], res);

}