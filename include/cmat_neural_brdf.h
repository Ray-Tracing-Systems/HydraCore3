#ifndef INCLUDE__CMAT_NEURAL_BRDF_H_
#define INCLUDE__CMAT_NEURAL_BRDF_H_
#include "../neural.h" 

static constexpr uint NBRDF_LATENT_CHANNELS = 0;
static constexpr uint NBRDF_INPUT_DIM = NBRDF_LATENT_CHANNELS + 6;
static constexpr uint NBRDF_HIDDEN_DIM = 64;

static constexpr uint NBRDF_OFFSET_BIAS0 = NBRDF_INPUT_DIM * NBRDF_HIDDEN_DIM;
static constexpr uint NBRDF_OFFSET1      = NBRDF_OFFSET_BIAS0 + NBRDF_HIDDEN_DIM;
static constexpr uint NBRDF_SIZE_MAT     = NBRDF_HIDDEN_DIM * NBRDF_HIDDEN_DIM;
static constexpr uint NBRDF_SIZE_LAYER   = NBRDF_SIZE_MAT + NBRDF_HIDDEN_DIM;

static constexpr uint NBRDF_SIZE_MAT6    = NBRDF_HIDDEN_DIM * 6;


static inline void neuralBrdfSampleAndEval(const Material* a_materials, const float *weights, 
                                            float3 vec, float *tex, BsdfSample* pRes)
{



}



static inline void neuralBrdfEval(const Material* a_materials, const float *weights,
                                    float3 v0, float3 v1, float *tex, BsdfEval *pRes)
{
    float x[NBRDF_HIDDEN_DIM];
    x[0] = v0.x;
    x[1] = v0.y;
    x[2] = v0.z;
    x[3] = v1.x;
    x[4] = v1.y;
    x[5] = v1.z;
    for(uint i = 0; i < NBRDF_LATENT_CHANNELS; ++i) {
        x[6 + i] = tex[i];
    }

    float buf[NBRDF_HIDDEN_DIM];


    //Layer0
    nn::Matmul(weights,
               x, buf, NBRDF_INPUT_DIM, NBRDF_HIDDEN_DIM);
    nn::Add(weights + NBRDF_OFFSET_BIAS0,
            buf, buf, NBRDF_HIDDEN_DIM);
    nn::ReLU(buf, buf, NBRDF_HIDDEN_DIM);


    uint32_t current_offset = NBRDF_OFFSET1;

    //Layer1
    nn::Matmul(weights + current_offset,
               buf, x, NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
    nn::Add(weights + current_offset + NBRDF_SIZE_MAT,
            x, x, NBRDF_HIDDEN_DIM);
    nn::ReLU(x, x, NBRDF_HIDDEN_DIM);
    current_offset += NBRDF_SIZE_LAYER;

    //Layer2
    nn::Matmul(weights + current_offset,
               x, buf, NBRDF_HIDDEN_DIM, NBRDF_HIDDEN_DIM);
    nn::Add(weights + current_offset + NBRDF_SIZE_MAT,
            buf, buf, NBRDF_HIDDEN_DIM);
    nn::ReLU(buf, buf, NBRDF_HIDDEN_DIM);
    current_offset += NBRDF_SIZE_LAYER;

    //Layer3
    nn::Matmul(weights + current_offset,
               buf, x, NBRDF_HIDDEN_DIM, 6);
    nn::Add(weights + current_offset + NBRDF_SIZE_MAT6,
            x, x, 6);
    nn::ReLU(x, x, 6);
    current_offset += NBRDF_SIZE_MAT6 + 6;

    //Layer4
    nn::Matmul(weights + current_offset,
               x, buf, 6, 3);
    nn::Add(weights + current_offset + 6 * 3,
            buf, buf, 3);

    pRes->val = float4(buf[0], buf[1], buf[2], 1.0f);
}

#endif