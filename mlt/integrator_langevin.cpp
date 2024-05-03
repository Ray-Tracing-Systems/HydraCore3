#include "utils.h" // for progress bar

#include <iomanip>
#include <omp.h>

#include "integrator_pt.h"

#include "include/cmaterial.h"
#include "include/cmat_gltf.h"
#include "include/cmat_conductor.h"
#include "include/cmat_glass.h"
#include "include/cmat_diffuse.h"
#include "include/cmat_plastic.h"
#include "include/cmat_dielectric.h"

class IntegratorLangevin : public Integrator {
public:
    IntegratorLangevin(int a_maxThreads, int a_spectral_mode, std::vector<uint32_t> a_features) : Integrator(a_maxThreads, a_spectral_mode, a_features) {}

    float GetRandomNumbersSpec(uint tid, RandomGen *a_gen, float *rands);
    float GetRandomNumbersTime(uint tid, RandomGen *a_gen, float *rands);
    float4 GetRandomNumbersLens(uint tid, RandomGen *a_gen, float *rands);
    float4 GetRandomNumbersMats(uint tid, RandomGen *a_gen, int a_bounce, float *rands);
    float4 GetRandomNumbersLgts(uint tid, RandomGen *a_gen, int a_bounce, float *rands);
    float GetRandomNumbersMatB(uint tid, RandomGen *a_gen, int a_bounce, int a_layer, float *rands);

    void kernel_InitEyeRay(
        uint tid,
        const uint *packedXY,
        float4 *rayPosAndNear,
        float4 *rayDirAndFar,
        float4 *wavelengths,
        float4 *accumColor,
        float4 *accumuThoroughput,
        RandomGen *gen,
        uint *rayFlags,
        MisData *misData,
        float *time,
        int *pX,
        int *pY,
        float *rands
    );

    void kernel_RayTrace(
        uint tid,
        uint bounce,
        const float4 *rayPosAndNear,
        const float4 *rayDirAndFar,
        const float *a_time,
        float4 *out_hit1,
        float4 *out_hit2,
        float4 *out_hit3,
        uint *out_instId,
        uint *rayFlags,
        float *rands
    );

    void kernel_SampleLightSource(
        uint tid,
        const float4 *rayPosAndNear,
        const float4 *rayDirAndFar,
        const float4 *wavelengths,
        const float4 *in_hitPart1,
        const float4 *in_hitPart2,
        const float4 *in_hitPart3,
        const uint *rayFlags,
        const float *a_time,
        uint bounce,
        RandomGen *a_gen,
        float4 *out_shadeColor,
        float *rands
    );

    void kernel_NextBounce(
        uint tid,
        uint bounce,
        const float4 *in_hitPart1,
        const float4 *in_hitPart2,
        const float4 *in_hitPart3,
        const uint *in_instId,
        const float4 *in_shadeColor,
        float4 *rayPosAndNear,
        float4 *rayDirAndFar,
        const float4 *wavelengths,
        float4 *accumColor,
        float4 *accumThoroughput,
        RandomGen *a_gen,
        MisData *misPrev,
        uint *rayFlags,
        float *rands
    );

    void kernel_HitEnvironment(
        uint tid,
        const uint *rayFlags,
        const float4 *rayDirAndFar,
        const MisData *a_prevMisData,
        const float4 *accumThoroughput,
        float4 *accumColor,
        float *rands
    );

    EyeRayData SampleCameraRay(
        RandomGen *pGen,
        uint tid,
        float *rands
    );

    BsdfSample MaterialSampleAndEval(
        uint a_materialId, uint tid, uint bounce, float4 wavelengths, 
        RandomGen *a_gen, float3 v, float3 n, float3 tan, float2 tc,
        MisData *a_misPrev, const uint a_currRayFlags, float *rands
    );

    uint32_t IntegratorLangevin::BlendSampleAndEval(
        uint a_materialId, uint tid, uint bounce, uint layer, float4 wavelengths, RandomGen* a_gen, float3 v, float3 n, float2 tc, 
        MisData* a_misPrev, BsdfSample* a_pRes, float *rands
    );

    void PathTraceBlock(uint pixelsNum, uint channels, float *out_color, uint a_passNum);

    float4 PathTraceF(
        uint tid,
        int *pX,
        int *pY,
        float *rands
    );

    static constexpr uint BOUNCE_START = 6;
    static constexpr uint LGHT_ID = 0;
    static constexpr uint MATS_ID = 4;
    static constexpr uint BLND_ID = 8;
    static constexpr uint PER_BOUNCE = 10;

    inline uint RandsPerThread() const {
        return PER_BOUNCE * m_traceDepth + BOUNCE_START;
    }

    std::vector<float> all_rands;
    uint m_randsPerThread = 0;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

constexpr float MUTATE_COEFF_SCREEN = 128.0f;
constexpr float MUTATE_COEFF_BSDF = 64.0f;

/**
\brief mutate random number in primary sample space -- interval [0,1]
\param valueX    - input original value
\param rands     - input pseudo random numbers
\param p2        - parameter of step size. The greater parameter is, the smaller step we gain. Default = 64.0f;
\param p1        - parameter of step size. The greater parameter is, the smaller step we gain. Default = 1024.0f;
\return mutated random float in range [0,1]

*/
static inline float MutateKelemen(float valueX, float2 rands, float p2, float p1) // mutate in primary space
{
    const float s1 = 1.0f / p1;
    const float s2 = 1.0f / p2;
    const float power = -std::log(s2 / s1);
    const float dv = std::max(s2 * (std::exp(power * std::sqrt(rands.x)) - std::exp(power)), 0.0f);

    if (rands.y < 0.5f) {
        valueX += dv;
        if (valueX > 1.0f)
            valueX -= 1.0f;
    } else {
        valueX -= dv;
        if (valueX < 0.0f)
            valueX += 1.0f;
    }

    return valueX;
}

float4 IntegratorLangevin::GetRandomNumbersLens(uint tid, RandomGen *a_gen, float *rands) {
    float *data = rands + tid * m_randsPerThread;
    return float4(data[0], data[1], data[2], data[3]);
}

float IntegratorLangevin::GetRandomNumbersSpec(uint tid, RandomGen *a_gen, float *rands) {
    float *data = rands + tid * m_randsPerThread;
    return data[4];
}

float IntegratorLangevin::GetRandomNumbersTime(uint tid, RandomGen *a_gen, float *rands) {
    float *data = rands + tid * m_randsPerThread;
    return data[5];
}

float4 IntegratorLangevin::GetRandomNumbersMats(uint tid, RandomGen *a_gen, int a_bounce, float *rands) {
    float *data = rands + tid * m_randsPerThread + BOUNCE_START + a_bounce * PER_BOUNCE + MATS_ID;
    return float4(data[0], data[1], data[2], data[3]);
}

float4 IntegratorLangevin::GetRandomNumbersLgts(uint tid, RandomGen *a_gen, int a_bounce, float *rands) {
    float *data = rands + tid * m_randsPerThread + BOUNCE_START + a_bounce * PER_BOUNCE + LGHT_ID;
    return float4(data[0], data[1], data[2], data[3]);
}

float IntegratorLangevin::GetRandomNumbersMatB(uint tid, RandomGen *a_gen, int a_bounce, int a_layer, float *rands) {
    float *data = rands + tid * m_randsPerThread + BOUNCE_START + a_bounce * PER_BOUNCE + BLND_ID;
    return data[a_layer];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void IntegratorLangevin::kernel_InitEyeRay(
    uint tid, const uint *packedXY, float4 *rayPosAndNear, float4 *rayDirAndFar,
    float4 *wavelengths, float4 *accumColor, float4 *accumuThoroughput, RandomGen *gen,
    uint *rayFlags, MisData *misData, float *time, int *pX, int *pY, float *rands) {
    if (tid >= m_maxThreadId)
        return;

    *accumColor = make_float4(0, 0, 0, 0);
    *accumuThoroughput = make_float4(1, 1, 1, 1);
    *rayFlags = 0;
    *misData = makeInitialMisData();

    RandomGen genLocal = m_randomGens[RandomGenId(tid)];

    EyeRayData r = SampleCameraRay(&genLocal, tid, rands);

    if (KSPEC_SPECTRAL_RENDERING != 0 && m_spectral_mode != 0)
        *wavelengths = SampleWavelengths(r.waveSam, LAMBDA_MIN, LAMBDA_MAX);
    else
        *wavelengths = float4(0.0f);

    *time = r.timeSam;

    transform_ray3f(m_worldViewInv, &r.rayPos, &r.rayDir);

    *rayPosAndNear = to_float4(r.rayPos, 0.0f);
    *rayDirAndFar = to_float4(r.rayDir, FLT_MAX);
    *gen = genLocal;

    (*pX) = int(r.x);
    (*pY) = int(r.y);
}

float4 IntegratorLangevin::PathTraceF(uint tid, int *pX, int *pY, float *rands) {
    float4 accumColor, accumThroughput;
    float4 rayPosAndNear, rayDirAndFar;
    float4 wavelengths;
    RandomGen gen;
    MisData mis;
    uint rayFlags;
    float time;

    kernel_InitEyeRay(
        tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar,
        &wavelengths, &accumColor, &accumThroughput, &gen, &rayFlags,
        &mis, &time, pX, pY, rands);

    for (uint depth = 0; depth < m_traceDepth; depth++) {
        float4 shadeColor, hitPart1, hitPart2, hitPart3;
        uint instId;
        kernel_RayTrace(
            tid, depth, &rayPosAndNear, &rayDirAndFar, &time,
            &hitPart1, &hitPart2, &hitPart3, &instId, &rayFlags, rands);
        if (isDeadRay(rayFlags))
            break;

        if (depth <= m_traceDepth - 1) {
            kernel_SampleLightSource(
                tid, &rayPosAndNear, &rayDirAndFar, &wavelengths, &hitPart1,
                &hitPart2, &hitPart3, &rayFlags, &time, depth, &gen, &shadeColor, rands);
        }

        kernel_NextBounce(
            tid, depth, &hitPart1, &hitPart2, &hitPart3, &instId, &shadeColor,
            &rayPosAndNear, &rayDirAndFar, &wavelengths, &accumColor, &accumThroughput,
            &gen, &mis, &rayFlags, rands);

        if (isDeadRay(rayFlags))
            break;
    }

    kernel_HitEnvironment(
        tid, &rayFlags, &rayDirAndFar, &mis,
        &accumThroughput, &accumColor, rands);

    return accumColor * m_exposureMult;
}

static inline float contribFunc(float4 color)
{
    return std::max(0.333334f * (color.x + color.y + color.z), 0.0f);
}

uint32_t AlignedSize(uint32_t a_size, uint32_t a_alignment);

bool SaveImage4fToBMP(const float *rgb, int width, int height, const char *outfilename, float a_normConst, float a_gamma);

void IntegratorLangevin::PathTraceBlock(uint pixelsNum, uint channels, float *out_color, uint a_passNum) {
    uint maxThreads = omp_get_max_threads();
    m_randsPerThread = AlignedSize(RandsPerThread(), uint32_t(16));
    m_maxThreadId = maxThreads;
    m_randomGens.resize(maxThreads);
    all_rands.resize(m_randsPerThread * maxThreads);
    const size_t samplesPerPass = (size_t(pixelsNum) * size_t(a_passNum)) / size_t(maxThreads);

    std::cout << "[IntegratorLangevin]: state size = " << m_randsPerThread << std::endl;

    ConsoleProgressBar progress(pixelsNum * a_passNum);
    progress.Start();
    auto start = std::chrono::high_resolution_clock::now();

    double avgBrightnessOut = 0.0f;
    float avgAcceptanceRate = 0.0f;
    const float plarge = 0.25f; // 25% of large step;

#ifndef _DEBUG
#pragma omp parallel default(shared)
#endif
    {
        int tid = omp_get_thread_num();

        RandomGen gen1 = RandomGenInit(tid * 7 + 1);
        RandomGen gen2 = RandomGenInit(tid);
        for (int i = 0; i < 10 + tid % 17; i++) {
            NextState(&gen1);
            NextState(&gen2);
        }

        // (1) Initial State

        int xScr = 0, yScr = 0;
        float *rands_new = all_rands.data() + tid * m_randsPerThread;
        std::vector<float> rands_cur(m_randsPerThread);
        for (size_t i = 0; i < rands_cur.size(); i++) {
            float rand = rndFloat1_Pseudo(&gen2);
            rands_cur[i] = rand;
            rands_new[i] = rand;
        }

        GetRandomNumbersLens(tid, &gen2, rands_new);

        float4 yColor = PathTraceF(tid, &xScr, &yScr, rands_new);
        float y = contribFunc(yColor);

        // (2) Markov Chain

        size_t accept = 0;
        size_t largeSteps = 0;
        double accumBrightness = 0.0;

        for (size_t i = 0; i < samplesPerPass; i++) {
            const bool isLargeStep = (rndFloat1_Pseudo(&gen1) < plarge);

            if (isLargeStep) {
                for (size_t i = 0; i < rands_cur.size(); ++i) {
                    rands_new[i] = rndFloat1_Pseudo(&gen2);
                }
            } else {
                const float4 r1 = rndFloat4_Pseudo(&gen2);
                rands_new[0] = MutateKelemen(rands_cur[0], float2(r1.x, r1.y), MUTATE_COEFF_SCREEN * 1.0f, 1024.0f); // screen
                rands_new[1] = MutateKelemen(rands_cur[1], float2(r1.z, r1.w), MUTATE_COEFF_SCREEN * 1.0f, 1024.0f); // screen

                for (size_t i = 2; i < rands_cur.size(); i += 2) {
                    const float4 r1 = rndFloat4_Pseudo(&gen2);
                    rands_new[i + 0] = MutateKelemen(rands_cur[i + 0], float2(r1.x, r1.y), MUTATE_COEFF_BSDF, 1024.0f);
                    rands_new[i + 1] = MutateKelemen(rands_cur[i + 1], float2(r1.z, r1.w), MUTATE_COEFF_BSDF, 1024.0f);
                }
            }

            float yOld = y;
            float4 yOldColor = yColor;

            int xScrOld = xScr, yScrOld = yScr;
            int xScrNew = 0, yScrNew = 0;

            float4 yNewColor = PathTraceF(tid, &xScrNew, &yScrNew, rands_new);
            float yNew = contribFunc(yNewColor);

            float a = (yOld == 0.0f) ? 1.0f : std::min(1.0f, yNew / yOld);
            float p = rndFloat1_Pseudo(&gen1);

            if (p <= a) // accept
            {
                for (size_t i = 0; i < rands_cur.size(); i++) {
                    rands_cur[i] = rands_new[i];
                }

                y = yNew;
                yColor = yNewColor;
                xScr = xScrNew;
                yScr = yScrNew;
                accept++;
            }

            if (isLargeStep) {
                accumBrightness += double(yNew);
                largeSteps++;
            }

            // (5) contrib to image
            //
            float w1 = 1.0f;

            float3 contribAtY = w1 * to_float3(yNewColor) * (1.0f / std::max(yNew, 1e-6f)) * a;
            float3 contribAtX = w1 * to_float3(yOldColor) * (1.0f / std::max(yOld, 1e-6f)) * (1.0f - a);

            if (dot(contribAtX, contribAtX) > 1e-12f) {
                const int offset = yScrOld * m_winWidth + xScrOld;
#pragma omp atomic
                out_color[offset * 4 + 0] += contribAtX.x;
#pragma omp atomic
                out_color[offset * 4 + 1] += contribAtX.y;
#pragma omp atomic
                out_color[offset * 4 + 2] += contribAtX.z;
            }

            if (dot(contribAtY, contribAtY) > 1e-12f) {
                const int offset = yScrNew * m_winWidth + xScrNew;
#pragma omp atomic
                out_color[offset * 4 + 0] += contribAtY.x;
#pragma omp atomic
                out_color[offset * 4 + 1] += contribAtY.y;
#pragma omp atomic
                out_color[offset * 4 + 2] += contribAtY.z;
            }

            progress.Update();
        }

#pragma omp atomic
        avgAcceptanceRate += float(accept);

        double avgBrightness = accumBrightness / double(largeSteps);
#pragma omp atomic
        avgBrightnessOut += avgBrightness;
    }

    progress.Done();

    avgBrightnessOut = avgBrightnessOut / double(m_maxThreadId);
    avgAcceptanceRate = avgAcceptanceRate / float(pixelsNum * a_passNum);

    std::cout << "[IntegratorLangevin]: average brightness      = " << std::fixed << std::setprecision(2) << avgBrightnessOut << std::endl;
    std::cout << "[IntegratorLangevin]: average acceptance rate = " << std::fixed << std::setprecision(2) << 100.0f * avgAcceptanceRate << "%" << std::endl;
    std::cout.flush();

    double actualBrightness = 0.0;
    {
        for (uint i = 0; i < pixelsNum; i++) {
            float4 color = float4(out_color[i * 4 + 0], out_color[i * 4 + 1], out_color[i * 4 + 2], out_color[i * 4 + 3]);
            actualBrightness += contribFunc(color);
        }
        actualBrightness /= double(pixelsNum);
    }

    const float normConst = float(a_passNum) * float(avgBrightnessOut / actualBrightness);
    for (uint i = 0; i < pixelsNum * channels; i++)
        out_color[i] = out_color[i] * normConst;

    shadowPtTime = float(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()) / 1000.f;
}

Integrator::EyeRayData IntegratorLangevin::SampleCameraRay(RandomGen *pGen, uint tid, float *rands) {
    const float4 pixelOffsets = GetRandomNumbersLens(tid, pGen, rands);
    float3 rayDir = EyeRayDirNormalized(pixelOffsets.x, pixelOffsets.y, m_projInv);
    float3 rayPos = float3(0, 0, 0);

    if (m_camLensRadius > 0.0f) {
        const float tFocus = m_camTargetDist / (-rayDir.z);
        const float3 focusPosition = rayPos + rayDir * tFocus;
        const float2 xy = m_camLensRadius * 2.0f * MapSamplesToDisc(float2(pixelOffsets.z - 0.5f, pixelOffsets.w - 0.5f));
        rayPos.x += xy.x;
        rayPos.y += xy.y;
        rayDir = normalize(focusPosition - rayPos);
    } else if (m_enableOpticSim != 0) // not nessesary part of QMC. Just implemented here for test cases, could be moved in main class further
    {
        const float2 xy = 0.25f * m_physSize * float2(2.0f * pixelOffsets.x - 1.0f, 2.0f * pixelOffsets.y - 1.0f);

        rayPos = float3(xy.x, xy.y, 0);

        const float2 rareSam = LensRearRadius() * 2.0f * MapSamplesToDisc(float2(pixelOffsets.z - 0.5f, pixelOffsets.w - 0.5f));
        const float3 shootTo = float3(rareSam.x, rareSam.y, LensRearZ());
        const float3 ray_dirF = normalize(shootTo - rayPos);

        float cosTheta = std::abs(ray_dirF.z);
        rayDir = ray_dirF;
        bool success = TraceLensesFromFilm(rayPos, rayDir);

        if (!success) {
            rayPos = float3(0, -10000000.0, 0.0); // shoot ray under the floor
            rayDir = float3(0, -1, 0);
        } else {
            rayDir = float3(-1, -1, -1) * normalize(rayDir);
            rayPos = float3(-1, -1, -1) * rayPos;
        }
    }

    /////////////////////////////////////////////////////
    uint x = uint(pixelOffsets.x * float(m_winWidth));
    uint y = uint(pixelOffsets.y * float(m_winHeight));
    if (x >= uint(m_winWidth - 1))
        x = uint(m_winWidth - 1);
    if (y >= uint(m_winHeight - 1))
        y = uint(m_winHeight - 1);
    /////////////////////////////////////////////////////

    EyeRayData res;
    {
        res.rayPos = rayPos;
        res.rayDir = rayDir;
        res.x = x;
        res.y = y;
        res.timeSam = 0.0f;
        res.waveSam = 1.0f;
        if (m_normMatrices2.size() != 0)
            res.timeSam = GetRandomNumbersTime(tid, pGen, rands);
        if (KSPEC_SPECTRAL_RENDERING != 0 && m_spectral_mode != 0)
            res.waveSam = GetRandomNumbersSpec(tid, pGen, rands);
        res.cosTheta = 1.0f;
    }

    RecordPixelRndIfNeeded(pixelOffsets, float2(res.waveSam, res.timeSam));

    return res;
}

void IntegratorLangevin::kernel_RayTrace(
    uint tid, uint bounce, const float4 *rayPosAndNear, const float4 *rayDirAndFar,
    const float *a_time, float4 *out_hit1, float4 *out_hit2, float4 *out_hit3,
    uint *out_instId, uint *rayFlags, float *rands
) {
    if (tid >= m_maxThreadId)
        return;
    uint currRayFlags = *rayFlags;
    if (isDeadRay(currRayFlags))
        return;

    const float4 rayPos = *rayPosAndNear;
    const float4 rayDir = *rayDirAndFar;
    const float time = *a_time;

    const CRT_Hit hit = m_pAccelStruct->RayQuery_NearestHit(rayPos, rayDir, time);
    RecordRayHitIfNeeded(bounce, hit);

    if (hit.geomId != uint32_t(-1)) {
        const float2 uv = float2(hit.coords[0], hit.coords[1]);

        // slightly undershoot the intersection to prevent self-intersection and other bugs
        const float3 hitPos = to_float3(rayPos) + hit.t * (1.f - 1e-6f) * to_float3(rayDir);
        // const float3 overHit  = to_float3(rayPos) + hit.t * (1.f + 1e-6f) * to_float3(rayDir);

        // alternative, you may consider Johannes Hanika solution from  Ray Tracing Gems2
        /////////////////////////////////////////////////////////////////////////////////
        // // get distance vectors from triangle vertices
        // vec3 tmpu = P - A, tmpv = P - B, tmpw = P - C
        // // project these onto the tangent planes
        // // defined by the shading normals
        // float dotu = min (0.0, dot(tmpu , nA))
        // float dotv = min (0.0, dot(tmpv , nB))
        // float dotw = min (0.0, dot(tmpw , nC))
        // tmpu -= dotu*nA
        // tmpv -= dotv*nB
        // tmpw -= dotw*nC
        // // finally P' is the barycentric mean of these three
        // vec3 Pp = P + u*tmpu + v*tmpv + w*tmpw
        /////////////////////////////////////////////////////////////////////////////////

        const uint triOffset = m_matIdOffsets[hit.geomId];
        const uint vertOffset = m_vertOffset[hit.geomId];

        const uint A = m_triIndices[(triOffset + hit.primId) * 3 + 0];
        const uint B = m_triIndices[(triOffset + hit.primId) * 3 + 1];
        const uint C = m_triIndices[(triOffset + hit.primId) * 3 + 2];

        const float4 data1 = (1.0f - uv.x - uv.y) * m_vNorm4f[A + vertOffset] + uv.y * m_vNorm4f[B + vertOffset] + uv.x * m_vNorm4f[C + vertOffset];
        const float4 data2 = (1.0f - uv.x - uv.y) * m_vTang4f[A + vertOffset] + uv.y * m_vTang4f[B + vertOffset] + uv.x * m_vTang4f[C + vertOffset];

        float3 hitNorm = to_float3(data1);
        float3 hitTang = to_float3(data2);
        float2 hitTexCoord = float2(data1.w, data2.w);

        // transform surface point with matrix and flip normal if needed
        //
        hitNorm = mul3x3(m_normMatrices[hit.instId], hitNorm);
        hitTang = mul3x3(m_normMatrices[hit.instId], hitTang);

        if (m_normMatrices2.size() > 0) {
            float3 hitNorm2 = mul3x3(m_normMatrices2[hit.instId], hitNorm);
            float3 hitTang2 = mul3x3(m_normMatrices2[hit.instId], hitTang);

            hitNorm = lerp(hitNorm, hitNorm2, time);
            hitTang = lerp(hitTang, hitTang2, time);
        }

        hitNorm = normalize(hitNorm);
        hitTang = normalize(hitTang);

        const float flipNorm = dot(to_float3(rayDir), hitNorm) > 0.001f ? -1.0f : 1.0f; // beware of transparent materials which use normal sign to identity "inside/outside" glass for example
        hitNorm = flipNorm * hitNorm;
        hitTang = flipNorm * hitTang; // do we need this ??

        if (flipNorm < 0.0f)
            currRayFlags |= RAY_FLAG_HAS_INV_NORMAL;
        else
            currRayFlags &= ~RAY_FLAG_HAS_INV_NORMAL;

        const uint midOriginal = m_matIdByPrimId[m_matIdOffsets[hit.geomId] + hit.primId];
        const uint midRemaped = RemapMaterialId(midOriginal, hit.instId);

        *rayFlags = packMatId(currRayFlags, midRemaped);
        *out_hit1 = to_float4(hitPos, hitTexCoord.x);
        *out_hit2 = to_float4(hitNorm, hitTexCoord.y);
        *out_hit3 = to_float4(hitTang, hit.t);
        *out_instId = hit.instId;
    } else {
        const uint flagsToAdd = (bounce == 0) ? (RAY_FLAG_PRIME_RAY_MISS | RAY_FLAG_IS_DEAD | RAY_FLAG_OUT_OF_SCENE) : (RAY_FLAG_IS_DEAD | RAY_FLAG_OUT_OF_SCENE);
        *rayFlags = currRayFlags | flagsToAdd;
    }
}

void IntegratorLangevin::kernel_SampleLightSource(
    uint tid, const float4 *rayPosAndNear, const float4 *rayDirAndFar, const float4 *wavelengths,
    const float4 *in_hitPart1, const float4 *in_hitPart2, const float4 *in_hitPart3,
    const uint *rayFlags, const float *a_time, uint bounce, RandomGen *a_gen, float4 *out_shadeColor,
    float *rands
) {
    if (tid >= m_maxThreadId)
        return;
    const uint currRayFlags = *rayFlags;
    if (isDeadRay(currRayFlags))
        return;

    const uint32_t matId = extractMatId(currRayFlags);
    const float3 ray_dir = to_float3(*rayDirAndFar);

    const float4 data1 = *in_hitPart1;
    const float4 data2 = *in_hitPart2;
    const float4 lambda = *wavelengths;

    SurfaceHit hit;
    hit.pos = to_float3(data1);
    hit.norm = to_float3(data2);
    hit.tang = to_float3(*in_hitPart3);
    hit.uv = float2(data1.w, data2.w);

    const int bounceTmp = int(bounce);
    const float4 acq_rands = GetRandomNumbersLgts(tid, a_gen, bounceTmp, rands);
    const int lightId = std::min(int(std::floor(acq_rands.w * float(m_lights.size()))), int(m_lights.size() - 1u));
    RecordLightRndIfNeeded(bounce, acq_rands);

    if (lightId < 0) // no lights or invalid light id
    {
        *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
        return;
    }

    const LightSample lSam = LightSampleRev(lightId, to_float3(acq_rands), hit.pos);
    const float hitDist = std::sqrt(dot(hit.pos - lSam.pos, hit.pos - lSam.pos));

    const float3 shadowRayDir = normalize(lSam.pos - hit.pos);                                 // explicitSam.direction;
    const float3 shadowRayPos = hit.pos + hit.norm * std::max(maxcomp(hit.pos), 1.0f) * 5e-6f; // TODO: see Ray Tracing Gems, also use flatNormal for offset

    float time = *a_time;
    const bool inIllumArea = (dot(shadowRayDir, lSam.norm) < 0.0f) || lSam.isOmni || lSam.hasIES;
    const bool needShade = inIllumArea && !m_pAccelStruct->RayQuery_AnyHit(to_float4(shadowRayPos, 0.0f), to_float4(shadowRayDir, hitDist * 0.9995f), time); /// (!!!) expression-way, RT pipeline bug work around, if change check test_213
    RecordShadowHitIfNeeded(bounce, needShade);

    if (needShade) /// (!!!) expression-way to compute 'needShade', RT pipeline bug work around, if change check test_213
    {
        const BsdfEval bsdfV = MaterialEval(matId, lambda, shadowRayDir, (-1.0f) * ray_dir, hit.norm, hit.tang, hit.uv);
        float cosThetaOut = std::max(dot(shadowRayDir, hit.norm), 0.0f);

        float lgtPdfW = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, shadowRayPos, shadowRayDir, lSam.pos, lSam.norm, lSam.pdf);
        float misWeight = (m_intergatorType == INTEGRATOR_MIS_PT) ? misWeightHeuristic(lgtPdfW, bsdfV.pdf) : 1.0f;
        const bool isDirect = (m_lights[lightId].geomType == LIGHT_GEOM_DIRECT);
        const bool isPoint = (m_lights[lightId].geomType == LIGHT_GEOM_POINT);

        if (isDirect) {
            misWeight = 1.0f;
            lgtPdfW = 1.0f;
        } else if (isPoint)
            misWeight = 1.0f;

        const bool isDirectLight = !hasNonSpecular(currRayFlags);
        if ((m_renderLayer == FB_DIRECT && !isDirectLight) ||
            (m_renderLayer == FB_INDIRECT && isDirectLight)) // skip some number of bounces if this is set
            misWeight = 0.0f;

        const float4 lightColor = LightIntensity(lightId, lambda, shadowRayPos, shadowRayDir);
        *out_shadeColor = (lightColor * bsdfV.val / lgtPdfW) * cosThetaOut * misWeight;
    } else
        *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
}

void IntegratorLangevin::kernel_NextBounce(
    uint tid, uint bounce, const float4 *in_hitPart1, const float4 *in_hitPart2, const float4 *in_hitPart3,
    const uint *in_instId, const float4 *in_shadeColor, float4 *rayPosAndNear, float4 *rayDirAndFar,
    const float4 *wavelengths, float4 *accumColor, float4 *accumThoroughput,
    RandomGen *a_gen, MisData *misPrev, uint *rayFlags, float *rands) {
    if (tid >= m_maxThreadId)
        return;
    const uint currRayFlags = *rayFlags;
    if (isDeadRay(currRayFlags))
        return;

    const uint32_t matId = extractMatId(currRayFlags);

    // process surface hit case
    //
    const float3 ray_dir = to_float3(*rayDirAndFar);
    const float3 ray_pos = to_float3(*rayPosAndNear);
    const float4 lambda = *wavelengths;

    const float4 data1 = *in_hitPart1;
    const float4 data2 = *in_hitPart2;

    SurfaceHit hit;
    hit.pos = to_float3(data1);
    hit.norm = to_float3(data2);
    hit.tang = to_float3(*in_hitPart3);
    hit.uv = float2(data1.w, data2.w);

    const float hitDist = in_hitPart3->w;

    const MisData prevBounce = *misPrev;
    const float prevPdfW = prevBounce.matSamplePdf;

    // process light hit case
    //
    if (m_materials[matId].mtype == MAT_TYPE_LIGHT_SOURCE) {
        const uint texId = m_materials[matId].texid[0];
        const float2 texCoordT = mulRows2x4(m_materials[matId].row0[0], m_materials[matId].row1[0], hit.uv);
        const float4 texColor = m_textures[texId]->sample(texCoordT);
        const uint lightId = m_instIdToLightInstId[*in_instId];

        const float4 emissColor = m_materials[matId].colors[EMISSION_COLOR];
        float4 lightIntensity = emissColor * texColor;

        if (lightId != 0xFFFFFFFF) {
            const float lightCos = dot(to_float3(*rayDirAndFar), to_float3(m_lights[lightId].norm));
            const float lightDirectionAtten = (lightCos < 0.0f || m_lights[lightId].geomType == LIGHT_GEOM_SPHERE) ? 1.0f : 0.0f;
            lightIntensity = LightIntensity(lightId, lambda, ray_pos, to_float3(*rayDirAndFar)) * lightDirectionAtten;
        }

        float misWeight = 1.0f;
        if (m_intergatorType == INTEGRATOR_MIS_PT) {
            if (bounce > 0 && lightId != 0xFFFFFFFF) {
                const float lgtPdf = LightPdfSelectRev(lightId) * LightEvalPDF(lightId, ray_pos, ray_dir, hit.pos, hit.norm, 1.0f);
                misWeight = misWeightHeuristic(prevPdfW, lgtPdf);
                if (prevPdfW <= 0.0f) // specular bounce
                    misWeight = 1.0f;
            }
        } else if (m_intergatorType == INTEGRATOR_SHADOW_PT && hasNonSpecular(currRayFlags))
            misWeight = 0.0f;

        const bool isDirectLight = !hasNonSpecular(currRayFlags);
        const bool isFirstNonSpec = (currRayFlags & RAY_FLAG_FIRST_NON_SPEC) != 0;
        if (m_renderLayer == FB_INDIRECT && (isDirectLight || isFirstNonSpec))
            misWeight = 0.0f;

        float4 currAccumColor = *accumColor;
        float4 currAccumThroughput = *accumThoroughput;

        currAccumColor += currAccumThroughput * lightIntensity * misWeight;

        *accumColor = currAccumColor;
        *rayFlags = currRayFlags | (RAY_FLAG_IS_DEAD | RAY_FLAG_HIT_LIGHT);
        return;
    }

    const uint bounceTmp = bounce;
    const BsdfSample matSam = MaterialSampleAndEval(matId, tid, bounceTmp, lambda, a_gen, (-1.0f) * ray_dir, hit.norm, hit.tang, hit.uv, misPrev, currRayFlags, rands);
    const float4 bxdfVal = matSam.val * (1.0f / std::max(matSam.pdf, 1e-20f));
    const float cosTheta = std::abs(dot(matSam.dir, hit.norm));

    MisData nextBounceData = *misPrev; // remember current pdfW for next bounce
    nextBounceData.matSamplePdf = (matSam.flags & RAY_EVENT_S) != 0 ? -1.0f : matSam.pdf;
    nextBounceData.cosTheta = cosTheta;
    *misPrev = nextBounceData;

    if (m_intergatorType == INTEGRATOR_STUPID_PT) {
        *accumThoroughput *= cosTheta * bxdfVal;
    } else if (m_intergatorType == INTEGRATOR_SHADOW_PT || m_intergatorType == INTEGRATOR_MIS_PT) {
        const float4 currThoroughput = *accumThoroughput;
        const float4 shadeColor = *in_shadeColor;
        float4 currAccumColor = *accumColor;

        currAccumColor += currThoroughput * shadeColor;
        *accumColor = currAccumColor;
        *accumThoroughput = currThoroughput * cosTheta * bxdfVal;
    }

    // compute point on the other side of the surface in case of transmission
    if ((matSam.flags & RAY_EVENT_T) != 0) {
        hit.pos = hit.pos + hitDist * ray_dir * 2 * 1e-6f;
    }

    *rayPosAndNear = to_float4(OffsRayPos(hit.pos, hit.norm, matSam.dir), 0.0f); // todo: use flatNormal for offset
    *rayDirAndFar = to_float4(matSam.dir, FLT_MAX);

    uint nextFlags = ((currRayFlags & ~RAY_FLAG_FIRST_NON_SPEC) | matSam.flags); // always force reset RAY_FLAG_FIRST_NON_SPEC;
    if (m_renderLayer == FB_DIRECT && hasNonSpecular(currRayFlags))              // NOTE: use currRayFlags for check, not nextFlags because of MIS: a ray may hit light source in next bounce
        nextFlags |= RAY_FLAG_IS_DEAD;                                           //       but if we already have non specular bounce previously, definitely can stop
    else if (!hasNonSpecular(currRayFlags) && hasNonSpecular(nextFlags))
        nextFlags |= RAY_FLAG_FIRST_NON_SPEC;
    *rayFlags = nextFlags;
}

void IntegratorLangevin::kernel_HitEnvironment(
    uint tid, const uint *rayFlags, const float4 *rayDirAndFar,
    const MisData *a_prevMisData, const float4 *accumThoroughput,
    float4 *accumColor, float *rands) {
    if (tid >= m_maxThreadId)
        return;
    const uint currRayFlags = *rayFlags;
    if (!isOutOfScene(currRayFlags))
        return;

    float envPdf = 1.0f;
    float4 envColor = EnvironmentColor(to_float3(*rayDirAndFar), envPdf);

    const auto misPrev = *a_prevMisData;
    const bool isSpec = isSpecular(&misPrev);
    const bool exitZero = (currRayFlags & RAY_FLAG_PRIME_RAY_MISS) != 0;

    if (m_intergatorType == INTEGRATOR_MIS_PT && m_envEnableSam != 0 && !isSpec && !exitZero) {
        float lgtPdf = LightPdfSelectRev(m_envLightId) * envPdf;
        float bsdfPdf = misPrev.matSamplePdf;
        float misWeight = misWeightHeuristic(bsdfPdf, lgtPdf); // (bsdfPdf*bsdfPdf) / (lgtPdf*lgtPdf + bsdfPdf*bsdfPdf);
        envColor *= misWeight;
    } else if (m_intergatorType == INTEGRATOR_SHADOW_PT && m_envEnableSam != 0) {
        envColor = float4(0.0f);
    }

    const uint camBackId = m_envCamBackId;
    if (exitZero && camBackId != uint(-1)) // apply camera back color to ray
    {
        const uint XY = m_packedXY[tid];
        const uint x = (XY & 0x0000FFFF);
        const uint y = (XY & 0xFFFF0000) >> 16;

        const float2 texCoord = float2((float(x) + 0.5f) / float(m_winWidth),
                                       (float(y) + 0.5f) / float(m_winHeight));

        envColor = m_textures[camBackId]->sample(texCoord);
    }

    if (m_intergatorType == INTEGRATOR_STUPID_PT) // todo: when explicit sampling will be added, disable contribution here for 'INTEGRATOR_SHADOW_PT'
        *accumColor = (*accumThoroughput) * envColor;
    else
        *accumColor += (*accumThoroughput) * envColor;
}

BsdfSample IntegratorLangevin::MaterialSampleAndEval(
    uint a_materialId, uint tid, uint bounce, float4 wavelengths, 
    RandomGen *a_gen, float3 v, float3 n, float3 tan, float2 tc,
    MisData *a_misPrev, const uint a_currRayFlags, float *rands
) {
    BsdfSample res;
    {
        res.val = float4(0, 0, 0, 0);
        res.pdf = 1.0f;
        res.dir = float3(0, 1, 0);
        res.ior = 1.0f;
        res.flags = a_currRayFlags;
        res.ior = 1.0f;
    }

    uint32_t currMatId = a_materialId;
    uint mtype = m_materials[currMatId].mtype;
    uint layer = 0;
    while (KSPEC_MAT_TYPE_BLEND != 0 && mtype == MAT_TYPE_BLEND) {
        currMatId = BlendSampleAndEval(currMatId, tid, bounce, layer, wavelengths, a_gen, v, n, tc, a_misPrev, &res, rands);
        mtype = m_materials[currMatId].mtype;
        layer++;
    }

    // BSDF is multiplied (outside) by cosThetaOut1.
    // When normal map is enables this becames wrong because normal is changed;
    // First : return cosThetaOut in sam;
    // Second: apply cos(theta2)/cos(theta1) to cos(theta1) to get cos(theta2)
    //
    const uint normalMapId = m_materials[currMatId].texid[1];
    const float3 geomNormal = n;
    float3 shadeNormal = n;

    if (KSPEC_BUMP_MAPPING != 0 && normalMapId != 0xFFFFFFFF)
        shadeNormal = BumpMapping(normalMapId, currMatId, geomNormal, tan, tc);

    const float2 texCoordT = mulRows2x4(m_materials[currMatId].row0[0], m_materials[currMatId].row1[0], tc);
    const uint texId = m_materials[currMatId].texid[0];
    const float4 texColor = m_textures[texId]->sample(texCoordT);
    const float4 acq_rands = GetRandomNumbersMats(tid, a_gen, int(bounce), rands);
    const uint cflags = m_materials[currMatId].cflags;
    RecordMatRndNeeded(bounce, acq_rands);

    float4 fourScalarMatParams = float4(1, 1, 1, 1);
    if (KSPEC_MAT_FOUR_TEXTURES != 0 && (cflags & FLAG_FOUR_TEXTURES) != 0) {
        const uint texId2 = m_materials[currMatId].texid[2];
        const uint texId3 = m_materials[currMatId].texid[3];

        const float2 texCoord2T = mulRows2x4(m_materials[currMatId].row0[2], m_materials[currMatId].row1[2], tc);
        const float2 texCoord3T = mulRows2x4(m_materials[currMatId].row0[3], m_materials[currMatId].row1[3], tc);

        const float4 color2 = m_textures[texId2]->sample(texCoord2T);
        const float4 color3 = m_textures[texId3]->sample(texCoord3T);

        if ((cflags & FLAG_PACK_FOUR_PARAMS_IN_TEXTURE) != 0)
            fourScalarMatParams = color2;
        else
            fourScalarMatParams = float4(color2.x, color3.x, 1, 1);
    }

    switch (mtype) {
    case MAT_TYPE_GLTF:
        if (KSPEC_MAT_TYPE_GLTF != 0) {
            const float4 color = m_materials[currMatId].colors[GLTF_COLOR_BASE] * texColor;
            gltfSampleAndEval(m_materials.data() + currMatId, acq_rands, v, shadeNormal, tc, color, fourScalarMatParams, &res);
        }
        break;
    case MAT_TYPE_GLASS:
        if (KSPEC_MAT_TYPE_GLASS != 0) {
            glassSampleAndEval(m_materials.data() + currMatId, acq_rands, v, geomNormal, tc, &res, a_misPrev);
        }
        break;
    case MAT_TYPE_CONDUCTOR:
        if (KSPEC_MAT_TYPE_CONDUCTOR != 0) {
            const float3 alphaTex = to_float3(texColor);
            const float2 alpha = float2(m_materials[currMatId].data[CONDUCTOR_ROUGH_V], m_materials[currMatId].data[CONDUCTOR_ROUGH_U]);
            const float4 etaSpec = SampleMatParamSpectrum(currMatId, wavelengths, CONDUCTOR_ETA, 0);
            const float4 kSpec = SampleMatParamSpectrum(currMatId, wavelengths, CONDUCTOR_K, 1);
            if (trEffectivelySmooth(alpha))
                conductorSmoothSampleAndEval(m_materials.data() + currMatId, etaSpec, kSpec, acq_rands, v, shadeNormal, tc, &res);
            else
                conductorRoughSampleAndEval(m_materials.data() + currMatId, etaSpec, kSpec, acq_rands, v, shadeNormal, tc, alphaTex, &res);
        }
        break;
    case MAT_TYPE_DIFFUSE:
        if (KSPEC_MAT_TYPE_DIFFUSE != 0) {
            const float4 color = texColor;
            float4 reflSpec = SampleMatColorSpectrumTexture(currMatId, wavelengths, DIFFUSE_COLOR, 0, tc);
            if (m_spectral_mode == 0)
                reflSpec *= color;

            diffuseSampleAndEval(m_materials.data() + currMatId, reflSpec, acq_rands, v, shadeNormal, tc, &res);
        }
        break;
    case MAT_TYPE_PLASTIC:
        if (KSPEC_MAT_TYPE_PLASTIC != 0) {
            const float4 color = texColor;
            float4 reflSpec = SampleMatColorSpectrumTexture(currMatId, wavelengths, PLASTIC_COLOR, 0, tc);
            // float4 reflSpec    = SampleMatColorParamSpectrum(currMatId, wavelengths, PLASTIC_COLOR, 0);
            if (m_spectral_mode == 0)
                reflSpec *= color;

            const uint precomp_id = m_materials[currMatId].datai[0];

            plasticSampleAndEval(m_materials.data() + currMatId, reflSpec, acq_rands, v, shadeNormal, tc, &res,
                                 m_precomp_coat_transmittance.data() + precomp_id * MI_ROUGH_TRANSMITTANCE_RES);
        }
        break;
    case MAT_TYPE_DIELECTRIC:
        if (KSPEC_MAT_TYPE_DIELECTRIC != 0) {
            const float4 intIORSpec = SampleMatParamSpectrum(currMatId, wavelengths, DIELECTRIC_ETA_INT, 0);
            const uint specId = m_materials[currMatId].spdid[0];
            dielectricSmoothSampleAndEval(m_materials.data() + currMatId, intIORSpec, a_misPrev->ior, acq_rands, v, shadeNormal, tc, &res);

            res.flags |= (specId < 0xFFFFFFFF) ? RAY_FLAG_WAVES_DIVERGED : 0;

            a_misPrev->ior = res.ior;
        }
        break;
    default:
        break;
    }

    // BSDF is multiplied (outside) by cosThetaOut1.
    // When normal map is enables this becames wrong because normal is changed;
    // First : return cosThetaOut in sam;
    // Second: apply cos(theta2)/cos(theta1) to cos(theta1) to get cos(theta2)
    //
    if (KSPEC_BUMP_MAPPING != 0 && normalMapId != 0xFFFFFFFF) {
        const float cosThetaOut1 = std::abs(dot(res.dir, geomNormal));
        const float cosThetaOut2 = std::abs(dot(res.dir, shadeNormal));
        res.val *= cosThetaOut2 / std::max(cosThetaOut1, 1e-10f);
    }

    return res;
}

uint32_t IntegratorLangevin::BlendSampleAndEval(
    uint a_materialId, uint tid, uint bounce, uint layer, float4 wavelengths, RandomGen* a_gen, float3 v, float3 n, float2 tc, 
    MisData* a_misPrev, BsdfSample* a_pRes, float *rands
) {
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const uint   texId     = m_materials[a_materialId].texid[0];
  const float4 weightDat = m_textures[texId]->sample(texCoordT);
  const float  weightTex = weightDat.x;
  const float  weight    = m_materials[a_materialId].data[BLEND_WEIGHT] * weightTex;

  const uint matId1 = m_materials[a_materialId].datai[0];
  const uint matId2 = m_materials[a_materialId].datai[1];

  uint32_t selectedMatId = matId1;
  const float select = GetRandomNumbersMatB(tid, a_gen, int(bounce), int(layer), rands);
  RecordBlendRndNeeded(bounce, layer, select);

  if(select < weight)
  {
    a_pRes->pdf *= weight;
    a_pRes->val *= weight;
    selectedMatId = matId2;
  }
  else
  {
    a_pRes->pdf *= 1.0f - weight;
    a_pRes->val *= 1.0f - weight;
    selectedMatId = matId1;
  }

  return selectedMatId;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Integrator> CreateIntegratorLangevin(int a_maxThreads = 1, int a_spectral_mode = 0, std::vector<uint32_t> a_features = {}) {
    return std::make_shared<IntegratorLangevin>(a_maxThreads, a_spectral_mode, a_features);
}
