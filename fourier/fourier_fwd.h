#ifndef FOURIER_FOURIER__FW_H_
#define FOURIER_FOURIER__FW_H_
#include "include/cglobals.h"

namespace fourier {
    static constexpr int _SIZE = 10;
    static constexpr int _SAMPLING_STEP = int(_SIZE / M_PI); 


    inline float to_phase(float wl, float start, float end)
    {
        return LiteMath::M_PI * (wl - start) / (end - start) - LiteMath::M_PI;
    }

    inline std::vector<float> wl_to_phases(std::vector<float> wavelenghts, float start = LAMBDA_MIN, float end = LAMBDA_MAX)
    {
        for(unsigned i = 0; i < wavelenghts.size(); ++i) {
            wavelenghts[i] = to_phase(wavelenghts[i], start, end);
        }
        return wavelenghts;
    }

    inline std::vector<float> make_sampling_phases()
    {
        const int count = int((LAMBDA_MAX - LAMBDA_MIN) / fourier::_SAMPLING_STEP);

        std::vector<float> data(count + LAMBDA_MIN + count * fourier::_SAMPLING_STEP != LAMBDA_MAX);
        for(int i = 0; i < count; ++i) {
            data[i] = fourier::to_phase(LAMBDA_MIN + float(i) * fourier::_SAMPLING_STEP, LAMBDA_MIN, LAMBDA_MAX);
        }
        if(LAMBDA_MIN + count * fourier::_SAMPLING_STEP != LAMBDA_MAX) {
            data.back() = fourier::to_phase(LAMBDA_MAX, LAMBDA_MIN, LAMBDA_MAX);
        }
        return data;
    }

}
#endif