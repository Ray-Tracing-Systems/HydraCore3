#include "fourier.h"
#include "levinson.h"
#include <cassert>
#include "LiteMath.h"

#include <iostream>

namespace fourier
{

    float to_phase(float wl, float start, float end)
    {
        return std::fma(LiteMath::M_PI, (wl - start) / (end - start), -LiteMath::M_PI);
    }

 
    fvec real_fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values, int n)
    {
        const unsigned N = unsigned(phases.size());
        assert(N == values.size());
        fvec moments;
        for(uint i = 0; i <= FOURIER_M; ++i) {
            float val = 0.0f;
            for(unsigned j = 0; j < N; ++j) {
                val += values[j] * std::cos(float(i) * phases[j]);
            }
            moments[i] = val / float(N) * 2.f;
        }
        moments[0] = moments[0] * 0.5f;
        return moments;
    }

    fvec_c fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values, int n)
    {
        const unsigned N = unsigned(phases.size());
        assert(N == values.size());
        fvec_c moments;
        for(uint i = 0; i <= FOURIER_M; ++i) {
            Complex val{0.0f, 0.0f};
            for(unsigned j = 0; j < N; ++j) {
                val += values[j] * std::exp(-I * float(i) * phases[j]);
            }
            moments[i] = val / float(N) * 2.f;
        }
        moments[0] = moments[0] * 0.5f;
        return moments;
    }


    fvec_c precompute_mese_coeffs(const fvec_c &gamma)
    {
        std::vector<float> e0(FOURIER_M + 1, 0.0f);
        e0[0] = 1.0f;
        std::vector<Complex> data(2 * FOURIER_M + 1); //Matrix values
        data[FOURIER_M] = LiteMath::INV_TWOPI * gamma[0];
        for(uint i = 1; i <= FOURIER_M; ++i) {
            data[FOURIER_M + i] = LiteMath::INV_TWOPI * gamma[i];
            data[FOURIER_M - i] = LiteMath::INV_TWOPI * std::conj(gamma[i]);
        }
        auto res0 = levinson<Complex>(data, e0);
        fvec_c res;
        std::copy(res0.begin(), res0.end(), res.v);
        return res;
    }

    fvec precompute_mese_coeffs(const fvec &gamma)
    {
        std::vector<float> e0(FOURIER_M + 1, 0.0f);
        e0[0] = 1.0f;
        std::vector<float> data(2 * FOURIER_M + 1); //Matrix values
        data[FOURIER_M] = LiteMath::INV_TWOPI * gamma[0];
        for(uint i = 1; i <= FOURIER_M; ++i) {
            data[FOURIER_M + i] = LiteMath::INV_TWOPI * gamma[i];
            data[FOURIER_M - i] = LiteMath::INV_TWOPI * gamma[i]; //[std::conj(gamma[i])] in fact;
        }
        auto res0 = levinson<float>(data, e0);
        fvec res;
        std::copy(res0.begin(), res0.end(), res.v);
        return res;
    }

    fvec_c exponential_moments(const fvec moments, Complex &gamma0)
    {
        gamma0 = 0.5f * LiteMath::INV_TWOPI * std::exp(LiteMath::M_PI * I * (moments[0] - 0.5f));
        fvec_c res;
        res[0] = 2.0f * gamma0.real();
        for(uint i = 1; i <= FOURIER_M; ++i) {
            Complex sm{0.0f, 0.0f};
            for(int j = 1; j <= int(i) - 1; ++j) {
                sm += float(i - j) * res[j] * moments[i - j]; 
            }

            res[i] = LiteMath::M_TWOPI * I * (float(i) * gamma0 * moments[i] + sm) / float(i);
        }
        return res;
    }

    fvec_c lagrange_multipliers(const fvec &moments)
    {      
        Complex gamma0;
        fvec_c gamma = exponential_moments(moments, gamma0);

        //Calculating q
        fvec_c q = precompute_mese_coeffs(gamma);

        fvec_c lambda;
        for(uint i = 0; i <= FOURIER_M; ++i) {
            Complex t{0.0f, 0.0f};

            for(uint k = 0; k <= FOURIER_M - i; ++k) {
                Complex ts{0.0f, 0.0f};
                for(uint j = 0; j <= FOURIER_M - k - i; ++j) {
                    ts += std::conj(q[j + k + i]) * q[j];
                }
                t += ts * (k == 0 ? gamma0 : gamma[k]);
            }
            lambda[i] = t / (LiteMath::M_PI * I * q[0]);
        }
        return lambda;
    }

    float bounded_mese_l(float phase, const fvec_c &lagrange_m)
    {
        float sum = std::real(lagrange_m[0]);
        for(uint i = 1; i <= FOURIER_M; ++i) {
            sum += 2.0f * std::real(lagrange_m[i] * std::exp(-I * float(i) * phase));
        } 

        return LiteMath::INV_PI * std::atan(sum) + 0.5f;
    }


    float mese_precomp(float phase, const fvec &q)
    {
        Complex t = 0.0f;
        for(uint i = 0; i <= FOURIER_M; ++i) t += LiteMath::INV_TWOPI * q[i] * std::exp(-I * float(i) * phase); 

        float div = std::fabs(t);
        div *= div;

        return (LiteMath::INV_TWOPI * std::real(q[0])) / div;
    }

}