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

    std::vector<float> real_fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values, int n)
    {
        const unsigned N = phases.size();
        assert(N == values.size());
        int M = n - 1;
        std::vector<float> moments(M + 1);
        for(int i = 0; i <= M; ++i) {
            Complex val{0.0f, 0.0f};
            for(unsigned j = 0; j < N; ++j) {
                val += values[j] * std::exp(-I * float(i) * phases[j]);
            }
            moments[i] = std::real(val) / float(N);
        }
        return moments;
    }

    std::vector<Complex> fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values, int n)
    {
        const unsigned N = phases.size();
        assert(N == values.size());
        int M = n - 1;
        std::vector<Complex> moments(M + 1);
        for(int i = 0; i <= M; ++i) {
            Complex val{0.0f, 0.0f};
            for(unsigned j = 0; j < N; ++j) {
                val += values[j] * std::exp(-I * float(i) * phases[j]);
            }
            moments[i] = val / float(N);
        }
        return moments;
    }

    static std::vector<Complex> exponential_moments(const std::vector<float> moments, Complex &gamma0)
    {
        gamma0 = 0.5f * LiteMath::INV_TWOPI * std::exp(LiteMath::M_PI * I * (moments[0] - 0.5f));
        const int M = moments.size() - 1;
        std::vector<Complex> res(M + 1);
        res[0] = 2.0f * gamma0.real();
        for(int i = 1; i <= M; ++i) {
            Complex sm{0.0f, 0.0f};
            for(int j = 1; j <= i - 1; ++j) {
                sm += float(i - j) * res[j] * moments[i - j]; 
            }

            res[i] = LiteMath::M_TWOPI * I * (float(i) * gamma0 * moments[i] + sm) / float(i);
        }
        return res;
    }


    std::vector<Complex> precompute_mese_coeffs(const std::vector<Complex> &gamma)
    {
        const int M = gamma.size() - 1;
        std::vector<float> e0(M + 1, 0.0f);
        e0[0] = 1.0f;
        std::vector<Complex> data(2 * M + 1); //Matrix values
        data[M] = LiteMath::INV_TWOPI * gamma[0];
        for(int i = 1; i <= M; ++i) {
            data[M + i] = LiteMath::INV_TWOPI * gamma[i];
            data[M - i] = LiteMath::INV_TWOPI * std::conj(gamma[i]);
        }
        return levinson<Complex>(data, e0);
    }

    std::vector<float> precompute_mese_coeffs(const std::vector<float> &gamma)
    {
        const int M = gamma.size() - 1;
        std::vector<float> e0(M + 1, 0.0f);
        e0[0] = 1.0f;
        std::vector<float> data(2 * M + 1); //Matrix values
        data[M] = LiteMath::INV_TWOPI * gamma[0];
        for(int i = 1; i <= M; ++i) {
            data[M + i] = LiteMath::INV_TWOPI * gamma[i];
            data[M - i] = LiteMath::INV_TWOPI * gamma[i]; //[std::conj(gamma[i])] in fact;
        }
        return levinson<float>(data, e0);
    }

    std::vector<Complex> lagrange_multipliers(const std::vector<float> &moments)
    {      
        const int M = moments.size() - 1;
        Complex gamma0;
        std::vector<Complex> gamma = exponential_moments(moments, gamma0);

        //Calculating q
        std::vector<Complex> q = precompute_mese_coeffs(gamma);

        std::vector<Complex> lambda(M + 1);
        for(int i = 0; i <= M; ++i) {
            Complex t{0.0f, 0.0f};

            for(int k = 0; k <= M - i; ++k) {
                Complex ts{0.0f, 0.0f};
                for(int j = 0; j <= M - k - i; ++j) {
                    ts += std::conj(q[j + k + i]) * q[j];
                }
                t += ts * (k == 0 ? gamma0 : gamma[k]);
            }
            lambda[i] = t / (LiteMath::M_PI * I * q[0]);
        }
        return lambda;
    }

    float bounded_mese_l(float phase, const std::vector<Complex> &lagrange_m)
    {
        const int M = lagrange_m.size() - 1;
        float sum = std::real(lagrange_m[0]);
        for(int i = 1; i <= M; ++i) {
            sum += 2.0f * std::real(lagrange_m[i] * std::exp(-I * float(i) * phase));
        } 

        return LiteMath::INV_PI * std::atan(sum) + 0.5f;
    }


    float mese_precomp(float phase, const std::vector<float> &q)
    {
        const int M = q.size() - 1;
        Complex t = 0.0f;
        for(int i = 0; i <= M; ++i) t += LiteMath::INV_TWOPI * q[i] * std::exp(-I * float(i) * phase); 

        float div = std::fabs(t);
        div *= div;

        return (LiteMath::INV_TWOPI * std::real(q[0])) / div;
    }

}