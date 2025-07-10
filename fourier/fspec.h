#ifndef FOURIER_FSPEC_H_
#define FOURIER_FSPEC_H_
#include <vector>
#include <cmath>
#include "fourier_fwd.h"

struct FourierSpec;

namespace fourier {
    extern FourierSpec real_fourier_moments_of(const std::vector<float> &phases, const std::vector<float> &values);
    extern std::vector<float> fourier_series(const std::vector<float> &phases, const FourierSpec &spec);
};

struct FourierSpec
{
    static constexpr int SIZE = fourier::_SIZE;
    static constexpr int M = SIZE - 1;
    static constexpr int SAMPLING_STEP = fourier::_SAMPLING_STEP; 
    static constexpr int SAMPLING_SIZE = int(LAMBDA_MAX - LAMBDA_MIN) / SAMPLING_STEP + int(LAMBDA_MAX - LAMBDA_MIN) % SAMPLING_STEP;
    static_assert(SAMPLING_STEP > 0);
    static inline std::vector<float> SAMPLING_PHASES = {}; 

    inline static bool unpack_on_multiply = false;


    float v[SIZE];

    FourierSpec() : v{} 
    {
        for(int i = 0; i < SIZE; ++i) {
            v[i] = 0.0f;
        }
    }

    explicit FourierSpec(const float *ptr, int max_size = SIZE) : v{} 
    { 
        int sz = max_size < SIZE ? max_size : SIZE;
        for(int i = 0; i < sz; ++i) {
            v[i] = ptr[i];
        }
    }

    FourierSpec(const FourierSpec &other) = default;

    explicit FourierSpec(float x) : v{} {
        v[0] = x; 
        for(int i = 1; i < SIZE; ++i) {
            v[i] = 0.0f;
        }
    }

    FourierSpec &operator=(const FourierSpec &other) = default;

    FourierSpec operator+(const FourierSpec &other) const 
    {
        FourierSpec c(*this);
        c += other;
        return c;
    }

    static FourierSpec  __attribute__((optimize("O3"))) convolve_freq(const FourierSpec &a, const FourierSpec &b)
    {
        FourierSpec res;
        
        for(int i = 0; i < SIZE; ++i) {
            for(int j = 0; j < SIZE; ++j) {
                if(i + j < SIZE) {
                    res[i + j] += a[i] * b[j];
                }
                if(std::abs(i - j) < SIZE) {
                    res[std::abs(i - j)] += a[i] * b[j];
                }
            }
        }

/*
        for(int i = 0; i < SIZE; ++i) {
            for(int j = 0; j <= i; ++j) {
                res[i] += a[j] * b[i - j];// + a[i - j] * b[j];
                if(i + j < SIZE)
                    res[i] += a[j] * b[i + j];// + a[i + j] * b[j];
            }
            for(int j = i + 1; j < SIZE; ++j) {
                if(i + j < SIZE)
                    res[i] += a[j] * b[i + j];// + a[i + j] * b[j];
                res[i] += a[j] * b[j - i];// + a[j - i] * b[j];
            } 
        }
*/
        return res * 0.5f;
    }

    FourierSpec operator-(const FourierSpec &other) const 
    {
        FourierSpec c(*this);
        c -= other;
        return c;
    }

    FourierSpec operator*(float f) const 
    {
        FourierSpec c(*this);
        c *= f;
        return c;
    }

    FourierSpec operator/(float f) const 
    {
        FourierSpec c(*this);
        c /= f;
        return c;
    }

    FourierSpec operator*(const FourierSpec &other) const 
    {
        if(unpack_on_multiply) {
            const std::vector<float> values1 = fourier::fourier_series(SAMPLING_PHASES, *this);
            std::vector<float> values2 = fourier::fourier_series(SAMPLING_PHASES, other);
            for(size_t i = 0; i < values1.size(); ++i) {
                values2[i] *= values1[i];
            }
            return fourier::real_fourier_moments_of(SAMPLING_PHASES, values2);
        }
        else return convolve_freq(*this, other);
    }

    //Convolution in time domain
    FourierSpec operator&(const FourierSpec &other) const
    {
        FourierSpec c(*this);
        c &= other;
        return c;
    }

    float &operator[](int i) 
    {
        return v[i];
    }

    float operator[](int i) const 
    {
        return v[i];
    }

    FourierSpec operator-() const 
    {
        FourierSpec res;
        for(uint i = 0; i < SIZE; ++i) {
            res.v[i] = -v[i];
        }
        return res;
    }

    FourierSpec &operator+=(const FourierSpec &other) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] += other.v[i];
        }
        return *this;
    }

    FourierSpec &operator-=(const FourierSpec &other) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] -= other.v[i];
        }
        return *this;
    }

    FourierSpec &operator*=(float f) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] *= f;
        }
        return *this;
    }

    FourierSpec &operator/=(float f) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] /= f;
        }
        return *this;
    }

    FourierSpec &operator*=(const FourierSpec &other) 
    {
        return *this = std::move((*this) * other);
    }

    //Convolution in time domain
    FourierSpec &operator&=(const FourierSpec &other) 
    {
        for(uint i = 0; i < SIZE; ++i) {
            v[i] *= other.v[i];
        }
        return *this;
    }

    bool operator==(const FourierSpec &other) const 
    {
        for(int i = 0; i < SIZE; ++i) {
            if(v[i] != other.v[i]) return false;
        }
        return true;
    }

    bool operator!=(const FourierSpec &other) const 
    {
        for(int i = 0; i < SIZE; ++i) {
            if(v[i] != other.v[i]) return true;
        }
        return false;
    }


};

inline FourierSpec operator*(float f, const FourierSpec &spec)
{
    return spec * f;
}


struct BsdfEvalF
{
  FourierSpec val;
  float  pdf; 
};

struct BsdfSampleF
{
  FourierSpec val;
  float3 dir;
  float  pdf; 
  uint   flags;
  float  ior;
};



#endif