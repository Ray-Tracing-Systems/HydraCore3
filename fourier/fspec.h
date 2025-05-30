#ifndef FOURIER_FSPEC_H_
#define FOURIER_FSPEC_H_
#include <algorithm>
#include <numeric>




struct FourierSpec
{
    static constexpr int SIZE = 8;
    static constexpr int M = SIZE - 1;

    float v[SIZE];

    FourierSpec() : v{} { std::fill_n(v, SIZE, 0.0f); }
    explicit FourierSpec(const float *ptr, int max_size = SIZE) : v{} { std::copy(ptr, ptr + std::min(SIZE, max_size), v); }
    FourierSpec(const FourierSpec &other) = default;
    explicit FourierSpec(float x) : v{} { v[0] = x; std::fill_n(v + 1, SIZE - 1, 0.0f); }

    FourierSpec &operator=(const FourierSpec &other) = default;

    FourierSpec operator+(const FourierSpec &other) const 
    {
        FourierSpec c(*this);
        c += other;
        return c;
    }

    static FourierSpec convolve_freq(const FourierSpec &a, const FourierSpec &b)
    {
        FourierSpec res;
        for(int i = 0; i < SIZE; ++i) {
            for(int j = 0; j < SIZE; ++j) {
                if(i + j < SIZE) {
                    res[i + j] += 0.5f * a[i] * b[j];
                }
                if(i - j >= 0 && i - j < SIZE) {
                    res[i - j] += 0.5f * a[i] * b[j];
                }
            }
        }
        return res;
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
        return convolve_freq(*this, other);
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
        return std::equal(v, v + SIZE, other.v);
    }

    bool operator!=(const FourierSpec &other) const 
    {
        return !(*this == other);
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