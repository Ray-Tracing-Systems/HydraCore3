#include <algorithm>
#include <numeric>

template<typename T, uint N>
struct vec
{
    static constexpr uint length = N;

    T v[N];

    vec() : v{} { std::fill_n(v, N, T(0)); }
    vec(const vec &other) = default;
    vec(T x)  : v{} { std::fill_n(v, N, x); }

    vec &operator=(const vec &other) = default;

    vec operator+(const vec &other) const 
    {
        vec c(*this);
        c += v;
        return c;
    }

    vec operator-(const vec &other) const 
    {
        vec c(*this);
        c -= v;
        return c;
    }

    vec operator*(T f) const 
    {
        vec c(*this);
        c *= f;
        return c;
    }

    vec operator/(T f) const 
    {
        vec c(*this);
        c /= f;
        return c;
    }

    vec operator*(const vec &other) const 
    {
        vec c(*this);
        c *= v;
        return c;
    }

    vec operator/(const vec &other) const 
    {
        vec c(*this);
        c /= v;
        return c;
    }

    T &operator[](int i) 
    {
        return v[i];
    }

    T operator[](int i) const 
    {
        return v[i];
    }

    vec operator-() const 
    {
        vec res;
        std::transform(v, v + N, [](T val) { return -val; });
        return res;
    }

    vec &operator+=(const vec &other) 
    {
        for(uint i = 0; i < N; ++i) {
            v[i] += other.v[i];
        }
        return *this;
    }

    vec &operator-=(const vec &other) 
    {
        for(uint i = 0; i < N; ++i) {
            v[i] -= other.v[i];
        }
        return *this;
    }

    vec &operator*=(T f) 
    {
        for(uint i = 0; i < N; ++i) {
            v[i] *= f;
        }
        return *this;
    }

    vec &operator/=(T f) 
    {
        const T invf = T(1.0) / f;
        for(uint i = 0; i < N; ++i) {
            v[i] *= invf;
        }
        return *this;
    }

    vec &operator*=(const vec &other) 
    {
        for(uint i = 0; i < N; ++i) {
            v[i] *= other.v[i];
        }
        return *this;
    }

    vec &operator/=(const vec &other) 
    {
        for(uint i = 0; i < N; ++i) {
            v[i] /= other.v[i];
        }
        return *this;
    }

    bool operator==(const vec &other) const 
    {
        return std::equal(v, v + N, other.v);
    }

    bool operator!=(const vec &other) const 
    {
        return !(*this == v);
    }


    bool operator<(const vec &other) const 
    {
        return std::lexicographical_compare(v, v + N, other.v, other.v + N);
    }

    bool operator<=(const vec &other) const 
    {
        return (*this == v) || (*this < v);
    }

    bool operator>(const vec &other) const 
    {
        return !(*this <= v);
    }

    bool operator>=(const vec &other) const 
    {
        return !(*this < v);
    }

    T sum() const 
    {
        return std::accumulate(v, v + N, T(0));
    }


    uint argmax() const 
    {
        const T *ptr = std::max_element(v, v + N);
        return uint(ptr - v);
    }

    uint argmin() const 
    {
        const T *ptr = std::min_element(v, v + N);
        return uint(ptr - v);
    }

    T max() const 
    {
        return v[argmax()];
    }

    T min() const 
    {
        return v[argmin()];
    }

    template<typename P>
    vec<P, N> cast() const
    {
        if constexpr(std::is_same_v<P, T>) {
            return *this;
        }
        else {
            vec<P, N> res;
            std::copy(v, v + N, res.v);
            return res;
        }
    }

    static T distance2(const vec &v1, const vec &v2) 
    {
        double res = 0.0;

        for(uint i = 0; i < N; ++i) {
            T dvi = v1.v[i] - v2.v[i];
            res += dvi * dvi;
        }

        return T(res);
    }

    static T distance(const vec &v1, const vec &v2)
    {
        return std::sqrt(distance2(v1, v2));
    }



};