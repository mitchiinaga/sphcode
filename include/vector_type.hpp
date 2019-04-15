#pragma once

#include <cmath>

#include "defines.hpp"

class vec_t {
    real vec[DIM];
public:
    // Constructor
#if DIM == 1
    vec_t(const real x = 0) {
        vec[0] = x;
    }
#elif DIM == 2
    vec_t(const real x = 0, const real y = 0) {
        vec[0] = x;
        vec[1] = y;
    }
#elif DIM == 3
    vec_t(const real x = 0, const real y = 0, const real z = 0) {
        vec[0] = x;
        vec[1] = y;
        vec[2] = z;
    }
#endif

    vec_t(const vec_t & a) {
        for(int i = 0; i < DIM; ++i) vec[i] = a[i];
    }

    vec_t(const real (&a)[DIM]) {
        for(int i = 0; i < DIM; ++i) vec[i] = a[i];
    }

    // Operator
    real & operator[](const int i) { return vec[i]; }
    const real & operator[](const int i) const { return vec[i]; }

    vec_t & operator=(const vec_t &a) {
        for(int i = 0; i < DIM; ++i) vec[i] = a[i];
        return *this;
    }

    vec_t & operator=(const real (&a)[DIM]) {
        for(int i = 0; i < DIM; ++i) vec[i] = a[i];
        return *this;
    }

    vec_t & operator=(const real a) {
        for(int i = 0; i < DIM; ++i) vec[i] = a;
        return *this;
    }

    const vec_t & operator+() const { return *this; }

    const vec_t operator-() const {
#if DIM == 1
        return vec_t(-vec[0]);
#elif DIM == 2
        return vec_t(-vec[0], -vec[1]);
#elif DIM == 3
        return vec_t(-vec[0], -vec[1], -vec[2]);
#endif
    }

    // +=
    vec_t & operator+=(const vec_t &a) {
        for(int i = 0; i < DIM; ++i) vec[i] += a[i];
        return *this;
    }

    vec_t & operator+=(const real (&a)[DIM]) {
        for(int i = 0; i < DIM; ++i) vec[i] += a[i];
        return *this;
    }

    vec_t & operator+=(const real a) {
        for(int i = 0; i < DIM; ++i) vec[i] += a;
        return *this;
    }

    // -=
    vec_t & operator-=(const vec_t &a) {
        for(int i = 0; i < DIM; ++i) vec[i] -= a[i];
        return *this;
    }

    vec_t & operator-=(const real (&a)[DIM]) {
        for(int i = 0; i < DIM; ++i) vec[i] -= a[i];
        return *this;
    }

    vec_t & operator-=(const real a) {
        for(int i = 0; i < DIM; ++i) vec[i] -= a;
        return *this;
    }

    vec_t & operator*=(const real a) {
        for(int i = 0; i < DIM; ++i) vec[i] *= a;
        return *this;
    }

    vec_t & operator/=(const real a) {
        for(int i = 0; i < DIM; ++i) vec[i] /= a;
        return *this;
    }

    // +
    vec_t operator+(const vec_t &a) const {
#if DIM == 1
        return vec_t(vec[0] + a[0]);
#elif DIM == 2
        return vec_t(vec[0] + a[0], vec[1] + a[1]);
#elif DIM == 3
        return vec_t(vec[0] + a[0], vec[1] + a[1], vec[2] + a[2]);
#endif
    }

    vec_t operator+(const real (&a)[DIM]) const {
#if DIM == 1
        return vec_t(vec[0] + a[0]);
#elif DIM == 2
        return vec_t(vec[0] + a[0], vec[1] + a[1]);
#elif DIM == 3
        return vec_t(vec[0] + a[0], vec[1] + a[1], vec[2] + a[2]);
#endif
    }

    vec_t operator+(const real a) const {
#if DIM == 1
        return vec_t(vec[0] + a);
#elif DIM == 2
        return vec_t(vec[0] + a, vec[1] + a);
#elif DIM == 3
        return vec_t(vec[0] + a, vec[1] + a, vec[2] + a);
#endif
    }

    // -
    vec_t operator-(const vec_t &a) const {
#if DIM == 1
        return vec_t(vec[0] - a[0]);
#elif DIM == 2
        return vec_t(vec[0] - a[0], vec[1] - a[1]);
#elif DIM == 3
        return vec_t(vec[0] - a[0], vec[1] - a[1], vec[2] - a[2]);
#endif
    }

    vec_t operator-(const real (&a)[DIM]) const {
#if DIM == 1
        return vec_t(vec[0] - a[0]);
#elif DIM == 2
        return vec_t(vec[0] - a[0], vec[1] - a[1]);
#elif DIM == 3
        return vec_t(vec[0] - a[0], vec[1] - a[1], vec[2] - a[2]);
#endif
    }

    vec_t operator-(const real a) const {
#if DIM == 1
        return vec_t(vec[0] - a);
#elif DIM == 2
        return vec_t(vec[0] - a, vec[1] - a);
#elif DIM == 3
        return vec_t(vec[0] - a, vec[1] - a, vec[2] - a);
#endif
    }

    vec_t operator*(const real a) const {
#if DIM == 1
        return vec_t(vec[0] * a);
#elif DIM == 2
        return vec_t(vec[0] * a, vec[1] * a);
#elif DIM == 3
        return vec_t(vec[0] * a, vec[1] * a, vec[2] * a);
#endif
    }

    vec_t operator/(const real a) const {
#if DIM == 1
        return vec_t(vec[0] / a);
#elif DIM == 2
        return vec_t(vec[0] / a, vec[1] / a);
#elif DIM == 3
        return vec_t(vec[0] / a, vec[1] / a, vec[2] / a);
#endif
    }

    const real *get_array() const { return vec; }
};

inline real inner_product(const vec_t &a, const vec_t &b)
{
#if DIM == 1
    return a[0] * b[0];
#elif DIM == 2
    return a[0] * b[0] + a[1] * b[1];
#elif DIM == 3
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
#endif
}

inline real inner_product(const vec_t &a, const real (&b)[DIM])
{
#if DIM == 1
    return a[0] * b[0];
#elif DIM == 2
    return a[0] * b[0] + a[1] * b[1];
#elif DIM == 3
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
#endif
}

inline real inner_product(const real (&a)[DIM], const vec_t &b)
{
#if DIM == 1
    return a[0] * b[0];
#elif DIM == 2
    return a[0] * b[0] + a[1] * b[1];
#elif DIM == 3
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
#endif
}

inline real abs2(const vec_t &a)
{
    return inner_product(a, a);
}

inline real abs(const vec_t &a)
{
    return std::sqrt(inner_product(a, a));
}

#if DIM == 2
inline real vector_product(const vec_t &a, const vec_t &b)
{
    return a[0] * b[1] - a[1] * b[0];
}

inline real vector_product(const vec_t &a, const real (&b)[DIM])
{
    return a[0] * b[1] - a[1] * b[0];
}

inline real vector_product(const real (&a)[DIM], const vec_t &b)
{
    return a[0] * b[1] - a[1] * b[0];
}
#elif DIM == 3
inline vec_t vector_product(const vec_t &a, const vec_t &b)
{
    return vec_t(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

inline vec_t vector_product(const vec_t &a, const real (&b)[DIM])
{
    return vec_t(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

inline vec_t vector_product(const real (&a)[DIM], const vec_t &b)
{
    return vec_t(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
#endif

inline real distance(const vec_t &a, const vec_t &b)
{
    return abs(a - b);
}
