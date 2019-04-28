#pragma once

#include <iostream>
#include <string>
#include <cstdlib>

#define DIM 1
typedef double real;

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751
#endif

inline real inner_product(const real (&v1)[DIM], const real (&v2)[DIM])
{
#if DIM == 1
    return v1[0] * v2[0];
#elif DIM == 2
    return v1[0] * v2[0] + v1[1] * v2[1];
#elif DIM == 3
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
#endif
}

inline real sqr(real x) { return x * x; }
inline real pow3(real x) { return x * x * x; }
inline real pow4(real x) { return x * x * x * x; }
inline real pow5(real x) { return x * x * x * x * x; }
inline real pow6(real x) { return x * x * x * x * x * x; }
