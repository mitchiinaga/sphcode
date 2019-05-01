#pragma once

#include "vector_type.hpp"

namespace sph
{

inline real powh(const real h) {
#if DIM == 1
    return h;
#elif DIM == 2
    return h * h;
#elif DIM == 3
    return h * h * h;
#endif
}

class KernelFunction {
public:
    virtual real w(const real r, const real h) = 0;                     // W(r,h)
    virtual vec_t dw(const vec_t &rij, const real r, const real h) = 0; // grad W(r,h)
    virtual real dhw(const real r, const real h) = 0;                   // dW(r,h)/dh
};

}
