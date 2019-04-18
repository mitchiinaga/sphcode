#pragma once

#include "vector_type.hpp"

namespace sph
{

class SPHParticle {
public:
    vec_t pos;    // position
    vec_t vel;    // velocity
    vec_t vel_i;  // velocity at t + dt/2
    vec_t acc;    // acceleration
    real mass;    // mass
    real dens;    // mass density
    real pres;    // pressure
    real ene;     // internal energy
    real ene_i;   // internal energy at t + dt/2
    real dene;    // du/dt
    real sml;     // smoothing length

    real balsara; // balsara switch
    real alpha;   // AV coefficient

    int id;
    SPHParticle *next = nullptr;
};

}