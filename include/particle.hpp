#pragma once

#include "vector_type.hpp"

namespace sph
{

class SPHParticle {
public:
    vec_t pos;    // position
    vec_t vel;    // velocity
    vec_t vel_p;  // velocity at t + dt/2
    vec_t acc;    // acceleration
    real mass;    // mass
    real dens;    // mass density
    real pres;    // pressure
    real ene;     // internal energy
    real ene_p;   // internal energy at t + dt/2
    real dene;    // du/dt
    real sml;     // smoothing length
    real sound;   // sound speed

    real balsara; // balsara switch
    real alpha;   // AV coefficient

    int id;
    int neighbor;
    SPHParticle *next = nullptr;
};

}