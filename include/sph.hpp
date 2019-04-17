#pragma once

#include <memory>

#include "defines.hpp"

namespace sph
{

struct SPHParameters;

class SPH {
public:
    SPH(std::shared_ptr<SPHParameters> & param);
    void output_particles(const real time);
    void output_energy(const real time);
    void integrate(real * time);
    int get_particle_num();
};

}