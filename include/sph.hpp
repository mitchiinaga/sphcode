#pragma once

#include <memory>

#include "defines.hpp"

namespace sph
{

class SPHParticle;
struct SPHParameters;

class SPH {
    std::unique_ptr<SPHParticle> m_particles;
    int m_particle_num;
public:
    SPH(std::shared_ptr<SPHParameters> & param);
    void integrate(real * time);
    int get_particle_num();
    const SPHParticle * get_particles();
};

}