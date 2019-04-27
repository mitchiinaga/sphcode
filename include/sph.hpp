#pragma once

#include <memory>
#include <unordered_map>

#include "defines.hpp"

#include "timestep.hpp"
#include "pre_interaction.hpp"

namespace sph
{

class SPHParticle;
struct SPHParameters;

class SPH {
    std::unique_ptr<SPHParticle>   m_particles;
    int                            m_particle_num;

    // modules
    TimeStep       m_timestep;
    PreInteraction m_pre;

    void initialize();
    void predict(const real dt);
    void correct(const real dt);
public:
    SPH(std::shared_ptr<SPHParameters> param);
    int get_particle_num();
    const SPHParticle * get_particles();
    void integrate(real * time);
};

}