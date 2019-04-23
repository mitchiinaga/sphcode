#include <cassert>

#include "sph.hpp"
#include "parameters.hpp"
#include "particle.hpp"

namespace sph
{

SPH::SPH(std::shared_ptr<SPHParameters> param)
{
    m_timestep.initialize(param);
}

void SPH::initialize()
{
    // calc_tree();
    // pre_interaction();
    // calc_force();
}

void SPH::integrate(real * time)
{
    real const dt = m_timestep.calculation(m_particles.get(), m_particle_num);
    // predict(dt);
    // calc_tree();
    // pre_interaction();
    // calc_force();
    // correct(dt);
}

int SPH::get_particle_num()
{
    return m_particle_num;
}

const SPHParticle * SPH::get_particles()
{
    const SPHParticle * p = m_particles.get();
    assert(p);
    return p;
}

}