#include <cassert>

#include "sph.hpp"
#include "parameters.hpp"
#include "particle.hpp"

namespace sph
{

SPH::SPH(std::shared_ptr<SPHParameters> & param)
{
}

void SPH::integrate(real * time)
{
    // dt = calc_dt();
    // predict(dt);
    // make_tree();
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