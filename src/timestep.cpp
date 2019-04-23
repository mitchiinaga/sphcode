#include "parameters.hpp"
#include "timestep.hpp"
#include "particle.hpp"

namespace sph
{

void TimeStep::initialize(std::shared_ptr<SPHParameters> param)
{
    c_sound = param->cfl.sound;
    c_force = param->cfl.force;
}

real TimeStep::calculation(SPHParticle * particles, int num)
{
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {

    }

    return 0;
}

}