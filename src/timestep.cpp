#include "parameters.hpp"
#include "timestep.hpp"
#include "particle.hpp"
#include "openmp.hpp"

namespace sph
{

void TimeStep::initialize(std::shared_ptr<SPHParameters> param)
{
    c_sound = param->cfl.sound;
    c_force = param->cfl.force;
}

real TimeStep::calculation(SPHParticle * particles, int num)
{
    omp_real dt_min;
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        const real dt_sound_i = c_sound * particles[i].sml / particles[i].v_sig;
        const real dt_force_i = c_force * std::sqrt(particles[i].sml / abs(particles[i].acc));
        const real dt_min_i = dt_sound_i < dt_force_i ? dt_sound_i : dt_force_i;
        if(dt_min_i < dt_min.get()) {
            dt_min.get() = dt_min_i;
        }
    }

    return dt_min.min();
}

}