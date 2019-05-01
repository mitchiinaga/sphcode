#include "parameters.hpp"
#include "timestep.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "openmp.hpp"

#include <algorithm>

namespace sph
{

void TimeStep::initialize(std::shared_ptr<SPHParameters> param)
{
    c_sound = param->cfl.sound;
    c_force = param->cfl.force;
}

void TimeStep::calculation(std::shared_ptr<Simulation> sim)
{
    auto * particles = sim->get_particles().get();
    auto * distance = sim->get_distance().get();
    const int num = sim->get_particle_num();

    omp_real dt_min;
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        const real dt_force_i = c_force * std::sqrt(particles[i].sml / abs(particles[i].acc));
        if(dt_force_i < dt_min.get()) {
            dt_min.get() = dt_force_i;
        }
    }

    const real dt_sound_i = c_sound * sim->get_h_per_v_sig_max();
    
    sim->set_dt(std::min(dt_sound_i, dt_min.min()));
}

}