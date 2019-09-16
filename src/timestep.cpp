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

namespace fixed
{
void Timestep::calculation(std::shared_ptr<Simulation> sim)
{
    auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();

    omp_real dt_min(std::numeric_limits<real>::max());
#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        const auto & p_i = particles[i];
        real dt_i = c_sound * p_i.sml / p_i.v_sig;

        const real acc_abs = std::abs(particles[i].acc);
        if(acc_abs > 0.0) {
            const real dt_force_i = c_force * std::sqrt(particles[i].sml / acc_abs);
            if(dt_force_i < dt_i) {
                dt_i = dt_force_i;
            }
        }

        if(dt_i < dt_min.get()) {
            dt_min.get() = dt_i;
        }
    }

    sim->set_dt(dt_min.min());
}
}

namespace indivisual
{
void Timestep::calculation(std::shared_ptr<Simulation> sim)
{
    const auto timeid = sim->get_timeid();
    if(timeid == 1) {
//        auto & particles = sim->get_particles();
//        const int num = sim->get_particle_num();
//
//        omp_real dt_min(std::numeric_limits<real>::max());
//#pragma omp parallel for
//        for(int i = 0; i < num; ++i) {
//            const real acc_abs = std::abs(particles[i].acc);
//            if(acc_abs > 0.0) {
//                const real dt_force_i = c_force * std::sqrt(particles[i].sml / acc_abs);
//                if(dt_force_i < dt_min.get()) {
//                    dt_min.get() = dt_force_i;
//                }
//            }
//        }
//
//        const real dt_sound_i = c_sound * sim->get_h_per_v_sig();
//        
//        sim->set_dt(std::min(dt_sound_i, dt_min.min()));
    } else {
        current++;
        sim->set_timeid(timeids[current]);
    }
}
}

}