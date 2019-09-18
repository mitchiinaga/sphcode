#include "parameters.hpp"
#include "timestep.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "openmp.hpp"

#include <algorithm>
#include <cassert>

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

    const real dt = dt_min.min();
    sim->set_dt(dt);
    sim->set_max_dt(dt);
}
}

namespace indivisual
{

void Timestep::calculation(std::shared_ptr<Simulation> sim)
{
    const auto timeid = sim->get_timeid();
    if(timeid == 1) {
        auto & particles = sim->get_particles();
        auto & dts = sim->get_scalar_array("dt");
        const int num = sim->get_particle_num();

        omp_real _dt_min(std::numeric_limits<real>::max());
        omp_real _dt_max(std::numeric_limits<real>::min());
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

            if(dt_i < _dt_min.get()) {
                _dt_min.get() = dt_i;
            }
            if(dt_i > _dt_max.get()) {
                _dt_max.get() = dt_i;
            }

            dts[i] = dt_i;
        }

        const real dt_min = _dt_min.min();
        const real dt_max = _dt_max.max();

        int level = 0;
        std::vector<real> dt_list;
        dt_list.push_back(dt_min);
        while(dt_list.back() * 2 <= dt_max) {
            level++;
            dt_list.push_back(dt_list.back() * 2);
        }

        // TODO: ‚È‚ñ‚Æ‚©‚·‚é1
        auto pow2 = [](int n) {
            int a = 1;
            for(int i = 0; i < n; ++i) {
                a <<= 1;
            }
            return a;
        };

        // timeids‚Ìì¬
        int id = pow2(level);
        if(timeids.size() != id) {
            timeids.resize(id);
            int interval = 2;
            int loop = 0;

            while(id > 0) {
                int i = pow2(loop) - 1;
                while(i < timeids.size()) {
                    timeids[i] = id;
                    i += interval;
                }

                loop++;
                id >>= 1;
                interval *= 2;
            }
        }

        int max_id = 1;
        for(int i = 0; i < level; ++i) {
            max_id <<= 1;
            ++max_id;
        }

        // TODO: ‚È‚ñ‚Æ‚©‚·‚é2
        auto calc_timeid = [&](const real dt) {
            int id_i = max_id;
            int i = 0;
            while(dt_list[i] * 2 < dt) {
                i++;
                id_i >>= 1;
                assert(id_i > 0);
            }
            return id_i;
        };

#pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < num; ++i) {
            auto & p_i = particles[i];
            const real dt_i = dts[i];
            p_i.timeid = calc_timeid(dt_i);
        }

        current = 0;
        sim->set_timeid(timeids[current]);
        sim->set_dt(dt_min);
        sim->set_max_dt(dt_min * pow2(level));
    } else {
        current++;
        sim->set_timeid(timeids[current]);
    }
}
}

}