#include <random>

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

namespace sph
{

// Evrard collapse (Evrard 1988)
void Solver::make_evrard()
{
#if DIM != 3
    THROW_ERROR("DIM != 3");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    auto & p = m_sim->get_particles();
    const real dx = 2.0 / N;

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            for(int k = 0; k < N; ++k) {
                vec_t r = {
                    (i + 0.5) * dx - 1.0,
                    (j + 0.5) * dx - 1.0,
                    (k + 0.5) * dx - 1.0
                };
                const real r_0 = std::abs(r);
                if(r_0 > 1.0) {
                    continue;
                }

                if(r_0 > 0.0) {
                    const real r_abs = std::pow(r_0, 1.5);
                    r *= r_abs / r_0;
                }

                SPHParticle p_i;
                p_i.pos = r;
                p.emplace_back(p_i);
            }
        }
    }

    const real mass = 1.0 / p.size();
    const real gamma = m_param->physics.gamma;
    const real G = m_param->gravity.constant;
    const real u = 0.05 * G;

    int i = 0;
    for(auto & p_i : p) {
        p_i.vel = 0.0;
        p_i.vel = 0.0;
        p_i.mass = mass;
        p_i.dens = 1.0 / (2.0 * M_PI * std::abs(p_i.pos));
        p_i.ene = u;
        p_i.pres = (gamma - 1.0) * p_i.dens * u;
        p_i.id = i;
        i++;
    }

    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

}
