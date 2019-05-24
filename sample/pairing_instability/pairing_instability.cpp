#include <random>

#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

namespace sph
{

void Solver::make_pairing_instability()
{
#if DIM != 2
    THROW_ERROR("DIM != 2");
#else
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real dx = 1.0 / N;

    const int num = N * N;
    std::vector<SPHParticle> p(num);

    real x = -0.5 + dx * 0.5;
    real y = -0.5 + dx * 0.5;
    const real mass = 1.0 / num;
    const real gamma = m_param->physics.gamma;

    std::mt19937 engine(1);
    std::uniform_real_distribution<real> dist(-dx * 0.05, dx * 0.05);

    for(int i = 0; i < num; ++i) {
        auto & p_i = p[i];

        p_i.pos[0] = x + dist(engine);
        p_i.pos[1] = y + dist(engine);
        p_i.vel = 0.0;
        p_i.dens = 1.0;
        p_i.pres = 1.0;
        p_i.mass = mass;
        p_i.ene = p_i.pres / ((gamma - 1.0) * p_i.dens);
        p_i.id = i;

        x += dx;
        if(x > 0.5) {
            x = -0.5 + dx * 0.5;
            y += dx;
        }
    }

    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

}
