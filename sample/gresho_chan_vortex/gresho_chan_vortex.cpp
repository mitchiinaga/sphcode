#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

namespace sph
{

real vortex_velocity(const real r)
{
    if(r < 0.2) {
        return 5.0 * r;
    } else if(r < 0.4) {
        return 2.0 - 5.0 * r;
    } else {
        return 0.0;
    }
}

real vortex_pressure(const real r)
{
    if(r < 0.2) {
        return 5.0 + 12.5 * r * r;
    } else if(r < 0.4) {
        return 9.0 + 12.5 * r * r - 20.0 * r + 4.0 * std::log(5.0 * r);
    } else {
        return 3.0 + 4.0 * std::log(2.0);
    }
}

void Solver::make_gresho_chan_vortex()
{
#if DIM == 2
    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real dx = 1.0 / N;

    const int num = N * N;
    std::vector<SPHParticle> p(num);

    real x = -0.5 + dx * 0.5;
    real y = -0.5 + dx * 0.5;
    const real mass = 1.0 / num;
    const real gamma = m_param->physics.gamma;

    for(int i = 0; i < num; ++i) {
        auto & p_i = p[i];

        p_i.pos[0] = x;
        p_i.pos[1] = y;
        const real r = std::abs(p_i.pos);
        const real vel = vortex_velocity(r);
        vec_t dir(-y, x);
        dir /= r;
        p_i.vel = dir * vel;
        p_i.dens = 1.0;
        p_i.pres = vortex_pressure(r);
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
#else
    THROW_ERROR("DIM != 2");
#endif
}

}
