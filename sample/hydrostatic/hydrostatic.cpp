#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

namespace sph
{

void Solver::make_hydrostatic()
{
    if(DIM != 2) {
        THROW_ERROR("DIM != 2");
    }

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const real dx1 = 0.5 / N;
    const real dx2 = dx1 * 2.0;
    const real mass = 1.0 / (N * N);
    const real gamma = m_param->physics.gamma;
    int i = 0;

    std::vector<SPHParticle> p;

    // dens region
    real x = -0.25 + dx1 * 0.5;
    real y = -0.25 + dx1 * 0.5;
    while(y < 0.25) {
        SPHParticle p_i;
        p_i.pos[0] = x;
        p_i.pos[1] = y;
        p_i.mass = mass;
        p_i.dens = 4.0;
        p_i.pres = 2.5;
        p_i.ene = p_i.pres / ((gamma - 1.0) * p_i.dens);
        p_i.id = i;
        p.push_back(p_i);

        x += dx1;
        if(x > 0.25) {
            x = -0.25 + dx1 * 0.5;
            y += dx1;
        }
        ++i;
    }

    // ambient
    x = -0.5 + dx2 * 0.5;
    y = -0.5 + dx2 * 0.5;
    while(y < 0.5) {
        SPHParticle p_i;
        p_i.pos[0] = x;
        p_i.pos[1] = y;
        p_i.mass = mass;
        p_i.dens = 1.0;
        p_i.pres = 2.5;
        p_i.ene = p_i.pres / ((gamma - 1.0) * p_i.dens);
        p_i.id = i;
        p.push_back(p_i);

        do {
            x += dx2;
            if(x > 0.5) {
                x = -0.5 + dx2 * 0.5;
                y += dx2;
            }
        } while(x > -0.25 && x < 0.25 && y > -0.25 && y < 0.25);
        ++i;
    }

    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
}

}
