#include "solver.hpp"
#include "simulation.hpp"
#include "particle.hpp"
#include "exception.hpp"
#include "parameters.hpp"

namespace sph
{

// Kelvin-Helmholtz instability (Springel 2010)

void Solver::make_khi()
{
#if DIM != 2
    THROW_ERROR("DIM != 2");
#else

    const int N = boost::any_cast<int>(m_sample_parameters["N"]);
    const int num = N * N * 3 / 4;
    const real dx = 1.0 / N;
    const real mass = 1.5 / num;
    const real gamma = m_param->physics.gamma;

    std::vector<SPHParticle> p(num);

    // dens region
    real x = dx * 0.5;
    real y = dx * 0.5;

    int region = 1;
    bool odd = true;

    auto vy = [](const real x, const real y) {
        constexpr real sigma2_inv = 2 / (0.05 * 0.05);
        return 0.1 * std::sin(4.0 * M_PI * x) * (std::exp(-sqr(y - 0.25) * 0.5 * sigma2_inv) + std::exp(-sqr(y - 0.75) * 0.5 * sigma2_inv));
    };

    for(int i = 0; i < num; ++i) {
        auto & p_i = p[i];
        p_i.pos[0] = x;
        p_i.pos[1] = y;
        p_i.vel[0] = region == 1 ? -0.5 : 0.5;
        p_i.vel[1] = vy(x, y);
        p_i.mass = mass;
        p_i.dens = static_cast<real>(region);
        p_i.pres = 2.5;
        p_i.ene = p_i.pres / ((gamma - 1.0) * p_i.dens);
        p_i.id = i;

        x += region == 1 ? 2.0 * dx : dx;
        
        if(x > 1.0) {
            y += dx;
            
            if(y > 0.25 && y < 0.75) {
                region = 2;
            } else {
                region = 1;
            }

            if(region == 1) {
                if(odd) {
                    odd = false;
                    x = dx * 1.5;
                } else {
                    odd = true;
                    x = dx * 0.5;
                }
            } else {
                x = dx * 0.5;
            }
        }
    }

    m_sim->set_particles(p);
    m_sim->set_particle_num(p.size());
#endif
}

}
