#pragma once

#include <vector>

#include "module.hpp"
#include "particle.hpp"

namespace sph
{
class Periodic;
class KernelFunction;

class PreInteraction : public Module {
protected:
    bool m_use_balsara_switch;
    bool m_use_time_dependent_av;
    real m_alpha_max;
    real m_alpha_min;
    real m_epsilon; // tau = h / (epsilon * c)
    real m_gamma;
    int  m_neighbor_number;
    real m_kernel_ratio;
    bool m_iteration;
    bool m_first;

    virtual real newton_raphson(
        const SPHParticle & p_i,
        const std::vector<SPHParticle> & particles,
        const std::vector<int> & neighbor_list,
        const int n_neighbor,
        const Periodic * periodic,
        const KernelFunction * kernel
    );
    void initial_smoothing(std::shared_ptr<Simulation> sim);

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    virtual void calculation(std::shared_ptr<Simulation> sim) override;
};
}
