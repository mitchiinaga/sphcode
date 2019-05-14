#pragma once

#include "module.hpp"

namespace sph
{
class SPHParticle;
class Periodic;

class PreInteraction : public Module {
    bool m_use_balsara_switch;
    bool m_use_time_dependent_av;
    real m_alpha_max;
    real m_alpha_min;
    real m_epsilon; // tau = h / (epsilon * c)
    real m_gamma;
    int  m_neighbor_number;
    real m_kernel_ratio;

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
    void initial_smoothing(std::shared_ptr<Simulation> sim);
};
}
