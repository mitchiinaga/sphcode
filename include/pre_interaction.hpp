#pragma once

#include <vector>

#include "module.hpp"

namespace sph
{
class SPHParticle;
class Periodic;

class PreInteraction : public Module {
    bool     m_use_balsara_switch;
    bool     m_use_time_dependent_av;
    real     m_gamma;
    int      m_neighbor_number;
    real     m_kernel_ratio;

    int exhaustive_search(
        SPHParticle & p_i,
        const real kernel_size,
        const std::vector<SPHParticle> & particles,
        const int num,
        std::vector<int> & neighbor_list,
        const int list_size,
        Periodic const * periodic
    ); // for debug

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
    void initial_smoothing(std::shared_ptr<Simulation> sim);
};
}
