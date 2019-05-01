#pragma once

#include "module.hpp"

namespace sph
{
class SPHParticle;
class Distance;

class PreInteraction : public Module {
    bool     m_use_balsara_switch;
    bool     m_use_time_dependent_av;
    real     m_gamma;
    int      m_neighbor_number;

    int exhaustive_search(
        SPHParticle & p_i,
        const real kernel_size,
        SPHParticle const * particles,
        const int num,
        std::shared_ptr<int[]> neighbor_list,
        const int list_size,
        Distance const * distance
    ); // for debug

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};
}
