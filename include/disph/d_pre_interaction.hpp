#pragma once

#include "pre_interaction.hpp"

namespace sph
{
namespace disph
{

class PreInteraction : public sph::PreInteraction {
    real newton_raphson (
        const SPHParticle & p_i,
        const std::vector<SPHParticle> & particles,
        const std::vector<int> & neighbor_list,
        const int n_neighbor,
        const Periodic * periodic,
        const KernelFunction * kernel
    ) override;

public:
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}
