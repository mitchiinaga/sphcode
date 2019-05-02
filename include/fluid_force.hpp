#pragma once

#include <vector>

#include "module.hpp"

namespace sph
{
class SPHParticle;
class Distance;

class FluidForce : public Module {
    int      m_neighbor_number;

    int exhaustive_search(
        SPHParticle & p_i,
        const real kernel_size,
        SPHParticle const * particles,
        const int num,
        std::vector<int> & neighbor_list,
        const int list_size,
        Distance const * distance
    ); // for debug

    real artificial_viscosity(const SPHParticle & p_i, const SPHParticle & p_j);

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};
}