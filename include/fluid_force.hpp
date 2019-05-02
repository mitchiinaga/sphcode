#pragma once

#include <vector>

#include "module.hpp"
#include "vector_type.hpp"

namespace sph
{
class SPHParticle;
class Periodic;

class FluidForce : public Module {
    int      m_neighbor_number;

    int exhaustive_search(
        SPHParticle & p_i,
        const real kernel_size,
        const std::vector<SPHParticle> & particles,
        const int num,
        std::vector<int> & neighbor_list,
        const int list_size,
        Periodic const * periodic
    ); // for debug

    real artificial_viscosity(const SPHParticle & p_i, const SPHParticle & p_j, const vec_t & r_ij);

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};
}