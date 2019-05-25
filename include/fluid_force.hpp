#pragma once

#include "module.hpp"
#include "vector_type.hpp"

namespace sph
{
class SPHParticle;
class Periodic;

class FluidForce : public Module {
protected:
    int  m_neighbor_number;
    bool m_use_ac;
    real m_alpha_ac;
    bool m_use_gravity;

    real artificial_viscosity(const SPHParticle & p_i, const SPHParticle & p_j, const vec_t & r_ij);
    real artificial_conductivity(const SPHParticle & p_i, const SPHParticle & p_j, const vec_t & r_ij, const vec_t & dw_ij);

public:
    virtual void initialize(std::shared_ptr<SPHParameters> param) override;
    virtual void calculation(std::shared_ptr<Simulation> sim) override;
};
}