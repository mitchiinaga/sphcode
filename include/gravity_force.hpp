#pragma once

#include "module.hpp"
#include "vector_type.hpp"

namespace sph
{
class SPHParticle;

class GravityForce : public Module {
    bool  m_is_valid;
    real m_constant;

public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};
}
