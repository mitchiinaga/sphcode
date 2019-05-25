#pragma once

#include "fluid_force.hpp"

namespace sph
{
namespace disph
{

class FluidForce : public sph::FluidForce {
    real m_gamma;
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}