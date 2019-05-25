#pragma once

#include "fluid_force.hpp"

namespace sph
{
namespace gsph
{

class FluidForce : public sph::FluidForce {
    bool m_is_2nd_order;
public:
    void initialize(std::shared_ptr<SPHParameters> param) override;
    void calculation(std::shared_ptr<Simulation> sim) override;
};

}
}