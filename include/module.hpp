#pragma once

#include <memory>

namespace sph
{
struct SPHParameters;
struct Simulation;

class Module {
public:
    virtual void initialize(std::shared_ptr<SPHParameters> param) = 0;
    virtual void calculation(std::shared_ptr<Simulation> sim) = 0;
};
}
