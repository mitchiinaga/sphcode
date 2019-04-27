#pragma once

#include <memory>

namespace sph
{
struct SPHParameters;

class Module {
public:
    virtual void initialize(std::shared_ptr<SPHParameters> param) = 0;
};
}
