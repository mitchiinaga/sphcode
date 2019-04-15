#pragma once

#include <string>
#include <memory>

#include "defines.hpp"
#include "vector_type.hpp"

namespace sph
{
struct SPHParameters;

class Core {
    std::shared_ptr<SPHParameters> param;
    std::string output_dir;

    void read_parameterfile(const char * filename);
public:
    Core(int argc, char * argv[]);
    void run();
};

}
