#pragma once

#include <string>
#include <memory>

#include "defines.hpp"
#include "vector_type.hpp"

namespace sph
{
struct SPHParameters;
class Output;

class Core {
    std::shared_ptr<SPHParameters> m_param;
    std::shared_ptr<Output>        m_output;
    std::string                    m_output_dir;

    void read_parameterfile(const char * filename);
public:
    Core(int argc, char * argv[]);
    void run();
};

}
