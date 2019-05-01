#pragma once

#include <memory>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "defines.hpp"
#include "timestep.hpp"
#include "pre_interaction.hpp"
#include "fluid_force.hpp"

namespace sph
{

struct SPHParameters;
class Simulation;
class Output;

class Solver {
    std::shared_ptr<SPHParameters> m_param;
    std::shared_ptr<Output>        m_output;
    std::string                    m_output_dir;
    std::shared_ptr<Simulation>    m_sim;

    // modules
    TimeStep       m_timestep;
    PreInteraction m_pre;
    FluidForce     m_fforce;

    void read_parameterfile(const char * filename);
    void initialize();
    void predict();
    void correct();
    void integrate();

public:
    Solver(int argc, char * argv[]);
    void run();
};

}