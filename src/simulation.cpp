#include "parameters.hpp"
#include "simulation.hpp"
#include "exception.hpp"
#include "periodic.hpp"
#include "bhtree.hpp"
#include "kernel/cubic_spline.hpp"
#include "kernel/wendland_kernel.hpp"

namespace sph
{

Simulation::Simulation(std::shared_ptr<SPHParameters> param)
{
    if(param->kernel == KernelType::CUBIC_SPLINE) {
        m_kernel = std::make_shared<Spline::Cubic>();
    } else if(param->kernel == KernelType::WENDLAND) {
        m_kernel = std::make_shared<Wendland::C4Kernel>();
    } else {
        THROW_ERROR("kernel is unknown.");
    }

    m_periodic = std::make_shared<Periodic>();
    m_periodic->initialize(param);

    m_tree = std::make_shared<BHTree>();
    m_tree->initialize(param);

    m_time = param->time.start;
}

void Simulation::update_time()
{
    m_time += m_dt;
}

void Simulation::make_tree()
{
    m_tree->make(m_particles, m_particle_num);
}

}