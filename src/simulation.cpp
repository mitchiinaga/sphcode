#include "parameters.hpp"
#include "simulation.hpp"
#include "exception.hpp"
#include "distance.hpp"
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

    m_distance = std::make_shared<Distance>();
    m_distance->initialize(param);

    m_time = param->time.start;
}

void Simulation::update_time()
{
    m_time += m_dt;
}

}