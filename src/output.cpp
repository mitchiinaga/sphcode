#include <ctime>

#include <boost/format.hpp>

#include "output.hpp"
#include "logger.hpp"
#include "defines.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "openmp.hpp"

namespace sph
{

inline void output(const SPHParticle & p, std::ofstream &out)
{
#define OUTPUT_SCALAR(name) do { out << p.name << ' '; } while(0)
#define OUTPUT_VECTOR(name) do { for(int i = 0; i < DIM; ++i) out << p.name[i] << ' '; } while(0)

    OUTPUT_VECTOR(pos);
    OUTPUT_VECTOR(vel);
    OUTPUT_VECTOR(acc);
    OUTPUT_SCALAR(mass);
    OUTPUT_SCALAR(dens);
    OUTPUT_SCALAR(pres);
    OUTPUT_SCALAR(ene);
    OUTPUT_SCALAR(sound);
    OUTPUT_SCALAR(sml);
    OUTPUT_SCALAR(id);
    OUTPUT_SCALAR(neighbor);
    out << std::endl;
}

Output::Output(int count)
{
    m_count = count;
    const std::string dir_name = Logger::get_dir_name();
    const std::string file_name = dir_name + "/energy.dat";
    m_out_energy.open(file_name);
    m_out_energy << "# time kinetic thermal potential total\n";
}

Output::~Output()
{
    m_out_energy.close();
}

void Output::output_particle(std::shared_ptr<Simulation> sim)
{
    const auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();
    const real time = sim->get_time();

    const std::string dir_name = Logger::get_dir_name();
    const std::string file_name = dir_name + (boost::format("/%05d.dat") % m_count).str();
    std::ofstream out(file_name);
    out << "# " << time << std::endl;

    for(int i = 0; i < num; ++i) {
        output(particles[i], out);
    }
    WRITE_LOG << "write " << file_name;
    ++m_count;
}

void Output::output_energy(std::shared_ptr<Simulation> sim)
{
    const auto & particles = sim->get_particles();
    const int num = sim->get_particle_num();
    const real time = sim->get_time();

    omp_real kinetic(0.0);
    omp_real thermal(0.0);
    omp_real potential(0.0);

#pragma omp parallel for
    for(int i = 0; i < num; ++i) {
        const auto & p_i = particles[i];
        kinetic.get() += 0.5 * p_i.mass * abs2(p_i.vel);
        thermal.get() += p_i.mass * p_i.ene;
        // potential;
    }

    const real e_k = kinetic.sum();
    const real e_t = thermal.sum();
    const real e_p = potential.sum();
    const real total = e_k + e_t + e_p;

    m_out_energy << time << " " << e_k << " " << e_t << " " << e_p << " " << total << std::endl;
}

}
