#include <ctime>

#include <boost/format.hpp>

#include "output.hpp"
#include "logger.hpp"
#include "defines.hpp"
#include "particle.hpp"
#include "simulation.hpp"

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
    OUTPUT_SCALAR(sml);
    OUTPUT_SCALAR(id);
    OUTPUT_SCALAR(neighbor);
    OUTPUT_SCALAR(alpha);
    OUTPUT_SCALAR(gradh);
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

    real kinetic = 0.0;
    real thermal = 0.0;
    real potential = 0.0;

#pragma omp parallel for reduction(+: kinetic, thermal, potential)
    for(int i = 0; i < num; ++i) {
        const auto & p_i = particles[i];
        kinetic += 0.5 * p_i.mass * abs2(p_i.vel);
        thermal += p_i.mass * p_i.ene;
        potential += 0.5 * p_i.mass * p_i.phi;
    }

    const real e_k = kinetic;
    const real e_t = thermal;
    const real e_p = potential;
    const real total = e_k + e_t + e_p;

    m_out_energy << time << " " << e_k << " " << e_t << " " << e_p << " " << total << std::endl;
}

}
