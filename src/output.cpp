#include <ctime>

#include <boost/format.hpp>

#include "output.hpp"
#include "logger.hpp"
#include "defines.hpp"
#include "particle.hpp"

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
    OUTPUT_SCALAR(id);
    out << std::endl;
}

Output::Output(int count)
{
    m_count = count;
    const std::string dir_name = Logger::get_dir_name();
    const std::string file_name = dir_name + "/energy.dat";
    m_out_energy.open(file_name);
}

Output::~Output()
{
    m_out_energy.close();
}

void Output::output_particle(const SPHParticle * particles, const int num, const real time)
{
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

void Output::output_energy(const SPHParticle * particles, const int num, const real time)
{
    // output energy
}

}
