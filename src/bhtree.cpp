#include "parameters.hpp"
#include "bhtree.hpp"

namespace sph
{

void BHTree::initialize(std::shared_ptr<SPHParameters> param)
{
    m_max_level         = param->tree.max_level;
    m_leaf_particle_num = param->tree.leaf_particle_num;
    m_periodic = param->periodic.is_valid;
    if(m_periodic) {
        m_range_max = param->periodic.range_max;
        m_range_min = param->periodic.range_min;
    }
}

void BHTree::assign(const int particle_num, const int tree_size)
{
    m_nodes.resize(particle_num * tree_size);

    const int n = m_nodes.size();
#pragma omp parallel for
    for(int i = 0; i < n; ++i) {
        m_nodes[i].clear();
    }
}

void BHTree::make(std::vector<SPHParticle> & particles, const int particle_num)
{
    if(!m_periodic) {
        // m_range_max, min‚ÌŒvŽZ
    }

    // center

#pragma omp parallel for
    for(int i = 0; i < particle_num - 1; ++i) {
        particles[i].next = &particles[i + 1];
    }
    particles[particle_num - 1].next = nullptr;
}

}