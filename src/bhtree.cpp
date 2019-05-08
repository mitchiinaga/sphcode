#include <cassert>

#include "parameters.hpp"
#include "bhtree.hpp"
#include "openmp.hpp"
#include "exception.hpp"

namespace sph
{

void BHTree::initialize(std::shared_ptr<SPHParameters> param)
{
    m_max_level         = param->tree.max_level;
    m_leaf_particle_num = param->tree.leaf_particle_num;
    m_root.clear();
    m_root.level = 1;
    m_periodic = param->periodic.is_valid;
    if(m_periodic) {
        m_range_max = param->periodic.range_max;
        m_range_min = param->periodic.range_min;
        m_root.center = (m_range_max + m_range_min) * 0.5;
        auto range = m_range_max - m_range_min;
        real l = 0.0;
        for(int i = 0; i < DIM; ++i) {
            if(l < range[i]) {
                l = range[i];
            }
        }
        m_root.edge = l;
    }
}

void BHTree::resize(const int particle_num, const int tree_size)
{
    assert(m_nodes.get() != nullptr);

    m_node_size = particle_num * tree_size;
    m_nodes = std::shared_ptr<BHNode>(new BHNode[m_node_size], std::default_delete<BHNode[]>());

#pragma omp parallel for
    for(int i = 0; i < m_node_size; ++i) {
        m_nodes.get()[i].clear();
    }
}

void BHTree::make(std::vector<SPHParticle> & particles, const int particle_num)
{
    if(!m_periodic) {
        omp_real r_min[DIM];
        omp_real r_max[DIM];
        for(int i = 0; i < DIM; ++i) {
            r_min[i].get() = std::numeric_limits<real>::max();
            r_max[i].get() = std::numeric_limits<real>::lowest();
        }

#pragma omp parallel for
        for(int i = 0; i < particle_num; ++i) {
            auto & r_i = particles[i].pos;
            for(int j = 0; j < DIM; ++j) {
                if(r_min[j].get() > r_i[j]) {
                    r_min[j].get() = r_i[j];
                }
                if(r_max[j].get() < r_i[j]) {
                    r_max[j].get() = r_i[j];
                }
            }
        }

        vec_t range_min, range_max;
        for(int i = 0; i < DIM; ++i) {
            range_min[i] = r_min[i].min();
            range_max[i] = r_max[i].max();
        }

        m_root.center = (range_max + range_min) * 0.5;
        auto range = range_max - range_min;
        real l = 0.0;
        for(int i = 0; i < DIM; ++i) {
            if(l < range[i]) {
                l = range[i];
            }
        }
        m_root.edge = l;
    }

#pragma omp parallel for
    for(int i = 0; i < particle_num - 1; ++i) {
        particles[i].next = &particles[i + 1];
    }
    particles[particle_num - 1].next = nullptr;
    m_root.first = &particles[0];

    int remaind = m_node_size;
    auto * p = m_nodes.get();
    m_root.create_tree(p, remaind, m_max_level, m_leaf_particle_num);
}

void BHTree::BHNode::create_tree(BHNode * & nodes, int & remaind, const int max_level, const int leaf_particle_num)
{
    for(int i = 0; i < NCHILD; ++i) {
        childs[i] = nullptr;
    }

    auto * pp = first;
    do {
        auto * pnext = pp->next;
        assign(pp, nodes, remaind);
        pp = pnext;
    } while(pp != nullptr);

    for(int i = 0; i < NCHILD; ++i) {
        auto * child = childs[i];
        if(child) {
            child->m_center /= child->mass;

            if(child->num > leaf_particle_num && level < max_level) {
                child->create_tree(nodes, remaind, max_level, leaf_particle_num);
            }
        }
    }
}

void BHTree::BHNode::assign(SPHParticle * particle, BHNode * & nodes, int & remaind)
{
    auto & p_i = *particle;
    const auto & pos = p_i.pos;

    int index = 0;
    int mask = 1;
    for(int i = 0; i < DIM; ++i) {
        if(pos[i] > center[i]) {
            index |= mask;
        }
        mask <<= 1;
    }

    auto * child = childs[index];
    if(!child) {
        if(remaind < 0) {
            THROW_ERROR("no more free node.");
        }
        child = nodes;
        ++nodes;
        --remaind;
        child->clear();
        child->level = level + 1;
        child->edge = edge * 0.5;

        int a = 1;
        real b = 2.0;
        for(int i = 0; i < DIM; ++i) {
            child->center[i] = center[i] + ((index & a) * b - 1.0) * edge * 0.25;
            a <<= 1;
            b *= 0.5;
        }
    }

    child->num++;
    child->mass += p_i.mass;
    child->m_center += pos * p_i.mass;
    p_i.next = child->first;
    child->first = particle;
}

}