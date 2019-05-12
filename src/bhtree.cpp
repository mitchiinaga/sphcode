#include <cassert>

#include "parameters.hpp"
#include "bhtree.hpp"
#include "openmp.hpp"
#include "exception.hpp"
#include "periodic.hpp"

namespace sph
{

void BHTree::initialize(std::shared_ptr<SPHParameters> param)
{
    m_max_level         = param->tree.max_level;
    m_leaf_particle_num = param->tree.leaf_particle_num;
    m_root.clear();
    m_root.level = 1;
    m_is_periodic = param->periodic.is_valid;
    if(m_is_periodic) {
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
    m_periodic = std::make_shared<Periodic>();
    m_periodic->initialize(param);
}

void BHTree::resize(const int particle_num, const int tree_size)
{
    assert(m_nodes.get() == nullptr);

    m_node_size = particle_num * tree_size;
    m_nodes = std::shared_ptr<BHNode>(new BHNode[m_node_size], std::default_delete<BHNode[]>());

#pragma omp parallel for
    for(int i = 0; i < m_node_size; ++i) {
        m_nodes.get()[i].clear();
    }
}

void BHTree::make(std::vector<SPHParticle> & particles, const int particle_num)
{
    if(!m_is_periodic) {
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

void BHTree::set_kernel()
{
    m_root.set_kernel();
}

int BHTree::neighbor_search(const SPHParticle & p_i, std::vector<int> & neighbor_list, const bool is_ij)
{
    int n_neighbor = 0;
    m_root.neighbor_search(p_i, neighbor_list, n_neighbor, is_ij, m_periodic.get());
    return n_neighbor;
}

// --------------------------------------------------------------------------------------------------------------- //

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

    int num_child = 0;
    for(int i = 0; i < NCHILD; ++i) {
        auto * child = childs[i];
        if(child) {
            ++num_child;
            child->m_center /= child->mass;

            if(child->num > leaf_particle_num && level < max_level) {
                child->create_tree(nodes, remaind, max_level, leaf_particle_num);
            }
        }
    }

    if(!num_child) {
        is_leaf = true;
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
        childs[index] = nodes;
        child = childs[index];
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

real BHTree::BHNode::set_kernel()
{
    real kernel = 0.0;
    if(is_leaf) {
        auto * p = first;
        while(p) {
            const real h = p->sml;
            if(h > kernel) {
                kernel = h;
            }
            p = p->next;
        }
    } else {
        for(int i = 0; i < NCHILD; ++i) {
            auto * child = childs[i];
            const real h = child->set_kernel();
            if(h > kernel) {
                kernel = h;
            }
        }
    }

    kernel_size = kernel;
    return kernel_size;
}

void BHTree::BHNode::neighbor_search(const SPHParticle & p_i, std::vector<int> & neighbor_list, int & n_neighbor, const bool is_ij, const Periodic * periodic)
{
    const vec_t & r_i = p_i.pos;
    const real h = is_ij ? std::max(p_i.sml, kernel_size) : p_i.sml;
    const real h2 = h * h;
    const real l2 = sqr(edge + h);
    const vec_t d = periodic->calc_r_ij(r_i, center);
    real dx2_max = sqr(d[0]);
    for(int i = 1; i < DIM; ++i) {
        const real dx2 = sqr(d[i]);
        if(dx2 > dx2_max) {
            dx2_max = dx2;
        }
    }

    if(dx2_max <= l2) {
        if(is_leaf) {
            auto * p = first;
            while(p) {
                const vec_t & r_j = p->pos;
                const vec_t r_ij = periodic->calc_r_ij(r_i, r_j);
                const real r2 = abs2(r_ij);
                if(r2 < h2) {
                    neighbor_list[n_neighbor] = p->id;
                    ++n_neighbor;
                }
                p = p->next;
            }
        } else {
            for(int i = 0; i < NCHILD; ++i) {
                if(childs[i]) {
                    childs[i]->neighbor_search(p_i, neighbor_list, n_neighbor, is_ij, periodic);
                }
            }
        }
    }
}

}