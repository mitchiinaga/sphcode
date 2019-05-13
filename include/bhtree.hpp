#pragma once

#include <memory>
#include <vector>
#include <algorithm>

#include "vector_type.hpp"
#include "particle.hpp"

namespace sph
{
struct SPHParameters;
class Periodic;

constexpr int NCHILD = DIM == 1 ? 2 : DIM == 2 ? 4 : 8;

class BHTree
{
    class BHNode
    {
    public:
        SPHParticle * first;
        real mass;
        int num;
        BHNode * childs[NCHILD];
        vec_t center;
        vec_t m_center; // mass center
        real edge;
        int level;
        real kernel_size;
        bool is_leaf;

        void clear() {
            first = nullptr;
            mass = 0.0;
            num = 0;
            std::fill(childs, childs + NCHILD, nullptr);
            center = 0.0;
            m_center = 0.0;
            edge = 0.0;
            level = 0;
            kernel_size = 0.0;
            is_leaf = false;
        }

        void root_clear() {
            first = nullptr;
            mass = 0.0;
            num = 0;
            m_center = 0.0;
            kernel_size = 0.0;
            is_leaf = false;
        }

        void create_tree(BHNode * & nodes, int & remaind, const int max_level, const int leaf_particle_num);
        void assign(SPHParticle * particle, BHNode * & nodes, int & remaind);
        real set_kernel();
        void neighbor_search(const SPHParticle & p_i, std::vector<int> & neighbor_list, int & n_neighbor, const bool is_ij, const Periodic * periodic);
    };

    int  m_max_level;
    int  m_leaf_particle_num;
    bool m_is_periodic;
    vec_t m_range_max;
    vec_t m_range_min;
    std::shared_ptr<Periodic> m_periodic;
    BHNode m_root;
    std::shared_ptr<BHNode> m_nodes;
    int m_node_size;
public:
    void initialize(std::shared_ptr<SPHParameters> param);
    void resize(const int particle_num, const int tree_size = 5);
    void make(std::vector<SPHParticle> & particles, const int particle_num);
    void set_kernel();
    int neighbor_search(const SPHParticle & p_i, std::vector<int> & neighbor_list, const bool is_ij = false);
};

}