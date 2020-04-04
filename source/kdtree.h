#pragma once
#include <algorithm>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
 

/**
 * C++ k-d tree implementation, based on the C version at rosettacode.org.
 */

class kdtree {
private:
    struct node
    {
        node(const Particle& pt) : m_point(pt), m_left(nullptr), m_right(nullptr) {
        }
        float get(size_t index) const {
            return m_point.position()[index];
        }
        double distance(const Particle& pt) const {
            return dist(m_point.position(), pt.position());
        }
        Particle m_point;
        node* m_left;
        node* m_right;
    };
    
    
    node* m_root;
    node* m_best;
    double m_bestdist;
    size_t m_visited;
    std::vector<node> m_nodes;
    std::vector<node> m_knearest;
 
    struct node_cmp
    {
        node_cmp(size_t index) : m_index(index) {}
        
        bool operator()(const node& n1, const node& n2) const {
            return n1.m_point.position()[m_index] < n2.m_point.position()[m_index];
        }
        size_t m_index;
    };
    
    // min heap comparator
    struct dist_cmp_max
    {
        dist_cmp_max(Particle m_target) : m_target(m_target){}
        
       bool operator()(const node& n1, const node& n2) {
           return n1.distance(m_target) < n2.distance(m_target);
       }
       Particle m_target;
    };
    
    struct dist_cmp_min
    {
        dist_cmp_min(Particle m_target) : m_target(m_target) { }
       bool operator()(const node& n1, const node& n2) {
           return n1.distance(m_target) > n2.distance(m_target);
       }
       Particle m_target;
    };

 
    node* make_tree(size_t begin, size_t end, size_t index) {
        if (end <= begin)
            return nullptr;
        size_t n = begin + (end - begin)/2;
        std::nth_element(&m_nodes[begin], &m_nodes[n], &m_nodes[end], node_cmp(index));
        index = (index + 1) % 3;
        m_nodes[n].m_left = make_tree(begin, n, index);
        m_nodes[n].m_right = make_tree(n + 1, end, index);
        return &m_nodes[n];
    }
 
    void nearest(node* root, const Particle& point, size_t index) {
        if (root == nullptr)
            return;
        ++m_visited;
        double d = root->distance(point);
        if (m_best == nullptr || d < m_bestdist) {
            m_bestdist = d;
            m_best = root;
        }
        if (m_bestdist == 0)
            return;
        double dx = root->get(index) - point.position()[index];
        index = (index + 1) % 3;
        nearest(dx > 0 ? root->m_left : root->m_right, point, index);
        if (dx * dx >= m_bestdist)
            return;
        nearest(dx > 0 ? root->m_right : root->m_left, point, index);
    }
    
    void knearest(node* root, const Particle& point, size_t index) {
        if (root == nullptr)
            return;
        
        ++m_visited;
        double d = root->distance(point);
        if (d < m_bestdist) {
            pop_heap (m_knearest.begin(),m_knearest.end(), dist_cmp_max(point));
            node n = m_knearest.front();
            m_knearest.pop_back();
           
            m_bestdist = n.distance(point);
            m_knearest.push_back(*root);
            push_heap (m_knearest.begin(),m_knearest.end(),  dist_cmp_max(point));
        }
        if (m_bestdist == 0)
            return;
        double dx = root->get(index) - point.position()[index];
        index = (index + 1) % 3;
        knearest(dx > 0 ? root->m_left : root->m_right, point, index);
        if (dx * dx >= m_bestdist)
            return;
        knearest(dx > 0 ? root->m_right : root->m_left, point, index);
    }
    
public:
    kdtree(const kdtree&) = delete;
    kdtree& operator=(const kdtree&) = delete;
    /**
     * Constructor taking a pair of iterators. Adds each
     * point in the range [begin, end) to the tree.
     *
     * @param begin start of range
     * @param end end of range
     */
    template<typename iterator>
    kdtree(iterator begin, iterator end) {
        m_best = nullptr;
        m_bestdist = 0;
        m_visited = 0;
        m_nodes.reserve(std::distance(begin, end));
        for (auto i = begin; i != end; ++i)
            m_nodes.emplace_back(*i);
        m_root = make_tree(0, m_nodes.size(), 0);
    }
 
    /**
     * Constructor taking a function object that generates
     * points. The function object will be called n times
     * to populate the tree.
     *
     * @param f function that returns a point
     * @param n number of points to add
     */
    template<typename func>
    kdtree(func&& f, size_t n) {
        m_best = nullptr;
        m_bestdist = 0;
        m_visited = 0;
        m_nodes.reserve(n);
        for (size_t i = 0; i < n; ++i)
            m_nodes.emplace_back(f());
        m_root = make_tree(0, m_nodes.size(), 0);
    }
 
    /**
     * Returns true if the tree is empty, false otherwise.
     */
    bool empty() const { return m_nodes.empty(); }
 
    /**
     * Returns the number of nodes visited by the last call
     * to nearest().
     */
    size_t visited() const { return m_visited; }
 
    /**
     * Returns the distance between the input point and return value
     * from the last call to nearest().
     */
    double distance() const { return m_bestdist; }
 
    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @param the nearest point in the tree to the given point
     */
    const Particle& nearest(const Particle& pt) {
        if (m_root == nullptr)
            throw std::logic_error("tree is empty");
        m_best = nullptr;
        m_visited = 0;
        m_bestdist = 0;
        nearest(m_root, pt, 0);
        return m_best->m_point;
    }
    
     void knearest(const Particle& pt, int k, vector<Particle>& result) {
        if (m_root == nullptr)
            throw std::logic_error("tree is empty");
        if (k > m_nodes.size())
            throw std::logic_error("k is greater than the number of nodes");
        m_knearest.clear();
        // initialize the vector
        for (int i =0; i< k; i++)
            m_knearest.push_back(m_nodes[i]);
            
        make_heap(m_knearest.begin(), m_knearest.end(), dist_cmp_max(pt));
        m_visited = 0;
        m_bestdist = m_knearest[0].distance(pt);
        knearest(m_root, pt, 0);
        
        sort_heap (m_knearest.begin(),m_knearest.end(), dist_cmp_max(pt));
        for (int i =0; i< k; i++)
            result.push_back(m_knearest[i].m_point);
    }
    
    
};
 
