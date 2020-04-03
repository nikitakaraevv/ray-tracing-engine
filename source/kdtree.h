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

class kdtree
{
private:
    struct node
    {
        node(const Particle& pt) : point_(pt), left_(nullptr), right_(nullptr)
        {
        }
        float get(size_t index) const
        {
            return point_.position()[index];
        }
        double distance(const Particle& pt) const
        {
            return dist(point_.position(), pt.position());
        }
        Particle point_;
        node* left_;
        node* right_;
    };
    
    
    node* root_;
    node* best_;
    double best_dist_;
    size_t visited_;
    std::vector<node> nodes_;
    std::vector<node> k_nearest_;
 
    struct node_cmp
    {
        node_cmp(size_t index) : index_(index)
        {
        }
        bool operator()(const node& n1, const node& n2) const
        {
            return n1.point_.position()[index_] < n2.point_.position()[index_];
        }
        size_t index_;
    };
    
    // min heap comparator
    struct dist_cmp_max
    {
        dist_cmp_max(Particle target) : target_(target)
        {
        }
       bool operator()(const node& n1, const node& n2)
       {
           return n1.distance(target_) < n2.distance(target_);
       }
       Particle target_;
    };
    
    struct dist_cmp_min
    {
        dist_cmp_min(Particle target) : target_(target)
        {
        }
       bool operator()(const node& n1, const node& n2)
       {
           return n1.distance(target_) > n2.distance(target_);
       }
       Particle target_;
    };

 
    node* make_tree(size_t begin, size_t end, size_t index)
    {
        if (end <= begin)
            return nullptr;
        size_t n = begin + (end - begin)/2;
        std::nth_element(&nodes_[begin], &nodes_[n], &nodes_[end], node_cmp(index));
        index = (index + 1) % 3;
        nodes_[n].left_ = make_tree(begin, n, index);
        nodes_[n].right_ = make_tree(n + 1, end, index);
        return &nodes_[n];
    }
 
    void nearest(node* root, const Particle& point, size_t index)
    {
        if (root == nullptr)
            return;
        ++visited_;
        double d = root->distance(point);
        if (best_ == nullptr || d < best_dist_)
        {
            best_dist_ = d;
            best_ = root;
        }
        if (best_dist_ == 0)
            return;
        double dx = root->get(index) - point.position()[index];
        index = (index + 1) % 3;
        nearest(dx > 0 ? root->left_ : root->right_, point, index);
        if (dx * dx >= best_dist_)
            return;
        nearest(dx > 0 ? root->right_ : root->left_, point, index);
    }
    
    void knearest(node* root, const Particle& point, size_t index)
    {
        if (root == nullptr)
            return;
        
        ++visited_;
        double d = root->distance(point);
        if (d < best_dist_)
        {
            pop_heap (k_nearest_.begin(),k_nearest_.end(), dist_cmp_max(point));
            node n = k_nearest_.front();
            k_nearest_.pop_back();
           
            best_dist_ = n.distance(point);
            k_nearest_.push_back(*root);
            push_heap (k_nearest_.begin(),k_nearest_.end(),  dist_cmp_max(point));
        }
        if (best_dist_ == 0)
            return;
        double dx = root->get(index) - point.position()[index];
        index = (index + 1) % 3;
        knearest(dx > 0 ? root->left_ : root->right_, point, index);
        if (dx * dx >= best_dist_)
            return;
        knearest(dx > 0 ? root->right_ : root->left_, point, index);
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
    kdtree(iterator begin, iterator end)
    {
        best_ = nullptr;
        best_dist_ = 0;
        visited_ = 0;
        nodes_.reserve(std::distance(begin, end));
        for (auto i = begin; i != end; ++i)
            nodes_.emplace_back(*i);
        root_ = make_tree(0, nodes_.size(), 0);
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
    kdtree(func&& f, size_t n)
    {
        best_ = nullptr;
        best_dist_ = 0;
        visited_ = 0;
        nodes_.reserve(n);
        for (size_t i = 0; i < n; ++i)
            nodes_.emplace_back(f());
        root_ = make_tree(0, nodes_.size(), 0);
    }
 
    /**
     * Returns true if the tree is empty, false otherwise.
     */
    bool empty() const
    {
        return nodes_.empty();
    }
 
    /**
     * Returns the number of nodes visited by the last call
     * to nearest().
     */
    size_t visited() const
    {
        return visited_;
    }
 
    /**
     * Returns the distance between the input point and return value
     * from the last call to nearest().
     */
    double distance() const
    {
        return best_dist_;
    }
 
    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @param the nearest point in the tree to the given point
     */
    const Particle& nearest(const Particle& pt)
    {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_ = 0;
        nearest(root_, pt, 0);
        return best_->point_;
    }
    
     void knearest(const Particle& pt, int k, vector<Particle>& result)
    {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        if (k > nodes_.size())
            throw std::logic_error("k is greater than the number of nodes");
        
        // initialize the vector
        for (int i =0; i< k; i++)
            k_nearest_.push_back(nodes_[i]);
            
        make_heap(k_nearest_.begin(), k_nearest_.end(), dist_cmp_max(pt));
        visited_ = 0;
        best_dist_ = k_nearest_[0].distance(pt);
        knearest(root_, pt, 0);
        
        sort_heap (k_nearest_.begin(),k_nearest_.end(), dist_cmp_max(pt));
        for (int i =0; i< k; i++)
            result.push_back(k_nearest_[i].point_);
    }
    
    
};
 
