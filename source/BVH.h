#pragma once

#include <algorithm>
#include <vector>
#include <functional>
#include <assert.h>
#include <unordered_map>
#include "Vec3.h"
#include "AABB.h"

using namespace std;


class BVH {
    public:
        AABB m_aabb;
        BVH *m_child_1;
        BVH *m_child_2;
        vector<int> m_triangles;
        int m_cut_axis;
        float m_cut;

        BVH() {
            this->m_aabb = AABB();
            this->m_child_1 = NULL;
            this->m_child_2 = NULL;
            this->m_triangles = vector<int>();
            this->m_cut_axis = 0;
            this->m_cut = 0;
        }

        // Deep copy
        BVH(const BVH &to_copy) {
            if (to_copy.m_child_1 != NULL)
                this->m_child_1 = new BVH(*to_copy.m_child_1);
            else
                this->m_child_1 = NULL;
            BVH *new_child_2;
            if (to_copy.m_child_2 != NULL)
                this->m_child_2 = new BVH(*to_copy.m_child_2);
            else
                this->m_child_2 = NULL;
            this->m_aabb = to_copy.m_aabb;
            this->m_triangles = to_copy.m_triangles;
            this->m_cut_axis = to_copy.m_cut_axis;
            this->m_cut = to_copy.m_cut;
        }

        ~BVH() {
            if (this->m_child_1 != NULL) {
                delete this->m_child_1;
            }
            if (this->m_child_2 != NULL) {
                delete this->m_child_2;
            }
        }

        // Deep copy
        BVH &operator=(const BVH &to_copy) {
            if (this != &to_copy) { // protect against invalid self-assignment
                // 1: allocate new memory and copy the elements
                BVH *new_child_1;
                if (to_copy.m_child_1 != NULL)
                    new_child_1 = new BVH(*to_copy.m_child_1);
                else
                    new_child_1 = NULL;
                BVH *new_child_2;
                if (to_copy.m_child_2 != NULL)
                    new_child_2 = new BVH(*to_copy.m_child_2);
                else
                    new_child_2 = NULL;
                // 2: deallocate old memory
                if (this->m_child_1 != NULL)
                    delete this->m_child_1;
                if (this->m_child_2 != NULL)
                    delete this->m_child_2;
                // 3: assign the new memory to the object
                this->m_child_1 = new_child_1;
                this->m_child_2 = new_child_2;
                this->m_aabb = to_copy.m_aabb;
                this->m_triangles = to_copy.m_triangles;
                this->m_cut_axis = to_copy.m_cut_axis;
                this->m_cut = to_copy.m_cut;
            }
            return *this;
        }

        void check_cut_axis() const {
            assert(0 <= this->m_cut_axis && this->m_cut_axis <= 2);
            if (m_child_1 != NULL) {
                m_child_1->check_cut_axis();
            }
            if (m_child_2 != NULL) {
                m_child_2->check_cut_axis();
            }
        }

        void from_triangles(vector<int>::const_iterator build_t_begin, vector<int>::const_iterator build_t_end, const vector<Vec3f> &vertices, const vector<Vec3i> &triangles) {
            // Get m_triangles
            this->m_triangles.assign(build_t_begin, build_t_end);
            // Calculate m_aabb
            float eps = 1e-3;
            Vec3f minCorner(m_aabb.minCorner()) ,
                maxCorner(m_aabb.maxCorner());
            
            for (int t_idx: m_triangles) {
                for (int v_idx=0; v_idx<3; v_idx++) {
                    for (int i=0; i<3; i++) {
                        float coord = vertices[triangles[t_idx][v_idx]][i];
                        minCorner[i] = fmin(minCorner[i], coord - __FLT_EPSILON__ * 2);
                        maxCorner[i] = fmax(maxCorner[i], coord + __FLT_EPSILON__ * 2);
                    }
                }
            }
            this->m_aabb.minCorner() = minCorner;
            this->m_aabb.maxCorner() = maxCorner;
            
            // Calculate children
            if (this->m_triangles.size() <= 1) {
                this->m_cut_axis = 0;
                this->m_child_1 = NULL;
                this->m_child_2 = NULL;
                assert(0 <= this->m_cut_axis && this->m_cut_axis <= 2);
                return;
            }
            else {
                // Calculate cut axis
                float longlength = 0;
                assert(0 <= this->m_cut_axis && this->m_cut_axis <= 2);
                for (int i=0; i<3; i++) {
                    float lengthi = this->m_aabb.maxCorner()[i] - this->m_aabb.minCorner()[i];
                    if (lengthi > longlength) {
                        longlength = lengthi;
                        this->m_cut_axis = i;
                    }
                }
                std::function<bool(int, int)> comparator = [this, &triangles, &vertices](int ti1, int ti2) {
                        return vertices[triangles[ti1][0]][m_cut_axis]
                             + vertices[triangles[ti1][1]][m_cut_axis]
                             + vertices[triangles[ti1][2]][m_cut_axis]
                             < vertices[triangles[ti2][0]][m_cut_axis]
                             + vertices[triangles[ti2][1]][m_cut_axis]
                             + vertices[triangles[ti2][2]][m_cut_axis];
                };
                sort(
                    this->m_triangles.begin(),
                    this->m_triangles.end(),
                    comparator
                );
                this->m_child_1 = new BVH();
                this->m_child_1->from_triangles(
                    m_triangles.begin(),
                    m_triangles.begin() + m_triangles.size() / 2,
                    vertices,
                    triangles
                );
                this->m_child_2 = new BVH();
                this->m_child_2->from_triangles(
                    m_triangles.begin() + m_triangles.size() / 2,
                    m_triangles.end(),
                    vertices,
                    triangles
                );
                assert(0 <= this->m_cut_axis && this->m_cut_axis <= 2);
            }
        }

        vector<float> intersection(Ray &r, Vec3i &t, const vector<Vec3f> &vertices, const vector<Vec3i> &triangles) const {
            Vec3f entry;
            Vec3f exit;
            if (this->m_aabb.intersectionTest(r, entry, exit)) {
                //cout << "AABB intersection found" << endl;
                if (this->m_child_1 == NULL) {
                    vector<float> best_inter = vector<float>();
                    for (auto t_idx : this->m_triangles) {
                        Vec3f p0 = vertices[triangles[t_idx][0]];
                        Vec3f p1 = vertices[triangles[t_idx][1]];
                        Vec3f p2 = vertices[triangles[t_idx][2]];
                        float u, v, t0;
                        bool intersect = r.triangleIntersect(p0, p1, p2, u, v, t0);
                        vector<float> this_inter = {1.f - u - v, u, v, t0};
                        if (intersect && (best_inter.size() == 0 || best_inter[3] > this_inter[3])) {
                            best_inter = this_inter;
                            t = triangles[t_idx];
                        }
                    }
                    return best_inter;
                }
                else {
                    Vec3i t2;
                    vector<float> inter1 = this->m_child_1->intersection(r, t, vertices, triangles);
                    vector<float> inter2 = this->m_child_2->intersection(r, t2, vertices, triangles);
                    if (inter1.size() == 0 || (inter2.size() > 0 && inter1[3] > inter2[3])) {
                        inter1 = inter2;
                        t = t2;
                    }
                    return inter1;
                }
            }
            else {
                return vector<float>();
            }
        }
};
