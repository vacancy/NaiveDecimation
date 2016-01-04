/**
 * File   : decimation.h
 * Author : Jiayuan Mao
 * Email  : mjy14@mails.tsinghua.edu.cn
 * Date   : 2016-01-02 18:23:12
 * This file is part of the school project Decimation of course
 * ``Advanced Computational Geometry''.
 **/

#ifndef _INCLUDE_DECIMATION_H_
#define _INCLUDE_DECIMATION_H_

#include "linalg.h"
#include "model.h"
#include "veckd.h"
#include <queue>

namespace decimation {

class Decimator {
public:
    struct VertexPair {
        union {
            struct { Vertex *a, *b; };
            struct { Vertex *first, *second; };
        };

        Vector3f final;
        double   error;

        int timestamp;

        VertexPair(void) : timestamp(0) { }
        VertexPair(Vertex *a, Vertex *b, int t = 0) : a(a), b(b), timestamp(t) { initialize(); }

        inline void initialize(void) {
            Matrix4f q = a->quad + b->quad;
            q[3][0] = q[3][1] = q[3][2] = 0; q[3][3] = 1;
            Matrix4f inv_q;
            bool invertible = inverse(q, inv_q);
            if (invertible) {
                final = (inv_q * affine<double>(0, 0, 0)).to3();
            } else {
                final = ((*a) + (*b)) * 0.5;
            }
            // std::cout << "invertible=" << invertible << " final=" << final;
            // std::cout << " compared with=" << ((*a) + (*b)) * 0.5 << std::endl;
            error = quadric(a->quad + b->quad, affine(final.x, final.y, final.z));
        }

        inline bool operator < (const VertexPair &rhs) const {
            return error > rhs.error;
        }
    };

    typedef std::priority_queue<VertexPair> queue_t;
    typedef VecKDTree<Vertex> kdtree_t;

    explicit Decimator(double ratio, double threshold) : _ratio(ratio), _threshold(threshold), _use_ratio(true) {
        std::cout << "Creating decimator using ratio=" << ratio << " threshold=" << threshold << std::endl;
    }

    explicit Decimator(int target, double threshold) : _target(target), _threshold(threshold), _use_ratio(false) {
        std::cout << "Creating decimator using target number=" << target << " threshold=" << threshold << std::endl;
    }

    TriangleMesh *decimate(TriangleMesh *&origin);

private:
    int _contract(const VertexPair &p, queue_t &pq, int *const timestamp, int &current_timestamp);

    union {
        double _ratio;
        int    _target;
    };
    bool   _use_ratio;
    double _threshold;
};

} // End namespace decimation

#endif
