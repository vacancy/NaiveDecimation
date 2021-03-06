/**
 * File   : vector.h
 * Author : Jiayuan Mao
 * Email  : mjy14@mails.tsinghua.edu.cn
 * Date   : 2016-01-04 10:34:59
 * This file is part of the school project RayTracing of course
 * ``Advanced Computational Geometry''.
 * This file is part of the school project Decimation of course
 * ``Advanced Computational Geometry''.
 **/

#ifndef _INCLUDE_VECTOR_KD_H_
#define _INCLUDE_VECTOR_KD_H_

#include "linalg.h"
#include <vector>

namespace decimation {

// {{{ Begin file: include/kd.h

/**
 * Using CRTP mode to inherit.
 **/
template <typename NodeType, typename ValueType>
struct KDNodeBase {
    typedef NodeType node_t;
    typedef ValueType value_t;

    KDNodeBase(void) : axis(-1), split(0), lson(NULL), rson(NULL) { }
    inline bool is_leaf(void) const { return lson == NULL && rson == NULL; }

    int axis;
    value_t split;
    node_t  *lson, *rson;
};

template <typename PtrType>
struct KDCompareBase {
    typedef PtrType ptr_t;

    virtual bool operator ()(const ptr_t &lhs, const ptr_t &rhs) const = 0;
};

// }}} End file: include/kd.h

template <typename T>
class VecKDTree {
public:
    typedef typename T::value_t value_t;
    typedef T vector_t;
    typedef T *ptr_vector_t;
    typedef std::vector<ptr_vector_t> list_t;
    typedef typename list_t::iterator list_iter_t;

    struct VecKDNode : KDNodeBase<VecKDNode, value_t> {
        VecKDNode(void) : KDNodeBase<VecKDNode, value_t>() { }
        list_t vectors;
    };

    struct VecKDCompare : KDCompareBase<ptr_vector_t> {
        VecKDCompare(int axis = 0) : KDCompareBase<ptr_vector_t>(), axis(axis) { }
        virtual inline bool operator ()(const typename KDCompareBase<ptr_vector_t>::ptr_t &lhs,
                                        const typename KDCompareBase<ptr_vector_t>::ptr_t &rhs) const {
            return lhs->_values[axis] < rhs->_values[axis];
        }

        int axis;
    };

    VecKDTree(const list_t &wrapper, int max_leaf_size = 30)
            : wrapper(wrapper), root(NULL), _max_leaf_size(max_leaf_size) {

        initialize();
    }

    inline void initialize() { _build(root, wrapper.begin(), wrapper.end(), 0); }

    inline list_t find_r(const vector_t &center, const value_t r) {
        list_t res;
        _traverse_r(res, root, center, r*r);
        return res;
    }

    inline list_t find_r_bf(const vector_t &center, const value_t r) {
        list_t res;
        value_t r2 = r*r;
        for (auto v : wrapper) {
            if ((*v - center).l2() < r2)
                res.push_back(v);
        }
        return res;
    }

    list_t wrapper;
    VecKDNode *root;

protected:
    void _build(VecKDNode *&root, list_iter_t begin, list_iter_t end, int current) {
        if (begin == end) {
            return ;
        }

        root = new VecKDNode();

        size_t n = end - begin;
        if (n <= _max_leaf_size) {
            root->vectors = list_t(begin, end);
        } else {
            std::nth_element(begin, begin + n/2, end, VecKDCompare(current));
            root->axis = current;
            root->split = (*(begin + n/2))->_values[current];

            int next = (current + 1) % 3;
            _build(root->lson, begin, begin + n/2, next);
            _build(root->rson, begin + n/2, end, next);
        }
    }

    void _traverse_r(list_t &res, const VecKDNode *root, const vector_t &center, const value_t r2) {
        if (root->is_leaf()) {
            for (auto v : root->vectors) {
                if ((*v - center).l2() < r2)
                    res.push_back(v);
            }
        } else {
            bool in_left = center[root->axis] < root->split;
            bool ignore = (center[root->axis] - root->split) * (center[root->axis] - root->split) > r2;
            if (in_left) {
                _traverse_r(res, root->lson, center, r2);
                if (!ignore) _traverse_r(res, root->rson, center, r2);
            } else {
                _traverse_r(res, root->rson, center, r2);
                if (!ignore) _traverse_r(res, root->lson, center, r2);
            }
        }
    }

private:
    int _max_leaf_size = 30;
};

} // end namespace decimation


#endif
