/**
 * File   : decimation.cpp
 * Author : Jiayuan Mao
 * Email  : mjy14@mails.tsinghua.edu.cn
 * Date   : 2016-01-03 10:59:10
 * This file is part of the school project Decimation of course
 * ``Advanced Computational Geometry''.
 **/

#include "../include/decimation.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>

using std::cin;
using std::cout;
using std::cerr;
using std::endl;

namespace decimation {

static const bool DEBUG = false;

int Decimator::_contract(const VertexPair &p, queue_t &pq, int *const timestamp, int &current_timestamp) {
    int count = 0;
    Vertex *a = p.a, *b = p.b;
    Vertex::adjtriangles_t &tria = a->adjtriangles;
    Vertex::adjtriangles_t &trib = b->adjtriangles;
    Vertex::adjvertexes_t  &vera = a->adjvertexes;
    Vertex::adjvertexes_t  &verb = b->adjvertexes;

    Vector3f posa = a->to_vector();
    Vector3f posb = b->to_vector();
    a->update(p.final), b->update(p.final);

    if (DEBUG) cout << "Successfully updated position, begin validating" << endl;

    bool has_flip = false;

    if (DEBUG) cout << "iter tria: " << tria.size() << endl;
    for (auto t: tria) {
        if (DEBUG) cout << "   itering: " << t->valid << " " << t->a->id << " " << t->b->id << " " << t->c->id << endl;
        if (t->valid) {
            if (!t->has_vertex(b) && t->has_flipped())
                has_flip = true;
        }
    }

    if (DEBUG) cout << "iter trib: " << trib.size() << endl;
    for (auto t: trib) {
        if (DEBUG) cout << "   itering: " << t->valid << " " << t->a->id << " " << t->b->id << " " << t->c->id << endl;
        if (t->valid) {
            if (!t->has_vertex(a) && t->has_flipped())
                has_flip = true;
        }
    }

    if (has_flip) {
        if (DEBUG) cout << "Validation failed, rollbacking" << endl;
         a->update(posa), b->update(posb);
        return 0;
    }

    a->quad += b->quad;
    timestamp[a->id] = ++current_timestamp;
    timestamp[b->id] = -1;

    if (DEBUG) cout << "Validation passed, begin updating adjacent vertex list" << endl;

    if (DEBUG) cout << "iter vertex-a: " << vera.size() << endl;
    if (DEBUG) for (Vertex *v : vera) {
        cout << "    itering: " << v->id << endl;
    }
    auto iter = std::find(vera.begin(), vera.end(), b);
    if (iter != vera.end()) vera.erase(iter);
    if (DEBUG) cout << "iter vertex-b: " << verb.size() << endl;
    for (Vertex *v : verb) {
        if (DEBUG) cout << "   itering: " << v->id << endl;
        if (v != a) {
            vera.push_back(v);
            iter = std::find(v->adjvertexes.begin(), v->adjvertexes.end(), a);
            if (iter != v->adjvertexes.end()) v->adjvertexes.erase(iter);
            iter = std::find(v->adjvertexes.begin(), v->adjvertexes.end(), b);
            if (iter != v->adjvertexes.end()) v->adjvertexes.erase(iter);
            v->adjvertexes.push_back(a);
        }
    }

    for (Vertex *v : vera) {
        assert(v != a);
        pq.push(VertexPair(a, v, current_timestamp));
    }

    if (DEBUG) cout << "Successfully updated adjacent vertex list, begin updating adjacent triangle list" << endl;

    Vertex::adjtriangles_t new_tria;
    if (DEBUG) cout << "iter triangle-a: " << tria.size() << endl;
    for (auto t: tria) {
        if (DEBUG) cout << "   itering: " << t->valid << " " << t->a->id << " " << t->b->id << " " << t->c->id << endl;
        if (t->valid) {
            if (!t->has_vertex(b)) t->update_character(), new_tria.push_back(t);
            else t->valid = false, ++count;
        }
    }
    if (DEBUG) cout << "iter triangle-b: " << trib.size() << endl;
    for (auto t: trib) {
        if (DEBUG) cout << "   itering: " << t->valid << " " << t->a->id << " " << t->b->id << " " << t->c->id << endl;
        if (t->valid) {
            if (!t->has_vertex(a)) {
                t->update_vertex(b, a);
                new_tria.push_back(t);
            }
        }
    }
    tria.clear(), trib.clear();
    a->adjtriangles = new_tria;

    return count;
}

static inline bool _is_valid(const Decimator::VertexPair &p, int *const timestamp) {
    if (timestamp[p.a->id] < 0 || timestamp[p.a->id] > p.timestamp)
        return false;
    if (timestamp[p.b->id] < 0 || timestamp[p.b->id] > p.timestamp)
        return false;
    return true;
}

/**
 * Decimation main loop.
 * The decimation algorithm is based on the following paper:
 * 
 * Garland, M., & Heckbert, P. S. (1997, August). 
 * ``Surface simplification using quadric error metrics. ''
 * In Proceedings of the 24th annual conference on Computer graphics and interactive techniques (pp. 209-216). 
 * ACM Press/Addison-Wesley Publishing Co..
 **/
TriangleMesh *Decimator::decimate(TriangleMesh *&origin) {
    queue_t pq;
    kdtree_t *kdtree = NULL;

    int n = origin->vertexes.size();
    int m = origin->triangles.size();

    bool use_threshold = _threshold > eps;
    if (use_threshold) {
        cout << "Will use threshold " << _threshold << " to find make pair" << endl;
        kdtree = new kdtree_t(origin->vertexes);
        cout << "KDTree building finished" << endl;
    }

    for (Vertex *u : origin->vertexes) {
        if (use_threshold) {
            kdtree_t::list_t neighbours = kdtree->find_r(*u, _threshold);
            for (Vertex *v : neighbours) {
                if (u != v) pq.push(VertexPair(u, v));
            }
        }
        for (Vertex *v : u->adjvertexes) {
            if (u->id < v->id)
                pq.push(VertexPair(u, v));
        }
        if (u->id % 100 == 0) cout << "\rInitialized pair with " << std::setw(9) << u->id
                              << " (" << std::setprecision(5) << double(u->id)/n*100. << "%)  ";
    }

    cout << endl;
    cout << "Initialization finished, heap size " << pq.size() << endl;

    int *timestamp = new int[n];
    int current_timestamp = 0;
    memset(timestamp, 0, sizeof(int) * n);

    int final_m = _use_ratio ? m * _ratio : _target;
    int target_dec = m - final_m;
    int decimated = 0;
    while (!pq.empty() && decimated < target_dec) {
        // cout << "Current m: " << m << endl;
        VertexPair pair = pq.top(); pq.pop();
        if (_is_valid(pair, timestamp)) {
            if (DEBUG) cout << "Contract pair: " << pair.a->id << " " << pair.b->id << endl;
            if (DEBUG) cout << "Error = " << pair.error << endl;
            decimated += _contract(pair, pq, timestamp, current_timestamp);
        }
        if (decimated % 100 == 0) cout << "\rContracted pair "
                                  << std::setw(7) << pair.a->id << " "
                                  << std::setw(7) << pair.b->id
                                  << " (" << std::setw(6) << std::setprecision(5)
                                  << double(decimated)/target_dec*100. << "%)  ";
    }
    cout << endl;
    cout << "Contraction finished heap size=" << pq.size() << " current m=" << m << endl;

    n = 0, m = 0;
    TriangleMesh *res = new TriangleMesh();
    for (Vertex *u : origin->vertexes) {
        if (timestamp[u->id] == -1)
            u->id = -1;
        else {
            u->id = n++;
            res->append_vertex(new Vertex(u->x, u->y, u->z, u->id));
        }
    }
    for (Triangle *t : origin->triangles) {
        if (t->valid) {
            assert(t->a->id != -1 && t->b->id != -1 && t->c->id != -1);
            res->append_triangle(new Triangle(
                    res->vertexes[t->a->id],
                    res->vertexes[t->b->id],
                    res->vertexes[t->c->id]
            ));
        }
    }
    cout << "Export finished num_vertexes=" << res->vertexes.size() << " num_triangles=" << res->triangles.size() << endl;

    delete []timestamp;
    delete origin;

    origin = res;
    return origin;
}

} // End namespace decimation
