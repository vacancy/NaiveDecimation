/**
 * File   : model.h
 * Author : Jiayuan Mao
 * Email  : mjy14@mails.tsinghua.edu.cn
 * Date   : 2016-01-02 18:23:12
 * This file is part of the school project Decimation of course
 * ``Advanced Computational Geometry''.
 **/

#ifndef _INCLUDE_MODEL_H_
#define _INCLUDE_MODEL_H_

#include "linalg.h"
#include <vector>
#include <list>
#include <iostream>

namespace decimation {

class Triangle;

class Vertex : public Vector3f {
public:
    typedef std::vector<Vertex *> adjvertexes_t;
    typedef std::vector<Triangle *> adjtriangles_t;

    inline Vertex(void) : id(-1) { }
    inline Vertex(double x, double y, double z, int id) : Vector3f(x, y, z), id(id) { }
    inline Vertex(const Vector3f &vec, int id=-1) : Vector3f(vec.x, vec.y, vec.z), id(id) { }

    inline Vector3f to_vector(void) {
        return Vector3f(x, y, z);
    }
    inline void update(const Vector3f &vec) {
        x = vec.x, y = vec.y, z = vec.z;
    }

    int id;
    adjvertexes_t  adjvertexes;
    adjtriangles_t adjtriangles;
    Matrix4f quad;
};

class Triangle {
public:
    union {
        struct { Vertex *a, *b, *c; };
        struct { Vertex *_vertexes[3]; };
    };
    Vector4f character;
    bool valid;

    inline Triangle(void) : valid(true) { }
    inline Triangle(Vertex *a, Vertex *b, Vertex *c) : a(a), b(b), c(c), valid(true) { initialize(); }
    inline bool has_vertex(Vertex *target) const {
        return (a == target) || (b == target) || (c == target);
    }
    inline void update_vertex(Vertex *from, Vertex *to) {
        if (a == from) a = to;
        if (b == from) b = to;
        if (c == from) c = to;
        initialize();
    }
    inline bool has_flipped(void) const {
        Vector3f new_norm = cross((*b - *a), (*c - *a)).normalize();
        return dot(character.to3(), new_norm) < 0;
    }

    inline void update_character(void) {
        Vector3f norm = cross((*b - *a), (*c - *a)).normalize();
        double d = - dot(*a, norm);
        character = Vector4f(norm.x, norm.y, norm.z, d);
    }
    inline void initialize(void) {
        update_character();
    }
};

class TriangleMesh {
public:
    typedef std::vector<Vertex *>   vertex_list_t;
    typedef std::vector<Triangle *> triangle_list_t;

    inline TriangleMesh(void) { }
    virtual inline ~TriangleMesh(void) {
        for (Vertex *v : vertexes) delete v;
        for (Triangle *t : triangles) delete t;
        vertexes.clear();
        triangles.clear();
    }

    inline void append_vertex(Vertex *ver) {
        vertexes.push_back(ver);
    }

    inline void append_triangle(Triangle *tri) {
        triangles.push_back(tri);
    }

    static TriangleMesh *from_stream(std::istream &);
    void to_stream(std::ostream &);

    vertex_list_t   vertexes;
    triangle_list_t triangles;

private:
};

} // End namespace decimation

#endif
