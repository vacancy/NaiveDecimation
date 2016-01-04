/**
 * File   : model.cpp
 * Author : Jiayuan Mao
 * Email  : mjy14@mails.tsinghua.edu.cn
 * Date   : 2015-11-15 15:51:20
 * This file is part of the school project Decimation of course
 * ``Advanced Computational Geometry''.
 **/

#include "../include/model.h"
#include <iostream>
#include <cstring>

using std::cout;
using std::endl;

namespace decimation {

TriangleMesh *TriangleMesh::from_stream(std::istream &f) {
    TriangleMesh *mesh = new TriangleMesh();

    int nvertexes = 0, ntriangles = 0;
    char buffer[128];
    while (f.getline(buffer, 128)) {
        if (buffer[0] == '#') {
            continue;
        } else if (buffer[0] == 'v') {
            char op;
            double x, y, z;
            sscanf(buffer, "%c %lf %lf %lf", &op, &x, &y, &z);
            Vertex *vertex = new Vertex(x, y, z, nvertexes++);
            mesh->append_vertex(vertex);
        } else if (buffer[0] == 'f') {
            char op;
            int x, y, z;
            sscanf(buffer, "%c %d %d %d", &op, &x, &y, &z);

            ++ntriangles;
            Triangle *triangle = new Triangle(mesh->vertexes[x-1],  mesh->vertexes[y-1], mesh->vertexes[z-1]);
            Matrix4f quad = MatrixFactory::get_quad_matrix(triangle->character);
            mesh->append_triangle(triangle);
#pragma unroll
            for (int i = 0; i < 3; ++i) {
                triangle->_vertexes[i]->adjtriangles.push_back(triangle);
                triangle->_vertexes[i]->adjvertexes.push_back(triangle->_vertexes[(i+1) % 3]);
                triangle->_vertexes[i]->adjvertexes.push_back(triangle->_vertexes[(i+2) % 3]);
                triangle->_vertexes[i]->quad += quad;
            }
        }
    }

    cout << "Read-in finished num_vertexes=" << nvertexes << " num_triangles=" << ntriangles << endl;

    return mesh;
}

void TriangleMesh::to_stream(std::ostream &f) {
    char buffer[128];

    sprintf(buffer, "# %d %d\n", vertexes.size(), triangles.size());
    f.write(buffer, strlen(buffer));
    for (Vertex *v : vertexes) {
        sprintf(buffer, "v %f %f %f\n", v->x, v->y, v->z);
        f.write(buffer, strlen(buffer));
    }
    for (Triangle *t : triangles) {
        sprintf(buffer, "f %d %d %d\n", t->a->id+1, t->b->id+1, t->c->id+1);
        f.write(buffer, strlen(buffer));
    }

    cout << "Output finished num_vertexes=" << vertexes.size() << " num_triangles=" << triangles.size() << endl;
}

}; // end namespace decimation
