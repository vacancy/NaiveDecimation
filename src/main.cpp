/**
 * File   : main.cpp
 * Author : Jiayuan Mao
 * Email  : mjy14@mails.tsinghua.edu.cn
 * Date   : 2016-01-03 10:18:31
 * This file is part of the school project Decimation of course
 * ``Advanced Computational Geometry''.
 **/

#include "../include/linalg.h"
#include "../include/model.h"
#include "../include/veckd.h"
#include "../include/decimation.h"
#include <iostream>
#include <fstream>

using namespace decimation;
using namespace std;


void test_matrix(void) {
    Matrix4f ori, inv;
    ori[0][0] = 2, ori[0][1] = 3, ori[0][2] = 5, ori[0][3] = 5;
    ori[1][0] = 1, ori[1][1] = 1, ori[1][2] = 5, ori[1][3] = 2;
    ori[2][0] = 2, ori[2][1] = 1, ori[2][2] = 3, ori[2][3] = 2;
    ori[3][0] = 1, ori[3][1] = 1, ori[3][2] = 3, ori[3][3] = 4;

    cout << ori;
    cout << inverse(ori, inv) << endl;
    cout << inv;
    cout << ori*inv << endl;
    getchar();
}

void test_kdtree(TriangleMesh *mesh) {
    typedef VecKDTree<Vertex>::list_t list_vec3f;
    VecKDTree<Vertex> *kdtree = new VecKDTree<Vertex>(mesh->vertexes);
    list_vec3f res1, res2;
    res1 = kdtree->find_r(Vertex(Vector3f(0, 0, 0)), 0.1);
    res2 = kdtree->find_r_bf(Vertex(Vector3f(0, 0, 0)), 0.1);
    cout << res1.size() << " " << res2.size() << endl;

    res1 = kdtree->find_r(Vertex(Vector3f(0, 0.5, 0.25)), 0.5);
    res2 = kdtree->find_r_bf(Vertex(Vector3f(0, 0.5, 0.25)), 0.5);
    cout << res1.size() << " " << res2.size() << endl;
    getchar();
}

int main(int argc, char *argv[]) {
    int input_pos = 0, output_pos = 0;
    int ratio_pos = 0, target_pos = 0, threshold_pos = 0;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--input") == 0) input_pos = ++i;
        if (strcmp(argv[i], "--output") == 0) output_pos = ++i;
        if (strcmp(argv[i], "--ratio") == 0) ratio_pos = ++i;
        if (strcmp(argv[i], "--target") == 0) target_pos = ++i;
        if (strcmp(argv[i], "--threshold") == 0) threshold_pos = ++i;
    }

    if (!(input_pos && output_pos && (ratio_pos ^ target_pos))) {
        cout << "Invalid parameter" << endl;
        return 1;
    }

    ifstream fin(argv[input_pos], ios::in);
    TriangleMesh *mesh = TriangleMesh::from_stream(fin);
    fin.close();

    double threshold = 0, ratio = 0; int target = 0;
    if (threshold_pos) sscanf(argv[threshold_pos], "%lf", &threshold);
    Decimator *decimator = NULL;
    if (ratio_pos) {
        sscanf(argv[ratio_pos], "%lf", &ratio);
        decimator = new Decimator(ratio, threshold);
    } else {
        sscanf(argv[target_pos], "%d", &target);
        decimator = new Decimator(target, threshold);
    }

    decimator->decimate(mesh);

    ofstream fout(argv[output_pos], ios::out);
    mesh->to_stream(fout);
    fout.close();

    return 0;
}
