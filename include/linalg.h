/**
 * File   : linalg.h
 * Author : Jiayuan Mao
 * Email  : mjy14@mails.tsinghua.edu.cn
 * Date   : 2016-01-02 11:18:17
 * This file is part of the school project Decimation of course
 * ``Advanced Computational Geometry''.
 **/

#ifndef _INCLUDE_LINALG_H_
#define _INCLUDE_LINALG_H_

#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

namespace decimation {

const double eps = 1e-7;

template<typename T>
inline int sgn(const T &a) {
    return a < -eps ? -1 : a > eps;
}

template<typename T, int N>
class Vector {
public:
    typedef T value_t;

    Vector(void) { memset(_values, 0, sizeof(_values)); }
    inline T &operator[](int index) { return _values[index]; }
    inline const T &operator[](int index) const { return _values[index]; }
    T _values[N];

private:

};

template <typename T>
class Vector<T, 3> {
public:
    typedef T value_t;

    Vector(void) { memset(_values, 0, sizeof(_values)); }
    Vector(const T &x, const T &y, const T &z) : x(x), y(y), z(z) { }
    inline T &operator[](int index) { return _values[index]; }
    inline const T &operator[](int index) const { return _values[index]; }
    inline Vector<T, 3> normalize(void) const {
        T norm = 1. / sqrt(x*x + y*y + z*z);
        return Vector<T, 3>(x*norm, y*norm, z*norm);
    };
    inline T l2(void) const {
        return x*x + y*y + z*z;
    };

    union {
        struct { T x, y, z; };
        struct { T _values[3]; };
    };
private:

};

template <typename T>
class Vector<T, 4> {
public:
    typedef T value_t;

    Vector(void) { memset(_values, 0, sizeof(_values)); }
    Vector(const T &x, const T &y, const T &z, const T &t = 1) : x(x), y(y), z(z), t(t) { }
    inline T &operator[](int index) { return _values[index]; }
    inline const T &operator[](int index) const { return _values[index]; }
    inline Vector<T, 3> to3(void) const {
        return Vector<T, 3>(x, y, z);
    }

    union {
        struct { T x, y, z, t; };
        struct { T _values[4]; };
    };
private:

};

template <typename T, int N>
inline Vector<T, N> operator +(const Vector<T, N> &a, const Vector<T, N> &b) {
    Vector<T, N> res;
#pragma unroll
    for (int i = 0; i < N; ++i)
        res._values[i] = a._values[i] + b._values[i];
    return res;
}

template <typename T, int N>
inline Vector<T, N> operator -(const Vector<T, N> &a, const Vector<T, N> &b) {
    Vector<T, N> res;
#pragma unroll
    for (int i = 0; i < N; ++i)
        res._values[i] = a._values[i] - b._values[i];
    return res;
}

template <typename T, int N>
inline Vector<T, N> operator *(const Vector<T, N> &a, const T &b) {
    Vector<T, N> res;
#pragma unroll
    for (int i = 0; i < N; ++i)
        res._values[i] = a._values[i] * b;
    return res;
}

typedef Vector<double, 3> Vector3f;
typedef Vector<double, 4> Vector4f;

template <typename T>
inline T dot(const Vector<T, 3> &a, const Vector<T, 3> &b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

template <typename T>
inline Vector<T, 3> cross(const Vector<T, 3> &op0, const Vector<T, 3> &op1) {
    return Vector<T, 3>(
        op0.y * op1.z - op0.z * op1.y,
        op0.z * op1.x - op0.x * op1.z,
        op0.x * op1.y - op0.y * op1.x
    );
}

template <typename T>
inline Vector<T, 4> affine(const T &x, const T &y, const T &z) {
    return Vector<T, 4>(x, y, z, 1);
}

template <typename T>
inline Vector<T, 4> affine(const Vector<T, 3> &v) {
    return Vector<T, 4>(v.x, v.y, v.z, 1);
}

template <typename T, int N, int M>
class Matrix {
public:
    typedef T value_t;
    typedef T row_t[M];

    Matrix(void) { memset(_values, 0, sizeof(_values)); }
    inline T &get(int x, int y) { return _values[x][y]; }
    inline void set(int x, int y, const T &v) { _values[x][y] = v; }
    inline row_t &operator[](int index) { return _values[index]; }
    inline const row_t &operator[](int index) const { return _values[index]; }

    row_t _values[N];
    int   rows = N, cols = M;

    inline Matrix<T, N, M> &operator +=(const Matrix<T, N, N> &rhs) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                _values[i][j] += rhs._values[i][j];
            }
        }
        return *this;
    }
private:

};

typedef Matrix<double, 4, 4> Matrix4f;

template <typename T, int N, int K, int M>
inline Matrix<T, N, M> operator *(const Matrix<T, N, K> &a, const Matrix<T, K, M> &b) {
    Matrix<T, N, M> c;
#pragma unroll
    for (int i = 0; i < N; ++i)
#pragma unroll
        for (int j = 0; j < M; ++j)
#pragma unroll
            for (int k = 0; k < K; ++k) {
                c[i][j] += a[i][k] * b[k][j];
            }
    return c;
}

template <typename T, int N, int M>
inline Matrix<T, N, M> operator +(const Matrix<T, N, M> &a, const Matrix<T, N, M> &b) {
    Matrix<T, N, M> c;
#pragma unroll
    for (int i = 0; i < N; ++i)
#pragma unroll
            for (int j = 0; j < M; ++j) {
                c[i][j] = a[i][j] + b[i][j];
            }
    return c;
}

template <typename T, int N>
inline T quadric(const Matrix<T, N, N> &a, const Vector<T, N> &x) {
    T res = T(0);
#pragma unroll
    for (int i = 0; i < N; ++i)
#pragma unroll
        for (int j = 0; j < N; ++j)
                res += a[i][j] * x[i] * x[j];
    return res;
}

template <typename T, int N, int M>
inline Vector<T, N> operator *(const Matrix<T, M, N> &a, const Vector<T, N> &x) {
    Vector<T, M> res;
#pragma unroll
    for (int i = 0; i < M; ++i)
#pragma unroll
            for (int j = 0; j < N; ++j)
                res[i] += a[i][j] * x[j];
    return res;
}

template <typename T, int N>
inline bool inverse(const Matrix<T, N, N> &ori, Matrix<T, N, N> &inv) {
    Matrix<T, N, N> mat = ori;

    for (int i = 0; i < N; ++i) {
        inv[i][i] = 1;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (sgn(mat[j][i]) == 0) continue;
            T lambda = mat[i][i] / mat[j][i];
            for (int k = 0; k < N; ++k) {
                mat[i][k] -= lambda * mat[j][k], std::swap(mat[i][k], mat[j][k]);
                inv[i][k] -= lambda * inv[j][k], std::swap(inv[i][k], inv[j][k]);
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        if (sgn(mat[i][i]) == 0) return false;
        for (int j =0; j < i; ++j) {
            T lambda = mat[j][i] / mat[i][i];
            for (int k = 0; k < N; ++k) {
                mat[j][k] -= lambda * mat[i][k];
                inv[j][k] -= lambda * inv[i][k];
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        T lambda = mat[i][i];
        for (int j = 0; j < N; ++j) {
            mat[i][j] /= lambda;
            inv[i][j] /= lambda;
        }
    }

    return true;
};

template <typename T, int N>
inline std::ostream &operator<<(std::ostream &os, const Vector<T, N> &vec) {
    os << "Vector(";

    bool is_first = true;
    for (int i = 0; i < N; ++i) {
        if (is_first) {
            os << vec._values[i];
            is_first = false;
        } else
            os << " " << vec._values[i];
    }

    os << ")";
    return os;
}

template <typename T, int N, int M>
inline std::ostream &operator<<(std::ostream &os, const Matrix<T, N, M> &mat) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0 ; j < M; ++j)
            os << mat[i][j] << " ";
        os << std::endl;
    }
    return os;
}

class MatrixFactory {
public:
    static inline Matrix4f get_quad_matrix(const Vector4f &character) {
        Matrix4f res;
#pragma unroll
        for (int i = 0; i < 4; ++i) {
            res[i][i] = character[i] * character[i];
#pragma unroll
            for (int j = 0; j < i; ++j)
                res[i][j] = res[j][i] = character[i] * character[j];
        }
        return res;
    }
};

} // End namespace decimation

#endif
