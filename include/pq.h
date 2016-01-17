/**
 * File   : pq.h
 * Author : Jiayuan Mao
 * Email  : mjy14@mails.tsinghua.edu.cn
 * Date   : 2016-01-02 17:55:23
 * This file is part of the school project Decimation of course
 * ``Advanced Computational Geometry''.
 **/


#ifndef _INCLUDE_PQ_H_
#define _INCLUDE_PQ_H_

#include <algorithm>
#include <iostream>

namespace decimation {

using std::make_heap;
using std::push_heap;
using std::pop_heap;

template <typename T>
class PriorityQueue {
public:
	PriorityQueue(void) : _heap(NULL), _size(0), _max_size(0) {
		_realloc(__INIT_SIZE__);
	}

	inline void push(const T &element) {
		_assert_size(_size + 1);
		_heap[_size++] = element;
		push_heap(_heap, _heap + _size);
	}
	inline void pop(void) {
		pop_heap(_heap, _heap + _size);
		_assert_size(--_size);
	}
	inline T top(void) const {
		return _heap[0];
	}
	inline int size(void) const {
		return _size;
	}
	inline bool empty(void) const {
		return _size == 0;
	}
private:
	inline void _realloc(int new_size) {
		T* new_heap = new T[new_size];
		if (_heap) {
			memcpy(new_heap, _heap, sizeof(T) * _max_size);
			delete []_heap;
		}
		_heap = new_heap;
		_max_size = new_size;
	}

	inline void _assert_size(int sz) {
		if (sz > _max_size) { 
			_realloc(_max_size * __INC_RATIO__);
		}
		// if (sz < _max_size * __DEC_RATIO__ && sz >= __INIT_SIZE__) {
		// 	_realloc(_max_size * __DEC_RATIO__);
		// }
	}

	const int __INIT_SIZE__ = 128;
	const double __INC_RATIO__ = 2;
	const double __DEC_RATIO__ = 0.5;
	T   *_heap;
	int _size, _max_size;
};

} // End namespace decimation

#endif
