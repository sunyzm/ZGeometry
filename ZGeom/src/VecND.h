#ifndef ZGEOM_VECND_H
#define ZGEOM_VECND_H

#include <cassert>

namespace ZGeom
{

template<typename T>
class VecND
{
public:
	VecND() : vec(nullptr), dim(0) {}
	VecND(int n) { resize(n); }
	~VecND() { delete []vec; }
	T& operator[] (int i) { return vec[i]; }
	void resize(int n);
	friend VecND<T> operator + (const VecND<T>& v1, const VecND<T>& v2);
	friend VecND<T> operator * (const VecND<T>& v1, double t);
	friend VecND<T> operator * (double t, const VecND<T>& v1);
private:
	T *vec;
	int dim;
};

template<typename T>
void VecND<T>::resize(int n)
{
	if (vec) delete []vec;
	this->dim = n;
	vec = new T[dim];
}

template<typename T>
VecND<T> operator + (const VecND<T>& v1, const VecND<T>& v2)
{
	assert(v1.dim == v2.dim);
	VecND v3(v1.dim);
	for (int k = 0; k < v3.dim; ++k) v3.vec[k] = v1.vec[k] + v2.vec[k];
	return v3;
}

template<typename T>
VecND<T> operator * (const VecND<T>& v1, double t)
{
	VecND v2(v1.dim);
	for (int k = 0; k < v2.dim; ++k) v2.vec[k] = v1.vec[k] * t;
	return v2;
}

template<typename T>
VecND<T> operator * (double t, const VecND<T>& v1)
{
	return v1 * t;
}

}// end of namespace



#endif