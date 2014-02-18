#ifndef ZGEOM_VEC4_H
#define ZGEOM_VEC4_H

namespace ZGeom
{
	template<typename T>
	class Vec4
	{
	public:
		Vec4() { for (int i = 0; i < 4; ++i) mVal[i] = 0.0; }
		Vec4(T a, T b, T c, T d) { mVal[0] = a; mVal[1] = b; mVal[2] = c; mVal[3] = d; }
		Vec4(T* data) { for (int i = 0; i < 4; ++i) mVal[i] = data[i]; }
		
		template<typename F> Vec4(const Vec4<F>& vec);		
		template<typename F> const Vec4<T>& operator = (const Vec4<F>& vec);

		T operator[] (int i) const { return mVal[i]; }
		T& operator[] (int idx) { return mVal[idx]; }
		T* c_ptr() const { return mVal; }

	protected:
		T mVal[4];
	};

	typedef Vec4<float> Vec4s;
	typedef Vec4<double> Vec4d;

	template<typename T>
	template<typename F>
	Vec4<T>::Vec4(const Vec4<F>& vec)
	{
		for (int i = 0; i < 4; ++i) 
			mVal[i] = static_cast<T>(vec.mVal[i]);
	}

	template<typename T>
	template<typename F>
	const Vec4<T>& Vec4<T>::operator = (const Vec4<F>& vec)
	{
		for (int i = 0; i < 4; ++i) 
			mVal[i] = static_cast<T>(vec.mVal[i]);
		return *this;
	}

} //end of namespace ZGeom

#endif