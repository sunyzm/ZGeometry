#include "arithmetic.h"
#include <stdexcept>
#include <amp.h>

double ZGeom::triArea( double e1, double e2, double e3 )
{
	if(!validAsTriangle(e1, e2, e3)) throw std::runtime_error("Invalid triangle side lengths");

	double p = (e1 + e2 + e3) / 2.0;
	return std::sqrt(p*(p-e1)*(p-e2)*(p-e3));
}

double ZGeom::cosTriSides( double e1, double e2, double e3 )
{
	if(!validAsTriangle(e1, e2, e3)) throw std::runtime_error("Invalid triangle side lengths");

	return (e1*e1 + e2*e2 - e3*e3) / (2.0*e1*e2);
}

void ZGeom::triangleCot( double a, double b, double c, double &cotan_a, double &cotan_c )
{
	if(!validAsTriangle(a, b, c)) throw std::runtime_error("Invalid triangle side lengths");

	double cosa = (b*b+c*c-a*a)/(2.0*b*c);
	double cosc = (b*b+a*a-c*c)/(2.0*b*a);
	cotan_a = cosa / std::sqrt(1-cosa*cosa);
	cotan_c = cosc / std::sqrt(1-cosc*cosc);
}

double ZGeom::calMixedTriArea( double a, double b, double c )
{
	double cosa = (b*b+c*c-a*a)/(2.0*b*c);
	double cosc = (b*b+a*a-c*c)/(2.0*b*a);
	double cotan_a = cosa / sqrt(1.0 - cosa*cosa);
	double cotan_c = cosc / sqrt(1.0 - cosc*cosc);

	if (a*a + c*c < b*b)
	{
		double s = (a+b+c)/2.0;
		return sqrt(s*(s-a)*(s-b)*(s-c))/2.0;
	}
	else if (a*a + b*b < c*c || b*b + c*c < a*a)
	{
		double s = (a+b+c)/2.0;
		return sqrt(s*(s-a)*(s-b)*(s-c))/4.0;
	}
	else 
	{
		return (a*a*cotan_a + c*c*cotan_c)/8.0;
	}
}

void ZGeom::quadricFormAMP(int dim1, int dim2, double* mat1, double* diag, double *matResult)
{
	using namespace Concurrency;
	// Y=X'*Q*X
	array_view<double, 2> X(dim2, dim1, mat1);
	array_view<double, 1> Q(dim2, diag);
	array_view<double, 2> Y(dim1, dim1, matResult);	
	Y.discard_data();

	parallel_for_each(Y.extent, [=](index<2> idx) restrict(amp) {
		int row = idx[0], col = idx[1];
		if (row >= col) {
			Y[idx] = 0;
			for (int k = 0; k < dim2; ++k) {
				Y[idx] += Q(k) * X(k, row) * X(k, col);
			}
			Y(col, row) = Y[idx];
		}	
	});
	Y.synchronize();
}

void ZGeom::matVecMulAMP( int dim1, int dim2, double *mat, double *vec, double *vResult )
{
	using namespace Concurrency;
	// Y = M*X
	array_view<double, 2> M(dim1, dim2, mat);
	array_view<double, 1> X(dim2, vec);
	array_view<double, 1> Y(dim1, vResult);
	Y.discard_data();

	parallel_for_each(Y.extent, [=](index<1> idx) restrict(amp) {
		int row = idx[0];
		Y[idx] = 0;
		for (int k = 0; k < dim2; ++k) {
			Y[idx] += M(row, k) * X(k);
		}
	});

	Y.synchronize();
}
