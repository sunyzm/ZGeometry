#include <cstdlib>
#include <iostream>
#include <random>
#include <numeric>
#include <fstream>
#include <ppl.h>
#include <ZGeom/VecN.h>
#include <ZUtil/timer.h>
#include <mkl.h>

using ZGeom::VecNd;

void test_dgemv()
{
	double a[]={1,2,3,9,8,7};
	double x[]={4,5,6};
	double y[2] = {0};
	char trans = 'N';
	double alpha = 1, beta = 0;
	int incx = 1, incy = 1;
	int m = 2, n = 3;
	//dgemv(&trans, &m, &n, &alpha, a, &m, x, &incx, &beta, y, &incy);	// default to column-major!
	cblas_dgemv(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, m, n, alpha, a, n, x, incx, beta, y, incy);
	std::cout << "\ny(:)" << y[0] << ' ' << y[1] << std::endl;
};

int main()
{
	CStopWatch timer;
	const int testSize(100000);
	VecNd v1(testSize), v2(testSize);
	
	double *va = new double[testSize], *vb = new double[testSize], *vc = new double[testSize];

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(1,2);	
	for (int i = 0; i < testSize; ++i) {
		va[i] = v1[i] = dis(gen);
		vb[i] = v2[i] = dis(gen);
	}

	double prod1(0), prod2(0), prod3(0);
	timer.startTimer();
	prod1 = std::inner_product(va, va + testSize, vb, 0.0);
	timer.stopTimer("Method 1: ");
	timer.startTimer();
	prod2 = std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
	timer.stopTimer("Method 2: ");
	timer.startTimer();
	prod3 = v1.dot(v2);
	timer.stopTimer("Method 4: ");

	std::cout << "inner product 1: " << prod1 << "; inner product 2: " << prod2 << std::endl;
	std::cout << "inner product 3: " << prod3 << std::endl;

	timer.startTimer();
	Concurrency::parallel_for(0, testSize, [&](int i){
		vc[i] = va[i] * vb[i];
	});
	timer.stopTimer("Product 1: ");
	timer.startTimer();
	for(int i = 0; i < testSize; ++i) {
		vc[i] = va[i] * vb[i];
	}
	timer.stopTimer("Product 2: ");

	test_dgemv();

	std::system("PAUSE");
	std::exit(0);
}