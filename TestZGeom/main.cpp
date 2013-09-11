#include <cstdlib>
#include <iostream>
#include <ZGeom/VecN.h>

using ZGeom::VecN;

int main()
{
	VecN<double> v1;
	v1.resize(10, 1.5);

	for (VecN<double>::iterator iter = v1.begin(); iter != v1.end(); ++iter) {
		std::cout << *iter << std::endl;
	}

	for(double& a : v1) a *= 2; 
	for(double a : v1) std::cout << a << std::endl;

	std::system("PAUSE");
	std::exit(0);
}