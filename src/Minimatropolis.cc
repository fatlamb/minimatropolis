#include "Walker.h"
#include "Function.h"

#include <vector>
#include <iostream>


int main()
{
	std::vector<double> pdf_params(3,1);
	pdf_params[0]=0;

	Walker falf(1,gaussian,pdf_params);
	std::vector<double> x(1,0);
	falf.initialize_position(&x);
	std::cout << "Position: " << x[0] << std::endl;
}
	
