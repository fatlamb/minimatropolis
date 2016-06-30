
#include <vector>
#include <iostream>
#include "Function.h"
#include "Walker.h"
#include "WriteOut.h"


double acc_pdf(std::vector<double>* x, std::vector<double>* params)

	//FIXME send pdf params as params
	std::vector<double> pdf_params(3,1);
    pdf_params[0]=0;
	int nsteps=1000;

    Walker texas_ranger(1,gaussian,&pdf_params,(*x)[0],false);
    std::vector<double> x(1,0);
    texas_ranger.initialize_position(&x);

    int n_acc=0;
    int step;
    for (step=0; step<nsteps; step++)
    {
        n_acc+=texas_ranger.jump(&x);
    }

	double acc_ratio = ((double)n_acc)/((double)nsteps);

    std::cout << "ACC RATIO: " << acc_ratio <<std::endl;

	double acc_diff = fabs(acc_ratio-(*params)[0]);

	return acc_diff;



int main()
{
	int nsteps=1000000;
	
	std::vector<double> pdf_params(3,1);
	pdf_params[0]=0;

	Walker texas_ranger(1,gaussian,&pdf_params,1);
	std::vector<double> x(1,0);
	texas_ranger.initialize_position(&x);
	std::cout << "Position: " << x[0] << std::endl;


	
	std::vector<double> positions(nsteps,0);

	double sum=0;
	double varsum=0;
	int n_acc=0;
	int step;
	for (step=0; step<nsteps; step++)
	{
		n_acc+=texas_ranger.jump(&x);
		positions[step]=x[0];
		sum+=pow((x[0]),4);
		//std::cout << "Position: " << x[0] << std::endl;
	}

	std::cout << "ACC RATIO: " << ((double)n_acc)/((double)nsteps)<<std::endl;

	double mean=sum/((double)nsteps);
	for (step=0; step<nsteps; step++)
	{
		varsum+=pow((pow(positions[step],4)-mean),2);
		//std::cout << "Position: " << x[0] << std::endl;
	}

	double sigma=sqrt(varsum/((double)(nsteps-1)));

	std::cout << "Integral Estimate: " << mean << std::endl;
	std::cout << "Error Estimate: " << sigma/(sqrt((double)nsteps)) << std::endl;
	std::cout << "Error: " << fabs(mean-3) << std::endl;

	std::vector<std::vector<double>> positions_container;
	positions_container.push_back(positions);

	write_vecs(positions_container,std::string("output/positions.txt"));

	
	texas_ranger.writeout_pdf(-5,5,100,std::string("output/pdfout.txt"));

}
	
