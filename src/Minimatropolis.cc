
#include <vector>
#include <iostream>
#include "Function.h"
#include "Walker.h"
#include "WriteOut.h"
#include <time.h>






int main()
{
	int nsteps=100;
	
	std::vector<double> metapdf_params(1,0.5);

	Walker meta_ranger(1,acc_pdf,&metapdf_params,0.25,true);
	std::vector<double> sigma(1,0);
	meta_ranger.initialize_position(&sigma);
	std::cout << "Position: " << sigma[0] << std::endl;

	double sum=0;
	double varsum=0;
	int n_acc=0;
	int step;
	for (step=0; step<nsteps; step++)
	{
		meta_ranger.jump(&sigma);
	}


    std::vector<double> pdf_params(3,1);
    pdf_params[0]=0;
	std::cout << "SIGMAFIN: " << sigma[0] << std::endl;

    nsteps=1000000;
    
    Walker texas_ranger(1,gaussian,&pdf_params,sigma[0],false);
    //Walker texas_ranger(1,gaussian,&pdf_params,1.0,false);
    std::vector<double> x(1,0);
    texas_ranger.initialize_position(&x);

	std::vector<double> positions(nsteps,0);

    clock_t start,end;
    double cpu_time_used;
    


 start=clock();

        
	sum=0;
	n_acc=0;
    for (step=0; step<nsteps; step++)
    { 
        n_acc+=texas_ranger.jump(&x);
		positions[step]=sigma[0];
		sum+=pow((x[0]),4);
		//std::cout << "Position: " << x[0] << std::endl;
    }  

    end=clock();
    cpu_time_used = ((double) (end-start))/CLOCKS_PER_SEC;
    //precalculate_cov(cov_mat_pre,sample_radii,k_cutoff);
    //precalculate_cov(cov_mat,sample_radii,k_cutoff);
    std::cout << "Chain Time:  " << cpu_time_used << "   (Seconds)" << std::endl;

	double mean=sum/((double)nsteps);
	for (step=0; step<nsteps; step++)
	{
		varsum+=pow((pow(positions[step],4)-mean),2);
		//std::cout << "Position: " << x[0] << std::endl;
	}

	double sigmaf=sqrt(varsum/((double)(nsteps-1)));

	std::cout << "Integral Estimate: " << mean << std::endl;
	std::cout << "Error Estimate: " << sigmaf/(sqrt((double)nsteps)) << std::endl;
	std::cout << "Error: " << fabs(mean-3) << std::endl;

    double acc_ratio = ((double)n_acc)/((double)nsteps);
        
    std::cout << "SIGMA: " << (sigma)[0] <<std::endl;
    std::cout << "ACC RATIO: " << acc_ratio <<std::endl;


	std::vector<std::vector<double>> positions_container;
	positions_container.push_back(positions);

	write_vecs(positions_container,std::string("output/positions.txt"));

	
	texas_ranger.writeout_pdf(-5,5,100,std::string("output/pdfout.txt"));
}
	
