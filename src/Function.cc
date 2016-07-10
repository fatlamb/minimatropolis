#include <gsl/gsl_matrix.h>
#include <vector>
#include <iostream>
#include "Function.h"
#include "Walker.h"

//const double pi=acos(-1.0);


double Q(gsl_vector* phi, gsl_vector* lambda, gsl_vector* sigma, gsl_matrix* eta)
{
	double T1,T2,T3,T4;
	int i,alpha;
	int d_fields=phi->size;
	int d_space = 3;

	double sigma0_sq = pow(gsl_vector_get(sigma,0),2);
	double sigma1_sq = pow(gsl_vector_get(sigma,1),2);
	double sigma2_sq = pow(gsl_vector_get(sigma,2),2);

	double gamma_sq = pow(sigma1_sq,2)/(sigma0_sq*sigma2_sq);


	T1 = 0;
	for (alpha=0;alpha<d_fields;alpha++)
	{
		T1+=pow(gsl_vector_get(phi,alpha),2);
	}
	T1/=sigma0_sq;

	T2 = 0;
	for (alpha=0;alpha<d_fields;alpha++)
	{
		for (i=0;i<d_space;i++)
		{
			T2+= pow(gsl_matrix_get(eta,alpha,i),2);
		}
	}
	T2*=(3.0/sigma1_sq);
	
	T3 = 0;
	for (i=0;i<d_space;i++)
	{
		T3+=pow(gsl_vector_get(lambda,i),2);
	}
	T3*=3.0;

	double sum_lam=0;
	for (i=0;i<d_space;i++)
	{
		sum_lam+=gsl_vector_get(lambda,i);
	}
	T3-=pow(sum_lam,2);	

	T3*=(5*gamma_sq*sigma0_sq)/(2*pow(sigma1_sq,2));


	double temp;
	T4=0;
	temp=sum_lam*sigma0_sq/sigma1_sq;
	temp+=gsl_vector_get(phi,0);
	T4=temp*gamma_sq/(sigma0_sq*(1-gamma_sq));

	return 0.5*(T1+T2+T3+T4);
}

double acc_pdf(std::vector<double>* sigma, std::vector<double>* params) 
{ 
    //FIXME send pdf params as params 
    std::vector<double> pdf_params(3,1); 
    pdf_params[0]=0; 
    int nsteps=1000; 
 
    Walker texas_ranger(1,gaussian,&pdf_params,(*sigma)[0],false); 
    std::vector<double> x(1,0); 
    texas_ranger.initialize_position(&x); 
 
    int n_acc=0; 
    int step; 
    for (step=0; step<nsteps; step++) 
    { 
        n_acc+=texas_ranger.jump(&x); 
    } 
 
    double acc_ratio = ((double)n_acc)/((double)nsteps); 
 
    std::cout << "SIGMA: " << (*sigma)[0] <<std::endl; 
    std::cout << "ACC RATIO: " << acc_ratio <<std::endl; 
 
    double acc_diff = fabs(acc_ratio-(*params)[0]); 
 
    return 1/acc_diff; 
} 



double gaussian(std::vector<double>* x, std::vector<double>* params)
{
    return (*params)[2]*(1/((*params)[1]*sqrt(2*pi)))*exp(-0.5*pow(((*x)[0]-(*params)[0])/(*params)[1],2));
}


double q_integrand(std::vector<double>* x, std::vector<double>* params)
{
	double nu,sig0,sig1,kap;
	nu=(*params)[0];
	sig0=(*params)[1];
	sig1=(*params)[2];
	kap=(*params)[3];

	double lamsum = std::accumulate(x->begin(),x->begin()+3,0);
	double lamsqsum = std::inner_product(x->begin(),x->begin()+3,x->begin(),0);
	std::vector<double> sqeta(3,0);
	unsigned int i;
	for (i=0; i<3; i++)
	{
		sqeta[i]=pow((*x)[i+3],2);
	}

	double sqetasum = std::accumulate(sqeta.begin(),sqeta.end(),0);
	double lamsqeta_sum = std::inner_product(sqeta.begin(),sqeta.end(),x->begin(),0);	

	double Q=0;

	Q += 0.5*pow(nu,2);
	Q += (3.0/(2.0*pow(sig1,2)))*squetasum;
	Q += ((5*pow(sig0,2))/(4*pow(kap,2)*pow(sig1,4))) * ( 3*(lamsqsum - (2/(nu*sig0))*lamsqeta_sum + pow((sqetasum/(nu*sig0)),2));
		- pow(lamsum - (1/(nu*sig0))*sqetasum,2));	
	Q += ((pow(sig0,2))/(2*(pow(kap,2)-1)*pow(sig1,4))) * pow( ((nu*pow(sig1,2))/sig0) + lamsum - (1/(nu*sig0))*sqetasum, 2);
	Q += (-1.0)*((pow(sig0,2))/(pow(kap,2)*pow(sig1,4)))*(((nu*pow(sig1,2))/sig0) + lamsum + (1/(nu*sig0))*sqetasum); 

	double detfac = ((*x)[0]-(*x)[1])*((*x)[0]-(*x)[2])* ((*x)[1]-(*x)[2])*(*x)[0]*(*x)[1]*(*x)[2];

	return detfac*exp((-1.0)*Q);

