#include <gsl/gsl_matrix.h>
#include <vector>
#include "Function.h"

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

double gaussian(std::vector<double>* x, std::vector<double>* params)
{
    return (*params)[2]*(1/((*params)[1]*sqrt(2*pi)))*exp(-0.5*pow(((*x)[0]-(*params)[0])/(*params)[1],2));
}


