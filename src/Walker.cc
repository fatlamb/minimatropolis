
#include "Walker.h"
#include "Function.h"

Walker::Walker(int mydim_parspace, double (*f)(std::vector<double>*, std::vector<double>*) , std::vector<double>* mypdf_params)
{
	dim_parspace = mydim_parspace;
	acceptor_pdf=mypdf;
	*acceptor_pdf_params=*mypdf_params;
	
	xt=new std::vector<double>(dim_parspace,0);
	par_sigmas=new std::vector<double>(dim_parspace,1.0);

    /* create a generator chosen by the 
       environment variable GSL_RNG_TYPE */

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

}

Walker::~Walker()
{
    gsl_rng_free (r);
	delete xt;
	delete par_sigmas;
	delete acceptor_pdf_params;
}


void Walker::jump(std::vector<double>* x)
{
	std::vector<double> xp(dim_parspace,0);
	jump_proposal(&xp,xt);
	bool accept=jump_acceptor(&xp,&xt);
	if (accept)
	{
		*x=*xp;
	}
	else
	{
		*x=*xt;
	}
}

void Walker::jump_proposal(std::vector<double>* xp, std::vector<double>* x)
{
	double x_i;
	int i;
	for (i = 0; i < dim_parspace; i++) 
	  {
	    x_i = gsl_ran_gaussian_ziggurat (r, par_sigmas[i]);
		*xp[i] = *x[i]+x_i;
	  }
}	

bool Walker::jump_acceptor(std::vector<double>* xp, std::vector<double>* x)
{
	double p = acceptor_pdf(xp,acceptor_pdf_params)
		/acceptor_pdf(x,acceptor_pdf_params);
	
	if (p>=1.0)
	{
		return true;
	}
	else
	{
		switch(gsl_ran_bernoulli(r,p))
		{
			case(1):
				return true;
			case(0):
				return false;
		}
		
	}
}

void Walker::initialize_position(std::vector<double>* x0)
{
	double x_i;
	int i;
	for (i = 0; i < dim_parspace; i++) 
	  {
	    x_i = gsl_ran_gaussian_ziggurat (r, init_sigma);
		*xt[i] = x_i;
	  }
	*x0 = *xt;
}	





