
#include "Walker.h"
#include "Function.h"
#include <iostream>

Walker::Walker(int mydim_parspace, double (*f)(std::vector<double>*, std::vector<double>*), std::vector<double>* mypdf_params, double propvar, bool myclimb)
{

	climb=myclimb;
	dim_parspace = mydim_parspace;
	acceptor_pdf=f;
	std::cout << "SHEEP " << std::endl;

	acceptor_pdf_params = new std::vector<double>;
	*acceptor_pdf_params=*mypdf_params;

	
	xt=std::vector<double>(dim_parspace,0);
	par_sigmas=std::vector<double>(dim_parspace,propvar);

    /* create a generator chosen by the 
       environment variable GSL_RNG_TYPE */

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

}

Walker::~Walker()
{
    gsl_rng_free (r);
	delete acceptor_pdf_params;
}


int Walker::jump(std::vector<double>* x)
{
	std::vector<double> xp(dim_parspace,0);
	jump_proposal(&xp,&xt);
	
	if (climb)
	{
		bool accept=jump_climb(&xp,&xt);
	}
	else
	{
		bool accept=jump_acceptor(&xp,&xt);
	}
	if (accept)
	{
		*x=xp;
		xt=*x;	
		return 1;
	}
	else
	{
		*x=xt;
		return 0;
	}
}

void Walker::jump_proposal(std::vector<double>* xp, std::vector<double>* x)
{
	double x_i;
	int i;
	for (i = 0; i < dim_parspace; i++) 
	  {
	    x_i = gsl_ran_gaussian_ziggurat (r, par_sigmas[i]);
		(*xp)[i] = (*x)[i]+x_i;
	  }
}	

bool Walker::jump_acceptor(std::vector<double>* xp, std::vector<double>* x)
{
	double p = acceptor_pdf(xp,acceptor_pdf_params)
		/acceptor_pdf(x,acceptor_pdf_params);
	
	bool ret;
	
	if (p>=1.0)
	{
		ret=true;
	}
	else
	{
		switch(gsl_ran_bernoulli(r,p))
		{
			case(1):
				ret = true;
				break;
			case(0):
				ret = false;
				break;
		}
		
	}
	return ret;
}

bool Walker::jump_climb(std::vector<double>* xp, std::vector<double>* x)
{
	double p = acceptor_pdf(xp,acceptor_pdf_params)
		/acceptor_pdf(x,acceptor_pdf_params);
	
	bool ret;
	
	if (p>=1.0)
	{
		ret=true;
	}
	else
	{
		ret=false;
	}
	return ret;
}

void Walker::initialize_position(std::vector<double>* x0)
{
	double x_i;
	int i;
	for (i = 0; i < dim_parspace; i++) 
	  {
	    x_i = gsl_ran_gaussian_ziggurat (r, init_sigma);
		xt[i] = x_i;
	  }
	*x0 = xt;
}	

void Walker::writeout_pdf(double first, double last, int npoints, std::string filename)
{

	std::vector<double> xvec(npoints,0);
	std::vector<double> yvec(npoints,0);

	double delta=last-first;

	std::vector<double> pdf_argvec(1,0);

    if (npoints==1)
    {
		xvec[0]=first;
    }
    else
    {
        int i;
        double rfrac;
        for(i=0; i<npoints; i++)
        {
            rfrac=((double)i/((double)npoints - 1));
            xvec[i]=first+delta*rfrac;
			pdf_argvec[0]=xvec[i];
			yvec[i]=acceptor_pdf(&pdf_argvec,acceptor_pdf_params);
        }
    }

	std::vector<std::vector<double>> write_arg;
	write_arg.push_back(xvec);
	write_arg.push_back(yvec);
	write_vecs(write_arg,filename);
}




