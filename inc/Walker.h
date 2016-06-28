#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


class Walker
{

	public:
		Walker(int mydim_parspace, double (*f)(std::vector<double>*, std::vector<double>*) ,std::vector<double>* mypdf_params);
		~Walker(void);
		void initialize_position(std::vector<double>* x0); 
		void jump(std::vector<double>* x);

	private:

		void jump_proposal(std::vector<double>* xp, std::vector<double>* x);
		//double proposal_pdf(std::vector<double>* xp std::vector<double> x, double* params);
		//double acceptor_pdf(std::vector<double>* xp std::vector<double> x, double* params);
		bool jump_acceptor(std::vector<double>* xp, std::vector<double>* x);
		double acceptor_pdf(std::vector<double>* x, std::vector<double>* params);

		std::vector<double>* xt;
		double init_sigma=1.0;
		std::vector<double>* par_sigmas;
		std::vector<double>* acceptor_pdf_params;		
		int dim_parspace;


	    const gsl_rng_type * T;
    	gsl_rng * r;


}
