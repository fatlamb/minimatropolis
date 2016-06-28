#ifndef __FUNCTION_H_INCLUDED__
#define __FUNCTION_H_INCLUDED__

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <vector>

const double pi=acos(-1.0);

double Q(gsl_vector* phi, gsl_vector* lambda, gsl_vector* sigma, gsl_matrix* eta);
double gaussian(std::vector<double>* x, std::vector<double>* params);


#endif // __FUNCTION_H_INCLUDED__
