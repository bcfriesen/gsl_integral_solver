#ifndef GET_FRED2_INTERP_H
#define GET_FRED2_INTERP_H

#include <gsl/gsl_vector.h>

void get_fred2_interp(gsl_vector* grid,
                      double      a,
                      double      b,
                      gsl_vector* t,
                      gsl_vector* f,
                      gsl_vector* w,
                      gsl_vector* new_result);

#endif
