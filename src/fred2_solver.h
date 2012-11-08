#ifndef FRED2_SOLVER_H
#define FRED2_SOLVER_H

#include <gsl/gsl_vector.h>

void fred2_solver(double      a,
                  double      b,
                  gsl_vector* t,
                  gsl_vector* f,
                  gsl_vector* w);

#endif
