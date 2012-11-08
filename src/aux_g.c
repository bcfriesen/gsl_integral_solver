#include <gsl/gsl_vector.h>
#include <params.h>
#include <math.h>

// g(x) = 2x^2 - 1 - (4/7)*(5^(7/2) - 2^(7/2)) - (2/3)*(5^(3/2) - 2^(3/2))
void aux_g(gsl_vector *g, gsl_vector *x)
{
  int i;

  for (i = 0; i < MAX; i++)
  {
    gsl_vector_set(g, i, 2.0*pow(gsl_vector_get(x, i), 2.0) - 1.0 -
                         (4.0/7.0)*(pow(5.0, (7.0/2.0)) - pow(2.0, (7.0/2.0))) -
			 (2.0/3.0)*(pow(5.0, (3.0/2.0)) - pow(2.0, (3.0/2.0))));
  }
  return;
}
