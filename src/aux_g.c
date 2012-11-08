#include <gsl/gsl_vector.h>
#include <params.h>
#include <math.h>

// g(x) = -34x - 1
void aux_g(gsl_vector* g, gsl_vector* x)
{
  int i;

  for (i = 0; i < MAX; i++)
  {
    gsl_vector_set(g, i, -34.0*gsl_vector_get(x, i) - 1.0);
  }
  return;
}
