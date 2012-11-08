#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <params.h>

// K(x,x') = K(x) = sqrt(x')
void kernel(gsl_matrix* ak,
            gsl_vector* x,
            gsl_vector* xp)
{
  int i, j;

  for (i = 0; i < MAX; i++)
  {
    for (j = 0; j < MAX; j++)
    {
      gsl_matrix_set(ak, i, j, sqrt(gsl_vector_get(xp, j)));
    }
  }
  return;
}
