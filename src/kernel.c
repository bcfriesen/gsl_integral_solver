#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <params.h>

void kernel(gsl_matrix* ak,
            gsl_vector* x,
            gsl_vector* xp)
{
  int i, j;

  for (i = 0; i < MAX; i++)
  {
    for (j = 0; j < MAX; j++)
    {
      // K(x, x') = 3x
      gsl_matrix_set(ak, i, j, 3.0*gsl_vector_get(x, i));
    }
  }
  return;
}
