#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <kernel.h>
#include <aux_g.h>
#include <params.h>

void get_fred2_interp(gsl_vector* grid,
                      double      a,
                      double      b,
                      gsl_vector* t,
                      gsl_vector* f,
                      gsl_vector* w,
                      gsl_vector* new_result)
{
  gsl_matrix* ak = gsl_matrix_alloc(MAX, MAX);
  gsl_vector* g = gsl_vector_alloc(MAX);
  int i, j;
  double sum;

  kernel(ak, grid, t);
  aux_g(g, grid);

  for (i = 0; i < MAX; i++)
  {
    sum = 0.0;
    for (j = 0; j < MAX; j++)
    {
      sum += gsl_matrix_get(ak, i, j) * gsl_vector_get(w, j) * gsl_vector_get(f, j);
    }
    gsl_vector_set(new_result, i, gsl_vector_get(g, i) + sum);
  }

  gsl_matrix_free(ak);
  gsl_vector_free(g);
  return;
}
