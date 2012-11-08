#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <params.h>
#include <kernel.h>
#include <fred2_solver.h>
#include <aux_g.h>
#include <fred2_interp.h>

const int MAX = 500;

int main (void)
{
  double a = 2.0, b = 5.0, i;
  gsl_vector *tmp        = gsl_vector_alloc(MAX);
  gsl_vector *grid       = gsl_vector_alloc(MAX);
  gsl_vector *t          = gsl_vector_alloc(MAX);
  gsl_vector *f          = gsl_vector_alloc(MAX);
  gsl_vector *w          = gsl_vector_alloc(MAX);
  gsl_vector *new_result = gsl_vector_alloc(MAX);

  // evenly spaced grid
  for (i = 0; i < MAX; i++)
  {
    gsl_vector_set(grid, i, a + (double)i*((b - a)/MAX));
  }

  // first solve integral equation at Gauss-Legendre quadrature points
  fred2_solver(a, b, t, f, w);
  // now use Nystrom interpolation to use the points WE want
  get_fred2_interp(grid, a, b, t, f, w, new_result);

  FILE *fp;
  fp = fopen("test2.out", "w");

  fprintf(fp, "%s %s %s\n", "x", "x (expected)", "x (calculated)");
  for (i = 0; i < MAX; i++)
  {
    fprintf(fp, "%18.5f %18.5f\n", gsl_vector_get(grid, i), gsl_vector_get(new_result, i));
  }

  fclose(fp);

  gsl_vector_free(tmp);
  gsl_vector_free(grid);
  gsl_vector_free(t);
  gsl_vector_free(f);
  gsl_vector_free(w);
  gsl_vector_free(new_result);
  return 0;
}

