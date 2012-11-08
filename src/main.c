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
  const double a = 0.0, b = 4.0;
  double i;
  gsl_vector* tmp        = gsl_vector_alloc(MAX);
  gsl_vector* grid       = gsl_vector_alloc(MAX);
  gsl_vector* t          = gsl_vector_alloc(MAX);
  gsl_vector* f          = gsl_vector_alloc(MAX);
  gsl_vector* w          = gsl_vector_alloc(MAX);
  gsl_vector* new_result = gsl_vector_alloc(MAX);
  gsl_vector* analytic   = gsl_vector_alloc(MAX);

  // evenly spaced grid
  for (i = 0; i < MAX; i++)
  {
    gsl_vector_set(grid, i, a + (double)i*((b - a)/(double)MAX));
    // fill in analytic result
    gsl_vector_set(analytic, i, 2.0*gsl_vector_get(grid, i) - 1.0);
  }

  // first solve integral equation at Gauss-Legendre quadrature points
  fred2_solver(a, b, t, f, w);
  // now use Nystrom interpolation to use the points WE want
  get_fred2_interp(grid, a, b, t, f, w, new_result);

  FILE *fp;
  fp = fopen("test.out", "w");

  fprintf(fp, "%25s %25s %25s\n", "x", "f(x) (expected)", "f(x) (calculated)");
  for (i = 0; i < MAX; i++)
  {
    fprintf(fp, "%25.9e %25.9e %25.9e\n",
            gsl_vector_get(grid, i),
            gsl_vector_get(analytic, i),
            gsl_vector_get(new_result, i));
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
