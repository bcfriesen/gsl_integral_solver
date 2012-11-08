#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <params.h>
#include <aux_g.h>
#include <kernel.h>

void fred2_solver(double a, double b, gsl_vector *t, gsl_vector *f, gsl_vector *w)
{
  gsl_integration_glfixed_table *glgrid = gsl_integration_glfixed_table_alloc(MAX);
  gsl_permutation *p = gsl_permutation_alloc(MAX);
  gsl_matrix *lhs    = gsl_matrix_alloc(MAX, MAX);
  gsl_matrix *ktilde = gsl_matrix_alloc(MAX, MAX);
  gsl_matrix *ak     = gsl_matrix_alloc(MAX, MAX);
  gsl_vector *g = gsl_vector_alloc(MAX);
  int i, j, error, s;
  double ptsi, wghtsi;

  // set Gauss-Legendre integration points and weights
  for (i = 0; i < MAX; i++)
  {
    error = gsl_integration_glfixed_point(a, b, i, &ptsi, &wghtsi, glgrid);
    gsl_vector_set(t, i, ptsi);
    gsl_vector_set(w, i, wghtsi);
  }

  kernel(ak, t, t);
  aux_g(g, t);

  // fill in unit matrix first
  gsl_matrix_set_identity(lhs);

  for (i = 0; i < MAX; i++)
  {
    for (j = 0; j < MAX; j++)
    {
      gsl_matrix_set(ktilde, i, j, gsl_matrix_get(ak, i, j)*gsl_vector_get(w, j));
    }
  }

  // set up LHS matrix
  error = gsl_matrix_sub(lhs, ktilde);

  gsl_linalg_LU_decomp(lhs, p, &s);
  gsl_linalg_LU_solve(lhs, p, g, f);

  gsl_integration_glfixed_table_free(glgrid);
  gsl_permutation_free(p);
  gsl_matrix_free(lhs);
  gsl_matrix_free(ktilde);
  gsl_vector_free(g);
  gsl_matrix_free(ak);
  return;
}


