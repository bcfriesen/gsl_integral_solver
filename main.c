#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

// size of vectors/matrices
#define MAX 500

void get_kernel(gsl_matrix *ak, gsl_vector *x, gsl_vector *xp);
void fred2_solver(double a, double b, gsl_vector *t, gsl_vector *f,
                  gsl_vector *w);
void get_g(gsl_vector *g, gsl_vector *eks);
void get_fred2_interp(gsl_vector *grid, double a, double b, gsl_vector *t,
                      gsl_vector *f, gsl_vector *w, gsl_vector *new_result);

int main (void)
{
  double a = 0.0, b = 4.0, i;
  gsl_vector *tmp        = gsl_vector_alloc(MAX);
  gsl_vector *grid       = gsl_vector_alloc(MAX);
  gsl_vector *t          = gsl_vector_alloc(MAX);
  gsl_vector *f          = gsl_vector_alloc(MAX);
  gsl_vector *w          = gsl_vector_alloc(MAX);
  gsl_vector *new_result = gsl_vector_alloc(MAX);

  for (i = 0; i < MAX; i++)
  {
    gsl_vector_set(grid, i, a + (double)i*((b - a)/MAX));
  }

  // first solve integral equation at Gauss-Legendre quadrature points
  fred2_solver(a, b, t, f, w);

  // now use Nystrom interpolation to use the points WE want
  get_fred2_interp(grid, a, b, t, f, w, new_result);

  for (i = 0; i < MAX; i++)
  {
    printf("%18.5f %18.5f\n", gsl_vector_get(grid, i), gsl_vector_get(new_result, i));
  }

  gsl_vector_free(tmp);
  gsl_vector_free(grid);
  gsl_vector_free(t);
  gsl_vector_free(f);
  gsl_vector_free(w);
  gsl_vector_free(new_result);

  return 0;
}


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

  get_kernel(ak, t, t);
  get_g(g, t);

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


void get_fred2_interp(gsl_vector *grid, double a, double b, gsl_vector *t,
                      gsl_vector *f, gsl_vector *w, gsl_vector *new_result)
{
  gsl_matrix *ak = gsl_matrix_alloc(MAX, MAX);
  gsl_vector *g = gsl_vector_alloc(MAX);
  int i, j;
  double sum;

  get_kernel(ak, grid, t);
  get_g(g, grid);

  for (i = 0; i < MAX; i++)
  {
    sum = 0.0;
    for (j = 0; j < MAX; j++)
    {
      sum += gsl_matrix_get(ak, i, j)*gsl_vector_get(w, j)*gsl_vector_get(f, j);
    }
    gsl_vector_set(new_result, i, gsl_vector_get(g, i) + sum);
  }

  gsl_matrix_free(ak);
  gsl_vector_free(g);
  return;
}

// K(x,x') = K(x) = 3x
void get_kernel(gsl_matrix *ak, gsl_vector *x, gsl_vector *xp)
{
  int i, j;
  
  for (i = 0; i < MAX; i++)
  {
    for (j = 0; j < MAX; j++)
    {
      gsl_matrix_set(ak, i, j, 3.0*gsl_vector_get(x, i));
    }
  }
  return;
}

// g(x) = -34x - 1
void get_g(gsl_vector *g, gsl_vector *eks)
{
  int i;

  for (i = 0; i < MAX; i++)
  {
    gsl_vector_set(g, i, -34.0*gsl_vector_get(eks, i) - 1.0);
  }
  return;
}
