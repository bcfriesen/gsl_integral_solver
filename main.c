#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

// size of vectors/matrices
#define MAX 100

int main (void)
{
  int a = 0, b = 4;
  int i;
  gsl_vector *tmp = gsl_vector_alloc(MAX);
  gsl_vector *grid = gsl_vector_alloc(MAX);

  for (i = 0; i < MAX; i++)
  {
    gsl_vector_set(grid, i, (double)i/10.0);
  }

  get_g(tmp, grid);
  
  for (i = 0; i < MAX; i++)
  {
    printf("%.18f %.18f\n", gsl_vector_get(grid, i), gsl_vector_get(tmp, i));
  }

  gsl_vector_free(grid);
  gsl_vector_free(tmp);
  return 0;
}



void fred2_solver(double a, double b, gsl_vector *t, gsl_vector *w, gsl_vector *g, gsl_matrix *ak)
{
  gsl_integration_glfixed_table *glgrid = gsl_integration_glfixed_table_alloc(MAX);
  gsl_matrix *lhs = gsl_matrix_alloc(MAX, MAX);
  gsl_matrix *ktilde = gsl_matrix_alloc(MAX, MAX);
  int i, j, error;
  double ptsi, wghtsi;

  // set Gauss-Legendre integration points and weights
  for (i = 0; i < MAX; i++)
  {
    error = gsl_integration_glfixed_point(a, b, i, &ptsi, &wghtsi, glgrid);
    gsl_vector_set(t, i, ptsi);
    gsl_vector_set(w, i, wghtsi);
  }
  // fill in unit matrix first
  gsl_matrix_set_identity(lhs);

  for (i = 0; i < MAX; i++)
  {
    for (j = 0; j < MAX; j++)
    {
      gsl_matrix_set(ktilde, i, j, gsl_matrix_get(ak, i, j)*gsl_vector_get(w, j));
    }
  }

  gsl_integration_glfixed_table_free(glgrid);
  gsl_matrix_free(lhs);
  return;
}




void get_kernel(gsl_matrix *ak, gsl_vector *eks, gsl_vector *eksp)
{
  int i, j;
  
  for (i = 0; i < MAX; i++)
  {
    for (j = 0; j < MAX; j++)
    {
      gsl_matrix_set(ak, i, j, 3.0*gsl_vector_get(eks, i));
    }
  }
  return;
}



void get_g(gsl_vector *g, gsl_vector *eks)
{
  int i;

  for (i = 0; i < MAX; i++)
  {
    gsl_vector_set(g, i, -34.0*gsl_vector_get(eks, i) - 1.0);
  }
  return;
}
