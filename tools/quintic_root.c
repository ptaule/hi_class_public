#include <math.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>

#include "common.h"

int rf_solve_quintic(const double coeffs[6], double roots[5][2], ErrorMsg errmsg) {

  class_test(coeffs == NULL, errmsg, "Input coefficient array is NULL.");

  const double a5 = coeffs[5];

  class_test(fabs(a5) < _EPSILON_, errmsg, "Leading coefficient (a5) must not be zero.");

  // Normalize to monic polynomial
  const double b0 = coeffs[0] / a5;
  const double b1 = coeffs[1] / a5;
  const double b2 = coeffs[2] / a5;
  const double b3 = coeffs[3] / a5;
  const double b4 = coeffs[4] / a5;

  gsl_matrix *companion = gsl_matrix_alloc(5, 5);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(5, 5);
  gsl_vector_complex *eval = gsl_vector_complex_alloc(5);
  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(5);

  class_test((companion == NULL) || (evec == NULL) || (eval == NULL) || (w == NULL),
             errmsg, "Memory allocation failed in rf_solve_quintic.");

  gsl_matrix_set_zero(companion);
  gsl_matrix_set(companion, 0, 0, -b4);
  gsl_matrix_set(companion, 0, 1, -b3);
  gsl_matrix_set(companion, 0, 2, -b2);
  gsl_matrix_set(companion, 0, 3, -b1);
  gsl_matrix_set(companion, 0, 4, -b0);

  for (int i = 1; i < 5; ++i) {
    gsl_matrix_set(companion, i, i - 1, 1.0);
  }

  int status = gsl_eigen_nonsymmv(companion, eval, evec, w);
  class_test(status != GSL_SUCCESS, errmsg,
             "GSL eigenvalue computation failed with code %d", status);

  // Extract roots
  for (int i = 0; i < 5; ++i) {
    gsl_complex z = gsl_vector_complex_get(eval, i);
    roots[i][0] = GSL_REAL(z);
    roots[i][1] = GSL_IMAG(z);
  }

  gsl_matrix_free(companion);
  gsl_matrix_complex_free(evec);
  gsl_vector_complex_free(eval);
  gsl_eigen_nonsymmv_free(w);

  return _SUCCESS_;
}

void rf_evaluate_poly_complex(
    const double coeffs[6],
    double z_real,
    double z_imag,
    double result[2]
) {
  gsl_complex z = gsl_complex_rect(z_real, z_imag);
  gsl_complex z_pow = gsl_complex_rect(1.0, 0.0);  // z^0 = 1
  gsl_complex sum = gsl_complex_rect(0.0, 0.0);

  for (int i = 0; i < 6; ++i) {
    gsl_complex term = gsl_complex_mul_real(z_pow, coeffs[i]);
    sum = gsl_complex_add(sum, term);
    z_pow = gsl_complex_mul(z_pow, z);
  }

  result[0] = GSL_REAL(sum);
  result[1] = GSL_IMAG(sum);
}
