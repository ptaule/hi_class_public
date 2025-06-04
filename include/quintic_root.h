#ifndef __QUINTIC_ROOT_H__
#define __QUINTIC_ROOT_H__

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

int rf_solve_quintic(
  const double coeffs[6],
  double roots[5][2],
  ErrorMsg errmsg);
void rf_evaluate_poly_complex(
  const double coeffs[6],
  double z_real,
  double z_imag,
  double result[2]
);

#ifdef __cplusplus
}
#endif

#endif
