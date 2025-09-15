/*
 * Solve y(u) = E(a) = H/H0 with u = -ln a (s.t. a = exp(-u)), IC y(0)=1, using GSL.
 * ODE:
 *   dy/du = (3/2) y * Om(a) * [ 1 - (1 - Om(a)) * S(a,B) ],
 *   Om(a) = Omega_m0 * a^{-3} / y^2,
 *   S(a,B) = 1 / ( 1 + (1 + B) * (1 - Om(a)) ).
 *
 * - Arrays are tabulated uniformly in u = -ln a, from u=0 to u_max = -ln(a_min),
 *   then stored in ascending a for spline construction (GSL requires increasing x).
 */

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>


struct gen_gal_params {
  double Om0;   /* Omega_m0 driving the RHS */
  double Or0;   /* Omega_m0 driving the RHS */
  double B;     /* model parameter */
};

/* ODE RHS: dy/du = f(u,y) with u = -ln a */
static int rhs(double u, const double y[], double dydu[], void *params) {
  const struct gen_gal_params *P = (const struct gen_gal_params*) params;
  const double Om   = P->Om0 * exp(3. * u) / (y[0]*y[0]);
  const double Or   = P->Or0 * exp(4. * u) / (y[0]*y[0]);
  const double Oext = Om + Or;
  const double O_smg = 1 - Oext;
  const double S    = 1.0 / (1.0 + (1.0 + P->B) * O_smg);
  const double fac  = 1.0 - O_smg * S;
  dydu[0] = 1.5 * y[0] * Oext * fac;
  return GSL_SUCCESS;
}

/* Simple guard against negative/NaN y during integration */
static int check_state(double y) {
  if (!isfinite(y) || y <= 0.0) return GSL_EFAILED;
  return GSL_SUCCESS;
}

/* Reverse arrays in-place */
static void reverse_double(double *arr, size_t n) {
  for (size_t i = 0; i < n/2; ++i) {
    double tmp = arr[i];
    arr[i] = arr[n-1-i];
    arr[n-1-i] = tmp;
  }
}

/* Build background table and splines */
int gen_gal_build_background(double H0,
                         double Omega_m0,
                         double Omega_r0,
                         double B,
                         gsl_interp_accel *s_acc,
                         gsl_spline *s_rho_smg,
                         gsl_spline *s_p_smg
                         )
{
  const double a_min = 1e-14; /* early-time limit */
  const size_t n_pts = 10000;
  const double rtol  = 1e-9;
  const double atol  = 1e-12;

  if (!(a_min > 0.0 && a_min < 1.0)) return GSL_EINVAL;
  if (n_pts < 5) return GSL_EINVAL;

  int status = GSL_SUCCESS;

  double* a_arr       = malloc(n_pts * sizeof(double));
  double* rho_smg_arr = malloc(n_pts * sizeof(double));
  double* p_smg_arr   = malloc(n_pts * sizeof(double));
  if (!a_arr || !rho_smg_arr || !p_smg_arr) { status = GSL_ENOMEM; goto fail; }

  /* Integration set-up */
  struct gen_gal_params P = { .Om0 = Omega_m0, .Or0 = Omega_r0, .B = B };
  gsl_odeiv2_system sys = { .function = rhs, .jacobian = NULL, .dimension = 1, .params = &P };

  gsl_odeiv2_driver *drv =
    gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-3, atol, rtol);
  if (!drv) { status = GSL_ENOMEM; goto fail; }

  /* Integrate uniformly in u from u=0 to u_max=-ln(a_min).
     We'll store temporary arrays in descending a, then reverse for splines. */
  const double u0 = 0.0;
  const double u1 = -log(a_min);
  const double du = (u1 - u0) / (double)(n_pts - 1); /* note: negative */

  double u = u0;
  double y[1];
  y[0] = 1.0;  /* E(a=1)=1 */

  for (size_t i = 0; i < n_pts; ++i) {
    const double a = exp(-u);
    /* Record state and analytic dEdu (from RHS) for this u,y */
    a_arr[i]    = a;
    const double E = y[0];

    const double Om    = Omega_m0 * pow(a,-3) / (E*E);
    const double Or    = Omega_r0 * pow(a,-4) / (E*E);
    const double Oext = Om + Or;

    rho_smg_arr[i] = 3 * H0*H0*E*E * (1 - Oext);
    if (rho_smg_arr[i] < 0) {
      rho_smg_arr[i] = 0;
    }

    const double S = 1.0 / (1.0 + (1.0 + B) * (1 - Oext));
    p_smg_arr[i] = - rho_smg_arr[i] * (1 + S * Oext);

    if (i + 1 < n_pts) {
      const double u_next = u0 + (double)(i + 1) * du;
      int s = gsl_odeiv2_driver_apply(drv, &u, u_next, y);
      if (s != GSL_SUCCESS) { status = s; break; }
      if ((status = check_state(y[0])) != GSL_SUCCESS) break;
    }
  }

  gsl_odeiv2_driver_free(drv);
  if (status != GSL_SUCCESS) goto fail;

  /* Arrays are descending in a (since x decreased). Reverse to ascending a. */
  reverse_double(a_arr,    n_pts);
  reverse_double(rho_smg_arr, n_pts);
  reverse_double(p_smg_arr, n_pts);

  /* Build splines on 'a' */
  s_acc   = gsl_interp_accel_alloc();
  s_rho_smg = gsl_spline_alloc(gsl_interp_cspline, n_pts);
  s_p_smg = gsl_spline_alloc(gsl_interp_cspline, n_pts);
  if (!s_acc || s_rho_smg || s_p_smg) { status = GSL_ENOMEM; goto fail; }

  if ((status = gsl_spline_init(s_rho_smg, a_arr, rho_smg_arr, n_pts)) != GSL_SUCCESS) goto fail;
  if ((status = gsl_spline_init(s_p_smg, a_arr, p_smg_arr, n_pts)) != GSL_SUCCESS) goto fail;

  free(a_arr);
  free(rho_smg_arr);
  free(p_smg_arr);

  return GSL_SUCCESS;

fail:
  if (s_rho_smg) gsl_spline_free(s_rho_smg);
  if (s_p_smg) gsl_spline_free(s_p_smg);
  if (s_acc)    gsl_interp_accel_free(s_acc);
  return status;
}
