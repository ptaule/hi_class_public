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
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include "solve_gen_galileon.h"

struct gen_gal_params {
  double Om0;   /* Omega_m0 driving the RHS */
  double Or0;   /* Omega_m0 driving the RHS */
  double B;     /* model parameter */
};

/* More stable version ? */
/* ODE RHS: dy/du = f(u,y) with u = -ln a, y[0] = E(a) */
static int rhs(double u, const double y[], double dydu[], void *params) {
  const struct gen_gal_params *P = (const struct gen_gal_params*) params;
  const double E = y[0];
  if (!isfinite(E) || E <= 0.0) return GSL_EFAILED;

  const double invE2 = 1.0 / (E*E);
  const double Om = P->Om0 * exp(3.0 * u) * invE2;   /* Ωm(a) */
  const double Or = P->Or0 * exp(4.0 * u) * invE2;   /* Ωr(a) */
  const double Oext = Om + Or;

  double num = 1.0 + P->B - P->B * Oext;
  double den = 2.0 + P->B - (1.0 + P->B) * Oext;

  /* Guard the RD/MD limit and pathological roundoff */
  double fac;
  if (!isfinite(num) || !isfinite(den) || den <= 0.0) {
    fac = 1.0;
  } else {
    fac = num / den;
    /* In RD/MD, Oext → 1 ⇒ fac → 1. Snap very close to 1 to kill noise. */
    if (Oext > 1.0 - 1e-14) fac = 1.0;
    /* Keep fac in a sane range */
    if (fac > 1.0) fac = 1.0;
    if (fac < 0.0) fac = 0.0;
  }

  dydu[0] = 1.5 * E * (Om + (4.0/3.0) * Or) * fac;
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
int gen_gal_build_background(struct background * pba)
{
  const double H0       = pba->H0;
  const double Omega_m0 = pba->Omega0_b + pba->Omega0_cdm;
  const double Omega_r0 = pba->Omega0_g + pba->Omega0_ur;
  const double B        = pba->parameters_smg[0];

  const double a_min = 1e-14; /* early-time limit */
  const size_t n_pts = 10000;
  const double rtol  = 1e-9;
  const double atol  = 1e-10;
  const double Omega_smg_match_th = 1e-4;  /* matching threshold */

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
  const double du = (u1 - u0) / (double)(n_pts - 1);

  double u = u0;
  double y[1];
  y[0] = 1.0;  /* E(a=1)=1 */

  /* Asymptotic matching state, using rho_smg ~ 1/H^2 */
  int use_asymp = 0;
  double C4 = 0.0; /* C4 = Omega_smg(a_match) * E(a_match)^4 */

  for (size_t i = 0; i < n_pts; ++i) {
    const double a = exp(-u);
    /* Record state and analytic dEdu (from RHS) for this u,y */
    a_arr[i]    = a;
    const double E = y[0];

    double dEdu_over_E;
    {
      double dy_temp[1];
      rhs(u,&E,dy_temp,&P);
      dEdu_over_E = dy_temp[0] / E;
    }

    const double Om    = Omega_m0 * pow(a,-3) / (E*E);
    const double Or    = Omega_r0 * pow(a,-4) / (E*E);

    const double denom = 1.5 * (Om + 4.0/3.0 * Or);
    double fac = dEdu_over_E / denom;

    const double Oext = Om + Or;
    if (Oext > 1.0 - 1e-14) fac = 1.0;
    if (fac < 0.0) fac = 0.0;
    if (fac > 1.0) fac = 1.0;

    double Omega_smg = (1.0 - fac) / (fac * (1.0 + B) - B);
    /* Clamp true negatives: */
    if (Oext > 1.0 - 1e-12 && Omega_smg < 0.0) Omega_smg = 0.0;

    /* Early-time asymptotic match: Ω_smg(a) = C4 / E(a)^4 once threshold crossed */
    if (!use_asymp && Omega_smg < Omega_smg_match_th) {
      use_asymp = 1;
      C4 = Omega_smg * E*E*E*E;  /* freeze C4 at the match point */
    }

    Omega_smg = use_asymp ? (C4 / (E*E*E*E)) : Omega_smg;

    rho_smg_arr[i] = 3.0 * H0*H0 * E*E * Omega_smg;

    const double S = 1.0 / (1.0 + (1.0 + B) * Omega_smg);
    p_smg_arr[i] = - rho_smg_arr[i] * (1.0 + S * (Om + 4.0/3.0 * Or));

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
  pba->s_acc   = gsl_interp_accel_alloc();
  pba->s_rho_smg = gsl_spline_alloc(gsl_interp_steffen, n_pts);
  pba->s_p_smg = gsl_spline_alloc(gsl_interp_steffen, n_pts);
  if (!pba->s_acc || !pba->s_rho_smg || !pba->s_p_smg) { status = GSL_ENOMEM; goto fail; }

  if ((status = gsl_spline_init(pba->s_rho_smg, a_arr, rho_smg_arr, n_pts)) != GSL_SUCCESS) goto fail;
  if ((status = gsl_spline_init(pba->s_p_smg, a_arr, p_smg_arr, n_pts)) != GSL_SUCCESS) goto fail;

  free(a_arr);
  free(rho_smg_arr);
  free(p_smg_arr);

  return GSL_SUCCESS;

fail:
  if (pba->s_rho_smg) gsl_spline_free(pba->s_rho_smg);
  if (pba->s_p_smg) gsl_spline_free(pba->s_p_smg);
  if (pba->s_acc)    gsl_interp_accel_free(pba->s_acc);
  class_test(_FALSE_,
             pba->error_message,
             "Something went wrong in gen_gal_build_background: %d", status
             );
  return status;
}
