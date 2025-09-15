#ifndef __SOLVE_GEN_GALILEON__
#define __SOLVE_GEN_GALILEON__

#include "common.h"
#include "background.h"

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int gen_gal_build_background(struct background * pba);

#ifdef __cplusplus
}
#endif

#endif
