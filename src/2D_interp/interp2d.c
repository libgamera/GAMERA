/*
 * This file is part of interp2d, a GSL-compatible two-dimensional
 * interpolation library. <http://www.ellipsix.net/interp2d.html>
 * 
 * Copyright 2012 David Zaslavsky
 * Portions based on GNU GSL interpolation code,
 *  copyright 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include "2D_interp/interp2d.h"

/**
 * Triggers a GSL error if the argument is not equal to GSL_SUCCESS.
 * If the argument is GSL_SUCCESS, this does nothing.
 */
#define DISCARD_STATUS(s) if ((s) != GSL_SUCCESS) { GSL_ERROR_VAL("interpolation error", (s),  GSL_NAN); }

interp2d* interp2d_alloc(const interp2d_type* T, size_t xsize, size_t ysize) {
    interp2d* interp;
    if (xsize < T->min_size || ysize < T->min_size) {
        GSL_ERROR_NULL("insufficient number of points for interpolation type", GSL_EINVAL);
    }
    interp = (interp2d*)malloc(sizeof(interp2d));
    if (interp == NULL) {
        GSL_ERROR_NULL("failed to allocate space for interp2d struct", GSL_ENOMEM);
    }
    interp->type = T;
    interp->xsize = xsize;
    interp->ysize = ysize;
    if (interp->type->alloc == NULL) {
        interp->state = NULL;
        return interp;
    }
    interp->state = interp->type->alloc(xsize, ysize);
    if (interp->state == NULL) {
        free(interp);
        GSL_ERROR_NULL("failed to allocate space for interp2d state", GSL_ENOMEM);
    }
    return interp;
}

int interp2d_init(interp2d* interp, const double xarr[], const double yarr[], const double zarr[], size_t xsize, size_t ysize) {
    size_t i;
    if (xsize != interp->xsize || ysize != interp->ysize) {
        GSL_ERROR("data must match size of interpolation object", GSL_EINVAL);
    }
    for (i = 1; i < xsize; i++) {
        if (xarr[i-1] >= xarr[i]) {
            GSL_ERROR("x values must be strictly increasing", GSL_EINVAL);
        }
    }
    for (i = 1; i < ysize; i++) {
        if (yarr[i-1] >= yarr[i]) {
            GSL_ERROR("y values must be strictly increasing", GSL_EINVAL);
        }
    }
    interp->xmin = xarr[0];
    interp->xmax = xarr[xsize - 1];
    interp->ymin = yarr[0];
    interp->ymax = yarr[ysize - 1];
    {
        int status = interp->type->init(interp->state, xarr, yarr, zarr, xsize, ysize);
        return status;
    }
}

void interp2d_free(interp2d* interp) {
    if (!interp) {
        return;
    }
    if (interp->type->free) {
        interp->type->free(interp->state);
    }
    free(interp);
}

/**
 * A wrapper function that checks boundary conditions, calls an evaluator
 * which implements the actual calculation of the function value or 
 * derivative etc., and checks the return status.
 */
static inline int interp2d_eval_impl(
    int (*evaluator)(const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z),
    const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya,
    double* result
) {
    if (x < interp->xmin || x > interp->xmax) {
        char errmsg[80];
        snprintf(errmsg, 80, "x value %g not in range %g to %g", x, interp->xmin, interp->xmax);
        GSL_ERROR(errmsg, GSL_EDOM);
    }
    if (y < interp->ymin || y > interp->ymax) {
        char errmsg[80];
        snprintf(errmsg, 80, "y value %g not in range %g to %g", y, interp->ymin, interp->ymax);
        GSL_ERROR(errmsg, GSL_EDOM);
    }
    return evaluator(interp->state, xarr, yarr, zarr, interp->xsize, interp->ysize, x, y, xa, ya, result);
}

/**
 * Another wrapper function that serves as a drop-in replacement for
 * interp2d_eval_impl but does not check the bounds. This can be used
 * for extrapolation.
 */
static inline int interp2d_eval_impl_no_boundary_check(
    int (*evaluator)(const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z),
    const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya,
    double* result
) {
    return evaluator(interp->state, xarr, yarr, zarr, interp->xsize, interp->ysize, x, y, xa, ya, result);
}

double interp2d_eval(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    double z;
    int status = interp2d_eval_impl(interp->type->eval, interp, xarr, yarr, zarr, x, y, xa, ya, &z);
    DISCARD_STATUS(status)
    return z;
}

double interp2d_eval_no_boundary_check(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    double z;
    int status = interp2d_eval_impl_no_boundary_check(interp->type->eval, interp, xarr, yarr, zarr, x, y, xa, ya, &z);
    DISCARD_STATUS(status)
    return z;
}

int interp2d_eval_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_impl(interp->type->eval, interp, xarr, yarr, zarr, x, y, xa, ya, z);
}

int interp2d_eval_e_no_boundary_check(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_impl_no_boundary_check(interp->type->eval, interp, xarr, yarr, zarr, x, y, xa, ya, z);
}

double interp2d_eval_deriv_x(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    double z;
    int status = interp2d_eval_impl(interp->type->eval_deriv_x, interp, xarr, yarr, zarr, x, y, xa, ya, &z);
    DISCARD_STATUS(status)
    return z;
}

int interp2d_eval_deriv_x_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_impl(interp->type->eval_deriv_x, interp, xarr, yarr, zarr, x, y, xa, ya, z);
}

double interp2d_eval_deriv_y(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    double z;
    int status = interp2d_eval_impl(interp->type->eval_deriv_y, interp, xarr, yarr, zarr, x, y, xa, ya, &z);
    DISCARD_STATUS(status)
    return z;
}

int interp2d_eval_deriv_y_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_impl(interp->type->eval_deriv_y, interp, xarr, yarr, zarr, x, y, xa, ya, z);
}

double interp2d_eval_deriv_xx(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    double z;
    int status = interp2d_eval_impl(interp->type->eval_deriv_xx, interp, xarr, yarr, zarr, x, y, xa, ya, &z);
    DISCARD_STATUS(status)
    return z;
}

int interp2d_eval_deriv_xx_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_impl(interp->type->eval_deriv_xx, interp, xarr, yarr, zarr, x, y, xa, ya, z);
}

double interp2d_eval_deriv_yy(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    double z;
    int status = interp2d_eval_impl(interp->type->eval_deriv_yy, interp, xarr, yarr, zarr, x, y, xa, ya, &z);
    DISCARD_STATUS(status)
    return z;
}

int interp2d_eval_deriv_yy_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_impl(interp->type->eval_deriv_yy, interp, xarr, yarr, zarr, x, y, xa, ya, z);
}

double interp2d_eval_deriv_xy(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    double z;
    int status = interp2d_eval_impl(interp->type->eval_deriv_xy, interp, xarr, yarr, zarr, x, y, xa, ya, &z);
    DISCARD_STATUS(status)
    return z;
}

int interp2d_eval_deriv_xy_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_impl(interp->type->eval_deriv_xy, interp, xarr, yarr, zarr, x, y, xa, ya, z);
}

size_t interp2d_type_min_size(const interp2d_type* T) {
    return T->min_size;
}

size_t interp2d_min_size(const interp2d* interp) {
    return interp->type->min_size;
}

const char* interp2d_name(const interp2d* interp) {
    return interp->type->name;
}
