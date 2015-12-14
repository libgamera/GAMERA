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

#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include "2D_interp/interp2d_spline.h"

interp2d_spline* interp2d_spline_alloc(const interp2d_type* T, size_t xsize, size_t ysize) {
    double* array_mem;
    interp2d_spline* interp;
    if (xsize < T->min_size || ysize < T->min_size) {
        GSL_ERROR_NULL("insufficient number of points for interpolation type", GSL_EINVAL);
    }
    interp = (interp2d_spline*)malloc(sizeof(interp2d_spline));
    if (interp == NULL) {
        GSL_ERROR_NULL("failed to allocate space for interp2d_spline struct", GSL_ENOMEM);
    }
    interp->interp_object.type = T;
    interp->interp_object.xsize = xsize;
    interp->interp_object.ysize = ysize;
    if (interp->interp_object.type->alloc == NULL) {
        interp->interp_object.state = NULL;
    }
    else {
        interp->interp_object.state = interp->interp_object.type->alloc(xsize, ysize);
        if (interp->interp_object.state == NULL) {
            free(interp);
            GSL_ERROR_NULL("failed to allocate space for interp2d_spline state", GSL_ENOMEM);
        }
    }
    // Use one contiguous block of memory for all three data arrays.
    // That way the code fails immediately if there isn't sufficient space for everything,
    // rather than allocating one or two and then having to free them.
    array_mem = (double*)calloc(xsize + ysize + xsize * ysize, sizeof(double));
    if (array_mem == NULL) {
        interp->interp_object.type->free(interp->interp_object.state);
        free(interp);
        GSL_ERROR_NULL("failed to allocate space for data arrays", GSL_ENOMEM);
    }
    interp->xarr = array_mem;
    interp->yarr = array_mem + xsize;
    interp->zarr = array_mem + xsize + ysize;
    return interp;
}

int interp2d_spline_init(interp2d_spline* interp, const double xarr[], const double yarr[], const double zarr[], size_t xsize, size_t ysize) {
    int status = interp2d_init(&(interp->interp_object), xarr, yarr, zarr, xsize, ysize);
    memcpy(interp->xarr, xarr, xsize * sizeof(double));
    memcpy(interp->yarr, yarr, ysize * sizeof(double));
    memcpy(interp->zarr, zarr, xsize * ysize * sizeof(double));
    return status;
}

void interp2d_spline_free(interp2d_spline* interp) {
    if (!interp) {
        return;
    }
    if (interp->interp_object.type->free) {
        interp->interp_object.type->free(interp->interp_object.state);
    }
    // interp->xarr points to the beginning of one contiguous block of memory
    // that holds interp->xarr, interp->yarr, and interp->zarr. So it all gets
    // freed with one call. cf. interp2d_spline_alloc() implementation
    free(interp->xarr);
    free(interp);
}

// All these methods just call their interp2d_eval* equivalents

double interp2d_spline_eval(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    return interp2d_eval(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya);
}

int interp2d_spline_eval_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_e(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya, z);
}

double interp2d_spline_eval_deriv_x(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    return interp2d_eval_deriv_x(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya);
}

int interp2d_spline_eval_deriv_x_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_deriv_x_e(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya, z);
}

double interp2d_spline_eval_deriv_y(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    return interp2d_eval_deriv_y(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya);
}

int interp2d_spline_eval_deriv_y_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_deriv_y_e(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya, z);
}

double interp2d_spline_eval_deriv_xx(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    return interp2d_eval_deriv_xx(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya);
}

int interp2d_spline_eval_deriv_xx_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_deriv_xx_e(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya, z);
}

double interp2d_spline_eval_deriv_yy(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    return interp2d_eval_deriv_yy(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya);
}

int interp2d_spline_eval_deriv_yy_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_deriv_yy_e(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya, z);
}

double interp2d_spline_eval_deriv_xy(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya) {
    return interp2d_eval_deriv_xy(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya);
}

int interp2d_spline_eval_deriv_xy_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z) {
    return interp2d_eval_deriv_xy_e(&(interp->interp_object), interp->xarr, interp->yarr, interp->zarr, x, y, xa, ya, z);
}

size_t interp2d_spline_min_size(const interp2d_spline* interp) {
    return interp2d_min_size(&(interp->interp_object));
}

const char* interp2d_spline_name(const interp2d_spline* interp) {
    return interp2d_name(&(interp->interp_object));
}
