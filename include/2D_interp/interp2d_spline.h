/*
 * This file is part of interp2d, a GSL-compatible two-dimensional
 * interpolation library. <http://www.ellipsix.net/interp2d.html>
 * 
 * Copyright 2013 David Zaslavsky
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
/** 
 * @file interp2d_spline.h
 * @brief Public interface for high-level 2D interpolation functions
 * 
 * This is the public interface to the high-level functions of the
 * `%interp2d` library, the ones that _do_ store the data arrays.
 * If you're using the high-level interface, you only use the functions in
 * this file, unless you're creating a new interpolation type.
 * 
 * The typical workflow is
 * 
 * 1. create an interpolation spline object using interp2d_spline_alloc()
 * 2. initialize it using interp2d_spline_init()
 * 3. evaluate the interpolating function or its derivatives using
 *    interp2d_spline_eval() or its counterparts, possibly many times
 * 4. free the memory using interp2d_spline_free()
 * 
 * @see interp2d.h
 */
#ifndef __INTERP_2D_SPLINE_H__
#define __INTERP_2D_SPLINE_H__

#include <gsl/gsl_interp.h>
#include "2D_interp/interp2d.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * A 2D interpolation object which stores the arrays defining the function.
 * In all other respects, this is just like an interp2d object.
 * Instances should be obtained from interp2d_spline_alloc().
 */
typedef struct {
    /** The low-level interpolation object */
    interp2d interp_object;
    /** The x data array */
    double* xarr;
    /** The y data array */
    double* yarr;
    /** The z data array */
    double* zarr;
} interp2d_spline;

/**
 * Allocate a new interpolation spline object. When you want to do an
 * interpolation using the high-level interface, you start by calling
 * this function. Don't forget to free the memory using
 * interp2d_spline_free() when you're done with it.
 * 
 * @param[in] T the `const struct` representing the interpolation algorithm
 *  you want to use. This should be one of the predefined types provided
 *  in this library, like interp2d_bilinear or interp2d_bicubic, unless you
 *  have your own interpolation type you want to use.
 * @param[in] xsize the size in the x direction of the grid that will
 *  specify the function being interpolated
 * @param[in] ysize the size in the y direction of the grid that will
 *  specify the function being interpolated
 * @return the newly allocated interpolation object
 */
interp2d_spline* interp2d_spline_alloc(const interp2d_type* T, size_t xsize, size_t ysize);

/**
 * Initialize an interpolation spline object with data points that define
 * the function being interpolated. This method stores the sizes of the
 * arrays and possibly other relevant data in the interp2d_spline object,
 * but unlike the low-level equivalent interp2d_init(), it also stores
 * the data arrays `xa`, `ya`, and `za`. The content of the arrays is copied,
 * so if you change the array elements after calling this, it won't affect
 * the interpolation spline object.
 * 
 * This completely resets the state of the object, so it is safe to reuse an
 * existing object to interpolate a new function by simply calling
 * interp2d_spline_init() on the object with the new arrays. The new arrays
 * must be the same size as the old arrays in both dimensions, otherwise this
 * function returns an error code.
 * 
 * @param[in] interp the interpolation object, previously initialized
 * @param[in] xa the x coordinates of the data, of length `xsize`
 * @param[in] ya the y coordinates of the data, of length `ysize`
 * @param[in] za the z coordinates of the data, of length `xsize*ysize`
 * @param[in] xsize the length of the array `xa`
 * @param[in] ysize the length of the array `ya`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_spline_init(interp2d_spline* interp, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize);

/**
 * Free the interpolation spline object.
 * 
 * @param[in] interp an interpolation spline object previously allocated
 *  with interp2d_spline_alloc()
 */
void interp2d_spline_free(interp2d_spline* interp);

/**
 * Evaluate the interpolating function at the point `(x,y)`.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the value of the interpolating function at `(x,y)`
 */
double interp2d_spline_eval(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the interpolating function at the point `(x,y)`.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the value of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_spline_eval_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the x derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the x derivative of the interpolating function at `(x,y)`
 */
double interp2d_spline_eval_deriv_x(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the x derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the x derivative of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_spline_eval_deriv_x_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the y derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the y derivative of the interpolating function at `(x,y)`
 */
double interp2d_spline_eval_deriv_y(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the y derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the y derivative of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_spline_eval_deriv_y_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the second x derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the second x derivative of the interpolating function at `(x,y)`
 */
double interp2d_spline_eval_deriv_xx(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the second x derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the second x derivative of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_spline_eval_deriv_xx_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the second y derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the second y derivative of the interpolating function at `(x,y)`
 */
double interp2d_spline_eval_deriv_yy(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the second y derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the second y derivative of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_spline_eval_deriv_yy_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the cross derivative of the interpolating function at `(x,y)`.
 * This is \f$\partial_{xy}f(x,y)\f$.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the cross derivative of the interpolating function at `(x,y)`
 */
double interp2d_spline_eval_deriv_xy(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the cross derivative of the interpolating function at `(x,y)`.
 * This is \f$\partial_{xy}f(x,y)\f$.
 * 
 * @param[in] interp the interpolation spline object, already initialized
 *  using interp2d_spline_init()
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the cross derivative of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_spline_eval_deriv_xy_e(const interp2d_spline* interp, const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Return the minimum number of points in each dimension needed by the
 * type of the given interpolation object. This just accesses the type
 * from `interp` and calls interp2d_type_min_size() on it.
 */
size_t interp2d_spline_min_size(const interp2d_spline* interp);

/**
 * Return the type name of the given interpolation. This just accesses
 * the type object from `interp` and returns its name.
 */
const char* interp2d_spline_name(const interp2d_spline* interp);

#ifdef __cplusplus
}
#endif

#endif
