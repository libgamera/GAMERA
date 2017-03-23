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
/** 
 * @file interp2d.h
 * @brief Public interface for low-level 2D interpolation functions
 * 
 * This is the public interface to the low-level functions of the
 * `%interp2d` library, i.e. the functions that don't store the data arrays.
 * If you're using the low-level interface, you only use the functions in
 * this file, unless you're creating a new interpolation type.
 * 
 * The typical workflow is
 * 
 * 1. create an interpolation object using interp2d_alloc()
 * 2. initialize it using interp2d_init()
 * 3. evaluate the interpolating function or its derivatives using
 *    interp2d_eval() or its counterparts, possibly many times
 * 4. free the memory using interp2d_free()
 * 
 * @see interp2d_spline.h
 */

#ifndef __INTERP_2D_H__
#define __INTERP_2D_H__

#include <gsl/gsl_interp.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * An algorithm for 2D interpolation.
 * 
 * There will be one instance of this `struct` representing each type
 * of interpolation that the library supports: bilinear, bicubic, etc.
 * You can use a custom interpolation type by declaring an instance of
 * this and populating it with pointers to the functions that implement
 * your interpolation, along with the other fields.
 * 
 * @see interp2d_bilinear
 *      interp2d_bicubic
 */
typedef struct {
    /** The name of the algorithm. Use interp2d_name() to access this field. */
    const char* name;
    /** The minimum number of points in each dimension required by the algorithm. Use interp2d_min_size() to access this field. */
    unsigned int min_size;
    /** The method that allocates memory for an interpolation object of this type. */
    void* (*alloc)(size_t xsize, size_t ysize);
    /** The method that initializes the interpolation type. */ 
    int (*init)(void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize);
    /** The method that evaluates the interpolation at the given point. */
    int (*eval)(const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z);
    /** The method that evaluates the x derivative of the interpolation at the given point. */
    int (*eval_deriv_x) (const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z_p);
    /** The method that evaluates the y derivative of the interpolation at the given point. */
    int (*eval_deriv_y) (const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z_p);
    /** The method that evaluates the second x derivative of the interpolation at the given point. */
    int (*eval_deriv_xx) (const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z_pp);
    /** The method that evaluates the second y derivative of the interpolation at the given point. */
    int (*eval_deriv_xy) (const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z_pp);
    /** The method that evaluates the cross derivative of the interpolation at the given point. */
    int (*eval_deriv_yy) (const void*, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize, double x, double y, gsl_interp_accel*, gsl_interp_accel*, double* z_pp);
    /** The method that frees the memory. */
    void (*free)(void*);
} interp2d_type;

/**
 * Represents a 2D interpolation object. This is the thing you create
 * to actually do an interpolation. Instances of this should be obtained
 * from interp2d_alloc(). (Even if you are using a custom interpolation
 * type.)
 */
typedef struct {
    /** The type object for the interpolation. */
    const interp2d_type* type;
    /** The minimum value of x for which data have been provided. */
    double xmin;
    /** The maximum value of x for which data have been provided. */
    double xmax;
    /** The minimum value of y for which data have been provided. */
    double ymin;
    /** The maximum value of y for which data have been provided. */
    double ymax;
    /** The number of x values provided. */
    size_t xsize;
    /** The number of y values provided. */
    size_t ysize;
    /** A state object. This is specific to the interpolation type. */
    void* state;
} interp2d;

/** The interpolation type for bilinear interpolation. */
extern const interp2d_type* interp2d_bilinear;
/** The interpolation type for bicubic interpolation. */
extern const interp2d_type* interp2d_bicubic;
// To be added:
// GSL_VAR const interp2d_type* interp2d_bipoly;

/**
 * Return the minimum number of points in each dimension needed by 
 * an interpolation type. For example bilinear interpolation needs
 * at least a 2 by 2 grid, so `interp2d_type_min_size(interp2d_bilinear)`
 * returns 2.
 */
size_t interp2d_type_min_size(const interp2d_type* T);
/**
 * Return the minimum number of points in each dimension needed by the
 * type of the given interpolation object. This just accesses the type
 * from `interp` and calls interp2d_type_min_size() on it.
 */
size_t interp2d_min_size(const interp2d* interp);
/**
 * Return the type name of the given interpolation. This just accesses
 * the type object from `interp` and returns its name.
 */
const char* interp2d_name(const interp2d* interp);

/**
 * Allocate a new interpolation object. When you want to do an interpolation,
 * you start by calling this function. Don't forget to free the memory using
 * interp2d_free() when you're done with it.
 * 
 * Implementation note: this just performs some checks and calls the allocator
 * method specified in the interp2d_type instance. So you can use it with
 * custom interpolation types.
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
interp2d* interp2d_alloc(const interp2d_type* T, size_t xsize, size_t ysize);

/**
 * Initialize an interpolation object with data points that define the
 * function being interpolated. This method stores the sizes of the arrays
 * and possibly other relevant data in the interp2d object, but it does not
 * actually store the arrays `xa`, `ya`, and `za`. So you need to pass the
 * same arrays to interp2d_eval() or whatever evaluation function you use.
 * 
 * This completely resets the state of the interpolation object, so it is safe
 * to reuse an existing interpolation object to interpolate a new function by
 * simply calling interp2d_init() on the interpolation object with the arrays
 * defining the new function.
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
int interp2d_init(interp2d* interp, const double xa[], const double ya[], const double za[], size_t xsize, size_t ysize);

/**
 * Free the interpolation object.
 * 
 * @param[in] interp an interpolation object previously allocated with
 *  interp2d_alloc()
 */
void interp2d_free(interp2d* interp);

/**
 * Evaluate the interpolating function at the point `(x,y)`.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the value of the interpolating function at `(x,y)`
 */
double interp2d_eval(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the interpolating function at the point `(x,y)`, without checking
 * whether `(x,y)` is within the bounds of the interpolation. This function
 * can be used for extrapolation - but do so at your own risk! Extrapolation
 * quickly loses its predictive value as you move away from the region in which
 * function values are known.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the value of the interpolating function at `(x,y)`
 */
double interp2d_eval_no_boundary_check(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the interpolating function at the point `(x,y)`.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the value of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_eval_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the interpolating function at the point `(x,y)`, without checking
 * whether `(x,y)` is within the bounds of the interpolation. This function
 * can be used for extrapolation - but do so at your own risk! Extrapolation
 * quickly loses its predictive value as you move away from the region in which
 * function values are known.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the value of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_eval_e_no_boundary_check(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the x derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the x derivative of the interpolating function at `(x,y)`
 */
double interp2d_eval_deriv_x(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the x derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the x derivative of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_eval_deriv_x_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the y derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the y derivative of the interpolating function at `(x,y)`
 */
double interp2d_eval_deriv_y(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the y derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the y derivative of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_eval_deriv_y_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the second x derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the second x derivative of the interpolating function at `(x,y)`
 */
double interp2d_eval_deriv_xx(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the second x derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the second x derivative of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_eval_deriv_xx_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the second y derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the second y derivative of the interpolating function at `(x,y)`
 */
double interp2d_eval_deriv_yy(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the second y derivative of the interpolating function at `(x,y)`.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the second y derivative of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_eval_deriv_yy_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Evaluate the cross derivative of the interpolating function at `(x,y)`.
 * This is \f$\partial_{xy}f(x,y)\f$.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @return the cross derivative of the interpolating function at `(x,y)`
 */
double interp2d_eval_deriv_xy(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya);

/**
 * Evaluate the cross derivative of the interpolating function at `(x,y)`.
 * This is \f$\partial_{xy}f(x,y)\f$.
 * 
 * @param[in] interp the interpolation object, which should have already
 *  been initialized using the arrays `xarr`, `yarr`, and `zarr`. This
 *  object stores the sizes of the arrays (among other values).
 * @param[in] xarr the x coordinates of the points defining the function.
 *  This should be the same x array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] yarr the y coordinates of the points defining the function
 *  This should be the same y array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] zarr the z coordinates of the points defining the function
 *  This should be the same z array that was passed to interp2d_init()
 *  when initializing `interp`.
 * @param[in] x the x coordinate at which to evaluate the interpolation
 * @param[in] y the y coordinate at which to evaluate the interpolation
 * @param[in] xa the accelerator object for the x direction (may be NULL)
 * @param[in] ya the accelerator object for the y direction (may be NULL)
 * @param[out] z the cross derivative of the interpolating function at `(x,y)`
 * @return a status code, either `GSL_SUCCESS` or an error code indicating
 *  what went wrong
 */
int interp2d_eval_deriv_xy_e(const interp2d* interp, const double xarr[], const double yarr[], const double zarr[], const double x, const double y, gsl_interp_accel* xa, gsl_interp_accel* ya, double* z);

/**
 * Compute the index into a 1D array for a given pair of `x,y` indices.
 * This implements row-major ordering.
 * 
 * In _theory_, if you are working with column-major arrays, you should
 * be able to redefine this macro as
 * 
 *     #define INDEX_2D(xi, yi, xsize, ysize) ((xi) * (ysize) + (yi))
 * 
 * and all `%interp2d` functions would then use column-major ordering.
 * But this is not guaranteed to work properly. It's probably better to
 * just transpose your arrays.
 * 
 * @param[in] xi the index of the desired grid position along the x dimension
 * @param[in] yi the index of the desired grid position along the y dimension
 * @param[in] xsize the size of the grid in the x dimension
 * @param[in] ysize the size of the grid in the y dimension
 */
#define INDEX_2D(xi, yi, xsize, ysize) ((yi) * (xsize) + (xi))

#ifdef __cplusplus
}
#endif

#endif
