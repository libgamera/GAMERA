%module gappa
%{
#include "../include/Radiation.h"
#include "../include/Particles.h"
#include "../include/Utils.h"
#include "../include/Astro.h"
#include "../include/2D_interp/interp2d_spline.h"
#include "../include/2D_interp/interp2d.h"
%}

%include "typemaps.i"
%include "std_vector.i"
%include "std_string.i"
%include "std_iostream.i"
namespace std
{
%template(OneIVector) vector<int>;
%template(OneDVector) vector<double>;
%template(TwoDVector) vector< vector<double> >;
}
%include "../include/Radiation.h"
%include "../include/Particles.h"
%include "../include/Utils.h"
%include "../include/Astro.h"
%include "../include/2D_interp/interp2d_spline.h"
%include "../include/2D_interp/interp2d.h"
