////////////////////////////////////////////////////////////////////////////////
// spline_interpolation.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Code for C2 cubic spline interpolation.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  11/18/2022 11:07:41
*///////////////////////////////////////////////////////////////////////////////
#ifndef SPLINE_INTERPOLATION_HH
#define SPLINE_INTERPOLATION_HH
#include "spline.hh"

Spline naturalCubicSplineInterpolant(std::vector<Eigen::Vector3f> &dataPts, Spline::ParametrizationType paramType);

#endif /* end of include guard: SPLINE_INTERPOLATION_HH */
