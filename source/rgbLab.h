/*************************************************************************************
Copyright (C) 2015 Qiyang Zhao, Beihang Univ.

This program is free software; you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation; either 
version 2 of the License, or (at your option) any later version.

Please cite the following paper if you want to use the source codes in your work:

Qiyang Zhao, Segmenting natural images with the least effort as humans, BMVC 2015.
**************************************************************************************/

#ifndef _RGBLAB
#define _RGBLAB

#include "math.h"
#include "misc.h"
using namespace std;

Lab rgb2Lab(double r, double g, double b)
{
	// designed for r, g, b in [0, 1], otherwise divided by 255 is necessary
	double thr = 0.008856;
	
	// RGB to XYZ
	double x = 0.412453 * r + 0.357580 * g + 0.180423 * b;
	double y = 0.212671 * r + 0.715160 * g + 0.072169 * b;
	double z = 0.019334 * r + 0.119193 * g + 0.950227 * b;


	// Normalize for D65 white point
	x = x / 0.950456;
	z = z / 1.088754;

	double y3 = pow(y, 1.0 / 3);
	
	double fx = x > thr ? pow(x, 1.0 / 3) : (7.787 * x + 16.0 / 116);
	double fy = y > thr ? y3 : (7.787 * y + 16.0 / 116);
	double fz = z > thr ? pow(z, 1.0 / 3) : (7.787 * z + 16.0 / 116);

	Lab cl;

	cl.L = y > thr ? (116 * y3 - 16.0) : (903.3 * y);
	cl.a = 500 * (fx - fy);
	cl.b = 200 * (fy - fz);

	return cl;
}

rgb Lab2rgb(double L, double a, double b)
{
	// Thresholds
	double thr1 = 0.008856;
	double thr2 = 0.206893;

	// Compute Y
	double fy0 = pow((double(L + 16) / 116), 3);
	double fy = (fy0 > thr1) ? fy0 : (L / 903.3);
	double y = fy;

	// Alter fY slightly for further calculations
	fy = (fy0 > thr1) ? pow(fy, 1.0 / 3) : (7.787 * fy + 16.0 / 116);

	// Compute X
	double fx = a / 500 + fy;
	double x = (fx > thr2) ? (fx * fx * fx) : ((fx - 16.0 / 116) / 7.787);

	// Compute Z
	double fz = fy - b / 200;
	double z = (fz > thr2) ? (fz * fz * fz) : ((fz - 16.0 / 116) / 7.787);

	x = x * 0.950456;
	z = z * 1.088754;

	rgb color;
	color.r = max(min(3.240479 * x - 1.537150 * y - 0.498535 * z, 1.0), 0.0) * 255;
	color.g = max(min(-0.969256 * x + 1.875992 * y + 0.041556 * z, 1.0), 0.0) * 255;
	color.b = max(min(0.055648 * x - 0.204043 * y + 1.057311 * z, 1.0), 0.0) * 255;

	return color;
}

#endif
