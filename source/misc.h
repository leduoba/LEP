/*************************************************************************************
Copyright (C) 2015 Qiyang Zhao, Beihang Univ.

This program is free software; you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation; either 
version 2 of the License, or (at your option) any later version.

Please cite the following paper if you want to use the source codes in your work:

Qiyang Zhao, Segmenting natural images with the least effort as humans, BMVC 2015.
**************************************************************************************/

#ifndef _MISC_H
#define _MISC_H

#include <cmath>

#ifndef SAVESEGIMG
#define SAVESEGIMG					0
#endif

#ifndef M_RGB_SIGMA2_RELATIVE_THR
#define M_RGB_SIGMA2_RELATIVE_THR	0.05
#endif

#ifndef M_LBP_SIGMA2_RELATIVE_THR
#define M_LBP_SIGMA2_RELATIVE_THR	0.05
#endif

#ifndef M_SIGMA2_RELATIVE_UNCHANGED
#define M_SIGMA2_RELATIVE_UNCHANGED 0.99
#endif

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#ifndef M_SQRT2_HALF
#define M_SQRT2_HALF 0.7071067811865475
#endif

#ifndef M_SAMPLE_SIZE_THR
#define M_SAMPLE_SIZE_THR 400
#endif

typedef unsigned char uchar;

typedef struct { double r, g, b; } rgb;
typedef struct { double L, a, b; } Lab;

inline bool operator==(const rgb &a, const rgb &b) {
  return ((a.r == b.r) && (a.g == b.g) && (a.b == b.b));
}

template <class T>
inline T abs(const T &x) { return (x > 0 ? x : -x); };

template <class T>
inline int sign(const T &x) { return (x >= 0 ? 1 : -1); };

template <class T>
inline T square(const T &x) { return x*x; };

template <class T>
inline T bound(const T &x, const T &min, const T &max) {
  return (x < min ? min : (x > max ? max : x));
}

template <class T>
inline bool check_bound(const T &x, const T&min, const T &max) {
  return ((x < min) || (x > max));
}

inline int vlib_round(float x) { return (int)(x + 0.5F); }

inline int vlib_round(double x) { return (int)(x + 0.5); }

inline double gaussian(double val, double sigma) {
  return exp(-square(val/sigma)/2)/(sqrt(2*M_PI)*sigma);
}

typedef struct {
	double v;
	int x, y;
} projections;

bool operator<(const projections &a, const projections &b) {
	return a.v < b.v;
}

// random color
rgb randomcolors(){
	rgb c = {(uchar)rand(), (uchar)rand(), (uchar)rand()};
	return c;
}

//double angleWeights[16] = {0.0175, 0.1385, 0.2981, 0.4269, 0.4269, 0.4269, 0.4269, 0.4269, 0.4269,
//	0.4269, 0.4269, 0.4269, 1.2429, 1.2429, 1.2429, 0.8134};

#endif
