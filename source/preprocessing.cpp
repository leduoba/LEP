/*************************************************************************************
Copyright (C) 2015 Qiyang Zhao, Beihang Univ.

This program is free software; you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation; either 
version 2 of the License, or (at your option) any later version.

Please cite the following paper if you want to use the source codes in your work:

Qiyang Zhao, Segmenting natural images with the least effort as humans, BMVC 2015.
**************************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <mex.h>
#include "image.h"
#include "misc.h"
#include "imgio.h"
#include "gmm_rgb.h"
#include "gmm_lbp.h"
#include "rgbLab.h"

// dissimilarity measure between pixels
static inline double colorDist(image<double> *r, image<double> *g, image<double> *b,
						 int x1, int y1, int x2, int y2)
{
	return sqrt(square(imRef(r, x1, y1) -imRef(r, x2, y2)) +
		square(imRef(g, x1, y1) - imRef(g, x2, y2)) +
		square(imRef(b, x1, y1) - imRef(b, x2, y2)));
}

/*
* Calculate sigma
*/
double calcSigma(image<rgb> *im)
{
	int width = im->width(), height = im->height();
	double sigma = 0;

	image<double> *r = new image<double>(width, height);
	image<double> *g = new image<double>(width, height);
	image<double> *b = new image<double>(width, height);

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			imRef(r, x, y) = imRef(im, x, y).r;
			imRef(g, x, y) = imRef(im, x, y).g;
			imRef(b, x, y) = imRef(im, x, y).b;
		}
	}

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (x < width-1) sigma += colorDist(r, g, b, x, y, x + 1, y);
			if (y < height-1) sigma += colorDist(r, g, b, x, y, x, y + 1);
			if ((x < width-1) && (y < height-1)) sigma += colorDist(r, g, b, x, y, x + 1, y + 1);
			if ((x < width-1) && (y > 0)) sigma += colorDist(r, g, b, x, y, x + 1, y - 1);
		}
	}

	delete r;
	delete g;
	delete b;

	return sigma / (height * width * 4 - 3 * height - 3 * width + 2);
}

// matlab调用格式：[rgb, lbp, edges] = features(imgfn, ck, tk)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	srand((unsigned int)(time(0)));

	char* imgfn = (char*)mxArrayToString(prhs[0]);
	int ock = (int)mxGetScalar(prhs[1]);
	int otk = (int)mxGetScalar(prhs[2]);
	double* outData;

	image<rgb> *im = loadimg(imgfn);
	int width = im->width(), height = im->height();
	image<double> * r = new image<double>(width, height);
	image<double> * g = new image<double>(width, height);
	image<double> * b = new image<double>(width, height);

//#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			Lab labcolor = rgb2Lab(imRef(im, x, y).r, imRef(im, x, y).g, imRef(im, x, y).b);

			imRef(im, x, y).r = labcolor.L; imRef(r, x, y) = labcolor.L;
			imRef(im, x, y).g = labcolor.a; imRef(g, x, y) = labcolor.a;
			imRef(im, x, y).b = labcolor.b; imRef(b, x, y) = labcolor.b;
		}
	}

	double sigma = calcSigma(im), sigma_inv = 1.0 / sigma;
	int * rgblabel = (int*)calloc(width * height, sizeof(int));
	int * lbplabel = (int*)calloc(width * height, sizeof(int));
	int ck = calcGMMs_rgb(im, ock, rgblabel);
	int tk = calcGMMs_lbp(im, otk, lbplabel);

	int * edge_a = new int[width * height * 4];
	int * edge_b = new int[width * height * 4];
	double * edge_w = new double[width * height * 4];
	int es1 = 0;
	int es2 = es1 + (width - 1) * height;
	int es3 = es2 + width * (height - 1);
	int es4 = es3 + (width - 1) * (height - 1);
	int ec1 = 0, ec2 = 0, ec3 = 0, ec4 = 0;

	int num = 0;

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (x < width-1) {
				int ap = y * width + x, bp = y * width + (x + 1);
				edge_a[es1 + ec1] = ap;
				edge_b[es1 + ec1] = bp;
				edge_w[es1 + ec1] = colorDist(r, g, b, x, y, x + 1, y) * sigma_inv;
				ec1++;
				num++;
			}

			if (y < height-1) {
				int ap = y * width + x, bp = (y + 1) * width + x;
				edge_a[es2 + ec2] = ap;
				edge_b[es2 + ec2] = bp;
				edge_w[es2 + ec2] = colorDist(r, g, b, x, y, x, y + 1) * sigma_inv;
				ec2++;
				num++;
			}

			if ((x < width-1) && (y < height-1)) {
				int ap = y * width + x, bp = (y + 1) * width + (x + 1);
				edge_a[es3 + ec3] = ap;
				edge_b[es3 + ec3] = bp;
				edge_w[es3 + ec3] = colorDist(r, g, b, x, y, x + 1, y + 1) * sigma_inv;
				ec3++;
				num++;
			}

			if ((x < width-1) && (y > 0)) {
				int ap = y * width + x, bp = (y - 1) * width + (x + 1);
				edge_a[es4 + ec4] = ap;
				edge_b[es4 + ec4] = bp;
				edge_w[es4 + ec4] = colorDist(r, g, b, x, y, x + 1, y - 1) * sigma_inv;
				ec4++;
				num++;
			}
		}
	}

	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); outData = mxGetPr(plhs[0]); outData[0] = height;
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL); outData = mxGetPr(plhs[1]); outData[0] = width;
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL); outData = mxGetPr(plhs[2]); outData[0] = ck;
	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL); outData = mxGetPr(plhs[3]); outData[0] = tk;
	plhs[4] = mxCreateDoubleMatrix(1, num, mxREAL); outData = mxGetPr(plhs[4]);
    for(int i = 0; i < num; i++) outData[i] = edge_a[i]; 
	plhs[5] = mxCreateDoubleMatrix(1, num, mxREAL); outData = mxGetPr(plhs[5]);
    for(int i = 0; i < num; i++) outData[i] = edge_b[i]; 
	plhs[6] = mxCreateDoubleMatrix(1, width * height, mxREAL); outData = mxGetPr(plhs[6]);
    for(int i = 0; i < width * height; i++) outData[i] = rgblabel[i];
	plhs[7] = mxCreateDoubleMatrix(1, width * height, mxREAL); outData = mxGetPr(plhs[7]);
    for(int i = 0; i < width * height; i++) outData[i] = lbplabel[i]; 
	plhs[8] = mxCreateDoubleMatrix(1, num, mxREAL); outData = mxGetPr(plhs[8]);
    for(int i = 0; i < num; i++) outData[i] = edge_w[i]; 

	delete r;
	delete g;
	delete b;
	delete im;
	delete[] edge_a;
	delete[] edge_b;
	delete[] edge_w;
	free(rgblabel);
	free(lbplabel);

	return;
}
