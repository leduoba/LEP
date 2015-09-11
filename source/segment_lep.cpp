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
#include <mex.h>
#include "image.h"
#include "misc.h"
#include "monotonic_merge.h"
#include "imgio.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	int amount_v, amount_e, ck, ock, tk, otk, height, width, thrn, *edge_a, *edge_b, *rgblabel, *lbplabel;
	double *edge_w, *thrs, *angleWeights, *ep_prior;
	
	amount_v = (int)mxGetScalar(prhs[0]);
	amount_e = (int)mxGetScalar(prhs[1]);
	edge_a = (int*) mxGetData(prhs[2]);
	edge_b = (int*) mxGetData(prhs[3]);
	edge_w = (double*) mxGetData(prhs[4]);
	rgblabel = (int*) mxGetData(prhs[5]);
	ck = (int)mxGetScalar(prhs[6]);
	lbplabel = (int*) mxGetData(prhs[7]);
	tk = (int)mxGetScalar(prhs[8]);
	height = (int)mxGetScalar(prhs[9]);
	width = (int)mxGetScalar(prhs[10]);
	angleWeights = (double*) mxGetData(prhs[11]);
	ep_prior = (double*) mxGetData(prhs[12]);
	thrs = (double*)mxGetData(prhs[13]);
	thrn = (int)mxGetScalar(prhs[14]);

	int* segrlt = new int[amount_v * thrn];
	segment_lep(amount_v, amount_e, edge_a, edge_b, edge_w, rgblabel, ck, lbplabel, tk,
		height, width, angleWeights, ep_prior, thrs, thrn, segrlt);
	
	plhs[0] = mxCreateDoubleMatrix(thrn, amount_v, mxREAL);

	double * outData = mxGetPr(plhs[0]);
	for(int i = 0; i < thrn; i++) 
        for(int j = 0; j < amount_v; j++) 
            outData[j * thrn + i] = segrlt[i * amount_v + j];

	delete [] segrlt;
}
