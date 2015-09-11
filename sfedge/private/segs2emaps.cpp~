/**************************************************************************************/
Copyright (C) 2015 Qiyang Zhao, Beihang Univ.

This program is free software; you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation; either 
version 2 of the License, or (at your option) any later version.

Please cite the following paper if you want to use the source codes in your work:

Qiyang Zhao, Segmenting natural images with the least effort as humans, BMVC 2015.
/**************************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	int amount_v, amount_e, ck, ock, tk, otk, height, width, thtnum, *edge_a, *edge_b, *rgblabel, *lbplabel;
	double *edge_w, alpha1, alpha2, mintht, thtstep, *angleWeights, *ep_prior;
	st_forestModel* model; 
	
	amount_v = (int)mxGetScalar(prhs[0]);
	amount_e = (int)mxGetScalar(prhs[1]);
	edge_a = (int*) mxGetData(prhs[2]);
	edge_b = (int*) mxGetData(prhs[3]);
	edge_w = (double*) mxGetData(prhs[4]);
	rgblabel = (int*) mxGetData(prhs[5]);
	ck = (int)mxGetScalar(prhs[6]);
	ock = (int)mxGetScalar(prhs[7]);
	lbplabel = (int*) mxGetData(prhs[8]);
	tk = (int)mxGetScalar(prhs[9]);
	otk = (int)mxGetScalar(prhs[10]);
	alpha1 = (double)mxGetScalar(prhs[11]);
	alpha2 = (double)mxGetScalar(prhs[12]);
	angleWeights = (double*) mxGetData(prhs[13]);
	ep_prior = (double*) mxGetData(prhs[14]);
	mintht = (double)mxGetScalar(prhs[15]);
	thtstep = (double)mxGetScalar(prhs[16]);
	thtnum = (int)mxGetScalar(prhs[17]);
	height = (int)mxGetScalar(prhs[18]);
	width = (int)mxGetScalar(prhs[19]);

    /* check empty field, proper data type, and data type consistency;
     * and get classID for each field. */
	st_forestModel* pModel = new st_forestModel;
    for(int i = 0; i < FORESTNUM; i++) {
        pModel->sizes[i] = (int)mxGetScalar(mxGetFieldByNumber(prhs[20], i, 0));
        pModel->nodenums[i] = (int)mxGetScalar(mxGetFieldByNumber(prhs[20], i, 1));
        pModel->pStartPos[i] = (int*)mxGetData(mxGetFieldByNumber(prhs[20], i, 2));
        pModel->pFids[i] = (int*)mxGetData(mxGetFieldByNumber(prhs[20], i, 3));
        pModel->pChild[i] = (int*)mxGetData(mxGetFieldByNumber(prhs[20], i, 4));
        pModel->pThrs[i] = (float*)mxGetData(mxGetFieldByNumber(prhs[20], i, 5));
        pModel->pDistr[i] = (float*)mxGetData(mxGetFieldByNumber(prhs[20], i, 6));
	}

	int* segrlt = new int[amount_v * thtnum];
	segment_lep(amount_v, amount_e, edge_a, edge_b, edge_w,
		rgblabel, ck, ock, lbplabel, tk, otk, alpha1, alpha2, angleWeights, ep_prior, 
		mintht, thtstep, thtnum, height, width, pModel, segrlt);
	
	plhs[0] = mxCreateDoubleMatrix(h1, w1 - 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(h1 - 1, w1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(h1 - 1, w1 - 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(h1 - 1, w1 - 1, mxREAL);

	double * outData = mxGetPr(plhs[0]);
	for(int i = 0; i < thtnum; i++) 
        for(int j = 0; j < amount_v; j++) 
            outData[j * thtnum + i] = segrlt[i * amount_v + j];

	delete [] segrlt;

function [emap1, emap2, emap3, emap4] = segs2EdgeMap_singleScale(segs, siz, tH, tW)

int h2 = siz[2], h1 = h2 * 2, w2 = siz[3], w1 = w2 * 2, Tn = siz[4];
double* emap1 = new double[h1 * (w1 - 1)], *emap2 = new double[(h1 - 1) * w1];
double* emap3 = new double[(h1 - 1) * (w1 - 1)], *emap4 = new double[(h1 - 1) * (w1 - 1)];
int* nmap1 = new int[h1 * (w1 - 1)], *nmap2 = new int[(h1 - 1) * w1];
int* nmap3 = new int[(h1 - 1) * (w1 - 1)], *nmap4 = new int[(h1 - 1) * (w1 - 1)];

memset(emap1, 0, h1 * (w1 - 1) * sizeof(double));
memset(emap2, 0, w1 * (h1 - 1) * sizeof(double));
memset(emap3, 0, (h1 - 1) * (w1 - 1) * sizeof(double));
memset(emap4, 0, (h1 - 1) * (w1 - 1) * sizeof(double));
memset(nmap1, 0, h1 * (w1 - 1) * sizeof(int));
memset(nmap2, 0, w1 * (h1 - 1) * sizeof(int));
memset(nmap3, 0, (h1 - 1) * (w1 - 1) * sizeof(int));
memset(nmap4, 0, (h1 - 1) * (w1 - 1) * sizeof(int));

for (int T = 0; T < Tn; T++) {
    for (int y = 0; y < h2; y++) {
        for (int x = 0; x < w2; x++) {
            // 乘以2的原因是detector是每隔1个计算localseg
            int ly1 = (y - 1) * 2 - 7, ry1 = (y - 1) * 2 + 8;
            int lx1 = (x - 1) * 2 - 7, rx1 = (x - 1) * 2 + 8;

            // local seg与图像重合的部分中在图像中的位置范围
            int ly2 = max(1, ly1), ry2 = min(h1, ry1);
            int lx2 = max(1, lx1), rx2 = min(w1, rx1);

            // local seg与图像重合的部分中在seg中的位置范围
            int ly3 = ly2 - ly1 + 1, ry3 = ry2 - ly1 + 1;        
            int lx3 = lx2 - lx1 + 1, rx3 = rx2 - lx1 + 1;

            localseg = reshape(segs(:, :, y, x, T), [16 16]);
            emap1(ly2 : ry2, lx2 : rx2 - 1) = emap1(ly2 : ry2, lx2 : rx2 - 1) + ...
                (localseg(ly3 : ry3, lx3 : rx3 - 1) ~= localseg(ly3 : ry3, lx3 + 1 : rx3));
            nmap1(ly2 : ry2, lx2 : rx2 - 1) = nmap1(ly2 : ry2, lx2 : rx2 - 1) + 1;
            
            emap2(ly2 : ry2 - 1, lx2 : rx2) = emap2(ly2 : ry2 - 1, lx2 : rx2) + ...
                (localseg(ly3 : ry3 - 1, lx3 : rx3) ~= localseg(ly3 + 1 : ry3, lx3 : rx3));
            nmap2(ly2 : ry2 - 1, lx2 : rx2) = nmap2(ly2 : ry2 - 1, lx2 : rx2) + 1;
            
            emap3(ly2 : ry2 - 1, lx2 : rx2 - 1) = emap3(ly2 : ry2 - 1, lx2 : rx2 - 1) + ...
                (localseg(ly3 : ry3 - 1, lx3 : rx3 - 1) ~= localseg(ly3 + 1 : ry3, lx3 + 1 : rx3));
            nmap3(ly2 : ry2 - 1, lx2 : rx2 - 1) = nmap3(ly2 : ry2 - 1, lx2 : rx2 - 1) + 1;
            
            emap4(ly2 : ry2 - 1, lx2 : rx2 - 1) = emap4(ly2 : ry2 - 1, lx2 : rx2 - 1) + ...
                (localseg(ly3 : ry3 - 1, lx3 + 1 : rx3) ~= localseg(ly3 + 1 : ry3, lx3 : rx3 - 1));
            nmap4(ly2 : ry2 - 1, lx2 : rx2 - 1) = nmap4(ly2 : ry2 - 1, lx2 : rx2 - 1) + 1;
        end
    end
end

t = 1.66;
emap1 = emap1 ./ nmap1; emap2 = emap2 ./ nmap2; emap3 = emap3 ./ nmap3; emap4 = emap4 ./ nmap4;
emap1 = emap1(1 : imgH, 1 : imgW - 1) * t; emap2 = emap2(1 : imgH - 1, 1 : imgW) * t;
emap3 = emap3(1 : imgH - 1, 1 : imgW - 1) * t; emap4 = emap4(1 : imgH - 1, 1 : imgW - 1) * t;
emap1 = convTri(single(emap1), 1); emap2 = convTri(single(emap2), 1);
emap3 = convTri(single(emap3), 1); emap4 = convTri(single(emap4), 1);

end
}
