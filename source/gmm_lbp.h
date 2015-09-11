/*************************************************************************************
Copyright (C) 2015 Qiyang Zhao, Beihang Univ.

This program is free software; you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation; either 
version 2 of the License, or (at your option) any later version.

Please cite the following paper if you want to use the source codes in your work:

Qiyang Zhao, Segmenting natural images with the least effort as humans, BMVC 2015.
**************************************************************************************/

#ifndef _GMM_LBP_H
#define _GMM_LBP_H

#include <algorithm>
#include "misc.h"
#include "image.h"
#include "cv.h"
#include "imgio.h"
#include "gmm_rgb.h"
#include <omp.h>

typedef struct { double v[8]; } lbp;

double eulardist2_lbp(lbp c1, lbp c2)
{
	double v = 0;
	for (int i = 0; i < 8; i++)
		v += (c1.v[i] - c2.v[i]) * (c1.v[i] - c2.v[i]);
	return v;
}

struct Gaussian_lbp
{
	lbp mu;					// mean of the gaussian
	double covariance[8][8];		// covariance matrix of the gaussian
	double pi;					// weighting of this gaussian in the GMM.

	// These are only needed during Orchard and Bouman clustering.
	double eigenvalues[8];		// eigenvalues of covariance matrix
	double eigenvectors[8][8];	// eigenvectors of   "          "
};

class GMM_lbp
{
public:

	// Initialize GMM with number of gaussians desired.
	GMM_lbp(unsigned int K);
	~GMM_lbp();

private:

	unsigned int m_K;		// number of gaussians
	Gaussian_lbp* m_gaussians;	// an array of K gaussians

	friend int calcGMMs_lbp(const image<rgb>* im, int tk, int* lbplabel);
};

// Helper class that fits a single Gaussian to color samples
class GaussianFitter_lbp
{
public:
	GaussianFitter_lbp();
	
	// Add a color sample
	void add(lbp c);
	
	// Build the gaussian out of all the added color samples
	void finalize(Gaussian_lbp& g, unsigned int totalCount) const;
	
private:

	lbp s;			// sum of r,g, and b
	double  p[8][8];		// matrix of products (i.e. r*r, r*g, r*b), some values are duplicated.

	unsigned int count;	// count of color samples added to the gaussian
};

GMM_lbp::GMM_lbp(unsigned int K) : m_K(K)
{
	m_gaussians = new Gaussian_lbp[m_K];
}

GMM_lbp::~GMM_lbp()
{
	if (m_gaussians)
		delete [] m_gaussians;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int calcGMMs_lbp(const image<rgb>* im, int tk, int* lbplabel)
{
	///////////////////////////////////////////////compute the lbp features of pixels
	GaussianFitter_rgb* colorFitters = new GaussianFitter_rgb;
	Gaussian_rgb colorgaussian;
	int height = im->height(), width = im->width();
	int totalCount = height * width;

//#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int ty = 0; ty < height; ++ty)	{
		for (int tx = 0; tx < width; ++tx) {
			rgb c;
			c.r = imRef(im, tx, ty).r; c.g = imRef(im, tx, ty).g; c.b = imRef(im, tx, ty).b;
			colorFitters->add(c);
		}
	}

	colorFitters->finalize(colorgaussian, totalCount);

	image<double> * tmpv = new image<double>(width, height);

#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int ty = 0; ty < height; ++ty)	{
		for (int tx = 0; tx < width; ++tx) {
			rgb c;

			c.r = imRef(im, tx, ty).r; c.g = imRef(im, tx, ty).g; c.b = imRef(im, tx, ty).b;
			imRef(tmpv, tx, ty) = colorgaussian.eigenvectors[0][0] * c.r + 
				colorgaussian.eigenvectors[1][0] * c.g + colorgaussian.eigenvectors[2][0] * c.b;
		}
	}

	image<lbp> * lbpim = new image<lbp>(width, height);

#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int ty = 1; ty < height - 1; ++ty)	{
#pragma omp parallel for num_threads(omp_get_max_threads()) 
		for (int tx = 1; tx < width - 1; ++tx) {
			double tv = imRef(tmpv, tx, ty);

			imRef(lbpim, tx, ty).v[0] = (double)(imRef(tmpv, tx - 1, ty - 1) - tv);  
			imRef(lbpim, tx, ty).v[1] = (double)(imRef(tmpv, tx, ty - 1) - tv);  
			imRef(lbpim, tx, ty).v[2] = (double)(imRef(tmpv, tx + 1, ty - 1) - tv);  
			imRef(lbpim, tx, ty).v[3] = (double)(imRef(tmpv, tx - 1, ty) - tv);  
			imRef(lbpim, tx, ty).v[4] = (double)(imRef(tmpv, tx + 1, ty) - tv);  
			imRef(lbpim, tx, ty).v[5] = (double)(imRef(tmpv, tx - 1, ty + 1) - tv);  
			imRef(lbpim, tx, ty).v[6] = (double)(imRef(tmpv, tx, ty + 1) - tv);  
			imRef(lbpim, tx, ty).v[7] = (double)(imRef(tmpv, tx + 1, ty + 1) - tv);
		}
	}

	/////////////////////////////////////////////////////////////////////////////Stage 1
	// Build GMMs on total image
	GMM_lbp * totalGMM = new GMM_lbp(tk);
	GaussianFitter_lbp* totalFitters = new GaussianFitter_lbp[tk];
	image<int>* components = new image<int>(width, height);

	// Initialize the first total clusters
//#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int ty = 0; ty < height; ++ty)	{
		for (int tx = 0; tx < width; ++tx) {
			lbp c;
			for (int ti = 0; ti < 8; ti++) c.v[ti] = imRef(lbpim, tx, ty).v[ti];
			imRef(components, tx, ty) = 0;
			totalFitters[0].add(c);
		}
	}

	totalFitters[0].finalize(totalGMM->m_gaussians[0], totalCount);
	double totalv = 0;
//#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int ti = 0; ti < 8; ti++) totalv += totalGMM->m_gaussians[0].eigenvalues[ti];
	totalv *= totalCount;

	int nTotal = 0;		// Which cluster will be split
	bool continuing = true;
	int tmptk = 1;
	
	// Compute clusters
	for (int tmpi = 1; tmpi < tk && continuing; tmpi++) {
		// Reset the fitters for the splitting clusters
		totalFitters[nTotal] = GaussianFitter_lbp();

		// For brevity, get references to the splitting Gaussians
		Gaussian_lbp& tt = totalGMM->m_gaussians[nTotal];

		double splitTotal = 0;
		for (int tmpii = 0; tmpii < 8; tmpii++) splitTotal += tt.eigenvectors[tmpii][0] * tt.mu.v[tmpii];

		// Split clusters nTotal, place split portion into cluster i
		for (int y = 0; y < height; ++y) {
			for(int x = 0; x < width; ++x) {
				lbp c;
				for (int tmpii = 0; tmpii < 8; tmpii++) c.v[tmpii] = imRef(lbpim, x, y).v[tmpii];

				// For each pixel
				if (tmpi < tk && imRef(components, x, y) == nTotal)
				{
					double tmpvv = 0;
					for (int tmpii = 0; tmpii < 8; tmpii++) tmpvv += tt.eigenvectors[tmpii][0] * c.v[tmpii];
					if (tmpvv > splitTotal)
					{
						imRef(components, x, y) = tmpi;
						totalFitters[tmpi].add(c);
					}
					else
					{
						totalFitters[nTotal].add(c);
					}
				}
			}
		}

		// Compute new split Gaussians
		totalFitters[nTotal].finalize(totalGMM->m_gaussians[nTotal], totalCount);
		if (tmpi < tk) totalFitters[tmpi].finalize(totalGMM->m_gaussians[tmpi], totalCount);
		tmptk++;

		// Whether continue splitting or not
		double v = 0;
		for (int j = 0; j <= tmpi; j++) {
			double v1 = 0;
			for (int ti = 0; ti < 8; ti++) v1 += totalGMM->m_gaussians[j].eigenvalues[ti];
			v1 *= totalGMM->m_gaussians[j].pi * totalCount;
			v += v1;
		}
		continuing = ((tmptk < tk) && (v > totalv * M_LBP_SIGMA2_RELATIVE_THR)
			&& (v < totalv * M_SIGMA2_RELATIVE_UNCHANGED));

		// Find clusters with highest eigenvalue
		nTotal = 0;
		while (totalGMM->m_gaussians[nTotal].pi * totalCount <= 10 && nTotal <= tmpi) nTotal++;

		double v1 = 0;
		for (int ti = 0; ti < 8; ti++) v1 += totalGMM->m_gaussians[nTotal].eigenvalues[ti];
		v1 *= totalGMM->m_gaussians[nTotal].pi;
		for (int j = nTotal + 1; j <= tmpi; j++) {
//			if (j < ck && totalGMM->m_gaussians[j].eigenvalues[0] > totalGMM->m_gaussians[nTotal].eigenvalues[0]
//			&& totalGMM->m_gaussians[j].pi * totalCount > 10)
//*
			double v2 = 0;
			for (int ti = 0; ti < 8; ti++) v2 += totalGMM->m_gaussians[j].eigenvalues[ti];
			v2 *= totalGMM->m_gaussians[j].pi;
			if (j < tk && v2 > v1 &&
				totalGMM->m_gaussians[j].pi * totalCount > 10) {
				v1 = v2; nTotal = j;
			}
		}
	}

	int otk = tk;
	tk = tmptk;
//#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			lbplabel[y * width + x] = imRef(components, x, y);
		}
	}

	delete[] totalFitters;
	delete totalGMM;
	delete colorFitters;
	delete components;
	delete tmpv;
	delete lbpim;

	return tk;
}

// GaussianFitter functions
GaussianFitter_lbp::GaussianFitter_lbp()
{
//#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int i = 0; i < 8; i++) s.v[i] = 0;
//#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			p[i][j] = 0;

	count = 0;
}

// Add a color sample
void GaussianFitter_lbp::add(lbp c)
{
//#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int i = 0; i < 8; i++) s.v[i] += c.v[i];
//#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			p[i][j] += c.v[i] * c.v[j];

	count++;
}

// Build the gaussian out of all the added colors
void GaussianFitter_lbp::finalize(Gaussian_lbp& g, unsigned int totalCount) const
{
	// Running into a singular covariance matrix is problematic. So we'll add a small epsilon
	// value to the diagonal elements to ensure a positive definite covariance matrix.
	const double Epsilon = (double)0.00000001;

	if (count==0)
	{
		g.pi = 0;
	}
	else
	{
		// Compute mean of gaussian
		for (int i = 0; i < 8; i++) g.mu.v[i] = s.v[i] / count;

		// Compute covariance matrix
		for (int i = 0; i < 8; i++)
			for (int j = 0; j < 8; j++)
				g.covariance[i][j] = p[i][j] / count - g.mu.v[i] * g.mu.v[j];

		for (int i = 0; i < 8; i++) g.covariance[i][i] += Epsilon;
		CvMat covmat = cvMat(8, 8, CV_64FC1, g.covariance);
		// The weight of the gaussian is the fraction of the number of pixels in this Gaussian to the number of 
		// pixels in all the gaussians of this GMM.
		g.pi = (double)count / totalCount;

		// Build OpenCV wrappers around our data.
		CvMat eval = cvMat(8, 1, CV_64FC1, g.eigenvalues);
		CvMat evec = cvMat(8, 8, CV_64FC1, g.eigenvectors);
			
		// Compute eigenvalues and vectors using SVD
		cvSVD( &covmat, &eval, &evec );
	}
}

#endif
