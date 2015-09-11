/*************************************************************************************
Copyright (C) 2015 Qiyang Zhao, Beihang Univ.

This program is free software; you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation; either 
version 2 of the License, or (at your option) any later version.

Please cite the following paper if you want to use the source codes in your work:

Qiyang Zhao, Segmenting natural images with the least effort as humans, BMVC 2015.
**************************************************************************************/

#ifndef _GMM_RGB_H
#define _GMM_RGB_H

#include <algorithm>
#include "misc.h"
#include "image.h"
#include "cv.h"
#include <omp.h>

double eulardist2_rgb(rgb c1, rgb c2)
{
	return (c1.r - c2.r) * (c1.r - c2.r) + (c1.g - c2.g) * (c1.g - c2.g) + (c1.b - c2.b) * (c1.b - c2.b);
}

struct Gaussian_rgb
{
	rgb mu;					// mean of the gaussian
	double covariance[3][3];		// covariance matrix of the gaussian
	double pi;					// weighting of this gaussian in the GMM.

	// These are only needed during Orchard and Bouman clustering.
	double eigenvalues[3];		// eigenvalues of covariance matrix
	double eigenvectors[3][3];	// eigenvectors of   "          "
};

class GMM_rgb
{
public:

	// Initialize GMM with number of gaussians desired.
	GMM_rgb(unsigned int K);
	~GMM_rgb();

private:

	unsigned int m_K;		// number of gaussians
	Gaussian_rgb* m_gaussians;	// an array of K gaussians

	friend int calcGMMs_rgb(const image<rgb>* im, int ck, int* rgblabel);
};

// Helper class that fits a single Gaussian to color samples
class GaussianFitter_rgb
{
public:
	GaussianFitter_rgb();
	
	// Add a color sample
	void add(rgb c);
	
	// Build the gaussian out of all the added color samples
	void finalize(Gaussian_rgb& g, unsigned int totalCount) const;
	
private:

	rgb s;			// sum of r,g, and b
	double  p[3][3];		// matrix of products (i.e. r*r, r*g, r*b), some values are duplicated.

	unsigned int count;	// count of color samples added to the gaussian
};

GMM_rgb::GMM_rgb(unsigned int K) : m_K(K)
{
	m_gaussians = new Gaussian_rgb[m_K];
}

GMM_rgb::~GMM_rgb()
{
	if (m_gaussians)
		delete [] m_gaussians;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int calcGMMs_rgb(const image<rgb>* im, int ck, int* rgblabel)
{
	/////////////////////////////////////////////////////////////////////////////Stage 1
	// Build GMMs on total image
	GMM_rgb * totalGMM = new GMM_rgb(ck);
	GaussianFitter_rgb* totalFitters = new GaussianFitter_rgb[ck];
	int height = im->height(), width = im->width();
	image<int>* components = new image<int>(width, height);
	unsigned int totalCount = height * width;

	// Initialize the first total clusters
//#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int ty = 0; ty < height; ++ty)	{
		for (int tx = 0; tx < width; ++tx) {
			rgb c;
			c.r = imRef(im, tx, ty).r; c.g = imRef(im, tx, ty).g; c.b = imRef(im, tx, ty).b;
			imRef(components, tx, ty) = 0;
			totalFitters[0].add(c);
		}
	}

	totalFitters[0].finalize(totalGMM->m_gaussians[0], totalCount);
	double totalv = totalGMM->m_gaussians[0].pi * totalCount *
			(totalGMM->m_gaussians[0].eigenvalues[0] + totalGMM->m_gaussians[0].eigenvalues[1] +
			totalGMM->m_gaussians[0].eigenvalues[2]);

	int nTotal = 0;		// Which cluster will be split
	bool continuing = true;
	int tmpck = 1;
	
	// Compute clusters
	for (int ti = 1; ti < ck && continuing; ti++) {

		// Reset the fitters for the splitting clusters
		totalFitters[nTotal] = GaussianFitter_rgb();

		// For brevity, get references to the splitting Gaussians
		Gaussian_rgb& tt = totalGMM->m_gaussians[nTotal];

		// Compute splitting points
		double splitTotal = tt.eigenvectors[0][0] * tt.mu.r + tt.eigenvectors[1][0] * tt.mu.g +
			tt.eigenvectors[2][0] * tt.mu.b;

		// Split clusters nTotal, place split portion into cluster i
		for (int y = 0; y < height; ++y) {
			for(int x = 0; x < width; ++x) {
				rgb c = {0, 0, 0};
				c.r = imRef(im, x, y).r; c.g = imRef(im, x, y).g; c.b = imRef(im, x, y).b;

				// For each pixel
				if (ti < ck && imRef(components, x, y) == nTotal)
				{
					if (tt.eigenvectors[0][0] * c.r + tt.eigenvectors[1][0] * c.g + tt.eigenvectors[2][0] * c.b > splitTotal)
					{
						imRef(components, x, y) = ti;
						totalFitters[ti].add(c);
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
		if (ti < ck) totalFitters[ti].finalize(totalGMM->m_gaussians[ti], totalCount);
		tmpck++;

		// Whether continue splitting or not
		double v = 0;
		for (int j = 0; j <= ti; j++) {
			v += totalGMM->m_gaussians[j].pi * totalCount *
			(totalGMM->m_gaussians[j].eigenvalues[0] + totalGMM->m_gaussians[j].eigenvalues[1] +
			totalGMM->m_gaussians[j].eigenvalues[2]);
		}
		continuing = ((tmpck < ck) && (v > totalv * M_LBP_SIGMA2_RELATIVE_THR)
			&& (v < totalv * M_SIGMA2_RELATIVE_UNCHANGED));

		// Find clusters with highest eigenvalue
		nTotal = 0;
		while (totalGMM->m_gaussians[nTotal].pi * totalCount <= 10 && nTotal <= ti) nTotal++;

		for (int j = nTotal + 1; j <= ti; j++) {
//			if (j < ck && totalGMM->m_gaussians[j].eigenvalues[0] > totalGMM->m_gaussians[nTotal].eigenvalues[0]
//			&& totalGMM->m_gaussians[j].pi * totalCount > 10)
//*
			if (j < ck && totalGMM->m_gaussians[j].eigenvalues[0] +
				totalGMM->m_gaussians[j].eigenvalues[1] +
				totalGMM->m_gaussians[j].eigenvalues[2] >
				totalGMM->m_gaussians[nTotal].eigenvalues[0] +
				totalGMM->m_gaussians[nTotal].eigenvalues[1] +
				totalGMM->m_gaussians[nTotal].eigenvalues[2] &&
				totalGMM->m_gaussians[j].pi * totalCount > 1000)
//*/
				nTotal = j;
		}
	}

#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			rgblabel[y * width + x] = imRef(components, x, y);
		}
	}
	ck = tmpck;

	delete [] totalFitters;
	delete totalGMM;
	delete components;

	return ck;
}

// GaussianFitter functions
GaussianFitter_rgb::GaussianFitter_rgb()
{
	s.r = 0; s.g = 0; s.b = 0;

	p[0][0] = 0; p[0][1] = 0; p[0][2] = 0;
	p[1][0] = 0; p[1][1] = 0; p[1][2] = 0;
	p[2][0] = 0; p[2][1] = 0; p[2][2] = 0;

	count = 0;
}

// Add a color sample
void GaussianFitter_rgb::add(rgb c)
{
	s.r += c.r; s.g += c.g; s.b += c.b;

	p[0][0] += c.r*c.r; p[0][1] += c.r*c.g; p[0][2] += c.r*c.b;
	p[1][0] += c.g*c.r; p[1][1] += c.g*c.g; p[1][2] += c.g*c.b;
	p[2][0] += c.b*c.r; p[2][1] += c.b*c.g; p[2][2] += c.b*c.b;

	count++;
}

// Build the gaussian out of all the added colors
void GaussianFitter_rgb::finalize(Gaussian_rgb& g, unsigned int totalCount) const
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
		g.mu.r = s.r/count;
		g.mu.g = s.g/count;
		g.mu.b = s.b/count;

		// Compute covariance matrix
		g.covariance[0][0] = p[0][0]/count-g.mu.r*g.mu.r + Epsilon; g.covariance[0][1] = p[0][1]/count-g.mu.r*g.mu.g; g.covariance[0][2] = p[0][2]/count-g.mu.r*g.mu.b;
		g.covariance[1][0] = p[1][0]/count-g.mu.g*g.mu.r; g.covariance[1][1] = p[1][1]/count-g.mu.g*g.mu.g + Epsilon; g.covariance[1][2] = p[1][2]/count-g.mu.g*g.mu.b;
		g.covariance[2][0] = p[2][0]/count-g.mu.b*g.mu.r; g.covariance[2][1] = p[2][1]/count-g.mu.b*g.mu.g; g.covariance[2][2] = p[2][2]/count-g.mu.b*g.mu.b + Epsilon;

		// The weight of the gaussian is the fraction of the number of pixels in this Gaussian to the number of 
		// pixels in all the gaussians of this GMM.
		g.pi = (double)count/totalCount;

		// Build OpenCV wrappers around our data.
		CvMat mat = cvMat(3, 3, CV_64FC1, g.covariance);
		CvMat eval = cvMat(3, 1, CV_64FC1, g.eigenvalues);
		CvMat evec = cvMat(3, 3, CV_64FC1, g.eigenvectors);

		// Compute eigenvalues and vectors using SVD
		cvSVD( &mat, &eval, &evec );
	}
}

#endif //GMM_RGB_H
