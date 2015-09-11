/*************************************************************************************
Copyright (C) 2015 Qiyang Zhao, Beihang Univ.

This program is free software; you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation; either 
version 2 of the License, or (at your option) any later version.

Please cite the following paper if you want to use the source codes in your work:

Qiyang Zhao, Segmenting natural images with the least effort as humans, BMVC 2015.
**************************************************************************************/

#ifndef _IMGIO_H
#define _IMGIO_H

#include <cv.h>
#include <string>
#include <highgui.h>
#include "image.h"
#include "misc.h"

image<rgb>* loadimg( std::string file_name ) //by zhaoqy on 2011-07-15 in NLSDE, BUAA
{
	IplImage *img;            
	img = cvLoadImage(file_name.c_str(), CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_ANYCOLOR);

	if (img) {
		unsigned int width, height, maxv;
		width = img->width;
		height = img->height;
		maxv = 255;
		
		image<rgb>* im = new image<rgb>(width, height);
		
		if (im)
		{
			for (unsigned int y = 0; y < height; ++y)
			{
				for (unsigned int x = 0; x < width; ++x)
				{
					char r, g, b;
					
					b = img->imageData[img->widthStep * y + img->nChannels * x + 0];
					g = img->imageData[img->widthStep * y + img->nChannels * x + 1];
					r = img->imageData[img->widthStep * y + img->nChannels * x + 2];

					imRef(im, x, y).r = (double)((unsigned char)r) / maxv;
					imRef(im, x, y).g = (double)((unsigned char)g) / maxv;
					imRef(im, x, y).b = (double)((unsigned char)b) / maxv;
				}
			}
		}
		return im;
	}

	return 0;
}

void saveimg( image<rgb>* im, std::string file_name ) //by zhaoqy on 2011-07-15 in NLSDE, BUAA
{
	if (im == NULL) return;
	unsigned int width, height;
	width = im->width(); height = im->height();
	CvSize sz = {width, height};
	IplImage *img = cvCreateImage(sz, IPL_DEPTH_8U, 3);

	if (img) {
		for (unsigned int y = 0; y < height; ++y) {
			for (unsigned int x = 0; x < width; ++x) {
				img->imageData[img->widthStep * y + img->nChannels * x + 0] = (unsigned char)imRef(im, x, y).b;
				img->imageData[img->widthStep * y + img->nChannels * x + 1] = (unsigned char)imRef(im, x, y).g;
				img->imageData[img->widthStep * y + img->nChannels * x + 2] = (unsigned char)imRef(im, x, y).r;
			}
		}
		cvSaveImage(file_name.c_str(), img);
		cvReleaseImage(&img);
	} else {
		printf("error saving U in images\n");
	}


	return;
}

#endif
