/*************************************************************************************
Copyright (C) 2015 Qiyang Zhao, Beihang Univ.

This program is free software; you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation; either 
version 2 of the License, or (at your option) any later version.

Please cite the following paper if you want to use the source codes in your work:

Qiyang Zhao, Segmenting natural images with the least effort as humans, BMVC 2015.
**************************************************************************************/

#ifndef _DISJOINT_SET_H
#define _DISJOINT_SET_H

typedef struct {
	int p;
	int size;
	double energy;
	double energy_t;
} uni_elt;

class universe {
public:
	universe(int elements);
	~universe();
	int find(int x);  
	int join(int x, int y);

	void setenergy(int x, double e) {elts[x].energy = e;}
	double getenergy(int x) {return elts[x].energy;}

	void setenergy_t(int x, double e) {elts[x].energy_t = e;}
	double getenergy_t(int x) {return elts[x].energy_t;}

	int size(int x) const { return elts[x].size;}
	int num_sets() const { return num; }

private:
	uni_elt *elts;
	int num;
	int unused;
};

universe::universe(int elements) {
	elts = new uni_elt[elements * 2];
	num = elements;
	unused = elements;
	for (int i = 0; i < elements * 2; i++) {
		elts[i].size = 1;
		elts[i].p = i;
		elts[i].energy = 0;
		elts[i].energy_t = 0;
	}
}

universe::~universe() {
	delete [] elts;
}

int universe::find(int x) {
	int label;
	if (x == elts[x].p) {
		label = x;
	} else {
		label = find(elts[x].p);
		elts[x].p = label;
	}

	return label;
}

int universe::join(int x, int y) {
	elts[x].p = unused;
	elts[y].p = unused;
	elts[unused].size = elts[x].size + elts[y].size;

	unused++;
	num--;

	return unused - 1;
}

#endif
