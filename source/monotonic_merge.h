/*************************************************************************************
Copyright (C) 2015 Qiyang Zhao, Beihang Univ.

This program is free software; you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation; either 
version 2 of the License, or (at your option) any later version.

Please cite the following paper if you want to use the source codes in your work:

Qiyang Zhao, Segmenting natural images with the least effort as humans, BMVC 2015.
**************************************************************************************/

#ifndef _MONOTONIC_MERGE_H
#define _MONOTONIC_MERGE_H

#include <algorithm>
#include <vector>
#include <cmath>
#include "disjoint-set.h"
#include "tablestruct.h"
#include "xlogxtable.h"
#include "acostable.h"
#include "normalizer.h"
#include "imgio.h"
#include <omp.h>

using namespace std;

///////////////////////////////////////////////////////////////
#define STABLE_STEP					0.0001
#define	STABLE_SCALE					10000
#define MAXSTRENGTH					100
#define EPS						1e-10
#define intlog(x)					log((double)(x))
#define entry_upperlimit(i)				((i) > 100 * STABLE_SCALE - 1 ? 100 * STABLE_SCALE - 1 : (i))
#define stable_entry(x)					((x) < 0 ? 0 : entry_upperlimit(int((x) * STABLE_SCALE)))
#define stable_value(i)					(STABLE_STEP * i)
#define mymax(x1, x2)					((x1) > (x2) ? (x1) : (x2))
#define mymin(x1, x2)					((x1) < (x2) ? (x1) : (x2))
#define normalise(v)					(1 / (1 + exp(- 5 * (v - 0.5))))
#define strength(w0, w1, w0_t, w1_t, w2, am, bl, ef, ep) \
	((normalise(w2 / am) * ((w1 - w0) + (w1_t - w0_t)) + w2) / (am + (1 / (1 + exp(- 0.1 * (am - 500)))) * ef) \
	+ (1 / (1 + exp(- 0.1 * (am - 500)))) * (ep / am))

/////////////////////////////////////////////////////////////// random color
rgb random_rgb(){
	rgb c = {(uchar)rand(), (uchar)rand(), (uchar)rand()};
	return c;
}

////////////////////////////////////////////////////////////////////////////
void saveUinImg(universe* u, int* oneseg, rgb *colors, int height, int width,
				const char * segimgpath, double thr)
{
#if SAVESEGIMG
	image<rgb> *segs = new image<rgb>(width, height);
#endif

	// pick random colors for each component
	int *label = (int*)calloc(width * height * 2, sizeof(int));

	int ind = 1;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int comp = u->find(y * width + x);
			if (!label[comp]) label[comp] = ind++;
			oneseg[x * height + y] = label[comp];
#if SAVESEGIMG
			imRef(segs, x, y) = colors[comp];
			imRef(segs, x, y).b = label[comp];
#endif
		}
	}

#if SAVESEGIMG
	char filename[128] = {0};
	sprintf_s(filename, "%s/%3.3f.bmp", segimgpath, thr);
	saveimg(segs, filename);
	delete segs;
#endif

	free(label);
}

bool notACorner(int x1, int y1, int x2, int y2, int x3, int y3)
{
	if (x1 == x2 && x2 == x3) return true;
	if (y1 == y2 && y2 == y3) return true;
	if ((y3 - y2) * (y2 - y1) + (x3 - x2) * (x2 - x1) < - EPS) return false;
	int a1 = x2 - x1, b1 = y2 - y1, a2 = x3 - x2, b2 = y3 - y2, a3 = x3 - x1, b3 = y3 - y1;
	int c1 = a1 * y2 - b1 * x2, c2 = a2 * y2 - b2 * x2, c3 = a3 * y3 - b3 * x3;
	int h1 = b1 * x3 - a1 * y3 + c1, h2 = b2 * x1 - a2 * y1 + c2, h3 = b3 * x2 - a3 * y2 + c3;
	if ((h1 >= 0 && h1 <= (abs(a1) + abs(b1))) ||
		(h2 >= 0 && h2 <= (abs(a2) + abs(b2))) ||
		(h3 >= 0 && h3 <= (abs(a3) + abs(b3)))) return true; else return false;
}

void merge2angleList(strengthNode* sida, strengthNode* sidb, gridNode* gtable, double* angleWeights)
{
	for (int i = 0; i < 16; i++) sida->angleHist[i] += sidb->angleHist[i];
	sida->effort += sidb->effort; sida->borderlen += sidb->borderlen;

	angleNode* aa = sida->angleList, * ba = sidb->angleList;
	angleNode* anext = NULL, * bnext = NULL, * lastaa = NULL, * lastba = NULL;
	int aid, bid;
	while (aa && ba) {
		anext = aa->next; bnext = ba->next; aid = aa->id; bid = ba->id;
		if (!aa->next) lastaa = aa;
		if (aid == bid) { // in both lists, so to merge
			removeANode(aa, sida); lastaa = aa->prev; // remove the common grid point in alist, but remains in blist
			int x1 = gtable[gtable[aid].dssEnds[aa->dir]].x;
			int y1 = gtable[gtable[aid].dssEnds[aa->dir]].y;
			int x2 = gtable[aid].x, y2 = gtable[aid].y;
			int x3 = gtable[gtable[bid].dssEnds[ba->dir]].x;
			int y3 = gtable[gtable[bid].dssEnds[ba->dir]].y;
			if (notACorner(x1, y1, x2, y2, x3, y3)) {
				gtable[gtable[aid].dssEnds[aa->dir]].dssEnds[aa->invdir] = gtable[bid].dssEnds[ba->dir];
				gtable[gtable[bid].dssEnds[ba->dir]].dssEnds[ba->invdir] = gtable[aid].dssEnds[aa->dir];
				sida->angleHist[0] += 1; sida->effort += angleWeights[0]; 
			} else {
#if 1
				int tmpind = (((y3 - y2) * (y2 - y1) + (x3 - x2) * (x2 - x1)) /
					sqrt((y3 - y2) * (y3 - y2) + (x3 - x2) * (x3 - x2)) /
					sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1))) * 10000 + 10000;
				tmpind = mymin(mymax(tmpind, 0), 20001);
				double ang = acostable[tmpind];
				int ind = floor((ang - M_PI / 32) * (16 / M_PI));
				ind = mymax(mymin(ind, 14), 0);
#else
				double tv1 = (y3 - y2) * (y2 - y1) + (x3 - x2) * (x2 - x1);
				double tv2 = (y3 - y2) * (y3 - y2) + (x3 - x2) * (x3 - x2);
				double tv3 = (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1);
				int ind = 0;
				if (tv1 > 0) ind = angletable[(int)(tv1 * tv1 / tv2 / tv3 * 100000 + 100000)];
				else ind = angletable[(int)(- tv1 * tv1 / tv2 / tv3 * 100000 + 100000)];
#endif
				sida->angleHist[ind] += 1; sida->effort += angleWeights[ind];
			}
			sida->angleHist[15] -= 2; sida->effort -= 2 * angleWeights[15] + angleWeights[16];
			sida->borderlen++;

			aa = anext;	ba = bnext;
		} else {
			if (aid > bid) { // angle point only in a list
				aa = anext;
			} else { // angle point only in b list
				insertANode(aa, ba, sida);
				ba = bnext;
			}
		}
	}

	if (ba) {
		if (lastaa) {
			ba->prev = lastaa; lastaa->next = ba;
		} else {
			ba->prev = NULL; sida->angleList = ba;
		}
	}
	sidb->angleList = NULL;
}

void merge2regions(universe* u,
						  strengthNode* stable, strengthNode** sheads,
						  neighborNode* ntable, neighborNode** nheads,
						  gridNode* gtable, angleNode* atable,
						  int* rgbcom, int* lbpcom, int ck, int tk,
						  double* angleWeights, strengthNode * si)
{
	int a = si->a, b = si->b;
	double ce = si->engy_c, te = si->engy_t;

	int ab = u->join(a, b);
	si->valid = false;

	register int * abp = (int*)(rgbcom + ab * ck);
	register int * ap = (int*)(rgbcom + a * ck);
 	register int * bp = (int*)(rgbcom + b * ck);
	for (int i = 0; i < ck; i++) *abp++ = *ap++ + *bp++;
	abp = (int*)(lbpcom + ab * tk);
	ap = (int*)(lbpcom + a * tk);
	bp = (int*)(lbpcom + b * tk);
	for (int i = 0; i < tk; i++) *abp++ = *ap++ + *bp++;
	u->setenergy(ab, ce); u->setenergy_t(ab, te);

	neighborNode * a_in_b, * b_in_a;
	neighborNode* na = *(nheads + a), * nb = *(nheads + b), * lastna = NULL, * nid;
	neighborNode** alist = nheads + a, ** blist = nheads + b;
	neighborNode* anext, * bnext;
	strengthNode* sida, * sidb, * sid_kept, * sid_dump;
	neighborNode** lp;
	neighborNode* tp, * ntp;
	int ca, cb;
	while (na && nb) {
		anext = na->next;
		bnext = nb->next;
		ca = na->id;
		cb = nb->id;
		if (!na->next) lastna = na;
		if (ca == cb) { // c is in both a list and b list
			sida = na->sid;
			sidb = nb->sid;
			if (sida->pb_entry < sidb->pb_entry) { // keep stable a-c
				sid_kept = sida;
				sid_dump = sidb;
			} else { // keep stable b-c
				sid_kept = sidb;
				sid_dump = sida;
			}
			sid_kept->a = ab;
			sid_kept->b = ca;
			sid_kept->nm += sid_dump->nm;
			sid_kept->ep_prior += sid_dump->ep_prior;
			sid_kept->nv += sid_dump->nv;
			sid_kept->updated = false;
			sid_dump->valid = false;

			merge2angleList(sid_kept, sid_dump, gtable, angleWeights);

			lp = nheads + ca;
			tp = *lp;
			if (a > b) { // keep a in c's neighbor list, and c in a
				nid = na->twin;
				nid->id = ab;
				nid->sid = sid_kept;
				if (tp != nid) moveNNode2First(nid, tp, lp);
				// remove b in c
				removeNNode(nb->twin);
			} else { // keep b in c's neighbor list, but keep c in a
				nid = nb->twin;
				nid->id = ab;
				nid->sid = sid_kept;
				if (tp != nid) moveNNode2First(nid, tp, lp);
				nid->twin = na;
				// remove a in c
				removeNNode(na->twin);
				na->twin = nid;	
			}
			na->sid = sid_kept;

			na = anext;
			nb = bnext;
		} else {
			if (ca > cb) { // c only in a list, keep a-c
				if (ca == b) {
					b_in_a = na;
				} else {				
					sida = na->sid;
					sida->a = ab;
					sida->b = ca;
					sida->updated = false;

					nid = na->twin;

					lp = nheads + ca;
					tp = *lp;

					nid->id = ab;
					if (tp != nid) moveNNode2First(nid, tp, lp);
				}
				na = anext;
			} else { // c only in b list
				if (cb == a) {
					a_in_b = nb;
				} else {
					sidb = nb->sid;
					sidb->a = ab;
					sidb->b = cb;
					sidb->updated = false;

					nid = nb->twin;

					lp = nheads + cb;
					tp = *lp;

					nid->id = ab;
					if (tp != nid) moveNNode2First(nid, tp, lp);
				}
				insertNNode(alist, na, nb);
				nb = bnext;
			}
		}
	}

	if (nb) {
		if (lastna) {
			nb->prev = lastna;
			lastna->next = nb;
		} else {
			nb->prev = NULL;
			*alist = nb;
		}
		while (nb) {
			bnext = nb->next;
			cb = nb->id;
			if (cb == a) {
				a_in_b = nb;
			} else {
				sidb = nb->sid;
				sidb->a = ab;
				sidb->b = cb;
				sidb->updated = false;

				lp = nheads + cb;
				tp = *lp;
				ntp = nb->twin;

				ntp->id = ab;
				if (tp != ntp) moveNNode2First(ntp, tp, lp);
			}
			nb = bnext;
		}
	}

	if (na) {
		while (na) {
			anext = na->next;
			ca = na->id;
			if (ca == b) {
				b_in_a = na;
			} else {
				sida = na->sid;
				sida->a = ab;
				sida->b = ca;
				sida->updated = false;

				lp = nheads + ca;
				tp = *lp;
				ntp = na->twin;

				ntp->id = ab;
				if (tp != ntp) moveNNode2First(ntp, tp, lp);
			}			
			na = anext;
		}
	}

	removeNNode(a_in_b, nheads + a); removeNNode(b_in_a, nheads + a);

	nheads[ab] = nheads[a];	nheads[a] = NULL; nheads[b] = NULL;
}

inline void constructNeighborLists(neighborNode* ntable, neighborNode** nheads, int height, int width)
{
	int es1 = 1;
	int es2 = es1 + (width - 1) * height * 2;
	int es3 = es2 + width * (height - 1) * 2;
	int es4 = es3 + (width - 1) * (height - 1) * 2;
	int nn = 0, nh = 0;

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			nh = y * width + x;
			if (x > 0 && y > 0) {
				nn = ((y - 1) * (width - 1) + (x - 1)) * 2 + 1 + es3;
				insertNNode(ntable + nn, nheads + nh);
			}
			if (y > 0) {
				nn = ((y - 1) * width + x) * 2 + 1 + es2;
				insertNNode(ntable + nn, nheads + nh);
			}
			if (x < width - 1 && y > 0) {
				nn = ((y - 1) * (width - 1) + x) * 2 + es4;
				insertNNode(ntable + nn, nheads + nh);
			}
			if (x > 0) {
				nn = (y * (width - 1) + (x - 1)) * 2 + 1 + es1;
				insertNNode(ntable + nn, nheads + nh);
			}
			if (x < width - 1) {
				nn = (y * (width - 1) + x) * 2 + es1;
				insertNNode(ntable + nn, nheads + nh);
			}
			if (x > 0 && y < height - 1) {
				nn = (y * (width - 1) + (x - 1)) * 2 + 1 + es4;
				insertNNode(ntable + nn, nheads + nh);
			}
			if (y < height - 1) {
				nn = (y * width + x) * 2 + es2;
				insertNNode(ntable + nn, nheads + nh);
			}
			if (x < width - 1 && y < height - 1) {
				nn = (y * (width - 1) + x) * 2 + es3;
				insertNNode(ntable + nn, nheads + nh);
			}
		}
	}
}

int segment_lep(int amount_v, int amount_e, int *edge_a, int *edge_b, double *edge_w,
				  int* rgblabel, int ck, int* lbplabel, int tk, int height, int width, 
				  double* angleWeights, double* ep_prior, double* thrs, int thrn, int* segrlt)
{
	universe*		u		= new universe(amount_v);
	strengthNode *	stable	= new strengthNode[amount_e + 1];
	neighborNode *	ntable	= new neighborNode[2 * amount_e + 1];
	strengthNode**	sheads	= (strengthNode**)calloc(MAXSTRENGTH * STABLE_SCALE, sizeof(strengthNode*));
	neighborNode**	nheads	= (neighborNode**)calloc(amount_v * 2, sizeof(neighborNode*));

	rgb *colors = new rgb[amount_v * 2];
#if SAVESEGIMG
#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int i = 0; i < amount_v * 2; i++) colors[i] = random_rgb();
#endif

	int * rgbcom = (int*)calloc(amount_v * ck * 2, sizeof(int));
	int * lbpcom = (int*)calloc(amount_v * tk * 2, sizeof(int));

	for (int i = 0; i < amount_v; i++) {
		u->setenergy(i, 0); u->setenergy_t(i, 0);
		rgbcom[i * ck + rgblabel[i]] = 1; lbpcom[i * tk + lbplabel[i]] = 1;
	}

#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int i = 0; i < amount_e; i++) {
		int a = edge_a[i], b = edge_b[i];
		double w0 = 0, w1 = 0, w2 = edge_w[i], sk = 0;
		double w0_t = 0, w1_t = 0, sk_t = 0, pb = 0;

		int nn = 1 + i * 2, sn = 1 + i;

		if (rgblabel[a] != rgblabel[b]) {w1 = 2 * log(2); pb += 2 * log(2);}
		if (lbplabel[a] != lbplabel[b]) {w1_t = 2 * log(2); pb += 2 * log(2);}

		pb *= normalise(w2); pb += w2;
		
		ntable[nn].id = b;
		ntable[nn].sid = stable + sn;
		ntable[nn].twin = ntable + nn + 1;
		ntable[nn].next = NULL;
		ntable[nn].prev = NULL;

		ntable[nn + 1].id = a;
		ntable[nn + 1].sid = stable + sn;
		ntable[nn + 1].twin = ntable + nn;
		ntable[nn + 1].next = NULL;
		ntable[nn + 1].prev = NULL;

		stable[sn].a = a;
		stable[sn].b = b;
		stable[sn].engy_c = w1;
		stable[sn].engy_t = w1_t;
		stable[sn].nm = 1;
		stable[sn].nv = w2;
		stable[sn].pb = pb;
		stable[sn].updated = true;
		stable[sn].valid = true;
		stable[sn].prev = 0;
		stable[sn].next = 0;
		stable[sn].pb_entry = stable_entry(pb);
		stable[sn].angleList = NULL;
		memset(stable[sn].angleHist, 0, sizeof(int) * 16);
		stable[sn].effort = 0;
		stable[sn].borderlen = 0;
		
		stable[sn].ep_prior = ep_prior[i];
	}
	for (int i = 0; i < amount_e; i++) insertSNode(stable + 1 + i, sheads + stable[1 + i].pb_entry);

	gridNode*		gtable	= new gridNode[(height + 1) * (width + 1)];
	angleNode*		atable	= new angleNode[((height - 1) * width + (width - 1) * height) * 2 + 1];
	int es1 = 0, es2 = es1 + (width - 1) * height, ec1 = 0, ec2 = 0;
#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int y = 0; y < height + 1; y++) {
		for (int x = 0; x < width + 1; x++) {
			int ind = y * (width + 1) + x;
			gtable[ind].id = y * (width + 1) + x;
			gtable[ind].x = x; gtable[ind].y = y;
			gtable[ind].dssEnds[0] = y * (width + 1) + (x - 1);
			gtable[ind].dssEnds[1] = (y - 1) * (width + 1) + x;
			gtable[ind].dssEnds[2] = y * (width + 1) + (x + 1);
			gtable[ind].dssEnds[3] = (y + 1) * (width + 1) + x;
		}
	}

	int ind = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (x < width-1) {
				int ap = y * (width + 1) + (x + 1), bp = (y + 1) * (width + 1) + (x + 1);
				atable[ind].id = ap;
				atable[ind].dir = 3;
				atable[ind].invdir = 1;
				atable[ind].next = NULL;
				atable[ind].prev = atable + ind + 1;
				
				atable[ind + 1].id = bp;
				atable[ind + 1].dir = 1;
				atable[ind + 1].invdir = 3;
				atable[ind + 1].next = atable + ind;
				atable[ind + 1].prev = NULL;

				stable[es1 + ec1 + 1].angleList = atable + ind + 1;
				memset(stable[es1 + ec1 + 1].angleHist, 0, sizeof(int) * 16);
				stable[es1 + ec1 + 1].angleHist[15] = 2;
				stable[es1 + ec1 + 1].borderlen = 0;
				stable[es1 + ec1 + 1].effort = 2 * angleWeights[15] + angleWeights[16];
				ec1++; ind += 2;
			}

			if (y < height-1) {
				int ap = (y + 1) * (width + 1) + x, bp = (y + 1) * (width + 1) + (x + 1);
				atable[ind].id = ap;
				atable[ind].dir = 2;
				atable[ind].invdir = 0;
				atable[ind].next = NULL;
				atable[ind].prev = atable + ind + 1;
				
				atable[ind + 1].id = bp;
				atable[ind + 1].dir = 0;
				atable[ind + 1].invdir = 2;
				atable[ind + 1].next = atable + ind;
				atable[ind + 1].prev = NULL;

				stable[es2 + ec2 + 1].angleList = atable + ind + 1;
				memset(stable[es2 + ec2 + 1].angleHist, 0, sizeof(int) * 16);
				stable[es2 + ec2 + 1].angleHist[15] = 2;
				stable[es2 + ec2 + 1].borderlen = 0;
				stable[es2 + ec2 + 1].effort = 2 * angleWeights[15] + angleWeights[16];
				ec2++; ind += 2;
			}
		}
	}

	constructNeighborLists(ntable, nheads, height, width);
	double curthr = thrs[0]; int segind = 0;
	for (int i = 0; i < MAXSTRENGTH * STABLE_SCALE; i++) {
		strengthNode* si = *(sheads + i);
		if (stable_value(i) > curthr) {
			saveUinImg(u, segrlt + amount_v * segind, colors, height, width, "./segimg", curthr);
			if (segind == thrn - 1) break;
			curthr = thrs[segind++];
		}
		if (si) {
			while (si) {
				strengthNode* sinext = si->next;
				if (si->valid) {
					if (si->updated) {
						merge2regions(u, stable, sheads, ntable, nheads, gtable, atable, 
								rgbcom, lbpcom, ck, tk, angleWeights, si);
					} else {
						int a = si->a, b = si->b;
						int nm = si->nm;
						double nv = si->nv;
						int sizea = u->size(a), sizeb = u->size(b);
						double w0 = 0, w1 = 0;
						double w0_t = 0, w1_t = 0;

						w0 = u->getenergy(a) + u->getenergy(b);
						w0_t = u->getenergy_t(a) + u->getenergy_t(b);
						double tlogt = xlogxs[sizea + sizeb];
						int *apos = (int*)(rgbcom + a * ck);
						int *bpos = (int*)(rgbcom + b * ck);

						for (int ti = 0; ti < ck; ti++) w1 -= xlogxs[apos[ti] + bpos[ti]];
						w1 += tlogt;
						apos = (int*)(lbpcom + a * tk);
						bpos = (int*)(lbpcom + b * tk);

						for (int ti = 0; ti < tk; ti++) w1_t -= xlogxs[apos[ti] + bpos[ti]];
						w1_t += tlogt;

						double pb = strength(w0, w1, w0_t, w1_t, nv, nm, si->borderlen, si->effort, si->ep_prior);

						si->engy_c = w1; si->engy_t = w1_t;
						if (stable_entry(pb) <= i) {
							merge2regions(u, stable, sheads, ntable, nheads, gtable, atable, 
								rgbcom, lbpcom, ck, tk, angleWeights, si);
						} else {
							si->pb = pb;
							si->updated = true;
							si->pb_entry = stable_entry(pb);
							insertSNode(si, sheads + stable_entry(pb));
						}
					}
				}
				si = sinext;
			}
		}
	}

	delete		u;
	delete[]	colors;
	delete[]	stable;
	delete[]	sheads;
	delete[]	ntable;
	delete[]	nheads;
	delete[]	gtable;
	delete[]	atable;
	free(rgbcom); free(lbpcom);

	return 0;
}

#endif
