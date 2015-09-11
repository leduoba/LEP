/*************************************************************************************
Copyright (C) 2015 Qiyang Zhao, Beihang Univ.

This program is free software; you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation; either 
version 2 of the License, or (at your option) any later version.

Please cite the following paper if you want to use the source codes in your work:

Qiyang Zhao, Segmenting natural images with the least effort as humans, BMVC 2015.
**************************************************************************************/
#ifndef TABLESTRUCT_H_
#define TABLESTRUCT_H_

#include "disjoint-set.h"

using namespace std;

////////////////////////////////////////// user effort
struct ANode{
	int					id;				// grid points
	int					dir;			// �ǵ���ʹ�õ�gridPoint�ķ���
	int					invdir;			// �ǵ�ĶԶ���ʹ�õ�gridPoint�ķ���
	struct ANode*		next;			// the next dssNode in the dssList
	struct ANode*		prev;			// the previous dssNode in the dssList
};

struct GNode{
	int					id;				// y * width + x
	int					x;
	int					y;
	int					dssEnds[4];		// ��0��1��2��3
};

typedef	struct ANode angleNode;
typedef struct GNode gridNode;

struct NNode{
	int					id;			// �ڽ�����ı��
	struct NNode*		twin;		// ǰһ����¼
	struct NNode*		prev;		// ǰһ����¼
	struct NNode*		next;		// ��һ����¼
	struct SNode*		sid;
};

struct SNode{
	int									a;
	int									b;
	int									nm;
	int									pb_entry;

	double								pb;			//α��������
	double								engy_c;		//��ɫ���������ϲ���
	double								engy_t;		//�������������ϲ���
	double								nv;

	struct SNode*						next;
	struct SNode*						prev;
	bool								valid;
	bool								updated;

	struct ANode*						angleList;
	int									angleHist[16];
	int									borderlen;
	double								effort;

	double								ep_prior;
};

typedef struct NNode neighborNode;
typedef struct SNode strengthNode;

inline void insertSNode(strengthNode* sn, strengthNode** sh)
{
	if (*sh) (*sh)->prev = sn;
	sn->next = *sh; sn->prev = NULL;
	*sh = sn;
}

inline void removeNNode(neighborNode* nn, neighborNode** nh)
{
	if (nn->prev) nn->prev->next = nn->next;
	if (nn->next) nn->next->prev = nn->prev;
	
	if (!nn->prev) *nh = nn->next;
}

inline void removeNNode(neighborNode* nn)
{
	// nn���ڵ��ڽ������еĽڵ�������Ȼ��С��2����nn��Ȼ���ǵ�һ���ڵ�
	nn->prev->next = nn->next;
	if (nn->next) nn->next->prev = nn->prev;
}

inline void insertNNode(neighborNode* nn, neighborNode** nh)
{
	nn->next = *nh;
	if (*nh) (*nh)->prev = nn;
	nn->prev = NULL;
	*nh = nn;
}

inline void moveNNode2First(neighborNode* nn, neighborNode* fn, neighborNode** nh)
{
	nn->prev->next = nn->next;
	if (nn->next) nn->next->prev = nn->prev;
	
	nn->next = fn;
	fn->prev = nn;
	nn->prev = NULL;
	*nh = nn;
}

inline void insertNNode(neighborNode** a, neighborNode* na, neighborNode* nb)
{
	// insert nb just before na
	nb->next = na;
	nb->prev = na->prev;
	if (na->prev) {
		na->prev->next = nb;
	} else {
		*a = nb;
	}
	na->prev = nb;
}

inline void removeANode(angleNode* an, strengthNode* si)
{
	if (an->prev) an->prev->next = an->next;
	if (an->next) an->next->prev = an->prev;
	
	if (!an->prev) si->angleList = an->next;
}

inline void insertANode(angleNode* an, angleNode* bn, strengthNode* si)
{
	// insert bn just before an
	bn->next = an;
	bn->prev = an->prev;
	if (an->prev) {
		an->prev->next = bn;
	} else {
		si->angleList = bn;
	}
	an->prev = bn;
}

#endif
