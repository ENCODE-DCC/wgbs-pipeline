#ifndef _SPARSE_H_
#define _SPARSE_H_

/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       March/April 1997                                   *
 *                                                                          *
 * sparse.h:                                                                *
 *                                                                          *
 * Structures for sparse matrix calculation and storage                     *
 *                                                                          *
 ****************************************************************************/

struct SparseMatRec
{
	double val;
	int x;
};

struct Off
{
	struct Off *Next;
	double val;
	int col;
};

struct Diag
{
	struct Off *First;
	double val;
  	int count;
};


#endif
