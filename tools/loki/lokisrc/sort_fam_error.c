/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                     Simon Heath - CNG, Evry                              *
 *                                                                          *
 *                            May 2003                                      *
 *                                                                          *
 * sort_fam_error.c:                                                        *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2003                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "utils.h"
#include "loki.h"
#include "loki_utils.h"
#include "ped_utils.h"
#include "gen_elim.h"

static nuc_fam *fam;
static int *errs;

static int cmp_err(const void *s1,const void *s2)
{
  return fam[*(int *)s2].n_err-fam[*(int *)s1].n_err;
}

static int cmp_err1(const void *s1,const void *s2)
{
  return errs[*(int *)s2]-errs[*(int *)s1];
}

void sort_fam_error(int *perm,int n,struct loki *loki)
{
  fam=loki->family->families;
  gnu_qsort(perm,n,sizeof(int),cmp_err);
}

void sort_ind_error(int *perm,int n,int *err_lst)
{
  errs=err_lst;
  gnu_qsort(perm,n,sizeof(int),cmp_err1);
}
