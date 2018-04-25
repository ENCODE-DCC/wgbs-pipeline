/*
 *  ibs_check.h
 *  loki
 *
 *  Created by Simon Heath on 10/13/05.
 *  Copyright 2005 Simon C. Heath.
 *  This is free software.  You can distribute it and/or modify it  
 *  under the terms of the Modified BSD license, see the file COPYING
 *
 */

struct result {
  double ibs_obs;
  double ibs_exp;
  double ibs_sd;
  int id,id1;
  int rel;
  int n,n1,n2;
  int n_disc;
};

void ibs_check(int,struct loki *);
