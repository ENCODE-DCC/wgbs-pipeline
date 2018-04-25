/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         August 2002                                      *
 *                                                                          *
 * count_bits.c:                                                            *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2002                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "fenris.h"

int count_bits(int *nb)
{
	int i,j,k,n,comp,cs,gkid;
	int *kids,nk,spouse,sp1,kid,*hap,fx=0;

	for(i=0;i<ped_size;i++) id_array[i].flag=(id_array[i].data || id_array[i].data1)?1:0;
	for(j=0;j<n_markers;j++) {
		hap=marker[j].haplo;
		for(i=0;i<ped_size;i++) if(hap[i]) id_array[i].flag|=2;
	}
	for(i=comp=0;comp<n_comp-singleton_flag;comp++) {
		cs=comp_size[comp];
		for(n=j=0;j<cs;j++,i++) {
			if(id_array[i].sire) n+=2;
			else {
				n--;
				if(!sex_map && !(id_array[i].flag&4)) {
					/* Is this founder a member of a symmetric founder couple? */
					id_array[i].flag|=4;
					kids=id_array[i].kids;
					nk=id_array[i].nkids;
					kid=kids[0];
					spouse=id_array[i].sex==2?id_array[kid].sire-1:id_array[kid].dam-1;
					id_array[spouse].flag|=4;
					if(!id_array[spouse].sire && id_array[spouse].nkids==nk) {
						for(k=1;k<nk;k++) {
							kid=kids[k];
							sp1=id_array[i].sex==2?id_array[kid].sire-1:id_array[kid].dam-1;
							if(sp1!=spouse) break;
						}
						if(k==nk && !((id_array[i].flag|id_array[spouse].flag)&3)) {
							gkid=-1;
							for(k=0;k<nk;k++) {
								kid=kids[k];
								if(id_array[kid].nkids) {
									gkid=id_array[kid].kids[0];
									id_array[gkid].flag|=(id_array[kid].sex==1?8:16);
									break;
								}
							}
							if(gkid>=0) {
								n--;
								fx++;
								fputs("Reflection with untyped founder couple ",stdout);
								print_orig_id(stdout,i+1);
								fputc(' ',stdout);
								print_orig_id(stdout,spouse+1);
								printf(" - fixing %s seg of grandchild ",(id_array[gkid].flag&8)?"paternal":"maternal");
								print_orig_id(stdout,gkid+1);
								fputc('\n',stdout);
							}
						}
					}
				}
			}
		}
		nb[comp]=n;
	}
	if(singleton_flag) nb[comp]=0;
	return fx;
}
