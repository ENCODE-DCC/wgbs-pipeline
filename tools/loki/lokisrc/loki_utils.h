#ifndef _LOKI_UTILS_H_
#define _LOKI_UTILS_H_

#ifndef _LOKI_PEEL_H
#include "loki_peel.h"
#endif
#ifndef _PED_UTILS_H_
#include "ped_utils.h"
#endif
#ifndef _STRING_UTILS_H_
#include "string_utils.h"
#endif

void init_utils(const struct loki *);
void print_marker_name(FILE *,const int);
size_t print_orig_id(FILE *,int);
size_t print_orig_id1(FILE *,int);
size_t print_orig_family(FILE *,const int,const int);
void print_orig_triple(FILE *,const int);
void print_orig_allele_id(FILE *,const int);
size_t get_max_idlen(void);
void print_allele_name(FILE *,const int,const int,const int);
void print_allele_type1(FILE *,const int,const int);
void print_allele_type(FILE *,const int,const int,const int);
int cmp_loci(const void *,const void *);
void get_locuslist(struct Locus **,const int,int *,int);
void print_orig_gtypes(FILE *,struct Id_Record *,struct Marker *);
void print_family_gtypes(FILE *,nuc_fam *,struct Marker *,int *,gtype **);
void print_derived_gtypes(FILE *,struct Id_Record *,struct Marker *,int *,gtype **);
void print_allele_type2(FILE *,struct Marker *,int,int);
void print_peel_elem(FILE *,void *,int,int);
void print_peelseq(FILE *,struct Peelseq_Head *);
string *str_print_orig_family(string *,const int);
string *str_print_orig_id(string *,int);

#endif
