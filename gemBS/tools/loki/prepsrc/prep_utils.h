#ifndef _PREP_UTILS_H_
#define _PREP_UTILS_H_

#define DataFileError(a) New_DFE(__FILE__,__LINE__,a)

extern void New_DFE(const char *,const int,const char *);
void print_orig_allele_id(FILE *,const int);
void add_to_list(const int,int *,int *);
int find_id_code(char *,int,int);
void cat_file(FILE *,FILE *,char *);
void blank_fam(const int,int *,const int);


#endif
