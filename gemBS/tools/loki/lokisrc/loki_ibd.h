#ifndef _LOKI_IBD_H_
#define _LOKI_IBD_H_

void set_n_ibd_rec(int);
int get_n_ibd_rec(void);
int SetupIBD(const struct loki *);
void Handle_IBD(const struct loki *);
double score_ibd(int,int *,int,int,int *,double *,const struct loki *);
void get_founders(unsigned long **,int **,int **,int **,const struct loki *);
void Output_Sample_IBD(int,int,const struct loki *);
void Output_Merlin_IBD(int,const struct loki *);
void Output_Solar_IBD(int,int *,const struct loki *);
void get_founder_params(unsigned long **,int **,int **,int **,const struct loki *);
void sample_segs(const struct loki *);
int write_ibd_dump(FILE *,int,const struct loki *);
int read_ibd_dump(FILE *,int *,char *,const struct loki *);

#define IBD_MIN_GRID_STEP .0001
#define IBD_MAX_GRID 10000

#define IBD_EST_DISCRETE 1
#define IBD_EST_MARKERS 2
#define IBD_EST_GRID 3

#define DEFAULT_IBD_MODE 0
#define MERLIN_IBD_MODE 1
#define SOLAR_IBD_MODE 2
#define IBD_MODE_MASK 3
#define IBD_SINGLE_POINT 4
#define COMPRESS_IBD 8

#endif
