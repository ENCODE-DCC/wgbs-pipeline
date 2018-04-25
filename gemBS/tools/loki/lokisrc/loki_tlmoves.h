#ifndef _LOKI_TLMOVES_H_
#define _LOKI_TLMOVES_H_

/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - Rockefeller University                         *
 *                                                                          *
 *                       November1997                                       *
 *                                                                          *
 * loki_tlmoves.h:                                                          *
 *                                                                          *
 * Defines for move probabilities                                           *
 *                                                                          *
 ****************************************************************************/

#define BIGMOVE_PROB .25
#define SMALLMOVE_P 0.5

#define BETA_A 5.75
#define BETA_B 1.25

#define BIRTH_STEP 0.5
#define DEATH_STEP 0.5

#define PROP_RATIO 0.1

void Sample_LinkageGroup(const int,struct loki *);
void Sample_TL_Position(const int,struct loki *);
void TL_Birth_Death(struct loki *);
void TL_Alloc(const struct loki *);
int get_tl_position(double *,const struct loki *);
void Flip_TL_Alleles(const int,struct loki *loki);
void Flip_TL_Mode(const int,struct loki *loki);
struct Locus **get_sorted_locuslist(const int,int *,int);
double safe_exp(double);
double calc_tl_like(struct Locus *,const int,struct loki *);

#endif
