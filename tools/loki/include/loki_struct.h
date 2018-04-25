#ifndef _LOKI_STRUCT_H_
#define _LOKI_STRUCT_H_

#ifndef _UTILS_H_
#include <utils.h>
#endif

/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *             Simon Heath - University of Washington                       *
 *                                                                          *
 *                       July 1997                                          *
 *                                                                          *
 * loki_struct.h:                                                           *
 *                                                                          *
 ****************************************************************************/

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#define ST_DATA 1
#define ST_ID 2
#define ST_SIRE 4
#define ST_DAM 8
#define ST_MARKER 0x10
#define ST_HAPLO 0x20
#define ST_LINKED 0x40
#define ST_FACTOR 0x80
#define ST_MODEL 0x100
#define ST_TRAITLOCUS 0x200
#define ST_RANDOM 0x400
#define ST_TRAIT 0x800
#define ST_RESTRICT 0x1000
#define ST_REQUIRED 0x4000
#define ST_STRING 0x8000
#define ST_REAL 0x10000
#define ST_INTEGER 0x20000
#define ST_REALTYPE 0x40000
#define ST_INTTYPE 0x80000
#define ST_ARRAY 0x100000
#define ST_SCALAR 0x200000
#define ST_CONSTANT 0x400000
#define ST_FLAG 0x800000
#define ST_CENSORED 0x1000000
#define ST_GROUP 0x2000000
#define ST_SEX 0x4000000
#define ST_FAMILY 0x8000000
#define ST_MULTIPLE 0x10000000
#define ST_LUMPED 0x20000000
#define ST_NOT_REALLY_REQUIRED 0x40000000
#define ST_PED (ST_ID|ST_SIRE|ST_DAM|ST_FAMILY)

#define LINK_AUTO 0
#define LINK_X 1
#define LINK_Y 2
#define LINK_MIT 3
#define LINK_Z 4
#define LINK_W 5
#define UNLINKED 6
#define N_LINK_TYPES 7
#define LINK_TYPES_MASK 7

#define LINK_MIRRORED 16 /* This linkage group will be mirrored */
#define LINK_PSEUDO 32 /* Pseudo chromosome (mirror of real chromosome) */

#define MAP_HALDANE 0
#define MAP_KOSAMBI 1

#define DEFAULT_ANALYSIS 0
#define AFFECTED_ANALYSIS 1
#define NULL_ANALYSIS 2
#define IBD_ANALYSIS 4
#define ESTIMATE_IBD 8

#define OUTPUT_TYPE_ORIGINAL 0
#define OUTPUT_VERSION_2_1 1
#define OUTPUT_VERSION_2_2 2
#define OUTPUT_VERSION_2_3 3
#define DEFAULT_OUTPUT_TYPE OUTPUT_VERSION_2_3

#define REMSIZE 256

#define NULL_STR "<NULL>"

struct remember {
     struct remember *next;
     void *mem[REMSIZE];
     int pos;
};

extern struct remember *AddRemem(void *p,struct remember *rblock);
extern void FreeRemem(struct remember *rblock);

#endif
