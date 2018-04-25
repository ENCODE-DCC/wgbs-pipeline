#ifndef _GET_PEELSEQ_H_
#define _GET_PEELSEQ_H_

void Get_Peel_Seq(struct loki *);
struct Peelseq_Head *get_peelseq(struct Locus *,struct loki *,int);

typedef struct rf_nd {
	struct rf_nd *next;
	int rf_idx;
} rf_node;

#define DONE_FND_MAT 0x100000
#define DONE_FND_PAT 0x200000
#define PEELED_MAT 0x400000
#define PEELED_PAT 0x800000
#define PEELED_FLAGS (PEELED_MAT|PEELED_PAT)
#define FIXED_MAT 0x1000000
#define FIXED_PAT 0x2000000
#define FIXED_FLAGS (FIXED_MAT|FIXED_PAT)
#define DONE_LINK 0x4000000
#define PAT_LINK 0x8000000
#define MAT_LINK 0x10000000
#define JNT_LINK 0x20000000
#define KID_LINK 0x40000000
#define PAR_LINK (PAT_LINK|MAT_LINK)
#define LINK_FLAGS (FIXED_MAT|FIXED_PAT|DONE_LINK|PAT_LINK|MAT_LINK|JNT_LINK|KID_LINK)

#endif
