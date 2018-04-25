#ifndef _SAMPLE_SEGS_H_
#define _SAMPLE_SEGS_H_

struct family {
	int ids,idd;
	int nkids;
	int *kids;
};

struct node {
	int type;
	int id;
	int locus;
};

struct sg_off {
	struct sg_off *next;
	int col;
};

#define MAT_SEG_IRR 0x20
#define PAT_SEG_IRR 0x40
#define MAT_SEG_FIXED_M 0x80
#define MAT_SEG_FIXED_P 0x100
#define PAT_SEG_FIXED_M 0x200
#define PAT_SEG_FIXED_P 0x400
#define NODE_MAT_SEG 0x800
#define NODE_PAT_SEG 0x1000
#define NODE_ORDER 0x2000



#endif
