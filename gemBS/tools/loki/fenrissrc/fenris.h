struct family {
	int *kids;
	int ids,idd;
	int nkids;	
};

void AllocFenrisStruct(void);
void FenrisSetup(void);
int count_bits(int *);
void calc_pen(int,int,int,int,double *,double *p,int *);

extern int trace_level;

#define TRACE_LEVEL_0 0
#define TRACE_LEVEL_1 1
#define TRACE_LEVEL_2 2
#define TRACE_LEVEL_3 3
#define TRACE_LEVEL_4 4
#define TRACE_MASK 7

#define CHK_TRACE(x) (((trace_level)&TRACE_MASK)>=(x))
