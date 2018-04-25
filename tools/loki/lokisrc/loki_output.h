#ifndef _LOKI_OUTPUT_H_
#define _LOKI_OUTPUT_H_

extern void OutputSample(FILE *,int,struct loki *);
extern void OutputFreqHeader(FILE *,struct loki *);
extern void OutputFreq(FILE *,int,struct loki *);
extern void OutputHeader(FILE *,struct loki *);
extern void Output_BV(FILE *,struct loki *);

#endif
