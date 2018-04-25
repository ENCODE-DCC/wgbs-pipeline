#ifndef _MEIOSIS_SCAN_H_
#define _MEIOSIS_SCAN_H_

#define MSCAN_INDIVIDUAL .1
#define MSCAN_HS_FAMILY .4
#define MSCAN_FS_FAMILY .25
#define MSCAN_GP_FAMILY .25 

extern void meiosis_scan(int,const struct loki *);
extern int set_mscan_probs(double *);

#define seg_pen1(loc,loc1,a,b,c,d) (loc==loc1?gen_pen(loc,a,b,c,d):seg_pen(loc,a,b,c,d))

#endif
