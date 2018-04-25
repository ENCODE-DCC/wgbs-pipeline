#ifndef _SAMPLE_RAND_H_
#define _SAMPLE_RAND_H_

#define RES_PRIOR_VC0 1.0
#define RES_PRIOR_SC0 1.0
#define RES_PRIOR_VA0 1.0
#define RES_PRIOR_SA0 1.0

void init_rand(struct loki *);
void sample_rand(struct loki *);
void sample_additive_var(struct loki *);

#endif
