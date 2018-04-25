#include <sys/types.h>

uint32_t hashword(const uint32_t *,size_t,uint32_t);
void hashword2(const uint32_t *,size_t,uint32_t *,uint32_t *);
uint32_t hashlittle(const void *,size_t,uint32_t);
void hashlittle2(const void *,size_t,uint32_t *,uint32_t *);
uint32_t hashbig(const void *,size_t,uint32_t);
