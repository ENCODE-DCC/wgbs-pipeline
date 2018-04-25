#ifndef _COUNT_DBR_H_
#define _COUNT_DBR_H_

int init_dbr_shm(void);
int init_dbr_count(void);
void free_dbr_count(void);
void zero_dbr_count(void);
void count_dbr(void);
void write_dbr_count(int);

#endif
