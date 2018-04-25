#ifndef _LOKI_DUMP_H_
#define _LOKI_DUMP_H_

extern int read_dump(int *,int *,int *,long *,int *,int,struct loki *);
extern void write_dump(int,int,int,long,int,int,struct loki *);
extern void write_xml_dump(int,int,int,long,int,int,struct loki *);

#endif
