#ifndef _LOCUS_H_
#define _LOCUS_H_

void alloc_marker(int,struct loki *);
void free_marker_data(struct Marker *);
void free_marker(struct Marker *,int);
void free_locus(struct Locus *);
void Init_Markers(int,struct loki *);
void init_marker_segs(struct Marker *,struct loki *);

#endif
