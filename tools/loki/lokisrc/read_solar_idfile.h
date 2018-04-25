#ifndef _READ_SOLAR_IDFILE_H_
#define _READ_SOLAR_IDFILE_H_

#include "bin_tree.h"

void read_solar_idfile(int *,struct loki *);

#define SOLAR_PEDINDEX_FILE "pedindex.cde"

struct rs_data {
	union arg_type *data;
	int id;
};


#endif
