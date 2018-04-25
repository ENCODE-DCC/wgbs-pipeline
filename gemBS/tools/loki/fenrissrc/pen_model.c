/****************************************************************************
 *                                                                          *
 *     Loki - Programs for genetic analysis of complex traits using MCMC    *
 *                                                                          *
 *                    Simon Heath - CNG, France                             *
 *                                                                          *
 *                         August 2002                                      *
 *                                                                          *
 * pen_model.c:                                                             *
 *                                                                          *
 * Copyright (C) Simon C. Heath 2002                                        *
 * This is free software.  You can distribute it and/or modify it           *
 * under the terms of the Modified BSD license, see the file COPYING        *
 *                                                                          *
 ****************************************************************************/

#include <stdio.h>

#include "utils.h"
#include "loki.h"
#include "loki_peel.h"
#include "fenris_peel.h"
#include "fenris.h"

fenris_pen_func *get_pen_model(int pen_type)
{
	fenris_pen_func *pen=0;
	
	switch(pen_type) {
	 case PEN_MODEL_EQUAL:
		pen=fpen_emodel_equal;
		break;
	 case PEN_MODEL_PROP:
		pen=fpen_emodel_prop;
		break;
	 case PEN_MODEL_EMP:
		pen=fpen_emodel_emp;
		break;
	 default:
		ABT_FUNC("Unknown penetrance model\n");
	}
	return pen;
}
