/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2013  For the list of authors, see file AUTHORS.
 *
 * This file is part of DENISE.
 * 
 * DENISE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * DENISE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with DENISE. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   Write program name etc to stdout                          
 *   last update 16/02/02
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void info(FILE *fp){

	fprintf(fp," ***********************************************************\n");
	fprintf(fp," This is program DENISE. Version 1.0                        \n");
	fprintf(fp," Parallel 2-D elastic Finite Difference FWT code            \n");
	fprintf(fp,"                                                            \n");
	fprintf(fp," FWT code written by D. Koehn                               \n");
	fprintf(fp," forward code written by  T. Bohlen                         \n");
	fprintf(fp," Institute of Geosciences, Kiel University, Germany         \n\n");
	fprintf(fp," See COPYING file for copying and redistribution conditions.\n");
	fprintf(fp," ***********************************************************\n");
	fprintf(fp,"\n");

}
