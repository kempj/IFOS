/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
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
 *   Write note to stdout                          
 *   last update 16/02/02
 *
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void note(FILE *fp){

extern char LOG_FILE[STRING_SIZE];
extern int MYID;

	fprintf(fp," Please note: \n");
	fprintf(fp," Each processing element (PE) is writing log information during program \n");
	fprintf(fp," execution to %s.PE .\n",LOG_FILE);
	fprintf(fp," See corresponding log-files for further information on program status.\n");
	fprintf(fp," Information about overall program execution \n");
	fprintf(fp," (numerical artefacts, accuracy, computing times etc) \n");
	fprintf(fp," will be written by PE 0 to ");
}
