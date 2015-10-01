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
 *   Read one single amplitude from file                                   
 *   last update 05/12/00
 *
 *  T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

/*
different data formats of output:
format=1  :  SU (IEEE)
format=2  :  ASCII
format=3  :  BINARY (IEEE)
*/


float readdsk(FILE *fp_in, int format){
	float amp;



	switch(format){
                case 1 : /* SU*/ 
                        err(" Sorry, SU-format for snapshots not implemented yet. \n");
                        break;
		case 2 :  /*ASCII*/
                        fscanf(fp_in,"%e\n", &amp); 
                        break;
		case 3 :   /* BINARY */

			fread(&amp, sizeof(float), 1, fp_in);
              		break;
	                
		default :
			printf(" Don't know the format for the snapshot-data !\n");
			err(" No output was written. ");
	}

	return amp;
}
