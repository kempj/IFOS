/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
 *
 * This file is part of IFOS.
 * 
 * IFOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.0 of the License only.
 * 
 * IFOS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with IFOS. See file COPYING and/or <http://www.gnu.org/licenses/gpl-2.0.html>.
-----------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
 *   Read source wavelet in su format                                  
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  inseis_source_wavelet(float *section, int ns, int ishot, int SH){

	/* declaration of extern variables */
	extern int MYID;
	extern char SIGNAL_FILE[STRING_SIZE];
    extern char SIGNAL_FILE_SH[STRING_SIZE];
	/* declaration of local variables */
	int j;
	float dump;
	segy tr;
	char data[STRING_SIZE];
    FILE *fpdata;
    
    if(SH==0) {
        sprintf(data,"%s.shot%d",SIGNAL_FILE,ishot);
    } else {
        sprintf(data,"%s.shot%d",SIGNAL_FILE_SH,ishot);
    }
    
	fpdata = fopen(data,"r");
	if (fpdata==NULL) err(" Source wavelet not found ");

	/* SEGY (without file-header) */
	fread(&tr,240,1,fpdata);
	fread(&tr.data,4,ns,fpdata);
			
	for(j=0;j<ns;j++){
	    dump=tr.data[j];
	    section[j+1]=dump;
	}
			  
	fclose(fpdata);
}
