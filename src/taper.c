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
 *   Write seismograms to disk                                  
 *   last update 19/01/02, T. Bohlen
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  taper(float **sectionpdiff, int ntr, int ns){

	/* declaration of extern variables */
	extern int MYID, TAPERLENGTH;
	extern FILE *FP;
	
	/* declaration of local variables */
	int i,j, h;
	int tracl1;
	float a;
	float damping, amp;
	float *window=NULL, *amp1=NULL;
	
	window = vector(1,ns);
        amp1 = vector(1,ns);
	
	/* "Cerjan"-Window */
        damping=100.0;
        amp=1.0-damping/100.0;
        a=sqrt(-log(amp)/((TAPERLENGTH-1)*(TAPERLENGTH-1)));
        for (i=1;i<=ns;i++){
		window[i]=1.0;
		amp1[i]=0.0;
	}
	if (MYID==0){
		fprintf(FP,"\n ntr: %d\n",ntr);
		fprintf(FP,"\n ns: %d\n",ns);
		fprintf(FP,"\n TAPERLENGTH: %d\n",TAPERLENGTH);
	}
        for (i=1;i<=TAPERLENGTH;i++){amp1[i]=exp(-(a*a*(TAPERLENGTH-i)*(TAPERLENGTH-i)));}
        for (i=1;i<=1;i++){window[i]=amp1[i];}
        h=1;
        for (i=ns;i>=(ns-TAPERLENGTH+3);i--){window[i]=amp1[h];h++;}
	
	/* Apply taper window */
	for(tracl1=1;tracl1<=ntr;tracl1++){ 
		
		for(j=1;j<=ns;j++){
			sectionpdiff[tracl1][j]*=window[j];
		}
	}
	free_vector(window,1,ns);
	free_vector(amp1,1,ns);
}
