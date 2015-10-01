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

/**/
/*------------------------------------------------------------------------
 *   Exchange FD-Parameters between PEs                         
 *   last update 29/06/2002
 *
 *  T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void exchange_mod_es(float ** matmod, int ncptot, int nparameter){


	/* definition of local variables */
	extern int MYID;
	int ntotal = ncptot * nparameter;
	int i,j,h; 
	float fdum1[ntotal];
	
	printf("MYID = %d \t ntotal = %d \n",MYID,ntotal);
	 
	if (MYID == 0){ 
        
	h=1;
	for(i=1;i<=ncptot;i++){
	    for(j=1;j<=nparameter;j++){
	        fdum1[h] = matmod[i][j];
		h++;
	    }
	}                                                                 
        
	} /** if (MYID == 0) **/
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&fdum1,ntotal,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	h=1;
	for(i=1;i<=ncptot;i++){
	    for(j=1;j<=nparameter;j++){
	        matmod[i][j] = fdum1[h] ;
		h++;
	    }
	}                                                                 
        
	
}
