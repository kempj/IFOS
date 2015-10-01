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
 *   calculate step length for material parameter update
 *   
 *   waveconv = conjugated gradient direction
 *   gradp    = preconditionated gradient
 *   
 *   Daniel Koehn
 *   last update 9.11.2007
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
float calc_opt_step_test(float *  L2t, float ** waveconv, float ** gradg, float * epst, int sws, float C_vp){

extern int NX, NY, IDX, IDY, MYID;
extern float EPSILON, EPSILON_u, EPSILON_rho;
int i, j, n;
float opteps, H1, H1sum, H2, H2sum, etas, H1max, H1maxsum;


/* calculate optimal step size after Tarantola (1986)*/
/*H1max=0.0;
H1=0.0;
H2=0.0;
for (i=1;i<=NX;i=i+IDX){
     for (j=1;j<=NY;j=j+IDY){
	 H1 += 1e+39 * waveconv[j][i]*gradg[j][i];
	 H2 += waveconv[j][i]*(1.0/C_vp)*waveconv[j][i];
	 
	 if(fabs(H1)>H1max){H1max=H1;}
	  
     }
  }

H1sum=0.0;
MPI_Allreduce(&H1,&H1sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

H2sum=0.0;
MPI_Allreduce(&H2,&H2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

H1maxsum=0.0;
MPI_Allreduce(&H1max,&H1maxsum,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);*/


etas = epst[sws];

/* calculate optimal step length for Vp update */
opteps =   etas  * L2t[sws]; 

if((opteps/etas) > 2.0){opteps = 2.0 * etas;}

printf("MYID = %d \t L2t[sws] = %e \t epst[sws] = %e \t sws = %d \t opteps = %e \t C_vp = %e \n ",MYID,L2t[sws],epst[sws],sws,opteps,C_vp);

return opteps;		
}

