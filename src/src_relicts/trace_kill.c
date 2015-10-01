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
 *   Apply trace killing
	15.12.11, S. Heider
 *  ----------------------------------------------------------------------*/
#include "fd.h"

void trace_kill(float **section, float *kill_vector, int nsrc_glob, int ns, int ntr, int ntr_glob, int ishot){

/* declaration of variables */
extern int REC1, REC2, MYID;
float ** kill_tmp;
int i, j, h;
int READ_KILL_TRACES;
char trace_kill_file[STRING_SIZE];

FILE *ftracekill;


READ_KILL_TRACES=1; /* read matrix from file trace_kill.dat; else define kill_tmp matrix on your self*/

if(READ_KILL_TRACES==1){

kill_tmp = matrix(1,ntr_glob,1,nsrc_glob);

sprintf(trace_kill_file,"./trace_kill/trace_kill.dat");
ftracekill=fopen(trace_kill_file,"r");

/* Some error checking */
if (ftracekill == NULL) err(" trace_kill.dat file could not be opened !");

for(i=1;i<=ntr_glob;i++){
	for(j=1;j<=nsrc_glob;j++){
		fscanf(ftracekill,"%f",&kill_tmp[i][j]);
	}
}

fclose(ftracekill);

/* distribute vector on CPUs */
int h=1;
for(i=REC1;i<=REC2;i++){
printf("kill_tmp[%i]\n",kill_tmp[i][ishot]);
if (kill_tmp[i][ishot] != 1){
if (kill_tmp[i][ishot] != 0){
	err(" Error! In trace_kill.dat only 0 and 1 are allowed!"); }}

   kill_vector[h] = kill_tmp[i][ishot];
   h++;
}

for(i=1;i<=ntr;i++){
	if(kill_vector[i]==1){
	      for(j=1;j<=ns;j++){
     	          section[i][j] = 0.0;
              }  
	}
}

free_matrix(kill_tmp,1,ntr_glob,1,nsrc_glob);

} /* end of if(READ_KILL_TRACES==1)*/

} /* end of function */
