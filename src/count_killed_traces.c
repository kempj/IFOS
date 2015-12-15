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
 *   Calculate number of killed traces                                  
 *   last update 11/07/12, L. Groos
 *  ----------------------------------------------------------------------*/
#include "fd.h"

void count_killed_traces(int ntr, int swstestshot, int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot, int* ptr_killed_traces, int* ptr_killed_traces_testshots){

/* declaration of variables */
extern int USE_WORKFLOW, WORKFLOW_STAGE;
extern char TRKILL_FILE[STRING_SIZE];
int i,j,h;


/* declaration of variables for trace killing */
int ** kill_tmp, *kill_vector;
char trace_kill_file[STRING_SIZE];	
FILE *ftracekill;
  

kill_tmp = imatrix(1,ntr_glob,1,nsrc_glob);
kill_vector = ivector(1,ntr);

    if(USE_WORKFLOW){
        sprintf(trace_kill_file,"%s_%i.dat",TRKILL_FILE,WORKFLOW_STAGE);
        ftracekill=fopen(trace_kill_file,"r");
        if (ftracekill==NULL){
            sprintf(trace_kill_file,"%s.dat",TRKILL_FILE);
            ftracekill=fopen(trace_kill_file,"r");
            if (ftracekill==NULL){
                err(" Trace kill file could not be opened!");
            }
        }
    }else{
        sprintf(trace_kill_file,"%s.dat",TRKILL_FILE);
        ftracekill=fopen(trace_kill_file,"r");
        if (ftracekill==NULL){
            err(" Trace kill file could not be opened!");
        }
    }

for(i=1;i<=ntr_glob;i++){
	for(j=1;j<=nsrc_glob;j++){
		fscanf(ftracekill,"%d",&kill_tmp[i][j]);
	}
}

fclose(ftracekill);

h=1;
for(i=1;i<=ntr;i++){
   kill_vector[h] = kill_tmp[recpos_loc[3][i]][ishot];
   if (kill_vector[h]==1) *ptr_killed_traces=*ptr_killed_traces+1;
   if ((kill_vector[h]==1) && (swstestshot==1)) *ptr_killed_traces_testshots=*ptr_killed_traces_testshots+1;
   h++;
}
}

