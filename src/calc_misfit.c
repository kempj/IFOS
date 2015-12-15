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
 *   Calculate Misfit                                  
 *   last update 18/04/11, L. Rehor
 *  ----------------------------------------------------------------------*/
#include "fd.h"

double calc_misfit(float **sectiondata, float **section, int ntr, int ns, int LNORM, float L2, int itest, int sws, int swstestshot,int ntr_glob, int **recpos_loc, int nsrc_glob, int ishot, int iter){

/* declaration of variables */
extern float DT;
extern int MYID, USE_WORKFLOW, WORKFLOW_STAGE;
extern int TRKILL, NORMALIZE, FC, TIMEWIN;
extern char TRKILL_FILE[STRING_SIZE];
extern int VELOCITY;
int i,j;
float l2;
int umax=0, h;
float abs_section, abs_sectiondata;
float intseis_s, intseis_sd;
float *picked_times=NULL;
float **intseis_section=NULL, **intseis_sectiondata=NULL;
float **intseis_sectiondata_envelope=NULL, **intseis_section_envelope=NULL;
if(LNORM==8){
	intseis_section_envelope = matrix(1,ntr,1,ns);
	intseis_sectiondata_envelope = matrix(1,ntr,1,ns);	
	}
intseis_section = matrix(1,ntr,1,ns);  /* declaration of variables for integration */
intseis_sectiondata = matrix(1,ntr,1,ns);
if(TIMEWIN) picked_times = vector(1,ntr); /* declaration of variables for TIMEWIN */

/* TRACE KILLING */
int ** kill_tmp, *kill_vector;	/* declaration of variables for trace killing */
char trace_kill_file[STRING_SIZE];	
FILE *ftracekill;

if(TRKILL){
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
	   h++;
	}
} /* end if(TRKILL)*/

/* Integration of measured and synthetic data  */
for(i=1;i<=ntr;i++){

	intseis_s=0;
	intseis_sd=0;
	if (VELOCITY==0){	/* only integtration if displacement seismograms are inverted */
      		for(j=1;j<=ns;j++){
			intseis_s += section[i][j];
			intseis_section[i][j]=intseis_s*DT;
			intseis_sd += sectiondata[i][j];
			intseis_sectiondata[i][j]=intseis_sd*DT;
		}
	}
	else{
		for(j=1;j<=ns;j++){
			intseis_section[i][j]=section[i][j];
			intseis_sectiondata[i][j]=sectiondata[i][j];
		}
	}
	
}

/* TIME WINDOWING */
if(TIMEWIN==1){
time_window(intseis_section, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
time_window(intseis_sectiondata, iter, ntr_glob,recpos_loc, ntr, ns, ishot);
}

/* NORMALIZE TRACES */
if(NORMALIZE==1){
normalize_data(intseis_section,ntr,ns);
normalize_data(intseis_sectiondata,ntr,ns);
}

/* calculate envelope for LNORM==8 */
if (LNORM==8){
calc_envelope(intseis_sectiondata,intseis_sectiondata_envelope,ns,ntr);
calc_envelope(intseis_section,intseis_section_envelope,ns,ntr);	
}
/* end of calculate envelope for LNORM==8 */

/* calculate kind of "energy" */

for(i=1;i<=ntr;i++){	
	
    	if((TRKILL==1)&&(kill_vector[i]==1))
    	continue;	
    
	if ((LNORM==5) || (LNORM==7)){
    		abs_sectiondata=0.0;
		abs_section=0.0;
	
		for(j=1;j<=ns;j++){
			abs_sectiondata+=intseis_sectiondata[i][j]*intseis_sectiondata[i][j];
			abs_section+=intseis_section[i][j]*intseis_section[i][j];
		}
		abs_sectiondata=sqrt(abs_sectiondata);
		abs_section=sqrt(abs_section);
	}
	
	/* calculate L2 residuals */
      	
	for(j=1;j<=ns;j++){

		if(LNORM==2){
		L2 += (intseis_section[i][j]-intseis_sectiondata[i][j]) * (intseis_section[i][j]-intseis_sectiondata[i][j]);}

		if((LNORM==5) || (LNORM==7)){
		L2 -= (intseis_sectiondata[i][j]*intseis_section[i][j]) / (abs_sectiondata*abs_section);}
		
		if(LNORM==8){
		L2 += (intseis_section_envelope[i][j]-intseis_sectiondata_envelope[i][j]) * (intseis_section_envelope[i][j]-intseis_sectiondata_envelope[i][j]);}
				                  
	}
}

l2=L2;

/*  printf("\n MYID = %i   IN CALC_MISFIT: L2 = %10.12f \n",MYID,l2); */

/* free memory for integrated seismograms */
free_matrix(intseis_section,1,ntr,1,ns);
free_matrix(intseis_sectiondata,1,ntr,1,ns);

/* free memory for time windowing and trace killing */
if(TIMEWIN==1) free_vector(picked_times,1,ntr);
if(TRKILL==1){
free_imatrix(kill_tmp,1,ntr_glob,1,nsrc_glob);
free_ivector(kill_vector,1,ntr);}
if(LNORM==8){
	free_matrix(intseis_sectiondata_envelope,1,ntr,1,ns);
	free_matrix(intseis_section_envelope,1,ntr,1,ns);
	}
return l2;
} /* end of function */
