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
 *   Apply time damping (after Brossier (2009))                                 
 *   last update 31/08/11, D.Koehn
 *   modified    02/02/12, S.Heider
 *  ----------------------------------------------------------------------*/
#include "fd.h"

void time_window(float **sectiondata, int iter, int ntr_glob, int **recpos_loc, int ntr, int ns, int ishot){

    /* declaration of variables */
    extern float DT;
    extern float GAMMA, TWLENGTH_PLUS, TWLENGTH_MINUS;
    extern int TW_IND, USE_WORKFLOW, WORKFLOW_STAGE;
    extern char PICKS_FILE[STRING_SIZE];
    char pickfile_char[STRING_SIZE];
    float time, dump, dump1, dump2, dump3, taper, taper1;
    float *pick_tmp, **pick_tmp_m;
    int i, j, h;

    float *picked_times=NULL, **picked_times_m=NULL;
    
    FILE *fptime;

    if(TW_IND){
        picked_times_m = matrix(1,3,1,ntr_glob);
        pick_tmp_m = matrix(1,3,1,ntr_glob);
    }else{
        picked_times = vector(1,ntr_glob);
        pick_tmp = vector(1,ntr_glob);
    }
    
    /* read picked first arrival times */
    if(USE_WORKFLOW){
        sprintf(pickfile_char,"%s_%i_%i.dat",PICKS_FILE,ishot,WORKFLOW_STAGE);
        fptime=fopen(pickfile_char,"r");
        if (fptime == NULL) {
            sprintf(pickfile_char,"%s_%i.dat",PICKS_FILE,ishot);
            fptime=fopen(pickfile_char,"r");
            if (fptime == NULL) {
                err(" picks_?.dat could not be opened !");
            }
        }
    }else{
        sprintf(pickfile_char,"%s_%i.dat",PICKS_FILE,ishot);
        fptime=fopen(pickfile_char,"r");
        if (fptime == NULL) {
            err(" picks_?.dat could not be opened !");
        }
    }
    
    if(TW_IND){
        for(i=1;i<=ntr_glob;i++){
            fscanf(fptime,"%f%f%f",&dump,&dump2,&dump3);
            pick_tmp_m[1][i] = dump;
            pick_tmp_m[2][i] = dump2;
            pick_tmp_m[3][i] = dump3;
        }
    }else{
        for(i=1;i<=ntr_glob;i++){
            fscanf(fptime,"%f",&dump);
            pick_tmp[i] = dump;
        }
    }
    
    fclose(fptime);
    
    /* distribute picks on CPUs */
    h=1;
    if(TW_IND){
        for(i=1;i<=ntr;i++){
            picked_times_m[1][h] = pick_tmp_m[1][recpos_loc[3][i]];
            picked_times_m[2][h] = pick_tmp_m[2][recpos_loc[3][i]];
            picked_times_m[3][h] = pick_tmp_m[3][recpos_loc[3][i]];
            
            h++;
        }
    }else{
        for(i=1;i<=ntr;i++){
            picked_times[h] = pick_tmp[recpos_loc[3][i]];
            
            h++;
        }
    }
    
    if(TW_IND)
        free_matrix(pick_tmp_m,1,3,1,ntr_glob);
    else
        free_vector(pick_tmp,1,ntr_glob);
    
    if(TW_IND){
        for(i=1;i<=ntr;i++){
        for(j=2;j<=ns;j++){
        
            time = (float)(j * DT);
            
            dump = (time-picked_times_m[1][i]-picked_times_m[2][i]);
            taper = exp(-GAMMA*dump*dump);
            
            dump1 = (time-picked_times_m[1][i]+picked_times_m[3][i]); 
            taper1 = exp(-GAMMA*dump1*dump1);
            
            if(time>=picked_times_m[1][i]+picked_times_m[2][i]){
            sectiondata[i][j] = sectiondata[i][j] * taper;}
            
            if(time<=picked_times_m[1][i]-picked_times_m[3][i]){
            sectiondata[i][j] = sectiondata[i][j] * taper1;}
            
            sectiondata[i][j] = sectiondata[i][j];
            
        }
        }
    }else{
        for(i=1;i<=ntr;i++){
        for(j=2;j<=ns;j++){
        
            time = (float)(j * DT);
            
            dump = (time-picked_times[i]-TWLENGTH_PLUS);
            taper = exp(-GAMMA*dump*dump);
            
            dump1 = (time-picked_times[i]+TWLENGTH_MINUS); 
            taper1 = exp(-GAMMA*dump1*dump1);
            
            if(time>=picked_times[i]+TWLENGTH_PLUS){
            sectiondata[i][j] = sectiondata[i][j] * taper;}
            
            if(time<=picked_times[i]-TWLENGTH_MINUS){
            sectiondata[i][j] = sectiondata[i][j] * taper1;}
            
            sectiondata[i][j] = sectiondata[i][j];
            
        }     
        }
    }
    
    if(TW_IND)
        free_matrix(picked_times_m,1,3,1,ntr_glob);
    else
        free_vector(picked_times,1,ntr_glob);
    
} /* end of function time_window.c */