/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2015  For the list of authors, see file AUTHORS.
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
 *  Restart workflow                                Oct 2015
 *
 *  written by Wittkamp, parts copied out of readmod.c
 *
 *  If you want to adjust the workflow, I think it is not necessary to
 *  modify this file. Have a look at apply_workflow.c and adjust
 *  WORKFLOW_MAX_VAR in fd.h.
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void restart_workflow(int iter,float  **  rho, float **  pi, float **  u,float ** taus, float **taup){
    
    extern int NX, NY, NXG, NYG,  POS[3], MYID, INVMAT1,L,WAVETYPE,ACOUSTIC;
    extern char  MFILE[STRING_SIZE],INV_MODELFILE[STRING_SIZE];
    extern FILE *FP;
    
    /* local variables */
    float rhov, muv, piv, vp, vs;
    float taupp, tauss;
    int i, j, ii, jj;
    FILE *fp_vs, *fp_vp, *fp_rho;
    FILE *fp_taus, *fp_taup;
    char modfile[STRING_SIZE];
    
    
    fprintf(FP,"\n Restart inversion from workflow \n");
    
    /* Check dependencies */
    if(INVMAT1!=1) err(" Restart of inversion only if INVMAT1=1! ");
    
    /* Read in old model files */
    if(INVMAT1==1){
        if(WAVETYPE==1||WAVETYPE==3){
            if(!ACOUSTIC){
                sprintf(modfile,"%s_vp_it%d.bin",INV_MODELFILE,iter);
                fprintf(FP,"\n trying to open %s",modfile);
                fp_vp=fopen(modfile,"r");
                if (fp_vp==NULL) {
                    fprintf(FP,"\n ****** FAILSAFE-MODE VP ******");
                    fprintf(FP,"\n NOT able to open %s",modfile);
                    sprintf(modfile,"%s_vp.bin",INV_MODELFILE);
                    fprintf(FP,"\n trying to open %s\n",modfile);
                    fp_vp=fopen(modfile,"r");
                    if (fp_vp==NULL) {
                        err(" Could not open model file for Vp to restart inversion! ");
                    }
                }
            }
        }
        sprintf(modfile,"%s_vs_it%d.bin",INV_MODELFILE,iter);
        fprintf(FP,"\n trying to open %s",modfile);
        fp_vs=fopen(modfile,"r");
        if (fp_vs==NULL) {
            fprintf(FP,"\n ****** FAILSAFE-MODE VS ******");
            fprintf(FP,"\n NOT able to open %s",modfile);
            sprintf(modfile,"%s_vs.bin",INV_MODELFILE);
            fprintf(FP,"\n trying to open %s\n",modfile);
            fp_vs=fopen(modfile,"r");
            if (fp_vs==NULL) {
                err(" Could not open model file for Vs to restart inversion! ");
            }
        }
        
        sprintf(modfile,"%s_rho_it%d.bin",INV_MODELFILE,iter);
        fprintf(FP,"\n trying to open %s",modfile);
        fp_rho=fopen(modfile,"r");
        if (fp_rho==NULL) {
            fprintf(FP,"\n ****** FAILSAFE-MODE RHO ******");
            fprintf(FP,"\n NOT able to open %s",modfile);
            sprintf(modfile,"%s_rho.bin",INV_MODELFILE);
            fprintf(FP,"\n trying to open %s\n",modfile);
            fp_rho=fopen(modfile,"r");
            if (fp_rho==NULL) {
                err(" Could not open model file for Rho to restart inversion! ");
            }
        }
        
        if(L) {
            
            /*******/
            iter=0; // Up to now no tau inversion, therefore we open at *it0.bin
            /*******/
            
            sprintf(modfile,"%s_taus_it%d.bin",INV_MODELFILE,iter);
            fprintf(FP,"\n trying to open %s",modfile);
            fp_taus=fopen(modfile,"r");
            if (fp_taus==NULL) err(" Could not open model file for Taus to restart inversion ! ");
            
            sprintf(modfile,"%s_taup_it%d.bin",INV_MODELFILE,iter);
            fprintf(FP,"\n trying to open %s",modfile);
            fp_taup=fopen(modfile,"r");
            if (fp_taup==NULL) err(" Could not open model file for Taup to restart inversion ! ");
        }
        
    }
    
    /* loop over global grid */
    for (i=1;i<=NXG;i++){
        for (j=1;j<=NYG;j++){
            if(WAVETYPE==1||WAVETYPE==3&&!ACOUSTIC) fread(&vp, sizeof(float), 1, fp_vp);
            fread(&vs, sizeof(float), 1, fp_vs);
            fread(&rhov, sizeof(float), 1, fp_rho);
            if(L) {
                fread(&tauss, sizeof(float), 1, fp_taus);
                fread(&taupp, sizeof(float), 1, fp_taup);
            }
            
            /* only the PE which belongs to the current global gridpoint
             is saving model parameters in his local arrays */
            if ((POS[1]==((i-1)/NX)) &&
                (POS[2]==((j-1)/NY))){
                ii=i-POS[1]*NX;
                jj=j-POS[2]*NY;
                
                u[jj][ii]=vs;
                rho[jj][ii]=rhov;
                if(WAVETYPE==1||WAVETYPE==3&&!ACOUSTIC) pi[jj][ii]=vp;
                
                if(L) {
                    taus[jj][ii]=tauss;
                    taup[jj][ii]=taupp;
                }
                
            }
        }
    }
    
    /* each PE writes his model to disk */
    if(WAVETYPE==1||WAVETYPE==3){
        if(!ACOUSTIC){
            sprintf(modfile,"%s.restart.vp",INV_MODELFILE);
            writemod(modfile,pi,3);
            MPI_Barrier(MPI_COMM_WORLD);
            if (MYID==0) mergemod(modfile,3);
            
            MPI_Barrier(MPI_COMM_WORLD);
            sprintf(modfile,"%s.restart.vp.%i.%i",INV_MODELFILE,POS[1],POS[2]);
            remove(modfile);
        }
    }
    /* each PE writes his model to disk */
    sprintf(modfile,"%s.restart.vs",INV_MODELFILE);
    writemod(modfile,u,3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYID==0) mergemod(modfile,3);
    
    MPI_Barrier(MPI_COMM_WORLD);
    sprintf(modfile,"%s.restart.vs.%i.%i",INV_MODELFILE,POS[1],POS[2]);
    remove(modfile);
    
    /* each PE writes his model to disk */
    sprintf(modfile,"%s.restart.rho",INV_MODELFILE);
    writemod(modfile,rho,3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (MYID==0) mergemod(modfile,3);
    
    MPI_Barrier(MPI_COMM_WORLD);
    sprintf(modfile,"%s.restart.rho.%i.%i",INV_MODELFILE,POS[1],POS[2]);
    remove(modfile);
    
    if(L) {
        /* each PE writes his model to disk */
        sprintf(modfile,"%s.restart.taup",INV_MODELFILE);
        writemod(modfile,taup,3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID==0) mergemod(modfile,3);
        
        MPI_Barrier(MPI_COMM_WORLD);
        sprintf(modfile,"%s.restart.taup.%i.%i",INV_MODELFILE,POS[1],POS[2]);
        remove(modfile);
        
        /* each PE writes his model to disk */
        sprintf(modfile,"%s.restart.taus",INV_MODELFILE);
        writemod(modfile,taus,3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID==0) mergemod(modfile,3);
        
        MPI_Barrier(MPI_COMM_WORLD);
        sprintf(modfile,"%s.restart.taus.%i.%i",INV_MODELFILE,POS[1],POS[2]);
        remove(modfile);
    }
}
