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
 *   Read elastic model properties (vp,vs,density) from files
 *
 *  Copyright (c)  T. Bohlen
 *  last update 29.06.2003
 *  ----------------------------------------------------------------------*/


/* This file contains function readmod, which has the purpose
 to read data from model-files for viscoelastic simulation */

#include "fd.h"

void readmod_elastic(float  **  rho, float **  pi, float **  u){
    
    extern int NX, NY, NXG, NYG,  POS[3], MYID, INVMAT1,WAVETYPE;
    extern char  MFILE[STRING_SIZE];
    extern FILE *FP;
    
    
    /* local variables */
    float rhov, muv, piv, vp, vs;
    int i, j, ii, jj;
    FILE *fp_vs, *fp_vp, *fp_rho;
    char filename[STRING_SIZE];
    
    
    
    
	   fprintf(FP,"\n...reading model information from modell-files...\n");
    
	   /* read density and seismic velocities */
	   /* ----------------------------------- */
	   if(INVMAT1==1){
           
           if(WAVETYPE==1||WAVETYPE==3){
               fprintf(FP,"\t Vp:\n\t %s.vp\n\n",MFILE);
               sprintf(filename,"%s.vp",MFILE);
               fp_vp=fopen(filename,"r");
               if (fp_vp==NULL) err(" Could not open model file for Vp ! ");
           }
           
           
           fprintf(FP,"\t Vs:\n\t %s.vs\n\n",MFILE);
           sprintf(filename,"%s.vs",MFILE);
           fp_vs=fopen(filename,"r");
           if (fp_vs==NULL) err(" Could not open model file for Vs ! ");
           
           fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
           sprintf(filename,"%s.rho",MFILE);
           fp_rho=fopen(filename,"r");
           if (fp_rho==NULL) err(" Could not open model file for densities ! ");
           
           
       }
	   
	   /* read density and Lame parameters */
	   /* ----------------------------------- */
	   if(INVMAT1==3){
           fprintf(FP,"\t Lame parameter lambda:\n\t %s.lam\n\n",MFILE);
           sprintf(filename,"%s.lam",MFILE);
           fp_vp=fopen(filename,"r");
           if (fp_vp==NULL) err(" Could not open model file for Lame parameter lambda ! ");
           
           
           fprintf(FP,"\t Lame parameter mu:\n\t %s.mu\n\n",MFILE);
           sprintf(filename,"%s.mu",MFILE);
           fp_vs=fopen(filename,"r");
           if (fp_vs==NULL) err(" Could not open model file for Lame parameter mu ! ");
           
           fprintf(FP,"\t Density:\n\t %s.rho\n\n",MFILE);
           sprintf(filename,"%s.rho",MFILE);
           fp_rho=fopen(filename,"r");
           if (fp_rho==NULL) err(" Could not open model file for densities ! ");
       }
	   
    
    /* loop over global grid */
    for (i=1;i<=NXG;i++){
        for (j=1;j<=NYG;j++){
            
            if(WAVETYPE==1||WAVETYPE==3){
                fread(&vp, sizeof(float), 1, fp_vp);
                
            }
            fread(&vs, sizeof(float), 1, fp_vs);
            fread(&rhov, sizeof(float), 1, fp_rho);
            
            
            /* only the PE which belongs to the current global gridpoint
             is saving model parameters in his local arrays */
            if ((POS[1]==((i-1)/NX)) &&
                (POS[2]==((j-1)/NY))){
                ii=i-POS[1]*NX;
                jj=j-POS[2]*NY;
                
                u[jj][ii]=vs;
                rho[jj][ii]=rhov;
                if(WAVETYPE==1||WAVETYPE==3){
                    pi[jj][ii]=vp;
                }
                
                
                
            }
        }
    }
    if(WAVETYPE==1||WAVETYPE==3){
        fclose(fp_vp);
    }
    fclose(fp_vs);
    fclose(fp_rho);
    
    
    
    /* each PE writes his model to disk */
    if(WAVETYPE==1||WAVETYPE==3){
        if(INVMAT1==1) sprintf(filename,"%s.out.vp",MFILE);
        if(INVMAT1==3) sprintf(filename,"%s.out.pi",MFILE);
        writemod(filename,pi,3);
        MPI_Barrier(MPI_COMM_WORLD);
        if (MYID==0) mergemod(filename,3);
        MPI_Barrier(MPI_COMM_WORLD);
        if(INVMAT1==1) sprintf(filename,"%s.out.vp.%i.%i",MFILE,POS[1],POS[2]);
        if(INVMAT1==3) sprintf(filename,"%s.out.pi.%i.%i",MFILE,POS[1],POS[2]);
        remove(filename);
    }
    
    if(INVMAT1==1) sprintf(filename,"%s.out.vs",MFILE);
    if(INVMAT1==3) sprintf(filename,"%s.out.mu",MFILE);
    writemod(filename,u,3);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (MYID==0) mergemod(filename,3);
    MPI_Barrier(MPI_COMM_WORLD);
    if(INVMAT1==1) sprintf(filename,"%s.out.vs.%i.%i",MFILE,POS[1],POS[2]);
    if(INVMAT1==3) sprintf(filename,"%s.out.mu.%i.%i",MFILE,POS[1],POS[2]);
    remove(filename);
    
    
    sprintf(filename,"%s.out.rho",MFILE);
    writemod(filename,rho,3);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (MYID==0) mergemod(filename,3);
    MPI_Barrier(MPI_COMM_WORLD);
    sprintf(filename,"%s.out.rho.%i.%i",MFILE,POS[1],POS[2]);
    remove(filename);
    
    
}




