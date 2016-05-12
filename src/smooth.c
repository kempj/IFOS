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

/*
 * Smoothing gradient / model with a median filter
 */

#include "fd.h"

void smooth(float ** mat, int sws, int filter, float Vs_avg, float F_LOW_PASS)
{
    
    /* extern variables */
    
    extern float DH, A;
    extern int FREE_SURF, NX, NY, NXG, NYG;
    extern int NPROCX, NPROCY, MYID, POS[3],WAVETYPE,VERBOSE;
    extern char JACOBIAN[STRING_SIZE], INV_MODELFILE[STRING_SIZE];
    extern FILE *FP;
    extern int FILT_SIZE_GRAD, FILT_SIZE, MODEL_FILTER, GRAD_FILTER, TIME_FILT,GRAD_FILT_WAVELENGTH;
    extern int VERBOSE;
    /* local variables */
    int i, j, ii, jj;
    int i1, j1, filtsize, hfs;
    
    float **model_tmp, **model_med, **filterpart, grad, normgauss, smooth_meter;
    
    char jac_tmp[STRING_SIZE];
    
    FILE *model;
    
    char modfile[STRING_SIZE];
    
    if(sws==1){
        sprintf(jac_tmp,"%s_g",JACOBIAN);
        write_matrix_disk(mat, jac_tmp);
    }
    
    if(sws==2){
        sprintf(jac_tmp,"%s_g_u",JACOBIAN);
        write_matrix_disk(mat, jac_tmp);
    }
    
    if(sws==3){
        sprintf(jac_tmp,"%s_g_rho",JACOBIAN);
        write_matrix_disk(mat, jac_tmp);
    }
    
    
    if(MYID==0){
        switch (filter){
            case 1:
                if((GRAD_FILT_WAVELENGTH==1)&&(TIME_FILT==1)){
                    if(VERBOSE) printf("\n -------------------------------------------------------------------------- \n");
                    if(VERBOSE) printf("\n Calculating a wavelength dependent filter size for smoothing the gradient: \n");
                    FILT_SIZE_GRAD = (int)(Vs_avg/F_LOW_PASS*A/DH);
                    if(VERBOSE) printf("\n FILT_SIZE_GRAD = Vs_avg = %4.2f m/s / F_LOW_PASS = %4.2f Hz * weighting factor A = %4.2f / grid spacing DH = %4.2f m  \n",Vs_avg,F_LOW_PASS,A,DH);
                    if(VERBOSE) printf("\n New FILT_SIZE_GRAD = %d (grid points) is used (-> %4.2f m).                \n",FILT_SIZE_GRAD,FILT_SIZE_GRAD*DH);
                }
                if (FILT_SIZE_GRAD==0)	return;
                if (!(FILT_SIZE_GRAD % 2)) {
                    if (FILT_SIZE_GRAD > 0)	FILT_SIZE_GRAD += 1;
                    else			FILT_SIZE_GRAD -= 1;
                }
                hfs = abs(FILT_SIZE_GRAD)/2;
                if(VERBOSE) printf("\n ----------------------------------------------------------------\n");
                if(VERBOSE) printf("\n Filter size is %d gridpoints, half filter size is %d gridpoints.\n",FILT_SIZE_GRAD,hfs);
                filterpart=matrix(1,abs(FILT_SIZE_GRAD),1,abs(FILT_SIZE_GRAD));
                model_tmp = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
                break;
                
            case 2:
                if (FILT_SIZE==0)	return;
                if (!(FILT_SIZE % 2)) {
                    if (FILT_SIZE > 0)	FILT_SIZE += 1;
                    else			FILT_SIZE -= 1;
                }
                hfs = abs(FILT_SIZE)/2;
                if(VERBOSE) printf("\n ----------------------------------------------------------------\n");
                if(VERBOSE) printf("\n Filter size is %d gridpoints, half filter size is %d gridpoints.\n",FILT_SIZE,hfs);
                filterpart=matrix(1,abs(FILT_SIZE),1,abs(FILT_SIZE));
                model_tmp = matrix(-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
                break;
        }
        model_med = matrix(1,NYG,1,NXG);
        
        if(sws==1){
            sprintf(jac_tmp,"%s_g.bin",JACOBIAN);
        }
        
        if(sws==2){
            sprintf(jac_tmp,"%s_g_u.bin",JACOBIAN);
        }
        
        if(sws==3){
            sprintf(jac_tmp,"%s_g_rho.bin",JACOBIAN);
        }
        
        if(sws==4){
            sprintf(jac_tmp,"%s_vp.bin",INV_MODELFILE);}
        
        if(sws==5){
            sprintf(jac_tmp,"%s_vs.bin",INV_MODELFILE);}
        
        if(sws==6){
            sprintf(jac_tmp,"%s_rho.bin",INV_MODELFILE);}
        
        model=fopen(jac_tmp,"rb");
        if (model==NULL) declare_error(" Could not open file !");
        
        /* load merged model */
        for (i=1;i<=NXG;i++){
            for (j=1;j<=NYG;j++){
                fread(&grad, sizeof(float), 1, model);
                model_tmp[j][i]=grad;
            }
        }
        
		      fclose(model);
        
		/* apply 2D-Gaussian filter on vp and vs model */
        /* extrapolate array */
        /* left/right boundary */
        for (j=1;j<=NYG;j++){
            for (i=-hfs+1;i<=0;i++){
                model_tmp[j][i] = model_tmp[j][1];}
            for (i=NXG+1;i<=NXG+hfs;i++){
                model_tmp[j][i] = model_tmp[j][NXG];}
        }
        /* top/bottom boundary incl. corners */
        for (j=-hfs+1;j<=0;j++){
            for (i=-hfs+1;i<=NXG+hfs;i++){
                model_tmp[j][i] = model_tmp[1][i];}
        }
        for (j=NYG+1;j<=NYG+hfs;j++){
            for (i=-hfs+1;i<=NXG+hfs;i++){
                model_tmp[j][i] = model_tmp[NYG][i];}
        }
        
        /* filter */
        for (j=1;j<=NYG;j++){
            for (i=1;i<=NXG;i++){
                /* create a filtersize x filtersize matrix */
                for (ii=-hfs;ii<=hfs;ii++){
                    for (jj=-hfs;jj<=hfs;jj++){
                        
                        filterpart[ii+hfs+1][jj+hfs+1] = model_tmp[j+jj][i+ii];
                    }
                }
                /* filter */
                switch (filter){
                    case 1:
                        model_med[j][i] = median2d(filterpart,abs(FILT_SIZE_GRAD),abs(FILT_SIZE_GRAD));
                        break;
                    case 2:
                        model_med[j][i] = median2d(filterpart,abs(FILT_SIZE),abs(FILT_SIZE));
                        break;
                }
                
            }
        }
        
        /* output of smoothed gradient / model*/
        if(sws==1){
            sprintf(jac_tmp,"%s_tmp_g.bin",JACOBIAN);}
        
        if(sws==2){
            sprintf(jac_tmp,"%s_tmp_g_u.bin",JACOBIAN);}
        
        if(sws==3){
            sprintf(jac_tmp,"%s_tmp_g_rho.bin",JACOBIAN);}
        
        if(sws==4){
            sprintf(jac_tmp,"%s_vp_tmp.bin",INV_MODELFILE);}
        
        if(sws==5){
            sprintf(jac_tmp,"%s_vs_tmp.bin",INV_MODELFILE);}
        
        if(sws==6){
            sprintf(jac_tmp,"%s_rho_tmp.bin",INV_MODELFILE);}
        
        model=fopen(jac_tmp,"wb");
        for (i=1;i<=NXG;i++){
            for (j=1;j<=NYG;j++){
                
                fwrite(&model_med[j][i],sizeof(float),1,model);
                
            }
        }
        fclose(model);
        
        switch (filter){
            case 1:
                free_matrix(model_tmp,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
                free_matrix(filterpart,1,abs(FILT_SIZE_GRAD),1,abs(FILT_SIZE_GRAD));
                break;
            case 2:
                free_matrix(model_tmp,-hfs+1,NYG+hfs,-hfs+1,NXG+hfs);
                free_matrix(filterpart,1,abs(FILT_SIZE),1,abs(FILT_SIZE));
                break;
        }
        free_matrix(model_med,1,NXG,1,NYG);
        
    }/* end of if(MYID==0)*/
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    switch (filter){
        case 1:
            smooth_meter=FILT_SIZE_GRAD*DH;
            if(VERBOSE) fprintf(FP,"\n Gradient %s is smoothed with filter length of %4.2f meter.\n",jac_tmp,smooth_meter);
            break;
        case 2:
            smooth_meter=FILT_SIZE*DH;
            if(VERBOSE) fprintf(FP,"\n Model %s is smoothed with filter length of %4.2f meter.\n",jac_tmp,smooth_meter);
            break ;
    }
    /* distribute smoothed gradient / model on computational nodes */
    if(sws==1){
        sprintf(jac_tmp,"%s_tmp_g.bin",JACOBIAN);}
    
    if(sws==2){
        sprintf(jac_tmp,"%s_tmp_g_u.bin",JACOBIAN);}
    
    if(sws==3){
        sprintf(jac_tmp,"%s_tmp_g_rho.bin",JACOBIAN);}
    
    if(sws==4){
        sprintf(jac_tmp,"%s_vp_tmp.bin",INV_MODELFILE);}   
    
    if(sws==5){
        sprintf(jac_tmp,"%s_vs_tmp.bin",INV_MODELFILE);}    
    
    if(sws==6){
        sprintf(jac_tmp,"%s_rho_tmp.bin",INV_MODELFILE);}
				
    model=fopen(jac_tmp,"rb");
    if (model==NULL) declare_error(" Could not open file! (distribute)");
    for (i=1;i<=NXG;i++){
        for (j=1;j<=NYG;j++){
            
            fread(&grad, sizeof(float), 1, model);
            
            if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
                ii=i-POS[1]*NX;
                jj=j-POS[2]*NY;
                
                mat[jj][ii]=grad;
                
            }
        }
    }
    
    fclose(model);
    
    if(VERBOSE) fprintf(FP,"\n Smoothed gradient / model is distributed on computational nodes.\n");
    if(VERBOSE) fprintf(FP,"\n ----------------------------------------------------------------\n");
}/* end of smoothing */
