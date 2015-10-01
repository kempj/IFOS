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

/*
 * taper model with a gaussian frame to damp inversion artefacts near the sources and receivers
 * M. Schaefer - more or less copied from spat_filt.c - May 2011
 */

#include "fd.h"

void smooth_model(float ** ppi, float ** pu, float ** prho, int iter)
{

	/* extern variables */

        extern float DH;
	extern int FREE_SURF, NX, NY, NXG, NYG;
	extern int NPROCX, NPROCY, MYID, POS[3];
	extern char INV_MODELFILE[STRING_SIZE];
	extern FILE *FP;
	extern int MODEL_FILTER, FILT_DEPTH, FILT_SIZE;
	
	/* local variables */
	int i, j, h, ifw, ii, jj, n, xb, yb, xe, ye, taperlength,taperlength2;
	int ijc, iy, ix, iii, jjj, xx, yy, srctaper_gridpt, i1, j1, filtsize;

	float **modeltmp_vp, **modeltmps_vp, **modeltmp_vs, **modeltmps_vs, grad, normgauss;
	char jac_vp_old[STRING_SIZE];
	char jac_vs_old[STRING_SIZE];
	char jac_rho_old[STRING_SIZE];
	
	char jac_vp[STRING_SIZE];
	char jac_vs[STRING_SIZE];
	char jac_rho[STRING_SIZE];
	
	FILE *model_vp;
	FILE *model_vs;
	FILE *model_rho;
	
	char modfile_vp[STRING_SIZE];
	char modfile_vs[STRING_SIZE];
	char modfile_rho[STRING_SIZE];
		
	/*int FILT_SIZE,FILT_DEPTH;
	FILT_SIZE=10;
	FILT_DEPTH=10;*/
	if(MODEL_FILTER>=1){
	
	modeltmp_vp = matrix(1,NYG,1,NXG);
	modeltmps_vp = matrix(1,NYG,1,NXG);

        modeltmp_vs = matrix(1,NYG,1,NXG);
	modeltmps_vs = matrix(1,NYG,1,NXG);
	
	sprintf(jac_vp_old,"%s_vp_it_%d.bin",INV_MODELFILE,iter);   
	sprintf(jac_vs_old,"%s_vs_it_%d.bin",INV_MODELFILE,iter);   
			
	sprintf(jac_vp,"%s_vp_smoothed_it_%d.bin",INV_MODELFILE,iter);
	sprintf(jac_vs,"%s_vs_smoothed_it_%d.bin",INV_MODELFILE,iter);
		
	/*fprintf(FP,"\n Read modelfile okay! \n");*/
	
	/*printf("\n Read modelfile okay! \n");*/
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/*printf("\n nach MPI_Barrier - Read modelfile okay! \n");*/
	
	if(MYID==0){
	model_vp=fopen(jac_vp_old,"rb");
	if (model_vp==NULL) err(" Could not open vp model file ! ");
	/* load merged model */
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
                        fread(&grad, sizeof(float), 1, model_vp);
			modeltmp_vp[j][i]=grad;
			modeltmps_vp[j][i]=grad;
		}
	}
	
	fclose(model_vp);
	
	/* apply spatial Gaussian filter on model */
	for (i=1+FILT_SIZE;i<=NXG-FILT_SIZE;i++){
	   for (j=1+FILT_SIZE;j<=NYG-FILT_SIZE;j++){

                     modeltmps_vp[j][i]=0.0;   
		     normgauss=0.0;
		     
		     for (i1=i-FILT_SIZE;i1<=i+FILT_SIZE;i1++){
	                 for (j1=j-FILT_SIZE;j1<=j+FILT_DEPTH;j1++){
			 
                         modeltmps_vp[j][i] += modeltmp_vp[j1][i1] * (1/(sqrt(2*PI)*FILT_SIZE)) * exp((-1/2)*((sqrt(((i-i1)*(i-i1))+((i-j1)*(i-j1))))*(sqrt(((i-i1)*(i-i1))+((i-j1)*(i-j1)))))/(FILT_SIZE*FILT_SIZE)); 
			 normgauss += (1/(sqrt(2*PI)*FILT_SIZE)) * exp((-1/2)*((sqrt(((i-i1)*(i-i1))+((i-j1)*(i-j1))))*(sqrt(((i-i1)*(i-i1))+((i-j1)*(i-j1)))))/(FILT_SIZE*FILT_SIZE)); 
			 
			 }
		     }
		     
		     modeltmps_vp[j][i]=modeltmps_vp[j][i]/normgauss;		  
			
	   }
	}
	
        model_vp=fopen(jac_vp,"wb");
        
	/* output of the preconditioned model */
        for (i=1;i<=NXG;i++){
           for (j=1;j<=NYG;j++){
            
           fwrite(&modeltmps_vp[j][i],sizeof(float),1,model_vp);

           }
        }

        fclose(model_vp);
	
	model_vs=fopen(jac_vs_old,"rb");
	if (model_vs==NULL) err(" Could not open vs model file ! ");
	
	/* load merged model */
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
                        fread(&grad, sizeof(float), 1, model_vs);
			modeltmp_vs[j][i]=grad;
			modeltmps_vs[j][i]=grad;
		}
	}
	
	fclose(model_vs);
	
	/* apply spatial Gaussian filter on model */
	for (i=1+FILT_SIZE;i<=NXG-FILT_SIZE;i++){
	   for (j=1+FILT_SIZE;j<=NYG-FILT_SIZE;j++){

                     modeltmps_vs[j][i]=0.0;   
		     normgauss=0.0;
		     
		     for (i1=i-FILT_SIZE;i1<=i+FILT_SIZE;i1++){
	                  
                         modeltmps_vs[j][i] += modeltmp_vs[j][i1] * (1/(sqrt(2*PI)*FILT_SIZE)) * exp((-1/2)*((sqrt(((i-i1)*(i-i1))+((i-i1)*(i-i1))))*(sqrt(((i-i1)*(i-i1))+((i-i1)*(i-i1)))))/(FILT_SIZE*FILT_SIZE)); 
			 normgauss += (1/(sqrt(2*PI)*FILT_SIZE)) * exp((-1/2)*((sqrt(((i-i1)*(i-i1))+((i-i1)*(i-i1))))*(sqrt(((i-i1)*(i-i1))+((i-i1)*(i-i1)))))/(FILT_SIZE*FILT_SIZE));  
		     }
		     
		     modeltmps_vs[j][i]=modeltmps_vs[j][i]/normgauss;		  
			
	   }
	}
	
        model_vs=fopen(jac_vs,"wb");

        /* output of the preconditioned model */
        for (i=1;i<=NXG;i++){
           for (j=1;j<=NYG;j++){
            
           fwrite(&modeltmps_vs[j][i],sizeof(float),1,model_vs);

           }
        }

        fclose(model_vs);
	
	}
	
	fprintf(FP,"\n \t ---- Models are smoothed up to a depth of %d gridpoints ---- \n",FILT_DEPTH); 
	fprintf(FP," \t ---- and for the smoothing the half filter length is %d gridpoints. ---- \n",FILT_SIZE);
	/*printf(" MYID = %d , Models are smoothed ... \n",MYID);*/
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	model_vp=fopen(jac_vp,"rb");
	if (model_vp==NULL) err(" Could not open vp model file ! ");
	/* distribute spatial filtered gradient on computational nodes */
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
	   
                        fread(&grad, sizeof(float), 1, model_vp);
			   			
			if ((POS[1]==((i-1)/NX)) && 
		   	 (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				ppi[jj][ii]=grad;

			}
		}
	}
		
	fclose(model_vp);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	model_vs=fopen(jac_vs,"rb");
	if (model_vs==NULL) err(" Could not open vs model file ! ");
	/* distribute spatial filtered gradient on computational nodes */
	for (i=1;i<=NXG;i++){
	   for (j=1;j<=NYG;j++){
	   
                        fread(&grad, sizeof(float), 1, model_vs);
			   			
			if ((POS[1]==((i-1)/NX)) && 
		   	 (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				pu[jj][ii]=grad;

			}
		}
	}
		
	fclose(model_vs);
	
	MPI_Barrier(MPI_COMM_WORLD); 

	fprintf(FP,"\n ---- Smoothed models are distributed on computational nodes ... ---- \n");
	/*printf(" MYID = %d , Smoothed models are distributed on computational nodes ... \n",MYID);*/

        sprintf(modfile_vp,"%s_vp_smoothed_it_%d.bin",INV_MODELFILE,iter);
	writemod(modfile_vp,ppi,3);
                                        
	MPI_Barrier(MPI_COMM_WORLD);
                                                
	if (MYID==0) mergemod(modfile_vp,3);
                                                        
	sprintf(modfile_vs,"%s_vs_smoothed_it_%d.bin",INV_MODELFILE,iter);
	writemod(modfile_vs,pu,3);
	
	MPI_Barrier(MPI_COMM_WORLD);
                                                                                                        
	if (MYID==0) mergemod(modfile_vs,3);                                                                                                                

	MPI_Barrier(MPI_COMM_WORLD);
	
	free_matrix(modeltmp_vs,1,NXG,1,NYG);
	free_matrix(modeltmps_vs,1,NXG,1,NYG);
	
	free_matrix(modeltmp_vp,1,NXG,1,NYG);
	free_matrix(modeltmps_vp,1,NXG,1,NYG);
	
	}/* end of application condition for the spatial filter */
} 
