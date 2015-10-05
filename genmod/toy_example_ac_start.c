/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2005  <name of author>
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

/* $Id: hh_el.c,v 2.3 2007/08/21 13:16:19 tbohlen Exp $*/
/*
 *   Model defined by flnode file. 
 */

#include "fd.h"
#include "../src/fd.h"

void model_acoustic(float **rho, float **pi){

	/* extern variables */
	extern int NX, NY, NXG, NYG, POS[3], MYID, FDORDER;
	extern float DH, DT;
	extern char MFILE[STRING_SIZE];
	extern char INV_MODELFILE[STRING_SIZE];
	extern char TAPER_FILE_NAME[STRING_SIZE], TAPER_FILE_NAME_RHO[STRING_SIZE];

	/* local variables */
	float piv, vp, rhov, taperv, **taper;
	int i, j, ii, jj, l, nd;
	char modfile[STRING_SIZE];
	
	FILE *flfile;
	int nodes;
	char cline[256];
	
	float *fldepth, *flrho, *flvp, *fltaper;
	
	nd = FDORDER/2 + 1;
	
	/*read FL nodes from File*/
	nodes=5;
	
	fldepth=vector(1,nodes);
	flrho=vector(1,nodes);
	flvp=vector(1,nodes);
	fltaper=vector(1,nodes);
	
	taper=matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	
	flfile=fopen("model_true/flnodes.toy_example_ac.start","r");
	if (flfile==NULL) err(" FL-file could not be opened !");
	
	
	for (l=1;l<=nodes;l++){
		fgets(cline,255,flfile);
		if (cline[0]!='#'){
			sscanf(cline,"%f%f%f%f",&fldepth[l], &flrho[l], &flvp[l], &fltaper[l]);
		}
		else l=l-1;
	
	}
	
	if(MYID==0){
	printf(" ------------------------------------------------------------------ \n\n");
	printf(" Information of FL nodes: \n\n");
	printf(" \t depth \t vp \n\n");
	
	for (l=1;l<=nodes;l++){
	printf(" \t %f \t %f\n\n",fldepth[l],flvp[l]);
	}
	printf(" ------------------------------------------------------------------ \n\n");
	}
	
	/* loop over global grid */
		for (i=1;i<=NXG;i++){
			for (l=1;l<nodes;l++){
				if (fldepth[l]==fldepth[l+1]){
					if ((i==1) && (MYID==0)){
					printf("depth: %f m: double node\n",fldepth[l]);}}
				else{
					for (j=(int)(fldepth[l]/DH)+1;j<=(int)(fldepth[l+1]/DH);j++){
						
						vp=0.0;
						
						vp=(DH*(j-1)-fldepth[l])*(flvp[l+1]-flvp[l])/(fldepth[l+1]-fldepth[l])+flvp[l];
						vp=vp*1000.0;
						
						piv=vp;
						
						/* only the PE which belongs to the current global gridpoint 
						is saving model parameters in his local arrays */
						if ((POS[1]==((i-1)/NX)) && 
						(POS[2]==((j-1)/NY))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;
						
						pi[jj][ii]=piv;
						}
					}
				}
			}
			
			for (j=(int)(fldepth[nodes]/DH)+1;j<=NYG;j++){
			  
				vp=0.0;
				vp=flvp[nodes]*1000.0;
				
				piv=vp;

				/* only the PE which belongs to the current global gridpoint 
				  		is saving model parameters in his local arrays */
						if ((POS[1]==((i-1)/NX)) && 
				    		(POS[2]==((j-1)/NY))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;

						pi[jj][ii]=piv;
						}
			}
		}
		
	/* loop over global grid - taper */
		for (i=1;i<=NXG;i++){
			for (l=1;l<nodes;l++){
				if (fldepth[l]==fldepth[l+1]){
					if ((i==1) && (MYID==0)){
					printf("depth: %f m: double node\n",fldepth[l]);}}
				else{
					for (j=(int)(fldepth[l]/DH)+1;j<=(int)(fldepth[l+1]/DH);j++){
						
						taperv=0.0;
						
						taperv=(DH*(j-1)-fldepth[l])*(fltaper[l+1]-fltaper[l])/(fldepth[l+1]-fldepth[l])+fltaper[l];
						
						/* only the PE which belongs to the current global gridpoint 
						is saving model parameters in his local arrays */
						if ((POS[1]==((i-1)/NX)) && 
						(POS[2]==((j-1)/NY))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;
						
						taper[jj][ii]=taperv;
						}
					}
				}
			}
			
			for (j=(int)(fldepth[nodes]/DH)+1;j<=NYG;j++){
			  
				taperv=0.0;
				taperv=fltaper[nodes];
				
				/* only the PE which belongs to the current global gridpoint 
						is saving model parameters in his local arrays */
						if ((POS[1]==((i-1)/NX)) && 
						(POS[2]==((j-1)/NY))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;

						taper[jj][ii]=taperv;
						}
			}
		}
		
	free_vector(fldepth,1,nodes);
	free_vector(flrho,1,nodes);
	free_vector(flvp,1,nodes);
	free_vector(fltaper,1,nodes);
	
	
	/**************************************************/
	/* creation of density models */
	/**************************************************/
	
	/*read FL nodes from File*/
	nodes=6;
	fldepth=vector(1,nodes);
	flrho=vector(1,nodes);
	flvp=vector(1,nodes);	
	flfile=fopen("model_true/flnodes.toy_example_ac","r");
	if (flfile==NULL) err(" FL-file could not be opened !");
	
	for (l=1;l<=nodes;l++){
		fgets(cline,255,flfile);
		if (cline[0]!='#'){
			sscanf(cline,"%f%f%f",&fldepth[l], &flrho[l], &flvp[l]);
		}
		else l=l-1;
	}
	
	
	if(MYID==0){
	printf(" ------------------------------------------------------------------ \n\n");
	printf(" Information of FL nodes: \n\n");
	printf(" \t depth \t vp \t rho \n\n");
	
	for (l=1;l<=nodes;l++){
	printf(" \t %f \t %f \t %f \n\n",fldepth[l],flvp[l],flrho[l]);
	}
	printf(" ------------------------------------------------------------------ \n\n");
	}
	
	/* loop over global grid */
	for (i=1;i<=NXG;i++){
		for (l=1;l<nodes;l++){
			if (fldepth[l]==fldepth[l+1]){
				if ((i==1) && (MYID==0)){
					printf("depth: %f m: double node\n",fldepth[l]);
				}
			}else{
				for (j=(int)(fldepth[l]/DH)+1;j<=(int)(fldepth[l+1]/DH);j++){
		
					rhov=0.0;
				
					rhov=(DH*(j-1)-fldepth[l])*(flrho[l+1]-flrho[l])/(fldepth[l+1]-fldepth[l])+flrho[l];
					rhov=rhov*1000.0;
					
					
					/* only the PE which belongs to the current global gridpoint 
					is saving model parameters in his local arrays */
					if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
						ii=i-POS[1]*NX;
						jj=j-POS[2]*NY;

						rho[jj][ii]=rhov;
					}
				}
			}
		}
		
		for (j=(int)(fldepth[nodes]/DH)+1;j<=NYG;j++){
			
			rhov=0.0;
			rhov=flrho[nodes]*1000.0;
			
			
			/* only the PE which belongs to the current global gridpoint 
			is saving model parameters in his local arrays */
			if ((POS[1]==((i-1)/NX)) && (POS[2]==((j-1)/NY))){
				ii=i-POS[1]*NX;
				jj=j-POS[2]*NY;

				rho[jj][ii]=rhov;
			}
		}
	}
	
	free_vector(fldepth,1,nodes);
	free_vector(flrho,1,nodes);
	free_vector(flvp,1,nodes);	
	
	
        sprintf(modfile,"%s_rho_it0.bin",INV_MODELFILE);
	writemod(modfile,rho,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_rho_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
	remove(modfile);
	
	sprintf(modfile,"%s_vp_it0.bin",INV_MODELFILE);
	writemod(modfile,pi,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s_vp_it0.bin.%i%i",INV_MODELFILE,POS[1],POS[2]);
	remove(modfile);
	
	sprintf(modfile,"%s",TAPER_FILE_NAME);
	writemod(modfile,taper,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s.%i%i",TAPER_FILE_NAME,POS[1],POS[2]);
	remove(modfile);
	
	sprintf(modfile,"%s",TAPER_FILE_NAME_RHO);
	writemod(modfile,taper,3);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0) mergemod(modfile,3);
	MPI_Barrier(MPI_COMM_WORLD); 
	sprintf(modfile,"%s.%i%i",TAPER_FILE_NAME_RHO,POS[1],POS[2]);
	remove(modfile);
	
	free_matrix(taper,-nd+1,NY+nd,-nd+1,NX+nd);
	
}
