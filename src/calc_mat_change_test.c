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
 *   calculate test step length for material parameter update
 *   
 *   Daniel Koehn
 *   last update 9.11.2007
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
void calc_mat_change_test(float  **  waveconv, float  **  waveconv_rho, float  **  waveconv_u, float  **  rho, float  **  rhonp1, float **  pi, float **  pinp1, float **  u, float **  unp1, int iter, 
                          int epstest, int INVMAT, float eps_scale, int itest, int nfstart, float ** u_start, float ** pi_start, float ** rho_start){


	/*--------------------------------------------------------------------------*/
	FILE *FP1;
	/* extern variables */
	extern float DH, DT, VP_VS_RATIO;
	extern float EPSILON, EPSILON_u, EPSILON_rho, MUN;
	extern int NX, NY, NXG, NYG,  POS[3], MYID, INVMAT1;
	
	extern int INV_RHO_ITER, INV_VP_ITER, INV_VS_ITER;
	
	extern char INV_MODELFILE[STRING_SIZE];
	
	extern float VPUPPERLIM, VPLOWERLIM, VSUPPERLIM, VSLOWERLIM, RHOUPPERLIM, RHOLOWERLIM;
	
	extern int S, ACOUSTIC;
	extern float S_VS, S_VP, S_RHO;

	/* local variables */

	float Rho, Vp, Vs, Vsnp1, Vpnp1, x, y, undf, r, pi0, K, mu, Zp, Zs, pro_vs, pro_vp, pro_rho;
	float dpi, pimax, rhomax, umax, gradmax, gradmax_rho, gradmax_u, epsilon1, pimaxr, gradmaxr, gradmaxr_u, umaxr, gradmaxr_rho, rhomaxr;
	int i, j, ii, jj, testuplow, pr=0;
	char modfile[STRING_SIZE];
	
	
	/*int RHO_VS=0;
	int f=0;*/
	
	
	/* invert for Zp and Zs */
	/* ------------------------------------------------------------------------------------ */	
	
	/* find maximum of Zp and gradient waveconv */	
	pimax = 0.0;
	gradmax = 0.0;
	
	for (i=1;i<=NX;i++){
	for (j=1;j<=NY;j++){
		
		Zp = pi[j][i];
		
		if(Zp>pimax)
			pimax=Zp;
		
		if((i*j == 1) || (gradmax == 0.0)) {
				gradmax = fabs(waveconv[j][i]);
		} else {
			if(fabs(waveconv[j][i]) > gradmax)
				gradmax = fabs(waveconv[j][i]);
		}
	}
	}
	
	/* find maximum of Zs and gradient waveconv_u */
	if(!ACOUSTIC){
		umax = 0.0;
		gradmax_u = 0.0;
		
		for (i=1;i<=NX;i++){
		for (j=1;j<=NY;j++){
			
			Zs = u[j][i];
			
			if(Zs>umax)
				umax=Zs;
			
			if((i*j == 1) || (gradmax_u == 0.0)) {
				gradmax_u = fabs(waveconv_u[j][i]);
			} else {
				if(fabs(waveconv_u[j][i]) > gradmax_u)
					gradmax_u = fabs(waveconv_u[j][i]);		                               
			}
		}
		}
	}
	
	/* find maximum of rho and gradient waveconv_rho */
	rhomax = 0.0;
	gradmax_rho = 0.0;
	
	for (i=1;i<=NX;i++){
	for (j=1;j<=NY;j++){
		
		if(rho[j][i]>rhomax){rhomax=rho[j][i];}
		
		if((i*j == 1) || (gradmax_rho == 0.0)) {
		        gradmax_rho = fabs(waveconv_rho[j][i]);
		} else {
			if(fabs(waveconv_rho[j][i]) > gradmax_rho)
				gradmax_rho = fabs(waveconv_rho[j][i]);
		}
		                               
	}
	}

		
	/* calculate scaling factor for the gradients */
	/* --------------------------------------------- */
	
	/* parameter 1 */
	MPI_Allreduce(&pimax,  &pimaxr,  1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&gradmax,&gradmaxr,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);

	EPSILON = eps_scale * (pimaxr/gradmaxr);
	if (iter<INV_VP_ITER){EPSILON = 0.0;}
	epsilon1=EPSILON;
	MPI_Allreduce(&EPSILON,&epsilon1,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);

	if (MYID==0)  EPSILON=epsilon1;

	exchange_par();
	
	if(!ACOUSTIC){
		/* parameter 2 */
		MPI_Allreduce(&umax,&umaxr,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
		MPI_Allreduce(&gradmax_u,&gradmaxr_u,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);

		EPSILON_u = eps_scale * (umaxr/gradmaxr_u);
		if (iter<INV_VS_ITER){EPSILON_u = 0.0;}
		epsilon1=EPSILON_u;
		MPI_Allreduce(&EPSILON_u,&epsilon1,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);

		if (MYID==0)  EPSILON_u=epsilon1;

		exchange_par();	
	}
	
	/* parameter 3 */
	MPI_Allreduce(&rhomax,&rhomaxr,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&gradmax_rho,&gradmaxr_rho,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);

	EPSILON_rho = eps_scale * (rhomaxr/gradmaxr_rho);
	if (iter<INV_RHO_ITER){EPSILON_rho = 0.0;}
	epsilon1=EPSILON_rho;
	MPI_Allreduce(&EPSILON_rho,&epsilon1,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);

	if (MYID==0)  EPSILON_rho=epsilon1;

	exchange_par();	
	
	
        if(MYID==0){
		printf("MYID = %d \t pimaxr = %e \t gradmaxr = %e \n",MYID,pimaxr,gradmaxr);
		printf("MYID = %d \t EPSILON = %e \n",MYID,EPSILON);
		if(!ACOUSTIC)
			printf("MYID = %d \t umaxr = %e \t gradmaxr_u = %e \n",MYID,umaxr,gradmaxr_u);
		printf("MYID = %d \t rhomaxr = %e \t gradmaxr_rho = %e \n",MYID,rhomaxr,gradmaxr_rho);
	}
	
	
	/* loop over local grid */
	for (i=1;i<=NX;i++){
	for (j=1;j<=NY;j++){		
		
		/* update lambda, mu, rho */
		if((INVMAT1==3) || (INVMAT1==1)){
			
			testuplow=0;
			
			pinp1[j][i] = pi[j][i] - EPSILON*waveconv[j][i];
			if(!ACOUSTIC)
				unp1[j][i] = u[j][i] - EPSILON_u*waveconv_u[j][i];
			rhonp1[j][i] = rho[j][i] - EPSILON_rho*waveconv_rho[j][i];
			
			/*rhonp1[j][i] = rho[j][i];
			pinp1[j][i] = pi[j][i];
			unp1[j][i] = u[j][i];*/
			
			
			if(pinp1[j][i]<VPLOWERLIM){
				pinp1[j][i] = VPLOWERLIM;
			}
			
			if(pinp1[j][i]>VPUPPERLIM){
				pinp1[j][i] = VPUPPERLIM;
			}
			
			
			if(!ACOUSTIC){
				if((unp1[j][i]<VSLOWERLIM)&&(unp1[j][i]>1e-6)){
					unp1[j][i] = VSLOWERLIM;
				}
				if(unp1[j][i]>VSUPPERLIM){
					unp1[j][i] = u[j][i];
				}
				
				/* checking poisson ratio */
				if(VP_VS_RATIO>1){
					if((pinp1[j][i]/unp1[j][i])<VP_VS_RATIO){
						pinp1[j][i]=unp1[j][i]*VP_VS_RATIO;
						pr=1;
					}
				}/* if(VP_VS_RATIO>1) */
			}
			
			if(S==1){
				if(!ACOUSTIC){
					/* limited update for Vs */
					pro_vs=((u_start[j][i]-unp1[j][i])/u_start[j][i])*100.0;
					if(((fabs(pro_vs))>S_VS)&&(pro_vs>0.0))
						unp1[j][i] = u_start[j][i] - u_start[j][i] * S_VS / 100.0;
					if(((fabs(pro_vs))>S_VS)&&(pro_vs<0.0))
						unp1[j][i] = u_start[j][i] + u_start[j][i] * S_VS / 100.0;
				}
				
				/* limited update for Vp */
				pro_vp=((pi_start[j][i]-pinp1[j][i])/pi_start[j][i])*100.0;
				if(((fabs(pro_vp))>S_VP)&&(pro_vp>0.0))
					pinp1[j][i] = pi_start[j][i] - pi_start[j][i] * S_VP / 100.0;
				if(((fabs(pro_vp))>S_VP)&&(pro_vp<0.0))
					pinp1[j][i] = pi_start[j][i] + pi_start[j][i] * S_VP / 100.0;
				
				/* limited update for Rho */
				pro_rho=((rho_start[j][i]-rhonp1[j][i])/rho_start[j][i])*100.0;
				if(((fabs(pro_rho))>S_RHO)&&(pro_rho>0.0))
					rhonp1[j][i] = rho_start[j][i] - rho_start[j][i] * S_RHO / 100.0;
				if(((fabs(pro_rho))>S_RHO)&&(pro_rho<0.0))
					rhonp1[j][i] = rho_start[j][i] + rho_start[j][i] * S_RHO / 100.0;
			}
			
			/*if(RHO_VS==1){
			/* Coupling Rho with Vs about Gardner Relationship */
			/*rhonp1[j][i] = pow(unp1[j][i],0.295)*0.3*1000; 
			f=1;
			}*/
			
			if((rhonp1[j][i]<RHOLOWERLIM) && (INV_RHO_ITER < iter)){
				rhonp1[j][i] = rho[j][i];
			}
			if((rhonp1[j][i]>RHOUPPERLIM) && (INV_RHO_ITER < iter)){
				rhonp1[j][i] = rho[j][i];
			}
			
			
			/* None of these parameters should be smaller than zero */
			if(pinp1[j][i]<0.0){
				pinp1[j][i] = pi[j][i];
			}
			if(!ACOUSTIC)
				if(unp1[j][i]<0.0){
					unp1[j][i] = u[j][i];
				} 
			if(rhonp1[j][i]<0.0){
				rhonp1[j][i] = rho[j][i];
			}
			
			if(itest==0){
				rho[j][i] = rhonp1[j][i];
				pi[j][i] = pinp1[j][i];
				if(!ACOUSTIC)
					u[j][i] = unp1[j][i];
				
			}
		}
	}
	}
	
	if(!ACOUSTIC)
		if((MYID==0)&&(pr==1))printf("\nThe Vp/Vs-ratio of %4.2f is violated. P-wave velocity will be increased that the Vp/Vs-ratio is at least %4.2f \n",VP_VS_RATIO,VP_VS_RATIO);
	/*if((MYID==0)&&(f==1))printf("\n RHO_VS == %d: Rho is coupled to Vp about Gardner Relationship \n",RHO_VS); */
	
	
	if(itest==0){
		/*sprintf(modfile,"model/waveform_test_model_vp_it_%d.bin",iter);*/
		/*sprintf(modfile,"model/waveform_test_model_vp.bin");*/
		sprintf(modfile,"%s_vp.bin",INV_MODELFILE);
		writemod(modfile,pinp1,3);
		
		MPI_Barrier(MPI_COMM_WORLD);

		if (MYID==0) mergemod(modfile,3);

		if(!ACOUSTIC){
			/*sprintf(modfile,"model/waveform_test_model_vs_it_%d.bin",iter);*/
			/*sprintf(modfile,"model/waveform_test_model_vs.bin");*/
			sprintf(modfile,"%s_vs.bin",INV_MODELFILE);
			
			writemod(modfile,unp1,3);
			
			MPI_Barrier(MPI_COMM_WORLD);

			if (MYID==0) mergemod(modfile,3);
		}

		/*sprintf(modfile,"model/waveform_test_model_rho_it_%d.bin",iter);*/
		/*sprintf(modfile,"model/waveform_test_model_rho.bin");*/
		sprintf(modfile,"%s_rho.bin",INV_MODELFILE);
		writemod(modfile,rho,3);
		
		MPI_Barrier(MPI_COMM_WORLD);

		if (MYID==0) mergemod(modfile,3);
	}
        
        if((itest==0)&&(iter==nfstart)){
		/*sprintf(modfile,"model/waveform_test_model_vp_it_%d.bin",iter);*/
		sprintf(modfile,"%s_vp_it%d.bin",INV_MODELFILE,iter);
		writemod(modfile,pinp1,3);
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (MYID==0) mergemod(modfile,3);
		MPI_Barrier(MPI_COMM_WORLD); 
		sprintf(modfile,"%s_vp_it%d.bin.%i.%i",INV_MODELFILE,iter,POS[1],POS[2]);
		remove(modfile);
		
		if(!ACOUSTIC){
			/*sprintf(modfile,"model/waveform_test_model_vs_it_%d.bin",iter);*/
			sprintf(modfile,"%s_vs_it%d.bin",INV_MODELFILE,iter);
			writemod(modfile,unp1,3);
			MPI_Barrier(MPI_COMM_WORLD);
			
			if (MYID==0) mergemod(modfile,3);
			MPI_Barrier(MPI_COMM_WORLD); 
			sprintf(modfile,"%s_vs_it%d.bin.%i.%i",INV_MODELFILE,iter,POS[1],POS[2]);
			remove(modfile);
		}
		
		/*sprintf(modfile,"model/waveform_test_model_rho_it_%d.bin",iter);*/
		sprintf(modfile,"%s_rho_it%d.bin",INV_MODELFILE,iter);
		writemod(modfile,rho,3);
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (MYID==0) mergemod(modfile,3);
		MPI_Barrier(MPI_COMM_WORLD); 
		sprintf(modfile,"%s_rho_it%d.bin.%i.%i",INV_MODELFILE,iter,POS[1],POS[2]);
		remove(modfile);
	}
}