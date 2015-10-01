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
 * Module for the Limited Memory - Broyden-Fletcher-Goldfarb-Shanno (L-BFGS)
 * method for the elastic multiparameter inversion of
 * vp, vs, rho and lambda, mu, rho respectively
 * 
 * Daniel Koehn
 * ----------------------------------------------------------------------*/

#include "fd.h"

void LBFGS(float ** waveconv, float ** taper_coeff, int nsrc, float ** srcpos, int ** recpos, int ntr_glob, int iter, float C_vp, float ** gradp, int nfstart_jac,
	     float ** waveconv_u, float C_vs, float ** gradp_u, float ** waveconv_rho, float C_rho, float ** gradp_rho, float * y_LBFGS_vp, float * s_LBFGS_vp, float * rho_LBFGS, float * alpha_LBFGS, 
		 float * y_LBFGS_vs, float * s_LBFGS_vs, float * y_LBFGS_rho, float * s_LBFGS_rho, float ** ppi, float ** pu, float ** prho, int nxnyi){

	extern int NX, NY, IDX, IDY, SPATFILTER;
	extern int HESSIAN, INVMAT, SWS_TAPER_GRAD_VERT, SWS_TAPER_GRAD_HOR, SWS_TAPER_GRAD_SOURCES, SWS_TAPER_FILE;
	extern int POS[3], MYID;
	extern char JACOBIAN[STRING_SIZE];
	
	char jac[225], jac1[225];
	int i, j, k;
	float betaz, betan, gradplastiter, gradclastiter, betar, beta;
	float gamma_LBFGS, beta_LBFGS, Vp_sum, Vs_sum;
    float LBFGSTMP, LBFGSTMP1, modellastiter, norm_fac, norm_fac_u, norm_fac_rho;
    int NLBFGS, ki, itershift, iter1;
	extern FILE *FP;
	FILE *FP3, *FP4, *FP6, *FP5, *FP7;
	
	NLBFGS = 200;
    itershift = 1;

/* =================================================================================================================================================== */
/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT ZP ================================================================================== */
/* ===================================================================================================================================================== */


if((INVMAT==1)||(INVMAT==0)){

/* Normalization of the gradient */
/* ------------------------------- */
        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
                 waveconv[j][i] = C_vp*waveconv[j][i];
           }
        }

/* apply Hessian^-1 and save in gradp*/
if (SWS_TAPER_FILE){ 
         taper_grad(waveconv,taper_coeff,srcpos,nsrc,recpos,ntr_glob,4);
}

/* Normalize gradient to maximum value */
norm_fac=norm(waveconv,iter,1);

if(MYID==0){printf("norm_fac=%e \n",norm_fac);}

for (i=1;i<=NX;i=i+IDX){
    for (j=1;j<=NY;j=j+IDY){
	    gradp[j][i] = waveconv[j][i];
    }
}

/* apply spatial wavelength filter */
if(SPATFILTER==1){
	if (MYID==0){
   	fprintf(FP,"\n Spatial filter is applied to gradient (written by PE %d)\n",MYID);}
spat_filt(waveconv,iter,1);}


/* save gradient for output as inversion result */
if(iter==nfstart_jac){
	sprintf(jac,"%s_p_it%d.old.%i.%i",JACOBIAN,iter,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

    for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
                	fwrite(&waveconv[j][i],sizeof(float),1,FP3);
        }
    }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_it%d.old",JACOBIAN,iter);
	if (MYID==0) mergemod(jac,3);
}

/* =================================================================================================================================================== */
/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT Zs ================================================================================== */
/* ===================================================================================================================================================== */

if((INVMAT==3)||(INVMAT==0)){
	
/* Normalization of the gradient   */
/* ------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv_u[j][i] = C_vs * waveconv_u[j][i];
   }
}

/* apply Hessian^-1 and save in gradp*/
if (SWS_TAPER_FILE){ 
  taper_grad(waveconv_u,taper_coeff,srcpos,nsrc,recpos,ntr_glob,5);
}

/* Normalize gradient to maximum value */
norm_fac_u=norm(waveconv_u,iter,2);
if(MYID==0){printf("norm_fac_u=%e \n",norm_fac_u);}

for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
	  gradp_u[j][i] = waveconv_u[j][i];
   }
}

/* apply spatial wavelength filter */
if(SPATFILTER==1){
	if (MYID==0){
   	fprintf(FP,"\n Spatial filter is applied to gradient (written by PE %d)\n",MYID);}
spat_filt(waveconv_u,iter,2);}

/* save gradient for output as inversion result */
if(iter==nfstart_jac){
	sprintf(jac,"%s_p_u_it%d.old.%i.%i",JACOBIAN,iter,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        	for (i=1;i<=NX;i=i+IDX){
           	for (j=1;j<=NY;j=j+IDY){
                	fwrite(&waveconv_u[j][i],sizeof(float),1,FP3);
           	}
        	}
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_u_it%d.old",JACOBIAN,iter);
	if (MYID==0) mergemod(jac,3);
}

}

/* ===================================================================================================================================================== */
/* ===================================================== GRADIENT rho ================================================================================== */
/* ===================================================================================================================================================== */

if((INVMAT==2)||(INVMAT==0)){

/* Normalization of the gradient   */
/* ------------------------------- */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
      waveconv_rho[j][i] = C_rho * waveconv_rho[j][i];
   }
}

/* apply Hessian^-1 and save in gradp*/
if (SWS_TAPER_FILE){ 
  taper_grad(waveconv_rho,taper_coeff,srcpos,nsrc,recpos,ntr_glob,6);
}

/* Normalize gradient to maximum value */
norm_fac_rho=norm(waveconv_rho,iter,3);

if(MYID==0){printf("norm_fac_rho=%e \n",norm_fac_rho);}

for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
	  gradp_rho[j][i] = waveconv_rho[j][i];
   }
} 

/* apply spatial wavelength filter */
if(SPATFILTER==1){
	if (MYID==0){
   	fprintf(FP,"\n Spatial filter is applied to gradient (written by PE %d)\n",MYID);}
spat_filt(waveconv_rho,iter,3);}

/* save gradient for output as inversion result */
if(iter==nfstart_jac){
	sprintf(jac,"%s_p_rho_it%d.old.%i.%i",JACOBIAN,iter,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        	for (i=1;i<=NX;i=i+IDX){
           	for (j=1;j<=NY;j=j+IDY){
                	fwrite(&waveconv_rho[j][i],sizeof(float),1,FP3);
           	}
        	}
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_rho_it%d.old",JACOBIAN,iter);
	if (MYID==0) mergemod(jac,3);
}
}

/* calculate H^-1 * waveconv, using the L-BFGS method, if iter > 1 */
/* --------------------------------------------------------------------- */

if(iter>1){
   
   /* load old models and gradients - Vp */
   /* ---------------------------------- */

   sprintf(jac,"%s_p.old.%i.%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");

   sprintf(jac1,"%s_p_vp.old.%i.%i",JACOBIAN,POS[1],POS[2]);
   FP7=fopen(jac1,"rb");

   iter1 = iter-itershift; /* shift iter counter by 1 because L-BFGS method starts at iter > 1 */

   LBFGSTMP = 0.0;
   k=1+((iter1-1)*nxnyi);
   
   /*printf("k = %d \t iter1 = %d \t MYID = %d \n",k,iter1, MYID);*/
   
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS_vp[k] = waveconv[j][i]-gradplastiter;

	      fread(&modellastiter,sizeof(float),1,FP7);
          s_LBFGS_vp[k] = ppi[j][i]-modellastiter; 

	  k++;
       }
     }
     
     fclose(FP6);
     fclose(FP7);

   /* load old models and gradients - Rho */
   /* ----------------------------------- */

   sprintf(jac,"%s_p_rho.old.%i.%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");
   
   sprintf(jac1,"%s_p_mrho.old.%i.%i",JACOBIAN,POS[1],POS[2]);
   FP7=fopen(jac1,"rb");

   iter1 = iter-itershift; /* shift iter counter by 1 because L-BFGS method starts at iter > 1 */
   
   k=1+((iter1-1)*nxnyi);
   
   /*printf("k = %d \t iter1 = %d \t MYID = %d \n",k,iter1, MYID);*/
   
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS_rho[k] = waveconv_rho[j][i]-gradplastiter;

	      fread(&modellastiter,sizeof(float),1,FP7);
          s_LBFGS_rho[k] = prho[j][i]-modellastiter;
 

	  k++;
       }
     }
     
     fclose(FP6);
     fclose(FP7);
   
   /* load old models and gradients - Vs */
   sprintf(jac,"%s_p_u.old.%i.%i",JACOBIAN,POS[1],POS[2]);
   FP6=fopen(jac,"rb");

   sprintf(jac1,"%s_p_vs.old.%i.%i",JACOBIAN,POS[1],POS[2]);
   FP7=fopen(jac1,"rb");

   iter1 = iter-itershift; /* shift iter counter by 1 because L-BFGS method starts at iter > 1 */

   k=1+((iter1-1)*nxnyi);
   
   /*printf("k = %d \t iter1 = %d \t MYID = %d \n",k,iter1, MYID);*/
   
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate and save y, s at iteration step iter */
          fread(&gradplastiter,sizeof(float),1,FP6);
          y_LBFGS_vs[k] = waveconv_u[j][i]-gradplastiter;

    	  fread(&modellastiter,sizeof(float),1,FP7);
          s_LBFGS_vs[k] = pu[j][i]-modellastiter;  

	  k++;
       }
     }
     
     fclose(FP6);
     fclose(FP7);
     
     /*printf("rho_LBFGS_vp[iter1] = %e \t LBFGSTMP = %e \t Vp_sum = %e \t iter1 = %d \t MYID = %d \n",rho_LBFGS_vp[iter1],LBFGSTMP,Vp_sum,iter1,MYID);*/
     
	 /* calculate q=waveconv[j][i] */
	 for(k=iter1;k>=1;k--){

            ki=1+((k-1)*nxnyi);
       
	   /* calculate rho at iteration step iter for all parameter classes*/
	   LBFGSTMP1 = 0.0;
	   for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
          
              LBFGSTMP1 += y_LBFGS_vp[k] * s_LBFGS_vp[k]; 
	          LBFGSTMP1 += y_LBFGS_vs[ki] * s_LBFGS_vs[ki];
              /*LBFGSTMP1 += y_LBFGS_rho[k] * s_LBFGS_rho[k];*/  /*XX no density update !!! */ 
         
	          ki++;

           }
	   }

           /* Sum over Rho of all CPUs */
           Vp_sum = 0.0;
           MPI_Allreduce(&LBFGSTMP1,&Vp_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
           rho_LBFGS[k] = 1.0/Vp_sum;

	   ki=1+((k-1)*nxnyi);

	   LBFGSTMP = 0.0;
	   for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
          
		  /* calculate alpha */
		  LBFGSTMP += s_LBFGS_vp[ki]*waveconv[j][i];  
		  LBFGSTMP += s_LBFGS_vs[ki]*waveconv_u[j][i];
		  /*LBFGSTMP += s_LBFGS_rho[ki]*waveconv_rho[j][i]; */  /*XX no density update !!! */
	      ki++;
       }
	 }

     /* Sum over alpha of all CPUs */
     Vp_sum = 0.0;
     MPI_Allreduce(&LBFGSTMP,&Vp_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     alpha_LBFGS[k] = rho_LBFGS[k] * Vp_sum;

	 ki=1+((k-1)*nxnyi);
	  for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               
		  /* update gradient */
		  waveconv[j][i] = waveconv[j][i] - alpha_LBFGS[k] * y_LBFGS_vp[ki];
		  waveconv_u[j][i] = waveconv_u[j][i] - alpha_LBFGS[k] * y_LBFGS_vs[ki];
		  /*waveconv_rho[j][i] = waveconv_rho[j][i] - alpha_LBFGS[k] * y_LBFGS_rho[ki];*/  /*XX no density update !!! */

	      ki++;
       }
	 }

	 }

     /* calculate approximated Hessian */
     LBFGSTMP = 0.0; 
     LBFGSTMP1 = 0.0;
     k=1+((iter1-1)*nxnyi);
     for (i=1;i<=NX;i=i+IDX){
       for (j=1;j<=NY;j=j+IDY){
   	  
          /* calculate nominator and denominator of gamma */
	  	  LBFGSTMP += y_LBFGS_vp[k] * s_LBFGS_vp[k];
		  LBFGSTMP += y_LBFGS_vs[k] * s_LBFGS_vs[k];
		  /*LBFGSTMP += y_LBFGS_rho[k] * s_LBFGS_rho[k];*/  /*XX no density update !!! */

		  LBFGSTMP1 += y_LBFGS_vp[k] * y_LBFGS_vp[k];  
		  LBFGSTMP1 += y_LBFGS_vs[k] * y_LBFGS_vs[k];
		  /*LBFGSTMP1 += y_LBFGS_rho[k] * y_LBFGS_rho[k];*/  /*XX no density update !!! */

	      k++;
       }
	 }

	 /* Sum nominator and denominator of gamma of all CPUs */
         Vp_sum = 0.0;
         MPI_Allreduce(&LBFGSTMP,&Vp_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     
	     Vs_sum = 0.0;
         MPI_Allreduce(&LBFGSTMP1,&Vs_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

	 gamma_LBFGS = Vp_sum/Vs_sum;

	 /* scale gradient with approximated Hessian and calculate r */
	 for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
	           waveconv[j][i] = gamma_LBFGS * waveconv[j][i];
	           waveconv_u[j][i] = gamma_LBFGS * waveconv_u[j][i];
		   /*waveconv_rho[j][i] = gamma_LBFGS * waveconv_rho[j][i];*/  /*XX no density update !!! */
		   }
	 }

	 /* calculate H^-1 * waveconv[j][i] */
	 for(k=1;k<=iter1;k++){

           ki=1+((k-1)*nxnyi);
	   LBFGSTMP = 0.0;
	   for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               
		  /* calculate beta */
		  LBFGSTMP += y_LBFGS_vp[ki]*waveconv[j][i];
		  LBFGSTMP += y_LBFGS_vs[ki]*waveconv_u[j][i];
		  /*LBFGSTMP += y_LBFGS_rho[ki]*waveconv_rho[j][i];*/  /*XX no density update !!! */

	      ki++;
       }
	 }

     /* Sum over beta of all CPUs */
     Vp_sum = 0.0;
     MPI_Allreduce(&LBFGSTMP,&Vp_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     beta_LBFGS = rho_LBFGS[k] * Vp_sum;

	 ki=1+((k-1)*nxnyi);
	  for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               
		  /* update gradient */
		  waveconv[j][i] = waveconv[j][i] + s_LBFGS_vp[ki] * (alpha_LBFGS[k]-beta_LBFGS);  
		  waveconv_u[j][i] = waveconv_u[j][i] + s_LBFGS_vs[ki] * (alpha_LBFGS[k]-beta_LBFGS);
		  /*waveconv_rho[j][i] = waveconv_rho[j][i] + s_LBFGS_rho[ki] * (alpha_LBFGS[k]-beta_LBFGS);*/  /*XX no density update !!! */

	      ki++;
          }
	  }

	 } /* end of calculate H^-1 * waveconv[j][i] */


	 /* Denormalize Gradients */
	 for (i=1;i<=NX;i=i+IDX){
         for (j=1;j<=NY;j=j+IDY){
               
		  /* update gradient */
		  waveconv[j][i] = waveconv[j][i] * norm_fac;  
		  waveconv_u[j][i] = waveconv_u[j][i] * norm_fac_u;
		  /*waveconv_rho[j][i] = waveconv_rho[j][i] * norm_fac_rho;*/  /*XX no density update !!! */

         }
	 }

}

}

    /* save old model - Vp */
	sprintf(jac,"%s_p_vp.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               fwrite(&ppi[j][i],sizeof(float),1,FP3);
           }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge model file */ 
	sprintf(jac,"%s_p_vp.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);

	/* save old gradient */
	sprintf(jac,"%s_p.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&gradp[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* save H^-1 * g */
        sprintf(jac,"%s_c.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_c.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
                         



/* save old models Vs */
/* ------------------ */

    /* save old model */
	sprintf(jac,"%s_p_vs.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               fwrite(&pu[j][i],sizeof(float),1,FP3);
           }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge model file */ 
	sprintf(jac,"%s_p_vs.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);

	/* save old gradient */
	sprintf(jac,"%s_p_u.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&gradp_u[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_u.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* save H^-1 * g */
        sprintf(jac,"%s_c_u.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_u[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_c_u.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);


/* save old models Rho */
/* ------------------ */

	sprintf(jac,"%s_p_mrho.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
               fwrite(&prho[j][i],sizeof(float),1,FP3);
           }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge model file */ 
	sprintf(jac,"%s_p_mrho.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);

	/* save old gradient */
	sprintf(jac,"%s_p_rho.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");

        for (i=1;i<=NX;i=i+IDX){
            for (j=1;j<=NY;j=j+IDY){
                	fwrite(&gradp_rho[j][i],sizeof(float),1,FP3);
            }
        }
	
	fclose(FP3);

	MPI_Barrier(MPI_COMM_WORLD);
          
	/* merge gradient file */ 
	sprintf(jac,"%s_p_rho.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
	/* save H^-1 * g_rho */
        sprintf(jac,"%s_c_rho.old.%i.%i",JACOBIAN,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	
	for (i=1;i<=NX;i=i+IDX){   
           for (j=1;j<=NY;j=j+IDY){
                 fwrite(&waveconv_rho[j][i],sizeof(float),1,FP3);
	   }
        }
        
	fclose(FP3);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* merge gradient file */ 
	sprintf(jac,"%s_c_rho.old",JACOBIAN);
	if (MYID==0) mergemod(jac,3);
	
}
