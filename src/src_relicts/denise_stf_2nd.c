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
 *  DENISE: 2D elastic time domain FWT Code 
 *
 * 
 *  ---------------------------------------------------------------------------------------*/

#include "fd.h"           /* general include file for viscoelastic FD programs */

#include "globvar.h"      /* definition of global variables  */
#include "cseife.h"

#include "stfinv/stfinv.h" /* libstfinv - inversion for source time function */

int main(int argc, char **argv){
/* variables in main */
int ns, nseismograms=0, nt, nd, fdo3, j, i, ii, jj, shotid, recid, k, nc, iter, h, infoout, SHOTINC, TIMEWIN, test_eps, lq, iq, jq, hin, hin1, s=0;
int DTINV, nxny, nxnyi, imat, imat1, imat2, IDXI, IDYI, hi, NTST, NTSTI, partest, FREQFILT;
int lsnap, nsnap=0, lsamp=0, buffsize, invtime, invtimer, sws, swstestshot, snapseis, snapseis1, PML;
int ntr=0, ntr_loc=0, ntr_glob=0, nsrc=0, nsrc_loc=0, nsrc_glob=0, ishot, irec, nshots=0, nshots1, Lcount, LNORM, itest, Lcountsum, itestshot;

float pum, ppim, ppim1, ppim2, thetaf, thetab, e33, e33b, e11, e11b, muss, lamss; 
float memdyn, memmodel, memseismograms, membuffer, memtotal, dngn, fphi, sum, avggrad, beta, betan, betaz, betaLog, betaVp, betaVs, betarho, eps_scale, L2old;
float fac1, fac2, wavefor, waverecipro, dump, dump1, epsilon, gradsign, mun, eps1, gradplastiter, gradglastiter, gradclastiter, betar, sig_max, sig_max1;
float signL1, RMS, opteps_vp, opteps_vs, opteps_rho, Vs, Vp, Vp_avg, C_vp, Vs_avg, C_vs, Cd, rho_avg, C_rho, Vs_sum, Vp_sum, rho_sum, Zp, Zs;
float freqshift, dfreqshift, memfwt, memfwt1, memfwtdata;
char *buff_addr, ext[10], *fileinp;
char wave_forward[225], wave_recipro[225], wave_conv[225], jac[225], jac2[225], jacsum[225], dwavelet[225], vyf[STRING_SIZE];

double time1, time2, time3, time4, time5, time6, time7, time8,
	time_av_v_update=0.0, time_av_s_update=0.0, time_av_v_exchange=0.0, 
	time_av_s_exchange=0.0, time_av_timestep=0.0;
	
float L2, L2sum, L2_all_shots, L2sum_all_shots, *L2t, alphanomsum, alphanom, alphadenomsum, alphadenom, scaleamp ,sdummy, lamr; 

float energy, energy_sum, energy_all_shots, energy_sum_all_shots;	

float  **  psxx, **  psxy, **  psyy, ** ux, ** uy, ** uxy, ** uyx, ** Vp0, ** uttx, ** utty, ** Vs0, ** Rho0;
float  **  pvx, **  pvy, **waveconv, **waveconv_lam, **waveconv_mu, **waveconv_rho, **waveconv_rho_s, **waveconv_u, **waveconvtmp, **wcpart, **wavejac;
float **waveconv_shot, **waveconv_u_shot, **waveconv_rho_shot;
float  **  pvxp1, **  pvyp1, **  pvxm1, **  pvym1;
float ** gradg, ** gradp,** gradg_rho, ** gradp_rho, ** gradg_u, ** gradp_u;
float  **  prho,**  prhonp1, **prip=NULL, **prjp=NULL, **pripnp1=NULL, **prjpnp1=NULL, **  ppi, **  pu, **  punp1, **  puipjp, **  ppinp1;
float  **  vpmat, *forward_prop_x, *forward_prop_y, *forward_prop_rho_x, *forward_prop_u, *forward_prop_rho_y;

float  ** sectionvx=NULL, ** sectionvy=NULL, ** sectionp=NULL, ** sectionpnp1=NULL,
	** sectioncurl=NULL, ** sectiondiv=NULL, ** sectionvxdata=NULL, ** sectionvxdiff=NULL, ** sectionvxdiffold=NULL, ** sectionvydiffold=NULL,
	** sectionvydiff=NULL, ** sectionvydata=NULL, ** sectionpn=NULL, ** sectionread=NULL, ** sectionvy_conv=NULL, ** sectionvy_obs=NULL,
	** sectionvy_syn=NULL, * source_time_function=NULL;
float  **  absorb_coeff, ** taper_coeff, * epst1, * epst2,  * epst3, * picked_times;
float  ** srcpos=NULL, **srcpos_loc=NULL, ** srcpos1=NULL, **srcpos_loc_back=NULL, ** signals=NULL, ** signals_rec=NULL, *hc=NULL, ** dsignals=NULL, ** signals_stf=NULL;
int   ** recpos=NULL, ** recpos_loc=NULL;
int   ** tracekill=NULL, TRKILL, DTRKILL;

float ** bufferlef_to_rig,  ** bufferrig_to_lef, ** buffertop_to_bot, ** bufferbot_to_top; 

/* PML variables */
float * d_x, * K_x, * alpha_prime_x, * a_x, * b_x, * d_x_half, * K_x_half, * alpha_prime_x_half, * a_x_half, * b_x_half, * d_y, * K_y, * alpha_prime_y, * a_y, * b_y, * d_y_half, * K_y_half, * alpha_prime_y_half, * a_y_half, * b_y_half;
float ** psi_sxx_x, ** psi_syy_y, ** psi_sxy_y, ** psi_sxy_x, ** psi_vxx, ** psi_vyy, ** psi_vxy, ** psi_vyx, ** psi_vxxs;

/* Variables for step length calculation */
int step1, step2, itests, iteste, stepmax, countstep;
float scalefac;

/* Variables for Pseudo-Hessian calculation */
int RECINC, ntr1;
float * jac_rho, * jac_u, * jac_lam_x, * jac_lam_y;
float * temp_TS, * temp_TS1, * temp_TS2, * temp_TS3, * temp_TS4, * temp_TS5, * temp_conv, * temp_conv1, * temp_conv2;
float TSHIFT_back, temp_hess, mulamratio;
float ** hessian, ** hessian_u, ** hessian_rho;
int QUELLART_OLD;

/* Variables of the L-BFGS method */
float *** y_LBFGS_vp, *** s_LBFGS_vp, * rho_LBFGS, * alpha_LBFGS; 
float *** y_LBFGS_vs, *** s_LBFGS_vs;
float *** y_LBFGS_rho, *** s_LBFGS_rho;
int NLBFGS;
float * rho_LBFGS_vp, * rho_LBFGS_vs, * alpha_LBFGS_vp, * alpha_LBFGS_vs;

int * recswitch=NULL;
float ** fulldata=NULL, ** fulldata_vx=NULL, ** fulldata_vy=NULL;

FILE *fprec, *FP2, *FP3, *FP4, *FP5, *FPL2, *FP6, *FP7;
	
MPI_Request *req_send, *req_rec;
MPI_Status  *send_statuses, *rec_statuses;

/* Initialize MPI environment */
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&NP);
MPI_Comm_rank(MPI_COMM_WORLD,&MYID);

setvbuf(stdout, NULL, _IONBF, 0);

if (MYID == 0){
   time1=MPI_Wtime(); 
   clock();
  }
		

/* print program name, version etc to stdout*/
if (MYID == 0) info(stdout);

  
/* read parameters from parameter-file (stdin) */
fileinp=argv[1];
FP=fopen(fileinp,"r");
read_par(FP);
exchange_par();


if (MYID == 0) note(stdout);


 
/* open log-file (each PE is using different file) */
/*	fp=stdout; */
sprintf(ext,".%i",MYID);  
strcat(LOG_FILE,ext);
if ((MYID==0) && (LOG==1)) FP=stdout;
else FP=fopen(LOG_FILE,"w");
fprintf(FP," This is the log-file generated by PE %d \n\n",MYID);
	
/* domain decomposition */
initproc();


NT=iround(TIME/DT);  	  /* number of timesteps */
ns=iround(NT/NDT);           /* number of samples per trace */
lsnap=iround(TSNAP1/DT);      /* first snapshot at this timestep */
lsamp=NDT;


/* output of parameters to log-file or stdout */
if (MYID==0)
   write_par(FP);



	
/* NXG, NYG denote size of the entire (global) grid */
NXG=NX;
NYG=NY;

/* In the following, NX and NY denote size of the local grid ! */
NX = IENDX;
NY = IENDY;


if (SEISMO){
   recpos=receiver(FP, &ntr);
   recswitch = ivector(1,ntr);
   recpos_loc = splitrec(recpos,&ntr_loc, ntr, recswitch);
   ntr_glob=ntr;
   ntr=ntr_loc;
}

fulldata = matrix(1,ntr_glob,1,NT);
fulldata_vx = matrix(1,ntr_glob,1,NT);
fulldata_vy = matrix(1,ntr_glob,1,NT);

/* estimate memory requirement of the variables in megabytes*/
	
switch (SEISMO){
case 1 : /* particle velocities only */
	nseismograms=2;	
	break;	
case 2 : /* pressure only */
	nseismograms=1;	
	break;	
case 3 : /* curl and div only */
	nseismograms=2;		
	break;	
case 4 : /* everything */
	nseismograms=5;		
	break;
}	
	
/* use only every DTINV time sample for the inversion */
DTINV=1;

/* save every IDXI and IDYI spatial point during the forward modelling */
IDXI=1;
IDYI=1;

/*allocate memory for dynamic, static and buffer arrays */
fac1=(NX+FDORDER)*(NY+FDORDER);
fac2=sizeof(float)*pow(2.0,-20.0);

nd = FDORDER/2 + 1;
fdo3 = 2*nd;

memdyn=5.0*fac1*fac2;
memmodel=6.0*fac1*fac2;
memseismograms=nseismograms*ntr*ns*fac2;

memfwt=5.0*((NX/IDXI)+FDORDER)*((NY/IDYI)+FDORDER)*(NT/DTINV)*fac2;
memfwt1=20.0*NX*NY*fac2;
memfwtdata=6.0*ntr*ns*fac2;

membuffer=2.0*fdo3*(NY+NX)*fac2;
buffsize=2.0*2.0*fdo3*(NX +NY)*sizeof(MPI_FLOAT);
memtotal=memdyn+memmodel+memseismograms+memfwt+memfwt1+memfwtdata+membuffer+(buffsize*pow(2.0,-20.0));


if (MYID==0){
   fprintf(FP,"\n **Message from main (printed by PE %d):\n",MYID);
   fprintf(FP," Size of local grids: NX=%d \t NY=%d\n",NX,NY);
   fprintf(FP," Each process is now trying to allocate memory for:\n");
   fprintf(FP," Dynamic variables: \t\t %6.2f MB\n", memdyn);
   fprintf(FP," Static variables: \t\t %6.2f MB\n", memmodel);
   fprintf(FP," Seismograms: \t\t\t %6.2f MB\n", memseismograms);
   fprintf(FP," Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
   fprintf(FP," Network Buffer for MPI_Bsend: \t %6.2f MB\n", buffsize*pow(2.0,-20.0));
   fprintf(FP," ------------------------------------------------ \n");
   fprintf(FP," Total memory required: \t %6.2f MB.\n\n", memtotal);
   }


/* allocate buffer for buffering messages */
buff_addr=malloc(buffsize);
if (!buff_addr) err("allocation failure for buffer for MPI_Bsend !");
MPI_Buffer_attach(buff_addr,buffsize);

/* allocation for request and status arrays */
req_send=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
req_rec=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
send_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));
rec_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));

/* memory allocation for dynamic (wavefield) arrays */
psxx =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
psxy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
psyy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvxp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvyp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvxm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pvym1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
ux   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
uy   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
uxy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
uyx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
uttx   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
utty   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
Vp0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
Vs0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
Rho0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

waveconv = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_lam = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconvtmp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
wcpart = matrix(1,3,1,3);
wavejac = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	
/* memory allocation for static (model) arrays */
prho =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
prhonp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
prip =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
prjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pripnp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
prjpnp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
ppi  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
ppinp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pu   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
punp1   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
vpmat   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
puipjp   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);

/* Time windowing  */
/* TIMEWIN==0 OFF */
/* TIMEWIN==1 ON  */
/* TIMEWIN==2 Time Window with STA/LTA Picker*/
TIMEWIN=0;

/* SPATFILTER = 1 activate spatial wavelengthfilter */
/* SPATFILTER = 0 deactivate spatial wavelengthfilter */
/*SPATFILTER=0;*/

/* FREQFILT = 1 activate frequency filter */
/* FREQFILT = 0 deactivate frequency filter */
FREQFILT=0;

lamr=0.01; /* Marquardt factor */

/*nf=4;
nfstart=4;*/

NTST=20;
NTSTI=NTST/DTINV;

nxny=NX*NY;
nxnyi=(NX/IDXI)*(NY/IDYI);

/* Parameters for step length calculations */
stepmax = 6; /* number of maximum misfit calculations/steplength 2/3*/
scalefac = 2.0; /* scale factor for the step length */

/* Variables for the L-BFGS method */
/*if(GRAD_METHOD==2){
  NLBFGS = 200;
  y_LBFGS_vp  =  vector(1,nxnyi*NLBFGS);
  s_LBFGS_vp  =  vector(1,nxnyi*NLBFGS);
  rho_LBFGS = vector(1,NLBFGS);
  alpha_LBFGS = vector(1,NLBFGS);

  y_LBFGS_vs  =  vector(1,nxnyi*NLBFGS);
  s_LBFGS_vs  =  vector(1,nxnyi*NLBFGS);

  y_LBFGS_rho  =  vector(1,nxnyi*NLBFGS);
  s_LBFGS_rho  =  vector(1,nxnyi*NLBFGS);
}*/

if(GRAD_METHOD==3){
  NLBFGS = 200;
  y_LBFGS_vp  =  f3tensor(1,NY,1,NX,1,NLBFGS);
  s_LBFGS_vp  =  f3tensor(1,NY,1,NX,1,NLBFGS);

  y_LBFGS_vs  =  f3tensor(1,NY,1,NX,1,NLBFGS);
  s_LBFGS_vs  =  f3tensor(1,NY,1,NX,1,NLBFGS);

  y_LBFGS_rho  =  f3tensor(1,NY,1,NX,1,NLBFGS);
  s_LBFGS_rho  =  f3tensor(1,NY,1,NX,1,NLBFGS);

  rho_LBFGS_vp = vector(1,NLBFGS);
  rho_LBFGS_vs = vector(1,NLBFGS);
  alpha_LBFGS_vp = vector(1,NLBFGS);
  alpha_LBFGS_vs = vector(1,NLBFGS);
  
}

if((INVMAT==1)||(INVMAT==0)){
/*forward_prop_x =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
forward_prop_y =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);*/

forward_prop_x =  vector(1,nxnyi*((NT/DTINV)));
forward_prop_y =  vector(1,nxnyi*((NT/DTINV)));

gradg = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
gradp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
}

if((INVMAT==2)||(INVMAT==0)){
/*forward_prop_rho_x =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
forward_prop_rho_y =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);*/

forward_prop_rho_x =  vector(1,nxnyi*((NT/DTINV)));
forward_prop_rho_y =  vector(1,nxnyi*((NT/DTINV)));

gradg_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
gradp_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_rho_s = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_rho_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
}

if((INVMAT==3)||(INVMAT==0)){
/*forward_prop_u =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);*/

forward_prop_u =  vector(1,nxnyi*((NT/DTINV)));

gradg_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
gradp_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
waveconv_u_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
}

if(HESSIAN){
	jac_rho =  vector(1,nxnyi*((NT/DTINV)));
	jac_u =  vector(1,nxnyi*((NT/DTINV)));
	jac_lam_x =  vector(1,nxnyi*((NT/DTINV)));
	jac_lam_y =  vector(1,nxnyi*((NT/DTINV)));
	temp_TS =  vector(1,(NT/DTINV));
	temp_TS1 =  vector(1,(NT/DTINV));
	temp_TS2 =  vector(1,(NT/DTINV));
	temp_TS3 =  vector(1,(NT/DTINV));
	temp_TS4 =  vector(1,(NT/DTINV));
	temp_TS5 =  vector(1,(NT/DTINV));
	temp_conv =  vector(1,(NT/DTINV));
	temp_conv1 =  vector(1,(NT/DTINV));
	temp_conv2 =  vector(1,(NT/DTINV));
	hessian = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	hessian_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	hessian_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
	TSHIFT_back = 0.0;
}

if(FW>0){

  d_x = vector(1,2*FW);
  K_x = vector(1,2*FW);
  alpha_prime_x = vector(1,2*FW);
  a_x = vector(1,2*FW);
  b_x = vector(1,2*FW);
  
  d_x_half = vector(1,2*FW);
  K_x_half = vector(1,2*FW);
  alpha_prime_x_half = vector(1,2*FW);
  a_x_half = vector(1,2*FW);
  b_x_half = vector(1,2*FW);

  d_y = vector(1,2*FW);
  K_y = vector(1,2*FW);
  alpha_prime_y = vector(1,2*FW);
  a_y = vector(1,2*FW);
  b_y = vector(1,2*FW);
  
  d_y_half = vector(1,2*FW);
  K_y_half = vector(1,2*FW);
  alpha_prime_y_half = vector(1,2*FW);
  a_y_half = vector(1,2*FW);
  b_y_half = vector(1,2*FW);

  psi_sxx_x =  matrix(1,NY,1,2*FW); 
  psi_syy_y =  matrix(1,2*FW,1,NX);
  psi_sxy_y =  matrix(1,2*FW,1,NX);
  psi_sxy_x =  matrix(1,NY,1,2*FW);
  psi_vxx   =  matrix(1,NY,1,2*FW);
  psi_vxxs  =  matrix(1,NY,1,2*FW); 
  psi_vyy   =  matrix(1,2*FW,1,NX);
  psi_vxy   =  matrix(1,2*FW,1,NX);
  psi_vyx   =  matrix(1,NY,1,2*FW);
   
}

taper_coeff=  matrix(1,NY,1,NX);



/* memory allocation for buffer arrays in which the wavefield
	   information which is exchanged between neighbouring PEs is stored */
bufferlef_to_rig = matrix(1,NY,1,fdo3);
bufferrig_to_lef = matrix(1,NY,1,fdo3);
buffertop_to_bot = matrix(1,NX,1,fdo3);
bufferbot_to_top = matrix(1,NX,1,fdo3);



if (ntr>0){
	switch (SEISMO){
	case 1 : /* particle velocities only */
		sectionvx=matrix(1,ntr,1,ns);
		sectionvy=matrix(1,ntr,1,ns);	
		break;	
	case 2 : /* pressure only */
		sectionp=matrix(1,ntr,1,ns);
		sectionpnp1=matrix(1,ntr,1,ns);
		sectionpn=matrix(1,ntr,1,ns);
		break;	
	case 3 : /* curl and div only */
		sectioncurl=matrix(1,ntr,1,ns);
		sectiondiv=matrix(1,ntr,1,ns);	
		break;	
	case 4 : /* everything */
		sectionvx=matrix(1,ntr,1,ns);
		sectionvy=matrix(1,ntr,1,ns);	
		sectioncurl=matrix(1,ntr,1,ns);
		sectiondiv=matrix(1,ntr,1,ns);		
		sectionp=matrix(1,ntr,1,ns);
		break;	
	
	}
}	

/* Memory for seismic data */
sectionvxdata=matrix(1,ntr,1,ns);
sectionread=matrix(1,ntr_glob,1,ns);
sectionvxdiff=matrix(1,ntr,1,ns);
sectionvydata=matrix(1,ntr,1,ns);
sectionvydiff=matrix(1,ntr,1,ns);
sectionvxdiffold=matrix(1,ntr,1,ns);
sectionvydiffold=matrix(1,ntr,1,ns);


/* Memory for inversion for source time function */
sectionvy_conv=matrix(1,ntr_glob,1,NT);
sectionvy_obs=matrix(1,ntr_glob,1,NT);
sectionvy_syn=matrix(1,ntr,1,NT);
source_time_function = vector(1,NT);

 /* memory for source position definition */
srcpos1=fmatrix(1,6,1,1);

/* memory of L2 norm */
L2t = vector(1,4);
epst1 = vector(1,3);
epst2 = vector(1,3);
epst3 = vector(1,3);
picked_times = vector(1,ntr);
	
fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);

		
/* Holberg coefficients for FD operators*/
hc = holbergcoeff();

MPI_Barrier(MPI_COMM_WORLD);


/* create model grids */



if (READMOD) readmod_elastic(prho,ppi,pu);
    else model_elastic(prho,ppi,pu);


/* check if the FD run will be stable and free of numerical dispersion */
checkfd_ssg_elastic(FP,prho,ppi,pu,hc);


/* calculate damping coefficients for CPMLs*/
if(FW>0){PML_pro(d_x, K_x, alpha_prime_x, a_x, b_x, d_x_half, K_x_half, alpha_prime_x_half, a_x_half, b_x_half, 
               d_y, K_y, alpha_prime_y, a_y, b_y, d_y_half, K_y_half, alpha_prime_y_half, a_y_half, b_y_half);}

MPI_Barrier(MPI_COMM_WORLD);

if (CHECKPTREAD){
	if (MYID==0){
		time3=MPI_Wtime();
 		fprintf(FP," Reading wavefield from check-point file %s \n",CHECKPTFILE);	
	}
	
	read_checkpoint(-1, NX+2, -1, NY+2, pvx, pvy, psxx, psyy, psxy);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0){
		time4=MPI_Wtime();
      		fprintf(FP," finished (real time: %4.2f s).\n",time4-time3);
	}
}
      

/* comunication initialisation for persistent communication */
comm_ini(bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec);


/* Reading source positions from SOURCE_FILE */ 	
srcpos=sources(&nsrc);
nsrc_glob=nsrc;
signals_stf=matrix(1,nsrc,1,(TIME/DT));
/* set residuals of traces defined in vector tracekill for each shot to zero */
TRKILL=0;
DTRKILL=3;

tracekill =  imatrix(1,nsrc,1,1);
h=1;
for(i=1;i<=nsrc;i++){
  tracekill[i][1]=h;
  h=h+8;
}

snapseis=1;
snapseis1=5;
SHOTINC=1;
RECINC=1;

/* Define used minimization norm */
/* LNORM==1 L1 Norm*/
/* LNORM==2 L2 Norm*/
/* LNORM==3 Cauchy*/
/* LNORM==4 SECH*/
LNORM=2;

if (INVMAT==4){
dsignals=fmatrix(1,nsrc,1,NT);}

QUELLART_OLD = QUELLART;

for(iter=1;iter<=ITERMAX;iter++){  /* fullwaveform iteration loop */	

if (MYID==0)
   {
   time2=MPI_Wtime();
   fprintf(FP,"\n\n\n ------------------------------------------------------------------\n");
   fprintf(FP,"\n\n\n                   TDFWI ITERATION %d \t of %d \n",iter,ITERMAX);
   fprintf(FP,"\n\n\n ------------------------------------------------------------------\n");
   }

/* For the calculation of the material parameters beteween gridpoints
   the have to be averaged. For this, values lying at 0 and NX+1,
   for example, are required on the local grid. These are now copied from the
   neighbouring grids */		
matcopy_elastic(prho, ppi, pu);
MPI_Barrier(MPI_COMM_WORLD);

av_mue(pu,puipjp,prho);
av_rho(prho,prip,prjp);

if(iter==1){
    for (i=1;i<=NX;i=i+IDX){ 
	for (j=1;j<=NY;j=j+IDY){
	
	if(INVMAT1==1){
	
	  Vp0[j][i] = ppi[j][i];
	  Vs0[j][i] = pu[j][i];
	  Rho0[j][i] = prho[j][i];}
	  
                 
		 
	if(INVMAT1==2){
        
	  Vp0[j][i] = sqrt((ppi[j][i]+2.0*pu[j][i])*prho[j][i]);
	  Vs0[j][i] = sqrt((pu[j][i])*prho[j][i]);
	  Rho0[j][i] = prho[j][i];
	
	}
	 
	if(INVMAT1==3){
        
	  Vp0[j][i] = ppi[j][i];
	  Vs0[j][i] = pu[j][i];
	  Rho0[j][i] = prho[j][i];
	
	}  
	
    }
    }

/* ----------------------------- */
/* calculate Covariance matrices */
/* ----------------------------- */

	 Lcount = 1;
	 Vp_avg = 0.0;
	 Vs_avg = 0.0;
	 rho_avg = 0.0;
	 
        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
	  
		 /* calculate average Vp, Vs */
                 Vp_avg+=ppi[j][i];
		 Vs_avg+=pu[j][i];
		 
		 /* calculate average rho */
		 rho_avg+=prho[j][i];
		 Lcount++;
           }
        }
		
        /* calculate average Vp, Vs and rho of all CPUs*/
        Lcountsum = 0.0;
        MPI_Allreduce(&Lcount,&Lcountsum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        Lcount=Lcountsum;
	
        Vp_sum = 0.0;
        MPI_Allreduce(&Vp_avg,&Vp_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        Vp_avg=Vp_sum;
	
	Vs_sum = 0.0;
        MPI_Allreduce(&Vs_avg,&Vs_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        Vs_avg=Vs_sum;
	
	rho_sum = 0.0;
        MPI_Allreduce(&rho_avg,&rho_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
        rho_avg=rho_sum;
	
	Vp_avg /=Lcount; 
	Vs_avg /=Lcount; 
	rho_avg /=Lcount;
	
	
        printf("MYID = %d \t Vp_avg = %e \t Vs_avg = %e \t rho_avg = %e \n ",MYID,Vp_avg,Vs_avg,rho_avg);	
	
	C_vp = Vp_avg*Vp_avg;
	C_vs = Vs_avg*Vs_avg;
	C_rho = rho_avg*rho_avg;


}

/* Open Log File for L2 norm */

if(MYID==0){
if(iter==1){
FPL2=fopen("L2_LOG.dat","w");}

if(iter>1){
FPL2=fopen("L2_LOG.dat","a");}
}

/* initialization of L2 calculation */
L2=0.0;
Lcount=0;
energy=0.0;
L2_all_shots=0.0;
energy_all_shots=0.0;

EPSILON=0.0;  /* test step length */
exchange_par();

/* initialize waveconv matrix*/
if((INVMAT==1)||(INVMAT==0)){
for (i=1;i<=NX;i=i+IDX){ 
	for (j=1;j<=NY;j=j+IDY){
	    waveconv[j][i]=0.0;    
        }
     }   
}

if((INVMAT==2)||(INVMAT==0)){
for (i=1;i<=NX;i=i+IDX){ 
	for (j=1;j<=NY;j=j+IDY){
	    waveconv_rho[j][i]=0.0;    
        }
     }   
}

if((INVMAT==3)||(INVMAT==0)){
for (i=1;i<=NX;i=i+IDX){ 
	for (j=1;j<=NY;j=j+IDY){
	    waveconv_u[j][i]=0.0;    
        }
     }   
}

if(HESSIAN){
  for (i=1;i<=NX;i=i+IDX){
      for (j=1;j<=NY;j=j+IDY){
         hessian[j][i]=0.0;
	 hessian_u[j][i]=0.0;
	 hessian_rho[j][i]=0.0;
      }
  }
}

itestshot=TESTSHOT_START;
swstestshot=0;

if(INVTYPE==2){ 
if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;
  
        for (ishot=1;ishot<=nshots;ishot+=SHOTINC){
        /*for (ishot=1;ishot<=1;ishot+=1){*/
 
        fprintf(FP,"\n==================================================================================\n");
        fprintf(FP,"\n MYID=%d *****  Starting simulation (forward model) for shot %d of %d  ********** \n",MYID,ishot,nshots);
	fprintf(FP,"\n==================================================================================\n\n");
		
                for (nt=1;nt<=6;nt++) srcpos1[nt][1]=srcpos[nt][ishot]; 
		
                if (RUN_MULTIPLE_SHOTS){

                                /* find this single source positions on subdomains */
                                if (nsrc_loc>0) free_matrix(srcpos_loc,1,6,1,1);
                                srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1);
		        }


                else{
                                /* Distribute multiple source positions on subdomains */
                               srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);
                }

QUELLART = QUELLART_OLD;
MPI_Barrier(MPI_COMM_WORLD);
/* calculate wavelet for each source point */
signals=NULL;
signals=wavelet(srcpos_loc,nsrc_loc);

/*char  source_signal_file[STRING_SIZE];
sprintf(source_signal_file,"source_signal.%d.su.shot%d.it%d",MYID,ishot,iter);
fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
output_source_signal(fopen(source_signal_file,"w"),signals,NT,3);

/* output source signal e.g. for cross-correlation of comparison with analytical solutions */

if(nsrc_loc>0){
        if(QUELLART==6){/*FFT_filt(signals,1.0,1,ns,1);}*/
			timedomain_filt(signals,FC_HESSIAN,ORDER_HESSIAN,nsrc_loc,ns,1);}
	
	if(HESSIAN==1){	/*calculation of Hessian matrix; in this case also the forward modelled wavefields
			must be filtered with the Butterworth filter applied to the delta impulse during
			the backpropagation*/
			timedomain_filt(signals,FC_HESSIAN,ORDER_HESSIAN,nsrc_loc,ns,1);}
			
            char  source_signal_file[STRING_SIZE];
	    sprintf(source_signal_file,"%s_source_signal.%d.su.shot%d", MFILE, MYID,ishot);
	    fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
	    output_source_signal(fopen(source_signal_file,"w"),signals,NT,1);
}                                
                                
MPI_Barrier(MPI_COMM_WORLD);

if((iter>1)&&(nsrc_loc>0)&&(INVMAT==4)){

   /* find maximum of dsignals */
   sig_max = 0.0;     
       
       for(iq=1;iq<=NT;iq++){
       
         if(fabs(dsignals[ishot][iq])>sig_max){
            sig_max = fabs(dsignals[ishot][iq]);
	 }
       }   
  
       /* update wavelet */
       for(iq=1;iq<=NT;iq++){
        signals[1][iq] -= 0.05 * (dsignals[ishot][iq]/sig_max);
	
       } 
          
}

/*if(FREQFILT==1){
FFT_filt(signals,freqshift,1,ns,1);}*/

/*if ((nsrc_loc>0)&&(INVMAT==4)){
     sprintf(dwavelet,"wavelet/dwavelet_shot_%d_iter_%d.dat",ishot,iter);

     FP7=fopen(dwavelet,"w");
       for(iq=1;iq<=NT;iq++){
        fprintf(FP7,"%e \n",signals[1][iq]);}
	
     fclose(FP7);	
}*/

/*if (nsrc_loc>0){
     sprintf(dwavelet,"wavelet/dwavelet_shot_%d_iter_%d.dat",ishot,iter);

     FP7=fopen(dwavelet,"w");
       for(iq=1;iq<=NT;iq++){
        fprintf(FP7,"%e \n",signals[1][iq]);}
	
     fclose(FP7);	
}*/

		    
/* initialize wavefield with zero */
zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);	


/*initialize gradient matrices for each shot with zeros*/
if((INVMAT==1)||(INVMAT==0)){
	for(i=1;i<=NX;i=i+IDX){
		for(j=1;j<=NY;j=j+IDY){
			waveconv_shot[j][i]=0.0;
		}
	}
}


if((INVMAT==2)||(INVMAT==0)){
	for(i=1;i<=NX;i=i+IDX){
		for(j=1;j<=NY;j=j+IDY){
			waveconv_rho_shot[j][i]=0.0;
		}
	}
}


if((INVMAT==3)||(INVMAT==0)){
	for(i=1;i<=NX;i=i+IDX){
		for(j=1;j<=NY;j=j+IDY){
			waveconv_u_shot[j][i]=0.0;
		}
	}
}

if(HESSIAN){
h=1;
 for (nt=1;nt<=NT;nt=nt+DTINV){
   for (i=1;i<=NX;i=i+IDX){ 
     	for (j=1;j<=NY;j=j+IDY){
	      
	      jac_lam_x[h]=0.0;
	      jac_lam_y[h]=0.0;
	      jac_u[h]=0.0;
	      jac_rho[h]=0.0;

	h++;

        }
   }
 }
}
                                                         
     
/*----------------------  loop over timesteps (forward model) ------------------*/

lsnap=iround(TSNAP1/DT);  
lsamp=NDT;
nsnap=0;

hin=1;
hin1=1;

imat=1;
imat1=1;
imat2=1;
hi=1;

for (nt=1;nt<=NT;nt++){     
                
	/* Check if simulation is still stable */
        /*if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");*/
	if (isnan(pvy[NY/2][NX/2])) {
	   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
	   err(" Simulation is unstable !");}

		
   infoout = !(nt%1000);

   if (MYID==0){
      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
      time3=MPI_Wtime();
   }

      /* update of particle velocities */
      /*update_v_hc(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,2);*/ 
      update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x,
                   b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);
                         
	if (MYID==0){
		time4=MPI_Wtime();
		time_av_v_update+=(time4-time3);
		if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
	}
                                                   
	/* exchange of particle velocities between PEs */
	exchange_v(pvx,pvy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec);
                                                               
	if (MYID==0){
	  time5=MPI_Wtime();
	  time_av_v_exchange+=(time5-time4);
	  if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
	}                                                                                         	

    /*update_s_elastic_hc(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout);*/
   update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout,
                        K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);  


    /* explosive source */
   if ((!CHECKPTREAD)&&(QUELLTYP==1)) 	
   psource(nt,psxx,psyy,srcpos_loc,signals,nsrc_loc);

   if ((FREE_SURF) && (POS[2]==0))
   surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy, ppi, pu, prho, hc, K_x, a_x, b_x, psi_vxxs);


   if (MYID==0){
      time6=MPI_Wtime();
	  time_av_s_update+=(time6-time5);
      if (infoout)  fprintf(FP," stress exchange between PEs ...");
      }


   /* stress exchange between PEs */
    exchange_s(psxx,psyy,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec);


 
   if (MYID==0){
      time7=MPI_Wtime();
 	  time_av_s_exchange+=(time7-time6);
     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
      }  

	/* store amplitudes at receivers in section-arrays */
	if ((SEISMO) && (nt==lsamp) && (nt<=NT)){
		seismo_ssg(lsamp, ntr, recpos_loc, sectionvx, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppi, pu, hc);
		lsamp+=NDT;
	}


if(nt==hin1){

    if((INVMAT==2)||(INVMAT==0)){
        for (i=1;i<=NX;i=i+IDXI){
	    for (j=1;j<=NY;j=j+IDYI){
		 forward_prop_rho_x[imat1]=pvxp1[j][i];
		 forward_prop_rho_y[imat1]=pvyp1[j][i];
                 imat1++;                                   
		 }
    }}   

    /* save snapshots from forward model */  
    if((INVMAT==1)||(INVMAT==0)){
    for (i=1;i<=NX;i=i+IDXI){ 
	for (j=1;j<=NY;j=j+IDYI){
	    forward_prop_x[imat]=psxx[j][i];
	    forward_prop_y[imat]=psyy[j][i];
	    
	    imat++;
        }
     }}
    
    if((INVMAT==3)||(INVMAT==0)){
    for (i=1;i<=NX;i=i+IDXI){ 
	for (j=1;j<=NY;j=j+IDYI){
	    forward_prop_u[imat2]=psxy[j][i];
	    imat2++;
        }
     } 
     
    hin++;
    hin1=hin1+DTINV;
                                                             
}

}

   /* WRITE SNAPSHOTS TO DISK */
   if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){

      snap(FP,nt,++nsnap,pvx,pvy,psxx,psyy,pu,ppi,hc);

      lsnap=lsnap+iround(TSNAPINC/DT);
      }

   
      
   if (MYID==0){
      time8=MPI_Wtime();
	  time_av_timestep+=(time8-time3);
      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
      }   		


   }/*--------------------  End  of loop over timesteps (forward model) ----------*/

	
	
if (SEISMO){

	catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
	catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
	
	if (MYID==0){
	  s=0;
	saveseis_glob(FP,fulldata_vx,fulldata_vy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,s);
	}
}


/*----------- Start of inversion for source time function -------------*/
if((INV_STF==1)&&(iter==iter)){
  
  N_STF=N_STF+N_STF;

catseis(sectionvy, fulldata, recswitch, ntr_glob, MPI_COMM_WORLD);

if (MYID==0){

inseis(fprec,ishot,sectionvy_obs,ntr_glob,ns,2);

stf(FP,fulldata,sectionvy_obs,sectionvy_conv,source_time_function,recpos,recpos_loc,ntr_glob,ntr,srcpos,ishot,ns,iter,nsrc_glob,signals_stf,signals);

	}

MPI_Barrier(MPI_COMM_WORLD);

/*----------------------  loop over timesteps (second forward model) due to source time function ------------------*/

fprintf(FP,"\n=================================================================================\n");
fprintf(FP,"\n MYID=%d *****  Starting second simulation (forward model) due to inversion *****\n",MYID);
fprintf(FP,"\n         *****  for source time function of shot %d of %d                  *****\n",ishot,nshots);
fprintf(FP,"\n=================================================================================\n\n");

/* calculate wavelet for each source point */
signals=NULL;
signals=wavelet_stf(nsrc_loc,ishot,signals_stf);

/*window_cos(signals,NT,nsrc_loc,0.0,0.001,0.19,TIME);*/

/*char  source_signal_file[STRING_SIZE];
sprintf(source_signal_file,"source_signal_stf.%d.su.shot%d.it%d",MYID,ishot,iter);
fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
output_source_signal(fopen(source_signal_file,"w"),signals,NT,3);*/
/*outseis_glob(FP,fopen(source_signal_file,"w"),1,signals,recpos,recpos_loc,ntr,srcpos,0,ns,SEIS_FORMAT,ishot);*/


lsnap=iround(TSNAP1/DT);  
lsamp=NDT;
nsnap=0;

hin=1;
hin1=1;

imat=1;
imat1=1;
imat2=1;
hi=1;

for (nt=1;nt<=NT;nt++){     
                
	/* Check if simulation is still stable */
        /*if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");*/
	if (isnan(pvy[NY/2][NX/2])) {
	   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
	   err(" Simulation is unstable !");}

		
   infoout = !(nt%1000);

   if (MYID==0){
      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
      time3=MPI_Wtime();
   }

      /* update of particle velocities */
      /*update_v_hc(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,2);*/ 
      update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x,
                   b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);
                         
	if (MYID==0){
		time4=MPI_Wtime();
		time_av_v_update+=(time4-time3);
		if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
	}
                                                   
	/* exchange of particle velocities between PEs */
	exchange_v(pvx,pvy, bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, req_send, req_rec);
                                                               
	if (MYID==0){
	  time5=MPI_Wtime();
	  time_av_v_exchange+=(time5-time4);
	  if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
	}                                                                                         	

    /*update_s_elastic_hc(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout);*/
   update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout,
                        K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);  


    /* explosive source */
   if ((!CHECKPTREAD)&&(QUELLTYP==1)) 	
   psource(nt,psxx,psyy,srcpos_loc,signals,nsrc_loc);

   if ((FREE_SURF) && (POS[2]==0))
   surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy, ppi, pu, prho, hc, K_x, a_x, b_x, psi_vxxs);


   if (MYID==0){
      time6=MPI_Wtime();
	  time_av_s_update+=(time6-time5);
      if (infoout)  fprintf(FP," stress exchange between PEs ...");
      }


   /* stress exchange between PEs */
    exchange_s(psxx,psyy,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec);


 
   if (MYID==0){
      time7=MPI_Wtime();
 	  time_av_s_exchange+=(time7-time6);
     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
      }  

	/* store amplitudes at receivers in section-arrays */
	if ((SEISMO) && (nt==lsamp) && (nt<=NT)){
		seismo_ssg(lsamp, ntr, recpos_loc, sectionvx, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppi, pu, hc);
		lsamp+=NDT;
	}


if(nt==hin1){

    if((INVMAT==2)||(INVMAT==0)){
        for (i=1;i<=NX;i=i+IDXI){
	    for (j=1;j<=NY;j=j+IDYI){
		 forward_prop_rho_x[imat1]=pvxp1[j][i];
		 forward_prop_rho_y[imat1]=pvyp1[j][i];
                 imat1++;                                   
		 }
    }}   

    /* save snapshots from forward model */  
    if((INVMAT==1)||(INVMAT==0)){
    for (i=1;i<=NX;i=i+IDXI){ 
	for (j=1;j<=NY;j=j+IDYI){
	    forward_prop_x[imat]=psxx[j][i];
	    forward_prop_y[imat]=psyy[j][i];
	    
	    imat++;
        }
     }}
    
    if((INVMAT==3)||(INVMAT==0)){
    for (i=1;i<=NX;i=i+IDXI){ 
	for (j=1;j<=NY;j=j+IDYI){
	    forward_prop_u[imat2]=psxy[j][i];
	    imat2++;
        }
     } 
     
    hin++;
    hin1=hin1+DTINV;
                                                             
}

}

   /* WRITE SNAPSHOTS TO DISK */
   if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){

      snap(FP,nt,++nsnap,pvx,pvy,psxx,psyy,pu,ppi,hc);

      lsnap=lsnap+iround(TSNAPINC/DT);
      }

   
      
   if (MYID==0){
      time8=MPI_Wtime();
	  time_av_timestep+=(time8-time3);
      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
      }   		


   }/*--------------------  End  of loop over timesteps (forward model) ----------*/

if (SEISMO){

	catseis(sectionvx, fulldata_vx, recswitch, ntr_glob, MPI_COMM_WORLD);
	catseis(sectionvy, fulldata_vy, recswitch, ntr_glob, MPI_COMM_WORLD);
	
	if (MYID==0){
	  s=1;
	saveseis_glob(FP,fulldata_vx,fulldata_vy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr_glob,srcpos,ishot,ns,iter,s);
	}
}
s=0;
}/*----------- End of inversion for source time function -------------*/


/*if ((ntr > 0) && (SEISMO) && (iter==snapseis) && (ishot==50)){
	saveseis(FP,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos1,ishot,ns,iter);
	snapseis = snapseis + snapseis1;
}*/


if(INVMAT==0){

/* define frequency shift */

freqshift=4.0;

/*if(iter>= 100){
freqshift=2.0;}*/

/*if(iter>= 200){
freqshift=4.0;}*/

if(iter>= 200){
FREQFILT=0;}

if (MYID==0){
printf("Calculate residuals  \n");
printf("-------------------  \n");
}

if ((ntr > 0)&&(HESSIAN==0)){

/* calculate L2-Norm and energy ? */
if((ishot==itestshot)&&(ishot<=TESTSHOT_END)){swstestshot=1;}

/* read seismic data from SU file vx */
/* --------------------------------- */
if((QUELLTYPB==1)||(QUELLTYPB==3)){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,1);

if((TRKILL==1)&&(ishot>1)){
  for(i=tracekill[ishot][1]-DTRKILL;i<=tracekill[ishot][1]+DTRKILL;i++){
     sectionread[i][1]=20000;
  }
}

/* assign input data to each PE */
h=1;
for(i=REC1;i<=REC2;i++){
   for(j=1;j<=ns;j++){
           sectionvxdata[h][j]=sectionread[i][j];
   }
   h++;
}
  L2=calc_res(sectionvxdata,sectionvx,sectionvxdiff,sectionvxdiffold,ntr,ns,LNORM,L2,0,1,swstestshot);
  if(swstestshot==1){energy=calc_energy(sectionvxdata,ntr,ns,energy);}
    /*fprintf(FP,"Energy vxdata for PE %d:   %f\n\n", MYID,energy);*/
  L2_all_shots=calc_misfit(sectionvxdata,sectionvx,ntr,ns,LNORM,L2_all_shots,0,1,1);
  energy_all_shots=calc_energy(sectionvxdata,ntr,ns,energy_all_shots);
  /*fprintf(FP,"Energy vxdata for PE %d:   %f\n\n", MYID,energy);*/


/* apply frequency filter on the residuals (x-component) */
if(FREQFILT==1){
FFT_filt(sectionvxdiff,freqshift,ntr,ns,1);}
} /* end QUELLTYPB */

/* read seismic data from SU file vy */
/* --------------------------------- */
if((QUELLTYPB==1)||(QUELLTYPB==2)){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,2);

if((TRKILL==1)&&(ishot>1)){
  for(i=tracekill[ishot][1]-DTRKILL;i<=tracekill[ishot][1]+DTRKILL;i++){
     sectionread[i][1]=20000;
  }
}

/* assign input data to each PE */
h=1;
for(i=REC1;i<=REC2;i++){
   for(j=1;j<=ns;j++){
           sectionvydata[h][j]=sectionread[i][j];
   }
   h++;
}


L2=calc_res(sectionvydata,sectionvy,sectionvydiff,sectionvydiffold,ntr,ns,LNORM,L2,0,1,swstestshot);
if(swstestshot==1){energy=calc_energy(sectionvydata,ntr,ns,energy);}
L2_all_shots=calc_misfit(sectionvydata,sectionvy,ntr,ns,LNORM,L2_all_shots,0,1,1);
energy_all_shots=calc_energy(sectionvydata,ntr,ns,energy_all_shots);
/*fprintf(FP,"Energy vydata for PE %d:   %f\n\n", MYID,energy);	*/	   	    


/* apply frequency filter on the residuals (y-component) */
if(FREQFILT==1){
FFT_filt(sectionvydiff,freqshift,ntr,ns,1);}

} /* end QUELLTYPB */

if((ishot==itestshot)&&(ishot<=TESTSHOT_END)){
       swstestshot=0;
       itestshot+=TESTSHOT_INCR;
}

if(TIMEWIN==1){
time_window(sectionvydiff,picked_times,iter,ntr,ns);
time_window(sectionvxdiff,picked_times,iter,ntr,ns);}

if(TIMEWIN==2){
stalta(sectionvy,ntr,ns,picked_times);
time_window(sectionvydiff,picked_times,iter,ntr,ns);

stalta(sectionvx,ntr,ns,picked_times);
time_window(sectionvxdiff,picked_times,iter,ntr,ns);}


if (SEISMO){
saveseis(FP,sectionvxdiff,sectionvydiff,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos1,ishot,ns,-1);
}

} /* end HESSIAN != 1 */

/*if ((ntr > 0)&&(TAPER)){*/
  /* taper seismic section at the boundaries to avoid artifical diffractions */
/*  taper(sectionvxdiff,ntr,ns);
}*/

	    		    
/* --------------------------------------------------------------------------------------------------------------------- */

/*srcpos_loc_back=splitsrc_back(recpos,&nsrc_loc,ntr_glob);*/

if (HESSIAN) nshots1=ntr_glob; else nshots1=1;
  
for (irec=1;irec<=nshots1;irec+=RECINC){ /* loop over shots at receiver positions */

/*for (irec=20;irec<=20;irec+=1){*/

hin=1;
hin1=1;

    if(MYID==0){
    printf("\n==================================================================================\n");
    printf("\n MYID=%d *****  Starting simulation (backward model) for shot %d of %d  ********** \n",MYID,irec,nshots1);
    printf("\n==================================================================================\n\n");}
    
    	
                if (HESSIAN){
                
                                /* find this single source positions on subdomains */
				srcpos_loc_back = matrix(1,6,1,1);
				ntr1=0;
				
				if((irec>=REC1)&&(irec<=REC2)){
				  srcpos_loc_back[1][1] = (recpos[1][irec]);
				  srcpos_loc_back[2][1] = (recpos[2][irec]);
				  srcpos_loc_back[4][1] = TSHIFT_back;
				  srcpos_loc_back[6][1] = 1.0;
				  ntr1=1;  
				}
				
				QUELLART=6;
	                        
	                        /*printf("MYID = %d \t REC1 = %d \t REC2 = %d \t ntr1 = %d x = %e \t y = %e \n",MYID,REC1,REC2,ntr1,srcpos_loc_back[1][1],srcpos_loc_back[2][1]);*/
	                        
	                        if(ntr1>0){
	                          /* calculate wavelet for each source point */
				  
	                          signals_rec=wavelet(srcpos_loc_back,ntr1);
	                        
	                          /* output source signal e.g. for cross-correlation of comparison with analytical solutions */
	                          if((QUELLART==6)&&(ntr1>0)){timedomain_filt(signals_rec,FC_HESSIAN,ORDER_HESSIAN,1,ns,1);}
	                          
			          /*char  source_signal_file[STRING_SIZE];
				  sprintf(source_signal_file,"%s_source_signal_back.%d.su", MFILE, MYID);
				  fprintf(stdout,"\n PE %d outputs source time function in SU format to %s \n ", MYID, source_signal_file);
				  output_source_signal(fopen(source_signal_file,"w"),signals_rec,NT,1);*/
	                                                                          
	                        }
	                        
	                        MPI_Barrier(MPI_COMM_WORLD);
	                                                
	                  if(ntr1 > 0){
	                  h=1;
	                   for(i=REC1;i<=REC2;i++){
				         for(j=1;j<=ns;j++){  
					     sectionvxdiff[h][j]=signals_rec[1][j];
				   	     sectionvydiff[h][j]=signals_rec[1][j];
				         }  
				      h++;
				      }   
	
	                                                                                
				  if (SEISMO){
				     saveseis(FP,sectionvxdiff,sectionvydiff,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos1,ishot,ns,-1);
				  }
			        }	   
				                                                                                                   
		}

		else{
                                /* Distribute multiple source positions on subdomains */
                                /* define source positions at the receivers */
                                srcpos_loc_back = matrix(1,6,1,ntr);
                                for (i=1;i<=ntr;i++){
                                  srcpos_loc_back[1][i] = (recpos_loc[1][i]);
                                  srcpos_loc_back[2][i] = (recpos_loc[2][i]);
                                }
                                ntr1=ntr;
                                    
                }

/* initialize wavefield with zero */
zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);	

/*----------------------  loop over timesteps (backpropagation) ------------------*/

lsnap=iround(TSNAP1/DT);  
lsamp=NDT;
nsnap=0;
if(HESSIAN==1){imat=1;}

/*invtimer=NT/DTINV;*/
for (nt=1;nt<=NT;nt++){     

	
	/* Check if simulation is still stable */
        /*if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");*/
	if (isnan(pvy[NY/2][NX/2])) {
	   fprintf(FP,"\n Time step: %d; pvy: %f \n",nt,pvy[NY/2][NX/2]);
	   err(" Simulation is unstable !");}

		
   infoout = !(nt%5000);

   if (MYID==0){
      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
      time3=MPI_Wtime();
   }

      /*if(MYID==0){  
      printf("Calculate zero lag x-correlation of forward and backpropagated wavefields \n");
      printf("------------------------------------------------------------------------- \n");
      }*/

   /* update of particle velocities */
   /*update_v_hc(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp,
   	       srcpos_loc_back,sectionvxdiff,sectionvydiff,ntr,absorb_coeff, hc, infoout,1);*/

   update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc_back, sectionvxdiff,sectionvydiff,ntr1,absorb_coeff,hc,infoout,1, K_x, a_x,
                   b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x); 

   if (MYID==0){
      time4=MPI_Wtime();
      time_av_v_update+=(time4-time3);
     if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
      }

   /* exchange of particle velocities between PEs */
   exchange_v(pvx,pvy, bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top, req_send, req_rec);

   if (MYID==0){
      time5=MPI_Wtime();
	  time_av_v_exchange+=(time5-time4);
      if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
   }

   /*update_s_elastic_hc(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp,absorb_coeff, prho, hc, infoout);*/

   update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppi, pu, puipjp, absorb_coeff, prho, hc, infoout,
                        K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);
			

     if ((FREE_SURF) && (POS[2]==0))
         surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy, ppi, pu, prho, hc, K_x, a_x, b_x, psi_vxxs);


   if (MYID==0){
      time6=MPI_Wtime();
	  time_av_s_update+=(time6-time5);
      if (infoout)  fprintf(FP," stress exchange between PEs ...");
   }


   /* stress exchange between PEs */
    exchange_s(psxx,psyy,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec);


 
   if (MYID==0){
      time7=MPI_Wtime();
 	time_av_s_exchange+=(time7-time6);
     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
   }

   /* calculate change of the source wavelet */
   if ((nsrc_loc>0)&&(INVMAT==4)){
      for (lq=1;lq<=nsrc_loc;lq++) {
      
                iq=(int)srcpos_loc[1][lq];
                jq=(int)srcpos_loc[2][lq];

        dsignals[ishot][invtimer] = (-psxx[jq][lq]-psyy[jq][lq])/2.0;}
	
      }


	/* store amplitudes at receivers in section-arrays */
	/*if ((SEISMO) && (nt==lsamp) && (nt<=NT)){
		seismo_ssg(lsamp, ntr, recpos_loc, sectionp, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppi, pu, hc);
		lsamp+=NDT;
	} */ 
      
   if(nt==hin1){
    
        if(HESSIAN==0){imat=((nxnyi*((NT/DTINV))) - hin*nxnyi)+1;}
              
        if((HESSIAN==0)&&(INVMAT==0)){
	    for (i=1;i<=NX;i=i+IDXI){   
	        for (j=1;j<=NY;j=j+IDYI){ 
                                             
		   waveconv_rho_shot[j][i]+=(pvxp1[j][i]*forward_prop_rho_x[imat])+(pvyp1[j][i]*forward_prop_rho_y[imat]);
		   waveconv_shot[j][i]+= (forward_prop_x[imat]+forward_prop_y[imat])*(psxx[j][i]+psyy[j][i]);  
		   
		   muss = prho[j][i] * pu[j][i] * pu[j][i];
	          lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;
		                
		   if(pu[j][i]>0.0){
		   waveconv_u_shot[j][i]+= ((1.0/(muss*muss))*(forward_prop_u[imat] * psxy[j][i])) 
                                  + ((1.0/4.0) * ((forward_prop_x[imat] + forward_prop_y[imat]) * (psxx[j][i] + psyy[j][i])) / ((lamss+muss)*(lamss+muss)))  
                                  + ((1.0/4.0) * ((forward_prop_x[imat] - forward_prop_y[imat]) * (psxx[j][i] - psyy[j][i])) / (muss*muss));}
                      		                                                                                                             
		   imat++;
		   }
	    }}

	 if((HESSIAN==1)&&(INVMAT==0)){
	    for (i=1;i<=NX;i=i+IDXI){   
	        for (j=1;j<=NY;j=j+IDYI){ 
                                             
		   /*jac_rho[imat]+=(pvxp1[j][i]*forward_prop_rho_x[imat])+(pvyp1[j][i]*forward_prop_rho_y[imat]);*/
		   jac_lam_x[imat] = psxx[j][i];
		   jac_lam_y[imat] = psyy[j][i];
		       jac_u[imat] = psxy[j][i];
		                         		                                                                                                             
		   imat++;
		   }
	    }}

                                                                                                                               
    hin1=hin1+DTINV;
    hin++;
    }
     
     /*fclose(FP2);*/
     
   /* WRITE SNAPSHOTS TO DISK */
   if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){

      snap(FP,nt,++nsnap,pvx,pvy,psxx,psyy,pu,ppi,hc);

      lsnap=lsnap+iround(TSNAPINC/DT);
      }
                                                                                                                                               
   
   if (MYID==0){
      time8=MPI_Wtime();
	time_av_timestep+=(time8-time3);
      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
      }   		


   }/*--------------------  End  of loop over timesteps (backpropagation)----------*/

if(HESSIAN){

  /* calculate Hessian for lambda from autocorrelation of the jacobian */
  /* ----------------------------------------------------------------- */
  imat=1;
  
  imat1=1;
  for (i=1;i<=NX;i=i+IDXI){   
	    for (j=1;j<=NY;j=j+IDYI){ 
             
			 h=1;
			 for (nt=1;nt<=NT;nt=nt+DTINV){
                             imat = imat1 + nxnyi*(nt-1);
		     
                             temp_TS[h]  =  jac_lam_x[imat] + jac_lam_y[imat];
			     temp_TS1[h] =  forward_prop_x[imat] + forward_prop_y[imat];

			     temp_TS2[h] =  jac_lam_x[imat] - jac_lam_y[imat];
			     temp_TS3[h] =  forward_prop_x[imat] - forward_prop_y[imat];

			     temp_TS4[h] =  jac_u[imat];
			     temp_TS5[h] =  forward_prop_u[imat];

			 h++;
	                 }

			 /* convolve time series in frequency domain */
			 conv_FD(temp_TS,temp_TS1,temp_conv,(NT/DTINV));
			 conv_FD(temp_TS2,temp_TS3,temp_conv1,(NT/DTINV));
			 conv_FD(temp_TS4,temp_TS5,temp_conv2,(NT/DTINV));

			 h=1;
			 for (nt=1;nt<=NT;nt=nt+DTINV){
				             
			     /* Hessian for lambda */
                             hessian[j][i] += temp_conv[h]*temp_conv[h];
							 
			     /* Hessian for Mu */
			     muss = prho[j][i] * pu[j][i] * pu[j][i];
			     lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 * muss;
		             
		             mulamratio = (muss * muss)/((lamss+muss)*(lamss+muss));
		                
			     if(pu[j][i]>0.0){
		                  
		                  /*temp_hess = ((1.0/(muss*muss))*(temp_conv2[h])) 
                                  + ((1.0/4.0) * (temp_conv[h]) / ((lamss+muss)*(lamss+muss)))  
                                  + ((1.0/4.0) * (temp_conv1[h]) / (muss*muss));*/
                                  
                                  temp_hess = temp_conv2[h] + (((mulamratio*temp_conv[h]) + temp_conv1[h])/4.0);
                                  
                                  hessian_u[j][i] += temp_hess*temp_hess;
                                  
			     }

				 h++;
			 }

			 imat1++;
         }
  }

  /* calculate Hessian for rho from autocorrelation of the jacobian */
  /* ----------------------------------------------------------------- */
  imat=1;

  for (nt=1;nt<=NT;nt=nt+DTINV){ 
     for (i=1;i<=NX;i=i+IDXI){   
	    for (j=1;j<=NY;j=j+IDYI){ 
                                               
		     hessian_rho[j][i] += jac_rho[imat]*jac_rho[imat];
                        		                                                                                                             
		     imat++;
	     }
     }
  }

sprintf(jac,"%s_hessian_rho_shot%i.%i%i",JACOBIAN,ishot,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
        fwrite(&hessian_rho[j][i],sizeof(float),1,FP4);
   }
}

fclose(FP4);

MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */ 
sprintf(jac,"%s_hessian_rho_shot%i",JACOBIAN,ishot); 
if (MYID==0) mergemod(jac,3);


} /* end HESSIAN */


} /* end of loop over shots */
	
/*if ((ntr > 0) && (SEISMO)){
	saveseis(FP,sectionpdiff,sectionvy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,srcpos1,ishot,ns,0);
}*/

/*if((INVMAT==1)||(INVMAT==0)){
if(MYID==0){
  printf("Output of zero lag x-correlation pi\n");
  printf("----------------------------------- \n");	 
}


/*if((INVMAT==2)||(INVMAT==0)){
if(MYID==0){
  printf("Output of zero lag x-correlation rho \n");
  printf("------------------------------------ \n");	 
}

sprintf(jac,"%s_rho.shot%i.%i%i",JACOBIAN,ishot,POS[1],POS[2]);*/
/*printf("%s\n",jac);*/
/*FP3=fopen(jac,"wb");*/

/* save Jacobian */

/*     for (i=1;i<=NX;i=i+IDX){ 
	for (j=1;j<=NY;j=j+IDY){
	    dump=waveconv_rho[j][i];
	    fwrite(&dump, sizeof(float), 1, FP3);
        }
     }
     
fclose(FP3);
}

if((INVMAT==3)||(INVMAT==0)){
if(MYID==0){
  printf("Output of zero lag x-correlation u \n");
  printf("------------------------------------ \n");	 
}

sprintf(jac,"%s_u.shot%i.%i%i",JACOBIAN,ishot,POS[1],POS[2]);*/
/*printf("%s\n",jac);*/
/*FP3=fopen(jac,"wb");*/

/* save Jacobian */

/*     for (i=1;i<=NX;i=i+IDX){ 
	for (j=1;j<=NY;j=j+IDY){
	    dump=waveconv_u[j][i];
	    fwrite(&dump, sizeof(float), 1, FP3);
        }
     }
     
fclose(FP3);
}*/


if (SWS_TAPER_CIRCULAR_PER_SHOT){    /* applying a circular taper at the source position to the gradient of each shot */
	
	/* saving the gradients for testing */
	/*sprintf(jac,"%s_shot%i.old.%i%i",JACOBIAN,ishot,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	for (i=1;i<=NX;i=i+IDX){
		for (j=1;j<=NY;j=j+IDY){
			fwrite(&waveconv_shot[j][i],sizeof(float),1,FP3);
		}
	}
	fclose(FP3);
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(jac,"%s_shot%i.old",JACOBIAN,ishot);
	if (MYID==0) mergemod(jac,3);


	sprintf(jac,"%s_u_shot%i.old.%i%i",JACOBIAN,ishot,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	for (i=1;i<=NX;i=i+IDX){
		for (j=1;j<=NY;j=j+IDY){
			fwrite(&waveconv_u_shot[j][i],sizeof(float),1,FP3);
		}
	}
	fclose(FP3);
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(jac,"%s_u_shot%i.old",JACOBIAN,ishot);
	if (MYID==0) mergemod(jac,3);

	sprintf(jac,"%s_rho_shot%i.old.%i%i",JACOBIAN,ishot,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	for (i=1;i<=NX;i=i+IDX){
		for (j=1;j<=NY;j=j+IDY){
			fwrite(&waveconv_rho_shot[j][i],sizeof(float),1,FP3);
		}
	}
	fclose(FP3);
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(jac,"%s_rho_shot%i.old",JACOBIAN,ishot);
	if (MYID==0) mergemod(jac,3); */
	
	
	
	/* applying the preconditioning */
	taper_grad_shot(waveconv_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot);
	taper_grad_shot(waveconv_rho_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot);
	taper_grad_shot(waveconv_u_shot,taper_coeff,srcpos,nsrc,recpos,ntr_glob,ishot);
	
	
	
	/* again saving the gradients for testing */
	/*sprintf(jac,"%s_precond_shot%i.old.%i%i",JACOBIAN,ishot,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	for (i=1;i<=NX;i=i+IDX){
		for (j=1;j<=NY;j=j+IDY){
			fwrite(&waveconv_shot[j][i],sizeof(float),1,FP3);
		}
	}
	fclose(FP3);
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(jac,"%s_precond_shot%i.old",JACOBIAN,ishot);
	if (MYID==0) mergemod(jac,3);


	sprintf(jac,"%s_u_precond_shot%i.old.%i%i",JACOBIAN,ishot,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	for (i=1;i<=NX;i=i+IDX){
		for (j=1;j<=NY;j=j+IDY){
			fwrite(&waveconv_u_shot[j][i],sizeof(float),1,FP3);
		}
	}
	fclose(FP3);
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(jac,"%s_u_precond_shot%i.old",JACOBIAN,ishot);
	if (MYID==0) mergemod(jac,3);

	sprintf(jac,"%s_rho_precond_shot%i.old.%i%i",JACOBIAN,ishot,POS[1],POS[2]);
	FP3=fopen(jac,"wb");
	for (i=1;i<=NX;i=i+IDX){
		for (j=1;j<=NY;j=j+IDY){
			fwrite(&waveconv_rho_shot[j][i],sizeof(float),1,FP3);
		}
	}
	fclose(FP3);
	MPI_Barrier(MPI_COMM_WORLD);
	sprintf(jac,"%s_rho_precond_shot%i.old",JACOBIAN,ishot);
	if (MYID==0) mergemod(jac,3);*/ 

} /* end of SWS_TAPER_CIRCULAR_PER_SHOT == 1 */

for(i=1;i<=NX;i=i+IDX){
	for(j=1;j<=NY;j=j+IDY){
		waveconv[j][i] += waveconv_shot[j][i];
		waveconv_rho[j][i] += waveconv_rho_shot[j][i];
		waveconv_u[j][i] += waveconv_u_shot[j][i];
	}
}

/* saving the last time for testing */
/*sprintf(jac,"%s_after_shot%i.old.%i%i",JACOBIAN,ishot,POS[1],POS[2]);
FP3=fopen(jac,"wb");
	for (i=1;i<=NX;i=i+IDX){
		for (j=1;j<=NY;j=j+IDY){
			fwrite(&waveconv[j][i],sizeof(float),1,FP3);
		}
	}
fclose(FP3);
MPI_Barrier(MPI_COMM_WORLD);
sprintf(jac,"%s_after_shot%i.old",JACOBIAN,ishot);
if (MYID==0) mergemod(jac,3);


sprintf(jac,"%s_u_after_shot%i.old.%i%i",JACOBIAN,ishot,POS[1],POS[2]);
FP3=fopen(jac,"wb");
	for (i=1;i<=NX;i=i+IDX){
		for (j=1;j<=NY;j=j+IDY){
			fwrite(&waveconv_u[j][i],sizeof(float),1,FP3);
		}
	}
fclose(FP3);
MPI_Barrier(MPI_COMM_WORLD);
sprintf(jac,"%s_u_after_shot%i.old",JACOBIAN,ishot);
if (MYID==0) mergemod(jac,3);



sprintf(jac,"%s_rho_after_shot%i.old.%i%i",JACOBIAN,ishot,POS[1],POS[2]);
FP3=fopen(jac,"wb");
	for (i=1;i<=NX;i=i+IDX){
		for (j=1;j<=NY;j=j+IDY){
			fwrite(&waveconv_rho[j][i],sizeof(float),1,FP3);
		}
	}
fclose(FP3);
MPI_Barrier(MPI_COMM_WORLD);
sprintf(jac,"%s_rho_after_shot%i.old",JACOBIAN,ishot);
if (MYID==0) mergemod(jac,3);*/




   
} /* end of invtype == 1*/

} /* end of invmat==10 */
nsrc_loc=0;
} /* end of loop over shots (forward and backpropagation) */   

if(HESSIAN){
/* save HESSIAN for lambda */
/* ----------------------- */
sprintf(jac,"%s_hessian.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
        fwrite(&hessian[j][i],sizeof(float),1,FP4);
   }
}

fclose(FP4);

MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */ 
sprintf(jac,"%s_hessian",JACOBIAN); 
if (MYID==0) mergemod(jac,3);

/* save HESSIAN for mu */
/* ----------------------- */
sprintf(jac,"%s_hessian_u.%i%i",JACOBIAN,POS[1],POS[2]);
FP4=fopen(jac,"wb");

/* output of the gradient */
for (i=1;i<=NX;i=i+IDX){
   for (j=1;j<=NY;j=j+IDY){
        fwrite(&hessian_u[j][i],sizeof(float),1,FP4);
   }
}

fclose(FP4);

MPI_Barrier(MPI_COMM_WORLD);

/* merge gradient file */ 
sprintf(jac,"%s_hessian_u",JACOBIAN); 
if (MYID==0) mergemod(jac,3);

} /* end HESSIAN */

/* calculate L2 norm of all CPUs*/
L2sum = 0.0;
MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
energy_sum = 0.0;
MPI_Allreduce(&energy,&energy_sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
L2sum_all_shots = 0.0;
MPI_Allreduce(&L2_all_shots,&L2sum_all_shots,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
energy_sum_all_shots = 0.0;
MPI_Allreduce(&energy_all_shots,&energy_sum_all_shots,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);



if(LNORM==2){
	L2t[1]=L2sum/energy_sum;
	L2t[4]=L2sum_all_shots/energy_sum_all_shots;}
else{L2t[1]=L2sum;
     L2t[4]=L2sum_all_shots;}

if(MYID==0){
	fprintf(FP,"L2sum: %f\n", L2sum);
	fprintf(FP,"energy_sum: %e\n\n", energy_sum);
	fprintf(FP,"L2sum_all_shots: %f\n", L2sum_all_shots);
	fprintf(FP,"energy_sum_all_shots: %e\n\n", energy_sum_all_shots);}
  
if((HESSIAN!=1)&&(INVMAT==0)){
/* calculate gradient direction pi */
/* ------------------------------- */

/* interpolate unknown values */
if((IDXI>1)||(IDYI>1)){
interpol(IDXI,IDYI,waveconv,1);
}

/* calculate complete gradient */
	
        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
                          
		 waveconv_lam[j][i] = - DT * waveconv[j][i]; 
		 
		 if(INVMAT1==1){
		 /* calculate Vp gradient */ 
		 muss = prho[j][i] * pu[j][i] * pu[j][i];
		 lamss = prho[j][i] * ppi[j][i] * ppi[j][i] - 2.0 *  muss;
		 waveconv_lam[j][i] = (1.0/(4.0 * (lamss+muss) * (lamss+muss))) * waveconv_lam[j][i];
		 waveconv[j][i] = 2.0 * ppi[j][i] * prho[j][i] * waveconv_lam[j][i];
		
		 }
		 
		 if(INVMAT1==2){
		 /* calculate Zp gradient */
		 
		 waveconv[j][i] = 2.0 * ppi[j][i] * waveconv_lam[j][i];}
		 
		 if(INVMAT1==3){
		 waveconv[j][i] = waveconv_lam[j][i];}
		 
		 
		 
           }
        }

}


if((HESSIAN!=1)&&(INVMAT==0)){
/* calculate gradient direction u */
/* -------------------------------- */

/* interpolate unknown values */
if((IDXI>1)||(IDYI>1)){
interpol(IDXI,IDYI,waveconv_u,1);
}

/* calculate complete gradient */

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
	         
		 
	         /* calculate mu gradient */ 
		 waveconv_mu[j][i] = - DT * waveconv_u[j][i];
		 
		 if(INVMAT1==1){
		 /* calculate Vs gradient */		 
		 waveconv_u[j][i] = (- 4.0 * prho[j][i] * pu[j][i] * waveconv_lam[j][i]) + 2.0 * prho[j][i] * pu[j][i] * waveconv_mu[j][i];
	 
		 }
		 
		 if(INVMAT1==2){
		 /* calculate Zs gradient */
		 waveconv_u[j][i] = (- 4.0 * pu[j][i] * waveconv_lam[j][i]) + (2.0 * pu[j][i] * waveconv_mu[j][i]);}
		 
		 if(INVMAT1==3){
		 /* calculate u gradient */
		 waveconv_u[j][i] = waveconv_mu[j][i];}
		 
		 
           }
        }

}

if((HESSIAN!=1)&&(INVMAT==0)){
/* calculate gradient direction rho */
/* -------------------------------- */

/* interpolate unknown values */
if((IDXI>1)||(IDYI>1)){
interpol(IDXI,IDYI,waveconv_rho,1);
}

/* calculate complete gradient */

        for (i=1;i<=NX;i=i+IDX){
           for (j=1;j<=NY;j=j+IDY){
	    
	         /* calculate density gradient rho' */
                 waveconv_rho_s[j][i]= - DT * waveconv_rho[j][i];
		 
		 if(INVMAT1==1){
		 /* calculate density gradient */
		 waveconv_rho[j][i] = ((((ppi[j][i] * ppi[j][i])-(2.0 * pu[j][i] * pu[j][i])) * waveconv_lam[j][i]) + (pu[j][i] * pu[j][i] * waveconv_mu[j][i]) + waveconv_rho_s[j][i]);		 
		 }
		 
		 
		 if(INVMAT1==3){
		 /* calculate density gradient */
		 waveconv_rho[j][i] = waveconv_rho_s[j][i];}
		 
           }
        }

}

if((HESSIAN==0)&&(GRAD_METHOD==1)){
	PCG(waveconv, taper_coeff, nsrc, srcpos, recpos, ntr_glob, iter, C_vp, gradp, nfstart_jac, waveconv_u, C_vs, gradp_u, waveconv_rho, C_rho, gradp_rho);
}

/*if((HESSIAN==0)&&(GRAD_METHOD==2)){
    LBFGS(waveconv, taper_coeff, nsrc, srcpos, recpos, ntr_glob, iter, C_vp, gradp, nfstart_jac, waveconv_u, C_vs, gradp_u, waveconv_rho, C_rho, gradp_rho, y_LBFGS_vp, s_LBFGS_vp, rho_LBFGS, alpha_LBFGS, y_LBFGS_vs, 
          s_LBFGS_vs, y_LBFGS_rho, s_LBFGS_rho, ppi, pu, prho, nxnyi);
}*/

if((HESSIAN==0)&&(GRAD_METHOD==3)){
    LBFGS1(waveconv, taper_coeff, nsrc, srcpos, recpos, ntr_glob, iter, C_vp, gradp, nfstart_jac, waveconv_u, C_vs, gradp_u, waveconv_rho, C_rho, gradp_rho, y_LBFGS_vp, s_LBFGS_vp, rho_LBFGS_vp, rho_LBFGS_vs, 
	       alpha_LBFGS_vp, alpha_LBFGS_vs, y_LBFGS_vs, s_LBFGS_vs, y_LBFGS_rho, s_LBFGS_rho, ppi, pu, prho, nxnyi);
}

opteps_vp=0.0;
opteps_vs=0.0;
opteps_rho=0.0;

/* ============================================================================================================================*/
/* =============================================== test loop L2 ===============================================================*/
/* ============================================================================================================================*/


if((INVMAT==0) && (HESSIAN!=1)){

step1=0;
step2=0;

/* start with first guess for step length alpha */
eps_scale=0.01; /* maximum model change = 1% of the maximum model value */
countstep=0; /* count number of forward calculations */

itests=2;
iteste=2;

while((step2!=1)||(step1!=1)){

for (itest=itests;itest<=iteste;itest++){ /* calculate 3 L2 values */

/* calculate change in the material parameters */
calc_mat_change_test(waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,pu,punp1,iter,1,INVMAT,eps_scale,1,nfstart);

char modfile[STRING_SIZE];

sprintf(modfile,"%s_vp_it_countstep%d.bin",INV_MODELFILE,countstep);
writemod(modfile,ppinp1,3);

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_vs_it_countstep%d.bin",INV_MODELFILE,countstep);

writemod(modfile,punp1,3);

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0) mergemod(modfile,3);

sprintf(modfile,"%s_rho_it_countstep%d.bin",INV_MODELFILE,countstep);
writemod(modfile,prhonp1,3);

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0) mergemod(modfile,3);


/*if(MYID==0){printf("EPSILON = %e \n",EPSILON);}*/

/* For the calculation of the material parameters beteween gridpoints
   the have to be averaged. For this, values lying at 0 and NX+1,
   for example, are required on the local grid. These are now copied from the
   neighbouring grids */		
matcopy_elastic(prhonp1, ppinp1, punp1);
MPI_Barrier(MPI_COMM_WORLD);

av_mue(punp1,puipjp,prhonp1);
av_rho(prhonp1,prip,prjp);


/* initialization of L2 calculation */
L2=0.0;

alphanom = 0.0;
alphadenom = 0.0;

exchange_par();
 
if (RUN_MULTIPLE_SHOTS) nshots=nsrc; else nshots=1;
    
        for (ishot=TESTSHOT_START;ishot<=TESTSHOT_END;ishot=ishot+TESTSHOT_INCR){
 
        fprintf(FP,"\n=================================================================================================\n");
        fprintf(FP,"\n MYID=%d *****  Starting simulation (test-forward model) no. %d for shot %d of %d (parameter %d) \n",MYID,itest,ishot,nshots,partest);
		fprintf(FP,"\n=================================================================================================\n\n");
		
        for (nt=1;nt<=6;nt++) srcpos1[nt][1]=srcpos[nt][ishot]; 
		
        if (RUN_MULTIPLE_SHOTS){

	    /* find this single source positions on subdomains */
           if (nsrc_loc>0) free_matrix(srcpos_loc,1,6,1,1);
           srcpos_loc=splitsrc(srcpos1,&nsrc_loc, 1);
		}
		           
		else{
		/* Distribute multiple source positions on subdomains */
		   srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);
		}

/* calculate wavelet for each source point */
signals=wavelet(srcpos_loc,nsrc_loc);
if (nsrc_loc){
	if(QUELLART==6){FFT_filt(signals,1.0,1,ns,1);}}
		    
/* initialize wavefield with zero */
zero_fdveps(-nd+1,NY+nd,-nd+1,NX+nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs);
     
/*----------------------  loop over timesteps (forward model) ------------------*/

lsnap=iround(TSNAP1/DT);  
lsamp=NDT;
nsnap=0;

for (nt=1;nt<=NT;nt++){     
                
	/* Check if simulation is still stable */
        if (isnan(pvy[NY/2][NX/2])) err(" Simulation is unstable !");

		
   infoout = !(nt%1000);

   if (MYID==0){
      if (infoout)  fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
      time3=MPI_Wtime();
   }

   /* update of particle velocities */
   /*update_v_hc(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp,
		srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0);*/

   update_v_PML(1, NX, 1, NY, nt, pvx, pvxp1, pvxm1, pvy, pvyp1, pvym1, uttx, utty, psxx, psyy, psxy, prip, prjp, srcpos_loc,signals,signals,nsrc_loc,absorb_coeff,hc,infoout,0, K_x, a_x,
                   b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_sxx_x, psi_syy_y, psi_sxy_y, psi_sxy_x);


   if (MYID==0){
      time4=MPI_Wtime();
      time_av_v_update+=(time4-time3);
     if (infoout)  fprintf(FP," particle velocity exchange between PEs ...");
      }

   /* exchange of particle velocities between PEs */
   exchange_v(pvx, pvy, bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top, req_send, req_rec);

   if (MYID==0){
      time5=MPI_Wtime();
	time_av_v_exchange+=(time5-time4);
      if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
      }

  
   /*update_s_elastic_hc(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppinp1, punp1, puipjp, absorb_coeff, prhonp1, hc, infoout);*/

   update_s_elastic_PML(1, NX, 1, NY, pvx, pvy, ux, uy, uxy, uyx, psxx, psyy, psxy, ppinp1, punp1, puipjp, absorb_coeff, prhonp1, hc, infoout,
                        K_x, a_x, b_x, K_x_half, a_x_half, b_x_half, K_y, a_y, b_y, K_y_half, a_y_half, b_y_half, psi_vxx, psi_vyy, psi_vxy, psi_vyx);

   if (MYID==0){
      time6=MPI_Wtime();
	time_av_s_update+=(time6-time5);
      if (infoout)  fprintf(FP," stress exchange between PEs ...");
      } 
 
    /* explosive source */
   if ((!CHECKPTREAD)&&(QUELLTYP==1)) 	
   psource(nt,psxx,psyy,srcpos_loc,signals,nsrc_loc);

   if ((FREE_SURF) && (POS[2]==0))
         surface_elastic_PML(1, pvx, pvy, psxx, psyy, psxy, ppinp1, punp1, prhonp1, hc, K_x, a_x, b_x, psi_vxxs);

   /* stress exchange between PEs */
    exchange_s(psxx,psyy,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec);

   if (MYID==0){
      time7=MPI_Wtime();
 	time_av_s_exchange+=(time7-time6);
     if (infoout)  fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
      }  

	/* store amplitudes at receivers in section-arrays */
	if ((SEISMO) && (nt==lsamp) && (nt<=NT)){
		seismo_ssg(lsamp, ntr, recpos_loc, sectionvx, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppinp1, punp1, hc);
		lsamp+=NDT;
	}

      
   if (MYID==0){
      time8=MPI_Wtime();
	time_av_timestep+=(time8-time3);
      if (infoout)  fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
      }   		


   }/*--------------------  End  of loop over timesteps (test forward model) ----------*/


if (MYID==0){
printf("Calculate residuals between test forward model m - mu * dm and actual model m \n");
printf("----------------------------------------------------------------------------- \n");
}

if (ntr > 0){

/* read seismic data from SU file vx */
/* --------------------------------- */

if((QUELLTYPB==1)||(QUELLTYPB==3)){ /* if QUELLTYPB */

inseis(fprec,ishot,sectionread,ntr_glob,ns,1);

/* assign input data to each PE */
h=1;
for(i=REC1;i<=REC2;i++){
   for(j=1;j<=ns;j++){
           sectionvxdata[h][j]=sectionread[i][j];
   }
   h++;
}

L2=calc_res(sectionvxdata,sectionvx,sectionvxdiff,sectionvxdiffold,ntr,ns,LNORM,L2,0,1,1);
} /* end QUELLTYPB*/

/* read seismic data from SU file vy */
/* --------------------------------- */

if((QUELLTYPB==1)||(QUELLTYPB==2)){ /* if QUELLTYPB */
inseis(fprec,ishot,sectionread,ntr_glob,ns,2);

/* assign input data to each PE */
h=1;
for(i=REC1;i<=REC2;i++){
   for(j=1;j<=ns;j++){
           sectionvydata[h][j]=sectionread[i][j];
   }
   h++;
}

L2=calc_res(sectionvydata,sectionvy,sectionvydiff,sectionvydiffold,ntr,ns,LNORM,L2,0,1,1);			   	    

} /* end QUELLTYPB */

}

} /* ===========================================================================================================================*/   
/* ==================================== end of loop over shots (test forward) ==================================================*/
/* =============================================================================================================================*/
epst1[itest]=eps_scale;
epst1[1] = 0.0;    

L2sum=0.0;
MPI_Allreduce(&L2,&L2sum,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
if(LNORM==2){   
    L2t[itest] = L2sum/energy_sum;}
else {L2t[itest] = L2sum;}
     
} /* end of L2 test */

/* Did not found a step size which reduces the misfit function */
if((step1==0)&&(L2t[1]<L2t[2])){
 eps_scale = eps_scale/scalefac; 
 countstep++;
}

/* Found a step size with L2t[2] < L2t[3]*/
if((step1==1)&&(L2t[2]<L2t[3])){
 epst1[3]=eps_scale;
 step2=1;
}

/* Could not found a step size with L2t[2] < L2t[3]*/
if((step1==1)&&(L2t[2]>L2t[3])){
 epst1[3]=eps_scale;
 /* increase step length to find  a larger misfit function than L2t[2]*/
 eps_scale = eps_scale + (eps_scale/scalefac);
 countstep++;                       
}         

/* found a step size which reduces the misfit function */
if((step1==0)&&(L2t[1]>L2t[2])){
 epst1[2]=eps_scale; 
 step1=1;
 iteste=3;
 itests=3;
 countstep=0;
 /* find a second step length with a larger misfit function than L2t[2]*/
 eps_scale = eps_scale + (eps_scale/scalefac);
}

if((step1==0)&&(countstep>stepmax)){
  err(" Steplength estimation failed!");
  step1=1;
  step2=1;
}

if((step1==1)&&(countstep>stepmax)){
  epst1[3]=eps_scale;
  if(MYID==0){
  printf("Could not found a proper 3rd step length which brackets the minimum\n");}
  step1=1;
  step2=1;
}

if(MYID==0){printf("iteste = %d \t itests = %d \t step1 = %d \t step2 = %d \t eps_scale = %e \t countstep = %d \t stepmax= %d \t scalefac = %e \t MYID = %d \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e \n",iteste,itests,step1,step2,eps_scale,countstep,stepmax,scalefac,MYID,L2t[1],L2t[2],L2t[3]);}

} /* end of while loop */

/* calculate optimal step length epsilon for Vp and Vs*/
opteps_vp=calc_opt_step(L2t,waveconv,gradg,epst1,1,C_vp);
eps_scale = opteps_vp;

if(MYID==0){
printf("MYID = %d \t opteps_vp = %e \t opteps_vs = %e \t opteps_rho = %e \n",MYID,opteps_vp,opteps_vs,opteps_rho);
printf("MYID = %d \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e \t L2t[4] = %e \n",MYID,L2t[1],L2t[2],L2t[3],L2t[4]);
printf("MYID = %d \t epst1[1] = %e \t epst1[2] = %e \t epst1[3] = %e \n",MYID,epst1[1],epst1[2],epst1[3]);}

if(MYID==0){
fprintf(FPL2,"%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n",opteps_vp,epst1[1],epst1[2],epst1[3],L2t[1],L2t[2],L2t[3],L2t[4]);}

/* calculate optimal change in the material parameters */
calc_mat_change_test(waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,pu,punp1,iter,1,INVMAT,eps_scale,0,nfstart);

} /* end of if(INVMAT!=4) */

/* smoothing the velocity models vp and vs */
smooth_model(ppi,pu,prho,iter);

if(MYID==0){	
/*	fprintf(FPL2,"=============================================================\n");
	fprintf(FPL2,"=============================================================\n");
	fprintf(FPL2,"STATISTICS FOR ITERATION STEP %d \n",iter);
	fprintf(FPL2,"=============================================================\n");
	fprintf(FPL2,"=============================================================\n");*/
/*	fprintf(FPL2,"Low-pass filter at %e Hz\n",freq);
	fprintf(FPL2,"----------------------------------------------\n");
*/	/*fprintf(FPL2,"L2 at iteration step n = %e \n",L2);*/
/*        fprintf(FPL2,"%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n",EPSILON,EPSILON_u,EPSILON_rho,L2t[4],betaVp,betaVs,betarho,sqrt(C_vp));*/

	/*fprintf(FPL2,"----------------------------------------------\n");*/
/*	fprintf(FPL2,"=============================================================\n");
	fprintf(FPL2,"=============================================================\n\n\n");*/
}


if(MYID==0){
  fclose(FPL2);
}

/*freqshift += dfreqshift;*/

if(iter==nfstart){
nfstart = nfstart + nf;
}

if(iter==nfstart_jac){
nfstart_jac = nfstart_jac + nf_jac;
}

/* ====================================== */
} /* end of fullwaveform iteration loop*/
/* ====================================== */

if (CHECKPTWRITE){
	if (MYID==0){
		time3=MPI_Wtime();
 		fprintf(FP," Saving wavefield to check-point file %s \n",CHECKPTFILE);	
	}
	
	save_checkpoint(-1, NX+2, -1, NY+2, pvx, pvy, psxx, psyy, psxy);
	MPI_Barrier(MPI_COMM_WORLD);
	if (MYID==0){
		time4=MPI_Wtime();
      		fprintf(FP," finished (real time: %4.2f s).\n",time4-time3);
	}
}

/* deallocation of memory */
free_matrix(psxx,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(psxy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(psyy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvx,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvxp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvyp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvxm1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvym1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(ux,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(uy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(uxy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(uyx,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(Vp0,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(Vs0,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(Rho0,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(uttx,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(utty,-nd+1,NY+nd,-nd+1,NX+nd);


free_matrix(prho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(prhonp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(prip,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(prjp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pripnp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(prjpnp1,-nd+1,NY+nd,-nd+1,NX+nd);

free_matrix(ppi,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(ppinp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pu,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(punp1,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(vpmat,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(puipjp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_lam,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconvtmp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_shot,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(wcpart,1,3,1,3);

if((INVMAT==1)||(INVMAT==0)){
free_matrix(gradg,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradp,-nd+1,NY+nd,-nd+1,NX+nd);
}

free_matrix(wavejac,-nd+1,NY+nd,-nd+1,NX+nd);

if(HESSIAN){
free_vector(jac_rho,1,nxnyi*((NT/DTINV)));
free_vector(jac_u,1,nxnyi*((NT/DTINV)));
free_vector(jac_lam_x,1,nxnyi*((NT/DTINV)));
free_vector(jac_lam_y,1,nxnyi*((NT/DTINV)));
free_vector(temp_TS,1,(NT/DTINV));
free_vector(temp_TS1,1,(NT/DTINV));
free_vector(temp_TS2,1,(NT/DTINV));
free_vector(temp_TS3,1,(NT/DTINV));
free_vector(temp_TS4,1,(NT/DTINV));
free_vector(temp_TS5,1,(NT/DTINV));
free_vector(temp_conv,1,(NT/DTINV));
free_vector(temp_conv1,1,(NT/DTINV));
free_vector(temp_conv2,1,(NT/DTINV));
free_matrix(hessian,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(hessian_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(hessian_rho,-nd+1,NY+nd,-nd+1,NX+nd);
}



if(FW>0){

free_vector(d_x,1,2*FW);
free_vector(K_x,1,2*FW);
free_vector(alpha_prime_x,1,2*FW);
free_vector(a_x,1,2*FW);
free_vector(b_x,1,2*FW);

free_vector(d_x_half,1,2*FW);
free_vector(K_x_half,1,2*FW);
free_vector(alpha_prime_x_half,1,2*FW);
free_vector(a_x_half,1,2*FW);
free_vector(b_x_half,1,2*FW);

free_vector(d_y,1,2*FW);
free_vector(K_y,1,2*FW);
free_vector(alpha_prime_y,1,2*FW);
free_vector(a_y,1,2*FW);
free_vector(b_y,1,2*FW);

free_vector(d_y_half,1,2*FW);
free_vector(K_y_half,1,2*FW);
free_vector(alpha_prime_y_half,1,2*FW);
free_vector(a_y_half,1,2*FW);
free_vector(b_y_half,1,2*FW);

free_matrix(psi_sxx_x,1,NY,1,2*FW);
free_matrix(psi_syy_y,1,2*FW,1,NX);
free_matrix(psi_sxy_x,1,NY,1,2*FW);
free_matrix(psi_sxy_y,1,2*FW,1,NX);
free_matrix(psi_vxx,1,NY,1,2*FW);
free_matrix(psi_vxxs,1,NY,1,2*FW);
free_matrix(psi_vyy,1,2*FW,1,NX);
free_matrix(psi_vxy,1,2*FW,1,NX);
free_matrix(psi_vyx,1,NY,1,2*FW);

}

/*absorb_coeff=  matrix(1,NY,1,NX);
free_matrix(absorb_coeff,1,NY,1,NX);*/
free_matrix(taper_coeff,1,NY,1,NX);

free_matrix(bufferlef_to_rig,1,NY,1,fdo3);
free_matrix(bufferrig_to_lef,1,NY,1,fdo3);
free_matrix(buffertop_to_bot,1,NX,1,fdo3);
free_matrix(bufferbot_to_top,1,NX,1,fdo3);

free_vector(hc,0,6);

if((INVMAT==1)||(INVMAT==0)){
free_vector(forward_prop_x,1,NY*NX*NT);
free_vector(forward_prop_y,1,NY*NX*NT);
}

if((INVMAT==2)||(INVMAT==0)){
free_vector(forward_prop_rho_x,1,NY*NX*NT);
free_vector(forward_prop_rho_y,1,NY*NX*NT);
free_matrix(gradg_rho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradp_rho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_rho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_rho_s,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_rho_shot,-nd+1,NY+nd,-nd+1,NX+nd);
}

if((INVMAT==3)||(INVMAT==0)){
free_vector(forward_prop_u,1,NY*NX*NT);
free_matrix(gradg_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(gradp_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_u,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_mu,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(waveconv_u_shot,-nd+1,NY+nd,-nd+1,NX+nd);
}

if (nsrc_loc>0){	
	free_matrix(signals,1,nsrc_loc,1,NT);
	/*dsignals=fmatrix(1,nsrc,1,NT);
	free_matrix(dsignals,1,nsrc_loc,1,NT);*/
	free_matrix(srcpos_loc,1,6,1,nsrc_loc);
	free_matrix(srcpos_loc_back,1,6,1,nsrc_loc);
}
		   

if (SEISMO) free_imatrix(recpos,1,3,1,ntr_glob);

/* free memory for global source positions */
free_matrix(srcpos,1,6,1,nsrc);

if ((ntr>0) && (SEISMO)){	

      free_imatrix(recpos_loc,1,3,1,ntr);
	switch (SEISMO){
	case 1 : /* particle velocities only */
		free_matrix(sectionvx,1,ntr,1,ns);
		free_matrix(sectionvy,1,ntr,1,ns);		
		break;	
	case 2 : /* pressure only */
		free_matrix(sectionp,1,ntr,1,ns);
		free_matrix(sectionpn,1,ntr,1,ns);
		free_matrix(sectionpnp1,1,ntr,1,ns);
		break;	
	case 3 : /* curl and div only */
		free_matrix(sectioncurl,1,ntr,1,ns);
		free_matrix(sectiondiv,1,ntr,1,ns);
		break;	
	case 4 : /* everything */
		free_matrix(sectionvx,1,ntr,1,ns);
		free_matrix(sectionvy,1,ntr,1,ns);
		free_matrix(sectionp,1,ntr,1,ns);
		free_matrix(sectioncurl,1,ntr,1,ns);
		free_matrix(sectiondiv,1,ntr,1,ns);		
		break;
	}	

}	

 /* free memory for source position definition */
 free_matrix(srcpos1,1,6,1,1);
 
 free_matrix(sectionvxdata,1,ntr,1,ns);
 free_matrix(sectionread,1,ntr_glob,1,ns);
 free_matrix(sectionvxdiff,1,ntr,1,ns);
 free_matrix(sectionvydata,1,ntr,1,ns);
 free_matrix(sectionvydiff,1,ntr,1,ns);	
 free_matrix(sectionvydiffold,1,ntr,1,ns);
 free_matrix(sectionvxdiffold,1,ntr,1,ns);		
 free_vector(L2t,1,4);
 free_vector(epst1,1,3);
 free_vector(epst2,1,3);
 free_vector(epst3,1,3); 
 free_vector(picked_times,1,ntr);
 
 /* free memory for inversion of source time function */
 free_matrix(sectionvy_conv,1,ntr_glob,1,ns);
 free_matrix(sectionvy_syn,1,ntr,1,ns);
 free_matrix(sectionvy_obs,1,ntr_glob,1,ns);
 free_vector(source_time_function,1,ns);
 free_matrix(signals_stf,1,nsrc_glob,1,TIME/DT);
 
 free_ivector(recswitch,1,ntr_glob);
 free_matrix(fulldata,1,ntr_glob,1,NT); 
 free_matrix(fulldata_vx,1,ntr_glob,1,NT);
 free_matrix(fulldata_vy,1,ntr_glob,1,NT);
 
 /* de-allocate buffer for messages */
MPI_Buffer_detach(buff_addr,&buffsize);
	
/* merge snapshot files created by the PEs into one file */
/* if ((SNAP) && (MYID==0)){ 
	snapmerge(nsnap);
}
*/

MPI_Barrier(MPI_COMM_WORLD);

if (MYID==0){
	fprintf(FP,"\n **Info from main (written by PE %d): \n",MYID);
	fprintf(FP," CPU time of program per PE: %li seconds.\n",clock()/CLOCKS_PER_SEC);
	time8=MPI_Wtime();
	fprintf(FP," Total real time of program: %4.2f seconds.\n",time8-time1);
	time_av_v_update=time_av_v_update/(double)NT;
	time_av_s_update=time_av_s_update/(double)NT;
	time_av_v_exchange=time_av_v_exchange/(double)NT;
	time_av_s_exchange=time_av_s_exchange/(double)NT;
	time_av_timestep=time_av_timestep/(double)NT;
	fprintf(FP," Average times for \n");
	fprintf(FP," velocity update:  \t %5.3f seconds  \n",time_av_v_update);
	fprintf(FP," stress update:  \t %5.3f seconds  \n",time_av_s_update);
	fprintf(FP," velocity exchange:  \t %5.3f seconds  \n",time_av_v_exchange);
	fprintf(FP," stress exchange:  \t %5.3f seconds  \n",time_av_s_exchange);
	fprintf(FP," timestep:  \t %5.3f seconds  \n",time_av_timestep);
		
}

fclose(FP);

MPI_Finalize();
return 0;	

}/*main*/
