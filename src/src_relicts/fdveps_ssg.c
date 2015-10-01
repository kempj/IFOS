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

/* $Id: fdveps_ssg.c,v 1.1.1.1 2007/11/21 22:44:52 koehn Exp $*/
/*------------------------------------------------------------------------
 *  Parallel 2-D-Viscoelastic Finite Difference Modelling 
 *  using the Standard Staggered Grid (SSG)
 *
 *
 *  
 *
 * 
 * In case of questions contact the author:                       
 *	Dr. Thomas Bohlen, Kiel University, Institute of Geoscience,
 *	Otto-Hahn-Platz 1, D-24098 Kiel, Germany, ph: +49 431 880 4648,
 *	fax: +49 431 880 4432, mailto:tbohlen@geophysik.uni-kiel.de,
 *	Homepage: http://www.geophysik.uni-kiel.de/~tbohlen
 *
 *
 *  Do not freely distribute the program. In case of interest of other
 *  people in obtaining the source code please refer them to the author.
 *  I would like to keep track where the program is used and modified.
 *  If you show modelling results in a paper or presentation please give a reference
 *  to the following paper:
 *  Bohlen, 2002, Parallel 3-D viscoelastic finite-difference seismic modelling,
 *  Computers and Geociences, 28, 887-899.
 *  
 *  Thank you for your co-operation, 
 *  Thomas Bohlen
 

	History of major modifications:

	Febr. 2002	Version 1.0	original parallel implementation
					T. Bohlen
					
	Febr. 2002  Version 1.1
					- fixed bugs which occured when reading
					  source signals from ASCII file
					- ported the coded to CRAY T3E system (changed outseis.c
					   to output correct SU format on CRAY T3E)
					- changed storage of density in field prho[j][i]:
					  it does now store density instead of 1/density ! 
					T. Bohlen

	April 2002  Version 1.5
					- dipping incident plane wave implemented
					T. Bohlen
					
	May 2002  Version 1.6
					- improved merging of snapshots (no memory allocation
					  required), output format of snapshot data was modified
					T. Bohlen

	July 2002  Version 1.7
					- snapshots files can now be merged to one data file using
					  program snapmerge. Usage: 'snapmerge < fdveps.inp' after
					  fdveps has finished. The merged data file can be visualized
					  with the SU program xmovie.
					- version for anisotropic media implemeted and tested.
					T. Bohlen

	June 2003 Version 1.9		- modification of  locations of wavefield and material parameters
	   				  on the standard staggered grid according to Vossen et al.,2002,
					  Geophysics, 67, No. 2, 618-624. In the new implementation only
					  density and shear slowness are averaged. This distribution has
					  higher accuracy in the presence of
					  fluid/solid contrasts (seafloor) compared to the
					  distribution I have used previously (Bohlen, Phd thesis,
					  1998)
	December, 2004, Version 2.0	- Since Revision 2.0 update with CVS

 *  ----------------------------------------------------------------------*/

#include "fd.h"           /* general include file for viscoelastic FD programs */

#include "globvar.h"      /* definition of global variables  */


int main(int argc, char **argv){
/* variables in main */
int ns, nseismograms=0, nt, nd, fdo3;
int lsnap, nsnap=0, lsamp=0, buffsize;
int ntr=0, ntr_loc=0, ntr_glob=0, nsrc=0, nsrc_loc=0;

float memdyn, memmodel, memseismograms, membuffer, memtotal;
float fac1, fac2;
char *buff_addr, ext[10], *fileinp;
double 	time1, time2, time3, time4, time5, time6, time7, time8,
	time_av_v_update=0.0, time_av_s_update=0.0, time_av_v_exchange=0.0, 
	time_av_s_exchange=0.0, time_av_timestep=0.0;


float  **  psxx, **  psxy, **  psyy;
float  **  pvx, **  pvy, ***  pr;
float  ***  pp, ***  pq;
float  **  prho, **  ppi, **  pu, **puipjp, **ptausipjp;
float  **  ptaus, **  ptaup, *  peta;
float * etaip, * etajm;

float  ** sectionvx=NULL, ** sectionvy=NULL, ** sectionp=NULL,
	** sectioncurl=NULL, ** sectiondiv=NULL;
float  **  absorb_coeff;
float  ** srcpos=NULL, **srcpos_loc=NULL, ** signals=NULL, *hc=NULL;
int   **recpos=NULL, ** recpos_loc=NULL;

float ** bufferlef_to_rig,  ** bufferrig_to_lef,
** buffertop_to_bot, ** bufferbot_to_top; 

	
MPI_Request *req_send, *req_rec;
MPI_Status  *send_statuses, *rec_statuses;



/* Initialize MPI environment */
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&NP);
MPI_Comm_rank(MPI_COMM_WORLD,&MYID);

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
/*fclose(FP);*/
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
   recpos_loc = splitrec(recpos,&ntr_loc, ntr);
   ntr_glob=ntr;
   ntr=ntr_loc;
}



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
	


/*allocate memory for dynamic, static and buffer arrays */
fac1=(NX+4.0)*(NY+4.0);
fac2=sizeof(float)*pow(2.0,-20.0);

nd = FDORDER/2 + 1;
fdo3 = 2*nd;

memdyn=(5.0+3.0*(float)L)*fac1*fac2;
memmodel=9.0*fac1*fac2;
memseismograms=nseismograms*ntr*ns*fac2;
membuffer=2.0*fdo3*(NY+NX)*fac2;
buffsize=2.0*2.0*fdo3*(NX +NY)*sizeof(MPI_FLOAT);
memtotal=memdyn+memmodel+memseismograms+membuffer+(buffsize*pow(2.0,-20.0));


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
pr   =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
pp   =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
pq   =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
	
/* memory allocation for static (model) arrays */
prho         =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
ppi          =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
pu           =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
puipjp       =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
ptaus        =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
ptausipjp    =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
ptaup        =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
peta         =  vector(1,L);
etaip        =  vector(1,L);
etajm        =  vector(1,L);
absorb_coeff =  matrix(1,NY,1,NX);

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


fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);


/* Holberg coefficients for FD operators*/
hc = holbergcoeff();

/* Reading source positions from SOURCE_FILE */ 	
srcpos=sources(&nsrc);

/* Distribute source positions on subdomains */
srcpos_loc = splitsrc(srcpos,&nsrc_loc, nsrc);


/* calculate wavelet for each source point */
signals=wavelet(srcpos_loc,nsrc_loc);

/* output source signal */
/*fp=fopen("source_signal.asc","w");
fprintf(fp,"%i\n",NT);
for (nt=1;nt<=NT;nt++){
	fprintf(fp,"%e ",signals[1][nt]);
}
fclose(fp);
*/



MPI_Barrier(MPI_COMM_WORLD);


/* create model grids */

if (READMOD) readmod(prho,ppi,pu,ptaus,ptaup,peta);
    else model(prho,ppi,pu,ptaus,ptaup,peta);


/* check if the FD run will be stable and free of numerical dispersion */
checkfd_hc(FP,prho,ppi,pu,ptaus,ptaup,peta,hc);

/* calculate 2-D array for exponential damping of reflections
   at the edges of the numerical mesh */
absorb(absorb_coeff);



/* For the calculation of the material parameters beteween gridpoints
	   the have to be averaged. For this, values lying at 0 and NX+1,
	for example, are required on the local grid. These are now copied from the
	neighbouring grids */		
matcopy(prho, ppi, pu, ptaus, ptaup);

av_mue(pu,puipjp);
av_tau(ptaus,ptausipjp);
 

MPI_Barrier(MPI_COMM_WORLD);


/* comunication initialisation for persistent communication */
comm_ini(bufferlef_to_rig, bufferrig_to_lef, buffertop_to_bot, bufferbot_to_top, 
	req_send, req_rec);


if (MYID==0)
   {
   time2=MPI_Wtime();
   fprintf(FP,"\n\n\n *********** STARTING TIME STEPPING ***************\n");
   fprintf(FP," real time before starting time loop: %4.2f s.\n",time2-time1);
   }

/*----------------------  loop over timesteps  ------------------*/

for (nt=1;nt<=NT;nt++){     
		
   if (MYID==0){
      fprintf(FP,"\n Computing timestep %d of %d \n",nt,NT);
      /*printf(" Computing timestep %d of %d \n",nt,NT);*/
      time3=MPI_Wtime();
   }

		
   /* update of particle velocities */

   update_v_hc(1, NX, 1, NY, nt, pvx, pvy, psxx, psyy, psxy, prho,
   	    srcpos_loc,signals,nsrc_loc,absorb_coeff,hc); 


   if (MYID==0){
      time4=MPI_Wtime();
      time_av_v_update+=(time4-time3);
     fprintf(FP," particle velocity exchange between PEs ...");
      }

   /* exchange of particle velocities between PEs */
   exchange_v(pvx,pvy, bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top, req_send, req_rec);

   if (MYID==0){
      time5=MPI_Wtime();
	time_av_v_exchange+=(time5-time4);
      fprintf(FP," finished (real time: %4.2f s).\n",time5-time4);
      }


   /* stress update */
   update_s_visc_hc(1, NX, 1, NY, pvx, pvy, psxx, psyy,
      psxy, pr, pp, pq, ppi, pu, puipjp, ptaup, ptaus, ptausipjp, etaip,
      etajm, peta, hc);
	
   
   /* explosive source */
   if (QUELLTYP==1) 	
   psource(nt,psxx,psyy,srcpos_loc,signals,nsrc_loc);


   if ((FREE_SURF) && (POS[2]==0))
      surface(1, pvx, pvy, psxx, psyy,
	 psxy, pp, pq, ppi, pu, ptaup, ptaus, etajm, peta);	


   if (MYID==0){
      time6=MPI_Wtime();
	time_av_s_update+=(time6-time5);
      fprintf(FP," stress exchange between PEs ...");
      }


   /* stress exchange between PEs */
    exchange_s(psxx,psyy,psxy, 
      bufferlef_to_rig, bufferrig_to_lef, 
      buffertop_to_bot, bufferbot_to_top,
      req_send, req_rec);

 
   if (MYID==0){
      time7=MPI_Wtime();
 	time_av_s_exchange+=(time7-time6);
     fprintf(FP," finished (real time: %4.2f s).\n",time7-time6);
      }


   /* store amplitudes at receivers in section-arrays */
   if ((SEISMO) && (nt==lsamp) && (nt<NT)){
		seismo_ssg(lsamp, ntr, recpos_loc, sectionvx, sectionvy, 
			sectionp, sectioncurl, sectiondiv, 
			pvx, pvy, psxx, psyy, ppi, pu, hc);
		lsamp+=NDT;
   }



   /* WRITE SNAPSHOTS TO DISK */
   if ((SNAP) && (nt==lsnap) && (nt<=TSNAP2/DT)){
      snap(FP,nt,++nsnap,pvx,pvy,psxx,psyy,pu,ppi,hc);
      lsnap=lsnap+iround(TSNAPINC/DT);
      }

   if (MYID==0){
      time8=MPI_Wtime();
	time_av_timestep+=(time8-time3);
      fprintf(FP," total real time for timestep %d : %4.2f s.\n",nt,time8-time3);
      }   		


   }/*--------------------  End  of loop over timesteps ----------*/

	
	
   if ((ntr > 0) && (SEISMO)){
	saveseis(FP,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,recpos,recpos_loc,ntr,
	          srcpos,nsrc,ns);
	}




/* deallocation of memory */
free_matrix(psxx,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(psxy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(psyy,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvx,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pvy,-nd+1,NY+nd,-nd+1,NX+nd);
free_f3tensor(pr,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
free_f3tensor(pp,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
free_f3tensor(pq,-nd+1,NY+nd,-nd+1,NX+nd,1,L);
free_matrix(prho,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(ppi,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(pu,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(puipjp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(ptaus,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(ptausipjp,-nd+1,NY+nd,-nd+1,NX+nd);
free_matrix(ptaup,-nd+1,NY+nd,-nd+1,NX+nd);
free_vector(peta,1,L);
free_vector(etaip,1,L);
free_vector(etajm,1,L);

free_matrix(absorb_coeff,1,NY,1,NX);

free_matrix(bufferlef_to_rig,1,NY,1,fdo3);
free_matrix(bufferrig_to_lef,1,NY,1,fdo3);
free_matrix(buffertop_to_bot,1,NX,1,fdo3);
free_matrix(bufferbot_to_top,1,NX,1,fdo3);
	
if (nsrc_loc>0){	
	free_matrix(signals,1,nsrc_loc,1,NT);
	free_matrix(srcpos_loc,1,6,1,nsrc_loc);
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
