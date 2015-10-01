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
 *   inversion for source time function 
 *   31. August 2011 L. Rehor, T. Forbriger, M. Sch�fer
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "stfinv/stfinv.h"
#include "segy.h"

void stf(FILE *fp, float **sectionvy, float ** sectionvy_obs, float ** sectionvy_conv, float * source_time_function, int  **recpos, int  **recpos_loc, 
int ntr_glob, int ntr,float ** srcpos, int ishot, int ns, int iter, int nshots, float ** signals_stf, float ** signals){ 

	/* declaration of global variables */
	extern float DT, DH;
	extern int SEIS_FORMAT, MYID, NT;
	extern char  SEIS_FILE_VY[STRING_SIZE], PARA[STRING_SIZE];
	
	/* --------------- declaration of variables --------------- */
	unsigned int nrec, nsamp, i, t;
	float dt;
	float xr=0.0, yr=0.0;
	float XS=0.0, YS=0.0;
	char conv[STRING_SIZE], qw[STRING_SIZE], qw_conv[STRING_SIZE];
	int nts, tt, k=1;
	float * signals_temp=NULL, * stf_conv=NULL, * stf_conv_temp=NULL;
	
	signals_temp=vector(1,NT);
	stf_conv=vector(1,NT);
	stf_conv_temp=vector(1,NT);
	
	printf("\n==================================================================================\n\n");
	printf("\n ***** Inversion of Source Time Function - shot: %d - it: %d ***** \n\n",ishot,iter);
				
	nrec=(unsigned int)ntr_glob;
	nsamp=(unsigned int)ns;
	dt=DT;
	
	
	/* source coordinates are written into trace header fields */
	XS=srcpos[1][ishot];	
	YS=srcpos[2][ishot];
		
	struct CTriples data;
	data.n=nrec;
	data.triples=(struct CWaveformTriple *)malloc(nrec*sizeof(struct CWaveformTriple));
	if (data.triples == NULL) {abort();}
	for (i=0;i<nrec;i++){
	
		xr=recpos[1][i+1]*DH;
		yr=recpos[2][i+1]*DH;
		
		data.triples[i].data=&sectionvy_obs[i+1][1];
		
		data.triples[i].synthetics=&sectionvy[i+1][1];
		
		data.triples[i].convolvedsynthetics=&sectionvy_conv[i+1][1];
		
		data.triples[i].header.sx=(unsigned int)iround(XS*1000.0);  /* X source coordinate */
		data.triples[i].header.sy=0.0;
		data.triples[i].header.sz=(unsigned int)iround(YS*1000.0);  /* source depth (positive) */
		data.triples[i].header.rx=(unsigned int)iround(xr*1000.0);  /* group coordinates */
		data.triples[i].header.ry=0.0;
		data.triples[i].header.rz=(unsigned int)iround(yr*1000.0);
		data.triples[i].header.sampling.n=nsamp;
		data.triples[i].header.sampling.dt=dt;
	}
	
	struct CWaveform stf;
	stf.series = &source_time_function[1];
	stf.sampling.n=nsamp;
	stf.sampling.dt=dt;
	
	
	/*char para[]="fbd:tshift=0.0"; /* parameter string */
	
	initstfinvengine(data, stf, PARA);
		
	runstfinvengine();
	
	/* --------------- writing out the convolved seismograms --------------- */
	
	sprintf(conv,"%s.conv.shot%d_it%d",SEIS_FILE_VY,ishot,iter);
	printf(" PE %d is writing %d convolved seismograms (vy) for shot = %d to\n\t %s \n",MYID,ntr_glob,ishot,conv);
	outseis_glob(fp,fopen(conv,"w"),1,sectionvy_conv,recpos,recpos_loc,ntr_glob,srcpos,0,ns,SEIS_FORMAT,ishot);

	/* --------------- writing out the source time function --------------- */
	
	sprintf(qw,"%s.stf.shot%d_it%d",SEIS_FILE_VY,ishot,iter);
	printf(" PE %d is writing source time function for shot = %d to\n\t %s \n",MYID,ishot,qw);
	outseis_vector(fp,fopen(qw,"w"),1,source_time_function,recpos,recpos_loc,ntr,srcpos,0,ns,SEIS_FORMAT,ishot);
	
		
	/*-------- convolving source time function with synthetic wavelet -------------*/
	
	for (i=0;i<NT;i++){ signals_temp[i]=signals[1][i]; }
				
	conv_FD(signals_temp,source_time_function,stf_conv_temp,NT);
	
	for (i=0;i<NT;i++){ stf_conv[i]=stf_conv_temp[i]/NT; }
	
	sprintf(qw_conv,"%s.stf_conv.shot%d_it%d",SEIS_FILE_VY,ishot,iter);
	printf(" PE %d is writing source time function for shot = %d to\n\t %s \n",MYID,ishot,qw_conv);
	outseis_vector(fp,fopen(qw_conv,"w"),1,stf_conv,recpos,recpos_loc,ntr,srcpos,0,ns,SEIS_FORMAT,ishot);
	
	
	/* ---------------------- filtering source time function cause of grid disperion due to high frequencies (is this helpful?????) --------*/
	
	/*timedomain_filt_vector(source_time_function,30,2,ntr_glob,ns,1);*/
	
	
	/* --------------- writing out the source time function in the same way as wavelet.c--------------- */
	
	for (tt=1;tt<=NT;tt++){
			signals_stf[ishot][tt]=stf_conv[tt];	
			        }	
				
	printf("\n\n==================================================================================\n");
	
	free_vector(signals_temp,1,NT);
	free_vector(stf_conv,1,NT);
	free_vector(stf_conv_temp,1,NT);

}
