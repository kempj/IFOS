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

/*------------------------------------------------------------------------
 *   inversion for source time function 
 *   31. August 2011 L. Rehor, T. Forbriger, M. Schaefer
 *  ----------------------------------------------------------------------*/

#include "fd.h"
#include "stfinv/stfinv.h"
#include "segy.h"

void stf(FILE *fp, float **sectionvy, float ** sectionvy_obs, float ** sectionvy_conv, float * source_time_function, int  **recpos, int  **recpos_loc, 
int ntr_glob,int ntr, float ** srcpos, int ishot, int ns, int iter, int nshots, float FC){ 

	/* declaration of global variables */
	extern float DT, DH;
	extern int SEIS_FORMAT, MYID, NT, SOURCE_SHAPE, TIME_FILT, TIMEWIN, TAPER_STF;
	extern char  PARA[STRING_SIZE], DATA_DIR[STRING_SIZE];
	extern int TRKILL_STF, NORMALIZE, USE_WORKFLOW, WORKFLOW_STAGE;
	extern char TRKILL_FILE_STF[STRING_SIZE];
	extern char SIGNAL_FILE[STRING_SIZE];
	
	/* declaration of variables for trace killing */
	int ** kill_tmp, *kill_vector, h, j;
	char trace_kill_file[STRING_SIZE];
	FILE *ftracekill;
	
	/* --------------- declaration of variables --------------- */
	unsigned int nrec, nsamp, i, npairs;
	float dt;
	float xr=0.0, yr=0.0;
	float XS=0.0, YS=0.0;
	char conv_y[STRING_SIZE], qw[STRING_SIZE], conv_y_tmp[STRING_SIZE], obs_y_tmp[STRING_SIZE], mod_y_tmp[STRING_SIZE];
	
	/* variables for wavelet */
	int nt, nts;
	float tshift, amp=0.0, fc, tau, t, ts, ag;
	float * wavelet, * stf_conv_wavelet, *psource=NULL;
	
	wavelet=vector(1,ns);
	stf_conv_wavelet=vector(1,ns);
	
	printf("\n================================================================================================\n\n");
	printf("\n ***** Inversion of Source Time Function - shot: %d - it: %d ***** \n\n",ishot,iter);
	
	/* if TRKILL_STF==1 a trace killing is applied */
	if(TRKILL_STF){
		kill_tmp = imatrix(1,ntr_glob,1,nshots);
		kill_vector = ivector(1,ntr_glob);
		
        if(USE_WORKFLOW){
            sprintf(trace_kill_file,"%s_%i.dat",TRKILL_FILE_STF,WORKFLOW_STAGE);
            ftracekill=fopen(trace_kill_file,"r");
            if (ftracekill==NULL){
                sprintf(trace_kill_file,"%s.dat",TRKILL_FILE_STF);
                ftracekill=fopen(trace_kill_file,"r");
                if (ftracekill==NULL){
                    err(" Trace kill file could not be opened!");
                }
            }
        }else{
            sprintf(trace_kill_file,"%s.dat",TRKILL_FILE_STF);
            ftracekill=fopen(trace_kill_file,"r");
            if (ftracekill==NULL){
                err(" Trace kill file could not be opened!");
            }
        }
		
		for(i=1;i<=ntr_glob;i++){
			for(j=1;j<=nshots;j++){
				fscanf(ftracekill,"%d",&kill_tmp[i][j]);
			}
		}
		
		fclose(ftracekill);
		
		h=1;
		for(i=1;i<=ntr_glob;i++){
			kill_vector[h] = kill_tmp[i][ishot];
			h++;
		}
	} /* end if(TRKILL_STF)*/
	
	if(TRKILL_STF){
		for(i=1;i<=ntr_glob;i++){
		if(i==1)printf("\n ***** Trace killing is applied for trace: ***** \n ***** \t");
			if(kill_vector[i]==1){
				printf("%d \t",i);
				for(j=1;j<=ns;j++){
					sectionvy[i][j]=0.0;
					sectionvy_obs[i][j]=0.0;
				}
			}
		if(i==ntr_glob)printf(" ***** \n\n");
		}
	}
	/* trace killing ends here */
	
	if(TIMEWIN==1){
		time_window_glob(sectionvy, iter, ntr_glob, ns, ishot);
		time_window_glob(sectionvy_obs, iter, ntr_glob, ns, ishot);
	}
	
	/* NORMALIZE TRACES */
	if(NORMALIZE==1){
		normalize_data(sectionvy,ntr_glob,ns);
		normalize_data(sectionvy_obs,ntr_glob,ns);
	}
	
	nrec=(unsigned int)ntr_glob;
	nsamp=(unsigned int)ns;
	dt=DT;
	npairs=nrec;
	
	/* source coordinates are written into trace header fields */
	XS=srcpos[1][ishot];	
	YS=srcpos[2][ishot];
	
	
	/* TF Software: see libstfinv */
	
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
	
	
	initstfinvengine(data, stf, PARA);
	
	runstfinvengine();
	
	/* END TF Software */
	
	
	psource=vector(1,ns);
	
	if (SOURCE_SHAPE==3) psource=rd_sour(&nts,fopen(SIGNAL_FILE,"r"));
	if (SOURCE_SHAPE==7){
		inseis_source_wavelet(psource,ns,ishot);
	}
	
	/* calculating wavelet SIN**3 for convoling with STF */
	tshift=srcpos[4][ishot];
	fc=srcpos[5][ishot];
	ts=1.0/fc;
	for (nt=1;nt<=ns;nt++){
		t=(float)nt*DT;
		switch (SOURCE_SHAPE){
					case 1 : 
						/* Old Ricker Wavelet */
						/* tau=PI*(t-ts-tshift)/(1.5*ts);
						amp=(((1.0-4.0*tau*tau)*exp(-2.0*tau*tau))); */ 

						/* New Ricker Wavelet, equal to SOFI2D */
						tau=PI*(t-1.5*ts-tshift)/(ts);
						amp=(((1.0-2.0*tau*tau)*exp(-tau*tau)));
					break;
					case 2 : 
						if ((t<tshift) || (t>(tshift+ts))) amp=0.0;
						else amp=((sin(2.0*PI*(t-tshift)*fc) 
			    				-0.5*sin(4.0*PI*(t-tshift)*fc)));

						/*amp=((sin(2.0*PI*(t+tshift)*fc) 
			    			-0.5*sin(4.0*PI*(t+tshift)*fc)));*/
					break;
					case 3 : 
						if (nt<=nts) amp=psource[nt]; 
						else amp=0.0;
// 						amp=psource[nt];
					break;  /* source wavelet from file SOURCE_FILE */
					case 4 : 
						/*tau=PI*(t-ts-tshift)/(1.5*ts);*/ /* Ricker */
						/*amp=((t-ts-tshift)*exp(-2.0*tau*tau));*/
						
						if ((t<tshift) || (t>(tshift+ts))) amp=0.0;
						else amp=pow(sin(PI*(t+tshift)/ts),3.0);
						break; /* sinus raised to the power of three */
						
					break; 
					case 5 : 
				                /* first derivative of a Gaussian */
				                ts=1.2/fc;
					        ag  = PI*PI*fc*fc;
			                        amp = - 2.0 * ag * (t-ts) * exp(-ag*(t-ts)*(t-ts));
					break;					
					case 6 : 
					        /* Bandlimited Spike */
						amp=0.0;
						if(nt==1+iround(tshift/DT)){
						amp = 1.0;}
					break; 
					case 7 :
					        /* source wavelet from file SOURCE_FILE */ 
						amp=psource[nt];
					break;  
					case 8 :
					        /* integral of sinus raised to the power of three */
						if (t<tshift) {
							 amp=0.0;}
						if ((t>=tshift) && (t<=(tshift+ts))){
							amp=(ts/(0.75*PI))*(0.5-0.75*cos(PI*(t-tshift)/ts)+0.25*pow(cos(PI*(t-tshift)/ts),3.0));}
						if (t>(tshift+ts))
							{amp=ts/(0.75*PI);}
						break;                                                                                                                                           	
					default : 
						err("Which source-wavelet ? ");
					
					
		}/* end of switch (SOURCE_SHAPE)  */
		wavelet[nt]=amp;
		
	}/*  end of for (nt=1;nt<=ns;nt++) */	
	
	
	/* convolving wavelet with STF */
	conv_FD(wavelet,source_time_function,stf_conv_wavelet,ns);
	
	if(TAPER_STF)
		taper(stf_conv_wavelet, ns, fc);
	
// 	/* --------------- writing out the observed seismograms --------------- */
// 	sprintf(obs_y_tmp,"%s.shot%d.obs",SEIS_FILE_P,ishot);
// 	printf(" PE %d is writing %d observed seismograms (vy) for shot = %d to\n\t %s \n",MYID,ntr_glob,ishot,obs_y_tmp);
// 	outseis_glob(fp,fopen(obs_y_tmp,"w"),1,sectionvy_obs,recpos,recpos_loc,ntr_glob,srcpos,0,ns,SEIS_FORMAT,ishot,0);
// 	
// 	/* --------------- writing out the modelled seismograms --------------- */
// 	sprintf(mod_y_tmp,"%s.shot%d.mod",SEIS_FILE_P,ishot);
// 	printf(" PE %d is writing %d modelled seismograms (vy) for shot = %d to\n\t %s \n",MYID,ntr_glob,ishot,mod_y_tmp);
// 	outseis_glob(fp,fopen(mod_y_tmp,"w"),1,sectionvy,recpos,recpos_loc,ntr_glob,srcpos,0,ns,SEIS_FORMAT,ishot,0);
	
	
	/* --------------- writing out the source time function --------------- */
	if((TIME_FILT==1)||(TIME_FILT==2)){
		sprintf(qw,"%s.shot%d_%dHz",SIGNAL_FILE,ishot,(int)FC);
		printf(" PE %d is writing source time function for shot = %d to\n\t %s \n",MYID,ishot,qw);
		outseis_vector(fp,fopen(qw,"w"),1,stf_conv_wavelet,recpos,recpos_loc,ntr,srcpos,0,ns,SEIS_FORMAT,ishot,0);
	}
	
	sprintf(qw,"%s.shot%d",SIGNAL_FILE,ishot);
	printf(" PE %d is writing source time function for shot = %d to\n\t %s \n",MYID,ishot,qw);
	outseis_vector(fp,fopen(qw,"w"),1,stf_conv_wavelet,recpos,recpos_loc,ntr,srcpos,0,ns,SEIS_FORMAT,ishot,0);
	
	/*sprintf(conv_y_tmp,"%s.shot%d.forward",SEIS_FILE_VY,ishot);
	printf(" PE %d is writing %d seismograms (vy) for shot = %d to\n\t %s \n",MYID,ntr_glob,ishot,conv_y_tmp);
	outseis_glob(fp,fopen(conv_y_tmp,"w"),1,sectionvy,recpos,recpos_loc,ntr_glob,srcpos,0,ns,SEIS_FORMAT,ishot,0);*/
	
	/*freestfinvengine();
	free(data.triples);*/
	
	free_vector(wavelet,1,ns);
	free_vector(stf_conv_wavelet,1,ns);
	free_vector(psource,1,ns);
	
	/* free memory for trace killing */
	if(TRKILL_STF){
		free_imatrix(kill_tmp,1,ntr_glob,1,nshots);
		free_ivector(kill_vector,1,ntr_glob);
	}
}