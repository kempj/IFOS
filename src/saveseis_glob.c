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
 *   write seismograms to files 
 *   last update Oct 2015 F.Wittkamp
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void saveseis_glob(FILE *fp, float **sectionvx, float **sectionvy,float **sectionvz,float **sectionp,
float **sectioncurl, float **sectiondiv, int  **recpos, int  **recpos_loc, 
int ntr, float ** srcpos, int ishot, int ns, int iter){ 
		
	extern int SEISMO, SEIS_FORMAT, RUN_MULTIPLE_SHOTS, WAVETYPE, VERBOSE;
	extern char  SEIS_FILE_VX[STRING_SIZE], SEIS_FILE_VY[STRING_SIZE], SEIS_FILE_VZ[STRING_SIZE];
	extern char  SEIS_FILE_CURL[STRING_SIZE], SEIS_FILE_DIV[STRING_SIZE], SEIS_FILE_P[STRING_SIZE];

        char vxf[STRING_SIZE], vyf[STRING_SIZE],vzf[STRING_SIZE], curlf[STRING_SIZE], divf[STRING_SIZE], pf[STRING_SIZE];
        int nsrc=1;
    
    if(iter>=0){
        sprintf(vxf,"%s.shot%d.it%d",SEIS_FILE_VX,ishot,iter);
        sprintf(vyf,"%s.shot%d.it%d",SEIS_FILE_VY,ishot,iter);
        if(WAVETYPE==2 || WAVETYPE==3) {
            sprintf(vzf,"%s.shot%d.it%d",SEIS_FILE_VZ,ishot,iter);
        }
        sprintf(pf,"%s.shot%d.it%d",SEIS_FILE_P,ishot,iter);
        sprintf(divf,"%s.shot%d.it%d",SEIS_FILE_DIV,ishot,iter);
        sprintf(curlf,"%s.shot%d.it%d",SEIS_FILE_CURL,ishot,iter);
    }
    
    if(iter==-1){
        sprintf(vxf,"%s.shot%d_adjoint_src",SEIS_FILE_VX,ishot);
        sprintf(vyf,"%s.shot%d_adjoint_src",SEIS_FILE_VY,ishot);
        if(WAVETYPE==2 || WAVETYPE==3) {
            sprintf(vzf,"%s.shot%d_adjoint_src",SEIS_FILE_VZ,ishot);
        }
        sprintf(pf,"%s.shot%d_adjoint_src",SEIS_FILE_P,ishot);
        sprintf(divf,"%s.shot%d_adjoint_src",SEIS_FILE_DIV,ishot);
        sprintf(curlf,"%s.shot%d_adjoint_src",SEIS_FILE_CURL,ishot);
    }
    
    
    switch (SEISMO){
        case 1 : /* particle velocities only */
            if (WAVETYPE==1 || WAVETYPE==3) {
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",0,ntr,vxf);
                outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",0,ntr,vyf);
                outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vz) to\n\t %s \n",0,ntr,vzf);
                outseis_glob(fp,fopen(vzf,"w"),2,sectionvz,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            break;
            
        case 2 : /* pressure only */
            if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of pressure to\n\t %s \n",0,ntr,pf);
            outseis_glob(fp,fopen(pf,"w"), 0, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            break;
            
        case 3 : /* curl and div only */
            if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of divergence to\n\t %s \n",0,ntr,divf);
            outseis_glob(fp,fopen(divf,"w"), 0, sectiondiv,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of curl to\n\t %s \n",0,ntr,curlf);
            outseis_glob(fp,fopen(curlf,"w"), 0, sectioncurl,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            break;
            
        case 4 : /* everything */
            if (WAVETYPE==1 || WAVETYPE==3) {
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",0,ntr,vxf);
                outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",0,ntr,vyf);
                outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of pressure to\n\t %s \n",0,ntr,pf);
                outseis_glob(fp,fopen(pf,"w"), 0, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of divergence to\n\t %s \n",0,ntr,divf);
                outseis_glob(fp,fopen(divf,"w"),0, sectiondiv,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of curl to\n\t %s \n",0,ntr,curlf);
                outseis_glob(fp,fopen(curlf,"w"),0, sectioncurl,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vz) to\n\t %s \n",0,ntr,vzf);
                outseis_glob(fp,fopen(vzf,"w"),2,sectionvz,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            break;
            
        case 5 : /* everything except curl and div */
            if (WAVETYPE==1 || WAVETYPE==3) {
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vx) to\n\t %s \n",0,ntr,vxf);
                outseis_glob(fp,fopen(vxf,"w"),1,sectionvx,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vy) to\n\t %s \n",0,ntr,vyf);
                outseis_glob(fp,fopen(vyf,"w"),2,sectionvy,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms of pressure to\n\t %s \n",0,ntr,pf);
                outseis_glob(fp,fopen(pf,"w"), 0, sectionp,recpos,recpos_loc,ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            if(WAVETYPE==2 || WAVETYPE==3) {
                if(VERBOSE==1) fprintf(fp," PE %d is writing %d seismograms (vz) to\n\t %s \n",0,ntr,vzf);
                outseis_glob(fp,fopen(vzf,"w"),2,sectionvz,recpos,recpos_loc, ntr,srcpos,nsrc,ns,SEIS_FORMAT,ishot,1);
            }
            break;
            
    }     
}
