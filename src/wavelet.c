/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2016  For the list of authors, see file AUTHORS.
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
*   Calculating source signal at different source positions with different
*   time-shift, centre frequency and amplitude (as specified in SOURCE_FILE).
*   Source signals are written to array signals 
*
*   last update 19/01/02, T. Bohlen
*  ----------------------------------------------------------------------*/

#include "fd.h"


float ** wavelet(float ** srcpos_loc, int nsrc, int ishot,int SH){
    
    
    /* extern variables */
    extern int QUELLART,QUELLART_SH, NT, MYID,VERBOSE;
    extern float  DT;
    extern char SIGNAL_FILE[STRING_SIZE];
    extern char SIGNAL_FILE_SH[STRING_SIZE];
    extern FILE *FP;
    
    /*local variables */
    int nts, nt, k;
    float *psource=NULL, tshift, amp=0.0, a, fc, tau, t, ts, ag;
    float ** signals;
    
    
    
    /* ---------------------------- */
    /*   P SV Source                */
    /* ---------------------------- */
    if (SH==0) {
        if (QUELLART==3) psource=rd_sour(&nts,fopen(SIGNAL_FILE,"r"));
        if (QUELLART==7){
            psource=vector(1,NT);
            inseis_source_wavelet(psource,NT,ishot);}
        
        signals=fmatrix(1,nsrc,1,NT);
        
        for (nt=1;nt<=NT;nt++){
            t=(float)nt*DT;
            
            for (k=1;k<=nsrc;k++) {
                tshift=srcpos_loc[4][k];
                fc=srcpos_loc[5][k];
                a=srcpos_loc[6][k];
                ts=1.0/fc;
                ag = PI*PI*fc*fc;
                
                switch (QUELLART){
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
                }
                
                
                signals[k][nt]=amp*a;
            }
        }
    
        if(VERBOSE){
            fprintf(FP," Message from function wavelet written by PE %d \n",MYID);
            fprintf(FP," %d source positions located in subdomain of PE %d \n",nsrc,MYID);
            fprintf(FP," have been assigned with a source signal. \n");
        }
        
        if (QUELLART==3 || QUELLART==7) free_vector(psource,1,NT);
        
        return signals;	
        
    }
    
    /* ---------------------------- */
    /*   SH Source                  */
    /* ---------------------------- */

    if (SH==1) {
        if (QUELLART_SH==3) psource=rd_sour(&nts,fopen(SIGNAL_FILE_SH,"r"));
        if (QUELLART_SH==7){
            psource=vector(1,NT);
            inseis_source_wavelet(psource,NT,ishot);}
        
        signals=fmatrix(1,nsrc,1,NT);
        
        for (nt=1;nt<=NT;nt++){
            t=(float)nt*DT;
            
            for (k=1;k<=nsrc;k++) {
                tshift=srcpos_loc[4][k];
                fc=srcpos_loc[5][k];
                a=srcpos_loc[6][k];
                ts=1.0/fc;
                ag = PI*PI*fc*fc;
                
                switch (QUELLART_SH){
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
                }
                
                
                signals[k][nt]=amp*a;
            }
        }
        if(VERBOSE){
            fprintf(FP,"\n Message from function wavelet SH written by PE %d \n",MYID);
            fprintf(FP," %d source positions located in subdomain of PE %d \n",nsrc,MYID);
            fprintf(FP," have been assigned with a source signal. \n");
        }
        
        if (QUELLART_SH==3 || QUELLART_SH==7) free_vector(psource,1,NT);
        
        return signals;	
        
    }
}
