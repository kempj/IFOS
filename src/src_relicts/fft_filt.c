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
 *   Write seismograms to disk                                  
 *   last update 19/01/02, T. Bohlen
 *  ----------------------------------------------------------------------*/
#include "fd.h"
#include "segy.h"

void  FFT_filt(float ** data, float freqshift, int ntr, int ns,int method){

	/* declaration of extern variables */
	extern int NDT, MYID, TAPERLENGTH;
	extern char DATA_DIR[STRING_SIZE];
	extern float  TIME, DH, DT, REFREC[4];
	const float xshift=800.0, yshift=800.0;
	
	/* declaration of local variables */
	int i,j, h, nfreq, fa, fb, npad, itr;
	int tracl1 ;
	float xr, yr, x, y, dump, a, times;
	float damping, amp, maximum;
	float *window, *amp1;
	float *trace_real, *trace_complex;
	double npadd;
	
	 /* declaration of local variables */
        int     if1, if2, if3, if4;
        float   f1, f2, f3, f4, df, fi;
	float   F1, F2, F3, F4;
	
	F1=0.0;
	F2=0.0;
	F3=380.0;
	F4=400.0;
	
	npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );  /* ns -> npad for usage in FFT*/
        npadd = (double)npad;
        trace_real = vector(1,npad);
        trace_complex = vector(1,npad);
	
	window = vector(1,(npad/2));
        
	/* calculate TIMEs after application of zero pad*/ 
        times=npad*(TIME/ns); 
          
        /* calculate frequency sample */
	df=1/times;
	
	
	 if (F1==F2)             method = 1;
        else if (F3==F4)         method = 3;
        else                     method = 2;

        f1 = F1 * freqshift;
        f2 = F2 * freqshift;
        f3 = F3 * freqshift;
        f4 = F4 * freqshift;

        if1 = (int)(floor(f1/df))+1;
        if2 = (int)(floor(f2/df))+1;
        if3 = (int)(floor(f3/df))+1;
        if4 = (int)(floor(f4/df))+1;


        /* damping coefficient */
        a = sqrt(-2*log(0.05));

	
	 for(itr=1;itr<=ntr;itr++){
                        for(j=1;j<=ns;j++){
                                trace_real[j] = data[itr][j];
                                trace_complex[j] = 0.0;
                        }
			
         /* apply fft */
         FFT(1,(unsigned long)(log(npadd)/log(2.0)),trace_real,trace_complex);
	
	
        if (method==1)  { /* low-pass */
                /* Gaussian window */  
                for (i=1;i<=(npad/2);i++){
                        fi = (float)(i)*df;
                        window[i] = 1.0;
                        if(i>=if3){
                                dump = a*(fi-f3)/(float)(f4-f3);
                                window[i] = exp(-(1.0/2.0)*dump*dump);
                        }
                }
        } else
        if (method==2)  { /* band-pass */
                /* Gaussian window */  
                for (i=1;i<=(npad/2);i++){
                        fi = (float)(i)*df;
                        window[i] = 1.0;
                        if(i<if2){
                                dump = a*(fi-f2)/(f2-f1);
                                window[i] = exp(-(1.0/2.0)*dump*dump);
                        }
                        if(i>=if3){
                                dump = a*(fi-f3)/(f4-f3);
                                window[i] = exp(-(1.0/2.0)*dump*dump);
                        }
                }
        }

	

	/* Apply taper window */
		for(i=1;i<=npad/2;i++){ 
                   trace_real[i]*=window[i];
		   trace_complex[i]*=window[i];
		   
		   trace_real[npad-i]*=window[i];
		   trace_complex[npad-i]*=window[i];
		}
		
	 /* apply ifft */
         FFT(-1,(unsigned long)(log(npadd)/log(2.0)),trace_real,trace_complex);	
	 
	 /* write output into data matrix */
        for(j=1;j<=ns;j++){
	    
	    data[itr][j] = trace_real[j];
	    
	}
		
        }
	
free_vector(trace_real,1,npad);
free_vector(trace_complex,1,npad);

free_vector(window,1,(npad/2));	

}
