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

void  conv_FD(float * temp_TS, float * temp_TS1, float * temp_conv, int ns){

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
	float *trace_real, *trace_complex, *trace_real1, *trace_complex1;
	float *trace_real_conv, *trace_complex_conv; 
	double npadd;
	
	 /* declaration of local variables */
        int     if1, if2, if3, if4;
    
	npad = (int)(pow(2.0, ceil(log((double)(ns))/log(2.0))+2.0) );  /* ns -> npad for usage in FFT*/
        npadd = (double)npad;
    
	trace_real = vector(1,npad);
        trace_complex = vector(1,npad);

	trace_real1 = vector(1,npad);
        trace_complex1 = vector(1,npad);

	trace_real_conv = vector(1,npad);
        trace_complex_conv = vector(1,npad);

	window = vector(1,(npad/2));
        
	/* calculate TIMEs after application of zero pad*/ 
        times=npad*(TIME/ns); 
          
        for(j=1;j<=ns;j++){
           trace_real[j] = temp_TS[j];
           trace_complex[j] = 0.0;

		   trace_real1[j] = temp_TS1[j];
           trace_complex1[j] = 0.0;

        }
			
         /* apply fft */
         FFT(1,(unsigned long)(log(npadd)/log(2.0)),trace_real,trace_complex);
	 FFT(1,(unsigned long)(log(npadd)/log(2.0)),trace_real1,trace_complex1);

	/* multiplication of complex vectors */
		for(i=1;i<=npad;i++){

           trace_real_conv[i]    = (trace_real[i] * trace_real1[i]) - (trace_complex[i] * trace_complex1[i]);
		   trace_complex_conv[i] = (trace_real[i] * trace_complex1[i]) + (trace_complex[i] * trace_real1[i]);
		   
		}
		
	 /* apply ifft */
         FFT(-1,(unsigned long)(log(npadd)/log(2.0)),trace_real_conv,trace_complex_conv);	
	 
	 /* write output into data matrix */
        for(j=1;j<=ns;j++){
	       temp_conv[j] = trace_real_conv[j];
	    }
		
	
free_vector(trace_real,1,npad);
free_vector(trace_complex,1,npad);

free_vector(trace_real_conv,1,npad);
free_vector(trace_complex_conv,1,npad);

free_vector(trace_real1,1,npad);
free_vector(trace_complex1,1,npad);

free_vector(window,1,(npad/2));	

}
