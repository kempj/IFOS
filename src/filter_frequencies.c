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

/* 
   Reading filter frequencies for the time domain filter from FREQ_FILE.
   
   T. Metz 02,2014 
*/

#include "fd.h"

float *filter_frequencies(int *nfrq){

extern char FREQ_FILE[STRING_SIZE];
extern FILE *FP;

FILE *freqf;
char cline[256];
int IT=1;
float *FC;


/*----------------------------------- open FREQ_FILE ------------------------------*/
fprintf(FP,"\n\n ------------------- FREQUENCY FILTERING --------------------------\n");

freqf=fopen(FREQ_FILE,"r");

if (freqf==NULL) err(" Freqency file could not be opened !");
*nfrq=0;

fprintf(FP,"\n Reading frequencies from %s \n",FREQ_FILE);

/*----------------------Read Number of requencies defined in FREQ_FILE-------------*/
fscanf(freqf,"%i",nfrq);

/*----------------------------alocate frequency array------------------------------*/
FC=vector(1,*nfrq);
//FC = malloc((*nfrq+1) * sizeof(float));

rewind (freqf);
fprintf(FP," Number of freqencies specified: %d \n",*nfrq);

/*----------------------Read frequencies from FREQ_FILE----------------------------*/
for (IT=1;IT<=(*nfrq);IT++){
      fgets(cline,255,freqf);
      fscanf(freqf,"%f",&FC[IT]);
      //Read only numbers from FREQ_FILE 
      if(FC[IT]==0){
      IT=IT-1;}else{
      fprintf(FP,"\n %d. %.2f Hz",IT,FC[IT]);
      }
}

fclose(freqf);

return(FC);
}

