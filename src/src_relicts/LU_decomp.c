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
 *   solve linear matrix equation using LU decomposition
 *   
 *   
 *   Daniel Koehn
 *   last update 9.11.2007
 *
 *  ---------------------------------------------------------------------*/

#include "fd.h"
double LU_decomp(double  **A, double *x, double *b,int n){

int i,j,k,s;
double suma, *y;

y = dvector(1,n);

for (j=1;j<=n;j++){
 for (i=1;i<=j;i++){
   if(i==1){
     A[i][j] = A[i][j];
   }else{
      suma = 0.0;
     for (k=1;k<=i-1;k++){
        suma = suma + A[i][k]*A[k][j];
     }
     A[i][j] = A[i][j] - suma;
   }
 } 
 
 if (j<n){ 
   for (s=1;s<=(n-j);s++){
     i = j + s;
     if (j==1){ 
       A[i][j] = A[i][j]/A[j][j];
     }else{
       suma = 0.0;
       for (k=1;k<=(j-1);k++){
         suma = suma + A[i][k]*A[k][j];
       }
       A[i][j] = (A[i][j] - suma)/A[j][j];
     }
   }
 }
} 

y[1] = b[1];
for (i=2;i<=n;i++){
 suma = 0.0;
 for (j=1;j<=(i-1);j++){
   suma = suma + A[i][j]*y[j];
 }
 y[i] = b[i] - suma;
}

i = n;
j = n;
x[i] = y[i]/A[i][j];
for (s=1;s<=(n-1);s++){
 i = n - s;
 suma = 0.0;
 for (j=i+1;j<=n;j++){
   suma = suma + A[i][j]*x[j];
 }
 x[i] = (y[i] - suma)/A[i][i];
}

free_dvector(y,1,n);

return 0;		
}

