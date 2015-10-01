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


#include "fd.h"

void av_mue(float ** u, float ** uipjp, float ** rho){

	extern int NX, NY, INVMAT1;
	int i, j;
	float u1, u2, u3, u4;
	
	if(INVMAT1==3){
	
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
	       
		       uipjp[j][i]=4.0/((1.0/u[j][i])+
			  	(1.0/u[j][i+1])+(1.0/u[j+1][i])+(1.0/u[j+1][i+1])); 
				
		       if((u[j][i]==0.0)||(u[j][i+1]==0.0)||(u[j+1][i]==0.0)||(u[j+1][i+1]==0.0)){
		           uipjp[j][i]=0.0;
		       }
		      	
	
		}
	}
	
	}
	
	if(INVMAT1==1){
	
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){
	       
	               u1 = rho[j][i] * u[j][i] * u[j][i];
		       u2 = rho[j][i+1] * u[j][i+1] * u[j][i+1];
		       u3 = rho[j+1][i] * u[j+1][i] * u[j+1][i];
		       u4 = rho[j+1][i+1] * u[j+1][i+1] * u[j+1][i+1];
		        
		       uipjp[j][i]=4.0/((1.0/u1)+
			  	(1.0/u2)+(1.0/u3)+(1.0/u4)); 
				
		       if((u1==0.0)||(u2==0.0)||(u3==0.0)||(u4==0.0)){
		           uipjp[j][i]=0.0;
		       }
		      	
	
		}
	}
	
	}
	
}
