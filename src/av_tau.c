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

/* $Id: av_tau.c,v 2.4 2007/08/21 13:16:19 tbohlen Exp $*/

#include "fd.h"

void av_tau(float **taus, float **tausipjp){

	extern int NX, NY;
	int i, j;
	for (j=1;j<=NY;j++){
		for (i=1;i<=NX;i++){

		       tausipjp[j][i] = 0.25*(taus[j][i]+taus[j][i+1]+taus[j+1][i]+taus[j+1][i+1]);
		}
	}
}
