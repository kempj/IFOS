/*------------------------------------------------------------------------
 *  read_parameters.c
 *
 *  This file was copied from the DENISE code and adjusted for its use
 *  in this test programme.
 *  For copyright information see the DENISE code.            
 *  ---------------------------------------------------------------------*/


#include "fd.h"

void read_parameters(FILE *fp_in){


/* declaration of extern variables */
extern int ntr;
extern char FILE_OBS[STRING_SIZE], FILE_SYN[STRING_SIZE], FILE_CONV[STRING_SIZE], FILE_STF[STRING_SIZE], FILE_ADD_CONV[STRING_SIZE], FILE_ADD[STRING_SIZE];
extern char stfinv_par[STRING_SIZE]; 

/* definition of local variables */
char s[74];
int  c=0, lineno=0;

while ((c=getc(fp_in)) != EOF){
      if ((c=='\n') && (getc(fp_in)!='#')){     
	 lineno++;
	 switch (lineno){
	 case 1 :
	    fscanf(fp_in,"%s =%s",s,FILE_OBS);
	    break;
	 case 2 :
	    fscanf(fp_in,"%s =%s",s,FILE_SYN);
	    break;     		
	 case 3 :
	    fscanf(fp_in,"%s =%i",s,&ntr);
	    break;
	 case 4 :
	    fscanf(fp_in,"%s =%s",s,FILE_CONV);
	    break;
	 case 5 :
	    fscanf(fp_in,"%s =%s",s,FILE_STF);
	    break;
	 case 6 :
	    fscanf(fp_in,"%s =%s",s,stfinv_par);
	    break;
	 case 7 :
	    fscanf(fp_in,"%s =%s",s,FILE_ADD);
	    break;
	 case 8 :
	    fscanf(fp_in,"%s =%s",s,FILE_ADD_CONV);  
	 }
	}
}

fclose(fp_in);

}



