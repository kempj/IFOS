/*------------------------------------------------------------------------
 *  fd.h - include file for FD programs
 *
 *  This file was copied from the DENISE code and adjusted for its use
 *  in this test programme.
 *  For copyright information see the DENISE code.            
 *  ---------------------------------------------------------------------*/

/* files to include */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>

#define iround(x) ((int)(floor)(x+0.5))
#define min(x,y) ((x<y)?x:y)    
#define max(x,y) ((x<y)?y:x)
#define fsign(x) ((x<0.0)?(-1):1)    

#define PI (3.141592653589793)
#define NPAR 50
#define STRING_SIZE 74
#define REQUEST_COUNT 4


/* declaration of functions */
void read_parameters(FILE *fp_in);



