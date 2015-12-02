/* This function is copied from DENISE Black Edition from D. Koehn  */

#include "fd.h"
void eprecond(float ** W, float ** vx, float ** vy){
    
    extern int NX, NY, IDX, IDY;
    int i, j;
    
    for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
            
            W[j][i]+=((vx[j][i]*vx[j][i])+(vy[j][i]*vy[j][i]));
            
        }
    }

}

void eprecond_SH(float ** W, float ** vz){
    
    extern int NX, NY, IDX, IDY;
    int i, j;
    
    for (i=1;i<=NX;i=i+IDX){
        for (j=1;j<=NY;j=j+IDY){
            
            W[j][i]+=((vz[j][i]*vz[j][i]));
            
        }
    }
    
}