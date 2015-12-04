/*-----------------------------------------------------------------------------------------
 * Copyright (C) 2015  For the list of authors, see file AUTHORS.
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
 *  Apply workflow
 *
 *  Written by Wittkamp Oct 2015
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void apply_workflow(float ** workflow,int workflow_lines,char workflow_header[STRING_SIZE],int *iter,float *FC,int wavetype_start, int * change_wavetype_iter, int * LBFGS_iter_start){
    
    /* local variables */
    int x;
    
    /* extern variables */
    extern int INV_RHO_ITER,INV_VS_ITER,INV_VP_ITER;
    extern int TIME_FILT,MYID;
    extern float PRO;
    extern int WAVETYPE;
    extern float JOINT_INVERSION_PSV_SH_ALPHA_VS;
    extern float JOINT_INVERSION_PSV_SH_ALPHA_RHO;
    extern int EPRECOND;
    extern float EPSILON_WE;
    extern int GRAD_METHOD;
    extern int WORKFLOW_STAGE;
    
    /******************/
    /* Apply Workflow */
    /******************/
    
    /* Print current workflow */
    if(MYID==0){
        printf("\n ---------- Applying Workflow -----------\n");
        printf(" %s ",workflow_header);
        for(x=1;x<=WORKFLOW_MAX_VAR;x++){
            printf("%.2f\t",workflow[WORKFLOW_STAGE][x]);
        }
    }
    
    /* Inversion of material parameter */
    if(workflow[WORKFLOW_STAGE][2]!=-1) {
        if(workflow[WORKFLOW_STAGE][2]==1) {
            if(INV_VS_ITER>*iter) INV_VS_ITER=*iter;
        } else {
            /* detect change and reset LBFGS */
            if(INV_VS_ITER<*iter) *LBFGS_iter_start=*iter;
            INV_VS_ITER=*iter+10;
        }
    }
    
    if(workflow[WORKFLOW_STAGE][3]!=-1){
        if(workflow[WORKFLOW_STAGE][3]==1) {
            if(INV_VP_ITER>*iter) INV_VP_ITER=*iter;
        } else {
            /* detect change and reset LBFGS */
            if(INV_VP_ITER<*iter) *LBFGS_iter_start=*iter;
            INV_VP_ITER=*iter+10;
        }
    }
    
    if(workflow[WORKFLOW_STAGE][4]!=-1){
        if(workflow[WORKFLOW_STAGE][4]==1) {
            if(INV_RHO_ITER>*iter) INV_RHO_ITER=*iter;
        } else {
            /* detect change and reset LBFGS */
            if(INV_RHO_ITER<*iter) *LBFGS_iter_start=*iter;
            INV_RHO_ITER=*iter+10;
        }
    }
    
    PRO=workflow[WORKFLOW_STAGE][5];
    
    /* Frequency filtering  */
    if(TIME_FILT==1) {
        TIME_FILT=workflow[WORKFLOW_STAGE][6];
        if(*FC>workflow[WORKFLOW_STAGE][7]&&(workflow[WORKFLOW_STAGE][6]>0)) {
            if(MYID==0)printf("\n Due to the abort criteriom FC is already higher than specified in workflow\n");
            if(MYID==0)printf(" therefore instead of %.2f HZ FC=%.2f HZ is used\n",workflow[WORKFLOW_STAGE][7],*FC);
        } else {
            if(*FC!=workflow[WORKFLOW_STAGE][7]) *LBFGS_iter_start=*iter;
            *FC=workflow[WORKFLOW_STAGE][7];
        }
    } else {
        if(MYID==0&&(workflow[WORKFLOW_STAGE][6]>0))printf("\n TIME_FILT cannot be activated due to it is not activated in the JSON File \n");
    }
    /* Change of wavetype  */
    if(wavetype_start!=3&&(WAVETYPE!=workflow[WORKFLOW_STAGE][8])){
        if(MYID==0)printf("\n Sorry, change of WAVETYPE with workflow only possible if WAVETYPE==3 in *.json");
        if(MYID==0)printf("\n WAVETYPE will remain unchanged %i",WAVETYPE);
    } else {
        /* detect change and reset some things */
        if(WAVETYPE!=workflow[WORKFLOW_STAGE][8]) {
            *change_wavetype_iter=*iter;
            *LBFGS_iter_start=*iter;
        }
        WAVETYPE=workflow[WORKFLOW_STAGE][8];
    }
    
    /* Joint inversion PSV and SH  */
    JOINT_INVERSION_PSV_SH_ALPHA_VS=workflow[WORKFLOW_STAGE][9];
    JOINT_INVERSION_PSV_SH_ALPHA_RHO=workflow[WORKFLOW_STAGE][10];
    
    /* Approx. Hessian  */
    EPRECOND=workflow[WORKFLOW_STAGE][11];
    EPSILON_WE=workflow[WORKFLOW_STAGE][12];
    
    if(*LBFGS_iter_start==*iter && GRAD_METHOD==2){
        if(MYID==0)printf("\n L-BFGS will be used from iteration %d on.",*LBFGS_iter_start+1);
    }
}