#include "fd.h"
#include "globvar.h"
#include "cseife.h"
#include "stfinv/stfinv.h"


void init_vars
{
    int ns, nseismograms=0, nt, nd, fdo3, j, i, iter, h, infoout, SHOTINC,  hin, hin1, do_stf=0;
    int NTDTINV, nxny, nxnyi, imat, imat1, imat2, IDXI, IDYI, hi, NTST, NTSTI;
    int lsnap, nsnap=0, lsamp=0, buffsize,  swstestshot, snapseis, snapseis1;
    int ntr=0, ntr_loc=0, ntr_glob=0, nsrc=0, nsrc_loc=0, nsrc_glob=0, ishot, irec, nshots=0, nshots1, Lcount, itest,
        itestshot;
    float muss, lamss;
    float memdyn, memmodel, memseismograms, membuffer, memtotal, eps_scale;
    float fac1, fac2;
    float opteps_vp, opteps_vs, opteps_rho, Vp_avg, C_vp, Vs_avg, C_vs, rho_avg, C_rho;
    float memfwt, memfwt1, memfwtdata;
    char *buff_addr, ext[10], *fileinp;
    char jac[225];
    double time1, time2, time3, time4, time5, time6, time7, time8,
           time_av_v_update=0.0, time_av_s_update=0.0, time_av_v_exchange=0.0,
           time_av_s_exchange=0.0, time_av_timestep=0.0;
    float L2, L2sum, L2_all_shots, L2sum_all_shots, *L2t, alphanom, alphadenom;
    int sum_killed_traces=0, sum_killed_traces_testshots=0, killed_traces=0, killed_traces_testshots=0;
    int *ptr_killed_traces=&killed_traces, *ptr_killed_traces_testshots=&killed_traces_testshots;
    float energy, energy_sum, energy_all_shots, energy_sum_all_shots = 0.0;
    float energy_SH, energy_sum_SH, energy_all_shots_SH, energy_sum_all_shots_SH;
    float L2_SH, L2sum_SH, L2_all_shots_SH, L2sum_all_shots_SH;
    // Pointer for dynamic wavefields:
    float ** psxx, **  psxy, **  psyy, **  psxz, **  psyz, **psp,
          ** ux, ** uy, ** uxy, ** uyx, ** u, ** Vp0, ** uttx, ** utty, ** Vs0, ** Rho0;
    float ** pvx, **  pvy, **  pvz, **waveconv, **waveconv_lam, **waveconv_mu,
          ** waveconv_rho, **waveconv_rho_s, **waveconv_u, **waveconvtmp, **wcpart,
          ** wavejac,**waveconv_rho_s_z,**waveconv_u_z,**waveconv_rho_z;
    float ** waveconv_shot, **waveconv_u_shot, **waveconv_rho_shot, **waveconv_u_shot_z, **waveconv_rho_shot_z;
    float ** pvxp1, **  pvyp1, **  pvzp1, **  pvxm1, **  pvym1, **  pvzm1;
    float ** gradg, ** gradp,** gradg_rho, ** gradp_rho, ** gradg_u, ** gradp_u, ** gradp_u_z,** gradp_rho_z;
    float ** prho,**  prhonp1, **prip=NULL, **prjp=NULL, **pripnp1=NULL, **prjpnp1=NULL,
                                 ** ppi, **  pu, **  punp1, **  puipjp, **  ppinp1;
    float ** vpmat;
    float ***forward_prop_x, ***forward_prop_y, ***forward_prop_rho_x, ***forward_prop_u, ***forward_prop_rho_y,
          ***forward_prop_p;
    float ***forward_prop_z_xz,***forward_prop_z_yz,***forward_prop_rho_z,**waveconv_mu_z;
    float ** uxz, ** uyz;
    float ** sectionvx=NULL, ** sectionvy=NULL, ** sectionvz=NULL, ** sectionp=NULL, ** sectionpnp1=NULL,
             ** sectioncurl=NULL, ** sectiondiv=NULL, ** sectionvxdata=NULL, ** sectionvydata=NULL,
                ** sectionvzdata=NULL, ** sectionvxdiff=NULL, ** sectionvzdiff=NULL, ** sectionvxdiffold=NULL,
                   ** sectionvydiffold=NULL, ** sectionvzdiffold=NULL,** sectionpdata=NULL, ** sectionpdiff=NULL,
                      ** sectionpdiffold=NULL, ** sectionvydiff=NULL, ** sectionpn=NULL, ** sectionread=NULL,
                         ** sectionvy_conv=NULL, ** sectionvy_obs=NULL, ** sectionvx_conv=NULL,** sectionvx_obs=NULL,
                            ** sectionvz_conv=NULL,** sectionvz_obs=NULL, ** sectionp_conv=NULL,** sectionp_obs=NULL,
                               * source_time_function=NULL;
    float ** absorb_coeff, ** taper_coeff, * epst1, * epst2,  * epst3, * picked_times;
    float ** srcpos=NULL, **srcpos_loc=NULL, ** srcpos1=NULL, **srcpos_loc_back=NULL, ** signals=NULL,** signals_SH=NULL,
             *hc=NULL;
    int   ** recpos=NULL, ** recpos_loc=NULL;
    /*int   ** tracekill=NULL, TRKILL, DTRKILL;*/
    int * DTINV_help;
    float ** bufferlef_to_rig,  ** bufferrig_to_lef, ** buffertop_to_bot, ** bufferbot_to_top;
    /* PML variables */
    float * d_x, * K_x, * alpha_prime_x, * a_x, * b_x, * d_x_half, * K_x_half,
          * alpha_prime_x_half, * a_x_half, * b_x_half, * d_y, * K_y,
          * * alpha_prime_y, * a_y, * b_y, * d_y_half, * K_y_half,
          * * alpha_prime_y_half, * a_y_half, * b_y_half;
    float ** psi_sxx_x, ** psi_syy_y, ** psi_sxy_y, ** psi_sxy_x, ** psi_vxx,
          ** psi_vyy, ** psi_vxy, ** psi_vyx, ** psi_vxxs;
    float ** psi_sxz_x, ** psi_syz_y, ** psi_vzx, ** psi_vzy;
    /* Variables for viscoelastic modeling */
    float **ptaus=NULL, **ptaup=NULL, *etaip=NULL, *etajm=NULL, *peta=NULL,
            **ptausipjp=NULL, **fipjp=NULL, ***dip=NULL, *bip=NULL, *bjm=NULL;
    float *cip=NULL, *cjm=NULL, ***d=NULL, ***e=NULL, ***pr=NULL, ***pp=NULL,
           ***pq=NULL, **f=NULL, **g=NULL;
    float ***pt=NULL, ***po=NULL; // SH Simulation
    /* Variables for step length calculation */
    int step1, step2, step3=0, itests, iteste, stepmax, countstep;
    float scalefac;
    int RECINC, ntr1;
    int SOURCE_SHAPE_OLD=0;
    int SOURCE_SHAPE_OLD_SH=0;
    /* Variables for conjungate gradient */
    int PCG_iter_start=1;
    /* Variables for L-BFGS */
    int LBFGS_NPAR=3;
    int LBFGS_iter_start=1;
    float **s_LBFGS,**y_LBFGS, *rho_LBFGS;
    int l=0;
    int m=0;
    /* Check wolfe */
    int steplength_search=0;
    int FWI_run=1;
    int gradient_optimization=1;
    float alpha_SL_min=0, alpha_SL_max=0, alpha_SL=1.0;
    float alpha_SL_old;
    float ** waveconv_old,** waveconv_u_old,** waveconv_rho_old;
    float ** waveconv_up,** waveconv_u_up,** waveconv_rho_up;
    float L2_SL_old=0, L2_SL_new=0;
    float c1_SL=1e-4, c2_SL=0.9;
    int wolfe_status;
    int wolfe_sum_FWI=0;
    int wolfe_found_lower_L2=0;
    float alpha_SL_FS;
    float L2_SL_FS;
    int use_wolfe_failsafe=0;
    int wolfe_SLS_failed=0;
    /* Variables for energy weighted gradient */
    float ** Ws, **Wr, **We;
    float ** Ws_SH, **Wr_SH, **We_SH;
    float ** We_sum,** We_sum_SH;
    float We_sum_max1;
    float We_max_SH,We_max;
    int * recswitch=NULL;
    float ** fulldata=NULL, ** fulldata_vx=NULL, ** fulldata_vy=NULL, ** fulldata_vz=NULL, ** fulldata_p=NULL,
             ** fulldata_curl=NULL, ** fulldata_div=NULL;
    /*vector for abort criterion*/
    float * L2_hist=NULL;
    /* help variable for MIN_ITER */
    int min_iter_help=0;
    float ** workflow=NULL;
    int workflow_lines;
    char workflow_header[STRING_SIZE];
    int change_wavetype_iter=-10; /* Have to be inialized negative */
    int wavetype_start; /* We need this due to MPI Comm */
    int buf1=0, buf2=0;
    WORKFLOW_STAGE=1;
    /* variable for time domain filtering */
    float F_LOW_PASS;
    float *F_LOW_PASS_EXT=NULL;
    int nfrq=0;
    int FREQ_NR=1;
    float JOINT_EQUAL_PSV=0.0, JOINT_EQUAL_SH=0.0;
    float JOINT_EQUAL_PSV_all=0.0, JOINT_EQUAL_SH_all=0.0;
    int JOINT_EQUAL_new_max=1;
    FILE *fprec, *FPL2;
    FILE *FPL2_JOINT;
    char L2_joint_log[STRING_SIZE];
    /* General parameters */
    int nt_out;
    MPI_Request *req_send, *req_rec;
    MPI_Status  *send_statuses, *rec_statuses;
    /* Initialize MPI environment */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&NP);
    MPI_Comm_rank(MPI_COMM_WORLD,&MYID);
    setvbuf(stdout, NULL, _IONBF, 0);
    if (MYID == 0) {
        time1=MPI_Wtime();
        clock();
    }
    /* print program name, version etc to stdout*/
    if (MYID == 0) info(stdout);
    /* read parameters from parameter-file (stdin) */
    fileinp=argv[1];
    FP=fopen(fileinp,"r");
    if(FP==NULL) {
        if (MYID == 0) {
            printf("\n==================================================================\n");
            printf(" Cannot open IFOS input file %s \n",fileinp);
            printf("\n==================================================================\n\n");
            declare_error(" --- ");
        }
    }
    /* read json formatted input file */
    read_par_json(stdout,fileinp);
    exchange_par();
    wavetype_start=WAVETYPE;
    if (MYID == 0) note(stdout);
    /* open log-file (each PE is using different file) */
    /*  fp=stdout; */
    sprintf(ext,".%i",MYID);
    strcat(LOG_FILE,ext);
    /* If Verbose==0, no PE will write a log file */
    if(!VERBOSE) sprintf(LOG_FILE,"/dev/null");
    if ((MYID==0)) FP=stdout;
    else {
        FP=fopen(LOG_FILE,"w");
    }
    fprintf(FP," This is the log-file generated by PE %d \n\n",MYID);
    /* domain decomposition */
    initproc();
    NT=iround(TIME/DT);       /* number of timesteps */
    /*ns=iround(NT/NDT);*/           /* number of samples per trace */
    ns=NT;  /* in a FWI one has to keep all samples of the forward modeled data
               at the receiver positions to calculate the adjoint sources and to do
               the backpropagation; look at function saveseis_glob.c to see that every
               NDT sample for the forward modeled wavefield is written to su files*/
    /* output of parameters to log-file or stdout */
    if (MYID==0) write_par(FP);
    /* NXG, NYG denote size of the entire (global) grid */
    NXG=NX;
    NYG=NY;
    /* In the following, NX and NY denote size of the local grid ! */
    NX = IENDX;
    NY = IENDY;
    /* Reading source positions from SOURCE_FILE */
    srcpos=sources(&nsrc);
    nsrc_glob=nsrc;
    ishot=0;
    if (SEISMO) {
        recpos=receiver(&ntr, srcpos, ishot);
        recswitch = ivector(1,ntr);
        recpos_loc = splitrec(recpos,&ntr_loc, ntr, recswitch);
        ntr_glob=ntr;
        ntr=ntr_loc;
    }
    /* memory allocation for abort criterion*/
    L2_hist = vector(1,1000);
    if(INV_STF)
        fulldata = matrix(1,ntr_glob,1,NT);
    /* estimate memory requirement of the variables in megabytes*/
    switch (SEISMO) {
    case 1 : /* particle velocities only */
        nseismograms=2;
        break;
    case 2 : /* pressure only */
        nseismograms=1;
        break;
    case 3 : /* curl and div only */
        nseismograms=2;
        break;
    case 4 : /* everything */
        nseismograms=5;
        break;
    case 5 : /* everything except curl and div */
        nseismograms=3;
        break;
    }
    /* use only every DTINV time sample for the inversion */
    /*DTINV=15;*/
    DTINV_help=ivector(1,NT);
    NTDTINV=ceil((float)NT/(float)DTINV);       /* round towards next higher integer value */
    /* save every IDXI and IDYI spatial point during the forward modelling */
    IDXI=1;
    IDYI=1;
    /*allocate memory for dynamic, static and buffer arrays */
    fac1=(NX+FDORDER)*(NY+FDORDER);
    fac2=sizeof(float)*pow(2.0,-20.0);
    nd = FDORDER/2 + 1;
    // decide how much space for exchange is needed
    switch (WAVETYPE) {
    case 1:
        fdo3 = 2*nd;
        break;
    case 2:
        fdo3 = 1*nd;
        break;
    case 3:
        fdo3 = 3*nd;
        break;
    default:
        fdo3 = 2*nd;
        break;
    }
    if (L) {
        memdyn=(5.0+3.0*(float)L)*fac1*fac2;
        memmodel=(12.0+3.0*(float)L)*fac1*fac2;
    } else {
        memdyn=5.0*fac1*fac2;
        memmodel=6.0*fac1*fac2;
    }
    memseismograms=nseismograms*ntr*ns*fac2;
    memfwt=5.0*((NX/IDXI)+FDORDER)*((NY/IDYI)+FDORDER)*NTDTINV*fac2;
    memfwt1=20.0*NX*NY*fac2;
    memfwtdata=6.0*ntr*ns*fac2;
    membuffer=2.0*fdo3*(NY+NX)*fac2;
    buffsize=2.0*2.0*fdo3*(NX+NY)*sizeof(MPI_FLOAT);
    memtotal=memdyn+memmodel+memseismograms+memfwt+memfwt1+memfwtdata+membuffer+(buffsize*pow(2.0,-20.0));
    if (MYID==0 && WAVETYPE == 1) {
        fprintf(FP,"\n **Message from main (printed by PE %d):\n",MYID);
        fprintf(FP," Size of local grids: NX=%d \t NY=%d\n",NX,NY);
        fprintf(FP," Each process is now trying to allocate memory for:\n");
        fprintf(FP," Dynamic variables: \t\t %6.2f MB\n", memdyn);
        fprintf(FP," Static variables: \t\t %6.2f MB\n", memmodel);
        fprintf(FP," Seismograms: \t\t\t %6.2f MB\n", memseismograms);
        fprintf(FP," Buffer arrays for grid exchange:%6.2f MB\n", membuffer);
        fprintf(FP," Network Buffer for MPI_Bsend: \t %6.2f MB\n", buffsize*pow(2.0,-20.0));
        fprintf(FP," ------------------------------------------------ \n");
        fprintf(FP," Total memory required: \t %6.2f MB.\n\n", memtotal);
    }
    /* allocate buffer for buffering messages */
    buff_addr=malloc(buffsize);
    if (!buff_addr) declare_error("allocation failure for buffer for MPI_Bsend !");
    MPI_Buffer_attach(buff_addr,buffsize);
    /* allocation for request and status arrays */
    req_send=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
    req_rec=(MPI_Request *)malloc(REQUEST_COUNT*sizeof(MPI_Request));
    send_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));
    rec_statuses=(MPI_Status *)malloc(REQUEST_COUNT*sizeof(MPI_Status));
    /* memory allocation for dynamic (wavefield) arrays */
    if(!ACOUSTIC) {
        switch (WAVETYPE) {
        case 1: // P and SV Waves
            psxx =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            psxy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            psyy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            break;
        case 2: // SH Waves
            psxz =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            psyz =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            break;
        case 3: // P, SH and SV Waves
            psxx =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            psxy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            psyy =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            psxz =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            psyz =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            break;
        }
    } else {
        psp  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    }
    if(GRAD_METHOD==2) {
        /* Allocate memory for L-BFGS */
        if(WAVETYPE==2) LBFGS_NPAR=2;
        s_LBFGS=fmatrix(1,N_LBFGS,1,LBFGS_NPAR*NX*NY);
        y_LBFGS=fmatrix(1,N_LBFGS,1,LBFGS_NPAR*NX*NY);
        rho_LBFGS=vector(1,N_LBFGS);
        for(l=1; l<=N_LBFGS; l++) {
            for(m=1; m<=LBFGS_NPAR*NX*NY; m++) {
                s_LBFGS[l][m]=0.0;
                y_LBFGS[l][m]=0.0;
            }
            rho_LBFGS[l]=0.0;
        }
    }
    if(!ACOUSTIC) {
        if(WAVETYPE==1||WAVETYPE==3) {
            ux   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            uy   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            uxy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            uyx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            uttx   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            utty   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        }
        if(WAVETYPE==2||WAVETYPE==3) {
            uxz   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            uyz   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        }
    } else {
        u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    }
    switch (WAVETYPE) {
    case 1: // P and SV Waves
        pvx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvxp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvyp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvxm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvym1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        break;
    case 2: // SH Waves
        pvz  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvzp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvzm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        break;
    case 3: // P and SV Waves
        pvx  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvy  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvxp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvyp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvxm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvym1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvz  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvzp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        pvzm1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        break;
    }
    Vp0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    if(!ACOUSTIC)
        Vs0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    Rho0  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    /* memory allocation for static (model) arrays */
    prho =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    prhonp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    prip =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    prjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    pripnp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    prjpnp1 =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    ppi  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    ppinp1  =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    if(!ACOUSTIC) {
        pu   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        punp1   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        puipjp   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    }
    vpmat   =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
    if((EPRECOND==1)||(EPRECOND==3)) {
        if(WAVETYPE==1 || WAVETYPE==3) {
            We_sum = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            Ws = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of the source wavefield */
            Wr = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of the receiver wavefield */
            We = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of source and receiver wavefield */
        }
        if(WAVETYPE==2 || WAVETYPE==3) {
            We_sum_SH = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            Ws_SH = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of the source wavefield */
            Wr_SH = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of the receiver wavefield */
            We_SH = matrix(-nd+1,NY+nd,-nd+1,NX+nd); /* total energy of source and receiver wavefield */
        }
    }
    if (L) {
        /* dynamic (wavefield) arrays for viscoelastic modeling */
        pr = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        pp = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        pq = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        /* memory allocation for static arrays for viscoelastic modeling */
        dip = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        d =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        e =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        ptaus =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        ptausipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        if(WAVETYPE==2 || WAVETYPE==3) {
            pt = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
            po = f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,L);
        }
        ptaup =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        fipjp =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        f =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        g =  matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        peta =  vector(1,L);
        etaip =  vector(1,L);
        etajm =  vector(1,L);
        bip =  vector(1,L);
        bjm =  vector(1,L);
        cip =  vector(1,L);
        cjm =  vector(1,L);
    }
    NTST=20;
    nxnyi=(NX/IDXI)*(NY/IDYI);
    /* Parameters for step length calculations */
    stepmax = STEPMAX; /* number of maximum misfit calculations/steplength 2/3*/
    scalefac = SCALEFAC; /* scale factor for the step length */
    if(FORWARD_ONLY==0) {
        waveconv = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconv_lam = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconv_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconvtmp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        wcpart = matrix(1,3,1,3);
        wavejac = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        if(!ACOUSTIC) {
            forward_prop_x =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            forward_prop_y =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        } else {
            forward_prop_p =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        }
        gradg = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        gradp = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        if(WAVETYPE==1 || WAVETYPE==3) {
            forward_prop_rho_x =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            forward_prop_rho_y =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
        }
        if(WAVETYPE==2 || WAVETYPE==3) {
            forward_prop_rho_z =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            forward_prop_z_xz =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            forward_prop_z_yz =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            waveconv_rho_shot_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_u_shot_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_mu_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_rho_s_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_u_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_rho_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            gradp_u_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            gradp_rho_z = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        }
        gradg_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        gradp_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconv_rho = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconv_rho_s = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        waveconv_rho_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        if(WOLFE_CONDITION) {
            c1_SL=WOLFE_C1_SL;
            c2_SL=WOLFE_C2_SL;
            waveconv_old= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            if(!ACOUSTIC) waveconv_u_old= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_rho_old= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_up= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            if(!ACOUSTIC) waveconv_u_up= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_rho_up= matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        }
        if(!ACOUSTIC) {
            forward_prop_u =  f3tensor(-nd+1,NY+nd,-nd+1,NX+nd,1,NT/DTINV);
            gradg_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            gradp_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_u = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_mu = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
            waveconv_u_shot = matrix(-nd+1,NY+nd,-nd+1,NX+nd);
        }
    }
    /* Allocate memory for boundary */
    if(FW>0) {
        d_x = vector(1,2*FW);
        K_x = vector(1,2*FW);
        alpha_prime_x = vector(1,2*FW);
        a_x = vector(1,2*FW);
        b_x = vector(1,2*FW);
        d_x_half = vector(1,2*FW);
        K_x_half = vector(1,2*FW);
        alpha_prime_x_half = vector(1,2*FW);
        a_x_half = vector(1,2*FW);
        b_x_half = vector(1,2*FW);
        d_y = vector(1,2*FW);
        K_y = vector(1,2*FW);
        alpha_prime_y = vector(1,2*FW);
        a_y = vector(1,2*FW);
        b_y = vector(1,2*FW);
        d_y_half = vector(1,2*FW);
        K_y_half = vector(1,2*FW);
        alpha_prime_y_half = vector(1,2*FW);
        a_y_half = vector(1,2*FW);
        b_y_half = vector(1,2*FW);
        if (WAVETYPE==1||WAVETYPE==3) {
            psi_sxx_x =  matrix(1,NY,1,2*FW);
            psi_syy_y =  matrix(1,2*FW,1,NX);
            psi_sxy_y =  matrix(1,2*FW,1,NX);
            psi_sxy_x =  matrix(1,NY,1,2*FW);
            psi_vxx   =  matrix(1,NY,1,2*FW);
            psi_vxxs  =  matrix(1,NY,1,2*FW);
            psi_vyy   =  matrix(1,2*FW,1,NX);
            psi_vxy   =  matrix(1,2*FW,1,NX);
            psi_vyx   =  matrix(1,NY,1,2*FW);
        }
        if(WAVETYPE==2||WAVETYPE == 3 ) {
            psi_sxz_x =  matrix(1,NY,1,2*FW);
            psi_syz_y =  matrix(1,2*FW,1,NX);
            psi_vzx   =  matrix(1,NY,1,2*FW);
            psi_vzy   =  matrix(1,2*FW,1,NX);
        }
    }
    taper_coeff=  matrix(1,NY,1,NX);
    /* memory allocation for buffer arrays in which the wavefield
       information which is exchanged between neighbouring PEs is stored */
    bufferlef_to_rig = matrix(1,NY,1,fdo3);
    bufferrig_to_lef = matrix(1,NY,1,fdo3);
    buffertop_to_bot = matrix(1,NX,1,fdo3);
    bufferbot_to_top = matrix(1,NX,1,fdo3);
    /* Allocate memory to save full seismograms */
    switch (SEISMO) {
    case 1 : /* particle velocities only */
        switch (WAVETYPE) {
        case 1:
            fulldata_vx = matrix(1,ntr_glob,1,NT);
            fulldata_vy = matrix(1,ntr_glob,1,NT);
            break;
        case 2:
            fulldata_vz = matrix(1,ntr_glob,1,NT);
            break;
        case 3:
            fulldata_vx = matrix(1,ntr_glob,1,NT);
            fulldata_vy = matrix(1,ntr_glob,1,NT);
            fulldata_vz = matrix(1,ntr_glob,1,NT);
            break;
        }
        break;
    case 2 : /* pressure only */
        fulldata_p = matrix(1,ntr_glob,1,NT);
        break;
    case 3 : /* curl and div only */
        fulldata_div = matrix(1,ntr_glob,1,NT);
        fulldata_curl = matrix(1,ntr_glob,1,NT);
        break;
    case 4 : /* everything */
        switch (WAVETYPE) {
        case 1:
            fulldata_vx = matrix(1,ntr_glob,1,NT);
            fulldata_vy = matrix(1,ntr_glob,1,NT);
            break;
        case 2:
            fulldata_vz = matrix(1,ntr_glob,1,NT);
            break;
        case 3:
            fulldata_vx = matrix(1,ntr_glob,1,NT);
            fulldata_vy = matrix(1,ntr_glob,1,NT);
            fulldata_vz = matrix(1,ntr_glob,1,NT);
            break;
        }
        fulldata_p = matrix(1,ntr_glob,1,NT);
        fulldata_div = matrix(1,ntr_glob,1,NT);
        fulldata_curl = matrix(1,ntr_glob,1,NT);
        break;
    case 5 : /* everything except curl and div*/
        switch (WAVETYPE) {
        case 1:
            fulldata_vx = matrix(1,ntr_glob,1,NT);
            fulldata_vy = matrix(1,ntr_glob,1,NT);
            break;
        case 2:
            fulldata_vz = matrix(1,ntr_glob,1,NT);
            break;
        case 3:
            fulldata_vx = matrix(1,ntr_glob,1,NT);
            fulldata_vy = matrix(1,ntr_glob,1,NT);
            fulldata_vz = matrix(1,ntr_glob,1,NT);
            break;
        }
        fulldata_p = matrix(1,ntr_glob,1,NT);
        break;
    }
    if (ntr>0) {
        alloc_sections(ntr,ns,&sectionvx,&sectionvy,&sectionvz,&sectionp,&sectionpnp1,&sectionpn,&sectioncurl,&sectiondiv,
                       &sectionpdata,&sectionpdiff,&sectionpdiffold,&sectionvxdata,&sectionvxdiff,&sectionvxdiffold,&sectionvydata,
                       &sectionvydiff,&sectionvydiffold,&sectionvzdata,&sectionvzdiff,&sectionvzdiffold);
    }
    /* Memory for seismic data */
    sectionread=matrix(1,ntr_glob,1,ns);
    /* Memory for inversion for source time function */
    if((INV_STF==1)||(TIME_FILT==1) || (TIME_FILT==2)) {
        sectionp_conv=matrix(1,ntr_glob,1,NT);
        sectionp_obs=matrix(1,ntr_glob,1,NT);
        source_time_function = vector(1,NT);
        switch (WAVETYPE) {
        case 1:
            sectionvy_conv=matrix(1,ntr_glob,1,NT);
            sectionvy_obs=matrix(1,ntr_glob,1,NT);
            sectionvx_conv=matrix(1,ntr_glob,1,NT);
            sectionvx_obs=matrix(1,ntr_glob,1,NT);
            break;
        case 2:
            sectionvz_conv=matrix(1,ntr_glob,1,NT);
            sectionvz_obs=matrix(1,ntr_glob,1,NT);
            break;
        case 3:
            sectionvy_conv=matrix(1,ntr_glob,1,NT);
            sectionvy_obs=matrix(1,ntr_glob,1,NT);
            sectionvx_conv=matrix(1,ntr_glob,1,NT);
            sectionvx_obs=matrix(1,ntr_glob,1,NT);
            sectionvz_conv=matrix(1,ntr_glob,1,NT);
            sectionvz_obs=matrix(1,ntr_glob,1,NT);
            break;
        }
    }
    /* memory for source position definition */
    srcpos1=fmatrix(1,8,1,1);
    /* memory of L2 norm */
    L2t = vector(1,4);
    epst1 = vector(1,3);
    epst2 = vector(1,3);
    epst3 = vector(1,3);
    picked_times = vector(1,ntr);
    fprintf(FP," ... memory allocation for PE %d was successfull.\n\n", MYID);
    /* Holberg coefficients for FD operators*/
    hc = holbergcoeff();
    MPI_Barrier(MPI_COMM_WORLD);
    if(FORWARD_ONLY==0&&USE_WORKFLOW) {
        read_workflow(FILE_WORKFLOW,&workflow, &workflow_lines,workflow_header);
    }
    /* create model grids */
    if(L) {
        if(!ACOUSTIC) {
            if (READMOD) {
                readmod(prho,ppi,pu,ptaus,ptaup,peta);
            } else {
                model(prho,ppi,pu,ptaus,ptaup,peta);
            }
        } else {
            if (READMOD) {
                readmod_viscac(prho,ppi,ptaup,peta);
            } else {
                model_viscac(prho,ppi,ptaup,peta);
            }
        }
    } else {
        if(!ACOUSTIC) {
            if (READMOD) {
                readmod_elastic(prho,ppi,pu);
            } else {
                model_elastic(prho,ppi,pu);
            }
        } else {
            if (READMOD) {
                readmod_acoustic(prho,ppi);
            } else {
                model_acoustic(prho,ppi);
            }
        }
    }
    /* check if the FD run will be stable and free of numerical dispersion */
    checkfd(FP, prho, ppi, pu, ptaus, ptaup, peta, hc, srcpos, nsrc, recpos, ntr_glob);
    /* calculate damping coefficients for CPMLs*/
    if(FW>0)
        PML_pro(d_x, K_x, alpha_prime_x, a_x, b_x, d_x_half, K_x_half, alpha_prime_x_half, a_x_half, b_x_half,
                d_y, K_y, alpha_prime_y, a_y, b_y, d_y_half, K_y_half, alpha_prime_y_half, a_y_half, b_y_half);
    MPI_Barrier(MPI_COMM_WORLD);
    SHOTINC=1;
    RECINC=1;
    switch(TIME_FILT) {
    case 1:
        F_LOW_PASS=F_LOW_PASS_START;
        break;
    /*read frequencies from file*/
    case 2:
        F_LOW_PASS_EXT=filter_frequencies(&nfrq);
        F_LOW_PASS=F_LOW_PASS_EXT[FREQ_NR];
        break;
    }
    /* Save old SOURCE_SHAPE, which is needed for STF */
    SOURCE_SHAPE_OLD = SOURCE_SHAPE;
    if(WAVETYPE==2 || WAVETYPE==3) SOURCE_SHAPE_OLD_SH=SOURCE_SHAPE_SH;
    nt_out=10000;
    if(!VERBOSE) nt_out=1e5;
}

void start_fullwaveform_iteration_loop() {
    for(iter=1; iter<=ITERMAX; iter++) {
        // While loop for Wolfe step length search            
        while(FWI_run || steplength_search || gradient_optimization) {
            // Calculate Misfit and gradient          
            if(FWI_run) {
                if(INVTYPE==2) {
                    // Start of loop over shots 
                    for (ishot=1; ishot<=nshots; ishot+=SHOTINC) {
                        // Start of inversion of source time function 
                        if(((INV_STF==1)&&( (iter==1) || (do_stf==1) )) && (gradient_optimization==1) ) {
                            // calculate wavelet */
                            // loop over timesteps ( STF ) 
                            // Inversion of source time function 
                        }
                        // determine source position on grid 
                        // Use STF wavelet  
                        // calculate wavelet 
                        // Time Domain Filtering 
                        // loop over timesteps (forward model) 
                        // start of inversion process 
                        if(FORWARD_ONLY==0) {
                            // Calculate residuals 
                            if ((ntr > 0)) {
                                if(WAVETYPE==1 || WAVETYPE==3) {
                                    // read seismic data from SU file vx 
                                    // read seismic data from SU file vy 
                                    // read seismic data from SU file p 
                                }
                                // read seismic data from SU file vz 
                            }
                            // start loop over shots at receiver positions (backward model) 
                            for (irec=1; irec<=nshots1; irec+=RECINC) { /* loop over shots at receiver positions */
                                // Start loop over timesteps (backpropagation) 
                                for (nt=1; nt<=NT; nt++) {
                                    // Check if simulation is still stable SH 
                                    // update of particle velocities 
                                    // exchange of particle velocities between PEs 
                                    // Applying free surface condition 
                                    // Calculate convolution for every DTINV time step 
                                }
                            }
                            // calculate gradient direction pi 
                            // calculate gradient direction u PSV 
                            // calculate gradient direction u SH 
                            // calculate gradient direction rho PSV 
                            // calculate gradient direction rho SH
                            // calculate and apply energy preconditioning   
                            // Apply taper SH  
                            // Apply taper PSV 
                            // Summing up the gradient for all shots PSV 
                            // Summing up the gradient for all shots SH 
                        }
                    }
                }
                if(FORWARD_ONLY==0) {
                    // Applying and output approx hessian  
                    // Set gradient to zero if no inversion     
                    // calculate L2 norm of all CPUs 
                }
            }
            // Gradient optimization with PCG or L-BFGS        
            if(gradient_optimization==1 && FORWARD_ONLY==0) {
                // Preconditioned Conjugate Gradient Method (PCG)  
                // Beginn Joint Inversion PSV and SH 
            }
            // Limited Memory - Broyden-Fletcher-Goldfarb-Shanno  
            if(GRAD_METHOD==2) {
                // TAPER       
            }
            //   Wolfe condition: Check and step length search     
        }

        // Wolfe condition: Model update
        
        // Step length estimation:
        if((FORWARD_ONLY==0)  && (!WOLFE_CONDITION || (WOLFE_CONDITION && GRAD_METHOD==2 && iter==LBFGS_iter_start))) {
            // Beginn Search three step lengths for parabolic search 
            while((step2!=1)||(step1!=1)) {
                // calculate 3 L2 values 
                for (itest=itests; itest<=iteste; itest++) {
                    // start of loop over shots (test forward) 
                    for (ishot=TESTSHOT_START; ishot<=TESTSHOT_END; ishot=ishot+TESTSHOT_INCR) {
                        // determine source position on grid 
                        // calculate wavelet 
                        // Start of Time Domain Filtering 
                        // End of Time Domain Filtering 
                        // loop over timesteps (forward model) step length
                        // end loop over timesteps (forward model) step length  
                        // Calculate residuals
                        if(ntr > 0) {
                            if(WAVETYPE==1 || WAVETYPE==3) {
                                // read seismic data from SU file vx
                                // read seismic data from SU file vy 
                                // read seismic data from SU file p 
                            }
                            // read seismic data from SU file vz 
                        }
                    }
                    // end of loop over shots (test forward) 
                }
                // end of L2 test     
            }
            // End Search three step lengths for parabolic search  
        } // end of if(FORWARD_ONLY!=4) 

        // smoothing the models vp, vs and rho 
        // Check abort criteriums 
        if (iter>min_iter_help) {
            // Check when NO workflow and TIME_FILT==0    
            // Check when Workflow is used          
            // Check when Workflow is NOT used and TIME_FILT==1  
            // Check when Workflow is NOT used and TIME_FILT==2 
        }
    }
}

int main(int argc, char **argv)
{
    init_vars();
    start_fullwaveform_iteration_loop();
    free_memory();
    return 0;
}
