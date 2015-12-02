# Makefile for DENISE

#--------------------------------------------------------
# edit here:

# source code for model generation

#MODEL = hh.c
MODEL = ../genmod/1D_linear_gradient_visc.c
MODEL_AC = ../genmod/1D_linear_gradient_ac.c
MODEL_EL = ../genmod/1D_linear_gradient_el.c
EXEC= ../bin

# Description:
# CC = Compiler
# LFLAGS = Linker flag
# CFLAGS = Compiler flag

# LINUX with OpenMPI / IntelMPI and INTEL Compiler
# Use icc whenever possible, this will be much faster than gcc
CC=mpiicc
LFLAGS=-lm -lcseife -lstfinv -laff -lfourierxx -lfftw3 -lstdc++
CFLAGS=-O3
SFLAGS=-L./../contrib/libcseife -L./../contrib/bin
IFLAGS=-I./../contrib/libcseife -I./../contrib/header -I.

# LINUX with OpenMPI / IntelMPI and GCC Compiler
#CC=mpicc
#LFLAGS=-lm -lcseife -lstfinv -laff -lfourierxx -lfftw3 -lstdc++
#CFLAGS=-O3
#SFLAGS=-L./../contrib/libcseife -L./../contrib/bin
#IFLAGS=-I./../contrib/libcseife -I./../contrib/header -I.


# after this line, no further editing should be necessary
# --------------------------------------------------------

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $< $(IFLAGS)

SNAPMERGE_SCR = \
	snapmerge.c \
	merge.c \
	json_parser.c \
	read_par_json.c \
	readdsk.c \
	writedsk.c \
	util.c 


DENISE= \
	denise.c \
	stf.c \
	window_cos.c \
	calc_mat_change_test.c \
        calc_res.c \
	calc_misfit.c \
        calc_opt_step.c \
	calc_energy.c \
	checkfd.c \
	checkfd_ssg_elastic.c \
	checkfd_ssg_visc.c \
	conv_FD.c \
	count_killed_traces.c \
	psource.c \
	holbergcoeff.c\
	comm_ini.c\
	exchange_v.c \
	exchange_s.c \
	exchange_p.c \
	snap_ssg.c \
    snap_ssg_SH.c \
	seismo_ssg.c \
	surface_elastic_PML.c \
	surface_acoustic_PML.c \
	surface_PML.c \
	update_v_ssg.c \
	update_v_PML.c \
	update_v_PML_SH.c \
	update_v_acoustic_PML.c \
	prepare_update_s.c \
	update_p_PML.c \
	update_s_elastic_ssg.c \
	update_s_elastic_PML.c \
    update_s_elastic_PML_SH.c \
	update_s_visc_PML.c \
    update_s_visc_PML_SH.c \
	av_mue.c \
	av_rho.c \
	av_tau.c \
	median2D.c \
	exchange_par.c \
	info.c \
	inseis.c \
	inseis_source_wavelet.c \
	initproc.c \
	interpol.c \
	json_parser.c \
	LBFGS.c \
	smooth.c \
	$(MODEL) \
	$(MODEL_AC) \
	$(MODEL_EL) \
	matcopy.c \
	matcopy_elastic.c \
	matcopy_acoustic.c \
	mergemod.c \
	max_grad.c \
	note.c \
	norm.c \
	outseis.c \
	outseis_vector.c \
	outseis_glob.c \
	catseis.c \
	output_source_signal.c \
	PCG.c \
	PCG_SH.c \
	PML_pro.c \
	readdsk.c \
	read_par_json.c \
	read_checkpoint.c\
	readmod.c \
	readmod_elastic.c \
	readmod_acoustic.c \
    	receiver.c \
	rd_sour.c \
    read_workflow.c \
    apply_workflow.c \
	saveseis.c \
	saveseis_glob.c \
	save_checkpoint.c\
	sources.c \
	solvelin.c \
	spat_filt.c \
	splitsrc.c \
	splitsrc_back.c \
	splitrec.c \
	stalta.c \
	taper.c \
	taper_grad.c \
	taper_grad_shot.c \
	timedomain_filt.c \
	timedomain_filt_vector.c \
	time_window.c \
	util.c \
	wavelet.c \
	wavelet_stf.c \
	writemod.c \
	write_par.c \
	writedsk.c \
	zero_fdveps.c \
	zero_fdveps_ac.c \
	zero_fdveps_visc.c \
	calc_envelope.c \
    joint_inversion.c \
    matrix_operations.c \
    wolfe_condition.c \
    calc_hilbert.c \
	eprecond.c \
	eprecond1.c \
        filter_frequencies.c
	
SNAPMERGE_OBJ = $(SNAPMERGE_SCR:%.c=%.o)

DENISE_OBJ = $(DENISE:%.c=%.o)

snapmerge:	$(SNAPMERGE_OBJ)
	$(CC) $(SFLAGS) $(SNAPMERGE_OBJ) -o ../bin/snapmerge $(LFLAGS)

denise:		$(DENISE_OBJ)
	$(CC) $(SFLAGS) $(DENISE_OBJ) -o $(EXEC)/denise $(LFLAGS)


clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find ../genmod -name "*.o" -exec rm {} \; 
