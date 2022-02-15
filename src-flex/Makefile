SHELL=/bin/sh

CC=gfortran -O3
# MPICC=mpif90 -O3
MPICC=/usr/lib64/mpich/bin/mpif90 -O3 -mcmodel=medium -fallow-argument-mismatch
LIBRARIES= -lfftw3 -lblas -llapack 

ifdef SYS
  include Make.sys/Make.$(SYS)
endif

OBJECTS_COMMON = \
	constants.o \
	bare_dispersion.o \
	h0.o \
	main.o \
	init_environ.o \
	readin.o \
	fft_3D_lattice.o \
	init_sigma1.o \
	init_pair_wave.o \
	psi_transform.o \
	gamma0_define.o \
	effective_field.o \
	calc_g_tau0.o \
	sigma_first.o \
	pair_wave.o \
	sigma1_out.o \
	pade_evaluate.o \
	sigma_input.o \
	kinetic_energy.o \
	analyze_psi.o \
	analyze_psi_exchange.o \
	analyze_psi_spin.o \
	analyze_psi_band.o \
	analyze_psi_1D.o \
	analyze_psi_2D_ortho.o \
	pair_wave_test.o \
	angular_matrices.o \
	hfock_energy.o \
	pair_energy.o
#	analyze_psi_3D_cubic.o \
#	analyze_psi_3D_tet.o \
#	sig_output.o \
#	chi_diag.o \
#	diag.o \



OBJECTS_SECOND = \
	generate_tau_eps_omega.o \
	analytic_functions.o \
	a_integrals.o \
	discontinuities.o \
	fft_3D.o \
	parameter_init.o \
	green_param_lat.o \
	dyson.o \
	calc_g_tau0_2nd.o \
	g_rtau.o \
	fft_4D.o \
	row_dist_to_col_dist.o \
	col_dist_to_row_dist.o \
	g_minus_tau.o \
	green_parameter.o \
	sigma_calc.o \
	chi_calc.o \
	t-generate.o \
	tmat_param.o \
	calc_new_alpha.o \
	convergence_test.o \
	pade.o \
	sigma_out.o \
        symmetrize.o \
	eval_tr_sigph_g.o \
	l_ints.o \

OBJECTS_THIRD = \
	transform_chi.o \
	t_transform.o \

ifeq ($(PREC),single)
        prec_flag = -DSINGLE_PREC         
	prec_tag = single
else
        prec_flag = -DDOUBLE_PREC
        prec_tag = double
endif

flex : defines = -DSECOND_ORDER -DTHIRD_ORDER -DFLEX $(prec_flag) 
flex : OBJECTS = $(OBJECTS_COMMON) $(OBJECTS_SECOND) $(OBJECTS_THIRD)
flex : $(OBJECTS_COMMON) $(OBJECTS_SECOND) $(OBJECTS_THIRD)
	$(CC) -o flex-$(prec_tag) $(OBJECTS) $(LIBRARIES) 

flex-mpi : defines = -DSECOND_ORDER -DTHIRD_ORDER -DFLEX -DUSE_MPI $(prec_flag) 
flex-mpi : CC=$(MPICC)
flex-mpi : OBJECTS = $(OBJECTS_COMMON) $(OBJECTS_SECOND) $(OBJECTS_THIRD)
flex-mpi : $(OBJECTS_COMMON) $(OBJECTS_SECOND) $(OBJECTS_THIRD)
	$(CC) -o flex-mpi-$(prec_tag) $(OBJECTS) $(LIBRARIES)

3rd-order : defines = -DSECOND_ORDER -DTHIRD_ORDER $(prec_flag) 
3rd-order : OBJECTS = $(OBJECTS_COMMON) $(OBJECTS_SECOND) $(OBJECTS_THIRD)
3rd-order :  $(OBJECTS_COMMON) $(OBJECTS_SECOND) $(OBJECTS_THIRD)
	$(CC) -o 3rd-order-$(prec_tag) $(OBJECTS) $(LIBRARIES) 

3rd-order-mpi : defines = -DSECOND_ORDER -DTHIRD_ORDER -DUSE_MPI $(prec_flag) 
3rd-order-mpi : CC=$(MPICC)
3rd-order-mpi : OBJECTS = $(OBJECTS_COMMON) $(OBJECTS_SECOND) $(OBJECTS_THIRD)
3rd-order-mpi :  $(OBJECTS_COMMON) $(OBJECTS_SECOND) $(OBJECTS_THIRD)
	$(MPICC) -o 3rd-order-mpi-$(prec_tag) $(OBJECTS) $(LIBRARIES)

2nd-order : defines = -DSECOND_ORDER $(prec_flag) 
2nd-order : OBJECTS = $(OBJECTS_COMMON) $(OBJECTS_SECOND)
2nd-order :  $(OBJECTS_COMMON) $(OBJECTS_SECOND)
	$(CC) -o 2nd-order-$(prec_tag) $(OBJECTS) $(LIBRARIES) 

2nd-order-mpi : defines = -DSECOND_ORDER -DUSE_MPI $(prec_flag) 
2nd-order-mpi : CC=$(MPICC)
2nd-order-mpi : OBJECTS = $(OBJECTS_COMMON) $(OBJECTS_SECOND)
2nd-order-mpi :  $(OBJECTS_COMMON) $(OBJECTS_SECOND)
	$(MPICC) -o 2nd-order-mpi-$(prec_tag) $(OBJECTS) $(LIBRARIES)

1st-order : defines = $(prec_flag) 
1st-order : OBJECTS = $(OBJECTS_COMMON) 
1st-order : $(OBJECTS_COMMON)
	$(CC) -o 1st-order-$(prec_tag) $(OBJECTS) $(LIBRARIES)

constants.o : constants.F90
	$(CC) -c constants.F90 $(defines)

main.o : main.F90 main_defs.F90 constants.F90
	$(CC) -c main.F90 $(defines)

init_environ.o : init/init_environ.F90 constants.F90
	$(CC) -c init/init_environ.F90 $(defines)

readin.o : init/readin.F90 constants.F90
	$(CC) -c init/readin.F90 $(defines)

bare_dispersion.o : init/bare_dispersion.F90 constants.F90
	$(CC) -c init/bare_dispersion.F90 $(defines)

fft_3D_lattice.o : ffts/fft_3D_lattice.F90 constants.F90
	$(CC) -c ffts/fft_3D_lattice.F90 $(defines)

generate_tau_eps_omega.o : init/generate_tau_eps_omega.F90 constants.F90
	$(CC) -c init/generate_tau_eps_omega.F90 $(defines)

analytic_functions.o : init/analytic_functions.F90 constants.F90
	$(CC) -c init/analytic_functions.F90 $(defines)

a_integrals.o : init/a_integrals.F90 constants.F90
	$(CC) -c init/a_integrals.F90 $(defines)

init_sigma1.o : init/init_sigma1.F90 constants.F90
	$(CC) -c init/init_sigma1.F90 $(defines)

init_pair_wave.o : init/init_pair_wave.F90 constants.F90
	$(CC) -c init/init_pair_wave.F90 $(defines)

gamma0_define.o : init/gamma0_define.F90 constants.F90
	$(CC) -c init/gamma0_define.F90 $(defines)

psi_transform.o : green_function/psi_transform.F90 constants.F90
	$(CC) -c green_function/psi_transform.F90 $(defines)

effective_field.o : green_function/effective_field.F90 constants.F90
	$(CC) -c green_function/effective_field.F90  $(defines)

h0.o : green_function/h0.F90 constants.F90
	$(CC) -c green_function/h0.F90  $(defines)

calc_g_tau0.o : green_function/calc_g_tau0.F90 constants.F90
	$(CC) -c green_function/calc_g_tau0.F90  $(defines)

sigma_first.o : self-energy-calculate/sigma_first.F90 constants.F90
	$(CC) -c self-energy-calculate/sigma_first.F90 $(defines)

pair_wave.o : green_function/pair_wave.F90 constants.F90
	$(CC) -c green_function/pair_wave.F90 $(defines)

kinetic_energy.o : thermodynamics/kinetic_energy.F90 constants.F90
	$(CC) -c thermodynamics/kinetic_energy.F90 $(defines)

analyze_psi.o : green_function/analyze_psi.F90 constants.F90
	$(CC) -c green_function/analyze_psi.F90 $(defines)

analyze_psi_exchange.o : green_function/analyze_psi_exchange.F90 constants.F90
	$(CC) -c green_function/analyze_psi_exchange.F90 $(defines)

analyze_psi_spin.o : green_function/analyze_psi_spin.F90 constants.F90
	$(CC) -c green_function/analyze_psi_spin.F90 $(defines)

analyze_psi_band.o : green_function/analyze_psi_band.F90 constants.F90
	$(CC) -c green_function/analyze_psi_band.F90 $(defines)

analyze_psi_1D.o : green_function/analyze_psi_1D.F90 constants.F90
	$(CC) -c green_function/analyze_psi_1D.F90 $(defines)

analyze_psi_2D_ortho.o : green_function/analyze_psi_2D_ortho.F90 constants.F90
	$(CC) -c green_function/analyze_psi_2D_ortho.F90 $(defines)

# analyze_psi_3D_cubic.o : green_function/analyze_psi_3D_cubic.F90 constants.F90
# 	$(CC) -c green_function/analyze_psi_3D_cubic.F90 $(defines)

# analyze_psi_3D_tet.o : green_function/analyze_psi_3D_tet.F90 constants.F90
# 	$(CC) -c green_function/analyze_psi_3D_tet.F90 $(defines)

sigma1_out.o : self-energy-read-write/sigma1_out.F90 constants.F90
	$(CC) -c self-energy-read-write/sigma1_out.F90 $(defines)

pade_evaluate.o : self-energy-read-write/pade_evaluate.F90 constants.F90
	$(CC) -c self-energy-read-write/pade_evaluate.F90 $(defines)

sigma_input.o : self-energy-read-write/sigma_input.F90 constants.F90
	$(CC) -c self-energy-read-write/sigma_input.F90 $(defines)

# discon_lat.o : green_function/discon_lat.F90 constants.F90
#	$(CC) -c green_function/discon_lat.F90 $(defines)

discontinuities.o : green_function/discontinuities.F90 constants.F90
	$(CC) -c green_function/discontinuities.F90 $(defines)

green_param_lat.o : green_function/green_param_lat.F90 constants.F90
	$(CC) -c green_function/green_param_lat.F90 $(defines)

green_parameter.o : green_function/green_parameter.F90 constants.F90
	$(CC) -c green_function/green_parameter.F90 $(defines)

dyson.o : green_function/dyson.F90 constants.F90
	$(CC) -c green_function/dyson.F90 $(defines)

calc_g_tau0_2nd.o : green_function/calc_g_tau0_2nd.F90 constants.F90
	$(CC) -c green_function/calc_g_tau0_2nd.F90 $(defines)

fft_3D.o : ffts/fft_3D.F90 constants.F90
	$(CC) -c ffts/fft_3D.F90 $(defines)

parameter_init.o : init/parameter_init.F90 constants.F90
	$(CC) -c init/parameter_init.F90 $(defines)

g_rtau.o : green_function/g_rtau.F90 constants.F90
	$(CC) -c green_function/g_rtau.F90 $(defines)

g_minus_tau.o : green_function/g_minus_tau.F90 constants.F90
	$(CC) -c green_function/g_minus_tau.F90 $(defines)

sigma_calc.o : self-energy-calculate/sigma_calc.F90 constants.F90
	$(CC) -c self-energy-calculate/sigma_calc.F90 $(defines)

calc_new_alpha.o : self-energy-calculate/calc_new_alpha.F90 constants.F90
	$(CC) -c self-energy-calculate/calc_new_alpha.F90 $(defines)

convergence_test.o : self-energy-calculate/convergence_test.F90 constants.F90
	$(CC) -c self-energy-calculate/convergence_test.F90 $(defines)

chi_calc.o : chi/chi_calc.F90 constants.F90
	$(CC) -c chi/chi_calc.F90 $(defines)

transform_chi.o : chi/transform_chi.F90 constants.F90
	$(CC) -c chi/transform_chi.F90 $(defines)

t-generate.o : t-matrix/t-generate.F90 constants.F90
	$(CC) -c t-matrix/t-generate.F90 $(defines)

t_transform.o : t-matrix/t_transform.F90 constants.F90
	$(CC) -c t-matrix/t_transform.F90 $(defines)

tmat_param.o : t-matrix/tmat_param.F90 constants.F90
	$(CC) -c t-matrix/tmat_param.F90 $(defines)

fft_4D.o : ffts/fft_4D.F90 constants.F90
	$(CC) -c ffts/fft_4D.F90 $(defines)

col_dist_to_row_dist.o : ffts/col_dist_to_row_dist.F90 constants.F90
	$(CC) -c ffts/col_dist_to_row_dist.F90 $(defines)

row_dist_to_col_dist.o : ffts/row_dist_to_col_dist.F90 constants.F90
	$(CC) -c ffts/row_dist_to_col_dist.F90 $(defines)

sigma_out.o : self-energy-read-write/sigma_out.F90 constants.F90
	$(CC) -c self-energy-read-write/sigma_out.F90 $(defines)

pade.o : self-energy-read-write/pade.F90 constants.F90
	$(CC) -c self-energy-read-write/pade.F90 $(defines)

pair_wave_test.o : green_function/pair_wave_test.F90 constants.F90
	$(CC) -c green_function/pair_wave_test.F90 $(defines)

symmetrize.o : green_function/symmetrize.F90 constants.F90
	$(CC) -c green_function/symmetrize.F90 $(defines)

angular_matrices.o : green_function/angular_matrices.F90 constants.F90
	$(CC) -c green_function/angular_matrices.F90 $(defines)

hfock_energy.o : thermodynamics/hfock_energy.F90 constants.F90
	$(CC) -c thermodynamics/hfock_energy.F90 $(defines)

eval_tr_sigph_g.o : thermodynamics/eval_tr_sigph_g.F90 constants.F90
	$(CC) -c thermodynamics/eval_tr_sigph_g.F90 $(defines)

l_ints.o : init/l_ints.F90 constants.F90
	$(CC) -c init/l_ints.F90 $(defines)

pair_energy.o : thermodynamics/pair_energy.F90 constants.F90
	$(CC) -c thermodynamics/pair_energy.F90 $(defines)

clean :
	rm $(OBJECTS_COMMON) $(OBJECTS_SECOND) $(OBJECTS_THIRD) constants.mod h_zero.mod pade_eval.mod bare_dispersion.mod green_param_lat.mod
