SHELL=/bin/sh

CC=gfortran -O3
LIBRARIES= -lfftw3 -lblas -llapack 

OBJECTS = \
	main.o \
	readin.o \
	bare_dispersion.o \
	fft_3D_lattice.o \
	init_sigma1.o \
	init_pair_wave.o \
	psi_transform.o \
	h0_generate.o \
	calc_g_tau0.o \
	sigma_first.o \
	pair_wave.o \
	sigma_input.o \
	pade_evaluate.o \
	fft_3D.o \
	dyson.o \
	g_rtau.o \
	sigma_out.o \

quasi.exe : $(OBJECTS)
	$(CC) -o quasi.exe $(OBJECTS) $(LIBRARIES)

main.o : main.F main_defs.F constants.F
	$(CC) -c main.F 

readin.o : readin.F constants.F
	$(CC) -c readin.F $(defines)

bare_dispersion.o : init/bare_dispersion.F constants.F
	$(CC) -c init/bare_dispersion.F $(defines)

fft_3D_lattice.o : ffts/fft_3D_lattice.F constants.F
	$(CC) -c ffts/fft_3D_lattice.F $(defines)

generate_tau_eps_omega.o : init/generate_tau_eps_omega.F constants.F
	$(CC) -c init/generate_tau_eps_omega.F $(defines)

analytic_functions.o : init/analytic_functions.F constants.F
	$(CC) -c init/analytic_functions.F $(defines)

a_integrals.o : init/a_integrals.F constants.F
	$(CC) -c init/a_integrals.F $(defines)

init_sigma1.o : init/init_sigma1.F constants.F
	$(CC) -c init/init_sigma1.F $(defines)

init_pair_wave.o : init/init_pair_wave.F constants.F
	$(CC) -c init/init_pair_wave.F $(defines)

gamma0_define.o : init/gamma0_define.F constants.F
	$(CC) -c init/gamma0_define.F $(defines)

psi_transform.o : green_function/psi_transform.F constants.F
	$(CC) -c green_function/psi_transform.F $(defines)

h0_generate.o : green_function/h0_generate.F constants.F
	$(CC) -c green_function/h0_generate.F  $(defines)

calc_g_tau0.o : green_function/calc_g_tau0.F constants.F
	$(CC) -c green_function/calc_g_tau0.F  $(defines)

sigma_first.o : self-energy-calculate/sigma_first.F constants.F
	$(CC) -c self-energy-calculate/sigma_first.F $(defines)

pair_wave.o : green_function/pair_wave.F constants.F
	$(CC) -c green_function/pair_wave.F $(defines)

kinetic_energy.o : thermodynamics/kinetic_energy.F constants.F
	$(CC) -c thermodynamics/kinetic_energy.F $(defines)

analyze_psi.o : green_function/analyze_psi.F constants.F
	$(CC) -c green_function/analyze_psi.F $(defines)

analyze_psi_1D.o : green_function/analyze_psi_1D.F constants.F
	$(CC) -c green_function/analyze_psi_1D.F $(defines)

analyze_psi_2D_tet.o : green_function/analyze_psi_2D_tet.F constants.F
	$(CC) -c green_function/analyze_psi_2D_tet.F $(defines)

analyze_psi_2D_ortho.o : green_function/analyze_psi_2D_ortho.F constants.F
	$(CC) -c green_function/analyze_psi_2D_ortho.F $(defines)

# analyze_psi_3D_cubic.o : green_function/analyze_psi_3D_cubic.F constants.F
# 	$(CC) -c green_function/analyze_psi_3D_cubic.F $(defines)

# analyze_psi_3D_tet.o : green_function/analyze_psi_3D_tet.F constants.F
# 	$(CC) -c green_function/analyze_psi_3D_tet.F $(defines)

sigma1_out.o : self-energy-read-write/sigma1_out.F constants.F
	$(CC) -c self-energy-read-write/sigma1_out.F $(defines)

sigma_input.o : self-energy-read-write/sigma_input.F constants.F
	$(CC) -c self-energy-read-write/sigma_input.F $(defines)

discon_lat.o : green_function/discon_lat.F constants.F
	$(CC) -c green_function/discon_lat.F $(defines)

discontinuities.o : green_function/discontinuities.F constants.F
	$(CC) -c green_function/discontinuities.F $(defines)

green_param_lat.o : green_function/green_param_lat.F constants.F
	$(CC) -c green_function/green_param_lat.F $(defines)

green_parameter.o : green_function/green_parameter.F constants.F
	$(CC) -c green_function/green_parameter.F $(defines)

dyson.o : green_function/dyson.F constants.F
	$(CC) -c green_function/dyson.F $(defines)

fft_3D.o : ffts/fft_3D.F constants.F
	$(CC) -c ffts/fft_3D.F $(defines)

parameter_init.o : init/parameter_init.F constants.F
	$(CC) -c init/parameter_init.F $(defines)

g_rtau.o : green_function/g_rtau.F constants.F
	$(CC) -c green_function/g_rtau.F $(defines)

g_minus_tau.o : green_function/g_minus_tau.F constants.F
	$(CC) -c green_function/g_minus_tau.F $(defines)

sigma_calc.o : self-energy-calculate/sigma_calc.F constants.F
	$(CC) -c self-energy-calculate/sigma_calc.F $(defines)

calc_new_alpha.o : self-energy-calculate/calc_new_alpha.F constants.F
	$(CC) -c self-energy-calculate/calc_new_alpha.F $(defines)

convergence_test.o : self-energy-calculate/convergence_test.F constants.F
	$(CC) -c self-energy-calculate/convergence_test.F $(defines)

chi_calc.o : chi/chi_calc.F constants.F
	$(CC) -c chi/chi_calc.F $(defines)

transform_chi.o : chi/transform_chi.F constants.F
	$(CC) -c chi/transform_chi.F $(defines)

t-generate.o : t-matrix/t-generate.F constants.F
	$(CC) -c t-matrix/t-generate.F $(defines)

t_transform.o : t-matrix/t_transform.F constants.F
	$(CC) -c t-matrix/t_transform.F $(defines)

tmat_param.o : t-matrix/tmat_param.F constants.F
	$(CC) -c t-matrix/tmat_param.F $(defines)

fft_4D.o : ffts/fft_4D.F constants.F
	$(CC) -c ffts/fft_4D.F $(defines)

col_dist_to_row_dist.o : ffts/col_dist_to_row_dist.F constants.F
	$(CC) -c ffts/col_dist_to_row_dist.F $(defines)

row_dist_to_col_dist.o : ffts/row_dist_to_col_dist.F constants.F
	$(CC) -c ffts/row_dist_to_col_dist.F $(defines)

sigma_out.o : self-energy-read-write/sigma_out.F constants.F
	$(CC) -c self-energy-read-write/sigma_out.F $(defines)

pade.o : self-energy-read-write/pade.F constants.F
	$(CC) -c self-energy-read-write/pade.F $(defines)

pade_evaluate.o : self-energy-read-write/pade_evaluate.F constants.F
	$(CC) -c self-energy-read-write/pade_evaluate.F $(defines)

clean :
	rm $(OBJECTS_COMMON) $(OBJECTS_SECOND) $(OBJECTS_THIRD)
