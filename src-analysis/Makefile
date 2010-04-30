SHELL=/bin/sh

CC=gfortran -O3
LIBRARIES= -lfftw3 -lblas -llapack 

OBJECTS = \
	main.o \
	readin.o \
	bare_dispersion.o \
	fft_3D_lattice.o \
	h0_generate.o \
	dyson.o \
	pade_evaluate.o

quasi : $(OBJECTS)
	$(CC) -o quasi.exe $(OBJECTS) $(LIBRARIES)

main.o : main.F constants.F
	$(CC) -c main.F 

readin.o : readin.F constants.F
	$(CC) -c readin.F 

bare_dispersion.o : bare_dispersion.F constants.F
	$(CC) -c bare_dispersion.F 

fft_3D_lattice.o : fft_3D_lattice.F constants.F fftw3.f
	$(CC) -c fft_3D_lattice.F

h0_generate.o : h0_generate.F constants.F
	$(CC) -c h0_generate.F

dyson.o : dyson.F constants.F
	$(CC) -c dyson.F

pade_evaluate.o : pade_evaluate.F constants.F
	$(CC) -c pade_evaluate.F

clean :
	rm $(OBJECTS) 
