#include "../convert.F"

      subroutine readin(t, flux, prfld, h,
     $   target_density, density_tol, mu, uu, up, uj, 
     $   ed, tij, prfld_pert, h_pert,
     $   read_input, sigma_input_file,
     $   write_output, sigma_output_file, max_pade_order,
     $   sigma_tol,  max_it, alpha, alpha_scheme)

#include "../constants.F"

#ifdef USE_MPI
      include 'mpif.h'
#endif

c Parameters to be read and returned
      REAL t
      REAL flux(1:3)
      REAL prfld
      REAL h(1:3)
      REAL target_density, density_tol
      REAL mu 

      REAL uu, up, uj 

      REAL ed(0:nb-1)
      COMPLEX tij(0:nb-1,0:nb-1,0:nl-1)
      REAL prfld_pert
      REAL h_pert(1:3)
      LOGICAL read_input, write_output
      CHARACTER*128 sigma_input_file, sigma_output_file
      INTEGER max_pade_order
      REAL sigma_tol
      INTEGER max_it
      REAL alpha
      INTEGER alpha_scheme

c MPI variables
      INTEGER rank
      INTEGER ierr

c     other variable
      INTEGER icode
      INTEGER i
      REAL tmp_hop(1:5)

      INTEGER ib, ibp, id, idp
      INTEGER ix, iy, iz, k
      INTEGER max_x, max_y, max_z

#ifdef USE_MPI
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
#else
      rank = 0
#endif /* USE_MPI */

      if (rank .eq. 0) then
        write(6,*) 'lcx = ', lcx
        write(6,*) 'lcy = ', lcy
        write(6,*) 'lcz = ', lcz
        write(6,*) 'm = ', m
        write(6,*) 'llx = ', llx
        write(6,*) 'lly = ', lly
        write(6,*) 'llz = ', llz
        write(6,*)
      endif

      if (rank .eq. 0) then

        read(5,*)
        read(5,*) t
        read(5,*) prfld
        read(5,*) h
        read(5,*) target_density, density_tol
        read(5,*) mu
        
        write(6,*) "Basic Thermodynamic Parameters" 
        write(6,*) "temperature in K = ", t  

        write(6,*) "external pairing field = ", prfld
        write(6,*) "uniform magnetic field, hx = ", h(1)
        write(6,*) "                        hy = ", h(2)
        write(6,*) "                        hz = ", h(3)       
        write(6,*) "target electron density = ", target_density
        write(6,*) "density tolerance = ", density_tol
        write(6,*) "default initial chemical potential = ", mu

        open (unit=20, file='converged_mu',status='old',
     $     iostat=icode)

        if (icode .ne. 0) then
          write (6,*) 'File converged_mu not found, ',
     $       'using default value for mu.'
        else 
          write (6,*) 'File converged_mu found.'
          read (20,*) mu
          write (6,*) 'Using mu = ',mu
        endif
        close (20)
        write(6,*) ' '

      endif

#ifdef USE_MPI
      call MPI_Bcast(t, 1, MPI_REAL, 0, 
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(prfld, 1, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(h, 3, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(target_density, 1, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(density_tol, 1, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(mu, 1, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
#endif

c----------------------Interaction Parameters ------------------------

      if (rank .eq. 0) then
        read(5,*)
        read(5,*)
        read(5,*) uu, up, uj 

        write(6,*) "Interaction parameters"
        write(6,*) "Intraorbital Coulomb = ", uu
        write(6,*) "Interorbital Coulomb = ", up
        write(6,*) "Interorbital exchange = ", uj

      endif

#ifdef USE_MPI
      call MPI_Bcast(uu, 1, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(up, 1, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(uj, 1, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
#endif 



c-----Perturbations: applied and relaxed over initial 10 iterations -------

      if (rank .eq. 0) then
        read(5,*)
        read(5,*)
        read(5,*) prfld_pert
        read(5,*) h_pert

        write(6,*) "Preturbations: applied and relaxed over",
     $   " initial 10 iterations"
        write(6,*) "Artificial pair field = ", prfld_pert
        write(6,*) "Artificial uniform magnetic field, hx = ", 
     $     h_pert(1)
        write(6,*) "Artificial uniform magnetic field, hy = ", 
     $     h_pert(2)
        write(6,*) "Artificial uniform magnetic field, hz = ", 
     $     h_pert(3)
      endif

#ifdef USE_MPI
      call MPI_Bcast(prfld_pert, 1, MPI_REAL, 0, 
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(h_pert, 3, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
#endif

c------------------------Sigma Input and Output --------------------------
      
      if (rank .eq. 0) then

        read(5,*)
        read(5,*)
        read(5,*) read_input, sigma_input_file
        read(5,*) write_output, sigma_output_file
        read(5,*) max_pade_order

        write(6,*) "Sigma input and output"
        write(6,*) "read sigma input file = ", read_input
        write(6,*) "Sigma input file = ", sigma_input_file
        write(6,*) "write sigma input file = ", write_output
        write(6,*) "Sigma output file = ", sigma_output_file
        write(6,*) "Maximum pade order for output sigma = ",
     $     max_pade_order

      endif

#ifdef USE_MPI
      call MPI_Bcast(read_input, 1, MPI_LOGICAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(write_output, 1, MPI_INTEGER, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(max_pade_order, 1, MPI_INTEGER, 0,
     $   MPI_COMM_WORLD, ierr)
#endif

c----------------------------FEA parameters ------------------------------

      if (rank .eq. 0) then

        read(5,*)
        read(5,*)
        read(5,*) sigma_tol,  max_it
        read(5,*) alpha, alpha_scheme

        write(6,*) "Flex parameters"
        write(6,*) "Tolerance in sigma = ", sigma_tol
        write(6,*) "Maximum iterations = ", max_it 
        write(6,*) "alpha = ", alpha, " alpha scheme = ",
     $     alpha_scheme

      endif

#ifdef USE_MPI
      call MPI_Bcast(sigma_tol, 1, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(max_it, 1, MPI_INTEGER, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(alpha, 1, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call  MPI_Bcast(alpha_scheme, 1, MPI_INTEGER, 0,
     $   MPI_COMM_WORLD, ierr)
#endif /* USE_MPI */



c-----Orbital energies and hopping matrix elements -------------------

      if (rank .eq. 0) then

        read(5,*) 
        read(5,*)
        read(5,*)
        read(5,*) flux
        write(6,*) "x, y and z flux = ", flux

        read(5,*)
        read(5,*)
        do ib = 0, nb-1
          read(5,*) id, ed(ib)
        enddo

        write(6,*) 'orbital energies '
        do ib = 0, nb-1
          write(6,*) ib, ed(ib)
        enddo
        write(6,*)
      endif

      tij = 0.0d0

      if (rank .eq. 0) then

        if (llx .gt. 2) then
          max_x = 2
        else
          max_x = llx - 1
        endif

        if (lly .gt. 2) then
          max_y = 2
        else
          max_y = lly - 1
        endif

        if (llz .gt. 2) then
          max_z = 2
        else
          max_z = llz - 1
        endif
          

        do ix = -max_x, max_x
          do iy = -max_y, max_y
            do iz = -max_z, max_z

              read(5,*)
              read(5,*)
              write(6,*) 
              write(6,200) ix, iy, iz
              
              k = mod(ix+llx,llx) + mod(iy+lly,lly)*llx +
     $           mod(iz+llz,llz)*llx*lly

              do ib = 0, nb-1
                do ibp = 0, nb-1
                  read(5,*) id, idp, tij(ib,ibp,k)
                  write(6,300) ib, ibp, real(tij(ib,ibp,k)),
     $               aimag(tij(ib,ibp,k))
                enddo
              enddo

            enddo
          enddo
        enddo

      endif

#ifdef USE_MPI
      call MPI_Bcast(flux, 3, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ed, nb, MPI_REAL, 0,
     $   MPI_COMM_WORLD, ierr)
      call MPI_Bcast(tij, nb*nb*nl, MPI_COMPLEX, 0,
     $   MPI_COMM_WORLD, ierr)
#endif 

 200  format('------------------- [',i3,',',i3,',',i3,'] hopping',
     $   '------------------------')
 300  format(i3,',',i3,'  ','(',D16.9,',',D16.9,')')


      return
      end
