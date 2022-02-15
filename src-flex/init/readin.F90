#include "../convert.F90"

subroutine readin(t, flux, prfld, h, target_density, density_tol, mu, &
     uu, up, uj, ed, tij, prfld_pert, h_pert, v_pert, h_so, &
     read_input, sigma_input_file, write_output, sigma_output_file, & 
     max_pade_order, sigma_tol,  max_it, alpha, alpha_scheme)

  USE CONSTANTS
  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif

  ! Parameters to be read and returned
  REAL t
  REAL flux(1:3)
  REAL prfld
  REAL h(0:nb-1,1:3)
  REAL target_density, density_tol
  REAL mu 

  REAL uu, up, uj 

  REAL ed(0:nb-1)
  COMPLEX tij(0:nb-1,0:nb-1,-2:2,-2:2,-2:2)
  COMPLEX h_so(0:2*nb-1, 0:2*nb-1)
  REAL prfld_pert
  REAL h_pert_amp(1:3)
  REAL h_pert(0:nb-1,1:3)
  REAL v_pert_amp
  REAL v_pert(0:nb-1)
  LOGICAL read_input, write_output
  CHARACTER*128 sigma_input_file, sigma_output_file
  INTEGER max_pade_order
  REAL sigma_tol
  INTEGER max_it
  REAL alpha
  INTEGER alpha_scheme

  REAL pi
  REAL theta_flux
  COMPLEX phi_flux

  ! MPI variables
  INTEGER rank
  INTEGER ierr

  !     other variable
  INTEGER icode
  INTEGER i
  REAL tmp_hop(1:5)

  INTEGER ib, ibp, id, idp
  INTEGER ix, iy, iz, k
  INTEGER max_x, max_y, max_z

  INTEGER ibs, ibsp
  LOGICAL so_flag
  REAL so_amp
  REAL delta_strain
  INTEGER ig

#ifdef USE_MPI
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
#else
  rank = 0
#endif /* USE_MPI */

  pi = acos(-1.0d0)
  
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
     do ib = 0, nb - 1
        read(5,*) h(ib,:)
     enddo
     read(5,*) target_density, density_tol
     read(5,*) mu
        
     write(6,*) "Basic Thermodynamic Parameters" 
     write(6,*) "temperature in K = ", t  
     
     write(6,*) "external pairing field = ", prfld
     write(6,*) "Band-dependent magnetic fields (Tesla)"
     do ib = 0, nb-1
        write(6,*) h(ib,:)
     enddo
     write(6,*) "target electron density = ", target_density
     write(6,*) "density tolerance = ", density_tol
     write(6,*) "default initial chemical potential = ", mu

     open (unit=20, file='converged_mu',status='old', iostat=icode)

     if (icode .ne. 0) then
        write (6,*) 'File converged_mu not found, using default value for mu.'
     else 
        write (6,*) 'File converged_mu found.'
        read (20,*) mu
        write (6,*) 'Using mu = ',mu
     endif
     close (20)
     write(6,*) ' '

  endif

#ifdef USE_MPI
  call MPI_Bcast(t, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(prfld, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(h, 3*nb, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(target_density, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(density_tol, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mu, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
#endif

  !----------------------Interaction Parameters ------------------------

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
  call MPI_Bcast(uu, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(up, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(uj, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
#endif 



  !-----Perturbations: applied and relaxed over initial 10 iterations -------

  if (rank .eq. 0) then
     read(5,*)
     read(5,*)
     read(5,*) prfld_pert
     read(5,*) h_pert_amp
     read(5,*) v_pert_amp

     write(6,*) "Preturbations: applied and relaxed over initial 10 iterations"
     write(6,*) "Artificial pair field = ", prfld_pert
     write(6,*) "Random orb.-dependent mag. field amp (T) ", h_pert_amp
     do ib = 0, nb-1
        h_pert(ib,1) = h_pert_amp(1) * (rand()-0.5d0)
        h_pert(ib,2) = h_pert_amp(2) * (rand()-0.5d0)
        h_pert(ib,3) = h_pert_amp(3) * (rand()-0.5d0)
     enddo
     write(6,*) "Random orbital potential amplitude  = ", v_pert_amp
     do ib = 0, nb-1
        v_pert(ib) = v_pert_amp * (rand() - 0.5d0)
     enddo
  endif

#ifdef USE_MPI
  call MPI_Bcast(prfld_pert, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(h_pert, 3*nb, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(v_pert, nb, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
#endif

  !------------------------Sigma Input and Output --------------------------
      
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
     write(6,*) "Maximum pade order for output sigma = ", max_pade_order
     
  endif

#ifdef USE_MPI
  call MPI_Bcast(read_input, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(write_output, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(max_pade_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

!----------------------------FEA parameters ------------------------------

  if (rank .eq. 0) then

     read(5,*)
     read(5,*)
     read(5,*) sigma_tol,  max_it
     read(5,*) alpha, alpha_scheme
     
     write(6,*) "Flex parameters"
     write(6,*) "Tolerance in sigma = ", sigma_tol
     write(6,*) "Maximum iterations = ", max_it 
     write(6,*) "alpha = ", alpha, " alpha scheme = ", alpha_scheme
     
  endif

#ifdef USE_MPI
  call MPI_Bcast(sigma_tol, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(max_it, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(alpha, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call  MPI_Bcast(alpha_scheme, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif /* USE_MPI */

  !-----Orbital energies and hopping matrix elements -------------------

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
  delta_strain = 0.000d0
  
  if (rank .eq. 0) then

     write(6,*) "delta_strain = ", delta_strain

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
              
!!$              k = mod(ix+llx,llx) + mod(iy+lly,lly)*llx + &
!!$                   mod(iz+llz,llz)*llx*lly

              theta_flux = 2.0d0*pi*(flux(1)*dfloat(ix)/dfloat(llx) + &
                        flux(2)*dfloat(iy)/dfloat(lly) + &
                        flux(3)*dfloat(iz)/dfloat(llz) )

              phi_flux = cexp(cmplx(0.0d0,1.0d0)*theta_flux)

              do ib = 0, nb-1
                 do ibp = 0, nb-1
                    read(5,*) id, idp, tij(ib,ibp,ix,iy,iz)

                    if ( (ib .eq. 2) .and. (ibp .eq. 2)) then
                      
                      ig = abs(ix) - abs(iy)
                      if (ig .gt. 0) then
                        tij(ib,ibp,ix,iy,iz) = tij(ib,ibp,ix,iy,iz)*(1.0+delta_strain)
                      endif

                      if (ig .lt. 0) then
                        tij(ib,ibp,ix,iy,iz) = tij(ib,ibp,ix,iy,iz)*(1.0-delta_strain)
                      endif 

                     endif


                    write(6,300) ib, ibp, real(tij(ib,ibp,ix,iy,iz)), &
                         aimag(tij(ib,ibp,ix,iy,iz))

                    tij(ib,ibp,ix,iy,iz) = tij(ib,ibp,ix,iy,iz)*phi_flux
                 enddo
              enddo

           enddo
        enddo
     enddo

  endif

#ifdef USE_MPI
  call MPI_Bcast(flux, 3, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ed, nb, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(tij, nb*nb*125, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif 

  if (rank .eq. 0) then
     read(5,*)
     read(5,*)
     read(5,*)
     read(5,*) so_flag, so_amp
     write(6,*)
     write(6,*) 'so_flag = ', so_flag
     write(6,*) 'so_amp = ', so_amp
     if (so_flag) then
        read(5,*)
        write(6,*) 'nus   nusp  so_matrix'
        do ibs = 0, 2*nb-1
           do ibsp = 0, 2*nb-1
              read(5,*) id, idp, h_so(ibs,ibsp)
              write(6,300) ibs, ibsp, real(h_so(ibs,ibsp)), imag(h_so(ibs,ibsp))
           enddo
        enddo
        h_so = so_amp*h_so
     else
        h_so = cmplx(0.0d0, 0.0d0)
     endif
  endif

#ifdef USE_MPI
  call MPI_Bcast(h_so,4*nb*nb, MPI_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif

200 format('------------------- [',i3,',',i3,',',i3,'] hopping', &
         '------------------------')
300 format(i3,',',i3,'  ','(',D16.9,',',D16.9,')')

  return
end subroutine readin
