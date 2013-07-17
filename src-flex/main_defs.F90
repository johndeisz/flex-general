  ! MPI variables
  INTEGER rank, size
  ! Timing variables
  Real start_time, end_time
  Real last_it_time, this_it_time

  ! Parameters read from input file
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
  COMPLEX Lvec(1:3,0:nb-1, 0:nb-1)
  COMPLEX Lorb(1:3)

  REAL prfld_pert
  REAL h_pert(0:nb-1,1:3)
  REAL v_pert(0:nb-1)
  REAL prfld_eff
  REAL h_eff(0:nb-1,1:3)
  REAL v_pert_eff(0:nb-1)
  LOGICAL read_input, write_output
  CHARACTER*128 sigma_input_file, sigma_output_file
  INTEGER max_pade_order
  REAL sigma_tol
  INTEGER max_it
  REAL alpha
  INTEGER alpha_scheme

  !     Dispersion and bare vertex
  REAL ek_min
  COMPLEX gamma0_ph(0:16*nb*nb-1, 0:16*nb*nb-1)

  ! Time and frequency arrays
#ifdef SECOND_ORDER
  REAL tau(0:mp1), epsilon(0:mp1), omega(0:mp1)
#endif

  ! Variables, arrays for analytic functions
#ifdef SECOND_ORDER
  REAL x(0:1,0:1), y(0:1,0:1)
  COMPLEX q_epsilon(0:1,0:1,0:mp1)
  COMPLEX r_omega(0:1,0:1,0:mp1)
  REAL q_tau(0:1,0:1,0:mp1), q_mtau(0:1,0:1,0:mp1)
  REAL r_tau(0:1,0:1,0:mp1)
  COMPLEX a_int(0:1,0:1,0:1,0:1,0:mp1)
#endif


  ! Self-energy arrays
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1) :: sigma1, & 
       sigma1_old, delta_sigma1, delta_sigma1_old
#ifdef SECOND_ORDER
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:mp-1,0:nc-1) :: sigma, sigma_old
  COMPLEX chi(0:16*nb*nb-1,0:16*nb*nb-1,0:mp1,0:nc1)
  COMPLEX delta_sigma_e0_old(0:4*nb-1,0:4*nb-1,0:nc1)
#endif

  !     Chemical potential, density, and interaction strength
  REAL mu_old, mu_current, dn_dmu, mu_min, mu_max
  REAL density, density_old, delta_mu
  INTEGER ib, jb, is, js, nu1, nu2
  LOGICAL density_converged, ft_sigma
  REAL m_psi, mag
  REAL denb(0:nb-1), magb(0:nb-1,1:3)

  !     Self-consistency loop variables
  INTEGER iteration, density_iteration
  LOGICAL sigma_converged

  ! Green function and related
  COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)

  COMPLEX g_tau0(0:4*nb-1,0:4*nb-1,0:nl-1)
  COMPLEX g_tau0_local(0:4*nb-1,0:4*nb-1)
#ifdef SECOND_ORDER
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:mp1,0:nc1) :: g, g_mtau
  COMPLEX, dimension (0:4*nb-1,0:4*nb-1,0:nc-1) :: &
       delta_g_r, delta_g_k, delta_gp_r, delta_gp_k
  COMPLEX c_r(0:1,0:1,0:4*nb-1,0:4*nb-1,0:nc-1)
  COMPLEX dominant_chi_eigenvector(0:16*nb*nb-1)
  REAL overall_eigenvalue_max
  INTEGER dominant_chi_index(0:1)
#endif

  ! FFT variables
  INTEGER isignv(0:3)

  ! Output
  CHARACTER*10 convergence_text

  ! Thermodynamic quantities
  REAL kinetic, e_hf, cor_energy, pair_energy
  REAL mag_energy, pot_energy
  COMPLEX so_energy, glocal
#ifdef SECOND_ORDER
  COMPLEX tr_sig_g
  COMPLEX d_r(0:1,0:1,0:16*nb*nb-1,0:16*nb*nb-1,0:nc1)
#endif /* SECOND_ORDER */

#ifdef USE_MPI
  INTEGER ierr
#endif /* USE_MPI */
