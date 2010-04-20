      program electronic_states

#include "constants.F"

      CHARACTER*128 sigma_input_file
      DOUBLE PRECISION t, mu, prfld
      DOUBLE PRECISION ed(0:nb-1), h(0:nb-1,1:3)
      DOUBLE COMPLEX tij(0:nb-1,0:nb-1,0:nl-1)
      DOUBLE COMPLEX psi(0:2*nb-1, 0:2*nb-1, 0:nl-1)
      DOUBLE COMPLEX sigma1(0:4*nb-1,0:4*nb-1)

      INTEGER ib
      DOUBLE COMPLEX ek(0:nb-1,0:nb-1,0:nl-1)

      write(6,*) "Self-energy file to be used"
      read(5,*) sigma_input_file

      call readin(sigma_input_file, t, mu, prfld, ed, h, 
     $   tij, psi, sigma1)

      write(6,*) "t = ", t
      write(6,*) "mu = ", mu
      write(6,*) "prfld = ", prfld
      do ib = 0, nb-1
        write(6,*) "ib, ed ", ib," ", ed(ib)
      enddo
      do ib = 0, nb-1
        write(6,*) "ib, h = ", ib, " ", h(ib,:)
      enddo

      t = kb*t

c     Bare bandstructure 
      call bare_dispersion(tij, ed, ek)

c$$$#ifdef SECOND_ORDER
c$$$      call sigma_input(sigma_input_file, uu, sigma1, psi,
c$$$     $   sigma, mu, epsilon)
c$$$#else
c$$$      call sigma_input(sigma_input_file, uu, sigma1, psi)
c$$$#endif
c$$$
c$$$:wq

c$$$c     Set the frequency limits and resolution
c$$$c     Determine if A(k,w) or N(w) is desired
c$$$
c$$$      call h0_generate(h0, ek, psi, h, prfld, mu, sigma1)
c$$$
c$$$      call dyson(rank, g, g_tau0, q_tau, q_epsilon, h0,
c$$$     $   sigma, epsilon, t, cl_k)


      stop
      end
