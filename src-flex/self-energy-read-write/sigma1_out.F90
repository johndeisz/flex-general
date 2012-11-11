#include "../convert.F90"

subroutine sigma1_out(t, mu, flux, prfld, ed, h, &
     tij, uu, up, uj, psi, h_so, sigma1)

  USE CONSTANTS
  IMPLICIT NONE

  REAL t, mu, flux(1:3)
  REAL prfld, ed(0:nb-1), h(0:nb-1,1:3)
  COMPLEX tij(0:nb-1,0:nb-1,-2:2,-2:2,-2:2)
  REAL uu, up, uj
  COMPLEX psi(0:2*nb-1,0:2*nb-1,0:nl-1)
  COMPLEX h_so(0:2*nb-1, 0:2*nb-1)
  COMPLEX sigma1(0:4*nb-1,0:4*nb-1)
  INTEGER k, nu1, nu2, is1, is2, ia1, ia2
  INTEGER ib, ibp
  INTEGER ix, iy, iz
  INTEGER max_x, max_y, max_z

#if defined (FLEX)
  write(40,*) 'Fluctuation exchange approximation self-energy'
#elif defined (THIRD_ORDER)
  write(40,*) 'Third order perturbation theory self-energy'
#elif defined (SECOND_ORDER)
  write(40,*) 'Second order perturbation theory self-energy'
#else
  write(40,*) 'First order perturbation theory self-energy'
#endif /* defined (FLEX) */
      
  write(40,*) "for the periodic Hubbard model."
  write(40,*) "------------------------------------------"
  write(40,*) "------Input Parameters---------"
  write(40,*) "Nb (number of bands)"
  write(40,*) nb
  write(40,*) "LCX LCY LCZ   M   LLX  LLY  LLZ "
  write(40,*) lcx, "  ", lcy, "  ", lcz, "  ", m, " ", llx, " ", lly, " ", llz
  write(40,*) "Temp(K)  mu   flux_x  flux_y  flux_z"
  write(40,*)  t/kb, " ",  mu, " ", flux(1), " ", flux(2), " ", flux(3)
  write(40,*) "prfld"
  write(40,*) prfld
  write(40,*) "bare orbital levels"
  do ib = 0, nb-1
     write(40,*) ib, ed(ib)
  enddo
  write(40,*) "orbital-dependent magnetic fields"
  do ib = 0, nb-1
     write(40,*) ib, h(ib,:)
  enddo
  write(40,*)  "uu, up, uj"
  write(40,*)  uu, up, uj
  write(40,*) "Tight-binding parameters"

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

           write(40,200) ix, iy, iz

           k = mod(ix+llx,llx) + mod(iy+lly,lly)*llx + &
                mod(iz+llz,llz)*llx*lly
            
           do ib = 0, nb-1
              do ibp = 0, nb-1
                 write(40,300) ib, ibp, real(tij(ib,ibp,ix,iy,iz)), &
                      aimag(tij(ib,ibp,ix,iy,iz))
              enddo
           enddo

        enddo
     enddo
  enddo

  write(40,*)  
  
  write(40,*) "--Pair field wave function (k-space)-------"
  write(40,*) "Print flag (0 not printed, 1 printed)"

  if (abs(prfld) < 1.0e-9) then

     write(40,*) 0

  else 

     write(40,*) 1
     write(40,*) "nu1 s1 nu2 s2   k    psi(nu1,s1,nu2,s2,k)"
     
     do nu1 = 0, nb-1
        do is1 = 0, 1
           do nu2 = 0,nb-1
              do is2 = 0, 1

                 do k = 0, nl-1

                    write(40,109) nu1, is1, nu2, is2, k, &
                         real(psi(2*nu1+is1,2*nu2+is2,k)), &
                         imag(psi(2*nu1+is1,2*nu2+is2,k))

                 enddo
              enddo
           enddo

        enddo
     enddo

  endif

  write(40,*)  
  write(40,*) "--Spin-orbit matrix -------"
  write(40,*) "Print flag (0 not printed, 1 printed)"
  write(40,*) 1
  write(40,*) "nu1 is1  nu2 is2 h_so"
  do nu1 = 0, nb-1
     do is1 = 0, 1
        do nu2 = 0, nb-1
           do is2 = 0, 1
              write(40, 119) nu1, is1, nu2, is2, &
                   real(h_so(2*nu1+is1,2*nu2+is2)), &
                   imag(h_so(2*nu1+is1,2*nu2+is2))
           enddo
        enddo
     enddo
  enddo

105 format(i3,2x,i1,2x,i3,2x,i1,2x,'(',d18.10,', ',d18.10,')')
109 format(i3,2x,i1,2x,i3,2x,i1,2x,i7,2x, '(',d18.10,', ',d18.10,')')
119 format(i3,2x,i1,2x,i3,2x,i1,3x, '(',d18.10,', ',d18.10,')')

  write(40,*)
  write(40,*) "------- First order sigma ------"
  write(40,*) "nu1 a1 nu2 a2 sigma1(i, j) "

  do nu1 = 0, nb-1
     do ia1 = 0, 3
        do nu2 = 0, nb-1
           do ia2 = 0, 3
              
              write(40,105) nu1,ia1, nu2, ia2, & 
                   real(sigma1(nu1*4+ia1, nu2*4+ia2)), &
                   imag(sigma1(nu1*4+ia1, nu2*4+ia2))

           enddo
        enddo

     enddo
  enddo

  write(40,*)

200 format('------------------- [',i3,',',i3,',',i3,'] hopping', &
         '------------------------')

300 format(i3,',',i3,'  ','(',D16.9,',',D16.9,')')

  return
end subroutine sigma1_out




