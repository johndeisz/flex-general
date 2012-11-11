#include "../convert.F90"

subroutine pade(nmax, a, z, p, q)

  USE CONSTANTS
  IMPLICIT NONE

  INTEGER nmax
  COMPLEX a(0:m/2-1)
  COMPLEX z(0:m/2-1)
  COMPLEX p(0:m/2-1)
  COMPLEX q(0:m/2-1)

  INTEGER row, col
  INTEGER l

  COMPLEX mat(0:2*nmax+1,0:2*nmax+1)
  COMPLEX c_vector(0:2*nmax+1)

  INTEGER lda,ldb, ipiv(2*nmax+2),info

  do row = 0, 2 * nmax + 1

     do col = 0, nmax
        mat(row,col) = a(row) * z(row)**col 
        mat(row,col + nmax + 1) = - z(row)**col 
     enddo
  
     c_vector(row) = -a(row) * z(row)**(nmax + 1)

  enddo

  lda = 2*nmax + 2
  ldb = 2*nmax + 2
  call cgetrf(2*nmax+2,2*nmax+2,mat,lda,ipiv,info)
  call cgetrs('N',2*nmax+2,1,mat,lda,ipiv,c_vector,ldb,info)

  do l = 0, nmax
     q(l) = c_vector(l)
     p(l) = c_vector(l+nmax+1)
  enddo

  return
end subroutine pade
