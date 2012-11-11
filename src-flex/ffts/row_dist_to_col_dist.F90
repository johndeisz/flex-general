#include "../convert.F90"

subroutine row_dist_to_col_dist(myid, a, b )

  USE CONSTANTS
  IMPLICIT NONE

#ifdef USE_MPI
  include 'mpif.h'
#endif /* USE_MPI */

  INTEGER myid
  COMPLEX b(0:m1,0:ncp1), a(0:mp1,0:nc1)
      
  COMPLEX buffer_send(0:mp*ncp-1)
  COMPLEX buffer_recv(0:mp*ncp-1)

  INTEGER i, pair_sum, partner
  INTEGER p, q, l

  INTEGER ierr
#ifdef USE_MPI
  INTEGER stat(MPI_STATUS_SIZE)
#endif /* USE_MPI */

  do i = 0, np - 1
    
     pair_sum = np - 1 - i
        
     if (myid .le. pair_sum) then
        partner = pair_sum - myid
     else  
        partner = pair_sum + np - myid
     endif

     l = 0

     do p = 0, ncp1
        do q = 0, mp1
	
           buffer_send(l) = a(q, p + partner * ncp)
           l = l+1

        enddo
     enddo

     if (myid .eq. partner) then

        do l = 0, mp*ncp - 1
           buffer_recv(l) = buffer_send(l)
        enddo

#ifdef USE_MPI
          
     else 

        if (myid .lt. partner) then

           call MPI_Send( buffer_send, mp * ncp, MPI_COMPLEX, partner, 0, &
                MPI_COMM_WORLD, ierr)
            
           call MPI_Recv( buffer_recv, mp * ncp, MPI_COMPLEX, partner, 0, &
                MPI_COMM_WORLD, stat, ierr)
         
        else 

           call MPI_Recv( buffer_recv, mp * ncp, MPI_COMPLEX, partner, 0, &
                MPI_COMM_WORLD, stat, ierr)

           call MPI_Send( buffer_send, mp * ncp, MPI_COMPLEX, partner, 0, &
                MPI_COMM_WORLD, ierr)

        endif
#endif /* USE_MPI */
        
     endif

     !     Unpack the received buffer, column after column

     l = 0

     do p = 0, ncp1
        do q = 0, mp1

           b(partner * mp + q,p) = buffer_recv(l)
           l = l + 1

        enddo
     enddo

  enddo

  return
end subroutine row_dist_to_col_dist
