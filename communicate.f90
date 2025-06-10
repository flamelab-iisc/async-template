subroutine communicate()
implicit none
include "mpif.h"
do nbr=1,tot_nbrs
call mpi_irecv(recv_buffer(1,neq,nbr),buff_size(nbr),mpi_double_precision,recv_rank(nbr),mpi_any_tag,mpi_comm_world,mpi_request_recv(nbr),ierr)
enddo
do nbr=1,tot_nbrs
call mpi_wait(mpi_request_recv(nbr),status_recv,ierr)
do pt=1,buff_pts(nbr)
call get_boundary_id(pt,i,j,k,blk)
do eq=1,neq
send_buffer(pt,eq,nbr) = u(i,j,k,eq,blk)
enddo
enddo
enddo
do nbr=1,tot_nbrs
call mpi_isend(send_buffer(1,neq,nbr),buff_size(nbr),mpi_double_precision,send_rank(nbr),mpi_any_tag,mpi_comm_world,mpi_request_send(nbr),ierr)
enddo
do nbr=1,tot_nbrs
call mpi_wait(mpi_request_recv(nbr),status_recv,ierr)
do pt=1,buff_pts(nbr)
call get_buffer_id(pt,i,j,k,blk)
do eq=1,neq
ubuff(i,j,k,eq,blk,nbr) = recv_buffer(pt,eq,nbr)
enddo
enddo
enddo
do nbr=1,tot_nbrs
call mpi_wait(mpi_request_send(nbr),status_send,ierr)
enddo
end subroutine

subroutine copy_aux_send()
implicit none
t_l=mod(iter,max_delay)
do nbr=1,tot_nbrs
do pt=1,buff_pts(nbr)
call get_boundary_id(pt,i,j,k,blk)
do eq=1,neq
send_buffer(pt,eq,nbr,t_l) = u(i,j,k,eq,blk)
enddo
enddo
enddo
end subroutine

subroutine communicate_aux()
implicit none
include "mpif.h"
do nbr=1,tot_nbrs
call mpi_irecv(recv_at_buffer(1,neq,nbr,1),buff_size_at(nbr),mpi_double_precision,recv_rank(nbr),mpi_any_tag,mpi_comm_world,mpi_request_recv(nbr),ierr)
enddo
do nbr=1,tot_nbrs
call mpi_isend(send_at_buffer(1,neq,nbr,1),buff_size_at(nbr),mpi_double_precision,send_rank(nbr),mpi_any_tag,mpi_comm_world,mpi_request_send(nbr),ierr)
enddo
do nbr=1,tot_nbrs
call mpi_wait(mpi_request_recv(nbr),status_recv_at,ierr)
do pt=1,buff_pts(nbr)
call get_buffer_id(pt,i,j,k,blk)
do t_l=1,c_time
do eq=1,neq
ubuff_at(i,j,k,eq,blk,nbr,t_l) = recv_at_buffer(pt,eq,nbr,t_l)
enddo
enddo
enddo
enddo
do nbr=1,tot_nbrs
call mpi_wait(mpi_request_send(nbr),status_send_at,ierr)
enddo
end subroutine
