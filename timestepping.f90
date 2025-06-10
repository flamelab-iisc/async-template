subroutine time_stepping()
implicit none

do stage=1,rk_stages
#ifdef ASYNC
if(stage.eq.1) then
call buffer_extrapolate()
else
call buffer_update(stage)
endif
#endif
call flux_calc()
call rhs_calc(stage)
call rk_update(stage)
#ifdef SYNC
call communicate()
#endif
call set_state()
enddo

#ifdef ASYNC
!------------------------------------------------------------------------------
!---------------------SEE LINES 15-24 OF ALGORITHM 2: CAA----------------------
!------------------------------------------------------------------------------
do nbr=1,tot_nbrs
delay(nbr) = delay(nbr) + 1
enddo
if((mod(iter,max_delay).ge.0).and.(mod(iter,max_delay).le.c_time-1)) then
call copy_aux_send()
endif
if(mod(iter,max_delay).eq.c_time) then
call communicate_aux()
do nbr=1,tot_nbrs
delay(nbr) = 0
enddo
endif
#endif

end subroutine

subroutine rhs_calc(stage)
implicit none
#ifdef ASYNC
if(stage.eq.1) then
call dphidx_at(f,fbuff_at,dfdx)
call dphidy_at(g,gbuff_at,dgdy)
call dphidz_at(h,hbuff_at,dhdz)
else
call dphidx(f,fbuff,dfdx)
call dphidy(g,gbuff,dgdy)
call dphidz(h,hbuff,dhdz)
endif
#elif SYNC
call dphidx(f,fbuff,dfdx)
call dphidy(g,gbuff,dgdy)
call dphidz(h,hbuff,dhdz)
#endif
rhs(:,:,:,:,:) = dfdx(:,:,:,:,:) + dgdy(:,:,:,:,:) + dhdz(:,:,:,:,:)
enddo

subroutine rk_update(stage)
implicit none
do blk=1,nblocks
do eq=1,neq
do k=1,nk(blk)
do j=1,nj(blk)
do i=1,ni(blk)
q(i,j,k,eq,blk) = a_rk(stage)*q(i,j,k,eq,blk) + dt*rhs(i,j,k,eq,blk)
u(i,j,k,eq,blk) = u(i,j,k,eq,blk) + b_rk(stage)*q(i,j,k,eq,blk)
enddo
enddo
enddo
enddo
enddo
end subroutine
