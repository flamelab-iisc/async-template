subroutine dphidx(phi,phibuff,dphidx)
implicit none

do blk=1,nblocks
do eq=1,neq
do k=1,nk(blk)
do j=1,nj(blk)
do i=1,ni(blk)
if(i.eq.1) then
call get_xbuffer_id(i,blk,ib,blk,nbr)
dphidx(i,j,k,eq,blk) = &
(phi(i+1,j,k,eq,blk)-phibuff(ib,j,k,eq,blk,nbr))/(2.*dx)
elseif(i.eq.ni(blk)) then
call get_xbuffer_id(i,blk,ib,blk,nbr)
dphidx(i,j,k,eq,blk) = &
(phibuff(ib,j,k,eq,blk,nbr)-phi(i-1,j,k,eq,blk))/(2.*dx)
else
dphidx(i,j,k,eq,blk) = (phi(i+1,j,k,eq,blk)-phi(i-1,j,k,eq,blk))/(2.*dx)
endif
enddo
enddo
enddo
enddo
enddo

end subroutine

subroutine dphidy(phi,phibuff,dphidy)
implicit none

do blk=1,nblocks
do eq=1,neq
do k=1,nk(blk)
do j=1,nj(blk)
do i=1,ni(blk)
if(j.eq.1) then
call get_ybuffer_id(j,blk,jb,blk,nbr)
dphidy(i,j,k,eq,blk) = &
(phi(i,j+1,k,eq,blk)-phibuff(i,jb,k,eq,blk,nbr))/(2.*dy)
elseif(j.eq.nj(blk)) then
call get_ybuffer_id(j,blk,jb,blk,nbr)
dphidy(i,j,k,eq,blk) = &
(phibuff(i,jb,k,eq,blk,nbr)-phi(i,j-1,k,eq,blk))/(2.*dy)
else
dphidy(i,j,k,eq,blk) = (phi(i,j+1,k,eq,blk)-phi(i,j-1,k,eq,blk))/(2.*dy)
endif
enddo
enddo
enddo
enddo
enddo

end subroutine

subroutine dphidz(phi,phibuff,dphidz)
implicit none

do blk=1,nblocks
do eq=1,neq
do k=1,nk(blk)
do j=1,nj(blk)
do i=1,ni(blk)
if(k.eq.1) then
call get_zbuffer_id(k,blk,kb,blk,nbr)
dphidz(i,j,k,eq,blk) = &
(phi(i,j,k+1,eq,blk)-phibuff(i,j,kb,eq,blk,nbr))/(2.*dz)
elseif(k.eq.nk(blk)) then
call get_zbuffer_id(k,blk,kb,blk,nbr)
dphidz(i,j,k,eq,blk) = &
(phibuff(i,j,kb,eq,blk,nbr)-phi(i,j,k-1,eq,blk))/(2.*dz)
else
dphidz(i,j,k,eq,blk) = (phi(i,j,k+1,eq,blk)-phi(i,j,k-1,eq,blk))/(2.*dz)
endif
enddo
enddo
enddo
enddo
enddo

end subroutine


subroutine buffer_extrapolate()
!------------------------------------------------------------------------------
!----------------SEE EQ.14: BUFFER EXTRAPOLATION SCHEME------------------------
!------------------------------------------------------------------------------
implicit none
do nbr=1,tot_nbrs
do pt=1,buff_pts(nbr)
call get_buffer_id(pt,i,j,k,blk)
do eq=1,neq
coeff_nmk   = 0.5*(delay(nbr)**2 + 3.*delay(nbr) + 2.)
coeff_nmkm1 = -(delay(nbr)**2 + 2.*delay(nbr))
coeff_nmkm2 = 0.5*(delay(nbr)**2 + delay(nbr))
ubuff(i,j,k,eq,blk,nbr) = coeff_nmk*ubuff_at(i,j,k,eq,blk,nbr,1) + coeff_nmkm1*ubuff_at(i,j,k,eq,blk,nbr,2) + coeff_nmkm2*ubuff_at(i,j,k,eq,blk,nbr,2)
enddo
enddo
enddo
end subroutine

subroutine buffer_update(stage)
!------------------------------------------------------------------------------
!----------------SEE APPENDIX C: BUFFER UPDATE SCHEMES-------------------------
!------------------------------------------------------------------------------
implicit none
select case(stage)
case(2)
do nbr=1,tot_nbrs
do pt=1,buff_pts(nbr)
call get_buffer_id(pt,i,j,k,blk)
do eq=1,neq
coeff_nmk   = (2.*delay(nbr) + 3.)/(2.*dt)
coeff_nmkm1 = -(2.*delay(nbr) + 2.)/dt
coeff_nmkm2 = (2.*delay(nbr) + 1.)/(2.*dt)
dubdt = coeff_nmk*ubuff_at(i,j,k,eq,blk,nbr,1) + coeff_nmkm1*ubuff_at(i,j,k,eq,blk,nbr,2) + coeff_nmkm2*ubuff_at(i,j,k,eq,blk,nbr,2)
qbuff(i,j,k,eq,blk,nbr) = a_rk(stage)*qbuff(i,j,k,eq,blk,nbr) + dt*dubdt
ubuff(i,j,k,eq,blk,nbr) = ubuff(i,j,k,eq,blk,nbr) + b_rk(stage)*qbuff(i,j,k,eq,blk,nbr)
enddo
enddo
enddo
end select
end subroutine

subroutine dphidx_at(phi,phibuff_at,dphidx)
!------------------------------------------------------------------------------
!-----------------SEE EQ.7: AT SCHEMES FOR PE BOUNDARY POINTS------------------
!------------------------------------------------------------------------------
implicit none
do blk=1,nblocks
do eq=1,neq
do k=1,nk(blk)
do j=1,nj(blk)
do i=1,ni(blk)
if(i.eq.1) then
call get_xbuffer_id(i,blk,ib,blk,nbr)
coeff_nmk   = 0.5*(delay(nbr)**2 + 3.*delay(nbr) + 2.)
coeff_nmkm1 = -(delay(nbr)**2 + 2.*delay(nbr))
coeff_nmkm2 = 0.5*(delay(nbr)**2 + delay(nbr))
dphidx(i,j,k,eq,blk) = &
+coeff_nmk  *(phi(i+1,j,k,eq,blk)-phibuff_at(ib,j,k,eq,blk,nbr,1))/(2.*dx)
+coeff_nmkm1*(phi(i+1,j,k,eq,blk)-phibuff_at(ib,j,k,eq,blk,nbr,2))/(2.*dx)
+coeff_nmkm2*(phi(i+1,j,k,eq,blk)-phibuff_at(ib,j,k,eq,blk,nbr,3))/(2.*dx)
elseif(i.eq.ni(blk)) then
call get_xbuffer_id(i,blk,ib,blk,nbr)
coeff_nmk   = 0.5*(delay(nbr)**2 + 3.*delay(nbr) + 2.)
coeff_nmkm1 = -(delay(nbr)**2 + 2.*delay(nbr))
coeff_nmkm2 = 0.5*(delay(nbr)**2 + delay(nbr))
dphidx(i,j,k,eq,blk) = &
+coeff_nmk  *(phibuff_at(ib,j,k,eq,blk,nbr,1)-phi(i-1,j,k,eq,blk))/(2.*dx)
+coeff_nmkm1*(phibuff_at(ib,j,k,eq,blk,nbr,2)-phi(i-1,j,k,eq,blk))/(2.*dx)
+coeff_nmkm2*(phibuff_at(ib,j,k,eq,blk,nbr,3)-phi(i-1,j,k,eq,blk))/(2.*dx)
else
dphidx(i,j,k,eq,blk) = (phi(i+1,j,k,eq,blk)-phi(i-1,j,k,eq,blk))/(2.*dx)
endif
enddo
enddo
enddo
enddo
enddo
end subroutine

subroutine dphidy_at(phi,phibuff_at,dphidy)
!------------------------------------------------------------------------------
!-----------------SEE EQ.7: AT SCHEMES FOR PE BOUNDARY POINTS------------------
!------------------------------------------------------------------------------
implicit none
do blk=1,nblocks
do eq=1,neq
do k=1,nk(blk)
do j=1,nj(blk)
do i=1,ni(blk)
if(j.eq.1) then
call get_ybuffer_id(j,blk,jb,blk,nbr)
coeff_nmk   = 0.5*(delay(nbr)**2 + 3.*delay(nbr) + 2.)
coeff_nmkm1 = -(delay(nbr)**2 + 2.*delay(nbr))
coeff_nmkm2 = 0.5*(delay(nbr)**2 + delay(nbr))
dphidy(i,j,k,eq,blk) = &
+coeff_nmk  *(phi(i,j+1,k,eq,blk)-phibuff_at(i,jb,k,eq,blk,nbr,1))/(2.*dy)
+coeff_nmkm1*(phi(i,j+1,k,eq,blk)-phibuff_at(i,jb,k,eq,blk,nbr,2))/(2.*dy)
+coeff_nmkm2*(phi(i,j+1,k,eq,blk)-phibuff_at(i,jb,k,eq,blk,nbr,3))/(2.*dy)
elseif(j.eq.nj(blk)) then
call get_ybuffer_id(j,blk,jb,blk,nbr)
coeff_nmk   = 0.5*(delay(nbr)**2 + 3.*delay(nbr) + 2.)
coeff_nmkm1 = -(delay(nbr)**2 + 2.*delay(nbr))
coeff_nmkm2 = 0.5*(delay(nbr)**2 + delay(nbr))
dphidy(i,j,k,eq,blk) = &
+coeff_nmk  *(phibuff_at(i,jb,k,eq,blk,nbr,1)-phi(i,j-1,k,eq,blk))/(2.*dy)
+coeff_nmkm1*(phibuff_at(i,jb,k,eq,blk,nbr,2)-phi(i,j-1,k,eq,blk))/(2.*dy)
+coeff_nmkm2*(phibuff_at(i,jb,k,eq,blk,nbr,3)-phi(i,j-1,k,eq,blk))/(2.*dy)
else
dphidy(i,j,k,eq,blk) = (phi(i,j+1,k,eq,blk)-phi(i,j-1,k,eq,blk))/(2.*dy)
endif
enddo
enddo
enddo
enddo
enddo
end subroutine

subroutine dphidz_at(phi,phibuff_at,dphidz)
!------------------------------------------------------------------------------
!-----------------SEE EQ.7: AT SCHEMES FOR PE BOUNDARY POINTS------------------
!------------------------------------------------------------------------------
implicit none
do blk=1,nblocks
do eq=1,neq
do k=1,nk(blk)
do j=1,nj(blk)
do i=1,ni(blk)
if(k.eq.1) then
call get_zbuffer_id(k,blk,kb,blk,nbr)
coeff_nmk   = 0.5*(delay(nbr)**2 + 3.*delay(nbr) + 2.)
coeff_nmkm1 = -(delay(nbr)**2 + 2.*delay(nbr))
coeff_nmkm2 = 0.5*(delay(nbr)**2 + delay(nbr))
dphidz(i,j,k,eq,blk) = &
+coeff_nmk  *(phi(i,j,k+1,eq,blk)-phibuff_at(i,j,kb,eq,blk,nbr,1))/(2.*dz)
+coeff_nmkm1*(phi(i,j,k+1,eq,blk)-phibuff_at(i,j,kb,eq,blk,nbr,2))/(2.*dz)
+coeff_nmkm2*(phi(i,j,k+1,eq,blk)-phibuff_at(i,j,kb,eq,blk,nbr,3))/(2.*dz)
elseif(k.eq.nk(blk)) then
call get_zbuffer_id(k,blk,kb,blk,nbr)
coeff_nmk   = 0.5*(delay(nbr)**2 + 3.*delay(nbr) + 2.)
coeff_nmkm1 = -(delay(nbr)**2 + 2.*delay(nbr))
coeff_nmkm2 = 0.5*(delay(nbr)**2 + delay(nbr))
dphidz(i,j,k,eq,blk) = &
+coeff_nmk  *(phibuff_at(i,j,kb,eq,blk,nbr,1)-phi(i,j,k-1,eq,blk))/(2.*dz)
+coeff_nmkm1*(phibuff_at(i,j,kb,eq,blk,nbr,2)-phi(i,j,k-1,eq,blk))/(2.*dz)
+coeff_nmkm2*(phibuff_at(i,j,kb,eq,blk,nbr,3)-phi(i,j,k-1,eq,blk))/(2.*dz)
else
dphidz(i,j,k,eq,blk) = (phi(i,j,k+1,eq,blk)-phi(i,j,k-1,eq,blk))/(2.*dz)
endif
enddo
enddo
enddo
enddo
enddo
end subroutine
