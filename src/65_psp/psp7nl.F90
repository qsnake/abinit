!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp7nl
!! NAME
!! psp7nl
!!
!! FUNCTION
!! Make paw projector form factors f_l(q) for each l
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  indlmn(6,lmnmax)= array giving l,m,n,lm,ln,s for i=lmn
!!  lmnmax=max number of (l,m,n) components
!!  lnmax=max number of (l,n) components
!!  mqgrid=number of grid points for q grid
!!  qgrid(mqgrid)=values at which form factors are returned
!!  radmesh <type(pawrad_type)>=data containing radial grid informations
!!  wfll(mmax,lnmax)=paw projector on radial grid
!!
!! OUTPUT
!!  ffspl(mqgrid,2,lnmax)= form factor f_l(q) and second derivative
!!
!! NOTES
!!  u_l(r) is the paw projector (input as wfll);
!!  j_l(q) is a spherical Bessel function;
!!  f_l(q) = $ \int_0^{rmax}[j_l(2\pi q r) u_l(r)  r dr]$
!!
!! PARENTS
!!      pawinit,psp7in
!!
!! CHILDREN
!!      copymesh,jbessel_4spline,leave_new,simp_gen,spline,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp7nl(ffspl,indlmn,lmnmax,lnmax,mqgrid,qgrid,radmesh,wfll)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_splines

 use m_special_funcs, only : jbessel_4spline
 use m_radmesh,       only : simp_gen, copymesh, compmesh

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp7nl'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmnmax,lnmax,mqgrid
 type(pawrad_type),intent(in) :: radmesh
!arrays
 integer,intent(in) :: indlmn(6,lmnmax)
 real(dp),intent(in) :: qgrid(mqgrid),wfll(radmesh%mesh_size,lnmax)
 real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)

!Local variables-------------------------------
!scalars
 integer :: ilmn,iln,iln0,iq,ir,ll,meshsz,mmax
 real(dp),parameter :: eps=tol14**4,TOLJ=0.001_dp
 real(dp) :: arg,argn,bes
 real(dp) :: besp,qr
 real(dp) :: yp1,ypn
 character(len=500) :: message
 type(pawrad_type) :: tmpmesh
!arrays
 real(dp),allocatable :: ff(:),gg(:),rr(:),rr2(:),rr2wf(:),rrwf(:),work(:)

!*************************************************************************

!Is mesh beginning with r=0 ?
 if (radmesh%rad(1)>tol10) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' psp7nl: BUG -',ch10,&
&   '  Radial mesh cannot begin with r<>0 !'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!Init. temporary arrays and variables
 call copymesh(radmesh,tmpmesh)
 meshsz=tmpmesh%mesh_size;mmax=meshsz
 ABI_ALLOCATE(ff,(meshsz))
 ABI_ALLOCATE(gg,(meshsz))
 ABI_ALLOCATE(rr,(meshsz))
 ABI_ALLOCATE(rr2,(meshsz))
 ABI_ALLOCATE(rrwf,(meshsz))
 ABI_ALLOCATE(rr2wf,(meshsz))
 ABI_ALLOCATE(work,(mqgrid))
 rr(:) =tmpmesh%rad(:)
 rr2(:)=two_pi*rr(:)*rr(:)
 argn=two_pi*qgrid(mqgrid)

!Loop on (l,n) projectors
 iln0=0
 do ilmn=1,lmnmax
   iln=indlmn(5,ilmn)
   if(iln>iln0) then
     iln0=iln;ll=indlmn(1,ilmn)

     ir=meshsz
     do while (abs(wfll(ir,iln))<eps)
       ir=ir-1
     end do
     ir=min(ir+1,meshsz)
     if (ir/=mmax) then
       mmax=ir;call compmesh(tmpmesh,rr(mmax))
     end if

     rrwf(:) =rr (:)*wfll(:,iln)
     rr2wf(:)=rr2(:)*wfll(:,iln)

!    1-Compute f_l(0<q<qmax)
     if (mqgrid>2) then
       do iq=2,mqgrid-1
         arg=two_pi*qgrid(iq)
         do ir=1,mmax
           qr=arg*rr(ir)
           call jbessel_4spline(bes,besp,ll,0,qr,TOLJ)
           ff(ir)=bes*rrwf(ir)
         end do
         call simp_gen(ffspl(iq,1,iln),ff,tmpmesh)
       end do
     end if

!    2-Compute f_l(q=0) and first derivative
     ffspl(1,1,iln)=zero;yp1=zero
     if (ll==0) call simp_gen(ffspl(1,1,iln),rrwf,tmpmesh)
     if (ll==1) then
       call simp_gen(yp1,rr2wf,tmpmesh)
       yp1=yp1*third
     end if

!    3-Compute f_l(q=qmax) and first derivative
     if (mqgrid>1) then
!      if (ll==0.or.ll==1) then
       do ir=1,mmax
         qr=argn*rr(ir)
         call jbessel_4spline(bes,besp,ll,1,qr,TOLJ)
         ff(ir)=bes*rrwf(ir) 
         gg(ir)=besp*rr2wf(ir)
       end do
!      else
!      do ir=1,mmax
!      qr=argn*rr(ir)
!      call jbessel(bes,besp,bespp,ll,1,qr)
!      ff(ir)=bes*rrwf(ir)
!      gg(ir)=besp*rr2wf(ir)
!      end do
!      end if
       call simp_gen(ffspl(mqgrid,1,iln),ff,tmpmesh)
       call simp_gen(ypn,gg,tmpmesh)
     else
       ypn=yp1
     end if

!    4-Compute second derivative of f_l(q)
     call spline(qgrid,ffspl(:,1,iln),mqgrid,yp1,ypn,ffspl(:,2,iln))

!    End loop on (l,n) projectors
   end if
 end do

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(gg)
 ABI_DEALLOCATE(rr)
 ABI_DEALLOCATE(rr2)
 ABI_DEALLOCATE(rrwf)
 ABI_DEALLOCATE(rr2wf)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(tmpmesh%rad)
 ABI_DEALLOCATE(tmpmesh%radfact)
 ABI_DEALLOCATE(tmpmesh%simfact)

end subroutine psp7nl
!!***
