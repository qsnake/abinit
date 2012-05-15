!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp7cg
!! NAME
!! psp7cg
!!
!! FUNCTION
!! Compute sine transform to transform from n(r) to n(q).
!! Computes integrals on (generalized) grid using corrected trapezoidal integration.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt.
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  radmesh <type(pawrad_type)>=data containing radial grid informations
!!  nr(radmesh%mesh_size)=n(r) on radial grid.
!!
!! OUTPUT
!!  dnqdq0= 1/q dn(q)/dq for q=0
!!{{\\ \begin{equation}
!!  nq(mqgrid)= n(q)
!!            = 4\pi\int[(\frac{\sin(2\pi q r)}{2\pi q r})(r^2 n(r))dr].
!!\end{equation} }}
!!  yp1,ypn=derivatives of n(q) wrt q at q=0 and q=qmax (needed for spline fitter).
!!
!! PARENTS
!!      psp7in
!!
!! CHILDREN
!!      simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp7cg(dnqdq0,mqgrid,qgrid,nq,radmesh,nr,yp1,ypn)

 use m_profiling

 use defs_basis
 use defs_datatypes

 use m_radmesh,   only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp7cg'
!End of the abilint section

 implicit none

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: mqgrid
 real(dp),intent(out) :: dnqdq0,yp1,ypn
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: nr(radmesh%mesh_size),qgrid(mqgrid)
 real(dp),intent(out) :: nq(mqgrid)

!Local variables-------------------------------
!scalars
 integer :: iq,ir
 real(dp) :: aexp,arg,bexp,dn,r0tor1,r1torm,rm,rmtoin
 logical :: begin_r0
!arrays
 real(dp),allocatable :: ff(:),rnr(:)

! *************************************************************************

 ABI_ALLOCATE(ff,(radmesh%mesh_size))
 ABI_ALLOCATE(rnr,(radmesh%mesh_size))
 ff=zero;rnr=zero

 do ir=1,radmesh%mesh_size
   rnr(ir)=radmesh%rad(ir)*nr(ir)
 end do

!Is mesh beginning with r=0 ?
 begin_r0=(radmesh%rad(1)<1.d-20)

!Adjustment of an exponentional at r_max (n_exp(r)=aexp*Exp[-bexp*r])
 rm=radmesh%rad(radmesh%mesh_size)
 dn=one/(12._dp*radmesh%stepint*radmesh%radfact(radmesh%mesh_size)) &
& *( 3._dp*nr(radmesh%mesh_size-4) &
& -16._dp*nr(radmesh%mesh_size-3) &
& +36._dp*nr(radmesh%mesh_size-2) &
& -48._dp*nr(radmesh%mesh_size-1) &
& +25._dp*nr(radmesh%mesh_size))
 if (dn<0._dp.or. &
& abs(radmesh%rad(radmesh%mesh_size)*nr(radmesh%mesh_size))>1.d-20) then
   bexp=-dn/nr(radmesh%mesh_size)
   aexp=nr(radmesh%mesh_size)*exp(bexp*rm)
 else
   bexp=0.001_dp
   aexp=zero
 end if

!===========================================
!=== Compute n(q) for q=0 separately
!===========================================

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero
 if (.not.begin_r0) r0tor1=(rnr(1)*radmesh%rad(1)**2)/3.d0

!Integral from r1 to rmax
 do ir=1,radmesh%mesh_size
   if (abs(rnr(ir))>1.d-20) ff(ir)=rnr(ir)*radmesh%rad(ir)
 end do
 call simp_gen(r1torm,ff,radmesh)

!Integral from rmax to infinity
!This part is approximated using an exponential density aexp*Exp[-bexp*r]
!(formulae obtained with mathematica)
 rmtoin=aexp*exp(-bexp*rm)/bexp**3*(two+two*bexp*rm+bexp*bexp*rm*rm)

!Some of the three parts
 nq(1)=four_pi*(r0tor1+r1torm+rmtoin)

!===========================================
!=== Compute n(q) for other q''s
!===========================================

!Loop over q values
 do iq=2,mqgrid
   arg=two_pi*qgrid(iq)

!  Integral from 0 to r1 (only if r1<>0)
   r0tor1=zero;if (.not.begin_r0) &
&   r0tor1=nr(1)*(sin(arg*radmesh%rad(1))/arg/arg&
&   -radmesh%rad(1)*cos(arg*radmesh%rad(1))/arg)

!  Integral from r1 to rmax
   do ir=1,radmesh%mesh_size
     if (abs(rnr(ir))>1.d-20) ff(ir)=sin(arg*radmesh%rad(ir))*rnr(ir)
   end do
   call simp_gen(r1torm,ff,radmesh)

!  Integral from rmax to infinity
!  This part is approximated using an exponential density aexp*Exp[-bexp*r]
!  (formulae obtained with mathematica)
   rmtoin=aexp*exp(-bexp*rm)/(arg**2+bexp**2)**2 &
&   *(arg*(two*bexp+arg**2*rm+bexp**2*rm)*cos(arg*rm) &
&   +(arg**2*(bexp*rm-one)+bexp**2*(bexp*rm+one))*sin(arg*rm))

!  Store q^2 v(q)
   nq(iq)=two/qgrid(iq)*(r0tor1+r1torm+rmtoin)
 end do

!===========================================
!=== Compute derivatives of n(q)
!=== at ends of interval
!===========================================

!yp(0)=zero
 yp1=zero

!yp(qmax)=$ 2\int_0^\infty[(-\sin(2\pi qmax r)+(2\pi qmax r)*\cos(2\pi qmax r) r n(r) dr]$
 arg=two_pi*qgrid(mqgrid)

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero;if (.not.begin_r0) &
& r0tor1=two_pi*nr(1)*(3.d0*radmesh%rad(1)/arg /arg*cos(arg*radmesh%rad(1))+ &
& (radmesh%rad(1)**2/arg-3.0d0/arg**3)*sin(arg*radmesh%rad(1)))

!Integral from r1 to rmax
 do ir=1,radmesh%mesh_size
   if (abs(rnr(ir))>1.d-20) ff(ir)=(two_pi*radmesh%rad(ir)*cos(arg*radmesh%rad(ir)) &
&   - sin(arg*radmesh%rad(ir))/qgrid(mqgrid)) *rnr(ir)
 end do
 call simp_gen(r1torm,ff,radmesh)

!Integral from rmax to infinity
!This part is approximated using an exponential density aexp*Exp[-bexp*r]
!(formulae obtained with mathematica)
 rmtoin=-one/(qgrid(mqgrid)*(arg**2+bexp**2)**3) &
& *aexp*exp(-bexp*rm) &
& *((arg**5*rm-two_pi*arg**4*qgrid(mqgrid)*rm*(bexp*rm-two) &
& +two*arg**3*bexp*(bexp*rm+one)+arg*bexp**3*(bexp*rm+two) &
& -four_pi*arg**2*bexp*qgrid(mqgrid)*(bexp**2*rm**2-three) &
& -two_pi*bexp**3*qgrid(mqgrid)*(bexp**2*rm**2+two*bexp*rm+two))*cos(arg*rm) &
& +(two*arg**2*bexp**3*rm+two_pi*arg**5*qgrid(mqgrid)*rm**2 &
& +arg**4*(bexp*rm-one)+bexp**4*(bexp*rm+one) &
& +four_pi*arg**3*qgrid(mqgrid)*(bexp**2*rm**2+two*bexp*rm-one) &
& +two_pi*arg*bexp**2*qgrid(mqgrid)*(bexp**2*rm**2+four*bexp*rm+six))*sin(arg*rm))

!Some of the three parts
 ypn=two/qgrid(mqgrid)*(r0tor1+r1torm+rmtoin)

!===========================================
!=== Compute 1/q dn(q)/dq at q=0
!===========================================

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero
 if (.not.begin_r0) r0tor1=(rnr(1)*radmesh%rad(1)**4)/5.d0

!Integral from r1 to rmax
 do ir=1,radmesh%mesh_size
   if (abs(rnr(ir))>1.d-20) ff(ir)=rnr(ir)*radmesh%rad(ir)**3
 end do
 call simp_gen(r1torm,ff,radmesh)

!Integral from rmax to infinity
!This part is approximated using an exponential density aexp*Exp[-bexp*r]
!(formulae obtained with mathematica)
 rmtoin=aexp*exp(-bexp*rm)/bexp**5 &
& *(24._dp+24._dp*bexp*rm+12._dp*bexp**2*rm**2+four*bexp**3*rm**3+bexp**4*rm**4)

!Some of the three parts
 dnqdq0=-(2.d0/3.d0)*two_pi**3*(r0tor1+r1torm+rmtoin)

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(rnr)

end subroutine psp7cg
!!***
