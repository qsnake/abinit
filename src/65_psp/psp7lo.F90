!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp7lo
!! NAME
!! psp7lo
!!
!! FUNCTION
!! Compute sine transform to transform from V(r) to q^2 V(q).
!! Computes integrals on (generalized) grid using corrected trapezoidal integration.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  radmesh <type(pawrad_type)>=data containing radial grid informations
!!  vloc(radmesh%mesh_size)=V(r) on radial grid.
!!  zion=nominal valence charge of atom.
!!
!! OUTPUT
!!  epsatm=$ 4\pi\int[r^2 (V(r)+\frac{Zv}{r}dr]$.
!!{{\\ \begin{equation}
!!  q2vq(mqgrid)
!!   =q^2 V(q)
!!   = -\frac{Zv}{\pi}
!!     + q^2 4\pi\int[(\frac{\sin(2\pi q r)}{2\pi q r})(r^2 V(r)+r Zv)dr].
!!\end{equation} }}
!!  yp1,ypn=derivatives of q^2 V(q) wrt q at q=0 and q=qmax (needed for spline fitter).
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

subroutine psp7lo(epsatm,mqgrid,qgrid,q2vq,radmesh,vloc,yp1,ypn,zion)

 use m_profiling

 use defs_basis
 use defs_datatypes

 use m_radmesh,  only : ifromr, simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp7lo'
!End of the abilint section

 implicit none

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: mqgrid
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm,yp1,ypn
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: qgrid(mqgrid),vloc(radmesh%mesh_size)
 real(dp),intent(out) :: q2vq(mqgrid)

!Local variables ------------------------------
!scalars
 integer :: iq,ir,irmax
 real(dp) :: arg,r0tor1,r1torm,rmtoin
 logical :: begin_r0
!arrays
 real(dp),allocatable :: ff(:),rvpz(:)

!************************************************************************

 irmax=ifromr(radmesh,min(20._dp,radmesh%rmax))

!Particular case of a zero potential
 if (maxval(abs(vloc(1:irmax)))<=1.e-20_dp) then
   q2vq=zero;yp1=zero;ypn=zero;epsatm=zero
   return
 end if

 ABI_ALLOCATE(ff,(radmesh%mesh_size))
 ABI_ALLOCATE(rvpz,(radmesh%mesh_size))
 ff=zero;rvpz=zero

!Is mesh beginning with r=0 ?
 begin_r0=(radmesh%rad(1)<1.e-20_dp)

!Store r.V+Z
 do ir=1,irmax
   rvpz(ir)=radmesh%rad(ir)*vloc(ir)+zion
 end do

!===========================================
!=== Compute q^2 v(q) for q=0 separately
!===========================================

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero;if (.not.begin_r0) &
& r0tor1=(zion*0.5_dp+radmesh%rad(1)*vloc(1)/3._dp)*radmesh%rad(1)**2

!Integral from r1 to rmax
 do ir=1,irmax
   if (abs(rvpz(ir))>1.e-20_dp) then
     ff(ir)=radmesh%rad(ir)*rvpz(ir)
   end if
 end do
 call simp_gen(r1torm,ff,radmesh)

!Integral from rmax to infinity
!This part is neglected... might be improved.
 rmtoin=zero

!Some of the three parts
 epsatm=four_pi*(r0tor1+r1torm+rmtoin)

 q2vq(1)=-zion/pi

!===========================================
!=== Compute q^2 v(q) for other q''s
!===========================================

!Loop over q values
 do iq=2,mqgrid
   arg=two_pi*qgrid(iq)

!  Integral from 0 to r1 (only if r1<>0)
   r0tor1=zero;if (.not.begin_r0) &
&   r0tor1=( vloc(1)/arg*sin(arg*radmesh%rad(1)) &
&   -rvpz(1)    *cos(arg*radmesh%rad(1)) +zion )/pi

!  Integral from r1 to rmax
   do ir=1,irmax
     if (abs(rvpz(ir))>1.e-20_dp) ff(ir)=sin(arg*radmesh%rad(ir))*rvpz(ir)
   end do
   call simp_gen(r1torm,ff,radmesh)

!  Integral from rmax to infinity
!  This part is neglected... might be improved.
   rmtoin=zero

!  Store q^2 v(q)
   q2vq(iq)=-zion/pi + two*qgrid(iq)*(r0tor1+r1torm+rmtoin)
 end do

!===========================================
!=== Compute derivatives of q^2 v(q)
!=== at ends of interval
!===========================================

!yp(0)=zero
 yp1=zero

!yp(qmax)=$ 2\int_0^\infty[(\sin(2\pi qmax r)+(2\pi qmax r)*\cos(2\pi qmax r)(r V(r)+Z) dr]$
 arg=two_pi*qgrid(mqgrid)

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero;if (.not.begin_r0) &
& r0tor1=zion*radmesh%rad(1)                  *sin(arg*radmesh%rad(1)) &
& +three*radmesh%rad(1)*vloc(1)/arg         *cos(arg*radmesh%rad(1)) &
& +(radmesh%rad(1)**2-one/arg**2)*vloc(1)*sin(arg*radmesh%rad(1))

!Integral from r1 to rmax
 do ir=1,irmax
   if (abs(rvpz(ir))>1.e-20_dp) ff(ir)=( arg*radmesh%rad(ir)*cos(arg*radmesh%rad(ir)) &
&   +                    sin(arg*radmesh%rad(ir))) *rvpz(ir)
 end do
 call simp_gen(r1torm,ff,radmesh)

!Integral from rmax to infinity
!This part is neglected... might be improved.
 rmtoin=zero

!Some of the three parts
 ypn=two*(r0tor1+r1torm+rmtoin)

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(rvpz)

end subroutine psp7lo
!!***
