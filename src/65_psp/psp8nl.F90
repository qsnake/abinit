!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp8nl
!! NAME
!! psp8nl
!!
!! FUNCTION
!! Make Kleinman-Bylander/Bloechl form factors f_ln(q) for each
!!  projector n for each angular momentum l excepting an l corresponding
!!  to the local potential.
!! Note that an arbitrary local potential can be used, so all l from
!!  0 to lmax may be represented.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, FrD, GZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  amesh=grid spacing for uniform (linear) radial grid
!!  lmax=maximum ang momentum for which nonlocal form factor is desired.
!!    lmax <= 2 allowed.
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mmax=number of radial grid points for atomic grid
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=number of grid points for q grid
!!  nproj(mpsang)=number of projection functions for each angular momentum
!!  qgrid(mqgrid)=values at which form factors are returned
!!  rad(mmax)=radial grid values
!!  vpspll(mmax,lnmax)=nonlocal projectors for each (l,n) on linear
!!   radial grid.  Here, these are the  product of the reference
!!   wave functions and (v(l,n)-vloc), calculated in the psp generation
!!   program and normalized so that integral(0,rc(l)) vpsll^2 dr = 1,
!!   which leads to the the usual convention for the energies ekb(l,n)
!!   also calculated in the psp generation program.
!!
!! OUTPUT
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_ln(q) and
!!   second derivative from spline fit for each (l,n).
!!
!! NOTES
!! u_l(r) is reference state wavefunction (input as wf);
!! j_l(q) is a spherical Bessel function;
!! dV_l(r) = vpsp_l(r)-vloc(r) for angular momentum l;
!! f_l(q) = $ \int_0^{rmax}[j_l(2\pi q r) u_l(r) dV_l(r) r dr]/\sqrt{dvms}$
!! where dvms = $\int_0^{rmax} [(u_l(r) dV_l(r))^2 dr]$ is the mean
!! square value of the nonlocal correction for angular momentum l.
!! Xavier Gonze s E_KB = $ dvms/\int_0^{rmax}[(u_l(r))^2 dV_l(r) dr]$.
!! This is the eigenvalue of the Kleinman-Bylander operator and sets
!! the energy scale of the nonlocal psp corrections.
!!
!! PARENTS
!!      psp8in
!!
!! CHILDREN
!!      ctrap,sbf8,spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp8nl(amesh,ffspl,lmax,lnmax,mmax,mpsang,mqgrid,nproj,qgrid,rad,vpspll)

 use m_profiling

 use defs_basis
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp8nl'
 use interfaces_32_util
 use interfaces_65_psp, except_this_one => psp8nl
!End of the abilint section

 implicit none

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: lmax,lnmax,mmax,mpsang,mqgrid
 real(dp),intent(in) :: amesh
!arrays
 integer,intent(in) :: nproj(mpsang)
 real(dp),intent(in) :: qgrid(mqgrid),rad(mmax),vpspll(mmax,lnmax)
 real(dp),intent(out) :: ffspl(mqgrid,2,lnmax)

!Local variables-------------------------------
!Following parameter controls accuracy of Fourier transform based on qmax
!and represents the minimun number of integration points in one period.
!scalars
 integer,parameter :: NPT_IN_2PI=200
 integer :: iln0,ipsang,iq,ir,irmu,irn,jj,mesh_mult,mmax_new,mvpspll
 real(dp) :: amesh_new,arg,c1,c2,c3,c4,dri,qmesh,result,tv,xp,xpm1,xpm2,xpp1
 real(dp) :: yp1,ypn
!arrays
 real(dp) :: sb_out(4)
 real(dp),allocatable :: rad_new(:),vpspll_new(:,:),work(:,:),work2(:)

! *************************************************************************

!Find r mesh spacing necessary for accurate integration at qmax
 amesh_new=2.d0*pi/(NPT_IN_2PI*qgrid(mqgrid))

!Choose submultiple of input mesh
 mesh_mult=int(amesh/amesh_new) + 1
 mmax_new=mesh_mult*(mmax-1)+1
 amesh_new=amesh/dble(mesh_mult)

 ABI_ALLOCATE(rad_new,(mmax_new))
 ABI_ALLOCATE(vpspll_new,(mmax_new,lnmax))

 if(mesh_mult==1) then
   rad_new(:)=rad(:)
 else
!  Set up new radial mesh
   irn=1
   do ir=1,mmax-1
     do irmu=0,mesh_mult-1
       rad_new(irn)=rad(ir)+dble(irmu)*amesh_new
       irn=irn+1
     end do
   end do
   rad_new(mmax_new)=rad(mmax)
 end if

!Interpolate projectors onto new grid if called for
!Cubic polynomial interpolation is used which is consistent
!with the original interpolation of these functions from
!a log grid to the input linear grid.
 dri = one/amesh
 do irn=1,mmax_new
!  index to find bracketing input mesh points
   if(mesh_mult>1) then
     ir = irn/mesh_mult + 1
     ir = max(ir,2)
     ir = min(ir,mmax-2)
!    interpolation coefficients
     xp = dri * (rad_new(irn) - rad(ir))
     xpp1 = xp + one
     xpm1 = xp - one
     xpm2 = xp - two
     c1 = -xp * xpm1 * xpm2 * sixth
     c2 = xpp1 * xpm1 * xpm2 * half
     c3 = - xp * xpp1 * xpm2 * half
     c4 = xp * xpp1 * xpm1 * sixth
!    Now do the interpolation on all projectors for this grid point
     iln0=0
     do ipsang=1,lmax+1
       if(nproj(ipsang)>0) then
         do jj=iln0+1,iln0+nproj(ipsang)
           tv =  c1 * vpspll(ir - 1, jj) &
&           + c2 * vpspll(ir    , jj) &
&           + c3 * vpspll(ir + 1, jj) &
&           + c4 * vpspll(ir + 2, jj)
           if(abs(tv)>tol10) then
             vpspll_new(irn,jj)=tv
             mvpspll=irn
           else
             vpspll_new(irn,jj)=zero
           end if
         end do
         iln0=iln0+nproj(ipsang)
       end if
     end do
   else
!    With no mesh multiplication, just copy projectors
     ir=irn
     iln0=0
     do ipsang=1,lmax+1
       if(nproj(ipsang)>0) then
         do jj=iln0+1,iln0+nproj(ipsang)
           tv = vpspll(ir,jj)
           if(abs(tv)>tol10) then
             vpspll_new(irn,jj)=tv
             mvpspll=irn
           else
             vpspll_new(irn,jj)=zero
           end if
         end do
         iln0=iln0+nproj(ipsang)
       end if
     end do
   end if
 end do

 ABI_ALLOCATE(work,(mvpspll,lnmax))

!Loop over q values
 do iq=1,mqgrid
   arg=2.d0*pi*qgrid(iq)

!  Set up integrands
   do  ir=1,mvpspll
     call sbf8(lmax+1,arg*rad_new(ir),sb_out)
     iln0=0
     do ipsang=1,lmax+1
       if(nproj(ipsang)>0) then
         do jj=iln0+1,iln0+nproj(ipsang)
           work(ir,jj)=sb_out(ipsang)*vpspll_new(ir,jj)*rad_new(ir)
         end do
         iln0=iln0+nproj(ipsang)
       end if
     end do
   end do

!  Do integral from zero to rad_new(mvpspll)
   iln0=0
   do ipsang=1,lmax+1
     if(nproj(ipsang)>0) then
       do jj=iln0+1,iln0+nproj(ipsang)
         call ctrap(mvpspll,work(1,jj),amesh_new,result)
         ffspl(iq,1,jj)=result
       end do
       iln0=iln0+nproj(ipsang)
     end if
   end do

!  End loop over q mesh
 end do


!Fit splines for form factors
 ABI_ALLOCATE(work2,(mqgrid))
 qmesh=qgrid(2)-qgrid(1)
 iln0=0
 do ipsang=1,lmax+1
   if(nproj(ipsang)>0) then
     do jj=iln0+1,iln0+nproj(ipsang)
!      Compute derivatives of form factors at ends of interval
       yp1=(-50.d0*ffspl(1,1,jj)+96.d0*ffspl(2,1,jj)-72.d0*ffspl(3,1,jj)&
&       +32.d0*ffspl(4,1,jj)- 6.d0*ffspl(5,1,jj))/(24.d0*qmesh)
       ypn=(6.d0*ffspl(mqgrid-4,1,jj)-32.d0*ffspl(mqgrid-3,1,jj)&
&       +72.d0*ffspl(mqgrid-2,1,jj)-96.d0*ffspl(mqgrid-1,1,jj)&
&       +50.d0*ffspl(mqgrid,1,jj))/(24.d0*qmesh)

       call spline(qgrid,ffspl(1,1,jj),mqgrid,yp1,ypn,ffspl(1,2,jj))

     end do
     iln0=iln0+nproj(ipsang)
   end if
 end do

 ABI_DEALLOCATE( rad_new)
 ABI_DEALLOCATE(vpspll_new)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(work2)

end subroutine psp8nl
!!***
