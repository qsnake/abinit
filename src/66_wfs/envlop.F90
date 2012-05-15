!{\src2tex{textfont=tt}}
!!****f* ABINIT/envlop
!!
!! NAME
!! envlop
!!
!! FUNCTION
!! Multiply random number values in cg by envelope function to lower
!! initial kinetic energy.
!! Envelope  $\left( 1-\left( G/G_{\max }\right) ^2\right) ^{power}$
!! for |G|<= Gmax.
!! Near G=0, little scaling, and goes to zero flatly near Gmax.Loop over perturbations
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cg(2,npw*nband)=initial random number wavefunctions
!! ecut=kinetic energy cutoff in Ha
!! gmet(3,3)=reciprocal space metric (bohr^-2)
!! icgmod=shift to be given to the location of data in cg
!! kg(3,npw)=reduced coordinates of G vectors in basis sphere
!! kpoint(3)=reduced coordinates of k point
!! mcg=maximum second dimension of cg (at least npw*nband*nspinor)
!! nband=number of bands being considered
!! npw=number of planewaves in basis sphere
!! nspinor=number of spinorial components of the wavefunctions
!!
!! OUTPUT
!!  cg(2,mcg)=revised values (not orthonormalized)
!!
!! PARENTS
!!      wfconv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine envlop(cg,ecut,gmet,icgmod,kg,kpoint,mcg,nband,npw,nspinor)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'envlop'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icgmod,mcg,nband,npw,nspinor
 real(dp),intent(in) :: ecut
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gmet(3,3),kpoint(3)
 real(dp),intent(inout) :: cg(2,mcg)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,power=12,re=1
 integer :: i1,i2,i3,ig,igs,ispinor,nn
 real(dp) :: cut,cutoff,fcut,gs,gsq,kpgsqc

! *************************************************************************

!$(k+G)^2$ evaluated using metric and kpoint
 gsq(i1,i2,i3)=gmet(1,1)*(kpoint(1)+dble(i1))**2+&
& gmet(2,2)*(kpoint(2)+dble(i2))**2+&
& gmet(3,3)*(kpoint(3)+dble(i3))**2+&
& 2.0_dp*(gmet(2,1)*(kpoint(2)+dble(i2))*(kpoint(1)+dble(i1))+&
& gmet(3,2)*(kpoint(3)+dble(i3))*(kpoint(2)+dble(i2))+&
& gmet(1,3)*(kpoint(1)+dble(i1))*(kpoint(3)+dble(i3)))

!cutoff function
 fcut(gs)=(1.0_dp-(gs/cutoff))**power

!$(k+G)^2$ cutoff from $(1/2)(2 Pi (k+G))^2 = ecut$
 kpgsqc=ecut/(2.0_dp*pi**2)
 cutoff = kpgsqc

!Run through G vectors in basis
 do ig=1,npw
   gs=gsq(kg(1,ig),kg(2,ig),kg(3,ig))
   if (gs>cutoff) then
     cut=0.0_dp
   else
     cut=fcut(gs)
   end if
!  Run through bands (real and imaginary components)
   do ispinor=1,nspinor
     igs=ig+(ispinor-1)*npw
     do nn=1,nband
       cg(re,igs+(nn-1)*npw*nspinor+icgmod)=cg(re,igs+(nn-1)*npw*nspinor+icgmod)*cut
       cg(im,igs+(nn-1)*npw*nspinor+icgmod)=cg(im,igs+(nn-1)*npw*nspinor+icgmod)*cut
     end do
   end do
 end do

end subroutine envlop
!!***
