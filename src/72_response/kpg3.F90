!{\src2tex{textfont=tt}}
!!****f* ABINIT/kpg3
!! NAME
!! kpg3
!!
!! FUNCTION
!! Compute elements of the derivative $d/dk$ of
!! the kinetic energy operator in reciprocal
!! space at given k point
!! $ E_kin = \frac{d}{dk(\textrm{idir})} ( (1/2) (2 \pi)^2 (k+G)^2 ) ) $
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass=effective mass for electrons (1. in common case)
!!  gmet(3,3) = reciprocal lattice metric tensor (Bohr**-2)
!!  idir      = gives the direction along which the derivative has to be taken
!!  kg(3,npw) = integer coordinates of planewaves in basis sphere.
!!  kpt(3)    = reduced coordinates of k point
!!  npw       = number of plane waves at kpt.
!!
!! OUTPUT
!!  dkinpw(npw)=d/dk(idir) ( (1/2)*(2 pi)**2 * (k+G)**2 )
!!
!!
!! NOTES
!!  (miscellaneous note, warning, etc., if any)
!!
!! TODO
!!  Should merge with mkkin.f , as well as a part of strnps.f
!!
!! PARENTS
!!      nstpaw3,rhofermi3,vtorho3
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine kpg3(dkinpw,ecut,ecutsm,effmass,gmet,idir,kg,kpt,npw)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kpg3'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: idir,npw
 real(dp),intent(in) :: ecut,ecutsm,effmass
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gmet(3,3),kpt(3)
 real(dp),intent(out) :: dkinpw(npw)

!Local variables -------------------------
!scalars
 integer :: ig
 real(dp),parameter :: break_symm=1.0d-11
 real(dp) :: dfsm,dkinetic,dkpg2,ecutsm_inv,fsm,gpk1,gpk2,gpk3,htpisq
! real(dp) :: d2fsm ! used in commented out section below
 real(dp) :: kpg2,xx
!arrays
 real(dp) :: gmet_break(3,3)

! *********************************************************************

!htpisq is (1/2) (2 Pi) **2:
 htpisq=0.5_dp*(two_pi)**2

 ecutsm_inv=0.0_dp
 if(ecutsm>1.0d-20)ecutsm_inv=1/ecutsm

 gmet_break(:,:)=gmet(:,:)
 gmet_break(1,1)=(1.0_dp+break_symm)*gmet(1,1)
 gmet_break(3,3)=(1.0_dp-break_symm)*gmet(3,3)

!$OMP PARALLEL DO PRIVATE(dkinetic,dkpg2,gpk1,gpk2,gpk3,ig,kpg2,xx) &
!$OMP SHARED(dkinpw,ecut,ecutsm,ecutsm_inv) &
!$OMP SHARED(gmet_break,htpisq,idir,kg,kpt,npw)
 do ig=1,npw
   gpk1=dble(kg(1,ig))+kpt(1)
   gpk2=dble(kg(2,ig))+kpt(2)
   gpk3=dble(kg(3,ig))+kpt(3)
   kpg2=htpisq*&
&   ( gmet_break(1,1)*gpk1**2+         &
&   gmet_break(2,2)*gpk2**2+         &
&   gmet_break(3,3)*gpk3**2          &
&   +2.0_dp*(gpk1*gmet_break(1,2)*gpk2+  &
&   gpk1*gmet_break(1,3)*gpk3+  &
&   gpk2*gmet_break(2,3)*gpk3 )  )
   dkpg2=htpisq*2.0_dp*&
&   (gmet_break(idir,1)*gpk1+gmet_break(idir,2)*gpk2+gmet_break(idir,3)*gpk3)
   dkinetic=dkpg2
   if(kpg2>ecut-ecutsm)then
     if(kpg2>ecut-tol12)then
!      The wavefunction has been filtered : no derivative
       dkinetic=0.0_dp
     else
       xx=(ecut-kpg2)*ecutsm_inv
!      This kinetic cutoff smoothing function and its xx derivatives
!      were produced with Mathematica and the fortran code has been
!      numerically checked against Mathematica.
       fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
       dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
!      d2fsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*&
!      &                        (-144+45*xx))))))*fsm**3
       dkinetic=dkpg2*(fsm-ecutsm_inv*kpg2*dfsm)
     end if
   end if
   dkinpw(ig)=dkinetic/effmass
 end do
!$OMP END PARALLEL DO

end subroutine kpg3
!!***
