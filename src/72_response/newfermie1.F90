!{\src2tex{textfont=tt}}
!!****f* ABINIT/newfermie1
!! NAME
!! newfermie1
!!
!! FUNCTION
!! This routine computes the derivative of the fermi energy wrt
!! the active perturbation for use in evaluating the edocc term
!! and active subspace contribution to the first-order wavefunctions
!! in the case of metals.  This is presently used only for the
!! strain perturbation, and only for Q = 0.
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  fe1fixed=fixed contribution to the first-order Fermi energy
!!  istep=index of the number of steps in the routine scfcv
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftot= total number of FFT grid points
!!  nspden=number of spin-density components
!!  occopt=option for occupancies
!!  prtvol=control print volume and debugging output
!!  rhorfermi(nfft,nspden)=array for fermi-level electron density
!!  ucvol=unit cell volume in bohr**3
!!  vtrial1(cplex*nfft,nspden)=INPUT RF Vtrial(r)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  fermie1=derivative of fermi energy wrt perturbation
!!   at input  : old value
!!   at output : updated value
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      dotprod_vn,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine newfermie1(cplex,fermie1,fe1fixed,istep,&
&  mpi_enreg,nfft,nfftot,nspden,occopt,&
&  prtvol,rhorfermi,ucvol,vtrial1)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'newfermie1'
 use interfaces_14_hidewrite
 use interfaces_53_spacepar
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,istep,nfft,nfftot,nspden,occopt,prtvol
 real(dp),intent(in) :: fe1fixed,ucvol
 real(dp),intent(inout) :: fermie1
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: rhorfermi(nfft,nspden),vtrial1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: option
 real(dp) :: doti,dotr,fermie1_new,fermie1rs
 character(len=500) :: message

! *********************************************************************

 if(occopt>=3 .and. occopt <=8) then

!  The product of the current trial potential and the so-called Fermi level
!  density is integrated to give the local potential contributions to the
!  first-order Fermi level.
   option=1
   call dotprod_vn(cplex,rhorfermi,dotr,doti,mpi_enreg,nfft,nfftot,nspden,option,&
&   vtrial1,ucvol)

!  The fixed contributions consisting of non-local potential and kinetic terms
!  are added
   fermie1_new=fe1fixed+dotr
   fermie1rs=(fermie1-fermie1_new)**2
   fermie1=fermie1_new


   if(prtvol>=10)then
     write(message, '(a,i5,es18.8,es18.8)' ) &
&     ' fermie1, residual squared',istep,fermie1,fermie1rs
     call wrtout(std_out,  message,'COLL')
   end if

 else
   fermie1=zero
 end if

end subroutine newfermie1
!!***
