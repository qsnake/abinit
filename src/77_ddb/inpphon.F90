!{\src2tex{textfont=tt}}
!!****f* ABINIT/inpphon
!!
!! NAME
!! inpphon
!!
!! FUNCTION
!! Interpolate phonon frequencies and eigenvectors at 1 qpt
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   phon_ds = datastructure containing all needed interatomic force constants etc...
!!      for interpolation
!!   qpt = qpoint coordinate we want to interpolate at
!!
!! OUTPUT
!!   displ_cart = phonon mode displacement vectors in Cartesian coordinates.
!!   pheigval = eigenvalues of modes
!!   pheigvec = eigenvectors of modes
!!   phfrq = frequencies of modes
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      interpolate_gkk,m_gamma,mka2f,mka2f_tr,mkph_linwid,read_gkk
!!
!! CHILDREN
!!      gtdyn9,phfrq3
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inpphon(displ_cart,pheigval,pheigvec,phfrq,phon_ds,qpt)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inpphon'
 use interfaces_72_response
 use interfaces_77_ddb, except_this_one => inpphon
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(phon_type),intent(inout) :: phon_ds
!arrays
 real(dp),intent(inout) :: qpt(3)
 real(dp),intent(out) :: displ_cart(2,3*phon_ds%natom,3*phon_ds%natom)
 real(dp),intent(out) :: pheigval(3*phon_ds%natom)
 real(dp),intent(out) :: pheigvec(2*3*phon_ds%natom*3*phon_ds%natom)
 real(dp),intent(out) :: phfrq(3*phon_ds%natom)

!Local variables-------------------------
!scalars
 integer :: symdynmat
 real(dp) :: qphnrm
!arrays
 real(dp) :: d2cart(2,3,phon_ds%mpert,3,phon_ds%mpert)

! *************************************************************************

!symdynmat = 0
 symdynmat = phon_ds%symdynmat
 qphnrm = one

!Generates the dynamical matrix from interatomic force constants and long-range electrostatic interactions.
 call gtdyn9(phon_ds%acell,phon_ds%atmfrc,phon_ds%dielt,phon_ds%dipdip,&
& phon_ds%dyewq0,d2cart,phon_ds%gmet,phon_ds%gprim,phon_ds%mpert,phon_ds%natom,&
& phon_ds%nrpt,qphnrm,qpt,phon_ds%rmet,phon_ds%rprim,phon_ds%rpt,&
& phon_ds%trans,phon_ds%ucvol,phon_ds%wghatm,phon_ds%xred,phon_ds%zeff)

!Diagonalize dynamical matrices getting the phonon frequencies and eigenvectors (as well as the corresponding displacements).
 call phfrq3(phon_ds%amu,displ_cart,d2cart,pheigval,pheigvec,phon_ds%indsym,&
& phon_ds%mpert,phon_ds%nsym,phon_ds%natom,phon_ds%nsym,phon_ds%ntypat,phfrq,&
& qphnrm,qpt,phon_ds%rprimd,symdynmat,phon_ds%symrel,phon_ds%typat,phon_ds%ucvol)


end subroutine inpphon
!!***
