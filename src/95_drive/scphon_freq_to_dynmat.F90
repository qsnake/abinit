!{\src2tex{textfont=tt}}
!!****f* ABINIT/scphon_freq_to_dynmat
!! NAME
!! scphon_freq_to_dynmat
!!
!! FUNCTION
!!  From the updated frequencies and constant, reference, phonon eigenvectors
!!  this routine recalculates the dynamical matrices, which will be used
!!  for phonon interpolation.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! natom_primitive_cell=number of atoms in primitive cell (not supercell used
!!   for SC phonon calculation)
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!! pcell = primitive cell information
!! phonon_eigval=phonon eigenfrequencies, updated inside SC phonon run
!! phonon_eigvec_ref=reference phonon eigenvectors, from the anaddb equil run
!!
!! OUTPUT
!! dynmat=dynamical matrix recalculated from phonon_eigval and phonon_eigvec_ref
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon_interpolate_phonon_and_dos,scphon_new_frequencies
!!
!! CHILDREN
!!      zgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_freq_to_dynmat(dynmat,natom_primitive_cell,&
&     nphononq,pcell,phonon_eigval,phonon_eigvec_ref)

 use m_profiling

 use defs_basis
 use m_primcell_ddb_info

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_freq_to_dynmat'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom_primitive_cell,nphononq
 type(primcell_ddb_info),intent(in) :: pcell
!arrays
 real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)
 real(dp),intent(in) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)
 real(dp),intent(out) :: dynmat(2,3,natom_primitive_cell,3,natom_primitive_cell,nphononq)

!Local variables-------------------------------
!scalars
 integer :: imode,iq,nmode
 integer :: iatom,iatom2
!arrays
 real(dp) :: tmpeigvec(2,3*natom_primitive_cell,3*natom_primitive_cell)
 real(dp) :: tmpmat(2,3*natom_primitive_cell,3*natom_primitive_cell)
 real(dp) :: tmpmat2(2,3*natom_primitive_cell,3*natom_primitive_cell),z_one(2)
 real(dp) :: z_zero(2)

! *************************************************************************

 nmode=3*natom_primitive_cell
 z_one=(/one,zero/)
 z_zero=zero

!loop over qpoints
 do iq=1,nphononq

   tmpeigvec=reshape(phonon_eigvec_ref(:,:,:,iq),&
&   (/2,3*natom_primitive_cell,3*natom_primitive_cell/))

!  square eigenvalues
   tmpmat=zero
   do imode=1,3*natom_primitive_cell
     tmpmat(1,imode,imode)=sign(one,phonon_eigval(imode,iq))*phonon_eigval(imode,iq)**2
   end do

!  ZGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!  multiply by eigvec left (no complex conjugation)
   tmpmat2=zero
   call zgemm ('N','N',nmode,nmode,nmode,z_one,&
&   tmpeigvec,nmode,&
&   tmpmat,nmode,&
&   z_zero,tmpmat2,nmode)

!  multiply by eigvec right, with trans conjugate
   tmpmat=zero
   call zgemm ('N','C',nmode,nmode,nmode,z_one,&
&   tmpmat2,nmode,&
&   tmpeigvec,nmode,&
&   z_zero,tmpmat,nmode)

   dynmat(:,:,:,:,:,iq)=reshape(tmpmat,(/2,3,natom_primitive_cell,3,natom_primitive_cell/))

!  de-compensate for sqrt of masses
   do iatom = 1, natom_primitive_cell
     do iatom2 = 1, natom_primitive_cell
       dynmat(:,:,iatom2,:,iatom,iq) = dynmat(:,:,iatom2,:,iatom,iq)  &
&       * sqrt(pcell%amu(pcell%typat(iatom))*pcell%amu(pcell%typat(iatom2))) &
&       * amu_emass
     end do
   end do

 end do ! iq


end subroutine scphon_freq_to_dynmat
!!***
