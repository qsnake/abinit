!{\src2tex{textfont=tt}}
!!****f* ABINIT/scphon_dynmat_to_freq2
!! NAME
!! scphon_dynmat_to_freq2
!!
!! FUNCTION
!!  From the dynamical matrices, calculate the corresponding frequencies (squared)
!!  and see if reference phonon eigenvectors are still eigenvectors.
!!  FIXME: Kind of duplicates phfrq3 - should be merged or eliminated
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dynmat=dynamical matrix recalculated from phonon_eigval and phonon_eigvec_ref
!! natom_primitive_cell=number of atoms in primitive cell (not supercell used
!!   for SC phonon calculation)
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!! pcell = primitive cell information
!! phonon_eigvec_ref=reference phonon eigenvectors, from the anaddb equil run
!!
!! OUTPUT
!! phonon_eigval=phonon eigenfrequencies, updated inside SC phonon run
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon_new_frequencies
!!
!! CHILDREN
!!      wrtout,zgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_dynmat_to_freq2(dynmat,natom_primitive_cell,&
&     nphononq,symfreq2,pcell,phonon_eigvec_ref)

 use m_profiling

 use defs_basis
 use m_primcell_ddb_info

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_dynmat_to_freq2'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------

 !scalars
 integer, intent(in) :: natom_primitive_cell,nphononq
 type(primcell_ddb_info),intent(in) :: pcell
 !arrays
 real(dp),intent(in) :: phonon_eigvec_ref(2,3*natom_primitive_cell,&
&                                           3*natom_primitive_cell,nphononq)
 real(dp),intent(in) :: dynmat(2,3,natom_primitive_cell,3,natom_primitive_cell,nphononq)
 real(dp),intent(out) :: symfreq2(3*natom_primitive_cell,nphononq)

!Local variables -------------------------

 integer :: iq,imode,nmode,imode2
 integer :: iatom, iatom2
 real(dp) :: freqsign
 real(dp) :: z_one(2), z_zero(2)
 real(dp) :: tmpmat(2,3*natom_primitive_cell,3*natom_primitive_cell)
 real(dp) :: tmpmat2(2,3*natom_primitive_cell,3*natom_primitive_cell)
 real(dp) :: tmpeigvec(2,3*natom_primitive_cell,3*natom_primitive_cell)
 character(len=500) :: message

! *********************************************************************

 nmode=3*natom_primitive_cell
 z_one=(/one,zero/)
 z_zero=zero

 write (message,'(2a)') 'Follow eventual off diagonal elements of evec*',&
& ' dynmat evec, if any are > 1.e-10'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!loop over qpoints
 do iq=1,nphononq

   tmpeigvec=reshape(phonon_eigvec_ref(:,:,:,iq),&
&   (/2,3*natom_primitive_cell,3*natom_primitive_cell/))

   tmpmat = zero
   do iatom = 1, natom_primitive_cell
     do iatom2 = 1, natom_primitive_cell
       tmpmat(:,(iatom2-1)*3+1:(iatom2-1)*3+3,(iatom-1)*3+1:(iatom-1)*3+3) = dynmat(:,:,iatom2,:,iatom,iq)  &
&       / sqrt(pcell%amu(pcell%typat(iatom))*pcell%amu(pcell%typat(iatom2))) &
&       / amu_emass
     end do
   end do

!  ZGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!  multiply by eigvec left with trans conjugate
   tmpmat2=zero
   call zgemm ('C','N',nmode,nmode,nmode,z_one,&
&   tmpeigvec,nmode,&
&   tmpmat,nmode,&
&   z_zero,tmpmat2,nmode)

!  multiply by eigvec right, (no complex conjugation)
   tmpmat=zero
   call zgemm ('N','N',nmode,nmode,nmode,z_one,&
&   tmpmat2,nmode,&
&   tmpeigvec,nmode,&
&   z_zero,tmpmat,nmode)

   do imode=1,nmode
!    should we check imaginary part too?
     freqsign = one
     symfreq2(imode,iq) = tmpmat(1,imode,imode)
     do imode2=1,nmode
       if (imode2==imode) cycle
       if (abs(tmpmat(1,imode,imode2)) > tol10) then
         write (message,'(2I6,2x,2E12.3)') imode,imode2,tmpmat(:,imode,imode2)
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')
       end if
     end do
   end do
 end do ! iq

end subroutine scphon_dynmat_to_freq2
!!***
