!!****f* ABINIT/scphon_new_frequencies
!! NAME
!! scphon_new_frequencies
!!
!! FUNCTION
!! Calculate new frequencies from forces on supercell atoms, then symmetrize
!! them and add them to the averaged frequencies
!!
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! istep= number of the present iteration of the SC phonon calculations, for
!!   printing purposes
!! natom_primitive_cell=number of atoms in primitive cell (not supercell used
!!   for SC phonon calculation)
!! normal_mode_displacements= calculated displacements of canonical coordinates
!!   (normal modes) of phonons at desired temperature.
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!! nsym_primitive_cell= number of symmetries in the primitive unit cell
!! pcell=container type with ancillary variables and dimensions from anaddb run
!! phonon_eigvec_ref=reference phonon eigenvectors, from the anaddb equil run
!! phononq= phonon q vectors used in anaddb run (reduced coordinates)
!! qsym_map= map of qpoints onto one another, by sym ops:
!!   q_{qmap(iq,isym)} = S_{isym} q_{iq}
!!
!! OUTPUT
!! phonon_eigval=phonon eigenfrequencies, updated inside SC phonon run
!!
!! SIDE EFFECTS
!! phonon_eigval2_averaged= phonon frequencies squared, averaged over all
!! iterations to date
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!      mkrdim,print_phonfreq,scphon_dynmat_to_freq2,scphon_freq_to_dynmat
!!      symdyma
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_new_frequencies(forces_on_atoms_ft,istep,natom_primitive_cell,&
&    normal_mode_displacements,nphononq,nsym_primitive_cell,pcell,phonon_eigvec_ref,&
&    phonon_eigval2_averaged,phonon_eigval,phononq,qsym_map)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_primcell_ddb_info

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_new_frequencies'
 use interfaces_42_geometry
 use interfaces_45_geomoptim
 use interfaces_72_response
 use interfaces_95_drive, except_this_one => scphon_new_frequencies
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(primcell_ddb_info),intent(inout) :: pcell
 integer,intent(in) :: istep,natom_primitive_cell,nphononq,nsym_primitive_cell
!arrays
 integer,intent(in) :: qsym_map(nphononq,nsym_primitive_cell,2)
 real(dp),intent(in) :: forces_on_atoms_ft(2,3*natom_primitive_cell,nphononq)
 real(dp),intent(in) :: normal_mode_displacements(3*natom_primitive_cell,nphononq)
 real(dp),intent(in) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)
 real(dp),intent(in) :: phononq(3,nphononq)
 real(dp),intent(inout) :: phonon_eigval2_averaged(3*natom_primitive_cell,nphononq)
 real(dp),intent(out) :: phonon_eigval(3*natom_primitive_cell,nphononq)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir,imode_primitive_cell,iq,iq_image,isym
 integer :: itimrev,msize,multiplicity_q
 real(dp) :: timrev_sign
!arrays
 integer :: symqpoint_flag(nphononq)

 real(dp) :: scalprod_eigvec_force(2)
 real(dp) :: symfreq(3*natom_primitive_cell)
 real(dp) :: symfreq2(3*natom_primitive_cell, nphononq)
 real(dp) :: freq2(3*natom_primitive_cell, nphononq)

 real(dp) :: rprimd(3,3)

 real(dp), allocatable :: dynmat(:,:,:,:,:,:)

!************************************************************************

 call mkrdim(pcell%acell,pcell%rprim,rprimd)

!calculate SQUARE of new frequencies from forces
!NOTE: only real part of scalar product is kept
 freq2 = zero
 do iq=1,nphononq
   do iatom=1,natom_primitive_cell
     do idir=1,3
       imode_primitive_cell=idir+3*(iatom-1)
       if (abs(normal_mode_displacements(imode_primitive_cell,iq)) < tol12) cycle

!      here sum is over all atoms and directions
       scalprod_eigvec_force(1)=sum(phonon_eigvec_ref(1,:,imode_primitive_cell,iq)*forces_on_atoms_ft(1,:,iq)) &
&       -sum(phonon_eigvec_ref(2,:,imode_primitive_cell,iq)*forces_on_atoms_ft(2,:,iq))
       scalprod_eigvec_force(2)=sum(phonon_eigvec_ref(1,:,imode_primitive_cell,iq)*forces_on_atoms_ft(2,:,iq)) &
&       +sum(phonon_eigvec_ref(2,:,imode_primitive_cell,iq)*forces_on_atoms_ft(1,:,iq))
       freq2(imode_primitive_cell,iq)= -scalprod_eigvec_force(1) &
&       / normal_mode_displacements(imode_primitive_cell,iq)
     end do
   end do
 end do

 write(std_out,*) '  new frequencies '
 write(std_out,'(6E20.10)')  sign(one,freq2)*sqrt(abs(freq2))
 write(std_out,*)

!Symmetrize new frequencies (still in recip space)
 ABI_ALLOCATE(dynmat,(2,3,natom_primitive_cell,3,natom_primitive_cell,1))
 symqpoint_flag = 0
 do iq=1,nphononq
!  if qpoint already symmetrized skip it
   if (symqpoint_flag(iq) == 1) cycle

   symfreq(:)=zero
   multiplicity_q=0
!  sum over images of qpoint (I think mode number is conserved btw equivalent
!  qpoints - the only alternative is for degenerate modes, so we do not care)
   do itimrev=1,2
     do isym=1,nsym_primitive_cell

       iq_image = qsym_map(iq,isym,itimrev)
       if (iq_image==0) cycle

!      add this contribution for squared frequencies
       symfreq(:) = symfreq(:) + freq2(:,iq_image)
       multiplicity_q=multiplicity_q+1
     end do
     timrev_sign=-one
   end do

!  average and take square root
   symfreq(:) = symfreq/dble(multiplicity_q)
   symfreq(:) = sign(one,symfreq)*sqrt(abs(symfreq))

   msize=3*pcell%mpert*3*pcell%mpert
!  cast the freq into dynamical matrix
!  only called for 1 qpoint, iq
   call scphon_freq_to_dynmat(dynmat,natom_primitive_cell,&
&   1,pcell,symfreq,phonon_eigvec_ref(:,:,:,iq))

!  symmetrize the dynamical matrix
   call symdyma(dynmat,pcell%indsym,natom_primitive_cell,pcell%nsym,&
&   phononq(:,iq),rprimd,pcell%symrel)

!  re-extract the frequencies from the dynamical matrix, diagonalizing
   call scphon_dynmat_to_freq2(dynmat,natom_primitive_cell,&
&   1,symfreq2(:,iq),pcell,phonon_eigvec_ref(:,:,:,iq))

!  copy symmetrized value to images and flag them as done
   do itimrev=1,2
     do isym=1,nsym_primitive_cell
       iq_image=qsym_map(iq,isym,itimrev)
       if (iq_image==0) cycle

       symfreq2(:,iq_image) = symfreq2(:,iq)
       symqpoint_flag(iq_image) = 1
     end do
   end do

 end do ! iq

 ABI_DEALLOCATE(dynmat)

!Update average frequencies (averaged over all iterations up to now, counting
!the first iteration with the equilibrium phonon freq in step 0)

 phonon_eigval2_averaged=(dble(istep)*phonon_eigval2_averaged + symfreq2)/dble(istep+1)

!return to ACTUAL frequencies
 phonon_eigval(:,:) = sign(one,phonon_eigval2_averaged(:,:))*sqrt( abs(phonon_eigval2_averaged(:,:)) )

 call print_phonfreq(istep,natom_primitive_cell,nphononq,phonon_eigval)

 write(std_out,*) '  averaged frequencies ', istep
 write(std_out,'(6E20.10)')  phonon_eigval
 write(std_out,*)

end subroutine scphon_new_frequencies
!!***




