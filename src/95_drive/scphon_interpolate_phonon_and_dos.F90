!!****f* ABINIT/scphon_interpolate_phonon_and_dos
!! NAME
!! scphon_interpolate_phonon_and_dos
!!
!! FUNCTION
!! Interpolate the phonon Density of States, from frequencies updated inside SC
!! phonon run
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
!! pcell=container type with ancillary variables and dimensions from anaddb run
!! phonon_eigval=phonon eigenfrequencies, updated inside SC phonon run
!! phonon_eigvec_ref=reference phonon eigenvectors, from the anaddb equil run
!! phononq= phonon q vectors used in anaddb run (reduced coordinates)
!! supercell_multiplicity=number of times the primitive unit cell is repeated
!!   along each axis, in the supercell
!!
!! OUTPUT
!! t_phonon_dos=type containing the DOS, partial DOS, and dimensions
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!      destroy_ddb_blk,ftifc_q2r,mkphbs,mkphdos,nullify_ddb_blk,print_phondos
!!      scphon_freq_to_dynmat
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_interpolate_phonon_and_dos (natom_primitive_cell,&
&    nphononq,pcell,t_phonon_dos,phonon_eigval,phonon_eigvec_ref,&
&    phononq,supercell_multiplicity)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_phdos
 use m_primcell_ddb_info
 use m_ddb_blk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_interpolate_phonon_and_dos'
 use interfaces_77_ddb
 use interfaces_95_drive, except_this_one => scphon_interpolate_phonon_and_dos
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom_primitive_cell,nphononq
 type(phonon_dos_type),intent(inout) :: t_phonon_dos
 type(primcell_ddb_info),intent(inout) :: pcell
!arrays
 integer,intent(in) :: supercell_multiplicity(3)
 real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)
 real(dp),intent(in) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)
 real(dp),intent(in) :: phononq(3,nphononq)

!Local variables-------------------------------
!scalars
 integer :: msize
 integer :: dummy_iodyn
 real(dp) :: tcpui,twalli
 type(anaddb_dataset_type) :: dummy_anaddb_dtset
 type(ddb_blk_type), pointer :: ddb_blk_dummy

!arrays
 integer :: ngqpt(3)
 real(dp),allocatable :: atmfrc(:,:,:,:,:,:),dynmat(:,:,:,:,:,:)

 real(dp), allocatable :: dummy_uinvers(:,:), dummy_vtinvers(:,:), dummy_singular(:), dummy_d2asr(:,:)

 character(len=fnlen) :: bsfilename
! *************************************************************************

!Source code

 write(std_out,*) ' in scphon_interpolate_phonon_and_dos'

!test if brav was used in initial anaddb run
 if(pcell%brav /= 1) stop ' scphon_interpolate_phonon_and_dos : Error - brav/=1 not coded '

 dummy_anaddb_dtset%dipdip = pcell%dipdip    ! include dipole dipole interaction term?
!
!The following initializations will be moved to a central routine later
!

!this is the fine grid on which the phonons are interpolated
 ngqpt(:) = 40

!fill defaults for anaddb_dtset and other inputs to mkphdos
 dummy_anaddb_dtset%prtdos = 1 ! use gaussian integ by default
 dummy_anaddb_dtset%dosdeltae = t_phonon_dos%omega_step
 dummy_anaddb_dtset%dossmear = t_phonon_dos%dossmear
 dummy_anaddb_dtset%brav = 1   ! do not use brav symmetrization for the moment
 dummy_anaddb_dtset%ngqpt = 0  ! for brav==1 the other 6 values should not be used
 dummy_anaddb_dtset%ngqpt(1:3) = supercell_multiplicity
 dummy_anaddb_dtset%nqshft = 1 ! force 1 shift to 000
 dummy_anaddb_dtset%q1shft = zero
 dummy_anaddb_dtset%ng2qpt = ngqpt
 dummy_anaddb_dtset%q2shft = zero
 dummy_anaddb_dtset%symdynmat = 1 ! force symmetrization of dynamical matrix
!END initializations indep of presently calculated frequencies

!reconstitute dynamical matrices on qpoints we know
 msize=3*pcell%mpert*3*pcell%mpert
 ABI_ALLOCATE(dynmat,(2,3,natom_primitive_cell,3,natom_primitive_cell,nphononq))
 call scphon_freq_to_dynmat(dynmat,natom_primitive_cell,&
& nphononq,pcell,phonon_eigval,phonon_eigvec_ref)

!calculate atomic force constants
 ABI_ALLOCATE(atmfrc,(2,3,natom_primitive_cell,3,natom_primitive_cell,pcell%nrpt))
 call ftifc_q2r(atmfrc,dynmat,pcell%gprim,pcell%natom,nphononq,&
& pcell%nrpt,pcell%rpt,phononq)
 ABI_DEALLOCATE(dynmat)

!interpolate the DOS (should also work for tetrahedron method)
 tcpui=zero
 twalli=zero

 call mkphdos(t_phonon_dos,dummy_anaddb_dtset%prtdos,dummy_anaddb_dtset%dosdeltae,dummy_anaddb_dtset%dossmear,&
& dummy_anaddb_dtset%dipdip,dummy_anaddb_dtset%symdynmat,&
& pcell%acell,pcell%amu,dummy_anaddb_dtset,atmfrc,pcell%dielt,pcell%dyewq0,&
& pcell%gmet,pcell%gprim,pcell%indsym,&
& pcell%mpert,pcell%msym,pcell%natom,pcell%nrpt,pcell%nsym,&
& t_phonon_dos%ntypat,pcell%rmet,pcell%rprim,&
& pcell%rpt,pcell%symrec,pcell%symrel,&
& pcell%trans,pcell%typat,pcell%ucvol,pcell%wghatm,pcell%xred,pcell%zeff)

 call print_phondos(t_phonon_dos,"PHDOS")

!make band structure as well:
 dummy_anaddb_dtset%nph1l = 0
 ABI_ALLOCATE(dummy_anaddb_dtset%qph1l,(3,1))
 
 dummy_anaddb_dtset%nqpath=5
 ABI_ALLOCATE(dummy_anaddb_dtset%qpath,(3,dummy_anaddb_dtset%nqpath))
 dummy_anaddb_dtset%qpath(:,1) = (/zero, zero, zero/)
 dummy_anaddb_dtset%qpath(:,2) = (/half, zero, zero/)
 dummy_anaddb_dtset%qpath(:,3) = (/half, half, zero/)
 dummy_anaddb_dtset%qpath(:,4) = (/zero, zero, zero/)
 dummy_anaddb_dtset%qpath(:,5) = (/half, half, half/)
 dummy_anaddb_dtset%eivec=3
 dummy_anaddb_dtset%outscphon = 0
 ABI_ALLOCATE(dummy_anaddb_dtset%qnrml1,(1))
 dummy_anaddb_dtset%qnrml1(1) = one
 dummy_anaddb_dtset%ifcflag = 1
 dummy_anaddb_dtset%rfmeth = 0
 dummy_anaddb_dtset%asr = 2
 dummy_anaddb_dtset%thmflag = 0
 dummy_anaddb_dtset%freeze_displ = zero
 dummy_anaddb_dtset%eivec = 1
 dummy_anaddb_dtset%enunit = 0


 ABI_ALLOCATE(dummy_uinvers,(1,1))
 ABI_ALLOCATE(dummy_vtinvers,(1,1))
 ABI_ALLOCATE(dummy_singular,(1))
 ABI_ALLOCATE(dummy_d2asr,(1,1))
 dummy_iodyn = 0

 ABI_ALLOCATE(ddb_blk_dummy,)
 call nullify_ddb_blk(ddb_blk_dummy)

 bsfilename = "abiscphon_interpolated_BS"
 call mkphbs(pcell%acell,pcell%amu,dummy_anaddb_dtset,atmfrc,ddb_blk_dummy,&
& dummy_d2asr,pcell%dielt,pcell%dyewq0,&
& bsfilename,pcell%gmet,pcell%gprim,pcell%indsym,&
& dummy_iodyn,pcell%mpert,msize,pcell%msym,pcell%natom,pcell%nrpt,pcell%nsym,&
& t_phonon_dos%ntypat,tol10,pcell%rmet,pcell%rprim,&
& pcell%rpt,dummy_singular,pcell%symrel,tcpui,pcell%trans,twalli,pcell%typat,pcell%ucvol,&
& dummy_uinvers,dummy_vtinvers,pcell%wghatm,pcell%xred,pcell%zeff)

 call destroy_ddb_blk(ddb_blk_dummy)
 ABI_DEALLOCATE(ddb_blk_dummy)

 ABI_DEALLOCATE(dummy_anaddb_dtset%qph1l)
 ABI_DEALLOCATE(dummy_anaddb_dtset%qpath)
 ABI_DEALLOCATE(dummy_anaddb_dtset%qnrml1)
 ABI_DEALLOCATE(dummy_uinvers)
 ABI_DEALLOCATE(dummy_vtinvers)
 ABI_DEALLOCATE(dummy_singular)
 ABI_DEALLOCATE(dummy_d2asr)

 ABI_DEALLOCATE(atmfrc)

end subroutine scphon_interpolate_phonon_and_dos
!!***
