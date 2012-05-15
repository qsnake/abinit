!{\src2tex{textfont=tt}}
!!****f* ABINIT/scphon
!! NAME
!! scphon
!!
!! FUNCTION
!!  Using a supercell, calculate a self consistent phonon structure
!!  as in PRL 100 095901 (2008). The initial phonon eigenvectors and
!!  eigenvalues are read in, and then atoms are displaced according
!!  to the normal modes populated at a given temperature until
!!  convergence of the vibrational free energy (or so I hope)
!!
!! Other references:
!!   Computational Materials Science 44 (2009) 888-894
!!   Phys Rev B 78, 184304 (2008)
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
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cpus= cpu time limit in seconds
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs for the "coarse" grid (see NOTES below)
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in cell.
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!
!! SIDE EFFECTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=updated wavefunctions;  if mkmem>=nkpt,
!!         these are kept in a disk file.
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  occ(mband*nkpt*nsppol)=occupation number for each band (often 2) at each k point
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!   (should be made a pure output quantity)
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in el./bohr**3
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  wffnew,wffnow=struct info for wf disk files.
!!  wvl <type(wvl_data)>=all wavelets data.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)= at input, previous reduced dimensionless atomic coordinates
!!                     at output, current xred is transferred to xred_old
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      destroy_primcell_ddb_info,init_phondos,print_phonfreq,prtxvf
!!      read_primcell_ddb_info,scfcv_new,scphon_build_qsym_map
!!      scphon_check_fcart,scphon_free_energy,scphon_freq_to_normmode
!!      scphon_ft_fcart,scphon_interpolate_phonon_and_dos
!!      scphon_make_phonon_dos,scphon_new_frequencies,scphon_phonon_init
!!      scphon_qpoint_init,scphon_supercell_vectors_init,scphon_update_xcart
!!      wrtout,xcomm_world,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon(amass,ab_scfcv_in,ab_scfcv_inout,&
& dtset, electronpositron,&
& paw_dmft,&
& rhog, rhor,&
& rprimd,&
& wffnew, wffnow,&
& xred, xred_old)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_wffile
 use defs_abitypes
 use defs_scftypes
 use m_electronpositron, only : electronpositron_type
 use m_primcell_ddb_info
 use m_phdos
 use m_paw_dmft, only: paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_45_geomoptim
 use interfaces_51_manage_mpi
 use interfaces_67_common
 use interfaces_95_drive, except_this_one => scphon
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ab_scfcv_args_in),intent(in) :: ab_scfcv_in
 type(ab_scfcv_args_inout),intent(inout) :: ab_scfcv_inout
 type(dataset_type),intent(inout) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(paw_dmft_type) :: paw_dmft
 type(wffile_type),intent(inout) :: wffnew,wffnow
 real(dp), intent(in) :: amass(dtset%natom)
 real(dp), intent(inout) :: rprimd(3,3)
 real(dp), pointer :: rhog(:,:),rhor(:,:)
 real(dp), intent(inout) :: xred(3,dtset%natom),xred_old(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: comm,iatom,iStep,iexit,me
 integer :: natom_primitive_cell,nphononq,nsym_primitive_cell
 real(dp) :: free_energy,old_free_energy,scphon_temp
 character(len=32) :: statusOut
 character(len=500) :: message
 character(len=fnlen) :: ddb_info_filename
 type(phonon_dos_type) :: t_phonon_dos
 type(primcell_ddb_info) :: pcell
!arrays
 integer :: supercell_multiplicity(3)
 integer,allocatable :: pcell_atom_in_supercell(:),qsym_map(:,:,:)
 integer,allocatable :: symrec_primitive_cell(:,:,:)
 real(dp),allocatable :: amass_pcell(:)
 real(dp),allocatable :: cartesian_displacements(:,:),forces_on_atoms_ft(:,:,:)
 real(dp),allocatable :: normal_mode_displacements(:,:)
 real(dp),allocatable :: normal_mode_displacements_old(:,:),phonon_eigval(:,:)
 real(dp),allocatable :: phonon_eigval2_averaged(:,:),phonon_eigval_ref(:,:)
 real(dp),allocatable :: phonon_eigvec_ref(:,:,:,:),phononq(:,:)
 real(dp),allocatable :: sqrt_amass_pcell(:),supercell_vectors(:,:),vel(:,:)
 real(dp),allocatable :: xcart(:,:),xcart0(:,:)
 

! ************************************************************************

!write(std_out,*) 'scphon 01'
!##########################################################
!### 01. Initial settings, only for first iteration

 write (message,'(a)') ' scphon: enter '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!the number of times the primitive unit cell has been replicated in each
!dimension
!FIXME: this could be inferred from the supercell size, avoiding possible mismatches
!of scphon_supercell with the real one...
 supercell_multiplicity = dtset%scphon_supercell
 write (message,'(a,3I6)') ' SC phonons: Found supercell multiplicity of', &
& supercell_multiplicity
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!this is the number of atoms for the primitive unit cell, in which the phonons
!were calculated
 natom_primitive_cell=nint(dble(dtset%natom)&
& /dble(supercell_multiplicity(1))&
& /dble(supercell_multiplicity(2))&
& /dble(supercell_multiplicity(3)))
 write (message,'(a,I6,a)') 'Deduced that there are ', natom_primitive_cell, &
& ' atoms in the primitive unit cell'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 ABI_ALLOCATE(amass_pcell,(natom_primitive_cell))
 ABI_ALLOCATE(sqrt_amass_pcell,(natom_primitive_cell))
 amass_pcell(:) = amass(1:natom_primitive_cell)
 sqrt_amass_pcell(:) = sqrt(amass_pcell(:))
 do iatom=1,natom_primitive_cell
   write (message,'(a,I6,a,E20.10)') ' mass of atom ', iatom, ' = ', amass_pcell(iatom)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end do

!the temperature we want the phonons to be at
 scphon_temp = dtset%scphon_temp

!set energies for DOS output and integration
!t_phonon_dos%dossmear  =0.00001   ! about 2. cm-1
!t_phonon_dos%omega_step=0.0000047 ! about 1. cm-1
!t_phonon_dos%nomega=1000
 call init_phondos(t_phonon_dos,ab_scfcv_in%psps%ntypat,dtset%natom,1,1000,1,1,&
& smallest_real,greatest_real,0.0000047_dp,0.00001_dp)
 write(std_out,*) 't_phonon_dos%nomega t_phonon_dos%dossmear, t_phonon_dos%omega_step ',&
& 't_phonon_dos%omega_min, t_phonon_dos%omega_max'
 write(std_out,*) t_phonon_dos%nomega,  t_phonon_dos%dossmear, t_phonon_dos%omega_step, &
& t_phonon_dos%omega_min, t_phonon_dos%omega_max


!END input variables

!write(std_out,*) 'scphon 02'
!##########################################################
!### 02. Several allocations also for first time

!
!Allocate working arrays
!
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 ABI_ALLOCATE(xcart0,(3, dtset%natom))
 ABI_ALLOCATE(vel,(3, dtset%natom))
!Transform xred to cartesian coordinates.
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)
!initial positions
 xcart0=xcart

!number of qpoints for phonon integration should be = natom
!could be reduced by symmetries for reciprocal space,
!but then you need to reconstruct the real space atomic displacements below
!This is not implemented yet.
 nphononq=dtset%natom/natom_primitive_cell

!for each atom, the vector to the unit cell it is in
 ABI_ALLOCATE(supercell_vectors,(3,dtset%natom))
 ABI_ALLOCATE(pcell_atom_in_supercell,(dtset%natom))

!phonon q vectors we will consider. Should correspond to the superlattice
!of atoms we are using, eg 4x4x4 supercell
 ABI_ALLOCATE(phononq,(3,nphononq))

!displacements of atoms from equilibrium positions
!written U_R in the PRL
!allocate(atom_displacements(3,dtset%natom))


!the following depend on natom_primitive_cell

!eigvec of reference phonon system
 ABI_ALLOCATE(phonon_eigvec_ref,(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq))

!eigval of reference phonon system
 ABI_ALLOCATE(phonon_eigval_ref,(3*natom_primitive_cell,nphononq))

!eigval of present phonon system, extracted from instantaneous forces
 ABI_ALLOCATE(phonon_eigval,(3*natom_primitive_cell,nphononq))

!eigval of present phonon system: averaged value of phonon_eigval
 ABI_ALLOCATE(phonon_eigval2_averaged,(3*natom_primitive_cell,nphononq))

!classical displacements along each of the 3*natom normal modes
!written A_ks in the reference PRL, but without the 1/sqrt(mass) factor
 ABI_ALLOCATE(normal_mode_displacements,(3*natom_primitive_cell,nphononq))
 ABI_ALLOCATE(normal_mode_displacements_old,(3*natom_primitive_cell,nphononq))

!the absolute displacement of the atoms, in cartesian coordinates, since
!the beginning of time
 ABI_ALLOCATE(cartesian_displacements,(3,dtset%natom))

!forces on each atom fcart are fourier transformed into the following array
!(in cartesian coordinates)
 ABI_ALLOCATE(forces_on_atoms_ft,(2,3*natom_primitive_cell,nphononq))


!write(std_out,*) 'scphon 03'
!##########################################################
!### 03. Initializations

!get generic phonon info from anaddb run, for later interpolation
 write(std_out,*) ' entering read_primcell'
 ddb_info_filename=trim(ab_scfcv_inout%dtfil%filnam_ds(3))//'_PCINFO'
 ddb_info_filename=trim(ddb_info_filename)
 call read_primcell_ddb_info(ddb_info_filename,pcell)

!initialize the supercell vector grid (pointing to each unit cell in the supercell)
 call scphon_supercell_vectors_init(dtset%natom,natom_primitive_cell,&
& pcell,pcell_atom_in_supercell,supercell_multiplicity,&
& supercell_vectors,xred)

 write (message,'(a)') '  xred = '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 do iatom=1,dtset%natom
   write (message,'(3E20.10)')  xred(:,iatom)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end do

!initialize the phonon q grid based on the supercell size
 call scphon_qpoint_init (nphononq,phononq,supercell_multiplicity)

!initialize the reference phonon eigenvectors and eigenvalues from file
 call scphon_phonon_init (ab_scfcv_inout%dtfil%fnameabi_phfrq,&
& ab_scfcv_inout%dtfil%fnameabi_phvec,natom_primitive_cell,&
& nphononq,phonon_eigvec_ref,phonon_eigval_ref)

!set the initial phonon frequency average
 phonon_eigval2_averaged(:,:)=sign(one,phonon_eigval_ref(:,:))*phonon_eigval_ref(:,:)**2

!get the symmetry operations for the primitive unit cells
 nsym_primitive_cell=pcell%nsym
 ABI_ALLOCATE(symrec_primitive_cell,(3,3,nsym_primitive_cell))
 symrec_primitive_cell(:,:,1:nsym_primitive_cell) = pcell%symrec

!mapping from 1 qpoint to another under symmetries
 ABI_ALLOCATE(qsym_map,(nphononq,nsym_primitive_cell,2))
!initialize the mapping from one qpoint to another under the symops
 call scphon_build_qsym_map(nphononq,nsym_primitive_cell,phononq,&
& qsym_map,symrec_primitive_cell)

!initialize the first normal mode displacements
 call scphon_freq_to_normmode (qsym_map(:,1,2),natom_primitive_cell,normal_mode_displacements,&
& nphononq,phonon_eigval_ref,scphon_temp)
 write (message,'(a)') ' Have initialized normal_mode_displacements'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!update xcart with difference between old and new normal mode displacements
!also initializes cartesian_displacements
 cartesian_displacements=zero
 call scphon_update_xcart (sqrt_amass_pcell,cartesian_displacements,dtset%natom,natom_primitive_cell,&
& normal_mode_displacements,&
& nphononq,pcell_atom_in_supercell,phonon_eigvec_ref,phononq,supercell_vectors,xcart,xcart0)

!in first iteration, apply full normal mode displacement
 normal_mode_displacements_old=normal_mode_displacements

 call print_phonfreq(-1,natom_primitive_cell,nphononq,phonon_eigval_ref)

!set old_free_energy to make sure at least 2 iterations are performed.
 old_free_energy=greatest_real

!write(std_out,*) 'scphon 04'
!##########################################################
!### 04. Begin scphon relaxation.

 do iStep = 1, dtset%ntime, 1
   write(message, '(a,a,i3,a)' ) ch10, ' SCPHON STEP NUMBER ', iStep,&
&   '  ---------------------------------------------------------'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')


!  transform to reduced coordinates for call to scfcv
   call xredxcart(dtset%natom, -1, rprimd, xcart, xred)
   write (message,'(a)') '  input xred = '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   do iatom=1,dtset%natom
     write (message,'(3E20.10)')  xred(:,iatom)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end do
   write (message,'(a)') ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

!  write(std_out,*) 'scphon 05'
!  ##########################################################
!  ### 05. SCFCV calculates electronic structure and forces

   call scfcv_new(ab_scfcv_in,ab_scfcv_inout,&
&   dtset,electronpositron,&
&   paw_dmft,&
&   rhog,rhor,&
&   rprimd,&
&   wffnew,wffnow,&
&   xred,xred_old)

!  write(std_out,*) 'scphon 06'
!  ##########################################################
!  ### 06. Output coordinates and forces

!  Output coordinates and forces (not velocities, prtvel = 0) and total energy
   call prtxvf(ab_scfcv_inout%results_gs%fcart,ab_scfcv_inout%results_gs%fred,dtset%iatfix, ab_out, dtset%natom, &
&   0, vel, xcart,xred)
   call prtxvf(ab_scfcv_inout%results_gs%fcart,ab_scfcv_inout%results_gs%fred,dtset%iatfix, 06 , dtset%natom, &
&   0, vel, xcart,xred)

!  write(std_out,*) 'scphon 07'
!  ##########################################################
!  ### 07. SCPHON Specialities

!  Check if fcart is in the opposite direction to cartesian_displacements for
!  each atom
   call scphon_check_fcart(cartesian_displacements,ab_scfcv_inout%results_gs%fcart,dtset%natom)

!  Fourier transform forces to recirocal space (inside BZ)
!  FIXME: could be optimized perhaps, make it an FFT, but the grid will never be
!  very big.
   call scphon_ft_fcart(sqrt_amass_pcell,ab_scfcv_inout%results_gs%fcart,dtset%natom,natom_primitive_cell,nphononq,phononq,&
&   forces_on_atoms_ft,pcell_atom_in_supercell,supercell_vectors)

!  Determine new frequencies according to forces
   call scphon_new_frequencies(forces_on_atoms_ft,istep,natom_primitive_cell,&
&   normal_mode_displacements,nphononq,nsym_primitive_cell,pcell,phonon_eigvec_ref,&
&   phonon_eigval2_averaged,phonon_eigval,phononq,qsym_map)

!  calculate the DOS
!  1. the following code is for the bare input phonon frequencies (no interpolation)
   t_phonon_dos%nomega=1000
   if (associated(t_phonon_dos%phdos))  then
     ABI_DEALLOCATE(t_phonon_dos%phdos)
   end if
   ABI_ALLOCATE(t_phonon_dos%phdos,(t_phonon_dos%nomega))
   call scphon_make_phonon_dos (t_phonon_dos%dossmear,natom_primitive_cell,&
&   t_phonon_dos%nomega,nphononq, t_phonon_dos%omega_max,t_phonon_dos%omega_min,&
&   t_phonon_dos%phdos,phonon_eigval)

!  2. interpolate phonon freq first, then calculate DOS
!  this will give nomega a new value!
   call scphon_interpolate_phonon_and_dos (natom_primitive_cell,&
&   nphononq,pcell,t_phonon_dos,phonon_eigval,phonon_eigvec_ref,&
&   phononq,supercell_multiplicity)

!  calculate vibrational free energy
   call scphon_free_energy(free_energy,istep,t_phonon_dos,scphon_temp)

!  Check whether Free Energy diff is below tolerance; if so, exit
!  from the istep loop
!  NOTE: tolmxf is usually in Ha/bohr, so we should probably have a conversion
!  factor like for the stress. F usually down to 1e-5 and energy down to 1e-10
!  or so...
   if (abs(free_energy-old_free_energy) < dtset%tolmxf) then
     write (message,'(a)') ' convergence in vibrational free energy has been reached: '
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write (message,'(E15.5,a,E15.5)') free_energy-old_free_energy, ' < ', dtset%tolmxf
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     statusOut = "OK"
     iexit = 1
     exit
   end if

   iexit = 0
   if (iStep == dtset%ntime) then
     iexit = 1
     statusOut = "Failed"
   end if

!  determine new classical normal mode displacements
   call scphon_freq_to_normmode (qsym_map(:,1,2),natom_primitive_cell,normal_mode_displacements,&
&   nphononq,phonon_eigval,scphon_temp)

!  update xcart with difference between old and new normal mode displacements
   call scphon_update_xcart (sqrt_amass_pcell,cartesian_displacements,&
&   dtset%natom,natom_primitive_cell,&
&   normal_mode_displacements-normal_mode_displacements_old,&
&   nphononq,pcell_atom_in_supercell,phonon_eigvec_ref,phononq,supercell_vectors,xcart,xcart0)

!  save present displacements for next step
   normal_mode_displacements_old=normal_mode_displacements

   old_free_energy=free_energy

 end do !                       end istep do

 call xredxcart(dtset%natom, -1, rprimd, xcart, xred)

!Free working arrays
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xcart0)
 ABI_DEALLOCATE(vel)
 ABI_DEALLOCATE(supercell_vectors)
 ABI_DEALLOCATE(phononq)
 ABI_DEALLOCATE(phonon_eigvec_ref)
 ABI_DEALLOCATE(phonon_eigval_ref)
 ABI_DEALLOCATE(phonon_eigval)
 ABI_DEALLOCATE(phonon_eigval2_averaged)
 ABI_DEALLOCATE(normal_mode_displacements)
 ABI_DEALLOCATE(normal_mode_displacements_old)
 ABI_DEALLOCATE(forces_on_atoms_ft)
 ABI_DEALLOCATE(qsym_map)

 call destroy_primcell_ddb_info(pcell)


!XML output of the status
 if (dtset%prtxml == 1) then
   call xcomm_world(ab_scfcv_inout%mpi_enreg,comm,myrank=me)
   if (me == 0 .and. dtset%prtxml == 1) then
     write(ab_xml_out, "(A)") '    <geometryMinimisation type="scphon">'
     write(ab_xml_out, "(A,A,A)") '      <status cvState="', trim(statusOut) , &
&     '" stop-criterion="tolmxf" />'
     write(ab_xml_out, "(A)") '    </geometryMinimisation>'
   end if
 end if

end subroutine scphon
!!***

