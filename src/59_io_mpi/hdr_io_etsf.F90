!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_io_etsf
!! NAME
!! hdr_io_etsf
!!
!! FUNCTION
!! This subroutine deals with the I/O of the hdr_type
!! structured variables (read/write/echo).
!! It handles variables according to the ETSF format, whenever
!! possible and uses new variables when not available in the ETSF
!! format.
!! According to the value of rdwr, it reads the header
!! of a file, writes it, or echo the value of the structured
!! variable to a file.
!! Note that, when reading, different records of hdr
!! are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated
!! correctly by a call to hdr_clean when hdr is not used anymore.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  rdwr= if 1, read the hdr structured variable from the header of the netCDF file,
!!        if 2, write the header to unformatted netCDF file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!        if 5, read the hdr without rewinding (unformatted), identical to 1 for netCDF
!!        if 6, read the hdr without rewinding (unformatted), identical to 2 for netCDF
!!  unitwff=the unit of the open NetCDF file.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!  hdr <type(hdr_type)>=the header structured variable
!!   if rdwr=1,5 : will be output
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!
!! NOTES
!!
!! PARENTS
!!      cut3d,inwffil,ioarr,kss2wfk,m_ebands,m_io_kss,m_io_screening
!!      m_sigma_results,m_wfs,outwf,pawmkaewf,wfk_read_ene
!!
!! CHILDREN
!!      etsf_io_basisdata_get,etsf_io_basisdata_put,etsf_io_dims_get
!!      etsf_io_electrons_get,etsf_io_electrons_put,etsf_io_geometry_get
!!      etsf_io_geometry_put,etsf_io_kpoints_get,etsf_io_kpoints_put
!!      etsf_io_low_read_dim,etsf_io_low_read_var,etsf_io_low_set_write_mode
!!      etsf_io_low_write_var,hdr_io_int,leave_new,strip,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine hdr_io_etsf(fform,hdr,rdwr,unitwff)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io_low_level
 use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_io_etsf'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_59_io_mpi, except_this_one => hdr_io_etsf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rdwr,unitwff
 integer,intent(inout) :: fform
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
#if defined HAVE_TRIO_ETSF_IO
 type(etsf_dims) :: dims
 type(etsf_kpoints) :: kpoints
 type(etsf_basisdata) :: basisdata
 type(etsf_geometry) :: geometry
 type(etsf_electrons) :: electrons
 type(etsf_io_low_error) :: error_data
 logical :: lstat
 character(len = 500)  :: message
 integer :: rhoijdim1, rhoijdim2, nresolution
 integer :: headform, iatom, itypat
 character(len=etsf_charlen), target :: basis_set
 real(dp), target :: ecut, fermie
 real(dp), target :: rprimd(3, 3)
!temp variables
! integer :: cplex, ilmn, irhoij, ispden, lmn2_size, nselect
! real(dp), allocatable :: rhoij(:,:,:)
#endif

! *************************************************************************

#if defined HAVE_TRIO_ETSF_IO
 write(message, '(A,I0)' ) ' hdr_io_etsf: accessing ABINIT specific data from unit ', unitwff
 call wrtout(std_out, message, 'COLL')

 if(rdwr==1 .or. rdwr==5)then
!  We switch off from define mode.
   call etsf_io_low_set_write_mode(unitwff, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

!  In case the file is just ETSF valid, we ignore the missing variables
!  and we use defualt values.
   hdr%codvsn   = "ETSF  "
   fform        = 1
   headform     = 57

!  First, we read the declaration of code, fform ...
!  We ignore errors, assuming that the file is at least ETSF valid.
   call etsf_io_low_read_var(unitwff, "codvsn", hdr%codvsn, 6, lstat, error_data = error_data)
   if (lstat) then ! We pad the returned string with " " instead of "\0"
     call strip(hdr%codvsn)
   end if

   call etsf_io_low_read_var(unitwff, "fform", fform, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_read_var(unitwff, "headform", hdr%headform, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   if (headform <= 42) then
     write(message,'(a,i0,3a)')&
&     '  headform is ',headform,', while it should be > 42.',ch10,&
&     '  Action : check the correctness of your file.'
     MSG_ERROR(message)
   end if

!  Then, read dimensions handled by ETSF
   call etsf_io_dims_get(unitwff, dims, lstat, error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

!  Copy dimensions to hdr structure
!  FIXME: don't handle k_dependent = 1
   hdr%bantot   = dims%max_number_of_states * dims%number_of_kpoints * dims%number_of_spins
   hdr%natom    = dims%number_of_atoms
   hdr%nkpt     = dims%number_of_kpoints
   hdr%nspden   = dims%number_of_components
   hdr%nspinor  = dims%number_of_spinor_components
   hdr%nsppol   = dims%number_of_spins
   hdr%nsym     = dims%number_of_symmetry_operations
   hdr%ntypat   = dims%number_of_atom_species
   hdr%ngfft(1) = dims%number_of_grid_points_vector1
   hdr%ngfft(2) = dims%number_of_grid_points_vector2
   hdr%ngfft(3) = dims%number_of_grid_points_vector3

!  We read other dimensions, not handled by ETSF format.
!  In case the file is just ETSF valid, we ignore the missing dimensions and we use default values.

   hdr%npsp    = hdr%ntypat
   call etsf_io_low_read_dim(unitwff, "npsp", hdr%npsp, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   rhoijdim1   = 1
   call etsf_io_low_read_dim(unitwff, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   rhoijdim2   = hdr%nspden
   call etsf_io_low_read_dim(unitwff, "rhoijdim2", rhoijdim2, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%usepaw  = 0
   call etsf_io_low_read_var(unitwff, "usepaw", hdr%usepaw, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%usewvl  = 0
   call etsf_io_low_read_var(unitwff, "usewvl", hdr%usewvl, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   nresolution=0
   if (hdr%usewvl == 1) then
!    This value must be 2...
     call etsf_io_low_read_dim(unitwff, "number_of_wavelet_resolutions",nresolution, lstat, error_data = error_data)
     ETSF_WARN(lstat,error_data)

!    We set the right ngfft, adding the padding space for wavelets.
     hdr%ngfft = hdr%ngfft + 31
   end if

!  Allocate all parts of hdr that need to be
   ABI_ALLOCATE(hdr%istwfk,(hdr%nkpt))
   ABI_ALLOCATE(hdr%lmn_size,(hdr%npsp))
   ABI_ALLOCATE(hdr%nband,(hdr%nkpt*hdr%nsppol))
   ABI_ALLOCATE(hdr%npwarr,(hdr%nkpt))
   ABI_ALLOCATE(hdr%pspcod,(hdr%npsp))
   ABI_ALLOCATE(hdr%pspdat,(hdr%npsp))
   ABI_ALLOCATE(hdr%pspso,(hdr%npsp))
   ABI_ALLOCATE(hdr%pspxc,(hdr%npsp))
   ABI_ALLOCATE(hdr%so_psp,(hdr%npsp))
   ABI_ALLOCATE(hdr%symafm,(hdr%nsym))
   ABI_ALLOCATE(hdr%symrel,(3,3,hdr%nsym))
   ABI_ALLOCATE(hdr%typat,(hdr%natom))
   ABI_ALLOCATE(hdr%kptns,(3,hdr%nkpt))
   ABI_ALLOCATE(hdr%occ,(hdr%bantot))
   ABI_ALLOCATE(hdr%tnons,(3,hdr%nsym))
   ABI_ALLOCATE(hdr%wtk,(hdr%nkpt))
   ABI_ALLOCATE(hdr%xred,(3,hdr%natom))
   ABI_ALLOCATE(hdr%znuclpsp,(hdr%npsp))
   ABI_ALLOCATE(hdr%znucltypat,(hdr%ntypat))
   ABI_ALLOCATE(hdr%zionpsp,(hdr%npsp))
   ABI_ALLOCATE(hdr%title,(hdr%npsp))
   if(hdr%usepaw==1)  then
     ABI_ALLOCATE(hdr%pawrhoij,(hdr%natom))
   end if

!  We get then all variables included in ETSF
   if (hdr%usewvl == 0) then
     basisdata%kinetic_energy_cutoff => ecut
     basisdata%number_of_coefficients => hdr%npwarr
     call etsf_io_basisdata_get(unitwff, basisdata, lstat, error_data)
     ETSF_CHECK_ERROR(lstat,error_data)
   else
     call etsf_io_low_read_var(unitwff, "number_of_wavelets", hdr%nwvlarr, lstat, error_data = error_data)
     ETSF_CHECK_ERROR(lstat,error_data)
   end if

   electrons%fermi_energy => fermie
   electrons%number_of_states%data1D => hdr%nband
   electrons%occupations%data1D => hdr%occ

   call etsf_io_electrons_get(unitwff, electrons, lstat, error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   geometry%primitive_vectors => rprimd
   geometry%reduced_symmetry_matrices => hdr%symrel
   geometry%atom_species => hdr%typat
   geometry%reduced_symmetry_translations => hdr%tnons
   geometry%reduced_atom_positions => hdr%xred
   geometry%atomic_numbers => hdr%znucltypat

   if (hdr%npsp == hdr%ntypat) then
     geometry%valence_charges => hdr%zionpsp
   end if

   call etsf_io_geometry_get(unitwff, geometry, lstat, error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   kpoints%reduced_coordinates_of_kpoints => hdr%kptns
   kpoints%kpoint_weights => hdr%wtk

   call etsf_io_kpoints_get(unitwff, kpoints, lstat, error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   hdr%fermie = fermie
   hdr%ecut   = ecut
   hdr%rprimd = rprimd
   hdr%znuclpsp(1:hdr%npsp) = hdr%znucltypat(1:hdr%npsp)

!  We get all other variables
!  In case the file is just ETSF valid, we ignore the missing variables and we use default values.

   hdr%date = 0
   call etsf_io_low_read_var(unitwff, "date", hdr%date, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%ecut_eff = hdr%ecut
   call etsf_io_low_read_var(unitwff, "ecut_eff", hdr%ecut_eff, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%ecutsm = zero
   call etsf_io_low_read_var(unitwff, "ecutsm", hdr%ecutsm, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%etot = zero 
   call etsf_io_low_read_var(unitwff, "etot", hdr%etot, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%intxc = 0
   call etsf_io_low_read_var(unitwff, "intxc", hdr%intxc, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%ixc = 1
   call etsf_io_low_read_var(unitwff, "ixc", hdr%ixc, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%occopt = 1
   call etsf_io_low_read_var(unitwff, "occopt", hdr%occopt, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%pertcase = 0
   call etsf_io_low_read_var(unitwff, "pertcase", hdr%pertcase, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%qptn(:) = 0
   call etsf_io_low_read_var(unitwff, "qptn", hdr%qptn, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%residm = zero 
   call etsf_io_low_read_var(unitwff, "residm", hdr%residm, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%stmbias = zero 
   call etsf_io_low_read_var(unitwff, "stmbias", hdr%stmbias, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%tphysel  = zero
   call etsf_io_low_read_var(unitwff, "tphysel", hdr%tphysel, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%tsmear = zero
   call etsf_io_low_read_var(unitwff, "tsmear", hdr%tsmear, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%ecutdg = hdr%ecut
   call etsf_io_low_read_var(unitwff, "ecutdg", hdr%ecutdg, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

!  test for old wavefunction style
   if(hdr%ecutsm>tol6 .and. headform<44 .and. &
&   .not.(fform==51.or.fform==52.or.fform==101.or.fform==102)) then
     write(message,'(a,es16.6,13a)' )&
&     '  The value of ecutsm is',hdr%ecutsm, &
&     ', while the file has been produced prior to v4.4 .',ch10,&
&     '  The definition of the smearing function has changed,', &
&     ' so that you are not allowed',ch10,&
&     '  to restart from a old wavefunction file. By contrast,', &
&     ' you can restart from an old',ch10,&
&     '  potential or density file, and perform a self-consistent', &
&     ' cycle with a new ABINIT version.',ch10,&
&     '  Action : produce a density or potential file using the old', &
&     ' version of ABINIT, and restart from it.'
     MSG_ERROR(message)
   end if

!  Multidimensional variables.
!  The case of istwfk is always 1, since ETSF don't use the time reversal symetry.
   hdr%istwfk(:)   = 1

   hdr%pspcod(:) = 0
   call etsf_io_low_read_var(unitwff, "pspcod", hdr%pspcod, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%pspdat(:) = 0
   call etsf_io_low_read_var(unitwff, "pspdat", hdr%pspdat, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%pspso(:) = 0
   call etsf_io_low_read_var(unitwff, "pspso", hdr%pspso, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%pspxc(:) = 0
   call etsf_io_low_read_var(unitwff, "pspxc", hdr%pspxc, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%so_psp(:) = 1
   call etsf_io_low_read_var(unitwff, "so_psp", hdr%so_psp, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%symafm(:) = 1
   call etsf_io_low_read_var(unitwff, "symafm", hdr%symafm, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%title(:) = ""
   call etsf_io_low_read_var(unitwff, "title", hdr%title, 132, lstat, error_data = error_data)
   if (lstat) then ! Pad the returned string with " " instead of "\0"
     do itypat = 1, size(hdr%title), 1
       call strip(hdr%title(itypat))
     end do
   end if

   call etsf_io_low_read_var(unitwff, "zionpsp", hdr%zionpsp, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   call etsf_io_low_read_var(unitwff, "znuclpsp", hdr%znuclpsp, lstat, error_data = error_data)
   ETSF_WARN(lstat,error_data)

   hdr%lmn_size = 1
   if (headform>=44) then ! Compared to 4.2, add lmn_size and
     call etsf_io_low_read_var(unitwff, "lmn_size", hdr%lmn_size, lstat, error_data = error_data)
!    DC: remove the PAW specific variables, they are not in accordance
!    with 59_io_mpi/hdr_io.F90. This part should be rewritten totally.
     if (hdr%usepaw==1) then
       write(message, '(12a)' ) ch10,&
&       ' hdr_io_etsf : ERROR -',ch10,&
&       '  The support for the internal variables of PAW are not yet',ch10,&
&       '  available with ETSF output. Restarting calculation from this',ch10,&
&       '  will not be possible.',ch10,&
&       '  Action : produce a density or potential file using the old',ch10,&
&       '  binary format of ABINIT, and restart from it.'
       call wrtout(std_out, message, 'COLL')
       call leave_new('COLL')
       do iatom=1,hdr%natom
         nullify(hdr%pawrhoij(iatom)%rhoijselect)
         nullify(hdr%pawrhoij(iatom)%rhoijp)
         hdr%pawrhoij(iatom)%ngrhoij = 0
         hdr%pawrhoij(iatom)%lmnmix_sz = 0
         hdr%pawrhoij(iatom)%use_rhoij_ = 0
         hdr%pawrhoij(iatom)%use_rhoijres = 0
       end do
!      !!    call etsf_io_low_read_var(unitwff, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
!      !!    allocate (rhoij(rhoijdim1,hdr%nspden,hdr%natom))
!      !!    call etsf_io_low_read_var(unitwff, "rhoij", rhoij, lstat, error_data = error_data)
!      !!    if (.not.lstat) goto 1000
!      !$
!      !!    cplex=1;if (rhoijdim1/=hdr%natom) cplex=2
!      !!    call rhoij_alloc(cplex,hdr%lmn_size,hdr%nspden,hdr%nspinor,hdr%nsppol,hdr%pawrhoij,hdr%typat)
!      !!    do iatom=1,hdr%natom
!      !!     itypat=hdr%typat(iatom)
!      !!     lmn2_size=hdr%lmn_size(itypat)*(hdr%lmn_size(itypat)+1)/2
!      !!     nselect=0
!      !!     if (cplex==1) then
!      !!      do ilmn=1,lmn2_size
!      !!       if (any(abs(rhoij(ilmn,:,iatom))>tol10)) then
!      !!        nselect=nselect+1
!      !!        hdr%pawrhoij(iatom)%rhoijselect(nselect)=ilmn
!      !!        do ispden=1,hdr%nspden
!      !!         hdr%pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoij(ilmn,ispden,iatom)
!      !!        end do
!      !!       end if
!      !!      end do
!      !!     else
!      !!      do ilmn=1,lmn2_size
!      !!       if (any(abs(rhoij(2*ilmn-1:2*ilmn,:,iatom))>tol10)) then
!      !!        nselect=nselect+1
!      !!        hdr%pawrhoij(iatom)%rhoijselect(nselect)=ilmn
!      !!        hdr%pawrhoij(iatom)%rhoijp(2*nselect-1,ispden)=rhoij(2*ilmn-1,ispden,iatom)
!      !!        hdr%pawrhoij(iatom)%rhoijp(2*nselect  ,ispden)=rhoij(2*ilmn  ,ispden,iatom)
!      !!       end if
!      !!      end do
!      !!     end if
!      !!     if (nselect<lmn2_size) then
!      !!      hdr%pawrhoij(iatom)%rhoijselect(nselect+1:lmn2_size)=0
!      !!      do ispden=1,hdr%nspden
!      !!       hdr%pawrhoij(iatom)%rhoijp(cplex*nselect+1:cplex*lmn2_size,ispden)=zero
!      !!      end do
!      !!     end if
!      !!     hdr%pawrhoij(iatom)%nrhoijsel=nselect
!      !!    end do
!      !!    deallocate(rhoij)
     end if
   end if

!  BigDFT private variables.
!  First implementation, we assume that the number of wavelet resolutions
!  is 2. Latter, we may add this value to hdr.

   lstat = .true.

!  -------------------------------------------------------------------------
!  Writing the header of an unformatted file
!  -------------------------------------------------------------------------
 else if(rdwr==2 .or. rdwr==6)then
!  We switch to write mode.
   call etsf_io_low_set_write_mode(unitwff, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

!  Associate and write values to ETSF groups.
   if (hdr%usewvl == 0) then  ! Plane wave case.
     ecut = hdr%ecut
     write(basis_set, "(A)") "plane_waves"
     basisdata%basis_set => basis_set
     basisdata%kinetic_energy_cutoff => ecut
     basisdata%number_of_coefficients => hdr%npwarr
   else  ! Wavelet case.
     write(basis_set, "(A)") "daubechies_wavelets"
     basisdata%basis_set => basis_set
!    Required variable than should enter the standard.
     call etsf_io_low_write_var(unitwff, "number_of_wavelets", hdr%nwvlarr, lstat, error_data = error_data)
     ETSF_CHECK_ERROR(lstat,error_data)
   end if

   call etsf_io_basisdata_put(unitwff, basisdata, lstat, error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   fermie = hdr%fermie
   electrons%fermi_energy => fermie
   electrons%number_of_states%data1D => hdr%nband
   electrons%occupations%data1D => hdr%occ

   call etsf_io_electrons_put(unitwff, electrons, lstat, error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   rprimd = hdr%rprimd
   geometry%primitive_vectors => rprimd
   geometry%reduced_symmetry_matrices => hdr%symrel
   geometry%atom_species => hdr%typat
   geometry%reduced_symmetry_translations => hdr%tnons
   geometry%reduced_atom_positions => hdr%xred
   geometry%atomic_numbers => hdr%znucltypat
   if (hdr%npsp == hdr%ntypat) then
     geometry%valence_charges => hdr%zionpsp
   end if

   call etsf_io_geometry_put(unitwff, geometry, lstat, error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   kpoints%reduced_coordinates_of_kpoints => hdr%kptns
   kpoints%kpoint_weights => hdr%wtk

   call etsf_io_kpoints_put(unitwff, kpoints, lstat, error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

!  Write non-ETSF variables.
   call etsf_io_low_write_var(unitwff, "date", hdr%date, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "codvsn", hdr%codvsn, 6, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "ecut_eff", hdr%ecut_eff, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "ecutsm", hdr%ecutsm, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "etot", hdr%etot, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "headform", 44, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "fform", fform, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "intxc", hdr%intxc, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "ixc", hdr%ixc, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "occopt", hdr%occopt, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "pertcase", hdr%pertcase, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "residm", hdr%residm, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "stmbias", hdr%stmbias, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "tphysel", hdr%tphysel, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "tsmear", hdr%tsmear, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

!  Version 44 add usepaw ecutdg
   call etsf_io_low_write_var(unitwff, "ecutdg", hdr%ecutdg, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "usepaw", hdr%usepaw, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

!  Array variables.
   call etsf_io_low_write_var(unitwff, "pspcod", hdr%pspcod, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "pspdat", hdr%pspdat, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "pspso", hdr%pspso, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "pspxc", hdr%pspxc, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "qptn", hdr%qptn, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "so_psp", hdr%so_psp, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "symafm", hdr%symafm, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "title", hdr%title, 132, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   call etsf_io_low_write_var(unitwff, "znuclpsp", hdr%znuclpsp, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

   if (hdr%npsp /= hdr%ntypat) then
     call etsf_io_low_write_var(unitwff, "zionpsp", hdr%zionpsp, lstat, error_data = error_data)
     ETSF_CHECK_ERROR(lstat,error_data)
   end if

!  Version 44 add lmn_size and rhoij
   call etsf_io_low_write_var(unitwff, "lmn_size", hdr%lmn_size, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

!  DC: remove the PAW specific variables, they are not in accordance
!  with 59_io_mpi/hdr_io.F90. This part should be rewritten totally.
   if (hdr%usepaw == 1) then
     write(message, '(12a)' ) ch10,&
&     ' hdr_io_etsf : WARNING -',ch10,&
&     '  The support for the internal variables of PAW are not yet',ch10,&
&     '  available with ETSF output. Restarting calculation from this',ch10,&
&     '  will not be possible.',ch10,&
&     '  Action : produce a density or potential file using the old',ch10,&
&     '  binary format of ABINIT, and restart from it.'
     call wrtout(std_out, message, 'COLL')
!    !!   rhoijdim1 = maxval(hdr%lmn_size)
!    !!   rhoijdim1 = hdr%pawrhoij(1)%cplex*rhoijdim1*(rhoijdim1+1)/2
!    !!   call etsf_io_low_write_var(unitwff, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
!    !!   allocate (rhoij(rhoijdim1,hdr%nspden,hdr%natom))
!    !!   do iatom=1,hdr%natom
!    !!    itypat=hdr%typat(iatom)
!    !!    lmn2_size = hdr%lmn_size(itypat)*(hdr%lmn_size(itypat)+1)/2
!    !!    cplex=hdr%pawrhoij(iatom)%cplex
!    !!    do ispden=1,hdr%nspden
!    !!     rhoij(1:cplex*lmn2_size,ispden,iatom)=zero
!    !!     if (cplex==1) then
!    !!      do irhoij=1,hdr%pawrhoij(iatom)%nrhoijsel
!    !!       ilmn=hdr%pawrhoij(iatom)%rhoijselect(irhoij)
!    !!       rhoij(ilmn,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(irhoij,ispden)
!    !!      end do
!    !!     else
!    !!      do irhoij=1,hdr%pawrhoij(iatom)%nrhoijsel
!    !!       ilmn=hdr%pawrhoij(iatom)%rhoijselect(irhoij)
!    !!       rhoij(2*ilmn-1,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(2*irhoij-1,ispden)
!    !!       rhoij(2*ilmn  ,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(2*irhoij  ,ispden)
!    !!      end do
!    !!     end if
!    !!    end do
!    !!   end do
!    !!   call etsf_io_low_write_var(unitwff, "rhoij", rhoij, lstat, error_data = error_data)
!    !!   if (.not.lstat) goto 1000
!    !!   deallocate (rhoij)
   end if
!  BigDFT variables.
   call etsf_io_low_write_var(unitwff, "usewvl", hdr%usewvl, lstat, error_data = error_data)
   ETSF_CHECK_ERROR(lstat,error_data)

 else if(rdwr==3 .or. rdwr==4)then

   call hdr_io_int(fform, hdr, rdwr, unitwff)
   lstat = .true.

 end if ! choice read/write/echo
#endif

 return
 fform=0 ; return   ! This is to allow treatment of old epsm1 format

end subroutine hdr_io_etsf
!!***
