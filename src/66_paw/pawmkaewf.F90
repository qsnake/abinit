!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmkaewf
!! NAME
!! pawmkaewf
!!
!! FUNCTION
!! Construct complete AE wave functions on the fine FFT grid adding onsite PAW corrections.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! dimcprj(natom)=array of dimensions of array cprj (not ordered)
!! natom=number of atoms in cell
!! ntypat=number of types of atoms in the cell
!! mpw=maximum dimensioned size of npw.
!! mband=maximum number of bands
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mkmem=number of k points which can fit in memory; set to 0 if use disk
!! nkpt=Total number of k-points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! unks=unit number for G vectors.
!! nband(nkpt*nsppol)=Number of bands for each k-point and spin.
!! istwfk(nkpt)=Storage mode at each k-point.
!! paral_kgb=Option for kgb parallelism
!! Pawfgrtab(natom) <type(pawfgrtab_type)> : data about the fine grid around each atom
!! Pawrad(ntypat) <type(pawrad_type)> : radial mesh data for each type of atom
!! Pawtab(ntypat) <type(pawtab_type)> : PAW functions around each type of atom
!! Psps <type(pseudopotential_type)> : basic pseudopotential data
!! Dtfil <type(datafiles_type)>=variables related to files
!! typat(natom) : list of atom types
!! cg(2,mcg)=planewave coefficients of wavefunctions.
!! Cprj(natom,nspinor*mband*mkmem*nsppol)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!   and each |p_lmn> non-local projector
!! Wffnow=struct info for current wf disk file
!! npwarr(nkpt)=Number of plane waves at each k-point
!! ngfftf(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  Note that ngfftf refers to the fine mesh.
!! kg(3,mpw*mkmem)=reduced planewave coordinates
!! eigen(mband*nkpt*nsppol)=eigenvalues (hartree) for all bands at each k point
!! occ(mband*nkpt*nsppol)=occupations for all bands at each k point
!! Hdr<hdr_type>=the header of wf, den and pot files
!! kpt(3,nkpt)=reduced coordinates of k points.
!!
!! OUTPUT
!!  ierr=Status error
!!  Main output is written on file (NETCDF file format).
!!
!! NOTES
!! In PAW calculations, the pseudized wavefunction us represented
!! on a relatively small plane wave basis set and is not normalized
!! as it does not include the on-site PAW contributions which is described
!! in terms of real spherical harmonics and radial functions.
!! For post-processing and proper visualization, it is necessary
!! to use the full electronic wave function, which is what this subroutine constructs.
!! Specifically, it computes the pseudo part by doing an FFT from G- to r-space
!! using the dense mesh defined by pawecutdg. The on-site PAW terms are also
!! computed in real space inside each sphere and added to the pseudo part.
!! Notice that this formula is expressed on the fine grid, and requires
!! interpolating the PAW radial functions onto this grid, as well as calling
!! initylmr in order to get the angular functions on the grid points.
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      abi_etsf_dims_init,abi_etsf_electrons_put,abi_etsf_geo_put,cprj_alloc
!!      cprj_diskinit_r,cprj_free,cprj_get,destroy_paw_pwaves_lmn
!!      etsf_io_data_init,etsf_io_electrons_put,etsf_io_low_close
!!      etsf_io_low_open_modify,etsf_io_low_set_write_mode,etsf_io_main_def
!!      etsf_io_main_put,flush_unit,fourwf,hdr_io_etsf,hdr_skip,ini_wf_etsf
!!      init_paw_pwaves_lmn,initmpi_seq,int2char4,nhatgrid,pawfgrtab_free
!!      pawfgrtab_init,pawfgrtab_print,printxsf,rdnpw,rwwf,sphereboundary
!!      wrap2_zero_one,wrtout,xbarrier_mpi,xcomm_init,xdefineoff,xmaster_init
!!      xmax_mpi,xme_init,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawmkaewf(Dtset,natom,mpw,mband,mcg,mcprj,nkpt,mkmem,nsppol,ntypat,nband,istwfk,npwarr,kpt,&
& paral_kgb,ngfftf,kg,dimcprj,Pawfgrtab,Pawrad,Pawtab,gmet,rprimd,ucvol,&
& Psps,Hdr,Dtfil,typat,eigen,occ,cg,Cprj,Wffnow,MPI_enreg,ierr,pseudo_norms,set_k,set_band)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes

 use m_xmpi
 use m_blas
 use m_splines
 use m_wffile
 use m_errors

 use m_io_tools,       only : get_unit, flush_unit
 use m_numeric_tools,  only : imax_loc
 use m_radmesh,        only : deducer0
 use m_crystal,        only : destroy_crystal, crystal_structure
 use m_crystal_io,     only : init_crystal_from_hdr
 use m_paw_toolbox,    only : pawfgrtab_init, pawfgrtab_free, pawfgrtab_print, &
&                             paw_pwaves_lmn_t, init_paw_pwaves_lmn, destroy_paw_pwaves_lmn

#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
 !use etsf_io_file,     only : etsf_io_file_merge
 use m_abi_etsf,       only : abi_etsf_dims_init
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawmkaewf'
 use interfaces_14_hidewrite
 use interfaces_27_toolbox_oop
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_59_io_mpi
 use interfaces_61_ionetcdf
 use interfaces_62_iowfdenpot
 use interfaces_66_paw, except_this_one => pawmkaewf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,mband,mcg,mcprj,mkmem,mpw,nsppol,paral_kgb,nkpt
 integer,intent(in),optional :: set_k,set_band
 integer,intent(out) :: ierr
 type(Datafiles_type),intent(in) :: Dtfil
 type(pseudopotential_type),intent(in) :: Psps
 type(MPI_type),intent(inout) :: MPI_enreg
 type(wffile_type),intent(inout) :: Wffnow
 type(hdr_type),intent(inout) :: Hdr
 type(dataset_type),intent(in) :: Dtset
!arrays
 integer,intent(in) :: typat(natom),nband(nkpt*nsppol),istwfk(nkpt),npwarr(nkpt),dimcprj(natom)
 integer,intent(in) :: ngfftf(18),kg(3,mpw*mkmem)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),target,intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),target,intent(in) :: occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: kpt(3,nkpt)
 real(dp),intent(in) :: gmet(3,3),rprimd(3,3),ucvol
 real(dp),optional,intent(out) :: pseudo_norms(nsppol,nkpt,mband)
 type(pawfgrtab_type),intent(in) :: Pawfgrtab(natom)
 type(pawrad_type),intent(in) :: Pawrad(ntypat)
 type(pawtab_type),intent(in) :: Pawtab(ntypat)
 type(cprj_type),intent(in) :: Cprj(natom,mcprj)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf0=0,tim_rwwf0=0
 integer :: bdtot_index,formeig,iband,icg,mgfftf,istat
 integer :: iatom,ifgd,ifftsph,ifft,itypat,ispinor,ipw,ndat,ii,i1,i2,i3
 integer :: jl,jm,jlmn
 integer :: nsp,max_nfgd,nfgd,ln_size,lmn_size,option
 integer :: iorder_cprj,spaceComm,rank,ibsp,ibg,isppol,ikpt,nband_k,cplex,master
 integer :: mcg_disk,n1,n2,n3,n4,n5,n6,ikg,npwout,istwf_k,npw_k
 integer :: indx,nfftot,my_spin,nprocs,tmp_unt
 integer :: optcut,optgr0,optgr1,optgr2,optrad,start_band,start_kpt,stop_kpt,stop_band
 real(dp),parameter :: weight1=one
 real(dp) :: phj,tphj,re_p,im_p,norm,norm_rerr,max_rerr,imur,reur,arg
 character(len=500) :: msg
 character(len=fnlen) :: xsf_fname
!arrays
 integer :: l_size_atm(natom)
 integer,allocatable :: gbound(:,:),kg_dum(:,:),kg_k(:,:)
 integer,allocatable :: atindx(:),atindx1(:),nattyp(:)
 integer,allocatable,target :: my_kpoints(:),my_spins(:)
 integer,allocatable :: my_kstable(:,:),my_nkpt(:)
 real(dp),allocatable :: r0shift(:,:,:),phk_atm(:,:,:)
 real(dp) :: red(3),shift(3)
 real(dp) :: rfft(3)
 real(dp) :: kpoint(3),cp_fact(2)
 real(dp),allocatable :: cg_disk(:,:),fofgin(:,:),fofgout(:,:)
 real(dp),allocatable :: eig_dum(:),occ_dum(:),denpot(:,:,:)
 real(dp),allocatable :: fofr(:,:,:,:),phkr(:,:)!,dummy_2d(:,:)
 real(dp),allocatable,target :: ur(:,:)
 real(dp),allocatable :: ur_ae_onsite(:,:),ur_ps_onsite(:,:)
 real(dp),allocatable :: ur_pw(:,:),ur_mask(:),xcart(:,:),dummy_1d(:)
 real(dp),allocatable :: rsph_red(:,:)
 real(dp),pointer :: rsph_cart(:,:)
 type(cprj_type),allocatable :: Cprj_k(:,:)
 type(pawfgrtab_type) :: local_pawfgrtab(natom)
 type(MPI_type) :: MPI_enreg_seq
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

#if defined HAVE_TRIO_ETSF_IO
 integer :: irank,fform,ncid,rdwr,var_main
 logical :: kdep,lstat
 character(len=4) :: tag
 character(len=80) :: file_title
 character(len=fnlen) :: my_fname,my_basename,out_file !,fname_part
 character(len=256),allocatable :: merge_files(:)
 type(etsf_dims) :: Dims
 type(etsf_groups_flags) :: Groups_flags
 !type(etsf_groups) :: Groups
 type(etsf_split) :: Split
 type(etsf_main) :: MainFolder
 type(etsf_io_low_error) :: Error_data
 type(etsf_electrons) :: Electrons_folder
 type(Wvl_wf_type) :: Dummy_wfs
#endif

! ************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 ABI_CHECK(paral_kgb==0,"paral_kgb/=0 not coded")
 ABI_CHECK(MPI_enreg%mode_para/='b',"mode_para=b not coded")
 ABI_CHECK(SIZE(dimcprj)>0,"dimcprj should be allocated")
 ABI_CHECK(mpi_enreg%paral_spin==0,"parallelisation over spinors not implemented")

 ABI_ALLOCATE(atindx,(natom))
 ABI_ALLOCATE(atindx1,(natom))
 ABI_ALLOCATE(nattyp,(ntypat))

 indx=1
 do itypat=1,ntypat
   nattyp(itypat)=0
   do iatom=1,natom
     if (typat(iatom)==itypat) then
       atindx (iatom )=indx
       atindx1(indx  )=iatom
       indx=indx+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do

 ABI_ALLOCATE(xcart,(3,Dtset%natom))
 call xredxcart(Dtset%natom,1,Hdr%rprimd,xcart,Hdr%xred)

!If collection of pseudo norms is enabled, make sure
!the array is initialised
 if (present(pseudo_norms)) pseudo_norms = zero

!use a local copy of pawfgrtab to make sure we use the correction in the paw spheres
!the usual pawfgrtab uses r_shape which may not be the same as r_paw
 do iatom = 1, natom
   l_size_atm(iatom) = pawtab(typat(iatom))%l_size
 end do
 call pawfgrtab_init(local_pawfgrtab,Pawfgrtab(1)%cplex,l_size_atm,Dtset%nspden)
 optcut = 1 ! use rpaw to construct local_pawfgrtab
 optgr0 = 0; optgr1 = 0; optgr2 = 0 ! dont need gY terms locally
 optrad = 1 ! do store r-R

 call initmpi_seq(MPI_enreg_seq)

 call nhatgrid(atindx1,gmet,MPI_enreg_seq,natom,natom,nattyp,ngfftf,ntypat,&
& optcut,optgr0,optgr1,optgr2,optrad,local_pawfgrtab,pawtab,rprimd,ucvol,Hdr%xred)
!now local_pawfgrtab is ready to use

 max_nfgd=MAXVAL(local_pawfgrtab(:)%nfgd) ! MAX no. of points in the fine grid for this PAW sphere
 ABI_ALLOCATE(r0shift,(3,max_nfgd,natom))
 ABI_ALLOCATE(phk_atm,(2,max_nfgd,natom))

 do iatom=1,natom
   nfgd=local_pawfgrtab(iatom)%nfgd ! no. of points in the fine grid for this PAW sphere
   ABI_ALLOCATE(rsph_red,(3,nfgd))
   ABI_ALLOCATE(rsph_cart,(3,nfgd))
!  rsph_cart => pawfgrtab(iatom)%rfgd
   do ifgd=1,nfgd
     rsph_cart(:,ifgd) = local_pawfgrtab(iatom)%rfgd(:,ifgd) + xcart(:,iatom)
   end do
   call xredxcart(nfgd,-1,rprimd,rsph_cart,rsph_red) ! we work in reduced coordinates.
   do ifgd=1,nfgd
     call wrap2_zero_one(rsph_red(1,ifgd),red(1),shift(1)) ! num = red + shift
     call wrap2_zero_one(rsph_red(2,ifgd),red(2),shift(2))
     call wrap2_zero_one(rsph_red(3,ifgd),red(3),shift(3))
     r0shift(:,ifgd,iatom) = shift
     if (ANY( ABS(shift) > tol12)) then
!      MSG_WARNING("rmR_red is outside the first unit cell.")
!      write(std_out,*)rsph_red(:,ifgd),shift
     end if
   end do
   ABI_DEALLOCATE(rsph_red)
   ABI_DEALLOCATE(rsph_cart)
 end do

 call pawfgrtab_print(local_pawfgrtab,unit=std_out,prtvol=Dtset%prtvol,mode_paral="COLL")

 ierr=0
#ifndef HAVE_TRIO_ETSF_IO
 ierr=-1
 write(msg,'(3a)')&
& " ETSF-IO support must be enabled in order to output AE PAW wavefunction. ",ch10,&
& " No output will be produced, use --enable-etsf-io at configure-time. "
 MSG_WARNING(msg)
 RETURN
!These statements are necessary to avoid the compiler complain about unused variables:
 ii=Dtset%usepaw;ii=Dtfil%unpaw;ii=Hdr%usepaw
#endif

!Init parallelism
 call xcomm_init(MPI_enreg,spaceComm)
 call xmaster_init(MPI_enreg,master)
 call xme_init(MPI_enreg,rank)
 nprocs = xcomm_size(spaceComm)

 ABI_CHECK(nprocs==1,"k spin parallelism not yet active")

 if (MPI_enreg%paral_compil_kpt==1)then
   call wrtout(std_out,'pawmkaewf: loop on k-points and spins done in parallel','COLL')
   call xbarrier_mpi(spaceComm)
 end if

 mgfftf=MAXVAL(ngfftf(1:3))

!=== Calculate my list of k-points and spin ===
!* my_kstable gives the sequential index for each k-point treated by rank.
!* cannot check for MPI_enreg%proc_distrb if nprocs ==1
 ABI_ALLOCATE(my_kstable,(nkpt,nsppol))
 my_kstable=0

 if (nprocs==1) then

   ii=0
   do ikpt=1,nkpt
     ii=ii+1; my_kstable(ikpt,:) = ii
   end do

   ABI_ALLOCATE(my_spins,(nsppol))
   do isppol=1,nsppol
     my_spins(isppol)=isppol
   end do

   ABI_ALLOCATE(my_kpoints,(nkpt))
   do ikpt=1,nkpt
     my_kpoints(ikpt) = ikpt
   end do

 else ! parallelism over k and spin.

   do isppol=1,nsppol
     ii=0
     do ikpt=1,nkpt
       nband_k = nband(ikpt+(isppol-1)*nkpt)
       if (ALL(MPI_enreg%proc_distrb(ikpt,1:nband_k,isppol)==rank)) then
         ii=ii+1
         my_kstable(ikpt,isppol)=ii
       end if
     end do
   end do

   ABI_ALLOCATE(my_nkpt,(nsppol))
   do isppol=1,nsppol
     my_nkpt(isppol) = COUNT(my_kstable(:,isppol)>0)
   end do

!  Each node has to deal with a single spin.
   if (nsppol>1 .and. ALL(my_nkpt>0)) then
     msg =' Non optimal distribution, some wave functions won''t be correctly initialized.'
     MSG_ERROR(msg)
   end if
   my_spin = imax_loc(my_nkpt)

   ABI_ALLOCATE(my_spins,(1))
   my_spins(1)=my_spin

   ABI_ALLOCATE(my_kpoints,(my_nkpt(my_spin)))
   ii=0
   do ikpt=1,nkpt
     if (my_kstable(ikpt,my_spin)/=0) then
       ii=ii+1
       my_kpoints(ii) = ikpt
     end if
   end do

   ABI_DEALLOCATE(my_nkpt)

 end if ! nprocs==1

#if defined HAVE_TRIO_ETSF_IO

!=== Initialize NETCDF files ===
 my_basename=dtfil%fnameabo_ae_wfk

!* For parallel case: the index of the processor must be appended.
!XG 100108 : One would better have it done in dtfil_init1 !!!!
 if (nprocs>1) then
   ABI_ALLOCATE(merge_files,(nprocs))
   do irank=1,nprocs
     call int2char4(irank,tag)
     merge_files(irank)=TRIM(Dtfil%filnam_ds(4))//'_P-'//tag//"_AE_WFK-etsf.nc"
     if (irank==rank+1) then
       my_basename=TRIM(Dtfil%filnam_ds(4))//'_P-'//tag//"_AE_WFK"
     end if
   end do
!  my_basename=TRIM(Dtfil%filnam_ds(4))//'_P-'//tag//"_AE_WFK"
!  my_basename=merge_files(rank+1)
 end  if

 my_fname=TRIM(my_basename)//"-etsf.nc"

 write(msg,'(2a)')' Opening file for AE PAW wave functions: ',TRIM(my_fname)
 call wrtout(std_out,msg,'PERS')
 call wrtout(ab_out,msg,'PERS')
 call xbarrier_mpi(spaceComm)

!Initialize Dims, remeber the hacking with with Dims
!call abi_etsf_dims_init(Dtset,2,Hdr%lmn_size,Psps,Dummy_wfs)
 call abi_etsf_dims_init(Dims,Dtset,2,Psps,Dummy_wfs)

!Change some values since we work in real space on the dense FFT mesh.
 Dims%max_number_of_coefficients      = etsf_no_dimension
 Dims%max_number_of_basis_grid_points = etsf_no_dimension
 Dims%number_of_localization_regions  = etsf_no_dimension
 Dims%real_or_complex_coefficients    = etsf_no_dimension

 Dims%real_or_complex_wavefunctions   = 2

 Dims%number_of_grid_points_vector1  = ngfftf(1)
 Dims%number_of_grid_points_vector2  = ngfftf(2)
 Dims%number_of_grid_points_vector3  = ngfftf(3)

!Dimensions for variables that can be splitted.
 Dims%my_number_of_kpoints = SIZE(my_kpoints) !etsf_no_dimension
 Dims%my_number_of_spins   = SIZE(my_spins)   !etsf_no_dimension

!Split data using k-points and spins
 if (nprocs>1) then
   Split%my_kpoints => my_kpoints(:)
   nullify(Split%my_spins)
   if (nsppol>1) Split%my_spins => my_spins(:)
 else
   nullify(Split%my_kpoints)
   nullify(Split%my_spins)
 end if

!=== Set-up the variables ===
!* These mandatory values are always written by the hdr_io_etsf() routine.
 Groups_flags%geometry  = etsf_geometry_all
 Groups_flags%electrons = etsf_electrons_all - etsf_electrons_x_functional - etsf_electrons_c_functional
 Groups_flags%kpoints   = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights

!These variables may be written depending on prt<something> input variables.
 Groups_flags%basisdata = etsf_basisdata_basis_set
 if (Dtset%usewvl==0) then
   Groups_flags%basisdata= Groups_flags%basisdata + &
&   etsf_basisdata_kin_cutoff + etsf_basisdata_n_coeff + etsf_basisdata_red_coord_pw
 else
   Groups_flags%basisdata= Groups_flags%basisdata + etsf_basisdata_coord_grid + etsf_basisdata_n_coeff_grid
 end if

!Groups_flags%basisdata = etsf_basisdata_none

!Wavefunctions in real space.
 Groups_flags%main = etsf_main_wfs_rsp

!=== Create the file ===
!* If the group contains main, we remove it for a while to be sure to
!add it at the end, after ABINIT private variables.
 var_main = Groups_flags%main
 Groups_flags%main = etsf_main_none
 write(file_title,'(a)')"PAW AE wavefunction given in real space"

 write(std_out,*)"Before  etsf_io_data_init"
 kdep=.TRUE.

 call etsf_io_data_init(my_fname,Groups_flags,Dims,file_title,'PAW AE Wavefunction File generated by ABINIT with ETSF_IO',&
& lstat,Error_data,k_dependent=kdep,overwrite=.TRUE.,Split_definition=Split)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 write(std_out,*)" my_number_of_kpoints ",Dims%my_number_of_kpoints
 write(std_out,*)" my_number_of_spins ",Dims%my_number_of_spins

!* Add the private ABINIT variables.
 call etsf_io_low_open_modify(ncid,my_fname,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call ini_wf_etsf(Dtset,Hdr%lmn_size,Psps%npsp,Psps%ntypat,ncid)

!Add the main part as last variables in the ETSF file.
 write(std_out,*)"Before  etsf_io_main_def"
 call etsf_io_main_def(ncid,lstat,Error_data,k_dependent=kdep,flags=var_main, Split=Split)
 ETSF_CHECK_MYERROR(lstat,Error_data)

!Close the file.
 call etsf_io_low_close(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

!Complete the geometry information with missing values from hdr_io().
 call abi_etsf_geo_put(Dtset,my_basename,Psps)

!To use the following statements, do not forget to declare:
!timrev(integer), Cryst(crystal_structure)
!timrev=2
!call init_crystal_from_hdr(Cryst,Hdr,timrev)
!call abi_crystal_put(Cryst,my_fname)
!call destroy_crystal(Cryst)

 call abi_etsf_electrons_put(Dtset,my_basename)

!We open again for further additions
 call etsf_io_low_open_modify(ncid,my_fname,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

!* Write the header.
!FIXME problem in the case of splitting over k-points due to SIZE mismatch
!in hdr%npwarr(number_of_kpoints) and number_of_coefficients(my_mkpt)

 fform=502; rdwr=2
 call hdr_io_etsf(fform,Hdr,rdwr,ncid)

 call etsf_io_low_close(ncid,lstat,Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

!=== Prepare the writing of the results ===
!
!1) Open file for writing
 call etsf_io_low_open_modify(ncid,my_fname,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

!2) Switch to write mode.
 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
#endif

!Init structure storing phi_{nlm} and tphi_(nlm} on the dense FFT points located in the PAW spheres.
 ABI_ALLOCATE(Paw_onsite,(natom))
 call init_paw_pwaves_lmn(Paw_onsite,natom,ntypat,typat,xcart,rprimd,Psps,Pawtab,Pawrad,local_pawfgrtab)

!FIXME check ordering in cprj and Eventually in external file
!why is iorder_cprj not stored in the file for crosschecking purpose?
!Here Im assuming cprj are not ordered!

 iorder_cprj=0
 call cprj_diskinit_r(atindx1,natom,iorder_cprj,mkmem,natom,0,dimcprj,dtset%nspinor,Dtfil%unpaw)

!mkmem==0 means wf and kg info on disk file
!* Skip header of Wffnow
!* Define offsets, in the case of MPI I/O.
 if (mkmem==0) then
   call hdr_skip(Wffnow,ierr)
   formeig=0
   call xdefineOff(formeig,Wffnow,MPI_enreg,nband,npwarr,Dtset%nspinor,nsppol,nkpt)
   mcg_disk=mpw*dtset%nspinor*mband
   ABI_ALLOCATE(cg_disk,(2,mcg_disk))
 end if

!n4,n5,n6 are FFT dimensions, modified to avoid cache trashing
 n1=ngfftf(1); n2=ngfftf(2); n3=ngfftf(3)
 n4=ngfftf(4); n5=ngfftf(5); n6=ngfftf(6)
 nfftot=PRODUCT(ngfftf(1:3))

!Loop over spin and k points
 bdtot_index=0; icg=0; ibg=0; norm_rerr=smallest_real
 do isppol=1,nsppol

   if (mkmem==0) rewind(Dtfil%unkg) ! Rewind the kpgsph data file on unit Dtfil%unkg
   ikg=0
   start_kpt=1
   stop_kpt=nkpt
!  Check if k-point was specified (only serial)
   if (present(set_k).AND.nprocs==1) then
     if (set_k/=0) then
       start_kpt = set_k
       stop_kpt = set_k
     end if
   end if

   do ikpt=start_kpt,stop_kpt

     kpoint  = kpt(:,ikpt)
     nband_k = nband(ikpt+(isppol-1)*nkpt)
     npw_k   = npwarr(ikpt)
     istwf_k = istwfk(ikpt)

     if (MPI_enreg%paral_compil_kpt==1)then
       if (MINVAL(ABS(MPI_enreg%proc_distrb(ikpt,1:nband_k,isppol)-rank))/=0) then
         bdtot_index=bdtot_index+nband_k; CYCLE
       end if
     end if

     ABI_ALLOCATE(phkr,(2,nfftot))
     do i3=0,n3-1
       rfft(3)=DBLE(i3)/n3
       do i2=0,n2-1
         rfft(2)=DBLE(i2)/n2
         do i1=0,n1-1
           rfft(1)=DBLE(i1)/n1
           ifft = 1 +i1 +i2*n1 +i3*n1*n2
           phkr(1,ifft) = COS(two_pi*dot_product(kpoint,rfft))
           phkr(2,ifft) = SIN(two_pi*dot_product(kpoint,rfft))
         end do
       end do
     end do
!    phkr(1,:)=one
!    phkr(2,:)=zero

!    Calculate the phase for the onsite PAW contributions.
     do iatom=1,natom
       nfgd=local_pawfgrtab(iatom)%nfgd ! no. of points in the fine grid for this PAW sphere
       do ifgd=1,nfgd
         arg = -two_pi* dot_product(r0shift(:,ifgd,iatom),kpoint)
         phk_atm(1,ifgd,iatom) = COS(arg)
         phk_atm(2,ifgd,iatom) = SIN(arg)
       end do
     end do

     ABI_ALLOCATE(Cprj_k,(natom,dtset%nspinor*nband_k))
     call cprj_alloc(Cprj_k,0,dimcprj)

!    Extract cprj for this k-point according to mkmem
     if (mkmem==0) then
       call cprj_get(atindx1,Cprj_k,Cprj,natom,1,ibg,ikpt,iorder_cprj,isppol,mband,mkmem,&
&       MPI_enreg,natom,nband_k,nband_k,dtset%nspinor,nsppol,Dtfil%unpaw)
     else
       ibsp=0
       do iband=1,nband_k
         do ispinor=1,dtset%nspinor
           ibsp=ibsp+1
           do iatom=1,natom
             Cprj_k(iatom,ibsp)%cp(:,:)=Cprj(iatom,ibsp+ibg)%cp(:,:)
           end do
         end do
       end do
     end if

     ABI_ALLOCATE(gbound,(2*mgfftf+8,2))
     ABI_ALLOCATE(kg_k,(3,npw_k))

!    Do i/o as needed
     if (mkmem==0) then
       nsp=dtset%nspinor
       call rdnpw(ikpt,isppol,nband_k,npw_k,nsp,0,Dtfil%unkg)
       read(Dtfil%unkg) kg_k(:,1:npw_k) ! Read k+g data
       call sphereboundary(gbound,istwf_k,kg_k,mgfftf,npw_k)
!      Read the wavefunction block for ikpt,isppol
       ABI_ALLOCATE(eig_dum,(mband))
       ABI_ALLOCATE(kg_dum,(3,0))
       ABI_ALLOCATE(occ_dum,(mband))
       call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,&
&       MPI_enreg,nband_k,nband_k,npw_k,dtset%nspinor,occ_dum,-2,0,tim_rwwf0,Wffnow)
       ABI_DEALLOCATE(eig_dum)
       ABI_DEALLOCATE(kg_dum)
       ABI_DEALLOCATE(occ_dum)
     else
       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       call sphereboundary(gbound,istwf_k,kg_k,mgfftf,npw_k)
     end if !mkmem==0

#if defined HAVE_TRIO_ETSF_IO
!    === Write eigenvalues and occupations ===
!    write(std_out,*)eigen(bdtot_index+1:bdtot_index+nband_k)
     Electrons_folder%eigenvalues__number_of_states = Dtset%mband
     Electrons_folder%eigenvalues%data1D         => eigen(bdtot_index+1:bdtot_index+nband_k)
     Electrons_folder%eigenvalues__kpoint_access = my_kstable(ikpt,isppol) !ikpt
     Electrons_folder%eigenvalues__spin_access   = isppol
     if (nprocs>1) Electrons_folder%eigenvalues__spin_access = 1

     Electrons_folder%occupations__number_of_states = Dtset%mband
     Electrons_folder%occupations%data1D         => occ(bdtot_index+1:bdtot_index+nband_k)
     Electrons_folder%occupations__kpoint_access = my_kstable(ikpt,isppol) !ikpt
     Electrons_folder%occupations__spin_access   = isppol
     if (nprocs>1) Electrons_folder%occupations__spin_access = 1

     write(std_out,*)"rank ",rank," about to write",my_kstable(ikpt,isppol)

     call etsf_io_electrons_put(ncid,Electrons_folder,lstat,Error_data)
     ETSF_CHECK_MYERROR(lstat,Error_data)
#endif

     start_band = 1
     stop_band = nband_k
!    If a single band is requested, neuter the loop (only serial)
     if (present(set_band).AND.nprocs==1) then
       if (set_band/=0) then
         start_band = set_band
         stop_band = set_band
       end if
     end if
     do iband=start_band,stop_band ! Loop over bands.

!      * Fourier transform on the real fft box of the smooth part.
       ndat=dtset%nspinor
       ABI_ALLOCATE(fofgin,(2,npw_k*ndat))
       ABI_ALLOCATE(fofr,(2,n4,n5,n6*ndat))

       if (mkmem/=0) then
         do ipw=1,npw_k*dtset%nspinor
           fofgin(:,ipw)=cg(:,ipw+(iband-1)*npw_k*dtset%nspinor+icg)
         end do
       else
         do ipw=1,npw_k*dtset%nspinor
           fofgin(:,ipw)=cg_disk(:,ipw+(iband-1)*npw_k*dtset%nspinor)
         end do
       end if

!      Complex can be set to 0 with this option(0) of fourwf
       option=0; cplex=0; npwout=1
       ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
       ABI_ALLOCATE(fofgout,(2,npwout*ndat))

       call fourwf(cplex,denpot,fofgin,fofgout,fofr,gbound,gbound,istwf_k,kg_k,kg_k,&
&       mgfftf,MPI_enreg,ndat,ngfftf,npw_k,npwout,n4,n5,n6,option,paral_kgb,tim_fourwf0,weight1,weight1,&
&       use_gpu_cuda=Dtset%use_gpu_cuda)
!      deallocate(fofgout)

!      Here I do not know if fourwf works in the case of spinors,
!      It seems that not all fftalg option support ndata! should check!
!      Do not forget to declare real(dp)::fofgin_down(:,:) to use the following statements
!      if (Dtset%nspinor==2) then
!      allocate(fofgin_down(2,npw_k))
!      fofgin_down(:,:)=fofgin(:,1+npw_k:2*npw_k)
!      ! Complex can be set to 0 with this option(0) of fourwf
!      cplex=1; option=1; npwout=1; ndat=1
!      call fourwf(cplex,denpot,fofgin_down,fofgout,fofr,gbound,gbound,istwf_k,kg_k,kg_k,&
!      &mgfftf,MPI_enreg,ndat,ngfftf,npw_k,npwout,n4,n5,n6,option,paral_kgb,tim_fourwf0,weight1,weight1)
!      deallocate(fofgin_down)
!      end if

       ABI_ALLOCATE(ur,(2,n1*n2*n3))
       ABI_ALLOCATE(ur_ae_onsite,(2,n1*n2*n3))
       ABI_ALLOCATE(ur_ps_onsite,(2,n1*n2*n3))
       ABI_ALLOCATE(ur_pw,(2,n1*n2*n3))
       ABI_ALLOCATE(ur_mask,(n1*n2*n3))
       ur=zero;ur_ae_onsite=zero;ur_ps_onsite=zero;ur_pw=zero;ur_mask=zero

!      * Add phase e^{ikr} since it is contained in cprj.
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             ii = i1 + n1*(i2-1)+ n1*n2*(i3-1)
             ur_pw(:,ii)=fofr(:,i1,i2,i3) ! Save pw part separately without the phase.
             ur(1,ii)= fofr(1,i1,i2,i3) * phkr(1,ii) - fofr(2,i1,i2,i3) * phkr(2,ii)
             ur(2,ii)= fofr(1,i1,i2,i3) * phkr(2,ii) + fofr(2,i1,i2,i3) * phkr(1,ii)
           end do
         end do
       end do
       ABI_DEALLOCATE(fofr)


!      === Add onsite term on the augmented FFT mesh ===
       do iatom=1,natom
         itypat  =typat(iatom)
         lmn_size=Pawtab(itypat)%lmn_size
         ln_size =Pawtab(itypat)%basis_size ! no. of nl elements in PAW basis
         nfgd    =local_pawfgrtab(iatom)%nfgd ! no. of points in the fine grid for this PAW sphere

         ibsp=(iband-1)*dtset%nspinor
         do ispinor=1,dtset%nspinor
           ibsp=ibsp+1
           do jlmn=1,lmn_size
             jl=Psps%indlmn(1,jlmn,itypat)
             jm=Psps%indlmn(2,jlmn,itypat)
             cp_fact(1) = Cprj_k(iatom,ibsp)%cp(1,jlmn) *sqrt(ucvol) ! Magic factor
             cp_fact(2) = Cprj_k(iatom,ibsp)%cp(2,jlmn) *sqrt(ucvol)

             do ifgd=1,nfgd ! loop over fine grid points in current PAW sphere.
               ifftsph = local_pawfgrtab(iatom)%ifftsph(ifgd) ! index of the point on the grid
               phj  = Paw_onsite(iatom)% phi(ifgd,jlmn)
               tphj = Paw_onsite(iatom)%tphi(ifgd,jlmn)
!              old code
!              re_p = cp_fact(1)
!              im_p = cp_fact(2)
!              apply the phase
               re_p = cp_fact(1) * phk_atm(1,ifgd,iatom) - cp_fact(2) * phk_atm(2,ifgd,iatom)
               im_p = cp_fact(1) * phk_atm(2,ifgd,iatom) + cp_fact(2) * phk_atm(1,ifgd,iatom)

               ur(1,ifftsph) = ur(1,ifftsph) + re_p * (phj-tphj)
               ur(2,ifftsph) = ur(2,ifftsph) + im_p * (phj-tphj)
               ur_ae_onsite(1,ifftsph) = ur_ae_onsite(1,ifftsph) + re_p * phj
               ur_ae_onsite(2,ifftsph) = ur_ae_onsite(2,ifftsph) + im_p * phj
               ur_ps_onsite(1,ifftsph) = ur_ps_onsite(1,ifftsph) + re_p * tphj
               ur_ps_onsite(2,ifftsph) = ur_ps_onsite(2,ifftsph) + im_p * tphj
               ur_mask(ifftsph) = one
             end do

           end do !jlmn
         end do !ispinor
       end do !iatom

!      * Remove the phase e^{ikr}, we store u(r).
#if 1
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             ii = i1 + n1*(i2-1)+ n1*n2*(i3-1)
             reur=ur(1,ii)
             imur=ur(2,ii)
             ur(1,ii)=  reur * phkr(1,ii) + imur * phkr(2,ii)
             ur(2,ii)= -reur * phkr(2,ii) + imur * phkr(1,ii)
             reur=ur_ae_onsite(1,ii)
             imur=ur_ae_onsite(2,ii)
             ur_ae_onsite(1,ii)=  reur * phkr(1,ii) + imur * phkr(2,ii)
             ur_ae_onsite(2,ii)= -reur * phkr(2,ii) + imur * phkr(1,ii)
             reur=ur_ps_onsite(1,ii)
             imur=ur_ps_onsite(2,ii)
             ur_ps_onsite(1,ii)=  reur * phkr(1,ii) + imur * phkr(2,ii)
             ur_ps_onsite(2,ii)= -reur * phkr(2,ii) + imur * phkr(1,ii)
           end do
         end do
       end do
#endif

       norm=zero
       do ifft=1,nfftot
         norm = norm + ur(1,ifft)**2+ur(2,ifft)**2
       end do
       norm=norm/nfftot

       norm_rerr = MAX((ABS(norm-one))*100,norm_rerr)
       write(std_out,*)"norm = g",norm
       call flush_unit(std_out)

!      MS: Various testing and debugging options
!      MG
       if (.TRUE..and.nprocs==1) then

!        Dump results to .xsf files if running in serial
         tmp_unt=get_unit()

         if (present(pseudo_norms)) then
!          Check the supposedly zero overlap |\tilde{Psi_n}-\tilde{Psi_n^1}|^2
           ABI_ALLOCATE(dummy_1d,(n1*n2*n3))
           dummy_1d=zero
           norm=zero
           do ifft=1,nfftot
             dummy_1d(ifft) = ((ur_pw(1,ifft)-ur_ps_onsite(1,ifft))**2 &
&             + (ur_pw(2,ifft)-ur_ps_onsite(2,ifft))**2) &
&             * ur_mask(ifft)
             norm = norm + dummy_1d(ifft)
           end do
           norm=norm/nfftot
           pseudo_norms(isppol,ikpt,iband) = norm
           if (Dtset%prtvol>9) then
             write(std_out,'(a,3(a,I0),a,F14.9,a)') ch10,' State sp',isppol,' kpt',ikpt,' bd',iband,&
&             ' |\tilde{Psi_n}-\tilde{Psi_n^1}|^2 norm:',norm,ch10
             write(xsf_fname,'(3(a,I0),a)') 'PAW_ps_norm_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
             open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
             call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&             Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
             close(tmp_unt)
           end if
           ABI_DEALLOCATE(dummy_1d)
         end if

!        if (present(delta_rhoij_norms)) then
!        Check the expression derived by G. Antonius:
!        Delta rho_ij =
!        (phi^t_i-psi^t_i)*(phi_j-phi^t_j)+(phi^t_j-psi^t_j)*(phi_i-phi^t_i)
!        allocate(dummy_2d(2,n1*n2*n3)); dummy_2d=zero
!        norm=zero
!        do ifft=1,nfftot
!        dummy_2d(:,ifft) = ((ur_pw(:,ifft)-ur_ps_onsite(1,ifft))**2 &
!        &             + (ur_pw(2,ifft)-ur_ps_onsite(2,ifft))**2) &
!        &             * ur_mask(ifft)
!        norm = norm + dummy_1d(ifft)
!        end do
!        norm=norm/nfftot
!        pseudo_norms(isppol,ikpt,iband) = norm
!        if (Dtset%prtvol>9) then
!        write(std_out,'(a,3(a,I0),a,F14.9,a)') ch10,' State sp',isppol,' kpt',ikpt,' bd',iband,&
!        &             ' |\tilde{Psi_n}-\tilde{Psi_n^1}|^2 norm:',norm,ch10
!        write(xsf_fname,'(3(a,I0),a)') 'PAW_ps_norm_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
!        open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
!        call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
!        &             Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
!        close(tmp_unt)
!        end if
!        deallocate(dummy_1d)
!        end if

         ABI_ALLOCATE(dummy_1d,(n1*n2*n3))
         dummy_1d=zero

         if (Dtset%prtvol>9) then
!          Onsite AE part
           dummy_1d = ur_ae_onsite(1,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_Re_ae_onsite_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)
           dummy_1d = ur_ae_onsite(2,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_Im_ae_onsite_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)
!          Calculate norm
           ur_ae_onsite(1,:) = ur_ae_onsite(1,:)**2 + ur_ae_onsite(2,:)**2
           dummy_1d = ur_ae_onsite(1,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_ABS_ae_onsite_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)


!          Onsite PS part-
           dummy_1d = ur_ps_onsite(1,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_Re_ps_onsite_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)
           dummy_1d = ur_ps_onsite(2,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_Im_ps_onsite_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)
!          Calculate norm
           ur_ps_onsite(1,:) = ur_ps_onsite(1,:)**2 + ur_ps_onsite(2,:)**2
           dummy_1d = ur_ps_onsite(1,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_ABS_ps_onsite_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)


!          PW part
           dummy_1d = ur_pw(1,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_Re_pw_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)
           dummy_1d = ur_pw(2,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_Im_pw_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)
!          Calculate norm
           ur_pw(1,:) = ur_pw(1,:)**2 + ur_pw(2,:)**2
           dummy_1d = ur_pw(1,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_ABS_pw_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)

!          Use ur_pw as dummy array
           ur_pw = ur

!          Full all-electron
           dummy_1d = ur_pw(1,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_Re_ae_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)
           dummy_1d = ur_pw(2,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_Im_ae_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)
!          Calculate norm
           ur_pw(1,:) = ur_pw(1,:)**2 + ur_pw(2,:)**2
           dummy_1d = ur_pw(1,1:nfftot)
           write(xsf_fname,'(3(a,I0),a)') 'PAW_ABS_ae_wfk_sp',isppol,'_kpt',ikpt,'_bd',iband,'.xsf'
           open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
           call printxsf(n1,n2,n3,dummy_1d,Hdr%rprimd,(/zero,zero,zero/),&
&           Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
           close(tmp_unt)

           if (iband==1.AND.ikpt==1) then ! This is printed only once
!            masking function - 1 in spheres, zero outside
             write(xsf_fname,'(a)') 'PAW_mask.xsf'
             open(tmp_unt,file=xsf_fname,status='unknown',form='formatted')
             call printxsf(n1,n2,n3,ur_mask,Hdr%rprimd,(/zero,zero,zero/),&
&             Dtset%natom,Dtset%ntypat,Dtset%typat,xcart,Dtset%znucl,tmp_unt,0)
             close(tmp_unt)
           end if

         end if ! prtvol>9

         ABI_DEALLOCATE(dummy_1d)

       else
         write(msg,'(5a)')&
&         " The option to print PAW all-electron wavefunctions is on, but execution ",ch10,&
&         " is in parallel on two or more processors. XcrysDen files with individual con-",ch10,&
&         " tributions will not be written. In order to enable this you must run in serial."
!        MG
         MSG_WARNING(msg)
       end if ! Check if serial run

#if defined HAVE_TRIO_ETSF_IO
       MainFolder%real_space_wavefunctions%data2D  => ur
       MainFolder%wfs_rsp__spin_access   =  isppol !this is wrong if para!
       if (nprocs>1) MainFolder%wfs_rsp__spin_access = 1
       MainFolder%wfs_rsp__kpoint_access = my_kstable(ikpt,isppol) !ikpt
       MainFolder%wfs_rsp__state_access  = iband
!      main_folder%wfs_coeff__number_of_coefficients = npw * dtset%nspinor

!      We use the group level write routine.
       call etsf_io_main_put(ncid,MainFolder,lstat,Error_data=Error_data)
       ETSF_CHECK_MYERROR(lstat,Error_data)
#endif

       ABI_DEALLOCATE(ur)
       ABI_DEALLOCATE(ur_ae_onsite)
       ABI_DEALLOCATE(ur_ps_onsite)
       ABI_DEALLOCATE(ur_pw)
       ABI_DEALLOCATE(ur_mask)

       ABI_DEALLOCATE(fofgin)
       ABI_DEALLOCATE(fofgout)
       ABI_DEALLOCATE(denpot)
       istat = ABI_ALLOC_STAT

     end do !nband_k

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
       ibg=ibg+dtset%nspinor*nband_k
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if

     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(phkr)

     call cprj_free(Cprj_k)
     ABI_DEALLOCATE(Cprj_k)

   end do !ikpt
 end do !nsppol

 if (mkmem==0)  then
   ABI_DEALLOCATE(cg_disk)
 end if

!* Free augmentation waves.
 call destroy_paw_pwaves_lmn(Paw_onsite)
 ABI_DEALLOCATE(Paw_onsite)

 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(nattyp)

 ABI_DEALLOCATE(my_kpoints)
 ABI_DEALLOCATE(my_spins)
 ABI_DEALLOCATE(my_kstable)

!* Maximum relative error over CPUs.
 call xmax_mpi(norm_rerr,max_rerr,spaceComm,ierr)
 write(std_out,*)"max_rerr=",max_rerr
 if (max_rerr>ten) then
   write(msg,'(7a)')&
&   " Inaccuracy on the normalization of the wave funtions exceeds 10%. ",ch10,&
&   " Likely due to the use of a too coarse FFT mesh or unconverged wavefunctions. ",ch10,&
&   " Numerical values inside the augmentation regions might be inaccurate. ",ch10,&
&   " Action: increase pawecutdg in your input file. "
   MSG_COMMENT(msg)
 end if

#if defined HAVE_TRIO_ETSF_IO
!=== Merge partial files ===
 if (nprocs>1) then
   call xbarrier_mpi(spaceComm)
   if (rank==master) then
     out_file="test_merge"
     write(msg,'(2a)')'Master node is merging NETCDF partial files into: ',TRIM(out_file)
     call wrtout(std_out, msg,'COLL')
!    call etsf_io_file_merge(out_file,merge_files,lstat,Error_data)
!    ETSF_CHECK_MYERROR(lstat,Error_data)
   end if
   call xbarrier_mpi(spaceComm)
 end if

 if (allocated(merge_files))  then
   ABI_DEALLOCATE(merge_files)
 end if
#endif

 ABI_DEALLOCATE(r0shift)
 ABI_DEALLOCATE(phk_atm)
 ABI_DEALLOCATE(xcart)
 call pawfgrtab_free(local_pawfgrtab)

 DBG_EXIT("COLL")

end subroutine pawmkaewf
!!***
