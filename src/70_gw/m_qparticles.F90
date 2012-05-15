!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_qparticles
!! NAME
!!  m_qparticles
!!
!! FUNCTION
!!  This module contains tools for the IO of the QP file and other procedures
!!  related to the calculation of the quasiparticle amplitudes represented in terms
!!  of KS states.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (FB, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_qparticles

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 !use m_xmpi
 use m_header
 use m_errors
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 use m_io_tools,       only : get_unit, file_exist, is_netcdf_file
 use m_numeric_tools,  only : linfit, c2r, set2unit
 use m_gwdefs,         only : sigma_parameters
 use m_crystal,        only : crystal_structure
 use m_crystal_io,     only : abi_crystal_put
 use m_bz_mesh,        only : bz_mesh_type
 use m_ebands,         only : get_valence_idx
 use m_sigma_results,  only : sigma_results

 implicit none

 private

 public :: wrqps           ! Write a QPS file.
 public :: rdqps           ! Read a QPS file.
 public :: show_QP         ! Report the components of a QP amplitude in terms of KS eigenstates.
 public :: rdgw            ! Read GW corrections from an external file.
 public :: updt_m_lda_to_qp! Updates the matrix of unitary transformation from
                           !   lda to qp states.

CONTAINS  !=========================================================================================================================
!!***

!!****f* m_qparticles/wrqps
!! NAME
!! wrqps
!!
!! FUNCTION
!!  Write the _QPS file containing information on the quasi-particles energies and wavefunctions.
!!
!! INPUTS
!!  fname=The name of the file
!!  Sigp<Sigma_parameters>=Parameters characterizing the self-energy calculation.
!!     %nsppol=1 for unpolarized, 2 for spin-polarized
!!     %nbnds=number of bands used for sigma
!!  Sr<Sigma_results>=Structure containing the results of the sigma run.
!!     %en_qp_diago(nbnds,nibz,nsppol)= NEW quasi-particle energies
!!     %eigvec_qp(nbnds,nbnds,nibz,nsppol)= NEW QP amplitudes in the KS basis set
!!      obtained by diagonalizing H0 + Herm(Sigma).
!!  m_lda_to_qp(nbnds,nbnds,nibz,nsppol)= expansion of the OLD QP amplitudes in terms of KS wavefunctions
!!  Kmesh<Bz_mesh_type>=information on the k-point sampling.
!!     %nibz=number of irreducible k-points
!!     %ibz(3,kibz)=reduced coordinates of the irreducible k-points
!!  nfftot=Total number of FFT points for density
!!  ngfftf(18)=Info on the FFT mesh for the density.
!!  nscf=Number of self consistent cycles performed
!!  nspden=number of spin-density components
!!  Cryst<crystal_structure>=Structure defining the crystal structure.
!!  Psps<type(pseudopotential_type)>=variables related to pseudopotentials.
!!  Pawrhoij(Cryst%natom*Psps%usepaw)<type(pawrhoij_type)>= rhoij datastructure.
!!  BSt<Bandstructure_type>=Structure containing the band structure energies (only used is nscf==-1)
!!
!! OUTPUT
!!  Only writing
!!
!! NOTES
!!  Old QPS fileformat:
!!   |
!!   | No. of QPSCF cycles already performed.
!!   | No. of k-points in the IBZ.
!!   | Total number of bands used to construct the Green's function (nbnds)
!!   | nsppol
!!   | For each spin and k-point in the IBZ:
!!   |   Reduced coordinates of the k-point.
!!   |   for each band:
!!   |     QP energies obtained by diagonalizing the QPSCGW Hamiltonian.
!!   |     <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}>$, ib=1,nbnds
!!   | FFT dimensions of the fine grid
!!   | QP density in real space.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrqps(fname,Sigp,Cryst,Kmesh,Psps,Pawtab,Pawrhoij,nspden,nscf,nfftot,ngfftf,Sr,Bst,m_lda_to_qp,rho_qp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrqps'
 use interfaces_14_hidewrite
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot,nscf,nspden
 character(len=fnlen),intent(in) :: fname
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Bandstructure_type),intent(in) :: BSt
 type(Sigma_parameters),intent(in) :: Sigp
 type(Sigma_results),intent(in) :: Sr
 type(Crystal_structure),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: rho_qp(nfftot,nspden)
 complex(dpc),intent(in) :: m_lda_to_qp(Sigp%nbnds,Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)
 type(Pawrhoij_type),intent(inout) :: Pawrhoij(Cryst%natom*Psps%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ib,ik,is,unqps,iatom,itypat
 character(len=500) :: msg
!arrays
 integer,allocatable :: nlmn_type(:)
 complex(dpc),allocatable :: mtmp(:,:)
#if defined HAVE_TRIO_ETSF_IO
! integer :: ncid, cplex
! logical :: lstat
! character(len=etsf_io_low_error_len) :: errmess
! character(len=etsf_charlen) :: file_title
! character(len=etsf_histlen) :: history
! character(len=fnlen) :: filetsf
! complex(dpc),allocatable :: full_mtmp(:,:,:,:)
! real(dp),allocatable :: rdata5(:,:,:,:,:)
! type(etsf_io_low_error) :: Error
! type(etsf_dims) :: Dims
 !type(etsf_groups_flags) :: Flags
#endif
! *************************************************************************

 DBG_ENTER("COLL")

 if (nscf>=0) then
  write(msg,'(3a)')ch10,' writing QP data on file : ',TRIM(fname)
  call wrtout(std_out,msg,'COLL')
  call wrtout(ab_out,msg,'COLL')
 end if


!#define MG_HAVE_TRIO_NETCDF_QPS 1
!#if ! defined MG_HAVE_TRIO_NETCDF_QPS

!if (.not.is_netcdf_file(fname)) then

 unqps=get_unit()
 open(unit=unqps,file=fname,form='formatted',status='unknown')

 write(unqps,*)nscf+1
 write(unqps,*)Kmesh%nibz
 write(unqps,*)Sigp%nbnds
 write(unqps,*)Sigp%nsppol

 ABI_ALLOCATE(mtmp,(Sigp%nbnds,Sigp%nbnds))

! call updt_m_lda_to_qp(Sigp,Kmesh,nscf,Sr,m_lda_to_qp) ! Calculate the new
!                                                       !  m_lda_to_qp

 if (nscf>=0) then ! Write the new m_lda_to_qp on file.
   do is=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       write(unqps,*)Kmesh%ibz(:,ik)
       !mtmp(:,:)=MATMUL(m_lda_to_qp(:,:,ik,is),Sr%eigvec_qp(:,:,ik,is))
       do ib=1,Sigp%nbnds
         write(unqps,*)Sr%en_qp_diago(ib,ik,is)
         !write(unqps,*)mtmp(:,ib)
         write(unqps,*)m_lda_to_qp(:,ib,ik,is)
       end do
     end do
   end do
 else if (nscf==-1) then ! Write fake QPS file with KS band structure (Mainly used for G0W)
   call set2unit(mtmp)
   do is=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       write(unqps,*)Kmesh%ibz(:,ik)
       do ib=1,Sigp%nbnds
         write(unqps,*)BSt%eig(ib,ik,is)
         write(unqps,*)mtmp(:,ib)
       end do
     end do
   end do
 else
   write(msg,'(a,i0)')" Wrong nscf= ",nscf
   MSG_ERROR(msg)
 end if

 ABI_DEALLOCATE(mtmp)
 !
 ! === Write FFT dimensions and QP density ===
 write(unqps,*)ngfftf(1:3)
 write(unqps,*)rho_qp(:,:)

 if (Psps%usepaw==1) then ! Write QP rhoij to be used for on-site density mixing.
   ABI_ALLOCATE(nlmn_type,(Cryst%ntypat))
   do itypat=1,Cryst%ntypat
     nlmn_type(itypat)=Pawtab(itypat)%lmn_size
   end do

   write(unqps,*) Cryst%natom, Cryst%ntypat
   write(unqps,*) (Cryst%typat(iatom), iatom=1,Cryst%natom)
   write(unqps,*) (nlmn_type(itypat), itypat=1,Cryst%ntypat)
   write(unqps,*) Pawrhoij(1)%nsppol, Pawrhoij(1)%nspden

   call rhoij_io(pawrhoij,unqps,Sigp%nsppol,Sigp%nspinor,nspden,nlmn_type,Cryst%typat,&
&                HDR_LATEST_HEADFORM,"Write",form="formatted")
   ABI_DEALLOCATE(nlmn_type)
 end if

 close(unqps)

#if 0
 ! Create the NetCDF file
 filetsf = "THIS_A_QP_TEST"
 file_title = "QPS file generated by ABINIT"; history = "BOO"

 call etsf_io_low_open_create(ncid, filetsf, etsf_file_format_version, lstat, file_title, &
&  history, Error, with_etsf_header=.TRUE., overwrite=.FALSE.)
 if (.not.lstat) goto 1000

 !call etsf_io_low_open_modify(ncid,filetsf,lstat,Error_data=Error)
 !if (.not.lstat) goto 1000

 ! =========================
 ! === Define dimensions ===
 ! =========================
 !call etsf_io_low_set_define_mode(ncid,lstat,Error)
 !if (.not.lstat) goto 1000

 Dims%max_number_of_angular_momenta  = Psps%mpsang
 Dims%max_number_of_projectors       = Psps%mproj !FIXME This does not work.
 Dims%max_number_of_states           = Sigp%nbnds
 Dims%number_of_atoms                = Cryst%natom
 Dims%number_of_atom_species         = Cryst%ntypat
 Dims%number_of_components           = nspden
 Dims%number_of_grid_points_vector1  = ngfftf(1)
 Dims%number_of_grid_points_vector2  = ngfftf(2)
 Dims%number_of_grid_points_vector3  = ngfftf(3)
 Dims%number_of_kpoints              = Kmesh%nibz
 !Dims%number_of_spinor_components    = nspinor  Dont think it is needed.
 Dims%number_of_spins                = Sigp%nsppol
 Dims%number_of_symmetry_operations  = Cryst%nsym
 Dims%real_or_complex_coefficients   = 2
 Dims%real_or_complex_density        = 1
 Dims%real_or_complex_gw_corrections = 2
 Dims%real_or_complex_wavefunctions  = 2

 call etsf_io_dims_def(ncid, Dims, lstat, Error)
 if (.not.lstat) goto 1000

 cplex=2
 call etsf_io_low_write_dim(ncid,'cplex',cplex,lstat,Error_data=Error)
 if (.not.lstat) goto 1000  ! Needed to store complex quantities

 !call etsf_io_low_write_dim(ncid,'max_number_of_scgw_iterations', NF90_UNLIMITED ,lstat,Error_data=Error)
 call etsf_io_low_write_dim(ncid,'max_number_of_scgw_iterations', 1 ,lstat,Error_data=Error)
 if (.not.lstat) goto 1000

 ! =======================
 ! == Define variables ===
 ! =======================

 call etsf_io_low_def_var(ncid,'usepaw',etsf_io_low_integer,lstat,Error_data=Error)
 if (.not.lstat) goto 1000

#if 0
 if (Psps%usepaw==1) then
 !version 44 add first dimension for rhoij = max(lmn_size)*(max(lmn_size)+1)/2
  rhoijdim1 = maxval(lmn_size)
  rhoijdim1 = rhoijdim1 * (rhoijdim1 + 1) / 2
 !impose rhoijdim1 >= 1 : if 0, it defaults to NF90_UNLIMITED
  rhoijdim1 = max(rhoijdim1, 1)
  call etsf_io_low_write_dim(unwff, "rhoijdim1", rhoijdim1, lstat, error_data = error_data)
  if (.not.lstat) goto 1000
 !Version 44 add lmn_size and rhoij
  call etsf_io_low_def_var(unwff, "lmn_size", etsf_io_low_integer, &
 & (/ "number_of_atom_species" /), lstat, error_data = error_data)
  if (.not.lstat) goto 1000
  call etsf_io_low_def_var(unwff, "rhoij", etsf_io_low_double, &
 & (/ pad("rhoijdim1"), pad("number_of_components"), &
 & pad("number_of_atoms") /), lstat, error_data = error_data)
  if (.not.lstat) goto 1000
 end if
#endif

 call etsf_io_low_def_var(ncid,'qp_density',etsf_io_low_double,&
& (/pad('real_or_complex_density'),&
&   pad('number_of_grid_points_vector1'),pad('number_of_grid_points_vector2'),pad('number_of_grid_points_vector3'),&
&   pad('number_of_components')/),lstat,Error_data=Error)
 if (.not.lstat) goto 1000

 call etsf_io_low_def_var(ncid,'m_lda_to_qp',etsf_io_low_double,&
& (/pad('cplex'),pad('max_number_of_states'),pad('max_number_of_states'),pad('number_of_kpoints'),pad('number_of_spins')/),&
& lstat,Error_data=Error)
 if (.not.lstat) goto 1000

 ! =====================
 ! === Start writing ===
 ! =====================

 ! Dump info on the crystal structure.
 !% call abi_crystal_put(Cryst,filetsf)

 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error)
 if (.not.lstat) goto 1000

 ! === Calculate the new m_lda_to_qp ===
 ABI_ALLOCATE(full_mtmp,(Sigp%nbnds,Sigp%nbnds,Kmesh%nibz,Sigp%nsppol))
 do is=1,Sigp%nsppol
  do ik=1,Kmesh%nibz
   full_mtmp(:,:,ik,is)=MATMUL(m_lda_to_qp(:,:,ik,is),Sr%eigvec_qp(:,:,ik,is))
  end do
 end do

 ABI_ALLOCATE(rdata5,(cplex,Sigp%nbnds,Sigp%nbnds,Kmesh%nibz,Sigp%nsppol))
 rdata5=c2r(full_mtmp)

 call etsf_io_low_write_var(ncid,'m_lda_to_qp',rdata5,lstat,Error_data=Error)
 if (.not.lstat) goto 1000

 ABI_DEALLOCATE(rdata5)
 ABI_DEALLOCATE(full_mtmp)

 call etsf_io_low_write_var(ncid,'qp_density',rho_qp,lstat,Error_data=Error)
 if (.not.lstat) goto 1000

 if (Psps%usepaw==1) then


 end if

 ! Close the file.
 call etsf_io_low_close(ncid,lstat,Error)
 if (.not.lstat) goto 1000

 1000 continue
 if (.not.lstat) then ! Handle the error.
  call etsf_io_low_error_to_str(errmess,Error)
  msg=errmess(1:MIN(500,LEN(errmess)))
  MSG_ERROR(msg)
 end if
#endif

 DBG_EXIT("COLL")

end subroutine wrqps
!!***

!----------------------------------------------------------------------

!!****f* m_qparticles/rdqps
!! NAME
!! rdqps
!!
!! FUNCTION
!!  Read a _QPS file containing the QP energies of the previous iteration, the coefficients
!!  defining the QP amplitudes in terms of the KS basis set and the QP density for mixing.
!!
!! INPUTS
!!  nfftot=Total number of FFT points for density
!!  ngfftf(18)=Info on the FFT mesh for the density.
!!  nspden=Number of SPin-DENsity components.
!!  usepaw=1 if we are using PAW.
!!  fname=Name of the file
!!  dimrho=1 if density has to be read, 0 otherwise
!!  BSt<Bandstructure_type>=Structure containing the initial band structure.
!!     %nsppol=1 for unpolarized, 2 for spin-polarized.
!!     %mband=Max number of bands used
!!     %nkpt=number of irreducible k-points.
!!     %kptns(3,nkpt)=reduced coordinates of each irreducible k-point.
!!  ucvol=Volume of the unit cell
!!
!! OUTPUT
!!  nbsc=number of bands used to describe the QP amplitudes
!!  nscf=number of iterations that have been performed (==0 if we start from a KS calculation)
!!  m_lda_to_qp(mband,mband,nibz,nsppol)=matrix giving the decomposition of the QP
!!   wavefunction in the mainfold generated by the KS wavefunctions
!!   (i.e. $ m_lda_to_qp(ib,jb,k,s) := <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}>$
!!  rhor_out(nfftot,nspden)=quasiparticle density
!!
!! SIDE EFFECTS
!!  BSt<Bandstructure_type>=Structure containing the initial band structure.
!!     %en_qp(mband,nkpt,nsppol)=QP energies at iteration nscf
!!
!! TODO
!!  The value of nspden is not reported in the QPS file thus we have a possible undetected error.
!!
!! PARENTS
!!      bethe_salpeter,mlwfovlp_qp,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine rdqps(BSt,fname,usepaw,nspden,dimrho,nscf,&
& nfftot,ngfftf,ucvol,paral_kgb,Cryst,Pawtab,MPI_enreg,nbsc,m_lda_to_qp,rhor_out,Pawrhoij)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rdqps'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot,nspden,usepaw,paral_kgb,dimrho
 integer,intent(out) :: nbsc,nscf
 real(dp),intent(in) :: ucvol
 character(len=fnlen),intent(in) :: fname
 type(crystal_structure),intent(in) :: Cryst
 type(Bandstructure_type),intent(inout) :: BSt
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(out) :: rhor_out(nfftot,nspden*dimrho)
 complex(dpc),intent(out) :: m_lda_to_qp(BSt%mband,BSt%mband,BSt%nkpt,BSt%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*usepaw)
 type(Pawrhoij_type),intent(inout) :: Pawrhoij(Cryst%natom*usepaw)

!Local variables-------------------------------
!scalars
 integer :: ib,ii,ik,isppol,nbandR,nkibzR,nsppolR,unqps,my_rank,ispden,master
 integer :: ifft,n1,n2,n3,ir1,ir2,ir3,ios
 integer :: cplex_fft,optin,optout,nfft_found
 integer :: iatom,natomR,nspdenR,ntypatR,itypat
 real(dp) :: uerr,rho_intp,nelect_qps,ratio
 logical,parameter :: use_FFT_interpolation=.TRUE.
 logical :: ltest
!logical :: lfile
 character(len=500) :: msg
!arrays
 integer :: ngfft_found(18)
 integer,allocatable :: typatR(:),nlmn_type(:)
 real(dp) :: kibz(3),rr(3),rhogdum(1,1)
 real(dp),allocatable :: en_tmp(:)
 real(dp),allocatable :: rhor_tmp(:,:)
 complex(dpc),allocatable :: mtmp(:,:),utest(:,:)

#if defined HAVE_TRIO_ETSF_IO
! integer :: cplex, ncid
! logical :: lstat
! character(len=etsf_io_low_error_len) :: errmess
!character(len=etsf_charlen) :: file_title
! character(len=etsf_histlen) :: history
! character(len=fnlen) :: filetsf
! real(dp),allocatable :: full_rmtmp(:,:,:,:,:)
! real(dp),allocatable ::  kptns(:,:)
! real(dp),allocatable :: rdata5(:,:,:,:,:)
! type(etsf_io_low_error) :: Error
! type(etsf_dims) :: Dims
 !type(etsf_groups_flags) :: Flags
#endif

! *************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(ALL(BSt%nband==BSt%nband(1)),"No. of bands must be constant")
 ABI_CHECK(dimrho==0.or.dimrho==1,'dimrho must be 0 or 1')

 ! This does not work in parallel !!?
 !% my_rank = xcomm_rank(MPI_enreg%spaceComm)
 call xme_init(MPI_enreg,my_rank)
 master=0
 !
 ! * Check whether file exists or not.
 write(msg,'(5a)')ch10,&
&  ' rdqps: reading QP wavefunctions of the previous step ',ch10,&
&  '        looking for file ',TRIM(fname)
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 if (.not.file_exist(fname)) then
   write(msg,'(2a)')' file not found, 1st iteration initialized with KS eigenelements ',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out, msg,'COLL')
   nscf=0; RETURN
 end if

 if (.not.is_netcdf_file(fname)) then
   unqps=get_unit()
   open(unit=unqps,file=fname,form='formatted',status='unknown')
   ! TODO the _QPS file should contain additional information

   read(unqps,*)nscf
   write(msg,'(a,i4,a)')' Number of iteration(s) already performed: ',nscf,ch10
   call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

   read(unqps,*)nkibzR
   if (nkibzR/=BSt%nkpt) then
     write(msg,'(2(a,i0))')&
&      ' Wrong number of k-points; Expected: ',BSt%nkpt,', Found: ',nkibzR
     MSG_ERROR(msg)
   end if

   read(unqps,*)nbandR
   nbsc=MIN(nbandR,BSt%mband)

   if (nbsc/=BSt%mband) then
     write(msg,'(3a,i4,a,i4)')&
&      ' QPS file contains less bands than that used in the present calculation ',ch10,&
&      ' Required: ',BSt%mband,', Found: ',nbandR
     MSG_WARNING(msg)
   end if

   if (nbsc/=nbandR) then
     write(msg,'(3a,i4,a)')&
&      ' The QPS file contains more bands than that used in the present calculation ',ch10,&
&      ' only the first ',nbandR,' bands will be read'
     MSG_COMMENT(msg)
   end if

   ABI_ALLOCATE(mtmp,(nbandR,nbandR))
   ABI_ALLOCATE(en_tmp,(nbandR))
   read(unqps,*)nsppolR

   ABI_CHECK(nsppolR==BSt%nsppol,"QPS generated with different nsppol")
   !
   ! === Read energies and transformation for each k-point and spin ===
   ! TODO: The format of the QPS file must be standardized !
   ! For example we might add the occupation numbers.
   do isppol=1,BSt%nsppol
     do ik=1,BSt%nkpt
       read(unqps,*)kibz(:)
       write(msg,'(a,i5,a,3(f6.3,1x),4x,a,i2)')' Reading ik ',ik,')  k = ',kibz(:),' is = ',isppol
       call wrtout(std_out,msg,'COLL')
       ltest=(ALL(ABS(kibz(:)-BSt%kptns(:,ik))<0.001))
       ABI_CHECK(ltest,'Wrong k-point read')
       do ib=1,nbandR
         read(unqps,*)en_tmp(ib)
         read(unqps,*)mtmp(:,ib)
       end do

       ! === Store transformation and update energies ===
       m_lda_to_qp(1:nbsc,1:nbsc,ik,isppol)=mtmp(1:nbsc,1:nbsc)
       BSt%eig(1:nbsc,ik,isppol)=en_tmp(1:nbsc)

       ! * Chech if matrix is unitary.
       ABI_ALLOCATE(utest,(nbsc,nbsc))
       utest=TRANSPOSE(mtmp(1:nbsc,1:nbsc)) !this is just for the buggy gfortran
       utest=MATMUL(CONJG(utest),mtmp(1:nbsc,1:nbsc))
       do ii=1,nbsc
         utest(ii,ii)=utest(ii,ii)-one
       end do
       uerr=MAXVAL(ABS(utest))
       if (uerr>tol6) then
         write(msg,'(a,es16.8)')' KS -> QP matrix is not unitary, MAX error = ',uerr
         MSG_WARNING(msg)
       end if
       ABI_DEALLOCATE(utest)
     end do !ik
   end do !isppol

   ABI_DEALLOCATE(mtmp)
   ABI_DEALLOCATE(en_tmp)
   !
   ! === Read the QP density ===
   ! * The two FFT grids might differ. In case perform an FFT interpolation to have rhor on the input mesh.
   if (dimrho==1) then
     read(unqps,*)n1,n2,n3

     if ( ALL(ngfftf(1:3)==(/n1,n2,n3/)) ) then
       read(unqps,*)rhor_out(:,:)
     else
       write(msg,'(2a,a,5(i3,a),i3)')&
&        ' FFT meshes differ. Performing Fourier interpolation. ',ch10,&
&        ' Found: ',n1,' x',n2,' x',n3,'; Expected: ',ngfftf(1),' x',ngfftf(2),' x',ngfftf(3)
       MSG_COMMENT(msg)

       ABI_ALLOCATE(rhor_tmp,(n1*n2*n3,nspden))
       read(unqps,*)rhor_tmp(:,:)

       if (use_FFT_interpolation) then
         ngfft_found(1:3)=(/n1,n2,n3/)
         ngfft_found(4)=2*(ngfft_found(1)/2)+1 ! 4:18 are not used, anyway!
         ngfft_found(5)=2*(ngfft_found(2)/2)+1
         ngfft_found(6)=ngfft_found(3)
         ngfft_found(7:18)=ngfftf(7:18)
         nfft_found=PRODUCT(ngfft_found(1:3)) !no FFT para

         cplex_fft =1 ! Real quantities.
         optin     =0 ! Input is taken from rhor.
         optout    =0 ! Output is only in real space.

         call fourier_interpol(cplex_fft,nspden,optin,optout,nfft_found,ngfft_found,nfftot,ngfftf,&
&          paral_kgb,MPI_enreg,rhor_tmp,rhor_out,rhogdum,rhogdum)
       else ! * Linear interpolation.
         do ispden=1,nspden
           do ir3=0,ngfftf(3)-1
             rr(3)=DBLE(ir3)/n3
             do ir2=0,ngfftf(2)-1
               rr(2)=DBLE(ir2)/n2
               do ir1=0,ngfftf(1)-1
                 rr(1)=DBLE(ir1)/n1
                 call interpol3d(rr,n1,n2,n3,rho_intp,rhor_tmp(:,ispden))
                 ifft = 1 +ir1 +ir2*ngfftf(1) +ir3*ngfftf(1)*ngfftf(2)
                 rhor_out(ifft,ispden)=rho_intp
               end do
             end do
           end do
         end do
       end if

       ABI_DEALLOCATE(rhor_tmp)
     end if
     !
     ! === Test the normalization of the QPS density ===
     ! * There might be errors due to the interpolation or the truncation of the G basis set
     ! * Density will be renormalized in the caller since for PAW we still have to add the onsite contribution.
     if (usepaw==0) then
       nelect_qps=SUM(rhor_out(:,1))*ucvol/nfftot; ratio=BSt%nelect/nelect_qps
       write(msg,'(3(a,f9.4))')&
&        ' Number of electrons calculated using the QPS density = ',nelect_qps,' Expected = ',BSt%nelect,' ratio = ',ratio
       call wrtout(std_out,msg,'COLL')
       !!rhor_out(:,:)=ratio*rhor_out(:,:)
     end if

     if (usepaw==1) then ! Write QP_rhoij for on-site density mixing.
       read(unqps,*,iostat=ios)natomR,ntypatR
       if (ios/=0) then
         msg="Old version of QPS file found. DO NOT USE rhoqpmix for this run."
         MSG_WARNING(msg)
         call wrtout(ab_out,msg,"COLL")
         ! Init dummy rhoij just to avoid problems in sigma when rhoij is freed.
         ABI_ALLOCATE(nlmn_type,(Cryst%ntypat))
         do itypat =1,Cryst%ntypat
           nlmn_type(itypat)=Pawtab(itypat)%lmn_size
         end do
         call rhoij_alloc(1,nlmn_type,nspden,BSt%nspinor,BSt%nsppol,Pawrhoij,Cryst%typat)
         ABI_DEALLOCATE(nlmn_type)
         close(unqps)
         RETURN
       end if

       ABI_CHECK(natomR ==Cryst%natom, "mismatch in natom")
       ABI_CHECK(ntypatR==Cryst%ntypat,"mismatch in ntypat")
       ABI_ALLOCATE(nlmn_type,(ntypatR))
       ABI_ALLOCATE(typatR,(ntypatR))

       read(unqps,*)(typatR(iatom), iatom=1,natomR)
       ABI_CHECK(ALL(Cryst%typat==typatR),"mismatch in typat")

       read(unqps,*)(nlmn_type(itypat), itypat=1,ntypatR)
       do itypat =1,Cryst%ntypat
         if (nlmn_type(itypat)/=Pawtab(itypat)%lmn_size) then
           MSG_ERROR("mismatch in nlmn_type, check QPS file")
         end if
       end do

       read(unqps,*) nsppolR,nspdenR
       ABI_CHECK(nsppolR==BSt%nsppol,"mismatch in nsppol")
       ABI_CHECK(nspdenR==nspden    ,"mismatch in nspden")

       call rhoij_io(pawrhoij,unqps,BSt%nsppol,BSt%nspinor,nspden,nlmn_type,Cryst%typat,&
&                    HDR_LATEST_HEADFORM,"Read",form="formatted")
       !% call rhoij_io(pawrhoij,std_out,BSt%nsppol,BSt%nspinor,nspden,nlmn_type,Cryst%typat,HDR_LATEST_HEADFORM,"Echo")

       ABI_DEALLOCATE(nlmn_type)
       ABI_DEALLOCATE(typatR)
     end if ! usepaw

   end if !dimrho=1

   close(unqps)

 else
!#if defined HAVE_TRIO_ETSF_IO
#if 0
  ! Open the file.
  call etsf_io_low_open_read(ncid,fname,lstat,Error_data=Error,with_etsf_header=.TRUE.)
  if (.not.lstat) goto 1000

  ! Read dimensions handled by ETSF.
  call etsf_io_dims_get(ncid,Dims,lstat,Error)
  if (.not.lstat) goto 1000

  !TODO read(unqps,*)nscf
  nkibzR = Dims%number_of_kpoints
  nbandR = Dims%max_number_of_states
  nsppolR = Dims%number_of_spins
  n1  = Dims%number_of_grid_points_vector1
  n2  = Dims%number_of_grid_points_vector2
  n3  = Dims%number_of_grid_points_vector3

  !call etsf_io_low_read_dim(ncid,'b1gw',Sr%b1gw,lstat,Error_data=Error)
  !if (.not.lstat) goto 1000

  ABI_ALLOCATE(kptns,(3,nkibzR))
  call etsf_io_low_read_var(ncid,'reduced_coordinates_of_kpoints',kptns,lstat,Error_data=Error)
  if (.not.lstat) goto 1000

  ltest=(ALL(ABS(kptns-BSt%kptns)<0.001))
  ABI_CHECK(ltest,'Wrong set of k-points read')
  ABI_DEALLOCATE(kptns)

  ABI_ALLOCATE(full_rmtmp,(2,nbandR,nbandR,nkibzR,nsppolR))
  call etsf_io_low_read_var(ncid,'m_lda_to_qp',full_rmtmp,lstat,Error_data=Error)
  if (.not.lstat) goto 1000

  !here be careful  !here be careful
  m_lda_to_qp(1:nbsc,1:nbsc,:,:)= DCMPLX(full_rmtmp(1,1:nbsc,1:nbsc,:,:),full_rmtmp(2,1:nbsc,1:nbsc,:,:))
  ABI_DEALLOCATE(full_rmtmp)

  ! Read energies.
  !% call etsf_io_low_read_var(ncid,'m_lda_to_qp',m_lda_to_qp,lstat,Error_data=Error)
  !% if (.not.lstat) goto 1000
  !% BSt%eig(1:nbsc,ik,isppol)=en_tmp(1:nbsc)

  if ( ALL(ngfftf(1:3)==(/n1,n2,n3/)) ) then
   call etsf_io_low_read_var(ncid,'qp_density',rhor_out,lstat,Error_data=Error)
   if (.not.lstat) goto 1000

  else
   write(msg,'(2a,a,5(i3,a),i3)')&
&   ' FFT meshes differ. Performing Fourier interpolation. ',ch10,&
&   ' Found: ',n1,' x',n2,' x',n3,'; Expected: ',ngfftf(1),' x',ngfftf(2),' x',ngfftf(3)
   MSG_COMMENT(msg)

   ABI_ALLOCATE(rhor_tmp,(n1*n2*n3,nspden))
   call etsf_io_low_read_var(ncid,'qp_density',rhor_tmp,lstat,Error_data=Error)
   if (.not.lstat) goto 1000

   if (use_FFT_interpolation) then
    ngfft_found(1:3)=(/n1,n2,n3/)
    ngfft_found(4)=2*(ngfft_found(1)/2)+1 ! 4:18 are not used, anyway!
    ngfft_found(5)=2*(ngfft_found(2)/2)+1
    ngfft_found(6)=ngfft_found(3)
    ngfft_found(7:18)=ngfftf(7:18)
    nfft_found=PRODUCT(ngfft_found(1:3)) !no FFT para

    cplex =1 ! Real quantities.
    optin =0 ! Input is taken from rhor.
    optout=0 ! Output is only in real space.

    call fourier_interpol(cplex,nspden,optin,optout,nfft_found,ngfft_found,nfftot,ngfftf,&
&    paral_kgb,MPI_enreg,rhor_tmp,rhor_out,rhogdum,rhogdum)

   else ! * Linear interpolation.
    do ispden=1,nspden
     do ir3=0,ngfftf(3)-1
      rr(3)=DBLE(ir3)/n3
      do ir2=0,ngfftf(2)-1
       rr(2)=DBLE(ir2)/n2
       do ir1=0,ngfftf(1)-1
        rr(1)=DBLE(ir1)/n1
        call interpol3d(rr,n1,n2,n3,rho_intp,rhor_tmp(:,ispden))
        ifft = 1 +ir1 +ir2*ngfftf(1) +ir3*ngfftf(1)*ngfftf(2)
        rhor_out(ifft,ispden)=rho_intp
       end do
      end do
     end do
    end do
   end if ! do FFT interpolation

   ABI_DEALLOCATE(rhor_tmp)
  end if

  ! Close the file.
  call etsf_io_low_close(ncid,lstat,Error)
  if (.not.lstat) goto 1000

1000 continue
  if (.not.lstat) then ! Handle the error.
   call etsf_io_low_error_to_str(errmess,Error)
   MSG_ERROR(errmess)
  end if

#else
  MSG_ERROR("ETSF-IO support is missing")
#endif
 end if

 DBG_EXIT("COLL")

end subroutine rdqps
!!***

!----------------------------------------------------------------------

!!****f* m_qparticles/show_QP
!! NAME
!! show_QP
!!
!! FUNCTION
!! Print in a nice format (?) the expansion coefficients of the quasiparticle
!! amplitudes in terms of the KS eigenvectors
!!
!! INPUTS
!!  Bst<Bandstructure_type>=Description of the band structure.
!!    %nsppol=1 for unpolarized, 2 for spin-polarized.
!!    %mband=Max number of bands (in GW doesn"t depend on k an spin)
!!    %nkpt=number of irreducible k-points.
!!    %eig(mband,nkpt,nsppol)= QP energies for each k-point, band and spin.
!!  m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)=matrix giving the decomposition of the QP
!!   amplitued in the mainfold generated by the KS wavefunctions
!!   (i.e $ m_lda_to_qp(ib,jb,k,s) := \langle \psi_{ib,k,s}^{KS}| \psi_{jb,k,s}^{QP}\rangle $
!!  fromb,tob=initial and final band index for QP, only states in this range are printed
!!  prtvol=Verbosity level (not used)
!!  unit=Unit number of the output file
!! tolmat[Optional]=Only components whose coefficient has modulus larger than tolmat are shown (default is 0.01)
!!
!! OUTPUT
!!  Only printing
!!
!! NOTES
!!  Only master node should call this routine.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine show_QP(Bst,m_lda_to_qp,fromb,tob,unit,prtvol,tolmat,kmask)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'show_QP'
 use interfaces_27_toolbox_oop
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: fromb,tob
 integer,optional,intent(in) :: prtvol,unit
 real(dp),optional,intent(in) :: tolmat
 type(Bandstructure_type),intent(in) :: Bst
!arrays
 logical,optional,intent(in) :: kmask(Bst%nkpt)
 complex(dpc),intent(in) :: m_lda_to_qp(Bst%mband,Bst%mband,Bst%nkpt,Bst%nsppol)

!Local variables-------------------------------
!scalars
 integer :: nbra=5
 integer :: ib_start,ib_stop,my_prtvol
 integer :: counter,ib_KS,ib_QP,ikibz,isp,nspace,my_unt,nband_k
 real(dp) :: my_tolmat
 character(len=10) :: bks,bqp,k_tag,spin_tag
 character(len=500) :: KS_row,KS_ket,tmpstr,QP_ket

! *********************************************************************

 my_unt   =std_out  ; if (PRESENT(unit  )) my_unt   =unit
 my_prtvol=0        ; if (PRESENT(prtvol)) my_prtvol=prtvol
 ib_start =1        ; if (PRESENT(fromb )) ib_start =fromb
 ib_stop  =Bst%mband; if (PRESENT(tob   )) ib_stop  =tob
 my_tolmat=0.001    ; if (PRESENT(tolmat)) my_tolmat=ABS(tolmat)

 ! * I suppose nband_k is constant thus the check is done here.
 if (ib_start<=0       ) ib_start=1
 if (ib_start>Bst%mband) ib_start=Bst%mband
 if (ib_stop<=0        ) ib_stop=1
 if (ib_stop>Bst%mband ) ib_stop=Bst%mband

 ! Have to follow rules 7.f.
 write(my_unt,'(/,a,/,a,/,a,f5.3,a,/,a)')&
   ' '//REPEAT('*',76),&
&  ' ***** QP amplitudes expressed as linear combination of KS eigenstates. *****',&
&  ' ***** Only KS components whose modulus is larger than ',my_tolmat,' are shown  ***** ',&
&  ' '//REPEAT('*',76)

 if (PRESENT(kmask)) then
   if (.not.ALL(kmask)) write(my_unt,'(/,a,i3,a)')' Only ',COUNT(kmask),' k-points are reported '
 end if

 do isp=1,Bst%nsppol
   call int2char(isp,spin_tag)
   write(my_unt,'(/,a,i2,a,/)')' >>>>> Begin block for spin ',isp,' <<<<< '

   do ikibz=1,Bst%nkpt
     if (PRESENT(kmask)) then
       if (.not.kmask(ikibz)) CYCLE
     end if
     call int2char(ikibz,k_tag)
     nband_k=Bst%nband(ikibz+(isp-1)*Bst%nkpt)
     write(my_unt,'(a,i4,a,3es16.8,a,f6.3,/)')' k-point: ',ikibz,') ',Bst%kptns(:,ikibz),'; wtk= ',Bst%wtk(ikibz)

     do ib_QP=ib_start,ib_stop
       call int2char(ib_QP,bqp)
       QP_ket=' |QP: b='//TRIM(bqp)//'; s='//TRIM(spin_tag)//'> = '
       write(my_unt,'(a)')TRIM(QP_ket)
       nspace=LEN(TRIM(QP_ket))

       counter=0 ; KS_row=REPEAT('',nspace+2)
       do ib_KS=1,Bst%mband
         if (ABS(m_lda_to_qp(ib_KS,ib_QP,ikibz,isp))<my_tolmat) CYCLE
         counter=counter+1
         call int2char(ib_KS,bks)
         write(tmpstr,'(3a)')' |',TRIM(bks),'>'
         write(KS_ket,'(1x,2f7.3,a,1x)')m_lda_to_qp(ib_KS,ib_QP,ikibz,isp),TRIM(tmpstr)
         KS_row=TRIM(KS_row)//TRIM(KS_ket)
         if (MOD(counter,nbra)==0) then  ! nbra KS kets per row
           write(my_unt,'(a)')TRIM(KS_row)
           KS_row=REPEAT('',nspace+2)
         end if
       end do

       if (MOD(counter,nbra)/=0) write(my_unt,'(a)')TRIM(KS_row) ! Last row, if any
       write(my_unt,'(a)')''
     end do !ib_QP

   end do !ikibz
 end do !isp

 write(my_unt,'(a,/)')' '//REPEAT('*',76)

end subroutine show_QP
!!***

!----------------------------------------------------------------------

!!****f* m_qparticles/rdgw
!! NAME
!! rdgw
!!
!! FUNCTION
!!  This subroutine reads the GW corrections from a _GW file.
!!
!! INPUTS
!!  [extrapolate]= if .TRUE., the routine extrapolates the
!!    GW corrections for the states that have not been explicitly evaluated (default).
!!    If .FALSE., only the GW states that have been calculated will be used to replace
!!    the input eigenvalues stored in Bst%eig
!!  Bst<Bandstructure_type>=type describing the Band structure.
!!    %nbnds=number of bands.
!!    %nkpt=number of irred k-points.
!!    %nsppol=number of spin
!!    %kptns(3,nkpt)=irreducible k-points
!!
!! SIDE EFFECTS
!!   Bst%eig(%mband,%nkpt,%nsppol)=Overwritten with GW energies according to extrapolate flag.
!!
!! OUTPUT
!!   igwene(Bst%mband,Bst%nkpt,Bst%nsppol)= The imaginary part of the QP energies.
!!
!! PARENTS
!!      mlwfovlp_qp,screening,setup_bse,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine rdgw(Bst,fname,igwene,extrapolate)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rdgw'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: fname
 type(Bandstructure_type),intent(inout) :: Bst
 logical,optional,intent(in) :: extrapolate
!arrays
 real(dp),intent(out) :: igwene(Bst%mband,Bst%nkpt,Bst%nsppol)

!Local variables ------------------------------
!scalars
 integer :: ib,ibr,ik,ikibz,ikr,ios,is,nn,nbandR,nkibzR,nsppolR,unt,nbv
 real(dp) :: alpha,beta,degw,egw_r,egw_i,smrt
 logical :: do_extrapolate
 character(len=500) :: msg
!arrays
 integer,allocatable :: vbik(:,:),seen(:)
 real(dp) :: kread(3)
 real(dp),allocatable :: gwcorr(:,:,:)

!************************************************************************

 ABI_CHECK(ALL(Bst%nband==Bst%mband),"nband must be constant")

 call wrtout(std_out,' Reading GW corrections from file '//TRIM(fname),'COLL')

 unt=get_unit()
 open(unt,file=fname,status='old',iostat=ios)
 if (ios/=0) then
   write(msg,'(3a)')' Opening file: ',TRIM(fname),' as old'
   MSG_ERROR(msg)
 end if

 read(unt,*)nkibzR,nsppolR

 ABI_CHECK(nsppolR==Bst%nsppol,"mismatch in nsppol")
 if (nkibzR/=Bst%nkpt) then
   write(msg,'(a,i4,a,i4,2a)')&
&   ' Found less k-points than that required ',nkibzR,'/',Bst%nkpt,ch10,&
&   ' Some k-points will be skipped. Continuing anyway '
   MSG_WARNING(msg)
 end if

 ABI_ALLOCATE(gwcorr,(Bst%mband,Bst%nkpt,Bst%nsppol))
 ABI_ALLOCATE(seen,(Bst%nkpt))
 gwcorr=zero
 igwene=zero

 do is=1,Bst%nsppol
   seen=0

   do ikr=1,nkibzR
     read(unt,*)kread(:)
     read(unt,*)nbandR
     ikibz=0
     do ik=1,Bst%nkpt
       if (ALL(ABS(kread(:)-Bst%kptns(:,ik))<0.0001)) then
         ikibz=ik
         seen(ik) = seen(ik) + 1
       end if
     end do
     do ib=1,nbandR
       read(unt,*)ibr,egw_r,degw,egw_i
       if (ibr<=Bst%mband .and. ikibz/=0) then
         gwcorr(ibr,ikibz,is)=degw/Ha_eV
         igwene(ibr,ikibz,is)=egw_i/Ha_eV
       end if
     end do
   end do

   if (ANY(seen/=1)) then
     do ik=1,Bst%nkpt
       if (seen(ik)/=1) then
         write(msg,'(a,3f8.3,a)')" k-point: ",Bst%kptns(:,ik)," not found in the GW file!"
         MSG_WARNING(msg)
       end if
     end do
   end if

 end do

 ABI_DEALLOCATE(seen)
 close(unt)

 do_extrapolate=.TRUE.; if (PRESENT(extrapolate)) do_extrapolate=extrapolate

 if (.not. do_extrapolate) then ! Only the bands calculated are updated.
   Bst%eig = Bst%eig + gwcorr

 else

   if (ANY(ABS(igwene)>tol6)) then
     write(msg,'(4a)')ch10,&
&      " The GW file contains QP energies with non-zero imaginary part",ch10,&
&      " Extrapolation not coded, change the source! "
     MSG_ERROR(msg)
   end if

   ABI_ALLOCATE(vbik,(BSt%nkpt,BSt%nsppol))
   vbik(:,:) = get_valence_idx(BSt)

   do is=1,Bst%nsppol
     do ik=1,Bst%nkpt

      nbv=vbik(ik,is) ! Index of the (valence band| Fermi band) for each spin
      nn=Bst%mband-nbv

      do ib=nbv+1,Bst%mband
        if ( ABS(gwcorr(ib,ik,is)) < tol16) then
          nn=ib-1-nbv
          if (nn>1) then
            call wrtout(std_out," Linear extrapolating (conduction) GW corrections beyond the read values","COLL")
            smrt=linfit(nn,Bst%eig(nbv+1:nbv+nn,ik,is),gwcorr(nbv+1:nbv+nn,ik,is),alpha,beta)
          else
            call wrtout(std_out," Assuming constant (conduction) GW corrections beyond the read values",'COLL')
            alpha=zero
            beta =gwcorr(nbv+nn,ik,is)
          end if
          EXIT !ib loop
        end if
      end do !ib

      do ib=nbv+nn+1,Bst%mband
        gwcorr(ib,ik,is)= alpha*Bst%eig(ib,ik,is) + beta
      end do

      nn=nbv
      do ib=nbv,1,-1
        if ( ABS(gwcorr(ib,ik,is)) < tol16) then
         nn=nbv-ib
         if (nn>1) then
           call wrtout(std_out,"Linear extrapolating (valence) GW corrections beyond the read values","COLL")
           smrt=linfit(nn,Bst%eig(nbv-nn+1:nbv,ik,is),gwcorr(nbv-nn+1:nbv,ik,is),alpha,beta)
         else
           call wrtout(std_out,"Assuming constant (valence) GW corrections beyond the read values","COLL")
           alpha=zero
           beta =gwcorr(nbv,ik,is)
         end if
         EXIT !ib
        end if
      end do !ib

      do ib=1,nbv-nn
        gwcorr(ib,ik,is)=alpha*Bst%eig(ib,ik,is) + beta
      end do

     end do !ik
   end do !is

   call wrtout(std_out,' k  s     GW corrections [eV] ','COLL')
   do is=1,Bst%nsppol
     do ik=1,Bst%nkpt
       write(msg,'(i3,1x,i3,10f7.2/50(10x,10f7.2/))')ik,is,(Ha_eV*gwcorr(ib,ik,is),ib=1,Bst%mband)
       call wrtout(std_out,msg,"COLL")
     end do
   end do
   Bst%eig = Bst%eig + gwcorr
   ABI_DEALLOCATE(vbik)
 end if

 call wrtout(std_out,' k   s    GW eigenvalues [eV]',"COLL")
 do is=1,Bst%nsppol
   do ik=1,Bst%nkpt
     write(std_out,'(2(i3,1x),7x,10f7.2/50(15x,10f7.2/))')ik,is,(Ha_eV*Bst%eig(ib,ik,is),ib=1,Bst%mband)
   end do
 end do

 ABI_DEALLOCATE(gwcorr)

end subroutine rdgw
!!***

!----------------------------------------------------------------------

!!****f* m_qparticles/updt_m_lda_to_qp
!! NAME
!! updt_m_lda_to_qp
!!
!! FUNCTION
!! Updates the matrix containing the unitary transformation from the lda states
!! to the quasiparticle states.
!!
!! INPUTS
!!  Sigp<Sigma_parameters>=Parameters characterizing the self-energy calculation.
!!     %nsppol=1 for unpolarized, 2 for spin-polarized
!!     %nbnds=number of bands used for sigma
!!  Sr<Sigma_results>=Structure containing the results of the sigma run.
!!     %en_qp_diago(nbnds,nibz,nsppol)= NEW quasi-particle energies
!!     %eigvec_qp(nbnds,nbnds,nibz,nsppol)= NEW QP amplitudes in the KS basis set
!!      obtained by diagonalizing H0 + Herm(Sigma).
!!  Kmesh<Bz_mesh_type>=information on the k-point sampling.
!!     %nibz=number of irreducible k-points
!!     %ibz(3,kibz)=reduced coordinates of the irreducible k-points
!!  nscf=Number of self consistent cycles performed
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  m_lda_to_qp(nbnds,nbnds,nibz,nsppol)= overwritten with the new QP amplitudes
!!                                        in terms of KS wavefunctions
!!
!! NOTES
!!  Only master node should call this routine.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine updt_m_lda_to_qp(Sigp,Kmesh,nscf,Sr,m_lda_to_qp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'updt_m_lda_to_qp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nscf
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Sigma_parameters),intent(in) :: Sigp
 type(Sigma_results),intent(in) :: Sr
!arrays
 complex(dpc),intent(inout) :: m_lda_to_qp(Sigp%nbnds,Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ik,is
!arrays
 complex(dpc),allocatable :: mtmp(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (nscf>=0) then ! Calculate the new m_lda_to_qp
   ABI_ALLOCATE(mtmp,(Sigp%nbnds,Sigp%nbnds))
   do is=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       mtmp(:,:)=m_lda_to_qp(:,:,ik,is)
       m_lda_to_qp(:,:,ik,is)=MATMUL(mtmp(:,:),Sr%eigvec_qp(:,:,ik,is))
     end do
   end do
 end if

 DBG_EXIT("COLL")

end subroutine updt_m_lda_to_qp

!----------------------------------------------------------------------

END MODULE m_qparticles
!!***
