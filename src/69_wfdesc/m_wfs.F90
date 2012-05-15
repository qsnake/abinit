!{\src2tex{textfont=tt}}
!****m* ABINIT/m_wfs
!! NAME
!!  m_wfs
!!
!! FUNCTION
!!  This module contains the declaration of the wfs_descriptor object and its methods.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG, FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_wfs

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_iterators

 use m_fstrings,       only : toupper, starts_with
 use m_io_tools,       only : get_unit, flush_unit
 use m_numeric_tools,  only : imin_loc
 use m_blas,           only : xcopy, xdotc
 use m_fft_mesh,       only : print_ngfft, rotate_fft_mesh, ceigr, ceikr, check_rot_fft, fft_check_rotrans
 use m_fftw3,          only : fftw3_fftpad_cplx, fftw3_fftpad
 use m_crystal,        only : crystal_structure
 use m_gsphere,        only : get_kg
 use m_header,         only : hdr_init_lowlvl, hdr_init, hdr_clean
 use m_bz_mesh,        only : bz_mesh_type, get_bz_item
 use m_paw_toolbox,    only : pawfgrtab_init, pawfgrtab_free, pawfgrtab_print,&
&                             init_paw_pwaves_lmn, destroy_paw_pwaves_lmn, paw_pwaves_lmn_t
 use m_ebands,         only : update_occ, bstruct_init, bstruct_clean
 use m_wffile,         only : wffile_type

 implicit none

 private
!!***

 ! Flags giving the status of the local %ug, %ur %cprj buffers.
 integer,public,parameter :: WFD_NOWAVE   =0
 integer,public,parameter :: WFD_ALLOCATED=1
 integer,public,parameter :: WFD_STORED   =2

 integer,public,parameter :: CPR_RANDOM   =1
 integer,public,parameter :: CPR_SORTED   =2

 ! ID used to identify different instances of wfs_descriptor
 integer,private,save :: WFD_ID=1

!----------------------------------------------------------------------

!!****t* m_wfs/kdata_t
!! NAME
!! kdata_t
!!
!! FUNCTION
!! Datatype storing k-dependent quantities and tables needed for performing the zero-padded FFT of wavefunctions.
!!
!! SOURCE

 type,public :: kdata_t

   integer :: istwfk
   ! Storage mode for this k point.

   integer :: npw
   ! Number of plane-waves for this k-point.

   integer :: useylm
   ! 1 if nonlocal part is applied using real spherical Harmonics. 0 for Legendre polynomial.

   integer :: has_ylm
   ! 0 if ylm is allocated.
   ! 1 if ylm is allocated.
   ! 2 if ylm is already computed.

   integer :: gbnds(3,2)
   ! gbnds(:,1)=Minval of kg_k.
   ! gbnds(:,2)=Maxval of kg_k.
   ! TODO: TO BE REMOVED when k-centered basis sets will be used.

   integer,pointer :: gc2kg(:,:,:)  SET2NULL
   ! gc2kg(ng1,ng2,ng3)
   ! The index of (g1,g2,g3) in the array kg_k
   ! Use the following indexing (N means ngfft of the adequate direction)
   ! 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= kg
   ! 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index
   ! The table is used to symmetrize the periodic part of the Bloch state in reciprocal space.

   integer,pointer :: ig1_inver(:,:)   SET2NULL
   ! ig1_inver(ng1,8)

   integer,pointer :: ig2_inver(:,:)   SET2NULL
   ! ig2_inver(ng2,8)

   integer,pointer :: ig3_inver(:,:)   SET2NULL
   ! ig3_inver(ng3,8)

   integer, pointer :: kg_k(:,:)  SET2NULL
   ! kg_k(3,npw)
   ! G vector coordinates in reduced cordinates.

   integer,pointer :: igfft0(:)   SET2NULL
   ! igfft0(npw)
   ! Index of the G-sphere in the FFT box.

   integer,pointer :: gbound(:,:)     SET2NULL
   ! gbound(2*mgfft+8,2))
   ! The boundary of the basis sphere of G vectors at a given k point.
   ! for use in improved zero padding of ffts in 3 dimensions.

   !% real(dp) :: kpoint(3)

   real(dp),pointer :: ph3d(:,:,:)   SET2NULL
   ! ph3d(2,npw,natom)
   ! 3-dim structure factors, for each atom and each plane wave.

   real(dp),pointer :: phkxred(:,:)     SET2NULL
   ! phkxred(2,natom))
   ! e^{ik.Ra} for each atom. Packed according to the atom type (atindx).

   real(dp),pointer :: fnl_dir0der0(:,:,:,:)    SET2NULL
   ! fnl_dir0der0(npw,1,lmnmax,ntypat)
   ! nonlocal form factors.
   ! fnl(k+G).ylm(k+G) if PAW
   ! f_ln(k+G)/|k+G|^l if NC

   real(dp),pointer :: ylm(:,:)        SET2NULL
   ! ylm(npw,mpsang**2*useylm)
   ! Real spherical harmonics for each k+G

   !% real(dp),pointer :: kinpw(:)
   ! kinpw(npw_k)
   ! compute elements of kinetic energy operator in reciprocal space.
   ! (1/2*effmass) (2 Pi)**2 (k+G)**2:

 end type kdata_t

 public :: kdata_init
 public :: kdata_nullify
 public :: kdata_free
 public :: kdata_copy

 interface kdata_nullify
   module procedure nullify_kdata_0D
   module procedure nullify_kdata_1D
 end interface kdata_nullify

 interface kdata_free
   module procedure destroy_kdata_0D
   module procedure destroy_kdata_1D
 end interface kdata_free

 interface kdata_copy
   module procedure copy_kdata_0D
   module procedure copy_kdata_1D
 end interface kdata_copy
!!***

!----------------------------------------------------------------------

!!****t* m_wfs/wave_t
!! NAME
!! wave_t
!!
!! FUNCTION
!!  Structure used to store a single wavefunction in reciprocal space and, optionally, its real space representation.
!!
!! SOURCE

 type,public :: wave_t

  !integer :: npw_k
  !integer :: nfft
  !integer :: nspinor
  !integer :: natom

  !! integer :: cplex
  ! 1 for real wavefunctions u(r)
  ! 2 for complex wavefunctions u(r).
  ! At gamma we always have real u(r) provided that time-reversal can be used.
  ! In systems with both time-reversal and spatial inversion, wavefunctions can be chosen to be real.
  ! One might use this to reduce memory in wave_t.

  integer :: has_ug=WFD_NOWAVE
  ! Flag giving the status of ug.

  integer :: has_ur=WFD_NOWAVE
  ! Flag giving the status of ur.

  integer :: has_cprj=WFD_NOWAVE
  ! Flag giving the status of cprj.

  integer :: cprj_order=CPR_RANDOM
  ! Flag defining whether cprj are sorted by atom type or ordered according
  ! to the typat variable used in the input file.

  complex(gwpc),pointer :: ug(:)  SET2NULL
  ! ug(npw_k*nspinor)
  ! The periodic part of the Bloch wavefunction in reciprocal space.

  complex(gwpc),pointer :: ur(:)  SET2NULL
  ! ur(nfft*nspinor)
  ! The periodic part of the Bloch wavefunction in real space.

  type(cprj_type),pointer :: Cprj(:,:)  SET2NULL
  ! Cprj(natom,nspinor)
  ! PAW projected wave function <Proj_i|Cnk> with all NL projectors.

 end type wave_t
!!***

 public :: wave_nullify
 public :: wave_init
 public :: wave_free
 public :: wave_copy

 interface wave_nullify
   module procedure nullify_wave_0D
   module procedure nullify_wave_3D
 end interface wave_nullify

 interface wave_init
   module procedure init_wave_0D
 end interface wave_init

 interface wave_free
   module procedure destroy_wave_0D
   module procedure destroy_wave_3D
 end interface wave_free
                                     
 interface wave_copy
   module procedure copy_wave_0D
   module procedure copy_wave_3D
 end interface wave_copy
!!***

!----------------------------------------------------------------------

!!****t* m_wfs/iskg_tabs_t
!! NAME
!! iskg_tabs_t
!!
!! FUNCTION
!!  This structure contains data and tables used to obtain u_{ISk}(G) from u_k(G)
!!  where S is one of the symrec operations and I is either the identity or the inversion
!!  (time reversal in reciprocal space).
!!
!! SOURCE

 type,public :: iskg_tabs_t

   integer :: istwf_k
   ! Storage mode for the rotated k-point.

   integer :: npw_k
   ! Number of planewaves

   integer :: isym
   ! Index of the operation S in reciprocal space that is used to get the rotated k-point.

   integer :: itim
   ! 2 is time-reversal is used. 1 otherwise.

   !% integer :: ngfft(18)
   !% integer :: mgfft

   real(dp) :: ecut
   ! Cutoff energy.

   real(dp) :: kpt(3)
   ! Rotated k-point in reduced coordinates. kpt = IS(k_ibz)

   integer,pointer :: gbound(:,:)
   ! gbound(2*mgfft+8,2))
   ! The boundary of the G-sphere used in improved zero padding of ffts in 3 dimensions.

   integer,pointer :: igfft(:)
   ! igfft(npw_k)
   ! Index of the rotated G-sphere in the FFT box.

   integer,pointer :: gt_k2isk(:)
   ! gt_k2isk(npw_k)
   ! Table with the correspondence G --> (IS)G
   ! Note that G belongs the sphere centered on k_ibz where we have u_k(G)
   ! whereas (IS)G belongs to the sphere centered on the rotated k-point.

   integer,pointer :: kg_k(:,:)
   ! kg_k(3,npw_k)
   ! G of the rotated sphere in reduced cordinates.

   complex(dpc),pointer :: ph_mskpgt(:)
   ! ph_mskpgt(npw_k)
   ! Non-symmorphic phase factor e^{-i S(k_ibz+G).t}.
   ! NOTE that time-reversal symmetry is not used here to allow for a more
   ! compact expression when time-reversal has to be applied.

 end type iskg_tabs_t
!!***

 public :: init_iskg_tabs
 public :: nullify_iskg_tabs
 public :: destroy_iskg_tabs

!----------------------------------------------------------------------

!!****t* m_wfs/wfs_descriptor
!! NAME
!! wfs_descriptor
!!
!! FUNCTION
!! Container gathering information on the set of wavefunctions treated by
!! this node as well as their distribution inside the MPI communicator.
!!
!! SOURCE

 type,public :: wfs_descriptor

  integer :: id                 ! Identifier.
  integer :: debug_level=0      ! Internal flag defining the debug level.
  integer :: lmnmax
  integer :: mband              ! MAX(nband)
  integer :: mgfft              ! Maximum size of 1D FFTs i.e. MAXVAL(ngfft(1:3)), used to dimension some arrays.
  !% integer :: mpsang
  integer :: natom
  integer :: nfft               ! Number of FFT points treated by this processor
  integer :: nfftot             ! Total number of points in the FFT grid
  integer :: nkibz              ! Number of irreducible k-points
  integer :: npwwfn             ! Number of G vectors for wavefunctions
  integer :: nspden             ! Number of independent spin-density components
  integer :: nspinor            ! Number of spinorial components
  integer :: nsppol             ! Number of independent spin polarizations
  integer :: ntypat
  integer :: paral_kgb          ! Option for kgb parallelism
  integer :: usepaw             ! 1 if PAW is used, 0 otherwise.
  !% integer :: usepawu           ! 1 if PAW+U is used, 0 otherwise.
  integer :: prtvol             ! Verbosity level.
  integer :: pawprtvol          ! Verbosity level for PAW.
  integer :: usewvl             ! 1 if BigDFT is used, 0 otherwise.
  !integer :: useylm            ! 1 if nonlocal part is applied using Ylm instead of Pl.

  integer :: comm               ! The MPI communicator.
  integer :: master             ! The rank of master node in comm.
  integer :: my_rank            ! The rank of my processor inside the MPI communicator comm.
  integer :: nproc              ! The number of processors in MPI comm.


  logical :: rfft_is_symok      ! .TRUE. if the real space FFT mesh is compatible with the rotational
                                ! part of the space group.

  real(dp) :: dilatmx

  real(dp) :: ecut
   ! Cutoff for plane wave basis set.

  real(dp) :: ecutsm
   ! ecutsm=smearing energy for plane wave kinetic energy (Ha)
   ! Cutoff for plane wave basis set.

  !% real(dp) :: pawecutdg=zero
   ! Cutoff for plane wave basis set.

  logical :: gamma_centered=.TRUE.
  !logical :: gamma_centered=.FALSE.
   ! .TRUE. if ug are given on the Gamma-centered G-sphere. Flag nedded to preserve the old Implementation.

  !% real(dp) :: effmass
  ! Effective mass for electrons

!arrays
  integer :: ngfft(18)
   ! Information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft

  integer :: nloalg(5)
   ! Governs the choice of the algorithm for nonlocal operator. See doc.

  integer,pointer :: gvec(:,:)  SET2NULL
  ! gvec(3,npwwfn)
  ! Reduced coordinates of the planewaves.
  ! TODO This is redundant and should be removed when k-centered G-sphere will be used.

  integer,pointer :: irottb(:,:)  SET2NULL
   ! irottb(nfftot,nsym)
   ! Index of $R^{-1}(r-\tau)$ in the FFT box.

  integer,pointer :: istwfk(:)   SET2NULL
   ! istwfk(nkibz)
   ! Storage mode for this k-point. At present only istwfk==1 is supported.

  integer,pointer :: nband(:,:)    SET2NULL
   ! nband(nkibz,nsppol)
   ! Number of bands at each k-point and spin.

  integer,pointer :: indlmn(:,:,:)  SET2NULL
   ! indlmn(6,lmnmax,ntypat)
   ! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
   !                                or i=lmn (if useylm=1)

  integer,pointer :: nlmn_atm(:)  SET2NULL
   ! nlmn_atm(natom)
   ! Number of (n,l,m) channels for each atom. Only for PAW

  integer,pointer :: nlmn_sort(:)  SET2NULL
   ! nlmn_sort(natom)
   ! Number of (n,l,m) channels for each atom (sorted by atom type). Only for PAW

  integer,pointer :: nlmn_type(:)  SET2NULL
   ! nlmn_type(ntypat)
   ! Number of (n,l,m) channels for each type of atom. Only for PAW.

  integer,pointer  :: npwarr(:)   SET2NULL
   ! npwarr(nkibz)
   ! Number of plane waves for this k-point.

  real(dp),pointer :: kibz(:,:)   SET2NULL
   ! kibz(3,nkibz)
   ! Reduced coordinates of the k-points in the IBZ.

  integer,pointer :: bks_tab(:,:,:,:)   SET2NULL
   ! bks_tab(mband,nkibz,nsppol,0:nproc-1)
   ! Global table used to keep trace of the distribution of the (b,k,s) states on each node inside Wfd%comm.
   ! 1 if the node has this state. 0 otherwise.
   ! A node owns a wavefunction if the corresponding ug is allocated AND computed.
   ! If a node owns ur but not ug, or ug is just allocated then its entry in the table is zero.

  integer,pointer :: bks_comm(:,:,:)  SET2NULL
   ! spin_comm(0:mband,0:nkibz,0:nsppol)
   ! MPI communicators.
   ! bks_comm(0,0,spin) MPI communicator for spin
   ! bks_comm(0,ik_ibz,spin)  MPI communicator for k-points.

  real(dp),pointer :: ph1d(:,:)    SET2NULL
   ! ph1d(2,3*(2*mgfft+1)*natom)
   ! 1-dim structure factor phase information.

  logical,pointer :: keep_ur(:,:,:) SET2NULL
   ! keep(mband,nkibz,nsppol)
   ! Storage strategy: keep or not keep calculated u(r) in memory.

  type(kdata_t),pointer :: Kdata(:)  SET2NULL
   ! Kdata(nkibz)
   ! datatype storing k-dependent quantities.

  type(wave_t),pointer :: Wave(:,:,:)  SET2NULL
   ! Wave(mband,nkibz,nsppol)
   ! Array of structures storing the periodic part of the wavefunctions in reciprocal- and real-space.

  type(MPI_type) :: MPI_enreg
   ! The MPI_type structured datatype gather different information about the MPI parallelisation :
   ! number of processors, the index of my processor, the different groups of processors, etc ...

 end type wfs_descriptor
!!***

 public :: wfd_init                ! Main creation method.
 public :: wfd_destroy             ! Destructor.
 public :: wfd_copy                ! Copy routine
 public :: wfd_reset_ur_cprj       ! Reinitialize memory storage of u(r) and <p_i|psi>
 public :: wfd_get_ur              ! Get one wavefunction in real space from its (b,k,s) indeces.
 public :: wfd_get_cprj            ! Get one PAW projection <Proj_i|Cnk> with all NL projectors from its (b,k,s) indeces.
 public :: wfd_change_ngfft        ! Reinitialize internal FFT tables.
 public :: wfd_nullify             ! Set all pointers to null()
 public :: wfd_print               ! Printout of basic info.
 public :: fft_onewfn              ! Helper function performing a single zero-padding FFT from G to R.
 public :: fft_ur                  ! Helper function performing a single zero-padding FFT from R to G.
 public :: wfd_mkall_ur            ! Calculate all ur owned by this node at once.
 public :: wfd_ug2cprj             ! Get PAW cprj from its (b,k,s) indeces.
 public :: wfd_ptr_ug              ! Return a pointer to ug from its (b,k,s) indeces. Use it carefully!
 public :: wfd_ptr_ur              ! Return a pointer to ur from its (b,k,s) indeces. Use it carefully!
 public :: wfd_wave_free           ! Free internal buffers used to store the wavefunctions.
 public :: wfd_push_ug
 public :: wfd_ihave_ug            ! .TRUE. if the node has this ug with the specified status.
 public :: wfd_ihave_ur            ! .TRUE. if the node has this ur with the specified status.
 public :: wfd_ihave_cprj          ! .TRUE. if the node has this cprj with the specified status.
 public :: wfd_mybands
 public :: wfd_distribute_bands
 public :: wfd_iterator_bks
 public :: wfd_bks_distrb
 public :: wfd_update_bkstab
 public :: wfd_set_mpicomm
 public :: wfd_rotate
 public :: wfd_sanity_check
 public :: wfd_gather_g2k
 public :: wfd_distribute_bbp
 public :: wfd_distribute_kb_kpbp
 public :: wfd_iam_master
 public :: wfd_test_ortho
 public :: wfd_barrier
 public :: wfd_sym_ur
 public :: wfd_paw_get_aeur
 public :: wfd_plot_ur
 public :: wfd_read_wfk
 public :: check_sym_ug

!----------------------------------------------------------------------

CONTAINS  !==============================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/nullify_kdata_0D
!! NAME
!!  nullify_kdata_0D
!!
!! FUNCTION
!!  Set all pointers to null.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine nullify_kdata_0D(Kdata)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_kdata_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kdata_t),intent(inout) :: Kdata

!************************************************************************

 !@kdata_t

!integer pointers
 nullify(Kdata%gc2kg )
 nullify(Kdata%ig1_inver)
 nullify(Kdata%ig2_inver)
 nullify(Kdata%ig3_inver)
 nullify(Kdata%kg_k  )
 nullify(Kdata%igfft0)
 nullify(Kdata%gbound)

!real pointers
 nullify(Kdata%ph3d   )
 nullify(Kdata%phkxred)
 nullify(Kdata%fnl_dir0der0)
 nullify(Kdata%ylm    )

end subroutine nullify_kdata_0D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/nullify_kdata_1D
!! NAME
!!  nullify_kdata_1D
!!
!! FUNCTION
!!  Set all pointers to null.
!!

subroutine nullify_kdata_1D(Kdata)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_kdata_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kdata_t),intent(inout) :: Kdata(:)

!Local variables ------------------------------
!scalars
 integer :: ik

!************************************************************************

 do ik=LBOUND(Kdata,DIM=1),UBOUND(Kdata,DIM=1)
   call nullify_kdata_0D(Kdata(ik))
 end do

end subroutine nullify_kdata_1D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/kdata_init
!! NAME
!!  kdata_init
!!
!! FUNCTION
!!  Main creation method for the kdata_t datatype.
!!
!! PARENTS
!!      debug_tools,m_shirley,m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine kdata_init(Kdata,Cryst,Psps,kpoint,istwfk,ngfft,MPI_enreg,ecut,kg_k)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kdata_init'
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk
 real(dp),optional,intent(in) :: ecut
 type(crystal_structure),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(kdata_t),intent(inout) :: Kdata
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,optional,target,intent(in) :: kg_k(:,:)
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: kpoint(3)

!Local variables ------------------------------
!scalars
 integer,parameter :: dum_unkg=0,dum_unylm=0,ider0=0,idir0=0
 integer :: mpw_,istat,npw_k,dimffnl,useylmgr,nkpg,iatom,ipw,ist
 integer :: mkmem_,nkpt_,optder,mgfft
 integer :: iatm,matblk,ng1,ng2,ng3,ig1,ig2,ig3
 real(dp) :: arg
 logical :: ltest
 !character(len=500) :: msg
!arrays
 integer :: nband_(1),npwarr_(1)
 real(dp),allocatable :: ylmgr_k(:,:,:),kpg_k(:,:)
 real(dp),allocatable :: ph1d(:,:)
 logical,allocatable :: kg_mask(:)

!************************************************************************

 !@kdata_t
 Kdata%istwfk  = istwfk
 Kdata%useylm  = Psps%useylm

 if (PRESENT(ecut)) then ! Calculate G-sphere from input ecut.
  ltest = (.not.associated(Kdata%kg_k))
  ABI_CHECK(ltest,"Kdata%kg_k is associated!")
  call get_kg(kpoint,istwfk,ecut,Cryst%gmet,npw_k,Kdata%kg_k)
 else if (PRESENT(kg_k)) then ! Use input g-vectors.
   npw_k = SIZE(kg_k,DIM=2)
   ABI_ALLOCATE(Kdata%kg_k,(3,npw_k))
   Kdata%kg_k = kg_k
 else
   MSG_ERROR("Either ecut or kg_k must be present")
 end if
 Kdata%npw = npw_k

 mgfft = MAXVAL(ngfft(1:3))
 !
 ! Finds the boundary of the basis sphere of G vectors (for this k point)
 ! for use in improved zero padding of ffts in 3 dimensions.
 ABI_ALLOCATE(Kdata%gbound,(2*mgfft+8,2))
 call sphereboundary(Kdata%gbound,istwfk,Kdata%kg_k,mgfft,npw_k)
 !
 ! Index of the G-sphere in the FFT box.
 ABI_ALLOCATE(Kdata%igfft0,(npw_k))

 ABI_ALLOCATE(kg_mask,(npw_k))
 call kgindex(Kdata%igfft0,Kdata%kg_k,kg_mask,MPI_enreg,ngfft,npw_k)

 ABI_CHECK(ALL(kg_mask),"FFT para not yet implemented")
 ABI_DEALLOCATE(kg_mask)

 ! Compute e^{ik.Ra} for each atom. Packed according to the atom type (atindx).
 ABI_ALLOCATE(Kdata%phkxred,(2,Cryst%natom))
 do iatom=1,Cryst%natom
   iatm=Cryst%atindx(iatom)
   arg=two_pi*(DOT_PRODUCT(kpoint,Cryst%xred(:,iatom)))
   Kdata%phkxred(1,iatm)=DCOS(arg)
   Kdata%phkxred(2,iatm)=DSIN(arg)
 end do
 !
 ! Calculate 1-dim structure factor phase information.
 mgfft = MAXVAL(ngfft(1:3))
 ABI_ALLOCATE(ph1d,(2,3*(2*mgfft+1)*Cryst%natom))
 call getph(Cryst%atindx,Cryst%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,Cryst%xred)

 matblk=Cryst%natom
 ABI_ALLOCATE(Kdata%ph3d,(2,npw_k,matblk))
 call ph1d3d(1,Cryst%natom,Kdata%kg_k,matblk,Cryst%natom,npw_k,ngfft(1),ngfft(2),ngfft(3),Kdata%phkxred,ph1d,Kdata%ph3d)
 ABI_DEALLOCATE(ph1d)
 !
 ! * Compute spherical harmonics if required.
 Kdata%has_ylm = 0
 ABI_ALLOCATE(Kdata%ylm,(npw_k,Psps%mpsang**2*Psps%useylm))
 useylmgr=0
 ABI_ALLOCATE(ylmgr_k,(npw_k,3,Psps%mpsang**2*useylmgr))

 if (Kdata%useylm==1) then
   mkmem_=1; mpw_=npw_k; nband_=0; nkpt_=1; npwarr_(1)=npw_k
   optder=0 ! only Ylm(K) are computed.

   call initylmg(Cryst%gprimd,Kdata%kg_k,kpoint,mkmem_,MPI_enreg,Psps%mpsang,mpw_,nband_,nkpt_,&
&    npwarr_,1,optder,Cryst%rprimd,dum_unkg,dum_unylm,Kdata%ylm,ylmgr_k)

   Kdata%has_ylm=2
 end if
 !
 ! * Compute (k+G) vectors.
 nkpg=0
 ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
 if (nkpg>0) call mkkpg(Kdata%kg_k,kpg_k,kpoint,nkpg,npw_k)
 !
 ! * Compute nonlocal form factors fnl_dir0der0 for all (k+G).
 dimffnl=1+3*ider0
 ABI_ALLOCATE(Kdata%fnl_dir0der0,(npw_k,dimffnl,Psps%lmnmax,Cryst%ntypat))

 call mkffnl(Psps%dimekb,dimffnl,Psps%ekb,Kdata%fnl_dir0der0,Psps%ffspl,&
&  Cryst%gmet,Cryst%gprimd,ider0,idir0,Psps%indlmn,Kdata%kg_k,kpg_k,kpoint,Psps%lmnmax,&
&  Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,Cryst%ntypat,&
&  Psps%pspso,Psps%qgrid_ff,Cryst%rmet,Psps%usepaw,Psps%useylm,Kdata%ylm,ylmgr_k)

 ABI_DEALLOCATE(kpg_k)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(ylmgr_k)
 !
 ! Setup of tables used to symmetrize u(g)
 ! TODO: Be careful here as FFT parallelism won't work.
 Kdata%gbnds(:,1) = MINVAL(Kdata%kg_k(:,:),DIM=2)
 Kdata%gbnds(:,2) = MAXVAL(Kdata%kg_k(:,:),DIM=2)

 ng1 = Kdata%gbnds(1,2) + ABS(Kdata%gbnds(1,1)) + 1
 ng2 = Kdata%gbnds(2,2) + ABS(Kdata%gbnds(2,1)) + 1
 ng3 = Kdata%gbnds(3,2) + ABS(Kdata%gbnds(3,1)) + 1
 !write(std_out,*)"ng1,ng2,ng3: ",ng1,ng2,ng3

 ! Build correspondence btw G coordinates and sequential index.
 ABI_ALLOCATE(Kdata%gc2kg,(ng1,ng2,ng3))
 Kdata%gc2kg=0
 do ipw=1,npw_k
   ig1=Kdata%kg_k(1,ipw); if (ig1<0) ig1=ig1+ng1; ig1=ig1+1
   ig2=Kdata%kg_k(2,ipw); if (ig2<0) ig2=ig2+ng2; ig2=ig2+1
   ig3=Kdata%kg_k(3,ipw); if (ig3<0) ig3=ig3+ng3; ig3=ig3+1
   Kdata%gc2kg(ig1,ig2,ig3)=ipw
 end do

 ABI_ALLOCATE(Kdata%ig1_inver,(ng1,8))
 ABI_ALLOCATE(Kdata%ig2_inver,(ng2,8))
 ABI_ALLOCATE(Kdata%ig3_inver,(ng3,8))

 do ist=1,8
   call make_istwfk_table(ist,ng1,ng2,ng3,Kdata%ig1_inver(:,ist),Kdata%ig2_inver(:,ist),Kdata%ig3_inver(:,ist))
 end do

end subroutine kdata_init
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/destroy_kdata_0D
!! NAME
!!  destroy_kdata_0D
!!
!! FUNCTION
!!  Deallocate memory
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine destroy_kdata_0D(Kdata)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_kdata_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kdata_t),intent(inout) :: Kdata

!************************************************************************

 !@kdata_t
 if (associated(Kdata%gc2kg    ))  then
   ABI_DEALLOCATE(Kdata%gc2kg)
 end if
 if (associated(Kdata%ig1_inver))  then
   ABI_DEALLOCATE(Kdata%ig1_inver)
 end if
 if (associated(Kdata%ig2_inver))  then
   ABI_DEALLOCATE(Kdata%ig2_inver)
 end if
 if (associated(Kdata%ig3_inver))  then
   ABI_DEALLOCATE(Kdata%ig3_inver)
 end if

 if (associated(Kdata%kg_k  ))  then
   ABI_DEALLOCATE(Kdata%kg_k)
 end if
 if (associated(Kdata%igfft0))  then
   ABI_DEALLOCATE(Kdata%igfft0)
 end if
 if (associated(Kdata%gbound))  then
   ABI_DEALLOCATE(Kdata%gbound)
 end if

 if (associated(Kdata%ph3d   ))  then
   ABI_DEALLOCATE(Kdata%ph3d)
 end if
 if (associated(Kdata%phkxred))  then
   ABI_DEALLOCATE(Kdata%phkxred)
 end if
 if (associated(Kdata%fnl_dir0der0 ))  then
   ABI_DEALLOCATE(Kdata%fnl_dir0der0)
 end if
 if (associated(Kdata%ylm    ))  then
   ABI_DEALLOCATE(Kdata%ylm)
 end if

end subroutine destroy_kdata_0D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/destroy_kdata_1D
!! NAME
!!  destroy_kdata_1D
!!
!! FUNCTION
!!   Deallocate memory.
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine destroy_kdata_1D(Kdata)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_kdata_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kdata_t),intent(inout) :: Kdata(:)

!Local variables ------------------------------
!scalars
 integer :: ik

!************************************************************************

 do ik=LBOUND(Kdata,DIM=1),UBOUND(Kdata,DIM=1)
   call destroy_kdata_0d(Kdata(ik))
 end do

end subroutine destroy_kdata_1D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/copy_kdata_0D
!! NAME
!!  copy_kdata_0D
!!
!! FUNCTION
!!  Deallocate memory
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine copy_kdata_0D(Kdata_in,Kdata_out)

 use defs_basis
 use m_copy

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_kdata_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kdata_t),intent(in) :: Kdata_in
 type(kdata_t),intent(inout) :: Kdata_out

!************************************************************************

 !@kdata_t

 call deep_copy(Kdata_in%istwfk,Kdata_out%istwfk)
 call deep_copy(Kdata_in%npw,Kdata_out%npw)
 call deep_copy(Kdata_in%useylm,Kdata_out%useylm)
 call deep_copy(Kdata_in%has_ylm,Kdata_out%has_ylm)
 call deep_copy(Kdata_in%gc2kg,Kdata_out%gc2kg)
 call deep_copy(Kdata_in%ig1_inver,Kdata_out%ig1_inver)
 call deep_copy(Kdata_in%ig2_inver,Kdata_out%ig2_inver)
 call deep_copy(Kdata_in%ig3_inver,Kdata_out%ig3_inver)
 call deep_copy(Kdata_in%kg_k,Kdata_out%kg_k)
 call deep_copy(Kdata_in%igfft0,Kdata_out%igfft0)
 call deep_copy(Kdata_in%gbound,Kdata_out%gbound)
 call deep_copy(Kdata_in%ph3d,Kdata_out%ph3d)
 call deep_copy(Kdata_in%phkxred,Kdata_out%phkxred)
 call deep_copy(Kdata_in%fnl_dir0der0,Kdata_out%fnl_dir0der0)
 call deep_copy(Kdata_in%ylm,Kdata_out%ylm)

 Kdata_out%gbnds=Kdata_in%gbnds
end subroutine copy_kdata_0D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/copy_kdata_1D
!! NAME
!!  copy_kdata_1D
!!
!! FUNCTION
!!   Deallocate memory.
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine copy_kdata_1D(Kdata_in,Kdata_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_kdata_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(kdata_t),intent(in) :: Kdata_in(:)
 type(kdata_t),intent(inout) :: Kdata_out(:)

!Local variables ------------------------------
!scalars
 integer :: ik

!************************************************************************

 if (size(Kdata_in,DIM=1) .ne. size(Kdata_out,DIM=1)) stop "Error in copy_kdata_1D: wrong sizes ! "

 do ik=LBOUND(Kdata_in,DIM=1),UBOUND(Kdata_in,DIM=1)
   call copy_kdata_0d(Kdata_in(ik),Kdata_out(ik))
 end do

end subroutine copy_kdata_1D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_init
!! NAME
!! wfd_init
!!
!! FUNCTION
!!  Initialize the object.
!!
!! INPUTS
!!  Cryst<crystal_structure>=Object defining the unit cell and its symmetries.
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=PAW tabulated starting data.
!!  Psps<Pseudopotential_type>=datatype storing data on the pseudopotentials.
!!  ngfft(18)=All needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkibz=Number of irreducible k-points.
!!  npwwfn=Number of plane waves for u(G).
!!  nsppol=Number of independent spin polarizations.
!!  nspden=Number of density components.
!!  nspinor=Number of spinorial components.
!!  ecutsm
!!  dilatmx
!!  mband
!!  nband(nkibz,nsppol)
!!  keep_ur(mband,nkibz,nsppol)=Option for memory storage of u(r).
!!  paral_kgb=Option for band-FFT parallelism (not yet available)
!!  gvec(3,npwwfn)=G-vectors in reduced coordinates.
!!  istwfk(nkibz)=Storage mode.
!!  kibz(3,nkibz)=Reduced coordinates of the k-points.
!!  nloalg(5)=Governs the choice of the algorithm for nonlocal operator. See doc.
!!  prtvol=Verbosity level.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Initialize the object with basic dimensions, allocate also memory for u(g) and u(r) according to keep_ur
!!    %ug in reciprocal space are always allocated.
!!    %ur in real space only if keep_ur.
!!
!! PARENTS
!!      bethe_salpeter,m_shirley,screening,sigma
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_init(Wfd,Cryst,Pawtab,Psps,keep_ur,paral_kgb,npwwfn,mband,nband,nkibz,nsppol,bks_mask,&
&  nspden,nspinor,ecutsm,dilatmx,istwfk,kibz,ngfft,gvec,nloalg,prtvol,pawprtvol,comm,opt_ecut)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_init'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: paral_kgb,mband,comm,prtvol,pawprtvol
 integer,intent(in) :: nkibz,npwwfn,nsppol,nspden,nspinor
 real(dp),optional,intent(in) :: opt_ecut
 real(dp),intent(in) :: ecutsm,dilatmx
 type(crystal_structure),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(out) :: Wfd
!array
 integer,intent(in) :: ngfft(18),istwfk(nkibz),nband(nkibz,nsppol)
 integer,intent(in) :: gvec(3,npwwfn),nloalg(5)
 real(dp),intent(in) :: kibz(3,nkibz)
 logical,intent(in) :: bks_mask(mband,nkibz,nsppol)
 logical,intent(in) :: keep_ur(mband,nkibz,nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: ider0=0,idir0=0,dum_unkg=0,dum_unylm=0,nfft0=0,mpw0=0,ikg0=0
 integer :: ig,ik_ibz,spin,band
 integer :: exchn2n3d,istwf_k,npw_k
 integer :: iatom,itypat,iat
 real(dp) :: ug_size,ur_size,cprj_size,gsq,g1,g2,g3
 logical :: iscompatibleFFT
 character(len=500) :: msg
 type(MPI_type) :: MPI_seq
!arrays
 integer :: dum_kg(3,0)
 real(dp) :: kpoint(3)

!************************************************************************

 DBG_ENTER("COLL")

 !@wfs_descriptor
 call wfd_nullify(Wfd)

 Wfd%gamma_centered=.TRUE.
 if (PRESENT(opt_ecut)) then
   if (opt_ecut > tol6) then
     Wfd%gamma_centered=.FALSE.
     MSG_WARNING("Using k-centered G-spheres.")
   end if
 end if

 Wfd%id=WFD_ID; WFD_ID=WFD_ID+1

 ! MPI info
 Wfd%comm    = comm
 Wfd%my_rank = xcomm_rank(Wfd%comm)
 Wfd%nproc   = xcomm_size(Wfd%comm)
 Wfd%master  = 0

 ABI_ALLOCATE(Wfd%bks_comm,(0:mband,0:nkibz,0:nsppol))
 Wfd%bks_comm  = xmpi_comm_null
 !
 ! Sequential MPI datatype to be passed to abinit routines.
 call initmpi_seq(Wfd%MPI_enreg)
 !
 ! === Basic dimensions ===
 Wfd%nkibz     = nkibz
 Wfd%nsppol    = nsppol
 Wfd%nspden    = nspden
 Wfd%nspinor   = nspinor
 Wfd%npwwfn    = npwwfn
 Wfd%paral_kgb = paral_kgb
 Wfd%nloalg    = nloalg

 Wfd%usepaw = Psps%usepaw
 Wfd%usewvl = 0 ! wavelets are not supported.
 Wfd%natom  = Cryst%natom
 Wfd%ntypat = Cryst%ntypat
 Wfd%lmnmax = Psps%lmnmax
 Wfd%prtvol = prtvol
 Wfd%pawprtvol = pawprtvol

 Wfd%ecutsm  = ecutsm
 Wfd%dilatmx = dilatmx

 ABI_ALLOCATE(Wfd%indlmn,(6,Wfd%lmnmax,Wfd%ntypat))
 Wfd%indlmn = Psps%indlmn

 if (Wfd%usepaw==1) then
   ABI_ALLOCATE(Wfd%nlmn_atm,(Cryst%natom))
   ABI_ALLOCATE(Wfd%nlmn_type,(Cryst%ntypat))
   do iatom=1,Cryst%natom
     Wfd%nlmn_atm(iatom)=Pawtab(Cryst%typat(iatom))%lmn_size
   end do

   do itypat=1,Cryst%ntypat
     Wfd%nlmn_type(itypat)=Pawtab(itypat)%lmn_size
   end do

   ABI_ALLOCATE(Wfd%nlmn_sort,(Cryst%natom))
   iat=0 ! nlmn dims sorted by atom type.
   do itypat=1,Cryst%ntypat
     Wfd%nlmn_sort(iat+1:iat+Cryst%nattyp(itypat))=Pawtab(itypat)%lmn_size
     iat=iat+Cryst%nattyp(itypat)
   end do
 end if

 ABI_ALLOCATE(Wfd%keep_ur,(mband,nkibz,nsppol))
 Wfd%keep_ur=keep_ur
 !
 ! Setup of the FFT mesh
 Wfd%ngfft  = ngfft
 Wfd%mgfft  = MAXVAL (Wfd%ngfft(1:3))
 Wfd%nfftot = PRODUCT(Wfd%ngfft(1:3))
 Wfd%nfft   = Wfd%nfftot ! At present no FFT parallelism.
 !
 ! Calculate ecut from input gvec.
 if (Wfd%gamma_centered) then
   Wfd%ecut=-one
   do ig=1,npwwfn
     g1=REAL(gvec(1,ig))
     g2=REAL(gvec(2,ig))
     g3=REAL(gvec(3,ig))
     gsq=      Cryst%gmet(1,1)*g1**2+Cryst%gmet(2,2)*g2**2+Cryst%gmet(3,3)*g3**2+ &
&         two*(Cryst%gmet(1,2)*g1*g2+Cryst%gmet(1,3)*g1*g3+Cryst%gmet(2,3)*g2*g3)
     Wfd%ecut=MAX(Wfd%ecut,gsq)
   end do
   Wfd%ecut=two*Wfd%ecut*pi**2
 else
   Wfd%ecut=opt_ecut
 end if
 !
 ! Precalculate the FFT index of $ R^{-1} (r-\tau) $ used to symmetrize u_Rk.
 ABI_ALLOCATE(Wfd%irottb,(Wfd%nfftot,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,Wfd%ngfft,Wfd%irottb,iscompatibleFFT)
 if (.not.iscompatibleFFT) then
   msg = "FFT mesh is not compatible with symmetries. Wavefunction symmetrization might be affected by large errors!"
   MSG_WARNING(msg)
 end if
 !
 ! Is the real space mesh compatible with the rotational part?
 Wfd%rfft_is_symok = check_rot_fft(Cryst%nsym,Cryst%symrel,Wfd%ngfft(1),Wfd%ngfft(2),Wfd%ngfft(3))
 !
 !
 ABI_ALLOCATE(Wfd%kibz,(3,Wfd%nkibz))
 Wfd%kibz=kibz
 ABI_ALLOCATE(Wfd%istwfk,(Wfd%nkibz))
 Wfd%istwfk=istwfk

 if (ANY(Wfd%istwfk/=1)) then
   if (.not.Wfd%gamma_centered) then
     MSG_ERROR("if (ANY(Wfd%istwfk/=1) then not Wfd%gamma_centered")
   end if
   MSG_WARNING("istwfk/=1 still under development!")
   write(std_out,*)Wfd%istwfk
 end if
 !
 ! * Get the number of planewaves npw_k
 ABI_ALLOCATE(Wfd%npwarr,(Wfd%nkibz))

 if (Wfd%gamma_centered) then
   Wfd%npwarr = npwwfn
 else
   ! TODO Here we should use ecut_eff instead of ecut
   exchn2n3d=0
   do ik_ibz=1,Wfd%nkibz
     istwf_k = Wfd%istwfk(ik_ibz)
     kpoint  = Wfd%kibz(:,ik_ibz)
     call kpgsph(Wfd%ecut,exchn2n3d,Cryst%gmet,ikg0,ik_ibz,istwf_k,dum_kg,kpoint,0,Wfd%MPI_enreg,mpw0,npw_k)
     Wfd%npwarr(ik_ibz)= npw_k
   end do
 end if

 ABI_ALLOCATE(Wfd%gvec,(3,npwwfn))
 Wfd%gvec=gvec  ! TODO For the time being, continue to use Gamma-centered basis set in Wfd%gvec.

 ABI_ALLOCATE(Wfd%nband,(nkibz,nsppol))
 Wfd%nband=nband

 Wfd%mband = mband
 ABI_CHECK(MAXVAL(Wfd%nband)==mband,"wrong mband")

 ! === Allocate u(g) and, if required, also u(r) ===
 ug_size=nspinor*npwwfn*COUNT(bks_mask)
 write(msg,'(a,f12.1,a)')' Memory needed for storing ug= ',2*gwpc*ug_size*b2Mb,' [Mb]'
 call wrtout(std_out,msg,'PERS')

 if (Wfd%usepaw==1) then
   cprj_size=nspinor*SUM(Wfd%nlmn_atm)*COUNT(bks_mask)
   write(msg,'(a,f12.1,a)')' Memory needed for storing Cprj= ',dp*cprj_size*b2Mb,' [Mb]'
   call wrtout(std_out,msg,'PERS')
 end if

 ur_size=nspinor*Wfd%nfft*COUNT(Wfd%keep_ur)
 write(msg,'(a,f12.1,a)')' Memory needed for storing ur= ',2*gwpc*ur_size*b2Mb,' [Mb]'
 call wrtout(std_out,msg,'PERS')

 ABI_ALLOCATE(Wfd%Wave,(Wfd%mband,Wfd%nkibz,Wfd%nsppol))
 call wave_nullify(Wfd%Wave)

 ! Allocate the wavefunctions in reciprocal space according to bks_mask.
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     npw_k = Wfd%npwarr(ik_ibz)
     do band=1,Wfd%nband(ik_ibz,spin)
       if (bks_mask(band,ik_ibz,spin)) then
         !if (Wfd%keep_ur(band,ik_ibz,spin)) then
         !  call init_wave_0D(Wfd%Wave(band,ik_ibz,spin),Wfd%usepaw,npw_k,Wfd%nfft,Wfd%nspinor,Wfd%natom,Wfd%nlmn_atm,CPR_RANDOM)
         !else
         call init_wave_0D(Wfd%Wave(band,ik_ibz,spin),Wfd%usepaw,npw_k,nfft0,Wfd%nspinor,Wfd%natom,Wfd%nlmn_atm,CPR_RANDOM)
         !end if
       end if
     end do
   end do
 end do

 ! Allocate the global table used to keep trace of the distribution, including a possible duplication.
 ABI_ALLOCATE(Wfd%bks_tab,(Wfd%mband,nkibz,nsppol,0:Wfd%nproc-1))
 Wfd%bks_tab=WFD_NOWAVE

 ! Update the kbs table storing the distribution of the ug.
 call wfd_update_bkstab(Wfd)
 !
 ! Initialize the MPI communicators.
 ! init MPI communicators:cannot be done here since waves are not stored yet.
 !call wfd_set_mpicomm(Wfd)
 !
 ! ===================================================
 ! ==== Precalculate nonlocal form factors for PAW ====
 ! ===================================================
 !
 ! Calculate 1-dim structure factor phase information.
 ABI_ALLOCATE(Wfd%ph1d,(2,3*(2*Wfd%mgfft+1)*Wfd%natom))
 call getph(Cryst%atindx,Wfd%natom,Wfd%ngfft(1),Wfd%ngfft(2),Wfd%ngfft(3),Wfd%ph1d,Cryst%xred)

 call initmpi_seq(MPI_seq)

 ABI_ALLOCATE(Wfd%Kdata,(Wfd%nkibz))
 call kdata_nullify(Wfd%Kdata)    ! Nullify the pointers defined in the datatype.

 do ik_ibz=1,Wfd%nkibz
   kpoint  = Wfd%kibz(:,ik_ibz)
   istwf_k = Wfd%istwfk(ik_ibz)
   npw_k   = Wfd%npwarr(ik_ibz)
   if (wfd_ihave_ug(Wfd,0,ik_ibz,0)) then
     if (Wfd%gamma_centered) then
       call kdata_init(Wfd%Kdata(ik_ibz),Cryst,Psps,kpoint,istwf_k,ngfft,Wfd%MPI_enreg,kg_k=Wfd%gvec)
     else
       call kdata_init(Wfd%Kdata(ik_ibz),Cryst,Psps,kpoint,istwf_k,ngfft,Wfd%MPI_enreg,ecut=Wfd%ecut)
     end if
   end if
 end do

 DBG_EXIT("COLL")

end subroutine wfd_init
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_destroy
!! NAME
!!  wfd_destroy
!!
!! FUNCTION
!!  Free the memory allocated in the wfs_descriptor data type.
!!
!! PARENTS
!!      bethe_salpeter,exc_interp_ham,m_shexc,screening,sigma
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_destroy(Wfd)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_destroy'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(inout) :: Wfd
!************************************************************************

 DBG_ENTER("COLL")

 !@wfs_descriptor
 !
 ! integer pointers.
 if (associated(Wfd%gvec     ))  then
   ABI_DEALLOCATE(Wfd%gvec)
 end if
 if (associated(Wfd%irottb   ))  then
   ABI_DEALLOCATE(Wfd%irottb)
 end if
 if (associated(Wfd%istwfk   ))  then
   ABI_DEALLOCATE(Wfd%istwfk)
 end if
 if (associated(Wfd%nband    ))  then
   ABI_DEALLOCATE(Wfd%nband)
 end if
 if (associated(Wfd%indlmn   ))  then
   ABI_DEALLOCATE(Wfd%indlmn)
 end if
 if (associated(Wfd%nlmn_atm ))  then
   ABI_DEALLOCATE(Wfd%nlmn_atm)
 end if
 if (associated(Wfd%nlmn_sort))  then
   ABI_DEALLOCATE(Wfd%nlmn_sort)
 end if
 if (associated(Wfd%nlmn_type))  then
   ABI_DEALLOCATE(Wfd%nlmn_type)
 end if
 if (associated(Wfd%npwarr   ))  then
   ABI_DEALLOCATE(Wfd%npwarr)
 end if
 if (associated(Wfd%bks_tab  ))  then
   ABI_DEALLOCATE(Wfd%bks_tab)
 end if

 ! Free the MPI communicators.
 call xcomm_free(Wfd%bks_comm)
 if (associated(Wfd%bks_comm ))  then
   ABI_DEALLOCATE(Wfd%bks_comm)
 end if

 ! real pointers.
 if (associated(Wfd%kibz))  then
   ABI_DEALLOCATE(Wfd%kibz)
 end if
 if (associated(Wfd%ph1d))  then
   ABI_DEALLOCATE(Wfd%ph1d)
 end if
 !
 ! logical pointers.
 if (associated(Wfd%keep_ur))  then
   ABI_DEALLOCATE(Wfd%keep_ur)
 end if
 !
 ! datatypes.
 if (associated(Wfd%Kdata)) then
   call kdata_free(Wfd%Kdata)
   ABI_DEALLOCATE(Wfd%Kdata)
 end if

 if (associated(Wfd%Wave)) then
   call wave_free(Wfd%Wave)
   ABI_DEALLOCATE(Wfd%Wave)
 end if

 call destroy_mpi_enreg(Wfd%MPI_enreg)

 DBG_EXIT("COLL")

end subroutine wfd_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_copy
!! NAME
!!  wfd_copy
!!
!! FUNCTION
!!  Duplicates a wfs_descriptor data type.
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_copy(Wfd_in,Wfd_out)

 use defs_basis
 use m_copy

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_copy'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wfs_descriptor),intent(inout) :: Wfd_in,Wfd_out

!************************************************************************

 DBG_ENTER("COLL")

 !@wfs_descriptor

 call deep_copy(Wfd_in%id             ,Wfd_out%id)
 call deep_copy(Wfd_in%debug_level    ,Wfd_out%debug_level)
 call deep_copy(Wfd_in%lmnmax         ,Wfd_out%lmnmax)
 call deep_copy(Wfd_in%mband          ,Wfd_out%mband)
 call deep_copy(Wfd_in%mgfft          ,Wfd_out%mgfft)
 call deep_copy(Wfd_in%natom          ,Wfd_out%natom)
 call deep_copy(Wfd_in%nfft           ,Wfd_out%nfft)
 call deep_copy(Wfd_in%nfftot         ,Wfd_out%nfftot)
 call deep_copy(Wfd_in%nkibz          ,Wfd_out%nkibz)
 call deep_copy(Wfd_in%npwwfn         ,Wfd_out%npwwfn)
 call deep_copy(Wfd_in%nspden         ,Wfd_out%nspden)
 call deep_copy(Wfd_in%nspinor        ,Wfd_out%nspinor)
 call deep_copy(Wfd_in%nsppol         ,Wfd_out%nsppol)
 call deep_copy(Wfd_in%ntypat         ,Wfd_out%ntypat)
 call deep_copy(Wfd_in%paral_kgb      ,Wfd_out%paral_kgb)
 call deep_copy(Wfd_in%usepaw         ,Wfd_out%usepaw)
 call deep_copy(Wfd_in%prtvol         ,Wfd_out%prtvol)
 call deep_copy(Wfd_in%pawprtvol      ,Wfd_out%pawprtvol)
 call deep_copy(Wfd_in%usewvl         ,Wfd_out%usewvl)
 call deep_copy(Wfd_in%comm           ,Wfd_out%comm)
 call deep_copy(Wfd_in%master         ,Wfd_out%master)
 call deep_copy(Wfd_in%my_rank        ,Wfd_out%my_rank)
 call deep_copy(Wfd_in%nproc          ,Wfd_out%nproc)
 call deep_copy(Wfd_in%rfft_is_symok  ,Wfd_out%rfft_is_symok)
 call deep_copy(Wfd_in%dilatmx        ,Wfd_out%dilatmx)
 call deep_copy(Wfd_in%ecut           ,Wfd_out%ecut)
 call deep_copy(Wfd_in%ecutsm         ,Wfd_out%ecutsm)
 call deep_copy(Wfd_in%gamma_centered ,Wfd_out%gamma_centered)
               Wfd_out%ngfft          =Wfd_in%ngfft
               Wfd_out%nloalg         =Wfd_in%nloalg
 call deep_copy(Wfd_in%gvec           ,Wfd_out%gvec)
 call deep_copy(Wfd_in%irottb         ,Wfd_out%irottb)
 call deep_copy(Wfd_in%istwfk         ,Wfd_out%istwfk)
 call deep_copy(Wfd_in%nband          ,Wfd_out%nband)
 call deep_copy(Wfd_in%indlmn         ,Wfd_out%indlmn)
 call deep_copy(Wfd_in%nlmn_atm       ,Wfd_out%nlmn_atm)
 call deep_copy(Wfd_in%nlmn_sort      ,Wfd_out%nlmn_sort)
 call deep_copy(Wfd_in%nlmn_type      ,Wfd_out%nlmn_type)
 call deep_copy(Wfd_in%npwarr         ,Wfd_out%npwarr)
 call deep_copy(Wfd_in%kibz           ,Wfd_out%kibz)
 call deep_copy(Wfd_in%bks_tab        ,Wfd_out%bks_tab)
 call deep_copy(Wfd_in%bks_comm       ,Wfd_out%bks_comm)
 call deep_copy(Wfd_in%ph1d           ,Wfd_out%ph1d)
 call deep_copy(Wfd_in%keep_ur        ,Wfd_out%keep_ur)

 
 if (size(Wfd_in%Kdata,DIM=1) .ne. size(Wfd_out%Kdata,DIM=1)) then
  if (associated(Wfd_out%Kdata))  then
    ABI_DEALLOCATE(Wfd_out%Kdata)
  end if
  ABI_ALLOCATE(Wfd_out%Kdata,(Wfd_out%nkibz))
 end if

 call kdata_copy(Wfd_in%Kdata,Wfd_out%Kdata)

 if (size(Wfd_in%Wave) .ne. size(Wfd_out%Wave)) then
  if (associated(Wfd_out%Wave))  then
    ABI_DEALLOCATE(Wfd_out%Wave)
  end if
  ABI_ALLOCATE(Wfd_out%Wave,(Wfd_out%mband,Wfd_out%nkibz,Wfd_out%nsppol))
 end if
 call wave_copy(Wfd_in%Wave,Wfd_out%Wave)

 call copy_mpi_enreg(Wfd_in%MPI_enreg,Wfd_out%MPI_enreg,1)

end subroutine wfd_copy
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_reset_ur_cprj
!! NAME
!!  wfd_reset_ur_cprj
!!
!! FUNCTION
!!  Reinitialize the storage mode of the ur treated by this node.
!!
!! PARENTS
!!      bethe_salpeter,m_shirley,sigma
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_reset_ur_cprj(Wfd)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_reset_ur_cprj'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(inout) :: Wfd

!************************************************************************

 where (Wfd%Wave(:,:,:)%has_ur == WFD_STORED)
   Wfd%Wave(:,:,:)%has_ur = WFD_ALLOCATED
 end where

 if (Wfd%usepaw==1) then
   where (Wfd%Wave(:,:,:)%has_cprj == WFD_STORED)
     Wfd%Wave(:,:,:)%has_cprj = WFD_ALLOCATED
   end where
 end if

end subroutine wfd_reset_ur_cprj
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_get_ur
!! NAME
!!  wfd_get_ur
!!
!! FUNCTION
!!  Get a wave function in real space, either by doing a G-->R FFT
!!  or by just retrieving the data already stored in Wfd.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=the wavefunction descriptor.
!!  band=Band index.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!
!! OUTPUT
!!  ur(Wfd%nfft*Wfd%nspinor)=The wavefunction in real space.
!!
!! PARENTS
!!      calc_sig_ppm_eet,calc_sigc_me,calc_sigx_me,calc_vhxc_me,cchi0,cchi0q0
!!      cchi0q0_intraband,check_completeness,classify_bands,cohsex_me
!!      exc_build_block,exc_build_ham,exc_den,m_oscillators,m_shirley,m_wfs
!!      wfd_mkrho
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_get_ur(Wfd,band,ik_ibz,spin,ur)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_get_ur'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 complex(gwpc),intent(out) :: ur(Wfd%nfft*Wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp=5,npw0=0
 integer :: npw_k,nfft,nspinor,has_this_ur
 character(len=500) :: msg
!arrays
 integer,pointer :: kg_k(:,:),gbound(:,:),igfft0(:)
 complex(gwpc),pointer :: wave_ug(:)
!************************************************************************

 has_this_ur = Wfd%Wave(band,ik_ibz,spin)%has_ur

 SELECT CASE (has_this_ur)

 CASE (WFD_NOWAVE, WFD_ALLOCATED) ! FFT is required.
   npw_k  = Wfd%npwarr(ik_ibz)
   nfft   = Wfd%nfft
   nspinor= Wfd%nspinor

   if (.not.wfd_ihave_ug(Wfd,band,ik_ibz,spin,"Stored")) then
     write(msg,'(a,3(i0,1x),a)')" ug for (band, ik_ibz, spin): ",band,ik_ibz,spin," is not stored in memory!"
     MSG_PERS_ERROR(msg)
   end if

   wave_ug => Wfd%Wave(band,ik_ibz,spin)%ug
   kg_k    => Wfd%Kdata(ik_ibz)%kg_k
   igfft0  => Wfd%Kdata(ik_ibz)%igfft0

   gbound  => Wfd%Kdata(ik_ibz)%gbound(:,:)

   call fft_onewfn(Wfd%paral_kgb,Wfd%istwfk(ik_ibz),nspinor,npw_k,nfft,Wfd%mgfft,Wfd%ngfft,&
&    wave_ug,ur,igfft0,kg_k,gbound,tim_fourdp,Wfd%MPI_enreg)

   if (Wfd%keep_ur(band,ik_ibz,spin)) then ! Store results
     if (has_this_ur==WFD_NOWAVE) then ! Alloc buffer for ur.
       call init_wave_0D(Wfd%Wave(band,ik_ibz,spin),Wfd%usepaw,npw0,Wfd%nfft,Wfd%nspinor,Wfd%natom,Wfd%nlmn_atm,CPR_RANDOM)
     end if
     call xcopy(Wfd%nfft*Wfd%nspinor,ur,1,Wfd%Wave(band,ik_ibz,spin)%ur,1)
     Wfd%Wave(band,ik_ibz,spin)%has_ur=WFD_STORED
   end if

 CASE (WFD_STORED) ! copy it back.
   call xcopy(Wfd%nfft*Wfd%nspinor,Wfd%Wave(band,ik_ibz,spin)%ur,1,ur,1)

 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong has_ur: ",has_this_ur
   MSG_PERS_BUG(msg)
 END SELECT

end subroutine wfd_get_ur
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_nullify
!! NAME
!!  wfd_nullify
!!
!! FUNCTION
!!  Nullify the pointers of the data structure.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_nullify(Wfd)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_nullify'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(inout) :: Wfd
!************************************************************************

 !@wfs_descriptor

! integer
 nullify(Wfd%irottb)
 nullify(Wfd%istwfk)
 nullify(Wfd%nband )
 nullify(Wfd%indlmn)
 nullify(Wfd%nlmn_atm )
 nullify(Wfd%nlmn_sort)
 nullify(Wfd%nlmn_type)
 nullify(Wfd%npwarr)
 nullify(Wfd%gvec  )

!integer arrays.
 nullify(Wfd%bks_tab)
 nullify(Wfd%bks_comm )

!real arrays
 nullify(Wfd%kibz)
 nullify(Wfd%ph1d)

!logical arrays
 nullify(Wfd%keep_ur)

! pointers to datatypes.
 nullify(Wfd%Kdata)
 nullify(Wfd%Wave)

! datatypes
 call nullify_mpi_enreg(Wfd%MPI_enreg)

end subroutine wfd_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_print
!! NAME
!! wfd_print
!!
!! FUNCTION
!!  Print the content of a wfs_descriptor datatype
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=The datatype.
!!  [header]=String to be printed as header for additional info.
!!  [unit]=Unit number for output
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS". Defaults to "COLL".
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!      bethe_salpeter,m_shirley,screening,sigma
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_print(Wfd,header,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(wfs_descriptor),intent(in) :: Wfd

!Local variables-------------------------------
!scalars
 integer :: my_prtvol,my_unt
 character(len=4) :: my_mode
 character(len=500) :: msg
! *************************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the Wfd% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(3(a,i5,a),a,i5,2a,f5.1)')&
&  '  Number of irreducible k-points ........ ',Wfd%nkibz,ch10,&
&  '  Number of spinorial components ........ ',Wfd%nspinor,ch10,&
&  '  Number of spin-density components ..... ',Wfd%nspden,ch10,&
&  '  Number of spin polarizations .......... ',Wfd%nsppol,ch10,&
&  '  Plane wave cutoff energy .............. ',Wfd%ecut
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(4(a,i0,a))')&
&  '  Max number of G-vectors ............... ',MAXVAL(Wfd%npwarr),ch10,&
&  '  Total number of FFT points ............ ',Wfd%nfftot,ch10,&
&  '  Number of FFT points treated by me .... ',Wfd%nfft,ch10,&
&  '  Parallelism over k-b-g (paral_kgb) .... ',Wfd%paral_kgb,ch10
 call wrtout(my_unt,msg,my_mode)

 call print_ngfft(Wfd%ngfft,'FFT mesh for wavefunctions',my_unt,my_mode,my_prtvol)

 !TODO
 ! Add addition info

end subroutine wfd_print
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/fft_onewfn
!! NAME
!! fft_onewfn
!!
!! FUNCTION
!! Low-level routine used to calculate ONE wavefunction in real space via FFT.
!!
!! INPUTS
!! nspinor=number of spinorial components
!! istwfk=Option describing the storage of the wavefunction. (at present must be 1)
!! igfft(npw_k)=index of each plane wave in FFT grid
!! mgfft=Max number of FFT divisions
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nkibz=number of k points
!! npw_k=number of plane waves for this k-point.
!! nsppol=number of independent spin polarizations
!! tim_fourdp=4 if called from within screening ; =5 if called from within sigma
!! ug(npw_k*nspinor)=wavefunctions in reciprocal space treated by this processor.
!! MPI_enreg= datatype containing information on parallelism to be passed to fourdp
!! gbound(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!!
!! OUTPUT
!!  ur(nfft*nspinor)=wavefunctions in real space.
!!
!! PARENTS
!!      calc_sig_ppm_eet,m_wfs,wfd_mkrho
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine fft_onewfn(paral_kgb,istwf_k,nspinor,npw_k,nfft,mgfft,ngfft,ug,ur,igfft,kg_k,gbound,tim_fourdp,MPI_enreg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_onewfn'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: paral_kgb,npw_k,nfft,tim_fourdp,nspinor,istwf_k,mgfft
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: igfft(npw_k)
 integer,intent(in) :: ngfft(18),gbound(2*mgfft+8,2)
 integer,intent(in) :: kg_k(3,npw_k)
 complex(gwpc),intent(in) :: ug(npw_k*nspinor)
 complex(gwpc),intent(out) :: ur(nfft*nspinor)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex0=0,option0=0,ndat1=1,iflag1=1,tim_fourwf0=0
 integer :: nx,ny,nz,ldx,ldy,ldz,fftalg,fftalga,fftalgc
 integer :: ispinor,ig,rspad,gspad
 integer :: ix,iy,iz,ifft
 real(dp),parameter :: weight1=one,xnorm1=one
 character(len=500) :: msg
!arrays
 integer :: shiftg(3),symm(3,3)
 integer :: dum_kg_kout(0,0)
 real(dp) :: dum_denpot(0,0,0),dum_fofgout(0,0)
 real(dp),allocatable :: fofgin(:,:),fofr(:,:,:,:),ftarr(:,:,:,:)
 complex(dpc),allocatable :: fg(:)

! *************************************************************************

 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=MOD(fftalg,10)
 nx=ngfft(1); ny=ngfft(2); nz=ngfft(3)

 SELECT CASE (fftalga)

 CASE (1,2,4) ! Goedecker routines or vendor provided FFT libraries.

   ldx=ngfft(4); ldy=ngfft(5); ldz=ngfft(6) ! Here augmentation is supported.
   ABI_ALLOCATE(fofgin,(2,npw_k))
   ABI_ALLOCATE(fofr,(2,ldx,ldy,ldz))

   if (istwf_k<=2) then

     do ispinor=1,nspinor
       gspad=(ispinor-1)*npw_k
       rspad=(ispinor-1)*nfft

       do ig=1,npw_k ! Have to convert from CPLX to REAL.
         fofgin(1,ig) = DBLE (ug(ig+gspad))
         fofgin(2,ig) = AIMAG(ug(ig+gspad))
       end do

#if 1
       call sg_fftrisc(cplex0,dum_denpot,fofgin,dum_fofgout,fofr,gbound,gbound,istwf_k,&
&        kg_k,dum_kg_kout,mgfft,ngfft,npw_k,0,ldx,ldy,ldz,option0,weight1)

#else
       ! TODO use fourwf but make sure that istwf_k>2 is supported with fftalg 112
       call fourwf(cplex0,dum_denpot,fofgin,dum_fofgout,fofr,gbound,gbound,istwf_k,&
&         kg_k,kg_k,mgfft,MPI_enreg,ndat1,ngfft,npw_k,npw_k,ldx,ldy,ldz,option0,&
&         paral_kgb,tim_fourwf0,weight1,weight1)
#endif

       do iz=1,nz     ! Fill the output array on the ngfft(1:3) mesh.
         do iy=1,ny
           do ix=1,nx
             ifft = ix + (iy-1)*nx + (iz-1)*nx*ny + rspad
             ur(ifft) = CMPLX(fofr(1,ix,iy,iz), fofr(2,ix,iy,iz), kind=gwpc)
           end do
         end do
       end do

       !MSG_WARNING("calling sg_fftpad! Results for R--G might be wrong")
       ! This call corresponds to fftalgc == 1
       !!!!!call sg_fftpad(ngfft(8),mgfft,ldx,ldy,ldz,nx,ny,nz,arr,ftarr,one,gbound)

       !call padded_fourwf_cplx(fg,ngfft,nx,ny,nz,ldx,ldy,ldz,mgfft,+1,gbound)
       !ur(rspad+1:rspad+nfft)=fg
     end do ! ispinor

   else  ! sg_fftpad does not accept istwf_k>2

     do ig=1,npw_k ! Have to convert from CPLX to REAL.
       fofgin(1,ig) = DBLE (ug(ig))
       fofgin(2,ig) = AIMAG(ug(ig))
     end do
     !
     ! Reconstruct full G-sphere.
     call sphere(fofgin,ndat1,npw_k,fofr,nx,ny,nz,ldx,ldy,ldz,kg_k,istwf_k,iflag1,MPI_enreg,shiftg,symm,xnorm1)

     MSG_WARNING("calling sg_fftpad! Results for R--G might be wrong")
     ABI_ALLOCATE(ftarr,(2,ldx,ldy,ldz))
     call sg_fftpad(ngfft(8),mgfft,ldx,ldy,ldz,nx,ny,nz,fofr,ftarr,one,gbound)

     do iz=1,nz     ! Fill the output array on the ngfft(1:3) mesh.
       do iy=1,ny
         do ix=1,nx
           ifft = ix + (iy-1)*nx + (iz-1)*nx*ny
           ur(ifft) = CMPLX(ftarr(1,ix,iy,iz), ftarr(2,ix,iy,iz), kind=gwpc)
         end do
       end do
     end do
     ABI_DEALLOCATE(ftarr)
   end if

   ABI_DEALLOCATE(fofgin)
   ABI_DEALLOCATE(fofr)
   RETURN

 CASE (3) ! FFTW3.
   !
   ! * FFT to give wfn in real space.
   ldx=nx; ldy=ny; ldz=nz ! TODO No augmentation at present

   if (istwf_k==1) then
     ABI_ALLOCATE(fg,(nfft))
     do ispinor=1,nspinor
       gspad=(ispinor-1)*npw_k
       rspad=(ispinor-1)*nfft
       !
       fg=czero ! Fill the FFT array from the PW array
       do ig=1,npw_k
         fg(igfft(ig))=ug(ig+gspad)
       end do

       !call fftw3_c2c_ip(nx,ny,nz,ldx,ldy,ldz,1,+1,fg) ! Version without Padding
       call fftw3_fftpad_cplx(fg,nx,ny,nz,ldx,ldy,ldz,mgfft,+1,gbound)
       ur(rspad+1:rspad+nfft)=fg
     end do ! ispinor
     ABI_DEALLOCATE(fg)

   else  ! NB: nspinor is always 1 in this case.

     MSG_ERROR("istwf_/=1 not tested with FFTW3")
     ! TODO: 312 is not correct --> one should use 311 for the zero-padded version.
     ABI_ALLOCATE(fofgin,(2,npw_k*ndat1))
     do ig=1,npw_k ! Convert from CPLX to REAL.
       fofgin(1,ig) = DBLE (ug(ig))
       fofgin(2,ig) = AIMAG(ug(ig))
     end do
     !
     ! Insert fofgin into fofr completing missing G.
     ABI_ALLOCATE(fofr,(2,ldx,ldy,ldz*ndat1))
     call sphere(fofgin,ndat1,npw_k,fofr,nx,ny,nz,ldx,ldy,ldz,kg_k,istwf_k,iflag1,MPI_enreg,shiftg,symm,xnorm1)
     ABI_DEALLOCATE(fofgin)

     call fftw3_fftpad(fofr,nx,ny,nz,ldx,ldy,ldz,mgfft,+1,gbound)

     do iz=1,nz     ! Fill the output array on the ngfft(1:3) mesh.
       do iy=1,ny
         do ix=1,nx
           ifft = ix + (iy-1)*nx + (iz-1)*nx*ny
           ur(ifft) = CMPLX(fofr(1,ix,iy,iz), fofr(2,ix,iy,iz), kind=gwpc)
         end do
       end do
     end do
     ABI_DEALLOCATE(fofr)
   end if

   RETURN

 CASE DEFAULT
   write(msg,"(a,i0,a)")" fftalga= ",fftalga," not available for GW calculations."
   MSG_ERROR(msg)
 END SELECT

 RETURN

 ABI_UNUSED((/paral_kgb,tim_fourdp/)) ! FIXME

end subroutine fft_onewfn
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/fft_ur
!! NAME
!! fft_ur
!!
!! FUNCTION
!! Low-level routine used to calculate one wavefunction in G-space space via FFT starting from its
!! real space representation.
!!
!! INPUTS
!! nspinor=number of spinorial components
!! istwfk=Option describing the storage of the wavefunction. (at present must be 1)
!! mgfft=Max number of FFT divisions
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nkibz=number of k points
!! npw_k=number of plane waves for this k-point.
!! nsppol=number of independent spin polarizations
!! tim_fourdp=4 if called from within screening ; =5 if called from within sigma
!!  ur(nfft*nspinor)=wavefunctions in real space.
!! MPI_enreg= datatype containing information on parallelism to be passed to fourdp
!! gbound(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!!
!! OUTPUT
!!  ug(npw_k*nspinor)=wavefunctions in reciprocal space treated by this processor.
!!
!! TODO
!!  Should merge fft_onewfn with fft_ur
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine fft_ur(nspinor,npw_k,istwf_k,paral_kgb,nfft,mgfft,ngfft,ur,gbound,kg_k,ug,MPI_enreg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_ur'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: paral_kgb,npw_k,nfft,nspinor,istwf_k,mgfft
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: ngfft(18),gbound(2*mgfft+8,2)
 integer,intent(in) :: kg_k(3,npw_k)
 complex(gwpc),intent(out) :: ug(npw_k*nspinor)
 complex(gwpc),intent(in) :: ur(nfft*nspinor)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex0=0,option3=3,ndat1=1,iflag_m1=-1,tim_fourwf0=0,npwin0=0
 integer :: nx,ny,nz,ldx,ldy,ldz,fftalg,fftalga,fftalgc
 integer :: ispinor,ig,rspad,gspad
 integer :: ix,iy,iz,ifft
 real(dp),parameter :: weight1=one,xnorm1=one
 character(len=500) :: msg
!arrays
 integer :: shiftg(3),symm(3,3),my_ngfft(18),dum_kg_kin(0,0)
 real(dp) :: dum_denpot(0,0,0),dum_fofgin(0,0)
 real(dp),allocatable :: fofgout(:,:),fofr(:,:,:,:),cg(:,:),fg(:,:,:,:)

! *************************************************************************

 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=MOD(fftalg,10)
 nx=ngfft(1); ny=ngfft(2); nz=ngfft(3)

 SELECT CASE (fftalga)

 CASE (1,2,4) ! Goedecker routines or vendor provided FFT libraries.

   ldx=ngfft(4); ldy=ngfft(5); ldz=ngfft(6)
   ABI_ALLOCATE(fofgout,(2,npw_k))
   ABI_ALLOCATE(fofr,(2,ldx,ldy,ldz))

   do ispinor=1,nspinor
     gspad=(ispinor-1)*npw_k
     rspad=(ispinor-1)*nfft

     do iz=1,nz  ! Fill fofr from input ur array.
       do iy=1,ny
         do ix=1,nx
           ifft = ix + (iy-1)*nx + (iz-1)*nx*ny + rspad
           fofr(1,ix,iy,iz) = REAL( ur(ifft))
           fofr(2,ix,iy,iz) = AIMAG(ur(ifft))
         end do
       end do
     end do
     !
     ! option=3 --> real space to reciprocal space.
     ! NOTE that in this case, fftalg=1x1 MUST be used. This may be changed in the future.
     my_ngfft=ngfft; my_ngfft(7)=111

     call fourwf(cplex0,dum_denpot,dum_fofgin,fofgout,fofr,gbound,gbound,istwf_k,&
&       dum_kg_kin,kg_k,mgfft,MPI_enreg,ndat1,my_ngfft,npwin0,npw_k,ldx,ldy,ldz,option3,&
&       paral_kgb,tim_fourwf0,weight1,weight1)

     do ig=1,npw_k ! Have to convert from REAL to CMPLX
       ug(ig+gspad) = DCMPLX( fofgout(1,ig), fofgout(2,ig))
     end do
   end do ! ispinor

   ABI_DEALLOCATE(fofgout)
   ABI_DEALLOCATE(fofr)
   RETURN

 CASE (3) ! FFTW3.
   !
   ldx=nx; ldy=ny; ldz=nz ! TODO No augmentation at present
   ABI_ALLOCATE(fg,(2,ldx,ldy,ldz))
   ABI_ALLOCATE(cg,(2,npw_k*ndat1))

   do ispinor=1,nspinor
     gspad=(ispinor-1)*npw_k
     rspad=(ispinor-1)*nfft
     !
     do iz=1,nz  ! Fill the fg array from input ur.
       do iy=1,ny
         do ix=1,nx
           ifft = ix + (iy-1)*nx + (iz-1)*nx*ny + rspad
           fg(1,ix,iy,iz) = REAL (ur(ifft))
           fg(2,ix,iy,iz) = AIMAG(ur(ifft))
         end do
       end do
     end do
     !
     ! In-place zero-padding FFT R-->G.
     call fftw3_fftpad(fg,nx,ny,nz,ldx,ldy,ldz,mgfft,-1,gbound)
     !
     ! iflag=-1 ==> extract the on-sphere cg from fg that given on the FFT box.
     call sphere(cg,ndat1,npw_k,fg,nx,ny,nz,ldx,ldy,ldz,kg_k,istwf_k,iflag_m1,MPI_enreg,shiftg,symm,xnorm1)

     do ig=1,npw_k ! Have to convert from REAL to CMPLX
       ug(ig+gspad) = DCMPLX(cg(1,ig), cg(2,ig))
     end do
   end do ! ispinor

   ABI_DEALLOCATE(fg)
   ABI_DEALLOCATE(cg)
   RETURN

 CASE DEFAULT
   write(msg,"(a,i0,a)")" fftalga= ",fftalga," not available for GW calculations."
   MSG_ERROR(msg)
 END SELECT

end subroutine fft_ur
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_mkall_ur
!! NAME
!! wfd_mkall_ur
!!
!! FUNCTION
!!  FFT transform from G to R the entire set of wavefunctions stored in the wfs_descriptor.
!!  Only those waves whose status is WFD_ALLOCATED are calculated unless force is used.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=Structure containing the wave functions for the GW.
!!
!! OUTPUT
!!  ncalc=Number of FFTs performed.
!!  [force]=If .TRUE. then force FFT even for waves whose status is WFD_STORED.
!!
!! SIDE EFFECTS
!!  %ur
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_mkall_ur(Wfd,ncalc,force)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_mkall_ur'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ncalc
 logical,optional,intent(in) :: force
 type(wfs_descriptor),intent(inout) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: spin,ik_ibz,band,npw_k
 logical :: do_fft
 integer,pointer :: kg_k(:,:),gbound(:,:),igfft0(:)

!************************************************************************

! TODO FFTs should be done in bunches.
!
 ncalc=0 !; if (.not.Wfd%keep_ur) RETURN

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     npw_k  =  Wfd%npwarr(ik_ibz)
     kg_k   => Wfd%Kdata(ik_ibz)%kg_k
     gbound => Wfd%Kdata(ik_ibz)%gbound
     igfft0 => Wfd%Kdata(ik_ibz)%igfft0

     do band=1,Wfd%nband(ik_ibz,spin)

       if (.not.Wfd%keep_ur(band,ik_ibz,spin)) CYCLE

       do_fft = wfd_ihave_ur(Wfd,band,ik_ibz,spin,"Allocated")
       if (PRESENT(force)) do_fft = (do_fft .or. wfd_ihave_ur(Wfd,band,ik_ibz,spin,"Stored"))

       if (do_fft) then
         call fft_onewfn(Wfd%paral_kgb,Wfd%istwfk(ik_ibz),Wfd%nspinor,npw_k,Wfd%nfft,Wfd%mgfft,Wfd%ngfft,&
&          Wfd%Wave(band,ik_ibz,spin)%ug,Wfd%Wave(band,ik_ibz,spin)%ur,igfft0,kg_k,gbound,0,Wfd%MPI_enreg)

         ncalc = ncalc + 1
         Wfd%Wave(band,ik_ibz,spin)%has_ur=WFD_STORED  ! Update the status
       end if

     end do
   end do
 end do

end subroutine wfd_mkall_ur
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_ug2cprj
!! NAME
!! wfd_ug2cprj
!!
!! FUNCTION
!!  Calculate the projected wave function <Proj_i|Cnk> with all NL projectors for a single
!!  k-point, band and spin.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=Structure containing the wave functions for the GW.
!!  ik_ibz=Index of the required k-point
!!  spin=Required spin index.
!!  choice=chooses possible output:
!!    In addition to projected wave function:
!!    choice=1 => nothing else
!!          =2 => 1st gradients with respect to atomic position(s)
!!          =3 => 1st gradients with respect to strain(s)
!!          =23=> 1st gradients with respect to atm. pos. and strain(s)
!!          =4 => 2nd derivatives with respect to atomic pos.
!!          =24=> 1st and 2nd derivatives with respect to atomic pos.
!!          =5 => 1st gradients with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!  idir=direction of the derivative, i.e. dir. of - atom to be moved  in the case choice=2
!!                                                 - strain component  in the case choice=3
!!                                                 - k point direction in the case choice=5
!!       Compatible only with choice=2,3,5; if idir=0, all derivatives are computed
!!  natom
!!  Cryst
!!  [sorted]=Logical flags defining if the output Cprj has to be sorted by atom type or not.
!!    By default, Cprj matrix elements are unsorted.
!!
!! OUTPUT
!!  cwaveprj
!!
!! PARENTS
!!      classify_bands,m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_ug2cprj(Wfd,band,ik_ibz,spin,choice,idir,natom,Cryst,cwaveprj,sorted)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_ug2cprj'
 use interfaces_44_abitypes_defs
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,idir,natom,band,ik_ibz,spin
 logical,optional,intent(in) :: sorted
 type(wfs_descriptor),intent(inout) :: Wfd
 type(crystal_structure),intent(in) :: Cryst
!arrays
 type(cprj_type),intent(inout) :: cwaveprj(natom,Wfd%nspinor)

!Local variables-------------------------------
!scalars
 integer :: cpopt,istwf_k,npw_k,nkpg,istat
 integer :: dimekb1,dimekb2,matblk,ia,iatm,dimffnl,itypat,iatom,isp
 logical :: want_sorted
!arrays
 integer,pointer :: kg_k(:,:)
 integer,allocatable :: dimcprj_srt(:)
 real(dp) :: dummy_ekb(0,0)
 real(dp) :: kpoint(3)
 real(dp),allocatable :: kpg(:,:)
 real(dp),pointer :: phkxred(:,:)
 real(dp),allocatable :: cwavef(:,:)
 !real(dp),allocatable :: ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),pointer :: ph3d(:,:,:)    ! ph3d(2,npw_k,matblk)
 real(dp),pointer :: ffnl(:,:,:,:)  ! ffnl(npw_k,dimffnl,lmnmax,ntypat)
 type(cprj_type),allocatable :: Cprj_srt(:,:)

! *********************************************************************

 ! different form factors have to be calculated and stored in Kdata.
 ABI_CHECK(choice==1,"choice/=1 not coded")

 dimffnl = 1
 npw_k   = Wfd%npwarr(ik_ibz)
 istwf_k = Wfd%istwfk(ik_ibz)
 kpoint  = Wfd%kibz(:,ik_ibz)

 kg_k    => Wfd%Kdata(ik_ibz)%kg_k
 ph3d    => Wfd%Kdata(ik_ibz)%ph3d
 ffnl    => Wfd%Kdata(ik_ibz)%fnl_dir0der0
 phkxred => Wfd%Kdata(ik_ibz)%phkxred

! Compute (k+G) vectors
 nkpg=0
 !% if (choice==3.or.choice==2.or.choice==23) nkpg=3*Wfd%nloalg(5)
 !% if (choice==4.or.choice==24) nkpg=9*Wfd%nloalg(5)
 ABI_ALLOCATE(kpg,(npw_k,nkpg))
 if (nkpg>0) call mkkpg(kg_k,kpg,kpoint,nkpg,npw_k)

 matblk = Cryst%natom
 !
 ! Copy wavefunction in reciprocal space.
 ABI_ALLOCATE(cwavef,(2,npw_k*Wfd%nspinor))
 cwavef(1,:) = DBLE (Wfd%Wave(band,ik_ibz,spin)%ug)
 cwavef(2,:) = AIMAG(Wfd%Wave(band,ik_ibz,spin)%ug)

 cpopt   = 0 ! Nothing is already calculated.
 dimekb1 = 0
 dimekb2 = 0

 want_sorted=.FALSE.; if (PRESENT(sorted)) want_sorted=sorted

 if (want_sorted) then ! Output cprj are sorted.

   call getcprj(choice,cpopt,cwavef,cwaveprj,dimekb1,dimekb2,dimffnl,dummy_ekb,ffnl,&
&    idir,Wfd%indlmn,istwf_k,kg_k,kpg,kpoint,Wfd%lmnmax,matblk,Wfd%mgfft,Wfd%MPI_enreg,&
&    Cryst%natom,Cryst%nattyp,Wfd%ngfft,nkpg,Wfd%nloalg,npw_k,Wfd%nspinor,Cryst%ntypat,&
&    phkxred,Wfd%ph1d,ph3d,Cryst%ucvol,Wfd%usepaw,1)

 else  ! Output cprj are unsorted.

   ABI_ALLOCATE(dimcprj_srt,(Cryst%natom))
   ia=0
   do itypat=1,Cryst%ntypat
     dimcprj_srt(ia+1:ia+Cryst%nattyp(itypat))=Wfd%nlmn_type(itypat)
     ia=ia+Cryst%nattyp(itypat)
   end do

   ABI_ALLOCATE(Cprj_srt,(natom,Wfd%nspinor))
   call cprj_alloc(Cprj_srt,0,dimcprj_srt)
   ABI_DEALLOCATE(dimcprj_srt)
   !
   ! Calculate sorted cprj.
   call getcprj(choice,cpopt,cwavef,Cprj_srt,dimekb1,dimekb2,dimffnl,dummy_ekb,ffnl,&
&    idir,Wfd%indlmn,istwf_k,kg_k,kpg,kpoint,Wfd%lmnmax,matblk,Wfd%mgfft,Wfd%MPI_enreg,&
&    Cryst%natom,Cryst%nattyp,Wfd%ngfft,nkpg,Wfd%nloalg,npw_k,Wfd%nspinor,Cryst%ntypat,&
&    phkxred,Wfd%ph1d,ph3d,Cryst%ucvol,Wfd%usepaw,1)
   !
   ! Reorder cprj (sorted --> unsorted)
   do iatom=1,Cryst%natom
     iatm=Cryst%atindx(iatom)
     do isp=1,Wfd%nspinor
       cwaveprj(iatom,isp)%cp=Cprj_srt(iatm,isp)%cp
     end do
   end do

   call cprj_free(Cprj_srt)
   ABI_DEALLOCATE(Cprj_srt)
 end if

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(kpg)
 istat = ABI_ALLOC_STAT

end subroutine wfd_ug2cprj
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/nullify_wave_0D
!! NAME
!!  nullify_wave_0D
!!
!! FUNCTION
!!  Set pointers to null.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine nullify_wave_0D(Wave)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_wave_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wave_t),intent(inout) :: Wave

!************************************************************************

 !@wave_t
 nullify(Wave%ug)
 nullify(Wave%ur)

! types
 nullify(Wave%Cprj)

end subroutine nullify_wave_0D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/nullify_wave_3D
!! NAME
!!  nullify_wave_3D
!!
!! FUNCTION
!!  Set all pointers to null.
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine nullify_wave_3D(Wave)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_wave_3D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wave_t),intent(inout) :: Wave(:,:,:)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3

!************************************************************************

 do i3=LBOUND(Wave,DIM=3),UBOUND(Wave,DIM=3)
   do i2=LBOUND(Wave,DIM=2),UBOUND(Wave,DIM=2)
     do i1=LBOUND(Wave,DIM=1),UBOUND(Wave,DIM=1)
       call nullify_wave_0D(Wave(i1,i2,i3))
     end do
   end do
 end do

end subroutine nullify_wave_3D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/init_wave_0D
!! NAME
!!  init_wave_0D
!!
!! FUNCTION
!!   Main creation method for the wave_t data type
!!
!! INPUTS
!!  usepaw=1 if PAW is used.
!!  npw =Number of plane-waves for ug
!!  nfft=Number of FFT points for the real space wavefunction.
!!  nspinor=Number of spinor components.
!!  natom=Number of atoms in cprj matrix elements.
!!  nlmn_size(natom)=Number of (n,l,m) channel for each atom.  Ordering of atoms depends on cprj_order
!!  cprj_order=Flag defining the ordering of the atoms in the cprj matrix elements (CPR_RANDOM|CPR_SORTED).
!!    Use to know if we have to reorder the matrix elements when wfd_get_cprj is called.
!!
!! OUTPUT
!!  Wave<wave_t>=The structure fully initialized.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine init_wave_0D(Wave,usepaw,npw,nfft,nspinor,natom,nlmn_size,cprj_order)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_wave_0D'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nfft,nspinor,usepaw,natom,cprj_order
 type(wave_t),intent(inout) :: Wave
!arrays
 integer,intent(in) :: nlmn_size(:)

!Local variables ------------------------------
!scalars
 integer,parameter :: ncpgr0=0  ! For the time being, no derivatives

!************************************************************************

 !Wave%npw_k   = npw
 !Wave%nfft    = nfft
 !Wave%nspinor = nspinor
 !Wave%natom   = natom

 !@wave_t
 if (npw >0) then
   ABI_ALLOCATE(Wave%ug,(npw*nspinor))
   Wave%has_ug=WFD_ALLOCATED
   Wave%ug=czero
   if (usepaw==1) then
     ABI_ALLOCATE(Wave%Cprj,(natom,nspinor))
     call cprj_alloc(Wave%Cprj,ncpgr0,nlmn_size)
     Wave%has_cprj=WFD_ALLOCATED
     Wave%cprj_order=cprj_order
   end if
 end if

 if (nfft>0) then
   ABI_ALLOCATE(Wave%ur,(nfft*nspinor))
   Wave%ur=czero
   Wave%has_ur=WFD_ALLOCATED
 end if

end subroutine init_wave_0D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/destroy_wave_0D
!! NAME
!!  destroy_wave_0D
!!
!! FUNCTION
!!  Main destruction method for the wave_t datatype.
!!
!! INPUTS
!!  [what]=String defining what has to be freed.
!!     "A" =Both ug and ur and Cprj. Default.
!!     "G" =Only ug.
!!     "R" =Only ur
!!     "C" =Only PAW Cprj.
!!
!! SIDE EFFECTS
!!  Memory in Wave is deallocated depending on what
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine destroy_wave_0D(Wave,what)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_wave_0D'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),optional,intent(in) :: what
 type(wave_t),intent(inout) :: Wave

!Local variables ------------------------------
!scalars
 character(len=10) :: my_what

!************************************************************************

 !@wave_t
 my_what="ALL"; if (PRESENT(what)) my_what=toupper(what)

 if ( .not.starts_with(my_what,(/"A", "G", "R", "C"/) )) then
   MSG_ERROR("unknow what"//TRIM(what))
 end if

 if ( starts_with(my_what,(/"A", "G"/) )) then
   if (associated(Wave%ug))  then
     ABI_DEALLOCATE(Wave%ug)
   end if
   Wave%has_ug=WFD_NOWAVE
 end if

 if ( starts_with(my_what,(/"A", "R"/) )) then
   if (associated(Wave%ur))  then
     ABI_DEALLOCATE(Wave%ur)
   end if
   Wave%has_ur=WFD_NOWAVE
 end if

 if ( starts_with(my_what,(/"A", "C"/) )) then
   if (associated(Wave%Cprj)) then
     call cprj_free(Wave%Cprj)
     ABI_DEALLOCATE(Wave%Cprj)
   end if
   Wave%has_cprj=WFD_NOWAVE
 end if

end subroutine destroy_wave_0D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/destroy_wave_3D
!! NAME
!!  destroy_wave_3D
!!
!! FUNCTION
!!   Destruction method used for a 3-D arrays of wave_t datatyps.
!!
!! INPUTS
!!  Wave(:,:,:)<wave_t>=The array of structures.
!!  [what]=String defining what has to be freed.
!!     "A"=Both ug and ur as PAW Cprj, if any. Default.
!!     "G"  =Only ug.
!!     "R"  =Only ur
!!     "C"  =Only PAW Cprj.
!!  [mask(:,:,:)]=Mask used to select the elements that have to be deallocated. All of them, if not specified.
!!
!! SIDE EFFECTS
!!  Memory in Wave is deallocated depending on what and mask.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine destroy_wave_3D(Wave,what,mask)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_wave_3D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),optional,intent(in) :: what
 type(wave_t),intent(inout) :: Wave(:,:,:)
!arrays
 logical,optional,intent(in) :: mask(:,:,:)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3
 logical :: do_free
 character(len=10) :: my_what

!************************************************************************

 my_what="ALL"; if (PRESENT(what)) my_what=toupper(what)

 do i3=LBOUND(Wave,DIM=3),UBOUND(Wave,DIM=3)
   do i2=LBOUND(Wave,DIM=2),UBOUND(Wave,DIM=2)
     do i1=LBOUND(Wave,DIM=1),UBOUND(Wave,DIM=1)

       do_free=.TRUE.; if (PRESENT(mask)) do_free=mask(i1,i2,i3)
       if (do_free) call destroy_wave_0D(Wave(i1,i2,i3),what=my_what)

     end do
   end do
 end do

end subroutine destroy_wave_3D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/copy_wave_0D
!! NAME
!!  copy_wave_0D
!!
!! FUNCTION
!!  Copy method for the wave_t datatype.
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine copy_wave_0D(Wave_in,Wave_out)

 use defs_basis
 use m_copy

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_wave_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wave_t),intent(in) :: Wave_in
 type(wave_t),intent(inout) :: Wave_out

!Local variables ------------------------------
 integer :: natom,nspinor
 integer :: iatom,ispinor

!************************************************************************

 !@wave_t
 !@cprj_type

 call deep_copy(Wave_in%has_ug,Wave_out%has_ug)
 call deep_copy(Wave_in%has_ur,Wave_out%has_ur)
 call deep_copy(Wave_in%has_cprj,Wave_out%has_cprj)
 call deep_copy(Wave_in%cprj_order,Wave_out%cprj_order)
 call deep_copy(Wave_in%ug,Wave_out%ug)
 call deep_copy(Wave_in%ur,Wave_out%ur)

 natom   = size(Wave_in%Cprj,dim=1)
 nspinor = size(Wave_in%Cprj,dim=2)
 if ((size(Wave_out%Cprj,dim=1) .ne. natom) .or. (size(Wave_out%Cprj,dim=2) .ne. nspinor)) then
   if (associated(Wave_out%Cprj))  then
     ABI_DEALLOCATE(Wave_out%Cprj)
   end if
   ABI_ALLOCATE(Wave_out%Cprj,(natom,nspinor))
 end if

 do ispinor=1,nspinor
   do iatom=1,natom
    Wave_out%Cprj(iatom,ispinor)%ncpgr=Wave_in%Cprj(iatom,ispinor)%ncpgr
    Wave_out%Cprj(iatom,ispinor)%nlmn=Wave_in%Cprj(iatom,ispinor)%nlmn
    call deep_copy(Wave_in%Cprj(iatom,ispinor)%cp,Wave_out%Cprj(iatom,ispinor)%cp)
    call deep_copy(Wave_in%Cprj(iatom,ispinor)%dcp,Wave_out%Cprj(iatom,ispinor)%dcp)
   end do
 end do

end subroutine copy_wave_0D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/copy_wave_3D
!! NAME
!!  copy_wave_3D
!!
!! FUNCTION
!!   Copy method used for a 3-D arrays of wave_t datatyps.
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine copy_wave_3D(Wave_in,Wave_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_wave_3D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wave_t),intent(in) :: Wave_in(:,:,:)
 type(wave_t),intent(inout) :: Wave_out(:,:,:)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3

!************************************************************************

 do i3=LBOUND(Wave_in,DIM=3),UBOUND(Wave_in,DIM=3)
   do i2=LBOUND(Wave_in,DIM=2),UBOUND(Wave_in,DIM=2)
     do i1=LBOUND(Wave_in,DIM=1),UBOUND(Wave_in,DIM=1)
       call copy_wave_0D(Wave_in(i1,i2,i3),Wave_out(i1,i2,i3))
     end do
   end do
 end do

end subroutine copy_wave_3D
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_push_ug
!! NAME
!!  wfd_push_ug
!!
!! FUNCTION
!!  This routine changes the status of the object by saving the wavefunction in the correct
!!  slot inside Wfd%Wave. It also set the corresponding has_ug flag to WFD_STORED.
!!  If the status of the corresponding ur is (WFD_STORED|WFD_ALLOCATED), then an G->R FFT transform
!!  is done (see also update_ur)
!!
!! INPUTS
!!   band=Band index.
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!   Cryst<crystal_structure>=Object defining the unit cell and its symmetries.
!!   ug(Wfd%npwwfn*Wfd%nspinor)=The ug to be saved.
!!   [update_ur]=If .FALSE.: no G-->R transform is done even if ur is (WFD_STORED|WFD_ALLOCATED) so be careful.
!!               Defaults to .TRUE.
!!   [update_cprj]=If .FALSE.: <C|p_i> matrix elements are not recalculatedd even if cprj is (WFD_STORED|WFD_ALLOCATED) so be careful.
!!               Defaults to .TRUE.
!!
!! SIDE EFFECTS
!!   Wfd<wfs_descriptor>=See above.
!!
!! PARENTS
!!      m_shirley,m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_push_ug(Wfd,band,ik_ibz,spin,Cryst,ug,update_ur,update_cprj)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_push_ug'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,band
 logical,optional,intent(in) :: update_ur,update_cprj
 type(wfs_descriptor),intent(inout) :: Wfd
 type(crystal_structure),intent(in) :: Cryst
!arrays
 complex(gwpc),intent(in) :: ug(:)

!Local variables ------------------------------
!scalars
 integer,parameter :: choice1=1,idir0=0,tim_fourdp=5
 integer :: npw_k
 logical :: do_update_ur,do_update_cprj,want_sorted
 character(len=500) :: msg
!arrays
 integer,pointer :: kg_k(:,:),gbound(:,:),igfft0(:)

!************************************************************************

 if (Wfd%debug_level>0) then
   if (.not.wfd_ihave_ug(Wfd,band,ik_ibz,spin)) then
     write(msg,'(a,i0,a,3(i0,1x))')" Node ",Wfd%my_rank," doesn't have ug for (band, ik_ibz, spin): ",band,ik_ibz,spin
     MSG_PERS_ERROR(msg)
   end if
 end if

 if (SIZE(ug)/=Wfd%npwarr(ik_ibz)*Wfd%nspinor) then
   MSG_ERROR("Wrong size in assumed shape array")
 end if

 !@wfs_descriptor
 Wfd%Wave(band,ik_ibz,spin)%ug = ug
 Wfd%Wave(band,ik_ibz,spin)%has_ug = WFD_STORED

 if (Wfd%usepaw==1.and.wfd_ihave_cprj(Wfd,band,ik_ibz,spin)) then ! Update the corresponding cprj if required.
   do_update_cprj=.TRUE.; if (PRESENT(update_cprj)) do_update_cprj=update_cprj
   if (do_update_cprj) then
     want_sorted = (Wfd%Wave(band,ik_ibz,spin)%cprj_order == CPR_SORTED)
     call wfd_ug2cprj(Wfd,band,ik_ibz,spin,choice1,idir0,Wfd%natom,Cryst,Wfd%Wave(band,ik_ibz,spin)%Cprj,sorted=want_sorted)
     Wfd%Wave(band,ik_ibz,spin)%has_cprj = WFD_STORED
   else
     Wfd%Wave(band,ik_ibz,spin)%has_cprj = WFD_ALLOCATED
   end if
 end if

 if (wfd_ihave_ur(Wfd,band,ik_ibz,spin)) then ! Update the corresponding ur if required.
   do_update_ur=.TRUE.; if (PRESENT(update_ur)) do_update_ur=update_ur

   if (do_update_ur) then
     npw_k  =  Wfd%npwarr(ik_ibz)
     kg_k   => Wfd%Kdata(ik_ibz)%kg_k
     gbound => Wfd%Kdata(ik_ibz)%gbound
     igfft0 => Wfd%Kdata(ik_ibz)%igfft0

     call fft_onewfn(Wfd%paral_kgb,Wfd%istwfk(ik_ibz),Wfd%nspinor,npw_k,Wfd%nfft,Wfd%mgfft,Wfd%ngfft,&
&      ug,Wfd%Wave(band,ik_ibz,spin)%ur,igfft0,kg_k,gbound,tim_fourdp,Wfd%MPI_enreg)

     Wfd%Wave(band,ik_ibz,spin)%has_ur = WFD_STORED
   else
     Wfd%Wave(band,ik_ibz,spin)%has_ur = WFD_ALLOCATED
   end if

 end if

end subroutine wfd_push_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_rank_has_ug
!! NAME
!!  wfd_rank_has_ug
!!
!! FUNCTION
!!  This function is used to ask a particular processor whether it has a particular ug and with which status.
!!
!! INPUTS
!!   rank=The MPI rank of the processor.
!!   band=Band index.
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!
!! NOTES
!!   A zero index can be used to inquire the status of a bunch of states.
!!   Thus (band,ik_ibz,spin) = (0,1,1) means: Do you have at least one band for the first k-point and the first spin.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wfd_rank_has_ug(Wfd,rank,band,ik_ibz,spin)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_rank_has_ug'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin,rank
 logical :: wfd_rank_has_ug
 type(wfs_descriptor),intent(in) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: nzeros,bks_flag
!arrays
 integer :: indeces(3)

!************************************************************************

 indeces = (/band,ik_ibz,spin/)
 bks_flag = WFD_STORED

 if ( ALL(indeces/=(/0,0,0/)) ) then
   wfd_rank_has_ug = (Wfd%bks_tab(band,ik_ibz,spin,rank) == bks_flag); RETURN
 else
   nzeros = COUNT(indeces==0)
   if (nzeros==3) MSG_ERROR("All indeces are zero!")

   if (band==0) then
     if (nzeros==1) wfd_rank_has_ug = ANY( Wfd%bks_tab(:,ik_ibz,spin,rank)==bks_flag)
     if (ik_ibz==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(:,:,spin,rank)     ==bks_flag)
     if (spin  ==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(:,ik_ibz,:,rank)   ==bks_flag)

   else if (ik_ibz==0) then
     if (nzeros==1) wfd_rank_has_ug = ANY( Wfd%bks_tab(band,:,spin,rank)==bks_flag)
     if (band  ==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(:,:,spin,rank)   ==bks_flag)
     if (spin  ==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(band,:,:,rank)   ==bks_flag)

   else
     if (nzeros==1) wfd_rank_has_ug = ANY( Wfd%bks_tab(band,ik_ibz,:,rank)==bks_flag)
     if (ik_ibz==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(band,:,:,rank)     ==bks_flag)
     if (band  ==0) wfd_rank_has_ug = ANY( Wfd%bks_tab(:,ik_ibz,:,rank)   ==bks_flag)
   end if
 end if

end function wfd_rank_has_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_ihave_ug
!! NAME
!!  wfd_ihave_ug
!!
!! FUNCTION
!!  This function is used to ask the processor whether it has a particular ug and with which status.
!!
!! INPUTS
!!   band=Band index.
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!   [how]=string defining which status is checked.
!!     Possible mutually exclusive values: "Allocated", "Stored".
!!     Only the first character is checked (no case-sensitive)
!!     By default the function returns .TRUE. if the wave is either WFD_ALLOCATED or WFD_STORED.
!!
!! NOTES
!!   A zero index can be used to inquire the status of a bunch of states.
!!   Thus (band,ik_ibz,spin) = (0,1,1) means: Do you have at least one band for the first k-point and the first spin.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wfd_ihave_ug(Wfd,band,ik_ibz,spin,how)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_ihave_ug'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 logical :: wfd_ihave_ug
 character(len=*),optional,intent(in) :: how
 type(wfs_descriptor),intent(in) :: Wfd

!************************************************************************

 if (PRESENT(how)) then
   wfd_ihave_ug = wfd_ihave(Wfd,"UG",band,ik_ibz,spin,how)
 else
   wfd_ihave_ug = wfd_ihave(Wfd,"UG",band,ik_ibz,spin)
 end if

end function wfd_ihave_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_ihave_ur
!! NAME
!!  wfd_ihave_ur
!!
!! FUNCTION
!!  This function is used to ask the processor whether it has a particular ur and with which status.
!!
!! INPUTS
!!   band=Band index.
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!   [how]=string defining which status is checked. By default the function returns
!!      .TRUE. if the wave is either WFD_ALLOCATED or WFD_STORED.
!!      Possible mutually exclusive values: "Allocated", "Stored".
!!      Only the first character is checked (no case-sensitive)
!!
!! NOTES
!!   A zero index can be used to inquire the status of a bunch of states.
!!   Thus (band,ik_ibz,spin) = (0,1,1) means: Do you have at least one band for the first k-point and the first spin.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wfd_ihave_ur(Wfd,band,ik_ibz,spin,how)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_ihave_ur'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 logical :: wfd_ihave_ur
 character(len=*),optional,intent(in) :: how
 type(wfs_descriptor),intent(in) :: Wfd

!************************************************************************

 if (PRESENT(how)) then
   wfd_ihave_ur = wfd_ihave(Wfd,"UR",band,ik_ibz,spin,how)
 else
   wfd_ihave_ur = wfd_ihave(Wfd,"UR",band,ik_ibz,spin)
 end if

end function wfd_ihave_ur
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_ihave_cprj
!! NAME
!!  wfd_ihave_cprj
!!
!! FUNCTION
!!  This function is used to ask the processor whether it has a particular cprj and with which status.
!!
!! INPUTS
!!   band=Band index.
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!   [how]=string defining which status is checked. By default the function returns
!!      .TRUE. if the wave is either WFD_ALLOCATED or WFD_STORED.
!!      Possible mutually exclusive values: "Allocated", "Stored".
!!      Only the first character is checked (no case-sensitive)
!!
!! NOTES
!!   A zero index can be used to inquire the status of a bunch of states.
!!   Thus (band,ik_ibz,spin) = (0,1,1) means: Do you have at least one band for the first k-point and the first spin.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wfd_ihave_cprj(Wfd,band,ik_ibz,spin,how)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_ihave_cprj'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 logical :: wfd_ihave_cprj
 character(len=*),optional,intent(in) :: how
 type(wfs_descriptor),intent(in) :: Wfd

!************************************************************************

 if (PRESENT(how)) then
   wfd_ihave_cprj = wfd_ihave(Wfd,"CPRJ",band,ik_ibz,spin,how)
 else
   wfd_ihave_cprj = wfd_ihave(Wfd,"CPRJ",band,ik_ibz,spin)
 end if

end function wfd_ihave_cprj
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_ihave
!! NAME
!!  wfd_ihave
!!
!! FUNCTION
!!  This function is used to ask the processor whether it has a particular (ug|ur|cprj) and with which status.
!!
!! INPUTS
!!   band=Band index.
!!   ik_ibz=k-point index
!!   spin=Spin index.
!!   what=String defining what has to be tested.
!!     ug
!!     ur
!!     cprj
!!   [how]=string defining which status is checked.
!!     Possible mutually exclusive values: "Allocated", "Stored".
!!     Only the first character is checked (no case-sensitive)
!!     By default the function returns .TRUE. if the wave is either WFD_ALLOCATED or WFD_STORED.
!!
!! NOTES
!!   A zero index can be used to inquire the status of a bunch of states.
!!   Thus (band,ik_ibz,spin) = (0,1,1) means: Do you have at least one band for the first k-point and the first spin.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wfd_ihave(Wfd,what,band,ik_ibz,spin,how)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_ihave'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 logical :: wfd_ihave
 character(len=*),intent(in) :: what
 character(len=*),optional,intent(in) :: how
 type(wfs_descriptor),intent(in) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: nzeros
 !character(len=500) :: msg
!arrays
 integer :: indeces(3),check(2)
 integer,pointer :: has_flags(:,:,:)

!************************************************************************

 check = (/WFD_ALLOCATED, WFD_STORED/)
 if (PRESENT(how)) then
   if (starts_with(how,(/"A","a"/))) check = (/WFD_ALLOCATED, WFD_ALLOCATED/)
   if (starts_with(how,(/"S","s"/))) check = (/WFD_STORED, WFD_STORED/)
 end if

 indeces = (/band,ik_ibz,spin/)

 select case (toupper(what))
 case ("UG")
   has_flags => Wfd%Wave(:,:,:)%has_ug
 case ("UR")
   has_flags => Wfd%Wave(:,:,:)%has_ur
 case ("CPRJ")
   has_flags => Wfd%Wave(:,:,:)%has_cprj
 case default
   MSG_ERROR("Wrong what"//TRIM(what))
 end select

 if ( ALL(indeces/=(/0,0,0/)) ) then
   wfd_ihave = ( ANY(has_flags(band,ik_ibz,spin) == check )); RETURN
 else
   nzeros = COUNT(indeces==0)
   if (nzeros==3) MSG_ERROR("All indeces are zero!")

   if (band==0) then
     if (nzeros==1) wfd_ihave = ANY( has_flags(:,ik_ibz,spin)==check(1) .or.&
&                                    has_flags(:,ik_ibz,spin)==check(2) )

     if (ik_ibz==0) wfd_ihave = ANY( has_flags(:,:,spin)==check(1) .or.&
&                                    has_flags(:,:,spin)==check(2) )

     if (spin  ==0) wfd_ihave = ANY( has_flags(:,ik_ibz,:)==check(1) .or.&
&                                    has_flags(:,ik_ibz,:)==check(2) )

   else if (ik_ibz==0) then
     if (nzeros==1) wfd_ihave = ANY( has_flags(band,:,spin)==check(1) .or.&
&                                    has_flags(band,:,spin)==check(2) )

     if (band  ==0) wfd_ihave = ANY( has_flags(:,:,spin)==check(1) .or.&
&                                    has_flags(:,:,spin)==check(2) )

     if (spin  ==0) wfd_ihave = ANY( has_flags(band,:,:)==check(1) .or.&
&                                    has_flags(band,:,:)==check(2) )
   else
     if (nzeros==1) wfd_ihave = ANY( has_flags(band,ik_ibz,:)==check(1) .or.&
&                                    has_flags(band,ik_ibz,:)==check(2) )

     if (ik_ibz==0) wfd_ihave = ANY( has_flags(band,:,:)==check(1) .or.&
&                                    has_flags(band,:,:)==check(2) )

     if (band  ==0) wfd_ihave = ANY( has_flags(:,ik_ibz,:)==check(1) .or.&
&                                    has_flags(:,ik_ibz,:)==check(2) )
   end if
 end if

end function wfd_ihave
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_mybands
!! NAME
!!  wfd_mybands
!!
!! FUNCTION
!!  Return the list of band indeces of the ug owned by this node at given (k,s).
!!
!! INPUTS
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!  [how]=string defining which status is checked.
!!    Possible mutually exclusive values: "Allocated", "Stored".
!!    Only the first character is checked (no case-sensitive)
!!    By default the list of bands whose status is either WFD_ALLOCATED or WFD_STORED is returned.
!!
!! OUTPUT
!!  how_manyb=The number of bands owned by this node
!!  my_band_list(Wfd%mband)=The first how_manyb values are the bands treated by this node.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_mybands(Wfd,ik_ibz,spin,how_manyb,my_band_list,how)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_mybands'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin
 integer,intent(out) :: how_manyb
 character(len=*),optional,intent(in) :: how
 type(wfs_descriptor),intent(in) :: Wfd
!arrays
 integer,intent(out) :: my_band_list(Wfd%mband)

!Local variables ------------------------------
!scalars
 integer :: band
 logical :: do_have

!************************************************************************

 how_manyb=0; my_band_list=-1
 do band=1,Wfd%nband(ik_ibz,spin)
   if (PRESENT(how)) then
     do_have = wfd_ihave_ug(Wfd,band,ik_ibz,spin,how=how)
   else
     do_have = wfd_ihave_ug(Wfd,band,ik_ibz,spin)
   end if
   if (do_have) then
     how_manyb = how_manyb +1
     my_band_list(how_manyb)=band
   end if
 end do

end subroutine wfd_mybands
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_bands_of_proc
!! NAME
!!  wfd_bands_of_proc
!!
!! FUNCTION
!!  Return the list of band index of the ug owned by a given processor at given (k,s).
!!
!! INPUTS
!!  Wfd
!!  rank=The MPI rank of the processor.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!
!! OUTPUT
!!  how_manyb=The number of bands owned by this node
!!  rank_band_list(Wfd%mband)=The first how_manyb values are the bands treated by the node.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wfd_bands_of_rank(Wfd,rank,ik_ibz,spin,how_manyb,rank_band_list)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_bands_of_rank'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,rank
 integer,intent(out) :: how_manyb
 type(wfs_descriptor),intent(in) :: Wfd
!arrays
 integer,intent(out) :: rank_band_list(Wfd%mband)

!Local variables ------------------------------
!scalars
 integer :: band
 logical :: it_has

!************************************************************************

 how_manyb=0; rank_band_list=-1
 do band=1,Wfd%nband(ik_ibz,spin)
   it_has = wfd_rank_has_ug(Wfd,rank,band,ik_ibz,spin)
   if (it_has) then
     how_manyb = how_manyb +1
     rank_band_list(how_manyb)=band
   end if
 end do

end subroutine wfd_bands_of_rank
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_get_ug
!! NAME
!!  wfd_get_ug
!!
!! FUNCTION
!!  Get a copy of a wave function in reciprocal space.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=the data type
!!  band=the index of the band.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!
!! OUTPUT
!!  ug(Wfd%npwwfn*Wfd%nspinor)=The required wavefunction in reciprocal space.
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_get_ug(Wfd,band,ik_ibz,spin,ug)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_get_ug'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 complex(gwpc),intent(out) :: ug(Wfd%npwarr(ik_ibz)*Wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer :: npw_k
 character(len=500) :: msg
!************************************************************************

 if (wfd_ihave_ug(Wfd,band,ik_ibz,spin,"Stored")) then
   npw_k = Wfd%npwarr(ik_ibz)
   call xcopy(npw_k*Wfd%nspinor,Wfd%Wave(band,ik_ibz,spin)%ug,1,ug,1)
 else
   write(msg,'(a,i0,a,3i0)')" Node ",Wfd%my_rank," doesn't have (band,ik_ibz,spin)=",band,ik_ibz,spin
   MSG_PERS_BUG(msg)
 end if

end subroutine wfd_get_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_ptr_ug
!! NAME
!!  wfd_ptr_ug
!!
!! FUNCTION
!!  Returns a pointer to ug
!!  WARNING: Do not use the returned pointer to modify the location of memory.
!!   The status of the object should always be modified via the appropriate method.
!!   Use the pointer only if you want to avoid a copy and you are not going to change the ug!
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=the data type
!!  band=the index of the band.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!
!! OUTPUT
!!  wfd_ptr_ug
!!  ierr=Status error.
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_ptr_ug(Wfd,band,ik_ibz,spin,ptr_ug,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_ptr_ug'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 integer,intent(out) :: ierr
 type(wfs_descriptor),intent(in) :: Wfd
!arrays
 complex(gwpc),pointer :: ptr_ug(:)

!Local variables ------------------------------
!scalars
 !character(len=500) :: msg

!************************************************************************

 if (wfd_ihave_ug(Wfd,band,ik_ibz,spin,how="Stored")) then
   ierr=0
   ptr_ug => Wfd%Wave(band,ik_ibz,spin)%ug
 else
   !write(msg,'(a,i0,a,3(i0,1x))')" Node ",Wfd%my_rank," doesn't have ug for (band, ik_ibz, spin): ",band,ik_ibz,spin
   !MSG_ERROR(msg)
   ierr=1
   nullify(ptr_ug)
 end if

end subroutine wfd_ptr_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_ptr_ur
!! NAME
!!  wfd_ptr_ur
!!
!! FUNCTION
!!  Returns a pointer to ur
!!  WARNING: Do not use the returned pointer to modify the location of memory.
!!   The status of the object should always be modified via the appropriate method.
!!   Use the pointer only if you want to avoid a copy and you are not going to change the ug!
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=the data type
!!  band=the index of the band.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!
!! OUTPUT
!!  wfd_ptr_ur
!!  ierr=Status error.
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_ptr_ur(Wfd,band,ik_ibz,spin,ptr_ur,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_ptr_ur'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 integer,intent(out) :: ierr
 type(wfs_descriptor),intent(in) :: Wfd
!arrays
 complex(gwpc),pointer :: ptr_ur(:)

!Local variables ------------------------------
!scalars
 !character(len=500) :: msg

!************************************************************************

 if (wfd_ihave_ur(Wfd,band,ik_ibz,spin,how="Stored")) then
   ptr_ur => Wfd%Wave(band,ik_ibz,spin)%ur
   ierr=0
 else
   !write(msg,'(a,i0,a,3(i0,1x))')" Node ",Wfd%my_rank," doesn't have ug for (band, ik_ibz, spin): ",band,ik_ibz,spin
   !MSG_ERROR(msg)
   ierr=1
   nullify(ptr_ur)
 end if

end subroutine wfd_ptr_ur
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_wave_free
!! NAME
!!  wfd_wave_free
!!
!! FUNCTION
!!  Free the set of waves specified by mask.
!!
!! INPUTS
!!  mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)=.TRUE. if the memory allocated for
!!    this state has to be freed
!!  [what]=String specifying which array have to be deallocated.
!!    Possible values (no case-sensitive).
!!      "All"= To free both ug and ur and PAW Cprj, if any. Default
!!      "G"  = Only ug
!!      "R"  = Only ur.
!!      "C"  = Only PAW Cprj.
!!
!! SIDE EFFECTS
!!  Wfd<wfs_descriptor>=See above.
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_wave_free(Wfd,what,bks_mask)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_wave_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(inout) :: Wfd
 character(len=*),optional,intent(in) :: what
!arrays
 logical,optional,intent(in) :: bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,band !,ierr,istat
 logical :: do_free
 character(len=10) :: my_what
!************************************************************************

 my_what="ALL"; if (PRESENT(what)) my_what=toupper(what)

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)

        do_free=.TRUE.; if (PRESENT(bks_mask)) do_free=bks_mask(band,ik_ibz,spin)
        if (do_free) then
          call destroy_wave_0D(Wfd%Wave(band,ik_ibz,spin),what=my_what)
          if ( starts_with(my_what,(/"A", "G"/) )) then ! Update the associated flags.
            Wfd%bks_tab(band,ik_ibz,spin,Wfd%my_rank) = WFD_NOWAVE
          end if
        end if

     end do
   end do
 end do
 !
 ! Reinit the MPI communicators.
 call wfd_set_mpicomm(Wfd)

end subroutine wfd_wave_free
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_who_has_ug
!! NAME
!!  wfd_who_has_ug
!!
!! FUNCTION
!!  Return the number of processors having a particular (b,k,s) state as well as their MPI rank.
!!  Warning: Wfd%bks_tab is supposed to be up-to-date (see wfd_update_bkstab).
!!
!! INPUTS
!!  band=the index of the band.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!
!! OUTPUT
!!  how_many=The number of nodes owing this ug state.
!!  proc_ranks(1:how_many)=Gives the MPI rank of the nodes owing the state.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_who_has_ug(Wfd,band,ik_ibz,spin,how_many,proc_ranks)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_who_has_ug'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 integer,intent(out) :: how_many
 type(wfs_descriptor),intent(in) :: Wfd
!arrays
 integer,intent(out) :: proc_ranks(Wfd%nproc)

!Local variables ------------------------------
!scalars
 integer :: irank
 logical :: bks_select,spin_select,kpt_select
 character(len=500) :: msg
!arrays

!************************************************************************

 bks_select  = (band/=0.and.ik_ibz/=0.and.spin/=0)
 spin_select = (band==0.and.ik_ibz==0.and.spin/=0)
 kpt_select = (band==0.and.ik_ibz/=0.and.spin/=0)

 how_many=0; proc_ranks=-1

 if (bks_select) then ! List the proc owining this (b,k,s) state.

   do irank=0,Wfd%nproc-1
     if (Wfd%bks_tab(band,ik_ibz,spin,irank)==WFD_STORED) then
       how_many = how_many +1
       proc_ranks(how_many)=irank
     end if
   end do

 else if (spin_select) then ! List the proc owining at least one state with this spin.

   do irank=0,Wfd%nproc-1
     if ( ANY(Wfd%bks_tab(:,:,spin,irank)==WFD_STORED) ) then
       how_many = how_many +1
       proc_ranks(how_many)=irank
     end if
   end do

 else if (kpt_select) then ! List the proc owining at least one state with this (k-point, spin).

   do irank=0,Wfd%nproc-1
     if ( ANY(Wfd%bks_tab(:,ik_ibz,spin,irank)==WFD_STORED) ) then
       how_many = how_many +1
       proc_ranks(how_many)=irank
     end if
   end do

 else
   write(msg,'(a,3(i0,1x))')" Wrong value for (b,k,s) ",band,ik_ibz,spin
   MSG_ERROR(msg)
 end if

end subroutine wfd_who_has_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_everybody_has_ug
!! NAME
!!  wfd_everybody_has_ug
!!
!! FUNCTION
!!  Return .TRUE. if all the nodes inside comm own the specified ug state.
!!
!! INPUTS
!!  band=the index of the band.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wfd_everybody_has_ug(Wfd,band,ik_ibz,spin) result(answer)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_everybody_has_ug'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 logical :: answer
 type(wfs_descriptor),intent(in) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: how_many,nzeros,ib
!arrays
 integer :: proc_ranks(Wfd%nproc),indeces(3)

!************************************************************************

 indeces = (/band,ik_ibz,spin/)

 if ( ALL(indeces/=(/0,0,0/)) ) then
   call wfd_who_has_ug(Wfd,band,ik_ibz,spin,how_many,proc_ranks)
   answer = (how_many==Wfd%nproc); RETURN
 else
   nzeros = COUNT(indeces==0)
   if (nzeros==3) MSG_ERROR("All indeces are zero!")

   answer=.TRUE.
   MSG_WARNING("Some cases are not coded!") ! TODO

   if (band==0) then

     if (nzeros==1)  then     ! All the bands for the given k-point and spin?
       ib=0
       do while(answer.and.ib<Wfd%nband(ik_ibz,spin))
         ib=ib+1
         call wfd_who_has_ug(Wfd,ib,ik_ibz,spin,how_many,proc_ranks)
         answer = (how_many==Wfd%nproc)
       end do; RETURN

     else if (ik_ibz==0) then ! All the bands and all the k-points for the the given spin?

     else if (spin==0) then   ! All the bands and all the spins for the given k-point?

     end if

   else if (ik_ibz==0) then
     if (nzeros==1) then     ! All the k-points for the the given band and spin?

     else if (band==0) then  ! All the k-points and all the bands for the the given spin?

     else if (spin==0) then  ! All the k-points and all the spins for the the given band?

     end if

   else
     if (nzeros==1) then      ! All the spins for the the given band and k-point?

     else if (ik_ibz==0) then ! All the spins and all the k-points for the the given band?

     else if (band==0) then   ! All the spins and all the bands for the the given k-point?

     end if
   end if

   MSG_ERROR("Not implemented error")
 end if

end function wfd_everybody_has_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_update_bkstab
!! NAME
!!  wfd_update_bkstab
!!
!! FUNCTION
!!  This routine should be called by all the nodes before any MPI operation involving the object.
!!  It updates the bks_tab storing information on the distribution of ug.
!!
!! SIDE EFFECTS
!!  Wfs%bks_tab
!!
!! PARENTS
!!      gw_tools,m_io_kss,m_wfs,wfd_mkrho
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_update_bkstab(Wfd)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_update_bkstab'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(inout) :: Wfd

!Local variables ------------------------------
!scalars
 integer :: ierr,nelem
 integer,allocatable :: my_vtab(:),gather_vtabs(:)

!************************************************************************

 ! Fill my slice of the global table.
 Wfd%bks_tab(:,:,:,Wfd%my_rank) = Wfd%Wave(:,:,:)%has_ug

 ! Gather flags on each node.
 nelem=Wfd%mband*Wfd%nkibz*Wfd%nsppol
 ABI_ALLOCATE(my_vtab,(nelem))
 my_vtab = RESHAPE(Wfd%bks_tab(:,:,:,Wfd%my_rank),(/nelem/))

 ABI_ALLOCATE(gather_vtabs,(nelem*Wfd%nproc))

 call xallgather_mpi(my_vtab,nelem,gather_vtabs,Wfd%comm,ierr)

 Wfd%bks_tab = RESHAPE(gather_vtabs,(/Wfd%mband,Wfd%nkibz,Wfd%nsppol,Wfd%nproc/))
 ABI_DEALLOCATE(my_vtab)
 ABI_DEALLOCATE(gather_vtabs)

end subroutine wfd_update_bkstab
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_set_mpicomm
!! NAME
!!  wfd_set_mpicomm
!!
!! FUNCTION
!!
!! PARENTS
!!      m_io_kss,m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_set_mpicomm(Wfd)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_set_mpicomm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(inout) :: Wfd

!Local variables ------------------------------
!scalars
#ifdef HAVE_MPI
 integer :: ik_ibz,spin,band,ierr,how_many,spin_comm,kpt_comm
 integer :: world_group,bks_group,kpt_group,spin_group
 character(len=500) :: msg
!arrays
 integer :: proc_ranks(Wfd%nproc)
#endif

!************************************************************************

#ifdef HAVE_MPI
 !
 ! First free the old communicators.
 call xcomm_free(Wfd%bks_comm)
 !
 ! Update the bks_tab.
 call wfd_update_bkstab(Wfd)

 call MPI_COMM_GROUP(Wfd%comm,world_group,ierr)
 !
 ! Init spin communicators.
 do spin=1,Wfd%nsppol
   !
   ! The list of procs owining at least one state with this spin.
   call wfd_who_has_ug(Wfd,0,0,spin,how_many,proc_ranks)

   if (how_many>0) then
     !write(std_out,*)" before incl ",how_many,proc_ranks(1:how_many)
     call MPI_GROUP_INCL(world_group,how_many,proc_ranks,spin_group,ierr)

     call MPI_COMM_CREATE(Wfd%comm,spin_group,spin_comm,ierr)

     Wfd%bks_comm(0,0,spin) = spin_comm

     !write(std_out,*)" after create ",world_group,spin_group,Wfd%spin_comm(spin)
     call MPI_GROUP_FREE(spin_group,ierr)
   else
     if (Wfd%debug_level>0) then
       write(msg,'(a,i0)')" Nobody has spin: ",spin
       MSG_WARNING(msg)
     end if
     Wfd%bks_comm(0,0,spin) = xmpi_comm_null
   end if

 end do
 !
 ! This section is problematic since several MPI implementations have a limit
 ! in the number of communicators that can be created at run-time.
 ! Disabled for the time being.
if (.FALSE.) then
 ! Init (k,s) communicators.
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     !
     ! The list of procs owining at least one state with this spin and this k-point.
     call wfd_who_has_ug(Wfd,0,ik_ibz,spin,how_many,proc_ranks)
     if (how_many>0) then
       call MPI_GROUP_INCL(world_group,how_many,proc_ranks,kpt_group,ierr)

       call MPI_COMM_CREATE(Wfd%comm,kpt_group,kpt_comm,ierr)

       Wfd%bks_comm(0,ik_ibz,spin) = kpt_comm

       call MPI_GROUP_FREE(kpt_group,ierr)
     else
       if (Wfd%debug_level>0) then
         write(msg,'(a,2i0)')" Nobody has (kpt, spin): ",ik_ibz,spin
         MSG_WARNING(msg)
       end if
       Wfd%bks_comm(0,ik_ibz,spin) = xmpi_comm_null
     end if
   end do
 end do
 !
 ! Init (b,k,s) communicators.
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)
       !
       ! The list of procs owining this (b,k,s) state.
       call wfd_who_has_ug(Wfd,band,ik_ibz,spin,how_many,proc_ranks)
       if (how_many>0) then
         call MPI_GROUP_INCL(world_group,how_many,proc_ranks,bks_group,ierr)

         call MPI_COMM_CREATE(Wfd%comm,bks_group,Wfd%bks_comm(band,ik_ibz,spin),ierr)

         call MPI_GROUP_FREE(bks_group,ierr)
       else
         if (Wfd%debug_level>0) then
           write(msg,'(a,3(i0,1x))')" Nobody has (band, kpt, spin): ",band,ik_ibz,spin
           MSG_WARNING(msg)
         end if
         Wfd%bks_comm(band,ik_ibz,spin) = xmpi_comm_null
       end if

     end do
   end do
 end do
end if

 call MPI_GROUP_FREE(world_group,ierr)
#endif

 RETURN
 ABI_UNUSED(Wfd%nkibz)

end subroutine wfd_set_mpicomm
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_distribute_bands
!! NAME
!!  wfd_distribute_bands
!!
!! FUNCTION
!!  This routines distributes a set of band indeces taking into account the
!!  distribution of the ug.
!!
!! INPUTS
!!  band=the index of the band.
!!  ik_ibz=Index of the k-point in the IBZ
!!  spin=spin index
!!  [got(Wfd%nproc)]=The number of tasks already assigned to the nodes.
!!  [bmask(Wfd%mband)]=The routine will raise an error if one band index
!!    is not treated by any processor. bmask can be used to select the subset of
!!    indeces that are expected to be available.
!!
!! OUTPUT
!!   my_nband=The number of bands that will be treated by this node.
!!   my_band_list(1:my_nband)=The band indeces for this node
!!
!! PARENTS
!!      cchi0q0_intraband,gw_tools,m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,got,bmask)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_distribute_bands'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin
 integer,intent(out) :: my_nband
 type(wfs_descriptor),intent(in) :: Wfd
!arrays
 integer,intent(out) :: my_band_list(Wfd%mband)
 integer,optional,intent(inout) :: got(Wfd%nproc)
 logical,optional,intent(in) :: bmask(Wfd%mband)

!Local variables ------------------------------
!scalars
 integer :: band,how_many,idle
 character(len=500) :: msg
!arrays
 integer :: proc_ranks(Wfd%nproc),get_more(Wfd%nproc)
 logical :: rank_mask(Wfd%nproc)

!************************************************************************

 my_nband=0; my_band_list=0
 get_more=0; if (PRESENT(got)) get_more = got

 do band=1,Wfd%nband(ik_ibz,spin)
   if (PRESENT(bmask)) then
     if (.not.bmask(band)) CYCLE
   end if

   call wfd_who_has_ug(Wfd,band,ik_ibz,spin,how_many,proc_ranks)

   if (how_many==1) then ! I am the only one owing this band. Add it to list.
     if (proc_ranks(1) == Wfd%my_rank) then
       my_nband=my_nband + 1
       my_band_list(my_nband) = band
     end if
   else if (how_many>1) then  ! This band is duplicated. Assign it trying to obtain a good load distribution.
     rank_mask=.FALSE.; rank_mask(proc_ranks(1:how_many)+1)=.TRUE.
     idle = imin_loc(get_more,mask=rank_mask)
     get_more(idle) = get_more(idle) + 1
     if (Wfd%my_rank==idle-1) then
       my_nband=my_nband + 1
       my_band_list(my_nband) = band
     end if
   else
     write(msg,'(a,3(i0,1x))')" No processor has (band, ik_ibz, spin): ",band,ik_ibz,spin
     MSG_PERS_ERROR(msg)
   end if
 end do

 if (PRESENT(got)) got = get_more

end subroutine wfd_distribute_bands
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_rotate
!! NAME
!! wfd_rotate
!!
!! FUNCTION
!!  This routine performs a linear trasformation of the wavefunctions stored in Wfd
!!  taking into account memory distribution. The transformation is done in reciprocal
!!  space therefore all the ug should be available. Wavefunctions in real space are then
!!  obtained via FFT. The implementation assumes that the matrix associated to the
!!  linear transformation is sparse (No BLAS-3 calls here).
!!
!! INPUTS
!!  Cryst<crystal_structure>=Object defining the unit cell and its symmetries.
!!  m_lda_to_qp(mband,mband,nkibz,nsppol)= expansion of the QP amplitudes in terms of KS wavefunctions.
!!  [bmask(mband,nkibz,nsppol)]=The routine will raise an error if one band index
!!    is not treated by any processor. bmask can be used to select the subset of
!!    indeces that are expected to be available.
!!
!! SIDE EFFECTS
!!   Wfd<wfs_descriptor>=See above.
!!
!! PARENTS
!!      bethe_salpeter,screening,sigma
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE
!!

subroutine wfd_rotate(Wfd,Cryst,m_lda_to_qp,bmask)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_rotate'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(inout) :: Wfd
 type(crystal_structure),intent(in) :: Cryst
!arrays
 complex(dpc),target,intent(in) :: m_lda_to_qp(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol)
 logical,optional,intent(in) :: bmask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)

!Local variables-------------------------------
!scalars
 integer :: band,ik_ibz,spin,ierr,icol,nnew,inew,my_nband,ib,npw_k
 !character(len=500) :: msg
!arrays
 integer :: new_list(Wfd%mband),my_band_list(Wfd%mband)
 complex(dpc),pointer :: umat_sk(:,:)
 complex(gwpc) :: mcol(Wfd%mband)
 complex(gwpc),allocatable :: new_ug(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ! Update the distribution table, first.
 call wfd_update_bkstab(Wfd)
 !
 ! === Calculate : $\Psi^{QP}_{r,b} = \sum_n \Psi^{KS}_{r,n} M_{n,b}$ ===
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz

     npw_k  = Wfd%npwarr(ik_ibz)
     umat_sk => m_lda_to_qp(:,:,ik_ibz,spin)
     !
     ! Select only those states that are mixed by the (sparse) m_lda_to_qp.
     nnew=0; new_list=0
     do icol=1,Wfd%nband(ik_ibz,spin)
       mcol = m_lda_to_qp(:,icol,ik_ibz,spin)
       mcol(icol) = mcol(icol) - cone
       if (ANY(ABS(mcol)>tol12)) then  ! Avoid a simple copy.
         nnew=nnew+1
         new_list(nnew)=icol
       end if
     end do
     if (nnew==0) CYCLE ! Nothing to do.
     !
     ! Retrieve the set of band indeces that have to be treated by
     ! this node taking into accoung a possible duplication.
     !
     if (PRESENT(bmask)) then
       call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,bmask=bmask(:,ik_ibz,spin))
     else
       call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list)
     end if

     !if (my_nband>0) then
     !  write(std_out,*)" At (ik_ibz,spin) ",ik_ibz,spin,&
     !  & ", rank ",Wfd%my_rank," will sum ",my_nband," bands, my_band_list: ",my_band_list(1:my_nband)
     !end if

     ABI_ALLOCATE(new_ug,(npw_k*Wfd%nspinor,nnew))
     new_ug=czero
     do inew=1,nnew
       icol = new_list(inew)
       do ib=1,my_nband
         band = my_band_list(ib)
         if (ABS(umat_sk(band,icol))>tol12) then
           new_ug(:,inew) = new_ug(:,inew) + umat_sk(band,icol)* Wfd%Wave(band,ik_ibz,spin)%ug
         end if
       end do
     end do

     call xsum_mpi(new_ug,Wfd%comm,ierr)
     call xbarrier_mpi(Wfd%comm)
     ! =======================================
     ! === Update the input wave functions ===
     ! =======================================
     do inew=1,nnew
       band = new_list(inew)
       if (wfd_ihave_ug(Wfd,band,ik_ibz,spin)) call wfd_push_ug(Wfd,band,ik_ibz,spin,Cryst,new_ug(:,inew))
     end do

     ABI_DEALLOCATE(new_ug)

   end do !ik_ibz
 end do !spin

 !This is needed only if FFTs are not done in wfd_push_ug. Do not know which one is faster.
 !call wfd_reset_ur_cprj(Wfd)

 call xbarrier_mpi(Wfd%comm)

 DBG_EXIT("COLL")

end subroutine wfd_rotate
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_iterator_bks
!! NAME
!!  wfd_iterator_bks
!!
!! FUNCTION
!!  This routines returns an iterator used to loop over bands, k-points and spin indeces
!!  taking into account the distribution of the ug.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=
!!  bks_mask(Wfd%mband.Wfd%nkibz,Wfd%nsppol)= mask used to select the (b,k,s) indeces.
!!
!! OUTPUT
!!  iter_bks<iter2_t>=Iterator over the bands treated by this node for each k-point and spin.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wfd_iterator_bks(Wfd,bks_mask) result(iter_bks)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_iterator_bks'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(in) :: Wfd
!arrays
 logical,optional,intent(in) :: bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)
 !type(coeffi1_type),intent(out) :: iter_bks(Wfd%nkibz,Wfd%nsppol)
 type(iter2_t) :: iter_bks

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,my_nband
 !character(len=500) :: msg
!arrays
 integer :: my_band_list(Wfd%mband)

!************************************************************************

 call iter_alloc(iter_bks,(/Wfd%nkibz,Wfd%nsppol/))

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz

     if (PRESENT(bks_mask)) then
       call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,bmask=bks_mask(:,ik_ibz,spin))
     else
       call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list)
     end if

     !iter_bks(ik_ibz,spin)%size = my_nband
     !allocate(iter_bks(ik_ibz,spin)%value(my_nband))
     !if (my_nband>0) iter_bks(ik_ibz,spin)%value = my_band_list(1:my_nband)

     call iter_push(iter_bks,ik_ibz,spin,my_band_list(1:my_nband))
   end do
 end do

end function wfd_iterator_bks
!!***

!----------------------------------------------------------------------


!!****f* m_wfs/wfd_bks_distrb
!! NAME
!!  wfd_bks_distrb
!!
!! FUNCTION
!!  This routines build a local logical table indexed by bands, k-points and spin that defines
!!  the distribution of the load inside the loops according to the availability of the ug.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=
!!  [bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)]=Mask used to skip selecter (b,k,s) entries.
!!  [got(Wfd%nproc)]=The number of tasks already assigned to the nodes.
!!
!! OUTPUT
!!  bks_distrbk(Wfd%mband,Wfd%nkibz,Wfd%nsppol)=Global table with the rank of the node treating (b,k,s)
!!
!! PARENTS
!!      wfd_pawrhoij
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_bks_distrb(Wfd,bks_distrb,got,bks_mask)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_bks_distrb'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(in) :: Wfd
!arrays
 integer,intent(out) :: bks_distrb(Wfd%mband,Wfd%nkibz,Wfd%nsppol)
 integer,optional,intent(inout) :: got(Wfd%nproc)
 logical,optional,intent(in) :: bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,band,how_many,idle
 character(len=500) :: msg
!arrays
 integer :: get_more(Wfd%nproc),proc_ranks(Wfd%nproc)
 logical :: rank_mask(Wfd%nproc)

!************************************************************************

 get_more=0; if (PRESENT(got)) get_more=got

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)

       if (PRESENT(bks_mask)) then
         if (.not.bks_mask(band,ik_ibz,spin)) CYCLE
       end if

       call wfd_who_has_ug(Wfd,band,ik_ibz,spin,how_many,proc_ranks)

       if (how_many==1) then ! I am the only one owing this band. Add it to list.
         bks_distrb(band,ik_ibz,spin) = proc_ranks(1)

       else if (how_many>1) then ! This band is duplicated. Assign it trying to obtain a good load distribution.
         rank_mask=.FALSE.; rank_mask(proc_ranks(1:how_many)+1)=.TRUE.
         idle = imin_loc(get_more,mask=rank_mask)
         get_more(idle) = get_more(idle) + 1
         bks_distrb(band,ik_ibz,spin) = proc_ranks(idle)

       else
         call wfd_dump_errinfo(Wfd)
         write(msg,'(a,3(i0,1x))')" Nobody has (band, ik_ibz, spin): ",band,ik_ibz,spin
         MSG_ERROR(msg)
       end if

     end do
   end do
 end do

 if (PRESENT(got)) got=get_more

end subroutine wfd_bks_distrb
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_sanity_check
!! NAME
!!  wfd_sanity_check
!!
!! FUNCTION
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_sanity_check(Wfd)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_sanity_check'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,band,mpi_ierr,ierr
 integer :: how_manyb,unt_dbg,irank
 character(len=500) :: msg
!arrays
 integer :: my_band_list(Wfd%mband)
 integer,pointer :: my_bkstab(:,:,:)

!************************************************************************

 call wfd_update_bkstab(Wfd)

 my_bkstab => Wfd%bks_tab(:,:,:,Wfd%my_rank)
 ierr=0

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
      do band=1,Wfd%nband(ik_ibz,spin)
        if (my_bkstab(band,ik_ibz,spin) == WFD_STORED .and. .not. wfd_ihave_ug(Wfd,band,ik_ibz,spin,"Stored") ) then
          write(msg,'(a,3(i0,1x))')" Found inconsistency in bks_tab for (band, ik_ibz, spin): ",band,ik_ibz,spin
          call wrtout(std_out,msg,"PERS")
          ierr=ierr+1
        end if
     end do
   end do
 end do

 call xsum_mpi(ierr,Wfd%comm,mpi_ierr)

 if (ierr/=0) then
   unt_dbg = get_unit()
   open(unit=unt_dbg,file="WFD_DEBUG")
   do irank=0,Wfd%nproc-1

     if (irank==Wfd%my_rank) then
       write(unt_dbg,*)" (k,b,s) states owned by rank: ",Wfd%my_rank

       do spin=1,Wfd%nsppol
         do ik_ibz=1,Wfd%nkibz
            write(unt_dbg,*)" (spin,ik_ibz) ",spin,ik_ibz
            call wfd_mybands(Wfd,ik_ibz,spin,how_manyb,my_band_list,"Stored")
            write(unt_dbg,*) (my_band_list(band),band=1,how_manyb)
          end do
       end do

     end if
   end do
   close(unt_dbg)
   call xbarrier_mpi(Wfd%comm)
   MSG_ERROR("Sanity check failed. Check WFD_DEBUG")
 end if

end subroutine wfd_sanity_check
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_dump_errinfo
!! NAME
!!  wfd_dump_errinfo
!!
!! FUNCTION
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=
!!
!! OUTPUT
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_dump_errinfo(Wfd,onfile)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_dump_errinfo'
 use interfaces_27_toolbox_oop
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,optional,intent(in) :: onfile
 type(wfs_descriptor),intent(in) :: Wfd
!arrays

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,band
 integer :: how_manyb,unt_dbg
 character(len=10) :: strank
 !character(len=500) :: msg
 character(len=fnlen) :: fname_dbg
!arrays
 integer :: my_band_list(Wfd%mband)

!************************************************************************

 unt_dbg=std_out

 if (PRESENT(onfile)) then
   if (onfile) then
     call int2char(Wfd%my_rank,strank)
     fname_dbg = "WFD_DEBUG_RANK"//TRIM(strank); unt_dbg = get_unit()
     open(unit=unt_dbg,file=fname_dbg)
   end if
 end if

 write(unt_dbg,*)" (k,b,s) states owned by rank: ",Wfd%my_rank
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
      write(unt_dbg,*)" ug stored at (ik_ibz, spin) ",ik_ibz,spin
      call wfd_mybands(Wfd,ik_ibz,spin,how_manyb,my_band_list,"Stored")
      write(unt_dbg,*) (my_band_list(band),band=1,how_manyb)
    end do
 end do

end subroutine wfd_dump_errinfo
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_gather_g2k
!! NAME
!! wfd_gather_g2k
!!
!! FUNCTION
!!  Converts a block of wavefunctions from the gamma-centered basis set to the k-centered
!!  basis set used in abinit. It works for a single block of states at fixed spin and k-point.
!!  gathering the data on each node.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=Structure containing the wave functions for the GW.
!!    %ecut=cutoff energy for GW wavefunctions.
!!    %kibz(3,nkibz)=K-points in reduced coordinates
!!    %istwfk(nkibz)=Storage mode for each k-point
!!    %npwwfn=Number of G-vectors in the gamma-centered basis set.
!!    %nspinor=Number of spinorial components.
!!    %gvec(3,npwwfn)=G-vectors in the gamma-centered basis set.
!!  ik_ibz=Index of the required k-point
!!  spin=Required spin index.
!!  ikg=Shift to be used in cg and. Mainly used to fill selected slices of the arrays.
!!  gmet(3,3)=Metric in reciprocal space.
!!
!! OUTPUT
!!  nmiss=Number of missing G-vectors, i.e. G-vectors in the k-centered basis set not contained
!!        in the gamma-centered basis set used to describe the GW wavefunctions.
!!  See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  cg(:,:)= cg(icg:) contains the states with Fourier components defined in the k-centered basis set.
!!  kg_k(:,:)= contains the G-vectors centered on the k-point.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_gather_g2k(Wfd,ik_ibz,spin,ikg,kg_k,cg,gmet,nmiss)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_gather_g2k'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,ikg
 integer,intent(out) :: nmiss
 type(wfs_descriptor),intent(in) :: Wfd
!arrays
 integer,pointer :: kg_k(:,:)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(inout)  :: cg(:,:)

!Local variables-------------------------------
!scalars
 integer :: band,ig,igp,gw_spad,cg_spad,ispinor,icg,igw,istwf_k,cg_bpad,npw_k
 integer :: comm,ierr,mb,my_nband
 real(dp) :: ecut
 logical :: found
!arrays
 integer :: gcur(3)
 integer :: my_band_list(Wfd%mband)
 integer,allocatable :: k2g(:)
 real(dp) :: kpt(3)

! *************************************************************************

 MSG_ERROR("Code has to be checked!")
 ABI_CHECK(ikg==0,"ikg/=0 not coded")

 comm    = Wfd%comm
 ecut    = Wfd%ecut
 kpt     = Wfd%kibz(:,ik_ibz)
 istwf_k = Wfd%istwfk(ik_ibz)

 ! Get the list of bands owened by this proc.
 call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list)

 if (Wfd%gamma_centered) then
   ! * Calculate set of G"s for this k-point (k-centered).
   call get_kg(kpt,istwf_k,ecut,gmet,npw_k,kg_k)
   !
   ! * Table with the correspondence btw the k-centered and the Gamma-centered basis set.
   ABI_ALLOCATE(k2g,(npw_k))
   nmiss=0
   do ig=1,npw_k
     gcur(:)=kg_k(:,ig)
     igp=0 ; found=.FALSE.
     do while ((.not.found) .and. igp<Wfd%npwwfn)
       igp=igp+1; found=ALL(gcur(:)==Wfd%gvec(:,igp))
     end do
     if (found) then ! Store it if found:
       k2g(ig) = igp
     else
       k2g(ig)=Wfd%npwwfn+1
       nmiss=nmiss+1
     end if
   end do

   cg=zero
   do mb=1,my_nband
     band=my_band_list(mb)
     cg_bpad= 0 !npw_k*Wfd%nspinor*(band-Wfd%my_minb)
     do ispinor=1,Wfd%nspinor
       cg_spad=(ispinor-1)*npw_k
       gw_spad=(ispinor-1)*Wfd%npwwfn
       do ig=1,npw_k
         icg = ig+cg_spad+cg_bpad+ikg
         igw = k2g(ig)+gw_spad
         if (k2g(ig)<Wfd%npwwfn+1) then
           cg(1,icg) = DBLE (Wfd%Wave(band,ik_ibz,spin)%ug(igw))
           cg(2,icg) = AIMAG(Wfd%Wave(band,ik_ibz,spin)%ug(igw))
         else  ! not in the gamma-centered basis set, set the Fourier component to zero.
           cg(:,icg) = zero
         end if
       end do
     end do
   end do

   ABI_DEALLOCATE(k2g)
   call xsum_mpi(cg,comm,ierr)

  else ! We are already using k-centered G-spheres
    nmiss=0
    npw_k = Wfd%npwarr(ik_ibz)
    ABI_ALLOCATE(kg_k,(3,npw_k))
    kg_k = Wfd%Kdata(ik_ibz)%kg_k

    cg=zero
    do mb=1,my_nband ! Fill my slice
      band=my_band_list(mb)
      cg_bpad= (band-1)*npw_k*Wfd%nspinor
      do ig=1,npw_k*Wfd%nspinor
        icg = ig+cg_bpad
        cg(1,icg) = DBLE (Wfd%Wave(band,ik_ibz,spin)%ug(ig))
        cg(2,icg) = AIMAG(Wfd%Wave(band,ik_ibz,spin)%ug(ig))
      end do
    end do
   call xsum_mpi(cg,comm,ierr) ! Collect on each node.
  end if

end subroutine wfd_gather_g2k
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_distribute_bbp
!! NAME
!!  wfd_distribute_bbp
!!
!! FUNCTION
!!  This routines distributes as set of (b,b') indeces taking into account the MPI distribution of the ug.
!!  It is used to calculate matrix elements of the form <b,k,s|O|b',k,s>
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=
!!  ik_ibz=The index of the k-point in the IBZ.
!!  spin=Spin index.
!!  allup=String used to select or not the upper triangle. Possible values:
!!    "All"  =Entire (b,b') matrix will be distributed.
!!    "Upper"=Only the upper triangle is distributed.
!!  [got(%nproc)]=The number of tasks already assigned to the nodes. Used to optimize the work load.
!!    Be careful when this routine is called inside several loops since each node should call the routine
!!    at each iteration with the same (local) copy of got so that bbp_distrb will assume the same value on each node.
!!  [bbp_mask(%mband,%mband)]= mask used to select a subset of (b,b') indeces.
!!
!! OUTPUT
!!  my_nbbp=The number of (b,b') indeces treated by this node.
!!  bbp_distrb(%mband%mband)=The rank of the node that will treat (b,b').
!!
!! PARENTS
!!      calc_optical_mels,calc_vhxc_me,cchi0q0
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_distribute_bbp(Wfd,ik_ibz,spin,allup,my_nbbp,bbp_distrb,got,bbp_mask)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_distribute_bbp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin
 integer,intent(out) :: my_nbbp
 type(wfs_descriptor),intent(in) :: Wfd
 character(len=*),intent(in) :: allup
!arrays
 integer,intent(out) :: bbp_distrb(Wfd%mband,Wfd%mband)
 integer,optional,intent(inout) :: got(Wfd%nproc)
 logical,optional,intent(in) :: bbp_mask(Wfd%mband,Wfd%mband)

!Local variables ------------------------------
!arrays
 integer :: loc_got(Wfd%nproc)

!************************************************************************

 ! Just a wrapper around wfd_distribute_kb_kpbp.
 loc_got=0; if (PRESENT(got)) loc_got = got

 if (PRESENT(bbp_mask)) then
   call wfd_distribute_kb_kpbp(Wfd,ik_ibz,ik_ibz,spin,allup,my_nbbp,bbp_distrb,loc_got,bbp_mask)
 else
   call wfd_distribute_kb_kpbp(Wfd,ik_ibz,ik_ibz,spin,allup,my_nbbp,bbp_distrb,loc_got)
 end if

end subroutine wfd_distribute_bbp
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_distribute_kb_kpbp
!! NAME
!!  wfd_distribute_kb_kpbp
!!
!! FUNCTION
!!  This routines distributes as set of (b,b') indeces taking into account the MPI distribution of the ug.
!!  It is used to calculate matrix elements of the form <b,k,s|O|b',k',s>
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=
!!  ik_ibz =The index of the k-point k  in the IBZ.
!!  ikp_ibz=The index of the k-point k' in the IBZ.
!!  spin=Spin index.
!!  allup=String used to select the upper triangle of the (b,b') matrix. Possible values:
!!    "All"  =Entire (b,b') matrix will be distributed.
!!    "Upper"=Only the upper triangle is distributed.
!!  [got(%nproc)]=The number of tasks already assigned to the nodes. Used to optimize the distribution of the tasks.
!!    Be careful when this routine is called inside several loops since each node should call the routine
!!    at each iteration with the same (local) copy of got so that bbp_distrb will assume the same value on each node.
!!  [bbp_mask(%mband,%mband)]= mask used to select a subset of (b,b') indeces.
!!
!! OUTPUT
!!  my_nbbp=The number of (b,b') indeces treated by this node.
!!  bbp_distrb(%mband%mband)=The rank of the node that will treat (b,b').
!!
!! PARENTS
!!      cchi0,check_completeness,m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_distribute_kb_kpbp(Wfd,ik_ibz,ikp_ibz,spin,allup,my_nbbp,bbp_distrb,got,bbp_mask)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_distribute_kb_kpbp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,ikp_ibz,spin
 integer,intent(out) :: my_nbbp
 type(wfs_descriptor),intent(in) :: Wfd
 character(len=*),intent(in) :: allup
!arrays
 integer,intent(out) :: bbp_distrb(Wfd%mband,Wfd%mband)
 integer,optional,intent(inout) :: got(Wfd%nproc)
 logical,optional,intent(in) :: bbp_mask(Wfd%mband,Wfd%mband)

!Local variables ------------------------------
!scalars
 integer :: my_nband,ib1,ib2,pcb2,pcb1,howmany_b,howmany_bp,workload_min
 integer :: rank,ncpus,idle,b1_stop,istat
 character(len=500) :: msg
!arrays
 integer :: rank_bandlist_k(Wfd%mband),rank_bandlist_kp(Wfd%mband)
 integer :: get_more(Wfd%nproc),my_band_list_k(Wfd%mband)
 integer,allocatable :: whocan_k(:,:),whocan_kp(:,:)
 logical :: b_mask(Wfd%mband)

!************************************************************************

 ABI_ALLOCATE(whocan_k ,(Wfd%mband,Wfd%nproc))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in whocan_k")
 ABI_ALLOCATE(whocan_kp,(Wfd%mband,Wfd%nproc))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in whocan_kp")
 whocan_k =0 !  Will be set to 1 if this node can calculate something containing (k,b)
 whocan_kp=0 !  Will be set to 1 if this node can calculate something containing (kp,bp)

 do rank=0,Wfd%nproc-1

   call wfd_bands_of_rank(Wfd,rank,ik_ibz ,spin,howmany_b, rank_bandlist_k )
   do pcb1=1,howmany_b
     ib1 = rank_bandlist_k(pcb1)
     whocan_k(ib1,rank+1) = 1
   end do

   call wfd_bands_of_rank(Wfd,rank,ikp_ibz,spin,howmany_bp,rank_bandlist_kp)
   do pcb2=1,howmany_bp
     ib2 = rank_bandlist_kp(pcb2)
     whocan_kp(ib2,rank+1) = 1
   end do

 end do

 get_more=0; if (PRESENT(got)) get_more=got
 b1_stop=Wfd%nband(ik_ibz,spin)

 bbp_distrb = xmpi_undefined_rank

 do ib2=1,Wfd%nband(ikp_ibz,spin)
   b_mask = .TRUE.; if (PRESENT(bbp_mask)) b_mask = bbp_mask(:,ib2)
   if (ANY(b_mask)) then
     my_nband=0; my_band_list_k=0
     if (starts_with(allup,(/"U","u"/))) b1_stop = MIN(ib2,Wfd%nband(ik_ibz,spin)) ! Only the upper triangle of the (b1,b2) matrix.

     do ib1=1,b1_stop
       if (b_mask(ib1)) then

         ! 
         ! find which CPUs can do the calculation (k,b)->(kp,bp) 
         ! find the one which is less busy
         ncpus=0
         workload_min=HUGE(0)
         do rank=0,Wfd%nproc-1
           if( whocan_k(ib1,rank+1)==1 .AND.  whocan_kp(ib2,rank+1)==1 ) then
             ncpus=ncpus+1
             if( get_more(rank+1) < workload_min ) then 
               idle=rank+1
               workload_min=get_more(idle)
             endif
             
           endif
         enddo

         if(ncpus>0) then
           bbp_distrb(ib1,ib2)=idle-1
           get_more(idle) = get_more(idle) + 1

         else
           call wfd_dump_errinfo(Wfd)
           write(msg,'(a,5(i0,1x))')" Nobody has (band1, ik_ibz) (band2, ikp_ibz) spin: ",ib1,ik_ibz,ib2,ikp_ibz,spin
           MSG_ERROR(msg)
         end if

       end if
     end do ! ib1
   end if
 end do ! ib2

 ABI_DEALLOCATE(whocan_k)
 ABI_DEALLOCATE(whocan_kp)

 my_nbbp = COUNT(bbp_distrb==Wfd%my_rank)
 if (PRESENT(got)) got=get_more

end subroutine wfd_distribute_kb_kpbp
!!***

!----------------------------------------------------------------------


!!****f* m_wfs/wfd_get_cprj
!! NAME
!!  wfd_get_cprj
!!
!! FUNCTION
!!  Return a copy of Cprj either by calculating it on-the-fly or by just retrieving the data already stored in the data type.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=the wavefunction descriptor.
!!  band=Band index.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!  sorted=.TRUE. if the output cprj matrix elements have to be sorted by atom type.
!!
!! OUTPUT
!!  Cprj_out(Wfd%natom,Wfd%nspinor) <type(cprj_type)>=Unsorted matrix elements.
!!
!! PARENTS
!!      calc_optical_mels,calc_sigc_me,calc_sigx_me,calc_vhxc_me,cchi0,cchi0q0
!!      cchi0q0_intraband,check_completeness,cohsex_me,debug_tools
!!      exc_build_block,exc_build_ham,m_shirley,m_wfs,sigma,wfd_pawrhoij
!!      wfd_vnlpsi
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_get_cprj(Wfd,band,ik_ibz,spin,Cryst,Cprj_out,sorted)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_get_cprj'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 logical,intent(in) :: sorted
 type(wfs_descriptor),intent(inout) :: Wfd
 type(crystal_structure),intent(in) :: Cryst
!arrays
 type(cprj_type),intent(out) :: Cprj_out(Wfd%natom,Wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer,parameter :: choice1=1,idir0=0
 integer :: want_order,iatom,sidx
 character(len=500) :: msg

!************************************************************************

 want_order=CPR_RANDOM; if (sorted) want_order=CPR_SORTED

 SELECT CASE (Wfd%Wave(band,ik_ibz,spin)%has_cprj)

 CASE (WFD_NOWAVE, WFD_ALLOCATED)  ! Have to calculate it!

   if (.not.wfd_ihave_ug(Wfd,band,ik_ibz,spin,"Stored")) then
     write(msg,'(a,3(i0,1x),a)')" ug for (band, ik_ibz, spin): ",band,ik_ibz,spin," is not stored in memory!"
     MSG_PERS_ERROR(msg)
   end if
   !
   ! Get cprj.
   call wfd_ug2cprj(Wfd,band,ik_ibz,spin,choice1,idir0,Wfd%natom,Cryst,Cprj_out,sorted=sorted)

   if (Wfd%Wave(band,ik_ibz,spin)%has_cprj==WFD_ALLOCATED) then ! Store it.

     if ( want_order == Wfd%Wave(band,ik_ibz,spin)%cprj_order) then
       call cprj_copy(Cprj_out,Wfd%Wave(band,ik_ibz,spin)%Cprj)
       Wfd%Wave(band,ik_ibz,spin)%has_cprj=WFD_STORED

     else ! Have to reorder cprj_out
       select case (want_order)

       case (CPR_SORTED)
         do iatom=1,Cryst%natom
           sidx = Cryst%atindx(iatom) ! random --> sorted table.
           call cprj_copy(Cprj_out(sidx:sidx,:),Wfd%Wave(band,ik_ibz,spin)%Cprj(iatom:iatom,:))
         end do

       case (CPR_RANDOM)
         do sidx=1,Cryst%natom
           iatom = Cryst%atindx1(sidx) ! sorted --> random table.
           call cprj_copy(Cprj_out(iatom:iatom,:),Wfd%Wave(band,ik_ibz,spin)%Cprj(sidx:sidx,:))
         end do

       case default
         write(msg,'(a,i0)')" Wrong value for want_order ",want_order
         MSG_PERS_ERROR(msg)
       end select

     end if
   end if

 CASE (WFD_STORED) ! copy it back.

   if (want_order == Wfd%Wave(band,ik_ibz,spin)%cprj_order) then
     call cprj_copy(Wfd%Wave(band,ik_ibz,spin)%Cprj,Cprj_out)

   else
     select case (want_order)

     case (CPR_SORTED)
       do iatom=1,Cryst%natom
         sidx = Cryst%atindx(iatom) ! random --> sorted table.
         call cprj_copy(Wfd%Wave(band,ik_ibz,spin)%Cprj(iatom:iatom,:),Cprj_out(sidx:sidx,:))
       end do

     case (CPR_RANDOM)
       do sidx=1,Cryst%natom
         iatom = Cryst%atindx1(sidx) ! sorted --> random table.
         call cprj_copy(Wfd%Wave(band,ik_ibz,spin)%Cprj(sidx:sidx,:),Cprj_out(iatom:iatom,:))
       end do

     case default
       write(msg,'(a,i0)')" Wrong value for want_order ",want_order
       MSG_PERS_ERROR(msg)
     end select

   end if

 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong has_cprj: ",Wfd%Wave(band,ik_ibz,spin)%has_cprj
   MSG_PERS_BUG(msg)
 END SELECT

end subroutine wfd_get_cprj
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_change_ngfft
!! NAME
!!  wfd_change_ngfft
!!
!! FUNCTION
!!   Reallocate and reinitialize internal tables for performing FFTs of wavefunctions.
!!
!! INPUTS
!!  Cryst<crystal_structure>=Info on unit cell.
!!  Psps<pseudopotential_type>=Pseudopotential info.
!!  new_ngfft(18)=FFT descriptor for the new FFT mesh.
!!
!!  SIDE EFFECTS
!!  Wfd<wfs_descriptor>=Wavefunction descriptor with new internal tables for FFT defined by new_ngfft.
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,calc_vhxc_me,cchi0,cchi0q0,cchi0q0_intraband
!!      check_completeness,classify_bands,cohsex_me,exc_build_block
!!      exc_build_ham,exc_plot,m_shirley,m_wfs,screening,sigma,wfd_mkrho
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_change_ngfft(Wfd,Cryst,Psps,new_ngfft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_change_ngfft'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: new_ngfft(18)
 type(crystal_structure),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays

!Local variables ------------------------------
!scalars
 integer,parameter :: npw0=0
 integer :: npw_k,ik_ibz,istwf_k
 logical :: iscompatibleFFT
 character(len=500) :: msg
!arrays
 integer,allocatable :: kg_k(:,:)

!************************************************************************

 !@wfs_descriptor
 if ( ALL(Wfd%ngfft(1:3) == new_ngfft(1:3)) ) RETURN ! Nothing to do.

 MSG_COMMENT("Changing FFT mesh")
 !
 ! Change FFT dimensions.
 Wfd%ngfft  = new_ngfft
 Wfd%mgfft  = MAXVAL(new_ngfft(1:3))
 Wfd%nfftot = PRODUCT(new_ngfft(1:3))
 Wfd%nfft   = Wfd%nfftot ! No FFT parallelism.

 if (associated(Wfd%ph1d))  then
   ABI_DEALLOCATE(Wfd%ph1d)
 end if
 ABI_ALLOCATE(Wfd%ph1d,(2,3*(2*Wfd%mgfft+1)*Cryst%natom))
 call getph(Cryst%atindx,Cryst%natom,Wfd%ngfft(1),Wfd%ngfft(2),Wfd%ngfft(3),Wfd%ph1d,Cryst%xred)
 !
 ! Recalculate FFT tables.
 ! Calculate the FFT index of $ R^{-1} (r-\tau) $ used to symmetrize u_Rk.
 if (associated(Wfd%irottb))  then
   ABI_DEALLOCATE(Wfd%irottb)
 end if
 ABI_ALLOCATE(Wfd%irottb,(Wfd%nfftot,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,Wfd%ngfft,Wfd%irottb,iscompatibleFFT)

 if (.not.iscompatibleFFT) then
   msg = " Real space FFT mesh not compatible with symmetries. Wavefunction symmetrization should not be done in real space!"
   MSG_WARNING(msg)
 end if
 !
 ! Is the new real space FFT mesh compatible with the rotational part?
 Wfd%rfft_is_symok = check_rot_fft(Cryst%nsym,Cryst%symrel,Wfd%ngfft(1),Wfd%ngfft(2),Wfd%ngfft(3))
 !
 ! Reallocate ur buffers with correct dimensions.
 call destroy_wave_3D(Wfd%Wave,"R")

 !do spin=1,Wfd%nsppol
 !  do ik_ibz=1,Wfd%nkibz
 !    do band=1,Wfd%nband(ik_ibz,spin)
 !      keep = keep_ur(band,ik_ibz,spin)
 !      if (wfd_ihave_ug(Wfd,band,ik_ibz,spin,"Stored") .and. keep) then
 !        call init_wave_0D(Wfd%Wave(band,ik_ibz,spin),Wfd%usepaw,npw0,Wfd%nfft,Wfd%nspinor,Wfd%natom,Wfd%nlmn_atm,CPR_RANDOM)
 !        Wfd%keep_ur(band,ik_ibz,spin) = .TRUE.
 !      else
 !        Wfd%keep_ur(band,ik_ibz,spin) = .FALSE.
 !      end if
 !    end do
 !  end do
 !end do
 !
 ! Reinit Kdata_t
 do ik_ibz=1,Wfd%nkibz
   if (wfd_ihave_ug(Wfd,0,ik_ibz,0)) then
     istwf_k = Wfd%istwfk(ik_ibz)
     npw_k   = Wfd%Kdata(ik_ibz)%npw
     ABI_ALLOCATE(kg_k,(3,npw_k))
     kg_k = Wfd%Kdata(ik_ibz)%kg_k
     call kdata_free(Wfd%Kdata(ik_ibz))
     call kdata_init(Wfd%Kdata(ik_ibz),Cryst,Psps,Wfd%kibz(:,ik_ibz),istwf_k,new_ngfft,Wfd%MPI_enreg,kg_k=kg_k)
     ABI_DEALLOCATE(kg_k)
   end if
 end do

end subroutine wfd_change_ngfft
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_iam_master
!! NAME
!!  wfd_iam_master
!!
!! FUNCTION
!!  Returns true if this rank is the master node. spin index can be specified.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wfd_iam_master(Wfd,bks_ids) result(ans)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_iam_master'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: bks_ids(3)
 type(wfs_descriptor),intent(in) :: Wfd
 logical :: ans

!Local variables ------------------------------
!scalars
! character(len=500) :: msg
!arrays

!************************************************************************

 if (.not.PRESENT(bks_ids)) then
   ans = (Wfd%my_rank == Wfd%master)
 else
   !FIXME
   MSG_WARNING(" spin optional argument not coded, have to introduce MPI communicators")
 end if

end function wfd_iam_master
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_test_ortho
!! NAME
!! wfd_test_ortho
!!
!! FUNCTION
!!  Test the orthonormalization of the wavefunctions stored in Wfd.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=wavefunction descriptor.
!!  Cryst<crystal_structure>=Object defining the unit cell and its symmetries.
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=PAW tabulated starting data.
!!
!! OUTPUT
!!   Only writing.
!!
!! PARENTS
!!      bethe_salpeter,m_shirley,screening,sigma
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_test_ortho(Wfd,Cryst,Pawtab,unit,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_test_ortho'
 use interfaces_14_hidewrite
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: unit
 character(len=4),optional,intent(in) :: mode_paral
 type(crystal_structure),intent(in) :: Cryst
 type(wfs_descriptor),intent(inout) :: Wfd
!array
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,band,band1,band2,ib,ib1,ib2,ierr,how_manyb,my_unt,npw_k,istwf_k
 real(dp) :: glob_cinf,my_cinf,glob_csup,my_csup,glob_einf,my_einf,glob_esup,my_esup
 complex(dpc) :: cdum
 logical :: bands_are_spread
 character(len=4) :: my_mode
 character(len=500) :: msg
!arrays
 integer :: my_bandlist(Wfd%mband)
 real(dp) :: pawovlp(2)
 complex(gwpc),pointer :: ug1(:),ug2(:)
 !complex(gwpc) :: ur(Wfd%nfft*Wfd%nspinor)
 character(len=6) :: tag_spin(2)
 type(Cprj_type),allocatable :: Cp1(:,:),Cp2(:,:)

!************************************************************************

 tag_spin(:)=(/'      ','      '/); if (Wfd%nsppol==2) tag_spin(:)=(/' UP   ',' DOWN '/)

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral
 !
 ! Update the kbs table storing the distribution of the ug.
 !call wfd_update_bkstab(Wfd)

 if (Wfd%usepaw==1) then
   ABI_ALLOCATE(Cp1,(Wfd%natom,Wfd%nspinor))
   call cprj_alloc(Cp1,0,Wfd%nlmn_atm)
   ABI_ALLOCATE(Cp2,(Wfd%natom,Wfd%nspinor))
   call cprj_alloc(Cp2,0,Wfd%nlmn_atm)
 end if

 bands_are_spread = .FALSE.

 do spin=1,Wfd%nsppol
   my_einf=greatest_real; my_esup=zero
   my_cinf=greatest_real; my_csup=zero
   !
   do ik_ibz=1,Wfd%nkibz
     istwf_k = Wfd%istwfk(ik_ibz)
     npw_k   = Wfd%npwarr(ik_ibz)
     !
     ! Select my band indeces.
     call wfd_mybands(Wfd,ik_ibz,spin,how_manyb,my_bandlist,"Stored")
     if (how_manyb/=Wfd%nband(ik_ibz,spin)) bands_are_spread = .TRUE.
     !
     ! 1) Normalization.
     do ib=1,how_manyb
       band = my_bandlist(ib)
       ug1 => Wfd%Wave(band,ik_ibz,spin)%ug
       cdum = xdotc(npw_k*Wfd%nspinor,ug1,1,ug1,1)
       if (istwf_k>1) then
         cdum=two*DBLE(cdum)
         if (istwf_k==2) cdum=cdum-CONJG(ug1(1))*ug1(1)
       end if
       if (Wfd%usepaw==1) then
         call wfd_get_cprj(Wfd,band,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
         pawovlp = paw_overlap(Cp1,Cp1,Cryst%typat,Pawtab,spinor_comm=Wfd%MPI_enreg%comm_spin)
         cdum = cdum + CMPLX(pawovlp(1),pawovlp(2))
       end if
       if (REAL(cdum)<my_einf) my_einf=REAL(cdum)
       if (REAL(cdum)>my_esup) my_esup=REAL(cdum)
     end do

     call xmin_mpi(my_einf,glob_einf,Wfd%comm,ierr) ! TODO should use the communicator for this spin
     call xmax_mpi(my_esup,glob_esup,Wfd%comm,ierr)
     !
     ! 2) Orthogonality of wavefunctions.
     do ib1=1,how_manyb
       band1 = my_bandlist(ib1)
       ug1 => Wfd%Wave(band1,ik_ibz,spin)%ug
       if (Wfd%usepaw==1) call wfd_get_cprj(Wfd,band1,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)

       do ib2=ib1+1,how_manyb
         band2 = my_bandlist(ib2)
         ug2 => Wfd%Wave(band2,ik_ibz,spin)%ug
         if (Wfd%usepaw==1) call wfd_get_cprj(Wfd,band2,ik_ibz,spin,Cryst,Cp2,sorted=.FALSE.)
         cdum = xdotc(npw_k*Wfd%nspinor,ug1,1,ug2,1)
         if (istwf_k>1) then
           cdum=two*DBLE(cdum)
           if (istwf_k==2) cdum=cdum-CONJG(ug1(1))*ug2(1)
         end if
         if (Wfd%usepaw==1) then
           pawovlp = paw_overlap(Cp1,Cp2,Cryst%typat,Pawtab,spinor_comm=Wfd%MPI_enreg%comm_spin)
           cdum = cdum + CMPLX(pawovlp(1),pawovlp(2))
         end if

         if (ABS(cdum)<my_cinf) my_cinf=ABS(cdum)
         if (ABS(cdum)>my_csup) my_csup=ABS(cdum)
         !if (ABS(cdum) > 0.1) write(std_out,*)" ib1,ib2,ABS_dotprod: ",ib1,ib2,ABS(cdum)
       end do !ib2
     end do !ib

     call xmin_mpi(my_cinf,glob_cinf,Wfd%comm,ierr) ! TODO should use the communicator for this spin
     call xmax_mpi(my_csup,glob_csup,Wfd%comm,ierr)
   end do ! ik_ibz
   !
   ! === Output results for this spin ===
   write(msg,'(2a)')ch10,' test on the normalization of the wavefunctions'
   if (Wfd%nsppol==2) write(msg,'(3a)')ch10,' test on the normalization of the wavefunctions with spin ',tag_spin(spin)
   call wrtout(my_unt,msg,mode_paral)
   write(msg,'(a,f9.6,a,a,f9.6)')&
&    ' min sum_G |a(n,k,G)| = ',glob_einf,ch10,&
&    ' max sum_G |a(n,k,G)| = ',glob_esup
   call wrtout(my_unt,msg,mode_paral)

   write(msg,'(a)')' test on the orthogonalization of the wavefunctions'
   if (Wfd%nsppol==2) write(msg,'(2a)')' test on the orthogonalization of the wavefunctions with spin ',tag_spin(spin)
   call wrtout(my_unt,msg,mode_paral)
   write(msg,'(a,f9.6,a,a,f9.6,a)')&
&    ' min sum_G a(n,k,G)* a(n",k,G) = ',glob_cinf,ch10,&
&    ' max sum_G a(n,k,G)* a(n",k,G) = ',glob_csup,ch10
   call wrtout(my_unt,msg,mode_paral)

 end do ! spin

 if (bands_are_spread) then
   write(msg,'(6a)')&
&    ' rdkss : COMMENT -',ch10,&
&    '  Note that the test on the orthogonalization is not complete ',ch10,&
&    '  since bands are spread among different processors',ch10
   call wrtout(my_unt,msg,mode_paral)
 end if

 if (Wfd%usepaw==1) then
   call cprj_free(Cp1)
   ABI_DEALLOCATE(Cp1)
   call cprj_free(Cp2)
   ABI_DEALLOCATE(Cp2)
 end if

 call flush_unit(my_unt)

end subroutine wfd_test_ortho
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_barrier
!! NAME
!!  wfd_barrier
!!
!! FUNCTION
!!
!! INPUTS
!!  Wfd<wfs_descriptor>
!!
!! PARENTS
!!      calc_sigc_me,cchi0,cchi0q0,cchi0q0_intraband,check_completeness
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_barrier(Wfd,bks_ids)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_barrier'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: bks_ids(3)
 type(wfs_descriptor),intent(in) :: Wfd

!************************************************************************

 if (.not.PRESENT(bks_ids)) then ! Synch all nodes in Wfd%comm.
   call xbarrier_mpi(Wfd%comm)
 else
   MSG_WARNING(" spin optional argument not coded, have to introduce MPI communicators")
   call xbarrier_mpi(Wfd%bks_comm(bks_ids(1),bks_ids(2),bks_ids(3)))
 end if

end subroutine wfd_barrier
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_sym_ur
!! NAME
!!  wfd_sym_ur
!!
!! FUNCTION
!!  Symmetrize a wave function in real space
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=the wavefunction descriptor.
!!  Cryst<crystal_structure>=Structure describing the crystal structure and its symmetries.
!!  Kmesh<bz_mesh_type>=Structure describing the BZ sampling
!!  band=Band index.
!!  ik_bz=Index of the k-point in the BZ.
!!  spin=Spin index
!!
!! OUTPUT
!!  ur_kbz(Wfd%nfft*Wfd%nspinor)=The symmetrized wavefunction in real space.
!!
!! PARENTS
!!      debug_tools,exc_plot,m_shirley
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_sym_ur(Wfd,Cryst,Kmesh,band,ik_bz,spin,ur_kbz)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_sym_ur'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_bz,spin
 type(crystal_structure),intent(in) :: Cryst
 type(bz_mesh_type),intent(in) :: Kmesh
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 complex(gwpc),intent(out) :: ur_kbz(Wfd%nfft*Wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,isym_k,itim_k,nr
 integer :: ispinor,spad,ir,ir2
 complex(dpc) :: ph_mkt,u2b,u2a
 logical :: isirred
 character(len=500) :: msg
!arrays
 integer :: umklp(3)
 integer,pointer :: tabr_k(:)
 real(dp) :: kbz(3)
 real(dp),pointer :: spinrot_k(:)
 complex(dpc) :: spinrot_mat(2,2)
 complex(dpc) :: eig0r(Wfd%nfft)
 complex(gwpc) :: ur_kibz(Wfd%nfft*Wfd%nspinor)

!************************************************************************

 ! k_bz =  S k - G0 ==> u_{k_bz} =  e^{iG0.r} u_{Sk}
 ! k_bz = -S k - G0 ==> u_{k_bz} =  e^{iG0.r} u_{Sk}^*

 ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz)
 !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
 !
 ! * Get ik_ibz, non-symmorphic phase, ph_mkt, and symmetries from ik_bz.
 call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k,ph_mkt,umklp,isirred)

 if (isirred) then ! Avoid symmetrization if this point is irreducible.
   call wfd_get_ur(Wfd,band,ik_ibz,spin,ur_kbz); RETURN
 end if
 !
 ! Reconstruct ur in the BZ from the corresponding wavefunction in IBZ.
 call wfd_get_ur(Wfd,band,ik_ibz,spin,ur_kibz)

 if (ANY(umklp/=0)) eig0r = ceigr(umklp,Wfd%nfft,Wfd%ngfft)

 tabr_k  => Wfd%irottb(:,isym_k) ! Table for rotated FFT points

 SELECT CASE (Wfd%nspinor)
 CASE (1)
   ur_kbz = ur_kibz(tabr_k)*ph_mkt
   if (itim_k==2) ur_kbz = CONJG(ur_kbz)
   if (ANY(umklp/=0)) ur_kbz = ur_kbz*eig0r

 CASE (2)
   MSG_ERROR("Implementation has to be tested")

   nr = Wfd%nfft
   spinrot_k => Cryst%spinrot(:,isym_k)
   !
   ! ==== Apply Time-reversal if required ====
   ! \psi_{-k}^1 =  (\psi_k^2)^*
   ! \psi_{-k}^2 = -(\psi_k^1)^*
   if (itim_k==1) then
     ur_kbz = ur_kibz
   else if (itim_k==2) then
     ur_kbz(1:nr)     = CONJG(ur_kibz(nr+1:2*nr))
     ur_kbz(nr+1:2*nr)=-CONJG(ur_kibz(1:nr))
   else
     MSG_ERROR('Wrong i2 in spinor')
   end if
   !
   ! Rotate wavefunctions in real space.
   do ispinor=1,Wfd%nspinor
     spad=(ispinor-1)*nr
     do ir=1,nr
       ir2=tabr_k(ir)
       ur_kbz(ir+spad) = ur_kbz(ir2+spad)*ph_mkt
     end do
   end do
   !
   ! Rotation in spinor space.
   spinrot_mat(1,1)= spinrot_k(1) + j_dpc*spinrot_k(4)
   spinrot_mat(1,2)= spinrot_k(3) + j_dpc*spinrot_k(2)
   spinrot_mat(2,1)=-spinrot_k(3) + j_dpc*spinrot_k(2)
   spinrot_mat(2,2)= spinrot_k(1) - j_dpc*spinrot_k(4)

   do ir=1,nr
     u2a=ur_kbz(ir)
     u2b=ur_kbz(ir+nr)
     ur_kbz(ir)   =spinrot_mat(1,1)*u2a+spinrot_mat(1,2)*u2b
     ur_kbz(ir+nr)=spinrot_mat(2,1)*u2a+spinrot_mat(2,2)*u2b
   end do

   if (ANY(umklp /=0)) then
     ur_kbz(1:Wfd%nfft) = ur_kbz(1:Wfd%nfft)  *eig0r
     ur_kbz(Wfd%nfft+1:) = ur_kbz(Wfd%nfft+1:)*eig0r
   end if

 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong value for nspinor: ",Wfd%nspinor
   MSG_ERROR(msg)
 END SELECT

end subroutine wfd_sym_ur
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_read_wfk
!! NAME
!! wfd_read_wfk
!!
!! FUNCTION
!!  This routine read a standard k-centered WFK file completing the initialization of the wavefunction
!!  descriptor used in the GW code.
!!
!! INPUTS
!!  wfd_fname=Name of the WFK file.
!!  accesswff=Option specifying the fileformat as well as the IO mode to be used.
!!
!! SIDE EFFECTS
!!  Wfd<wfs_descriptor>=All the states owned by this node whose status is (STORED|ALLOCATED) read.
!!
!! PARENTS
!!      bethe_salpeter,screening,sigma
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_read_wfk(Wfd,wfk_fname,accesswff)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_read_wfk'
 use interfaces_14_hidewrite
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff
 character(len=*),intent(in) :: wfk_fname
 type(wfs_descriptor),intent(inout) :: Wfd

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_rwwf0=0,headform0=0,icg0=0,formeig0=0,optkg1=1
 integer :: wfk_unt,npw_k,nmiss,ig,igp,ierr
 integer :: comm,master,my_rank,rdwr,fform
 integer :: spin,ik_ibz,option
 integer :: mcg,nband_k,nband_disk,mband,band
 integer :: spinor,cg_spad,gw_spad,icg,igw,cg_bpad,ib
 character(len=500) :: msg
 logical :: found
 type(Wffile_type) :: Wff
 type(Hdr_type) :: Hdr
!arrays
 integer :: gcur(3)
 integer,allocatable :: k2g(:),kg_k(:,:)
 real(dp),allocatable :: eig_k(:),occ_k(:),cg_k(:,:)
 logical,allocatable :: my_readmask(:,:,:)
 character(len=6) :: tag_spin(2)

!************************************************************************

 DBG_ENTER("COLL")

 if (ANY(accesswff == (/IO_MODE_NETCDF, IO_MODE_FORTRAN_MASTER/) )) then
   write(msg,'(a,i0)')" Unsupported value for accesswff ",accesswff
   MSG_ERROR(msg)
 end if

 comm    = Wfd%comm
 master  = Wfd%master
 my_rank = Wfd%my_rank

 wfk_unt = get_unit()
 call wrtout(std_out," wfd_read_wfk : about to read "//TRIM(wfk_fname),"COLL")

 tag_spin(:)=(/'      ','      '/); if (Wfd%nsppol==2) tag_spin(:)=(/' UP   ',' DOWN '/)
 !
 ! * Init Wff structure.
 call WffOpen(accesswff,comm,wfk_fname,ierr,Wff,master,my_rank,wfk_unt)
 !
 ! * Read Header.
 rdwr=1; fform=2
 if (ANY(Wff%accesswff == (/IO_MODE_FORTRAN, IO_MODE_FORTRAN_MASTER, IO_MODE_MPI/) )) then
   call hdr_io(fform,Hdr,rdwr,Wff)
   call WffKg(wff,optkg1)  ! 1 for reading the G-vectors.
 else if (Wff%accesswff==IO_MODE_ETSF) then
   call hdr_io_etsf(fform,Hdr,rdwr,Wff%unwff)
 end if
 !
 ! * Output the header of the GS wavefunction file.
 if (Wfd%prtvol>0) call hdr_io(fform,Hdr,4,std_out)
 !
 ! TODO: Perform consistency check btw Hdr and Wfd.
 !
 ! For each spin and k-point, do:
 !  1) Convert Wfd from gamma-centered to k-centered basis set
 !  2) Write G vectors, energies, occ and u(G) on file.
 !
 ! Each node will read the waves whose status if (WFD_ALLOCATED|WFD_STORED).
 ABI_ALLOCATE(my_readmask,(Wfd%mband,Wfd%nkibz,Wfd%nsppol))
 my_readmask=.FALSE.
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do band=1,Wfd%nband(ik_ibz,spin)

       if (wfd_ihave_ug(Wfd,band,ik_ibz,spin)) then
         my_readmask(band,ik_ibz,spin) = .TRUE.
         if (wfd_ihave_ug(Wfd,band,ik_ibz,spin,"Stored")) then
           MSG_WARNING("Wavefunction is already stored!")
         end if
       end if

     end do
   end do
 end do

 write(msg,'(a,i0,a)')" wfd_read_wfk: will read ",COUNT(my_readmask)," (b,k,s) states"
 call wrtout(std_out,msg,"PERS")

 ! TODO, to be removed. Needed only to have the same output as the one given by wfd_rdkss.
 write(msg,'(2a)')ch10,' k       eigenvalues [eV]'
 call wrtout(ab_out,msg,'COLL')

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz

     npw_k     = Hdr%npwarr(ik_ibz)
     !nband_k   = Hdr%nband(ik_ibz+(spin-1)*Hdr%nkpt)

     nband_disk= Hdr%nband(ik_ibz+(spin-1)*Hdr%nkpt)
     nband_k   = Wfd%nband(ik_ibz,spin)
     if ( nband_k > nband_disk ) then
       write(msg,'(a,2(i0,1x))')&
&       " nband_k to be read cannot be greater than nband_disk while: ",nband_k,nband_disk
       MSG_ERROR(msg)
     end if
     mband     = nband_k
     mcg       = npw_k*Wfd%nspinor*mband

     ABI_ALLOCATE(eig_k,((2*mband)**formeig0*mband))
     ABI_ALLOCATE(occ_k,(mband))
     option=1 ! for reading cg and eigen,

     ABI_ALLOCATE(cg_k,(2,mcg))
     ABI_ALLOCATE(kg_k,(3,optkg1*npw_k))
     !
     ! Read the block of bands for this (k,s).
     call rwwf(cg_k,eig_k,formeig0,headform0,icg0,ik_ibz,spin,kg_k,mband,mcg,Wfd%MPI_enreg,nband_k,&
&       nband_disk,npw_k,Wfd%nspinor,occ_k,option,optkg1,tim_rwwf0,Wff)

     ! TODO, to be removed. Needed only to have the same output as the one given by wfd_rdkss.
     if (Wfd%my_rank==Wfd%master) then
       if (Wfd%nsppol==2) then
         write(ab_out,'(i3,a,10f7.2/50(10x,10f7.2/))') ik_ibz,tag_spin(spin),(eig_k(ib)*Ha_eV,ib=1,nband_k)
       else
         write(ab_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik_ibz,(eig_k(ib)*Ha_eV,ib=1,nband_k)
       end if
     end if

     if (Wfd%gamma_centered) then
       !
       ! * Table with the correspondence btw the k-centered and the Gamma-centered basis set.
       ABI_ALLOCATE(k2g,(npw_k))
       nmiss=0
       do ig=1,npw_k
         gcur=kg_k(:,ig)
         igp=0; found=.FALSE.
         do while (.not.found .and. igp<Wfd%npwwfn)
           igp=igp+1; found=ALL(gcur(:)==Wfd%gvec(:,igp))
         end do
         if (found) then ! Store it if found:
           k2g(ig) = igp
         else
           k2g(ig)=Wfd%npwwfn+1
           nmiss=nmiss+1
         end if
       end do

       if (nmiss/=0) then
         write(msg,'(a,2(1x,i0),a,i0)')" For (k,s) ",ik_ibz,spin," the number of missing G is ",nmiss
         MSG_WARNING(msg)
       end if
       !
       ! * Conversion of the basis set.
       do band=1,Wfd%nband(ik_ibz,spin)
         if (my_readmask(band,ik_ibz,spin)) then

           Wfd%Wave(band,ik_ibz,spin)%ug = czero
           cg_bpad=npw_k*Wfd%nspinor*(band-1)
           do spinor=1,Wfd%nspinor
             cg_spad=(spinor-1)*npw_k
             gw_spad=(spinor-1)*Wfd%npwwfn
             do ig=1,npw_k
               icg = ig+cg_spad+cg_bpad
               igw = k2g(ig)+gw_spad
               if (k2g(ig)<Wfd%npwwfn+1) then
                 Wfd%Wave(band,ik_ibz,spin)%ug(igw) = DCMPLX(cg_k(1,icg),cg_k(2,icg))
               !else  ! not in the gamma-centered basis set, the Fourier component is set to zero.
               !  Wfd%Wave(band,ik_ibz,spin)%ug(igw) = czero
               end if
             end do
           end do
           Wfd%Wave(band,ik_ibz,spin)%has_ug = WFD_STORED
         end if
       end do

     else
       !
       ! * Table with the correspondence btw the k-centered sphere of the WFK file
       !   and the one used in Wfd possibly smaller due to ecutwfn.
       ABI_ALLOCATE(k2g,(npw_k))
       nmiss=0
       do ig=1,npw_k
         gcur=kg_k(:,ig)
         igp=0; found=.FALSE.
         do while (.not.found .and. igp<Wfd%npwarr(ik_ibz))
           igp=igp+1; found=ALL(gcur(:)==Wfd%Kdata(ik_ibz)%kg_k(:,igp))
         end do
         if (found) then ! Store it if found.
           k2g(ig) = igp
         else
           k2g(ig)=Wfd%npwarr(ik_ibz)+1
           nmiss=nmiss+1
         end if
       end do

       if (nmiss/=0) then
         write(msg,'(a,2(1x,i0),a,i0)')" For (k,s) ",ik_ibz,spin," the number of missing G is ",nmiss
         MSG_WARNING(msg)
       end if
       !
       ! * Conversion of the basis set.
       do band=1,Wfd%nband(ik_ibz,spin)

         if (my_readmask(band,ik_ibz,spin)) then

           Wfd%Wave(band,ik_ibz,spin)%ug = czero
           cg_bpad=npw_k*Wfd%nspinor*(band-1)
           do spinor=1,Wfd%nspinor
             cg_spad=(spinor-1)*npw_k
             gw_spad=(spinor-1)*Wfd%npwarr(ik_ibz)
             do ig=1,npw_k
               icg = ig+cg_spad+cg_bpad
               igw = k2g(ig)+gw_spad
               if (k2g(ig)<Wfd%npwarr(ik_ibz)+1) then
                 Wfd%Wave(band,ik_ibz,spin)%ug(igw) = CMPLX(cg_k(1,icg),cg_k(2,icg))
               !else  ! not in the gamma-centered basis set, the Fourier component is set to zero.
               !  Wfd%Wave(band,ik_ibz,spin)%ug(igw) = czero
               end if
             end do
           end do
           Wfd%Wave(band,ik_ibz,spin)%has_ug = WFD_STORED
         end if

       end do
     end if ! Gamma centered

     ABI_DEALLOCATE(eig_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(cg_k)
     ABI_DEALLOCATE(k2g)

   end do !ik_ibz
 end do !spin
 !
 ! * Close the wavefunction file (and do NOT delete it !)
 call WffClose(Wff,ierr)
 !
 ! * Free local memory.
 ABI_DEALLOCATE(my_readmask)
 call hdr_clean(Hdr)

 call wfd_update_bkstab(Wfd)

 DBG_EXIT("COLL")

end subroutine wfd_read_wfk
!!***

!----------------------------------------------------------------------


!!****f* m_wfs/wfd_write_wfk
!! NAME
!! wfd_write_wfk
!!
!! FUNCTION
!!  Produces a standard k-centered WFK file containing QP energies, occupancies
!!  and wavefunction. The WFK file will contain a modified header with the
!!  basic parameters used for the GW calculation and the QP onsite occupancies
!!  in the case of PAW.
!!  **NOTE** that Bandstructure_type and Pawrhoij are supposed to have been updated outside the routine.
!!
!! INPUTS
!!  codvsn=code version
!!  etot=total energy (Hartree)
!!  residm=maximal residual
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_write_wfk(Wfd,Cryst,Psps,Pawtab,Pawrhoij,BSt,Wvl,accesswff,wfk_fname,codvsn,&
&  etot,residm,intxc,ixc,pawcpxocc,pawecutdg,stmbias,qptn,so_psp,ngfftdg)

 use defs_basis
 use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_write_wfk'
 use interfaces_14_hidewrite
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff,intxc,ixc,pawcpxocc
 real(dp),intent(in) :: etot,residm,stmbias,pawecutdg
 character(len=fnlen),intent(in) :: wfk_fname
 character(len=6),intent(in) :: codvsn
 type(wfs_descriptor),intent(inout) :: Wfd
 type(Bandstructure_type),intent(in) :: BSt
 type(Crystal_structure),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
 type(wvl_internal_type),intent(in) :: wvl
!arrays
 integer,intent(in) ::  so_psp(Psps%npsp),ngfftdg(18)
 real(dp),intent(in) :: qptn(3)
 type(Pawrhoij_type),intent(in) :: Pawrhoij(Cryst%natom*Psps%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_rwwf0=0,headform0=0,icg0=0
 integer :: wfk_unt,formeig,npw_k,ierr
 integer :: spaceComm,master,my_rank,rdwr,fform,imiss,ikg
 integer :: spin,ik_ibz,option,optkg
 integer :: bantot,ib,jj,mcg,nband_k,nband_disk,mband,pertcase
 character(len=500) :: msg
 type(Wffile_type) :: Wff
 type(Hdr_type) :: my_Hdr
!arrays
 real(dp),allocatable :: eig_k(:),occ_k(:)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:)
 integer,pointer :: kg_k(:,:)
 real(dp),allocatable  :: cg_k(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 MSG_ERROR("Recheck the implementation as Wfd has been changed")

 if (ANY(accesswff == (/IO_MODE_NETCDF, IO_MODE_MPI/) )) then
   write(msg,'(a,i0)')" Unsupported value for accesswff ",accesswff
   MSG_ERROR(msg)
 end if

 spaceComm = Wfd%comm
 master    = Wfd%master
 my_rank   = Wfd%my_rank

 call wfd_update_bkstab(Wfd)

 call wfd_set_mpicomm(Wfd)

 wfk_unt = get_unit()
 call wrtout(std_out,"wfd_write_wfk : about to write "//TRIM(wfk_fname),"COLL")

 ! === Initialize a local band structure with QP results ===
 ! * Copy QP energies and occupations

 bantot=SUM(Wfd%nband)
 ABI_ALLOCATE(doccde,(bantot))
 ABI_ALLOCATE(eigen,(bantot))
 ABI_ALLOCATE(occfact,(bantot))
 doccde=zero; eigen=zero; occfact=zero

 jj=0
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do ib=1,Wfd%nband(ik_ibz,spin)
       if (ib<=BSt%nband(ik_ibz+(spin-1)*BSt%nkpt)) then
         jj=jj+1
         occfact(jj)=BSt%occ(ib,ik_ibz,spin)
         eigen  (jj)=BSt%eig(ib,ik_ibz,spin)
       else
         MSG_ERROR("I was about to SIGFAULT")
       end if
     end do
   end do
 end do

 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen)

 ! === Init new header ===
 pertcase=0

 call hdr_init_lowlvl(my_Hdr,BSt,Psps,Pawtab,Wvl,&
&  codvsn,pertcase,Cryst%natom,Cryst%nsym,Wfd%nspden,Wfd%ecut,pawecutdg,Wfd%ecutsm,Wfd%dilatmx,&
&  intxc,ixc,stmbias,Wfd%usewvl,pawcpxocc,Wfd%ngfft,ngfftdg,so_psp,qptn,&
&  Cryst%rprimd,Cryst%xred,Cryst%symrel,Cryst%tnons,Cryst%symafm,Cryst%typat)

 ! Copy some quantities in Hdr: Fermi level, Pawrhoij for PAW
 !
 call hdr_update(bantot,etot,BSt%fermie,my_Hdr,Cryst%natom,residm,Cryst%rprimd,occfact,&
&  Wfd%MPI_enreg,Pawrhoij,Psps%usepaw,Cryst%xred)

 ABI_DEALLOCATE(occfact)
 !
 ! * Init Wff structure.
 call WffOpen(accesswff,spaceComm,wfk_fname,ierr,Wff,master,my_rank,wfk_unt)

 !call wfd_barrier(Wfd)

 ! * Write Header to unformatted file
 rdwr=2; fform=2

 if (ANY(Wff%accesswff == (/IO_MODE_FORTRAN, IO_MODE_FORTRAN_MASTER, IO_MODE_MPI/) )) then
   call hdr_io(fform,my_Hdr,rdwr,Wff)
   call WffKg(Wff,1)
 else if (Wff%accesswff==IO_MODE_ETSF .and. my_rank==Wff%master) then
   call hdr_io_etsf(fform,my_Hdr,rdwr,Wff%unwff)
 end if
 !
 ! For each spin and k-point, do:
 !  1) Convert Wfd from gamma-centered to k-centered basis set
 !  2) Write G vectors, energies, occ and u(G) on file.
 !
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz

     nband_k   = Wfd%nband(ik_ibz,spin)
     mband     = nband_k
     nband_disk= nband_k
     !
     ! * Perform conversion of the basis set for a bloch of bands taking into account the distribution of the ug.
     !   cg_k is on each node when wfd_gather_g2k returns.
     ikg=0
     call wfd_gather_g2k(Wfd,ik_ibz,spin,ikg,kg_k,cg_k,Cryst%gmet,imiss)
     if (imiss>0) then
       write(msg,'(2(a,i4))')' Missing ',imiss,' components for ik_ibz= ',ik_ibz
       MSG_ERROR(msg)
     end if

     npw_k = SIZE(kg_k,DIM=2)
     mcg   =npw_k*Wfd%nspinor*mband

     ABI_CHECK(SIZE(cg_k,DIM=2)==mcg,"Size mismatch")

     formeig=0
     ABI_ALLOCATE(eig_k,((2*mband)**formeig * mband))
     ABI_ALLOCATE(occ_k,(mband))
     occ_k = BSt%occ(1:mband,ik_ibz,spin)
     eig_k = BSt%eig(1:mband,ik_ibz,spin)
     option=2; optkg=1

     if (wfd_iam_master(Wfd)) then ! Write set of bands for this (k,s).
       call rwwf(cg_k,eig_k,formeig,headform0,icg0,ik_ibz,spin,kg_k,mband,mcg,Wfd%MPI_enreg,nband_k,&
&        nband_disk,npw_k,Wfd%nspinor,occ_k,option,optkg,tim_rwwf0,Wff)
     end if

     ABI_DEALLOCATE(eig_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(cg_k)

   end do !ik_ibz
 end do !spin
 !
 ! * Close the wavefunction file (and do NOT delete it !)
 call WffClose(Wff,ierr)
 !
 ! * Free local memory.
 call hdr_clean(my_Hdr)

 DBG_EXIT("COLL")

end subroutine wfd_write_wfk
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_paw_get_aeur
!! NAME
!! wfd_paw_get_aeur
!!
!! FUNCTION
!!
!! INPUTS
!!   band,ik_ibz,spin=indeces specifying the band, the k-point and the spin.
!!   Psps<pseudopotential_type>=variables related to pseudopotentials
!!   Cryst<Crystal_structure>= data type gathering info on symmetries and unit cell.
!!   Wfd<wfs_descriptor>=wavefunction descriptor.
!!   Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data.
!!   Pawfgrtab(natom)<pawfgrtab_type>=atomic data given on fine rectangular grid.
!!     NB: rpaw should be used in nhatgrid to initialize the datatype (optcut=1 option) instead of the radius for the
!!     shape functions (rpaw /= rshp).
!!   Paw_onsite(natom)<paw_pwaves_lmn_t>=3D PAW partial waves in real space for each FFT point in the PAW spheres.
!!
!! OUTPUT
!! ur_ae(Wfd%nfft*Wfd%nspinor)=AE PAW wavefunction in real space.
!! [ur_ae_onsite(Wfd%nfft*Wfd%nspinor)]
!! [ur_ps_onsite(Wfd%nfft*Wfd%nspinor)]
!!
!! NOTES
!!  (1) The true wavefunction integrates in real space to the unit cell volume.
!!      The definition of the cprj matrix elements includes the term 1/SQRT(ucvol) that comes
!!      from the use of a normalized planewave e^(iG.r)/SQRT(omega) in the FFT transform G-->R (see e.g. opernla_ylm)
!!      On the contrary, the convention for the G-->R transform employed in the FFT routines used in abinit is
!!      u(r) = sum_G u(G) e^(iG.r); u(G) = one/omega \int u(r) e^(-iG.r)dr.
!!      Hence we have to multiply the onsite part by SQRT(uvol) before adding the smooth FFT part in real space.
!!
!!  (2) Care has to be taken in the calculation of the onsite contribution when the FFT point belongs to the PAW
!!      sphere of a periodically repeated atom. In this case one evaluates the onsite term associated to the
!!      atom in the first unit cell then the contribution has to be multiplied by a k- dependent
!!      phase factor to account for the wrapping of the real-space point in the first unit cell.
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,check_completeness
!!      classify_bands,m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_paw_get_aeur(Wfd,band,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae,ur_ae_onsite,ur_ps_onsite)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_paw_get_aeur'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 type(pseudopotential_type),intent(in) :: Psps
 type(crystal_structure),intent(in) :: Cryst
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)
 type(pawfgrtab_type),intent(in) :: Pawfgrtab(Cryst%natom)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
 complex(gwpc),intent(out) :: ur_ae(Wfd%nfft*Wfd%nspinor)
 complex(gwpc),optional,intent(out) :: ur_ae_onsite(Wfd%nfft*Wfd%nspinor)
 complex(gwpc),optional,intent(out) :: ur_ps_onsite(Wfd%nfft*Wfd%nspinor)

!Local variables-------------------------------
!scalars
 integer :: itypat,ln_size,lmn_size,iatom,spinor
 integer :: nfgd,ifgd,jlmn,jl,jm,ifftsph
 real(dp) :: phj,tphj,arg,re_cp,im_cp
 complex(dpc) :: cp,cnorm
!arrays
 real(dp) :: kpoint(3)
 complex(dpc) :: eikr(Wfd%nfftot)
 complex(dpc),allocatable :: phk_atm(:)
 type(cprj_type),allocatable :: Cp1(:,:)

! *************************************************************************

 ! TODO ngfft should be included in pawfgrtab_type
 !% if (ANY(Wfd%ngfft(1:3)/=Pawfgrtab%ngfft(1:3)) then
 !!  MSG_ERROR("Wfd%ngfft(1:3)/=Pawfgrtab%ngfft(1:3)")
 !% end if

 call wfd_get_ur(Wfd,band,ik_ibz,spin,ur_ae)

 kpoint = Wfd%kibz(:,ik_ibz)

 eikr = ceikr(kpoint,Wfd%nfftot,Wfd%ngfft)
 ur_ae = ur_ae * eikr

 ABI_ALLOCATE(Cp1,(Wfd%natom,Wfd%nspinor))
 call cprj_alloc(Cp1,0,Wfd%nlmn_atm)

 call wfd_get_cprj(Wfd,band,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
 !
 ! === Add onsite term on the augmented FFT mesh ===
 if (PRESENT(ur_ae_onsite)) ur_ae_onsite = czero
 if (PRESENT(ur_ps_onsite)) ur_ps_onsite = czero

 ABI_CHECK(Wfd%nspinor==1,"nspinor==1 not coded")

 do iatom=1,Cryst%natom
   itypat  =Cryst%typat(iatom)
   lmn_size=Pawtab(itypat)%lmn_size
   ln_size =Pawtab(itypat)%basis_size   ! no. of nl elements in PAW basis.
   nfgd    =Pawfgrtab(iatom)%nfgd       ! no. of points in the fine grid for this PAW sphere.

   ABI_ALLOCATE(phk_atm,(nfgd))
   do ifgd=1,nfgd
     arg = -two_pi* DOT_PRODUCT(Paw_onsite(iatom)%r0shift(:,ifgd),kpoint)
     phk_atm(ifgd) = DCMPLX(COS(arg),SIN(arg))
   end do

   do spinor=1,Wfd%nspinor
     do jlmn=1,lmn_size
       jl=Psps%indlmn(1,jlmn,itypat)
       jm=Psps%indlmn(2,jlmn,itypat)
       re_cp = Cp1(iatom,spinor)%cp(1,jlmn)
       im_cp = Cp1(iatom,spinor)%cp(2,jlmn)
       cp = DCMPLX(re_cp, im_cp) * SQRT(Cryst%ucvol) ! Pay attention here. see (1).

       do ifgd=1,nfgd ! loop over fine grid points in current PAW sphere.
         ifftsph = Pawfgrtab(iatom)%ifftsph(ifgd) ! index of the point on the grid
         phj  = Paw_onsite(iatom)% phi(ifgd,jlmn)
         tphj = Paw_onsite(iatom)%tphi(ifgd,jlmn)
         ur_ae(ifftsph)           = ur_ae(ifftsph) + cp * (phj-tphj) * phk_atm(ifgd)
         if (PRESENT(ur_ae_onsite)) ur_ae_onsite(ifftsph) = ur_ae_onsite(ifftsph) + cp *  phj * phk_atm(ifgd)
         if (PRESENT(ur_ps_onsite)) ur_ps_onsite(ifftsph) = ur_ps_onsite(ifftsph) + cp * tphj * phk_atm(ifgd)
       end do
     end do !jlmn
   end do !spinor

   ABI_DEALLOCATE(phk_atm)
 end do !iatom
 !
 ! * Remove the phase e^{ikr}, u(r) is returned.
 ur_ae = ur_ae * CONJG(eikr)
 cnorm = xdotc(Wfd%nfft*Wfd%nspinor,ur_ae,1,ur_ae,1)/Wfd%nfft

 !write(std_out,*)" AE PAW norm: (b,k,s)",band,ik_ibz,spin,REAL(cnorm)

 if (PRESENT(ur_ae_onsite)) ur_ae_onsite = ur_ae_onsite * CONJG(eikr)
 if (PRESENT(ur_ps_onsite)) ur_ps_onsite = ur_ps_onsite * CONJG(eikr)

 call cprj_free(Cp1)
 ABI_DEALLOCATE(Cp1)

end subroutine wfd_paw_get_aeur
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_plot_ur
!! NAME
!! wfd_plot_ur
!!
!! FUNCTION
!!  This routine writes the squared modulus of the wavefunctions in real space
!!  on external files, one for each (k,b,s). File are written in the XSF format (Xcrysden).
!!  A subset of (b,k,s) states can be specified via the bks_mask. The routine is MPI parallelized.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=Wavefunction descriptor.
!!  Cryst<Crystal_structure>= Information on symmetries and unit cell.
!!  Psps<pseudopotential_type>=Pseudopotential info.
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=PAW tabulated starting data.
!!  Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data.
!!  ngfftf(18)=The FFT mesh used for plotting |wfr|**2, it can differ from the one internally used in Wfd.
!!    For example, PAW wavefunctions should be plotted on a much finer FFT mesh.
!!  bks_mask(mband,nkibz,nsppol)=logical mask used to select the states to be plotted.
!!
!! OUTPUT
!!  Output is written on file.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_plot_ur(Wfd,Cryst,Psps,Pawtab,Pawrad,ngfftf,bks_mask)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_plot_ur'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfftf(18)
 logical,target,intent(in) :: bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: spin,band,ik_ibz,iatom,optcut,optgr0,optgr1,optgr2,optrad
 integer :: n1,n2,n3,my_nplots,plot,funt,my_nband,cplex
 !character(len=500) :: msg
 character(len=fnlen) :: xsf_fname
!arrays
 integer :: got(Wfd%nproc)
 integer,allocatable :: l_size_atm(:),my_plot_list(:,:)
 integer :: my_band_list(Wfd%mband)
 real(dp),allocatable :: data_plot(:)
 logical,pointer :: bmask(:)
 complex(gwpc),allocatable :: ur_ae(:),nc_ur(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

!************************************************************************

 if (ALL(.not.bks_mask)) RETURN

 DBG_ENTER("COLL")

 call wrtout(std_out," Plotting |wfs|^2 ...","COLL")
 !
 ! Change the FFT mesh if needed because we want u(r) on the ngfftf mesh (pawecutd for PAW).
 call wfd_change_ngfft(Wfd,Cryst,Psps,ngfftf)
 n1 = ngfftf(1)
 n2 = ngfftf(2)
 n3 = ngfftf(3)
 !
 ! Distribute the plots among the nodes taking into account the distribution of the waves.
 ! my_plot_list gives the list of (b,k,s) states plotted by this node.
 ABI_ALLOCATE(my_plot_list,(3,Wfd%mband*Wfd%nkibz*Wfd%nsppol))

 my_nplots=0; got=0
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     bmask => bks_mask(:,ik_ibz,spin)
     call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,got,bmask)

     if (my_nband>0) then
       my_plot_list(1,my_nplots+1:my_nplots+my_nband) = my_band_list(1:my_nband)
       my_plot_list(2,my_nplots+1:my_nplots+my_nband) = ik_ibz
       my_plot_list(3,my_nplots+1:my_nplots+my_nband) = spin
       my_nplots = my_nplots + my_nband
     end if
   end do
 end do
 funt = get_unit()

 if (Wfd%usepaw==1) then
   MSG_WARNING("Testing the calculation of AE PAW wavefunctions.")
   ! Use a local pawfgrtab to make sure we use the correction in the paw spheres
   ! the usual pawfgrtab uses r_shape which may not be the same as r_paw.
   ABI_ALLOCATE(Pawfgrtab,(Cryst%natom))
   ABI_ALLOCATE(l_size_atm,(Cryst%natom))
   do iatom=1,Cryst%natom
     l_size_atm(iatom) = Pawtab(Cryst%typat(iatom))%l_size
   end do

   cplex=1
   call pawfgrtab_init(Pawfgrtab,cplex,l_size_atm,Wfd%nspden)
   ABI_DEALLOCATE(l_size_atm)

   optcut=1                     ! use rpaw to construct Pawfgrtab.
   optgr0=0; optgr1=0; optgr2=0 ! dont need gY terms locally.
   optrad=1                     ! do store r-R.

   call nhatgrid(Cryst%atindx1,Cryst%gmet,Wfd%MPI_enreg,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
&    optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%ucvol,Cryst%xred)
   !Pawfgrtab is ready to use

   if (Wfd%pawprtvol>0) then
     call pawfgrtab_print(Pawfgrtab,unit=std_out,prtvol=Wfd%pawprtvol,mode_paral="COLL")
   end if

   ABI_ALLOCATE(Paw_onsite,(Cryst%natom))
   call init_paw_pwaves_lmn(Paw_onsite,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%rprimd,&
&    Psps,Pawtab,Pawrad,Pawfgrtab)

   ABI_ALLOCATE(ur_ae,(Wfd%nfft*Wfd%nspinor))
   ABI_ALLOCATE(data_plot,(Wfd%nfft))

   do plot=1,my_nplots
     band  =my_plot_list(1,plot)
     ik_ibz=my_plot_list(2,plot)
     spin  =my_plot_list(3,plot)

     call wfd_paw_get_aeur(Wfd,band,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae)

     data_plot = DBLE(ur_ae(1:Wfd%nfft)*CONJG(ur_ae(1:Wfd%nfft)))/Cryst%ucvol
     if (Wfd%nspinor==2) &
&      data_plot = data_plot + DBLE(ur_ae(Wfd%nfft+1:)*CONJG(ur_ae(Wfd%nfft+1:)))/Cryst%ucvol

     write(xsf_fname,'(3(a,i0),a)') 'PAW_AE_wfk2_sp',spin,'_kpt',ik_ibz,'_bd',band,'.xsf'
     open(funt,file=xsf_fname,status='unknown',form='formatted')

     call printxsf(n1,n2,n3,data_plot,Cryst%rprimd,(/zero,zero,zero/),&
&      Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,funt,0)

     close(funt)

   end do

   ABI_DEALLOCATE(ur_ae)
   ABI_DEALLOCATE(data_plot)

   call pawfgrtab_free(Pawfgrtab)
   ABI_DEALLOCATE(Pawfgrtab)
   call destroy_paw_pwaves_lmn(Paw_onsite)
   ABI_DEALLOCATE(Paw_onsite)

 else  ! NC case. Just a simple FFT G-->R and then dump the results.

   ABI_ALLOCATE(nc_ur,(Wfd%nfft*Wfd%nspinor))
   ABI_ALLOCATE(data_plot,(Wfd%nfft))

   do plot=1,my_nplots
     band  =my_plot_list(1,plot)
     ik_ibz=my_plot_list(2,plot)
     spin  =my_plot_list(3,plot)

     call wfd_get_ur(Wfd,band,ik_ibz,spin,nc_ur)

     data_plot = DBLE(nc_ur(1:Wfd%nfft)*CONJG(nc_ur(1:Wfd%nfft)))/Cryst%ucvol
     if (Wfd%nspinor==2) &
&      data_plot = data_plot + DBLE(nc_ur(Wfd%nfft+1:)*CONJG(nc_ur(Wfd%nfft+1:)))/Cryst%ucvol

     write(xsf_fname,'(3(a,i0),a)') 'NC_wfk2_sp',spin,'_kpt',ik_ibz,'_bd',band,'.xsf'
     open(funt,file=xsf_fname,status='unknown',form='formatted')

     call printxsf(n1,n2,n3,data_plot,Cryst%rprimd,(/zero,zero,zero/),&
&      Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,funt,0)

     close(funt)
   end do

   ABI_DEALLOCATE(nc_ur)
   ABI_DEALLOCATE(data_plot)
 end if

 ABI_DEALLOCATE(my_plot_list)

 DBG_EXIT("COLL")

end subroutine wfd_plot_ur
!!***

!----------------------------------------------------------------------

#if 0

!!****f* m_wfs/wfd_iterator_bbks
!! NAME
!!  wfd_iterator_bbks
!!
!! FUNCTION
!!  This routines returns an iterator used to loop over the composite index (b,b',k,s)
!!  taking into account the distribution of the ug.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=
!!  [got(Wfd%nproc)]=The number of tasks already assigned to the nodes.
!!  [bbks_mask(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol)]= mask used to select the (b,b',k,s) indeces.
!!
!! OUTPUT
!!  Iter_bbks<iter3_t>=Iterator over the composite index (b,b',k,s) treated by this node.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function wfd_iterator_bbks(Wfd,allup,got,bbks_mask) result(Iter_bbks)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_iterator_bbks'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wfs_descriptor),intent(in) :: Wfd
 character(len=*),intent(in) :: allup
!arrays
 integer,optional,intent(inout) :: got(Wfd%nproc)
 logical,optional,intent(in) :: bbks_mask(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol)
 type(iter3_t) :: Iter_bbks

!Local variables ------------------------------
!scalars
 integer :: ik_ibz,spin,my_nband,ib1,ib2,pcb,how_manyb
 integer :: rank,ncpus,idle,rich_rank
 character(len=500) :: msg
!arrays
 integer :: rank_bandlist(Wfd%mband),get_more(Wfd%nproc)
 integer :: my_band_list(Wfd%mband)
 integer :: whocan(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol,Wfd%nproc)
 logical :: rank_mask(Wfd%nproc)

!************************************************************************

 call iter_alloc(Iter_bbks,(/Wfd%mband,Wfd%nkibz,Wfd%nsppol/))

 whocan=0 ! 1 means that this node can calculate (b,b',k,s)

 do rank=0,Wfd%nproc-1

   do spin=1,Wfd%nsppol
     do ik_ibz=1,Wfd%nkibz

       call wfd_bands_of_rank(Wfd,rank,ik_ibz,spin,how_manyb,rank_bandlist,how="Stored")

       if (how_manyb>0) then

         do pcb=1,how_manyb
           ib2 = rank_bandlist(pcb)
           whocan(rank_bandlist,ib2,ik_ibz,spin,rank+1) = 1
         end do

         !if (starts_with(allup,(/"U","u"/))) then ! only upper triangle of b1,b2 matrix.
         !  do ib2=1,Wfd%nband(ik_ibz,spin)
         !    do ib1=1,ib2-1
         !      whocan(ib1,ib2,ik_ibz,spin) = 0
         !    end do
         !  end do
         !end if

         !if (PRESENT(bbks_mask)) then
         !  where (.not.bbks_mask(:,:,ik_ibz,spin))
         !    whocan(:,:,ik_ibz,spin) = 0
         !  end where
         !end if
       end if
       !if (PRESENT(bbks_mask)) then
       !call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,rank_bandlist,bmask=bbks_mask(:,ik_ibz,spin))
       !else
       !call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,rank_bandlist)
       !end if
       !call iter_push(Iter_bbks,ik_ibz,spin,rank_bandlist(1:my_nband))
     end do
   end do

 end do !rank

 get_more=0; if (PRESENT(got)) get_more=got

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz

     do ib2=1,Wfd%nband(ik_ibz,spin)

       my_nband=0; my_band_list=0
       do ib1=1,ib2-1
         ncpus = COUNT(whocan(ib1,ib2,ik_ibz,spin,:)/=0)

         if (ncpus==1) then
           ! Only one node have the data required. This section should not interfere with the one below.
           rich_rank = imin_loc(ABS(whocan(ib1,ib2,ik_ibz,spin,:)-1)) - 1
           if (Wfd%my_rank==rich_rank) then
             my_nband=my_nband+1
             my_band_list(my_nband)=ib1
           end if

         else if (ncpus>1) then
           ! More nodes might calculate this element. Assign it trying to obtain a good load distribution.
           rank_mask = (whocan(ib1,ib2,ik_ibz,spin,:)==1)
           idle = imin_loc(get_more,mask=rank_mask)
           get_more(idle) = get_more(idle) + 1
           if (Wfd%my_rank==idle-1) then
             my_nband=my_nband + 1
             my_band_list(my_nband) = band
           end if

         else
           call wfd_dump_errinfo(Wfd)
           write(msg,'(a,4(i0,1x))')" Nobody has (band1, band2, ik_ibz, spin): ",ib1,ib2,ik_ibz,spin
           MSG_ERROR(msg)
         end if
       end do ! ib1

       call iter_push(Iter_bbks,ib2,ik_ibz,spin,my_bandlist(1:my_nband))
     end do ! ib2

   end do ! ik_ibz
 end do ! spin

 if (PRESENT(got)) got=get_more

end function wfd_iterator_bbks
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_wfd/make_istwk_table
!! NAME
!! make_istwfk_table
!!
!! FUNCTION
!!
!! INPUTS
!!  ng1,ng2,ng3
!!
!! OUTPUT
!!
!! NOTES
!!   Useful relations:
!!     u_k(G) = u_{k+G0}(G-G0); u_{-k}(G) = u_k(G)^*
!!   and therefore:
!!     u_{G0/2}(G) = u_{G0/2}(-G-G0)^*.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine make_istwfk_table(istwf_k,ng1,ng2,ng3,ig1_inver,ig2_inver,ig3_inver)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_istwfk_table'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng1,ng2,ng3,istwf_k
!arrays
 integer,intent(out) :: ig1_inver(ng1),ig2_inver(ng2),ig3_inver(ng3)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3
 character(len=500) :: msg

!************************************************************************

! Initialize the inverse coordinates
 select case (istwf_k)

 case (1)
   ig1_inver(1)=1
   do i1=2,ng1
     ig1_inver(i1)=ng1+2-i1
   end do
   ig2_inver(1)=1
   do i2=2,ng2
     ig2_inver(i2)=ng2+2-i2
   end do
   ig3_inver(1)=1
   do i3=2,ng3
     ig3_inver(i3)=ng3+2-i3
   end do

 case (2:8)
   if (istwf_k==2 .or. istwf_k==4 .or. istwf_k==6 .or. istwf_k==8) then
     ig1_inver(1)=1
     do i1=2,ng1
       ig1_inver(i1)=ng1+2-i1
     end do
   else
     do i1=1,ng1
       ig1_inver(i1)=ng1+1-i1
     end do
   end if
   if (istwf_k>=2 .and. istwf_k<=5) then
     ig2_inver(1)=1
     do i2=2,ng2
       ig2_inver(i2)=ng2+2-i2
     end do
   else
     do i2=1,ng2
       ig2_inver(i2)=ng2+1-i2
     end do
   end if
   if (istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7) then
     ig3_inver(1)=1
     do i3=2,ng3
       ig3_inver(i3)=ng3+2-i3
     end do
   else
     do i3=1,ng3
       ig3_inver(i3)=ng3+1-i3
     end do
   end if

 case default
   write(msg,'(a,i0)')" Wrong value for istwf_k: ",istwf_k
   MSG_ERROR(msg)
 end select

end subroutine make_istwfk_table
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/check_sym_ug
!! NAME
!! check_sym_ug
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine check_sym_ug(Wfd,Cryst,Kmesh,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_sym_ug'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays

!Local variables ------------------------------
!scalars
 integer :: band,spin,ierr
 integer :: ik_bz,ik_ibz,npw_k
 integer :: k_sym,k_tim,istwf_k
 real(dp) :: err
 logical :: k_isirred
 !character(len=500) :: msg
 type(iskg_tabs_t) :: ISkg
!arrays
 integer,pointer :: kg_k(:,:)
 integer :: k_umklp(3)
 real(dp) :: k_bz(3),k_ibz(3)
 complex(dpc) :: k_eimkt
 complex(gwpc),allocatable :: ug_kbz(:)
 complex(gwpc) :: ur1_kbz(Wfd%nfft*Wfd%nspinor),ur2_kbz(Wfd%nfft*Wfd%nspinor)

!************************************************************************

 ierr = 0
 ABI_CHECK(Wfd%rfft_is_symok," not rfft_is_symok ")

 do spin=1,Wfd%nsppol
   do ik_bz=1,Kmesh%nbz

     call get_BZ_item(Kmesh,ik_bz,k_bz,ik_ibz,k_sym,k_tim,k_eimkt,k_umklp,k_isirred)

     k_ibz   = Kmesh%ibz(:,ik_ibz)

     npw_k   =  Wfd%Kdata(ik_ibz)%npw
     istwf_k =  Wfd%istwfk(ik_ibz)
     kg_k    => Wfd%Kdata(ik_ibz)%kg_k

     call init_iskg_tabs(ISkg,Cryst,k_tim,k_sym,Wfd%ecut,istwf_k,npw_k,k_ibz,kg_k,Cryst%gmet,Wfd%ngfft,ierr)
     ABI_CHECK(ierr==0," init_iskg_tabs returned ierr/=0")

     ABI_ALLOCATE(ug_kbz,(npw_k*Wfd%nspinor))

     do band=1,Wfd%nband(ik_ibz,spin)

       call wfd_get_ur_isk(Wfd,ISkg,Cryst,band,ik_ibz,spin,ur1_kbz,rspace=.TRUE.)
       call wfd_get_ur_isk(Wfd,ISkg,Cryst,band,ik_ibz,spin,ur2_kbz,rspace=.FALSE.)

       err = MAXVAL(ABS(ur2_kbz-ur1_kbz))
       if (err > tol12) then
         write(std_out,*)" (b,k,s) ",band,ik_bz,spin," MAX Error: ",err
         ierr = ierr+1
       end if
     end do

     call destroy_iskg_tabs(ISkg)
     !
     ABI_DEALLOCATE(ug_kbz)
   end do
 end do

end subroutine check_sym_ug
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/init_iskg_tabs
!! NAME
!! init_iskg_tabs
!!
!! FUNCTION
!!  Main creation method for the iskg_tabs_t datatype. This routine constructs
!!  the G-sphere of cutoff energy ecut centered on (IS)k. Then it initializes tables
!!  that are used to reconstruct u_{ISk}(G) from u_k(G) by symmetry.
!!
!! INPUTS
!! itim=2 if time-reversal symmetry is used. 1 otherwise.
!! isym=Index of the symmetry operation S in the array symrec.
!! ecut=cutoff enery.
!! istwf_k=Time-reversal Storage mode for kpt
!! npw_k=Number of planewaves in the sphere.
!! kpt(3)=The k-point that will be rotated via (itim,isym)
!! kg_k(3,npw_k)=Coordinates of the plane waves in the G-sphere centered on kpt
!! gmet(3,3)=Metric in G-space.
!! ngfft(18)=Info on FFT mesh.
!!
!! OUTPUT
!!  ISk<iskg_tabs_t>=Datatype completely initialized.
!!  ierr=Status error.
!!
!! NOTES:
!!  I is either the identity or the inversion (time reversal in reciprocal space).
!!  S is one of the symmetry operation in reciprocal space belonging to the space group of the crystal.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine init_iskg_tabs(ISkg,Cryst,itim,isym,ecut,istwf_k,npw_k,kpt,kg_k,gmet,ngfft,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_iskg_tabs'
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itim,isym,istwf_k,npw_k
 integer,intent(out) :: ierr
 real(dp),intent(in) :: ecut
 type(crystal_structure),intent(in) :: Cryst
 type(iskg_tabs_t),intent(inout) :: ISkg
!arrays
 integer,intent(in) :: kg_k(3,npw_k),ngfft(18)
 real(dp),intent(in) :: kpt(3),gmet(3,3)

!Local variables ------------------------------
!scalars
 integer :: ig1,ig2,isg1_idx,ii,jj,kk,mgfft
 real(dp) :: msk1pgt
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: isg1(3),g2bd(3,2),symrec(3,3)
 integer,allocatable :: box2(:,:,:)
 real(dp) :: tnons(3)
 logical,allocatable :: kg_mask(:)

!************************************************************************

 !Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)

 symrec = Cryst%symrec(:,:,isym)
 tnons  = Cryst%tnons(:,isym)

 ISkg%isym    = isym
 ISkg%itim    = itim

 ISkg%ecut    = ecut
 ISkg%kpt     = (3-2*itim)*MATMUL(symrec, kpt) ! (IS)k
 ISkg%istwf_k = 1
 ! Use time-reversal only if the input k-point already uses it.
 if (istwf_k>1) ISkg%istwf_k = set_istwfk(ISkg%kpt)

 ! TODO this cannot be used yet.
 if (istwf_k/=1 .or. ISkg%istwf_k/=1) then
   MSG_ERROR("time reversal symmetry not yet coded")
 end if
 !
 ! The G-sphere centered in (IS)k.
 call get_kg(ISkg%kpt,ISkg%istwf_k,ISkg%ecut,gmet,ISkg%npw_k,ISkg%kg_k)
 !
 ! Make sure that the two spheres have the same number of PWs.
 if (npw_k/=ISkg%npw_k) then
   !do ig1=1,MIN(npw_k,ISkg%npw_k); write(std_out,*)kg_k(:,ig1),ISkg%kg_k(:,ig1); end do
   write(msg,'(a,2(1x,i0))')" npw_k/=ISkg%npw_k: ",npw_k,ISkg%npw_k
   MSG_ERROR(msg)
 end if
 !
 ! Place the rotated G-sphere in a box to speed up the setup of the table.
 ! box2 gives the correspondence btw the coordinates of ISkg%kg_k and its sequential index.
 ! Set to 0 if the point is outside the sphere
 g2bd(:,1) = MINVAL(ISkg%kg_k,DIM=2)
 g2bd(:,2) = MAXVAL(ISkg%kg_k,DIM=2)

 ABI_ALLOCATE(box2,(g2bd(1,1):g2bd(1,2), g2bd(2,1):g2bd(2,2), g2bd(3,1):g2bd(3,2)))
 box2=0
 do ig2=1,ISkg%npw_k
   ii = ISkg%kg_k(1,ig2)
   jj = ISkg%kg_k(2,ig2)
   kk = ISkg%kg_k(3,ig2)
   box2(ii,jj,kk) = ig2
 end do
 !
 ! For each G in the first sphere:
 !  1) Form SG and check whether it is contained in box2
 !  2) If in box2 and in ISkg%kg_k then store the mapping,
 !     negative indeces are used to signal that cc has to be taken. (useful when istwf_k>2)
 !
 ABI_ALLOCATE(ISkg%gt_k2isk,(npw_k))
 ISkg%gt_k2isk=0
 ABI_ALLOCATE(ISkg%ph_mskpgt,(npw_k))

 ierr=0
 do ig1=1,npw_k
   ! Form (IS)G1.
   isg1 = (3-2*itim)*MATMUL(symrec, kg_k(:,ig1))
   !
   ! Precompute the non-symmorphic phase. Note tha time-reversal is not included here.
   msk1pgt = -(3-2*itim) * two_pi * DOT_PRODUCT(ISkg%kpt+isg1, tnons)
   ISkg%ph_mskpgt(ig1) = DCMPLX( COS(msk1pgt), SIN(msk1pgt) )
   !
   if ( ALL(isg1>=g2bd(:,1)) .and. ALL(isg1<=g2bd(:,2)) ) then
      isg1_idx = box2(isg1(1),isg1(2),isg1(3))
      if (isg1_idx>0) then
        ISkg%gt_k2isk(ig1) = isg1_idx !* (3-2*itim)  ! negative sign means we have to take the cc.
      else
        ierr=ierr+1
        write(std_out,*)" out-of-sphere: ",isg1(:)
      end if
   else
     ierr=ierr+1
     write(std_out,*)" out-of-box: ",isg1(:)
   end if
   !
 end do

 ABI_DEALLOCATE(box2)
 !
 ! Tables for zero-padding FFT.
 mgfft = MAXVAL(ngfft(1:3))

 ABI_ALLOCATE(ISkg%igfft,(ISkg%npw_k))
 ABI_ALLOCATE(kg_mask,(ISkg%npw_k))

 call kgindex(ISkg%igfft,ISkg%kg_k,kg_mask,MPI_enreg_seq,ngfft,ISkg%npw_k)

 !ABI_CHECK(ALL(kg_mask),"FFT para not yet implemented")
 ABI_DEALLOCATE(kg_mask)

 ABI_ALLOCATE(ISkg%gbound,(2*mgfft+8,2))
 call sphereboundary(ISkg%gbound,ISkg%istwf_k,ISkg%kg_k,mgfft,ISkg%npw_k)

end subroutine init_iskg_tabs
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/nullify_iskg_tabs
!! NAME
!! nullify_iskg_tabs
!!
!! FUNCTION
!!  Nullify all pointers defined in iskg_tabs_t datatype.
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine nullify_iskg_tabs(ISkg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_iskg_tabs'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(iskg_tabs_t),intent(inout) :: ISkg

!************************************************************************

!@iskg_tabs_t

! integer pointers.
 nullify(ISkg%gbound  )
 nullify(ISkg%igfft   )
 nullify(ISkg%gt_k2isk)
 nullify(ISkg%kg_k    )

! complex pointers.
 nullify(ISkg%ph_mskpgt)

end subroutine nullify_iskg_tabs
!!***

!----------------------------------------------------------------------

!!****f* m_wfd/destroy_iskg_tabs
!! NAME
!! destroy_iskg_tabs
!!
!! FUNCTION
!!  Deallocate all pointers defined in iskg_tabs_t datatype.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine destroy_iskg_tabs(ISkg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_iskg_tabs'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(iskg_tabs_t),intent(inout) :: ISkg

!************************************************************************

!@iskg_tabs_t

! integer pointers.
 if (associated(ISkg%gbound  ))   then
   ABI_DEALLOCATE(ISkg%gbound)
 end if
 if (associated(ISkg%igfft   ))   then
   ABI_DEALLOCATE(ISkg%igfft)
 end if
 if (associated(ISkg%gt_k2isk))   then
   ABI_DEALLOCATE(ISkg%gt_k2isk)
 end if
 if (associated(ISkg%kg_k     ))  then
   ABI_DEALLOCATE(ISkg%kg_k)
 end if

! complex pointers.
 if (associated(ISkg%ph_mskpgt))  then
   ABI_DEALLOCATE(ISkg%ph_mskpgt)
 end if

end subroutine destroy_iskg_tabs
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_get_ur_isk
!! NAME
!!  wfd_get_ur_isk
!!
!! FUNCTION
!!  Symmetrize a wave function in real space
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=the wavefunction descriptor.
!!  ISkg<iskg_tabs_t>=Info on the symmetry operation to be applied and the rotated G-sphere.
!!  Cryst<crystal_structure>=Structure describing the crystal structure and its symmetries.
!!  band=Band index.
!!  ik_bz=Index of the k-point in the BZ.
!!  spin=Spin index
!!
!! OUTPUT
!!  ur_isk(Wfd%nfft*Wfd%nspinor)=The symmetrized wavefunction in real space.
!!
!! PARENTS
!!      m_wfs
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_get_ur_isk(Wfd,ISkg,Cryst,band,ik_ibz,spin,ur_isk,rspace)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_get_ur_isk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 type(crystal_structure),intent(in) :: Cryst
 type(wfs_descriptor),intent(inout) :: Wfd
 type(iskg_tabs_t),intent(in) :: ISkg
 logical,optional,intent(in) :: rspace
!arrays
 complex(gwpc),intent(out) :: ur_isk(Wfd%nfft*Wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0
 integer :: isym,npw_k,ig1,isg,npw_isk
 complex(dpc) :: ph_mkt
 logical :: use_rspace
 character(len=500) :: msg
!arrays
 integer,pointer :: tabr_k(:)
 real(dp) :: rm1t(3)
 complex(gwpc) :: ur_kibz(Wfd%nfft*Wfd%nspinor)
 complex(gwpc),allocatable :: ug_isk(:)

!************************************************************************

 isym = ISkg%isym

 if (isym==1) then ! Avoid the rotation.
   call wfd_get_ur(Wfd,band,ik_ibz,spin,ur_isk)
   if (ISkg%itim==2) ur_isk = CONJG(ur_isk)
   RETURN
 end if
 !
 ! Reconstruct ur_ISk from ur_k.
 ! Symmetrization is done in real space provided that:
 !   1) real space FFT mesh is compatible with the rotational part of the space group.
 !   2) u_k(r) are stored in memeory so that FFT G--R can be skipped.
 !
 !use_rspace = Wfd%rfft_is_symok .and. wfd_ihave_ur(Wfd,band,ik_ibz,spin,"Stored")
 use_rspace = Wfd%rfft_is_symok .and. wfd_ihave_ur(Wfd,band,ik_ibz,spin)

 if (PRESENT(rspace)) use_rspace = rspace

 npw_k = Wfd%npwarr(ik_ibz)

 SELECT CASE (Wfd%nspinor)

 CASE (1)
   if (use_rspace) then
     ! Symmetrization in done in R-space.
     !  1) FFT G-->R if needed
     !  2) Rotation in R space.
     ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz)
     !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal

     call wfd_get_ur(Wfd,band,ik_ibz,spin,ur_kibz)
     ! TODO recheck this
     rm1t=MATMUL(TRANSPOSE(Cryst%symrec(:,:,isym)),Cryst%tnons(:,isym))
     ph_mkt = EXP(-(0.,1.)*two_pi*DOT_PRODUCT(Wfd%kibz(:,ik_ibz),rm1t))
     tabr_k  => Wfd%irottb(:,isym)  ! Table for rotated FFT points

     ur_isk = ur_kibz(tabr_k) * ph_mkt
     if (ISkg%itim==2) ur_isk = CONJG(ur_isk)
   else
     ! Symmetrization in done in G-space:
     !  1) Rotation in G space
     !  2) FFT G-->R
     ! u_{Sk}(SG) = e^{-iS(k+G).\tau}  u_k(G).
     ! u_{-k}(G)  = u_k(-G)^* for time-reversal symmetry.
     ! u_{-Sk}(-SG) = e^{iS(k+G).\tau}  u_k^*(G).
     ! TODO recheck time-reveral symmetry and istwfk case.
     npw_isk = ISkg%npw_k
     ABI_ALLOCATE(ug_isk,(npw_isk))
     do ig1=1,npw_k
       isg = ISkg%gt_k2isk(ig1)  ! index of (IS)G1 in kg_sym.
       !if (isg==0) stop "isg==0"
       ug_isk(isg) = Wfd%Wave(band,ik_ibz,spin)%ug(ig1) * ISkg%ph_mskpgt(ig1)
       if (ISkg%itim==2) ug_isk(isg) = CONJG(ug_isk(isg))
     end do
     !
     ! zero-padding FFT on the G-sphere centered on (IS)k.
     call fft_onewfn(Wfd%paral_kgb,ISkg%istwf_k,Wfd%nspinor,ISkg%npw_k,Wfd%nfft,Wfd%mgfft,Wfd%ngfft,&
&       ug_isk,ur_isk,ISkg%igfft,ISkg%kg_k,ISkg%gbound,tim_fourdp0,Wfd%MPI_enreg)
     ABI_DEALLOCATE(ug_isk)
   end if

 CASE (2)
   MSG_ERROR("nspinor==2 not implemented.")

 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong value for nspinor: ",Wfd%nspinor
   MSG_ERROR(msg)
 END SELECT

end subroutine wfd_get_ur_isk
!!***

!----------------------------------------------------------------------

!!****f* m_wfs/wfd_grad_ur
!! NAME
!!  wfd_grad_ur
!!
!! FUNCTION
!!   This routine calculates <r|\nabla|u_k>. Two options are coded depending of the
!!   presence of the optional argument ISkg.
!!     (1) If ISkg is NOT present:
!!         <r|\nabla|u_k> is returned where k is the IBZ point corresponding to ik_ibz.
!!     (2) If ISkg is in input:
!!         <r|\nabla|u_ISk> is returned where k is the IBZ point corresponding to ik_ibz.
!!         (I,S) are the set of symmetries used to rotate the point.
!!
!! INPUTS
!!  Wfd<wfs_descriptor>=Wavefunction descriptor.
!!  Cryst<crystal_structure>=Structure describing the crystal structure and its symmetries.
!!  band=Band index.
!!  ik_ibz=Index of the k-point in the IBZ.
!!  spin=Spin index
!!  [ ISkg<iskg_tabs_t> ] =Data needed to symmetrize the wavefunction at the rotated k-point.
!!
!! OUTPUT
!!  grad_ur(Wfd%nfft*Wfd%nspinor,3)=
!!     If (1) ==>  <r|\nabla|u_k>.
!!     If (2) ==>  <r|\nabla|u_ISk>.
!!
!! PARENTS
!!
!! CHILDREN
!!      fft_onewfn
!!
!! SOURCE

subroutine wfd_grad_ur(Wfd,Cryst,band,ik_ibz,spin,grad_ur,ISkg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_grad_ur'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin
 type(crystal_structure),intent(in) :: Cryst
 type(wfs_descriptor),intent(inout) :: Wfd
 type(iskg_tabs_t),optional,intent(in) :: ISkg
!arrays
 complex(gwpc),intent(out) :: grad_ur(Wfd%nfft*Wfd%nspinor,3)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0
 integer :: dir,npw_k,ig,ig_isk,npw_isk,ig_sp,istwf_k,spinor,spad
 logical :: do_sym
 character(len=500) :: msg
!arrays
 integer,pointer :: kg_k(:,:),igfft_k(:),gbound_k(:,:)
 complex(dpc),allocatable :: ug_isk(:)
 complex(gwpc),allocatable :: grad_dir(:)

!************************************************************************

 do_sym = .FALSE.
 if (PRESENT(ISkg)) then
   do_sym = .TRUE.
   !! TODO do_sym = (ISkg%isym/=1 .or. ISkg%itim/=1)
 end if

 if (.not. do_sym) then
   npw_k    = Wfd%npwarr(ik_ibz)
   istwf_k  = Wfd%istwfk(ik_ibz)
   if (istwf_k/=1) then
     ! symmetries of \nabla u are not the same as the wavefunctions.
     ! f_i(G) = -f_i(-G-G0)^* - G0_i u(G)
     ! where f_i(G) =  G_i u(G)
     MSG_ERROR("Not programmed for istwfk/=1")
   end if

   gbound_k => Wfd%Kdata(ik_ibz)%gbound
   kg_k     => Wfd%Kdata(ik_ibz)%kg_k
   igfft_k  => Wfd%Kdata(ik_ibz)%igfft0

   ABI_ALLOCATE(grad_dir,(npw_k*Wfd%nspinor))
   do dir=1,3
     ! Fill array with i G_dir u(G).
     do spinor=1,Wfd%nspinor
       spad = (spinor-1)*npw_k
       do ig=1,npw_k
         ig_sp = ig + spad
         grad_dir(ig_sp)= j_dpc * kg_k(dir,ig) * Wfd%Wave(band,ik_ibz,spin)%ug(ig_sp)
       end do
     end do

     call fft_onewfn(Wfd%paral_kgb,istwf_k,Wfd%nspinor,npw_k,Wfd%nfft,Wfd%mgfft,Wfd%ngfft,&
&       grad_dir,grad_ur(:,dir),igfft_k,kg_k,gbound_k,tim_fourdp0,Wfd%MPI_enreg)
   end do

   ABI_DEALLOCATE(grad_dir)

 else
  !
  ! Reconstruct ur in BZ from the corresponding wavefunction in IBZ.
  SELECT CASE (Wfd%nspinor)
  CASE (1)

    if (Wfd%istwfk(ik_ibz)/=1 .or. ISkg%istwf_k/=1) then
      ! symmetries of \nabla u are not the same as the wavefunctions.
      MSG_ERROR("Not programmed for istwfk/=1")
    end if

    ! symmetrization in done in G-space: rotation + FFT G-->R
    ! u_{Sk}(SG) = e^{-iS(k+G).\tau}  u_k(G).
    ! u_{-k}(G)  = u_k(-G)^* for time-reversal symmetry.
    ! u_{-Sk}(-SG) = e^{iS(k+G).\tau}  u_k^*(G).
    ! TODO recheck time-reveral symmetry and istwfk case.
    npw_k = Wfd%npwarr(ik_ibz)
    npw_isk = ISkg%npw_k
    ABI_ALLOCATE(ug_isk,(npw_isk))

    do ig=1,npw_k
      ig_isk = ISkg%gt_k2isk(ig)  ! index of (IS)G1 in kg_sym.
      !if (ig_isk==0) stop "ig_isk==0"
      ug_isk(ig_isk) = Wfd%Wave(band,ik_ibz,spin)%ug(ig) * ISkg%ph_mskpgt(ig)
      if (ISkg%itim==2) ug_isk = CONJG(ug_isk)
    end do

    ABI_ALLOCATE(grad_dir,(npw_isk))

    do dir=1,3
      do ig=1,npw_isk
        grad_dir(ig)= j_dpc * ISkg%kg_k(dir,ig) * ug_isk(ig)
      end do
      !
      ! zero-padding FFT on the G-sphere centered on ISk.
      call fft_onewfn(Wfd%paral_kgb,ISkg%istwf_k,Wfd%nspinor,npw_isk,Wfd%nfft,Wfd%mgfft,Wfd%ngfft,&
&        grad_dir,grad_ur(:,dir),ISkg%igfft,ISkg%kg_k,ISkg%gbound,tim_fourdp0,Wfd%MPI_enreg)
    end do
    ABI_DEALLOCATE(ug_isk)
    ABI_DEALLOCATE(grad_dir)

  CASE (2)
    MSG_ERROR("Not implemented error")

  CASE DEFAULT
    write(msg,'(a,i0)')" Wrong value for nspinor: ",Wfd%nspinor
    MSG_ERROR(msg)
  END SELECT

 end if

 RETURN
 ABI_UNUSED(Cryst%nsym)

end subroutine wfd_grad_ur
!!***

END MODULE m_wfs
!!***

