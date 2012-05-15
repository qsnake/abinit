#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_screen

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_splines
 use m_lebedev

 use m_gwdefs,         only : GW_TOLQ0, czero_gw, cone_gw
 use m_fstrings,       only : starts_with
 use m_io_tools,       only : get_unit
 use m_blas,           only : xgemv
 use m_numeric_tools,  only : print_arr
 use m_geometry,       only : normv
 use m_crystal,        only : crystal_structure, destroy_crystal
 use m_bz_mesh,        only : bz_mesh_type, get_BZ_item, box_len, has_bz_item, find_qmesh, destroy_bz_mesh_type
 use m_gsphere,        only : gvectors_type, init_gsphere, destroy_gsphere
 use m_vcoul,          only : vcoul_t
 use m_io_screening,   only : free_scrhdr, nullify_HScr, scr_hdr_io, read_screening, write_screening, &
&                             copy_scrhdr, scrhdr_type, scrhdr_comm
 use m_spectra,        only : spectra_type, destroy_spectra, repr_dielconst
 use m_ppmodel,        only : ppmodel_type, ppm_init, ppm_free, ppm_nullify, PPM_NONE, new_setup_ppmodel, ppm_symmetrizer
 use m_screening

 implicit none

 private 

 ! Flags defining the content of the %mat buffer in the fgg_t type. 
 integer,public,parameter :: MAT_NOTYPE         = 0
 integer,public,parameter :: MAT_CHI0           = 1
 integer,public,parameter :: MAT_CHI            = 2
 integer,public,parameter :: MAT_EPSILON        = 3
 integer,public,parameter :: MAT_INV_EPSILON    = 4
 integer,public,parameter :: MAT_INV_EPSILON_M1 = 5
 integer,public,parameter :: MAT_W              = 6
 integer,public,parameter :: MAT_W_M1           = 7
 !
 ! Family vertex.
 integer,public,parameter :: VTX_FAMILY_NONE  = 0   ! No vertex correction.
 integer,public,parameter :: VTX_FAMILY_TDDFT = 1   ! TDDFT-based vertex.
 integer,public,parameter :: VTX_FAMILY_ADA   = 2   ! ADA vertex.
 !
 ! Test charge or test particle.
 integer,public,parameter :: VTX_TEST_CHARGE   = 0
 integer,public,parameter :: VTX_TEST_PARTICLE = 1   
 !
 ! Named constants for the frequency mesh.
 integer,private,parameter :: WMESH_LINEAR    = 1
 integer,private,parameter :: WMESH_GAUSS_LEG = 2
 integer,private,parameter :: WMESH_TAN_GRID  = 3
 !
 ! Method used for the frequency integration.
 integer,public,parameter ::  WINT_NONE    = 0
 integer,public,parameter ::  WINT_PPMODEL = 1
 integer,public,parameter ::  WINT_CONTOUR = 2
 integer,public,parameter ::  WINT_AC      = 3
 !
 ! Parameters used for the model dielectric function.
 integer,public,parameter :: MDL_NONE      = 0
 integer,public,parameter :: MDL_BECHSTEDT = 1
 !
 ! Flags giving the status of the local buffers defined in fgg_t.
 integer,private,parameter :: MAT_NODATA    = 0
 integer,private,parameter :: MAT_ALLOCATED = 1 
 integer,private,parameter :: MAT_STORED    = 2
 !
 ! Flags giving the status of the local buffers defined in fgg_t.
 integer,private,parameter :: FGG_QBZ_ISPOINTER  =1 ! Fgg_qbz is used to store the address in memory.
 integer,private,parameter :: FGG_QBZ_ISALLOCATED=2 ! Fgg_qbz is used as an allocable array.
!!***

!----------------------------------------------------------------------

!!****t* m_screen/screen_info_t
!! NAME
!! screen_info_t
!!
!! FUNCTION
!!  Container storing the parameters used to initialize a screen_t datatype or to 
!!  calculate a new SCR file from the SUSC file containing the independent-particle
!!  polarizability.
!!
!! NOTES
!!  The list of parameters passed to screening_init is copied in W%Info.
!!  At present there is no need to provide a copy method since the structure does 
!!  not contain pointers but such a method must be defined and used if the 
!!  dynamic entities are added to the datatype.
!!
!! SOURCE
                          
type,public :: screen_info_t

  integer :: mat_type = MAT_NOTYPE
  ! Matrix identifier. See MAT_* flags.

  integer :: vtx_family = VTX_FAMILY_NONE
  ! Vertex correction family.

  integer :: ixc=0
  ! XC functional used for the TDDFT-based vertex.

  integer :: use_ada=0
  ! >0 if ADA vertex is used.

  integer :: use_mdf = MDL_NONE
  ! >0 if model dielectric function is used.

  integer :: use_ppm = PPM_NONE
  ! >0 if ppmodel is used.

  integer :: vtx_test = VTX_TEST_CHARGE
  ! test charge or test particle.

  integer :: wint_method = WINT_NONE
  ! Defines the frequency integration technique. See WIN_ flags.
  ! NOTE that this flag can be changed at run time. For example 
  ! one can switch from the CD to the PPm if the ppmodel parameters are in memory 

  real(dp) :: ada_kappa = 2.1_dp
  ! Inverse smearing length used for ADA.

  real(dp) :: eps_inf = 12.0_dp           
  ! Dielectric constant used for the model dielectric function.

  real(dp) :: drude_plsmf = zero
  ! Drude plasma frequency used for PPmodel 1.

end type screen_info_t
!!***

 public :: screen_info_print

!----------------------------------------------------------------------

!!****t* m_screen/fgg_t
!! NAME
!! fgg_t
!!
!! FUNCTION
!!  Object used to store F(G,G')(w) for a given q-point.
!!
!! SOURCE

 type,public :: fgg_t

  integer :: nomega
  ! Number of frequencies.

  integer :: npw       
  ! Number of G vectors.

  integer :: nqlwl
  ! Number of points for the treatment of the long wave-length limit.

  integer :: has_mat = MAT_NODATA
  ! Flag giving the status of mat.
                                    
  integer :: has_wings = MAT_NODATA
  ! Flag giving the status of ur.

  !arrays
  complex(gwpc),pointer :: mat(:,:,:)  SET2NULL
  ! mat(npw,npw,nomega)
  ! The component of the two-point function $F_{G,G',w}$ for a given q.

  complex(dpc),pointer :: lwing(:,:,:)  SET2NULL
  ! lwing(npw,nqlwl,nomega)
  ! Lower wings for the different q --> 0.

  complex(dpc),pointer :: uwing(:,:,:)  SET2NULL
  ! uwing(npw,nqlwl,nomega)
  ! Upper wings for the different q --> 0.

  complex(dpc),pointer :: head(:,:,:)  SET2NULL
  ! uwing(nqlwl,nqlwl,nomega)
  ! Head for the different q --> 0.

 end type fgg_t   
!!***

 public :: fgg_nullify    ! Nullify all pointers before use.
 public :: fgg_init       ! Creation method.
 public :: fgg_free       ! Free all associated pointers.

 interface fgg_nullify
   module procedure fgg_nullify_0D
 end interface fgg_nullify

 interface fgg_free
   module procedure fgg_free_0D
   module procedure fgg_free_1D
 end interface fgg_free

!----------------------------------------------------------------------

!!****t* m_screen/screen_t
!! NAME
!! screen_t
!!
!! FUNCTION
!!
!! SOURCE

 type,public :: screen_t

  integer :: accesswff            ! Flag defining the IO mode.
  integer :: debug_level=0        ! Internal Flag defining the debug level.
  integer :: mqmem                ! =0 for out-of-core solution, =nqibz if entire matrix is stored in memory.
  integer :: nI,nJ                ! Number of components (rows,columns) in chi|eps^-1. (1,1) if collinear.
  integer :: nqibz                ! Number of q-points in the IBZ used.
  integer :: nqlwl                ! Number of points used for the treatment of the long wave-length limit.
  integer :: nomega               ! Total Number of frequencies used.
  integer :: nomega_i             ! Number of purely imaginary frequencies used.
  integer :: nomega_r             ! Number of real frequencies used.
  integer :: npw                  ! Number of G vectors.
  integer :: prtvol

  integer :: has_ppmodel 
  ! 1 if PPmodel tables are stored. 

  integer :: has_fgg
  ! 1 if Fgg tables are stored. 

  !%integer :: wmesh_type
  !%real(dp) :: ecut

  integer :: nfftf_tot
  integer :: nspden

  integer :: ngfftf(18)          ! Info on the FFT mesh used for ae_rhor (used for the model dielectric function)

  real(dp),pointer :: ae_rhor(:,:)
  ! ae_rhor(nfft,nspden)
  ! Density in real space used to construct the TDDFT kernel or the model dielectric function.
  ! NOTE that ae_rhor is given on the dense mesh as it contains the PAW onsite contribution.

  character(len=fnlen) :: fname=ABI_NOFILE  ! Name of the file used for the out-of-core solution.

!arrays
  real(dp),pointer :: qibz(:,:)  SET2NULL
  ! qibz(3,nqibz)
  ! q-points in reduced coordinates

  !type(bz_mesh_type) :: Qmesh
  ! Info the q-point sampling.

  real(dp),pointer :: qlwl(:,:) SET2NULL
  ! qlwl(3,nqlwl)
  ! q-points used for the long wave-length limit treatment.

  complex(dpc),pointer :: omega(:)  SET2NULL
  ! omega(nomega)
  ! List of frequencies. Real frequencies are packed first.

  integer,pointer :: gvec(:,:)   SET2NULL
  ! gvec(3,npw)
  ! G-vectors used to describe the two-point function (r.l.u.).

  !type(gvectors_type) :: Gsphere
  ! Info on the G-sphere. Note that the basis set does not depend on q.
  ! The G-vectors are ordered in shells to facilitate the application of symmetry operations.
  ! See m_gsphere.F90.

  !type(vcoul_t) :: Vcp

  logical,pointer :: keep_q(:)  SET2NULL   
   ! keep_q(nqibz)
   ! Storage strategy: keep or not keep Em1(q) in memory.

  type(fgg_t),pointer :: Fgg(:)
  ! Fgg(nqibz)
  ! F_{G,G'}(q,w) for q in the IBZ.

  integer :: fgg_qbz_stat=FGG_QBZ_ISPOINTER
  ! Status of Fgg_qbz

  integer :: fgg_qbz_idx=0
  ! The index of the q-point in BZ pointed by Fgg_qbz. Used for debugging purpose.

  type(fgg_t),pointer :: Fgg_qbz     SET2NULL
  ! Buffer used for storing F_GG' at the point q_bz in the BZ
  ! If q_bz is in the IBZ, Fgg_qbz points to Fgg(iq_ibz)
  ! If q_bz is not in the IBZ, Fgg_qbz is allocated and used to store the symmetrized matrix.

  type(ppmodel_type) :: PPm
  ! Structure storing the plasmon-pole paramenters.

  type(screen_info_t) :: Info
  ! Paramenters used to construct the screening.

 end type screen_t
!!***

 public :: screen_nullify       ! Nullify all pointers before use.
 public :: screen_init          ! Creation method.
 public :: screen_free          ! Free all associated pointers.
 public :: screen_symmetrizer   ! Prepare the object for applying W_qbz.
 public :: screen_w0gemv        ! Matrix vector multiplication \sum_{G'} F_{G,G') |u(G')>.
 public :: screen_times_ket     
 public :: screen_mdielf

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_info_print
!! NAME
!!  screen_info_print
!!
!! FUNCTION
!!  Printout of screen_info_t.
!!
!! INPUTS
!!  W_info<screen_info_t>=The structure.
!!  [unit]=Unit number for output
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS"
!!  [header]=String to be printed as header for additional info.
!!
!! OUTPUT
!!  Only printing 
!!
!! PARENTS
!!      m_screen
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine screen_info_print(W_info,header,unit,mode_paral,prtvol) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_info_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral 
 character(len=*),optional,intent(in) :: header
 type(screen_info_t),intent(in) :: W_info

!Local variables-------------------------------
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg      
! ********************************************************************* 

 !@screen_info_t
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol 
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the screen_info_t% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(std_out,*)" mat_type    ",W_info%mat_type
 write(std_out,*)" vtx_family  ",W_info%vtx_family
 write(std_out,*)" ixc         ",W_info%ixc
 write(std_out,*)" use_ada     ",W_info%use_ada
 write(std_out,*)" use_mdf     ",W_info%use_mdf
 write(std_out,*)" use_ppm     ",W_info%use_ppm
 write(std_out,*)" vtx_test    ",W_info%vtx_test
 write(std_out,*)" wint_method ",W_info%wint_method

 write(std_out,*)" ada_kappa   ",W_info%ada_kappa
 write(std_out,*)" eps_inf     ",W_info%eps_inf
 write(std_out,*)" drude_plsmf ",W_info%drude_plsmf

end subroutine screen_info_print
!!***

!----------------------------------------------------------------------

!!****f* m_screen/fgg_nullify_0D
!! NAME
!! fgg_nullify_0D
!!
!! FUNCTION
!! Initialize the pointers to null()
!!
!! INPUTS
!! Fgg<fgg_t>=The data structure to be nullified.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine fgg_nullify_0D(Fgg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fgg_nullify_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(fgg_t),intent(inout) :: Fgg
! *************************************************************************

 !@fgg_t
 !complex
 nullify(Fgg%mat  ) 
 nullify(Fgg%lwing)
 nullify(Fgg%uwing)       
 nullify(Fgg%head )

end subroutine fgg_nullify_0D
!!***

!----------------------------------------------------------------------

!!****f* m_screen/fgg_free_0D
!! NAME
!! fgg_free_0D
!!
!! FUNCTION
!! Deallocate all the pointers in the structure that result to be associated.
!!
!! PARENTS
!!      m_screen
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine fgg_free_0D(Fgg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fgg_free_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(fgg_t),intent(inout) :: Fgg
! *************************************************************************

 !@fgg_t
 !complex
 if (associated(Fgg%mat))  then
   ABI_DEALLOCATE(Fgg%mat)
 end if
 Fgg%has_mat = MAT_NODATA

 if (associated(Fgg%lwing))  then
   ABI_DEALLOCATE(Fgg%lwing)
 end if
 if (associated(Fgg%uwing))  then
   ABI_DEALLOCATE(Fgg%uwing)
 end if
 if (associated(Fgg%head ))  then
   ABI_DEALLOCATE(Fgg%head)
 end if
 Fgg%has_wings = MAT_NODATA

end subroutine fgg_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_screen/fgg_free_1D
!! NAME
!! fgg_free_1D
!!
!! FUNCTION
!! Deallocate all the pointers in the structure that result to be associated.
!!
!! INPUT
!!  [keep_q(:)]=Optional logical mask used to select the q-points that are deallocated.
!!
!! PARENTS
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine fgg_free_1D(Fgg,keep_q)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fgg_free_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(fgg_t),intent(inout) :: Fgg(:)
 logical,optional,intent(in) :: keep_q(:)

!Local variables ------------------------------
!scalars
 integer :: iq_ibz
 logical :: keep_it
! *************************************************************************

 do iq_ibz=LBOUND(Fgg,DIM=1),UBOUND(Fgg,DIM=1)
   keep_it = .FALSE.; if (PRESENT(keep_q)) keep_it = keep_q(iq_ibz)
   if (.not.keep_it) call fgg_free_0D(Fgg(iq_ibz))
 end do

end subroutine fgg_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_screen/fgg_init
!! NAME
!! fgg_init
!!
!! FUNCTION
!! Initialize the structure allocating the memory and initializing the internal variables.
!!
!! INPUT
!!  npw 
!!  nqlwl
!!
!! SIDE EFFECTS
!!  Fgg, See function description
!!
!! PARENTS
!!      m_screen
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine fgg_init(Fgg,npw,nomega,nqlwl)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fgg_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nqlwl,nomega
 type(fgg_t),intent(inout) :: Fgg

!Local variables ------------------------------
!scalars
 integer :: istat
! *************************************************************************

 !@fgg_t
 Fgg%nomega= nomega
 Fgg%npw   = npw
 Fgg%nqlwl = nqlwl

 if (npw>0.and.nomega>0) then
   ABI_ALLOCATE(Fgg%mat,(npw,npw,nomega))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"out of memory in Fgg%mat")
   Fgg%has_mat = MAT_ALLOCATED
 end if

 if (nqlwl>0.and.nomega>0) then
   ABI_ALLOCATE(Fgg%head,(nqlwl,nqlwl,nomega))
   istat = ABI_ALLOC_STAT
   ABI_ALLOCATE(Fgg%lwing,(npw,nqlwl,nomega))
   istat = ABI_ALLOC_STAT
   ABI_ALLOCATE(Fgg%uwing,(npw,nqlwl,nomega))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"out of memory in Fgg%wings")
   Fgg%has_wings = MAT_ALLOCATED
 end if

end subroutine fgg_init
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_fgg_qbz_set
!! NAME
!!  screen_fgg_qbz_set
!!
!! FUNCTION
!!  Helper function used to perform the setup W%Fgg_qbz setting also the internal 
!!  flag that defines its status.
!!
!! INPUTS
!!  W<screen_t>=Object describing the screened interaction.
!!  iq_bz=Index of the q-point in the BZ.
!!  nqlwl=Number of wings wanted.
!!  how= "Pointer" is a true pointer is wanted.
!!       "Allocated" if memory has to be allocated.
!!
!! NOTES
!!  iq_bz and nqlwl are not used if how="Pointer".
!!
!! PARENTS
!!      m_screen
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine screen_fgg_qbz_set(W,iq_bz,nqlwl,how) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_fgg_qbz_set'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nqlwl
 character(len=*),intent(in) :: how
 type(screen_t),intent(inout) :: W

!Local variables ------------------------------
!scalars
 character(len=500) :: msg

!************************************************************************

 W%fgg_qbz_idx = iq_bz ! For Debugging.

 if (starts_with(how,(/"P"/)) ) then ! We want a pointer.
   !
   select case (W%fgg_qbz_stat)

   case (FGG_QBZ_ISALLOCATED)
     call fgg_free_0D(W%Fgg_qbz)
     ABI_DEALLOCATE(W%Fgg_qbz)
     nullify(W%Fgg_qbz)
     W%fgg_qbz_stat = FGG_QBZ_ISPOINTER

   case (FGG_QBZ_ISPOINTER)
     nullify(W%Fgg_qbz) ! Return null().

   case default
     write(msg,'(a,i0)')" Wrong status =",W%fgg_qbz_stat
     MSG_ERROR(msg)
   end select
   !
 else if (starts_with(how,(/"A"/)) ) then ! We want an allocatable array.
   !
   select case (W%fgg_qbz_stat)

   case (FGG_QBZ_ISPOINTER) ! Allocate space.ABI_ALLOCATE(FGG_QBZ_ISPOINTER,)
     nullify(W%Fgg_qbz)
     ABI_ALLOCATE(W%Fgg_qbz,)
     
     call fgg_init(W%Fgg_qbz,W%npw,W%nomega,nqlwl) 
     W%fgg_qbz_stat = FGG_QBZ_ISALLOCATED

   case (FGG_QBZ_ISALLOCATED)
     !W%Fgg_qbz%mat = czero_gw
     W%Fgg_qbz%has_mat = MAT_ALLOCATED  ! STORED --> ALLOCATED
     if (nqlwl>0 .and. W%Fgg_qbz%has_wings == MAT_NODATA) then ! Unlikely but oh well.
       call fgg_free_0D(W%Fgg_qbz)
       call fgg_init(W%Fgg_qbz,W%npw,W%nomega,nqlwl)
       W%fgg_qbz_stat = FGG_QBZ_ISALLOCATED
     end if

   case default
     write(msg,'(a,i0)')" Wrong status =",W%fgg_qbz_stat
     MSG_ERROR(msg)
   end select
   !
 else 
   MSG_BUG("Wrong how: "//how)
 end if

end subroutine screen_fgg_qbz_set
!!***                                                                                                

!----------------------------------------------------------------------

!!****f* m_screen/screen_ihave_fgg
!! NAME
!!  screen_ihave_fgg
!!
!! FUNCTION
!!  Inquire the processor whether it has a particular F_{GG')(q) and with which status.
!!
!! INPUTS
!!  W<screen_t>=Object describing the screened interaction.
!!  iq_ibz=k-point index
!!  [how]=string defining which status is checked. By default the function returns
!!     .TRUE. if the wave is either MAT_ALLOCATED or MAT_STORED.
!!     Possible mutually exclusive values: "Allocated", "Stored". 
!!     Only the first character is checked (no case-sensitive)
!!
!! NOTES
!!   A zero index can be used to inquire the status of the full set of q-points.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function screen_ihave_fgg(W,iq_ibz,how) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_ihave_fgg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_ibz
 logical :: screen_ihave_fgg
 character(len=*),optional,intent(in) :: how
 type(screen_t),intent(in) :: W

!Local variables ------------------------------
!scalars
 integer :: ii
 !character(len=500) :: msg
!arrays
 integer :: check(2) 

!************************************************************************

 check = (/MAT_ALLOCATED, MAT_STORED/)
 if (PRESENT(how)) then
   if (starts_with(how,(/"A","a"/))) check = (/MAT_ALLOCATED, MAT_ALLOCATED/)
   if (starts_with(how,(/"S","s"/))) check = (/MAT_STORED, MAT_STORED/)
 end if

 if (iq_ibz>0) then
   screen_ihave_fgg = ( W%Fgg(iq_ibz)%has_mat==check(1) .or.&
&                    W%Fgg(iq_ibz)%has_mat==check(2) )
 else 
   ! check the status of the full set of q-tables.
   screen_ihave_fgg=.TRUE.
   do ii=1,W%nqibz
     screen_ihave_fgg = screen_ihave_fgg .and. &
&                    ( W%Fgg(ii)%has_mat==check(1) .or.&
&                      W%Fgg(ii)%has_mat==check(2) )
   end do
 end if

end function screen_ihave_fgg
!!***                                                                                                

!----------------------------------------------------------------------

!!****f* m_screen/screen_nullify
!! NAME
!! screen_nullify
!!
!! FUNCTION
!! Initialize the pointers to null()
!!
!! SIDE EFFECTS
!! W<screen_t>=The data structure to be nullified.
!!
!! PARENTS
!!      bethe_salpeter,m_screen
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine screen_nullify(W)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(screen_t),intent(inout) :: W
! *************************************************************************

 !@screen_t
 !
 ! integer
 nullify(W%gvec)
 !
 !real
 nullify(W%ae_rhor)
 nullify(W%qibz)
 nullify(W%qlwl)
 !
 !logical
 nullify(W%keep_q)
 !
 !complex
 nullify(W%omega)       
 !
 !types
 nullify(W%Fgg_qbz); W%fgg_qbz_stat=FGG_QBZ_ISPOINTER ! Nedded since the initial status is undefined.
 nullify(W%Fgg)

 call ppm_nullify(W%PPm)

end subroutine screen_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_free
!! NAME
!! screen_free
!!
!! FUNCTION
!! Free the memory allocate in the datatype.
!!
!! SIDE EFFECTS
!! W<screen_t>=The data structure to be nullified.
!!
!! OUTPUT
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine screen_free(W)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(screen_t),intent(inout) :: W

! *************************************************************************

 !@screen_t
 !
 ! integer
 if (associated(W%gvec))  then
   ABI_DEALLOCATE(W%gvec)
 end if
 !
 !real
 if (associated(W%ae_rhor))  then
   ABI_DEALLOCATE(W%ae_rhor)
 end if
 if (associated(W%qibz   ))  then
   ABI_DEALLOCATE(W%qibz)
 end if
 if (associated(W%qlwl   ))  then
   ABI_DEALLOCATE(W%qlwl)
 end if
 !
 !complex
 if (associated(W%omega))  then
   ABI_DEALLOCATE(W%omega)
 end if

 ! logical
 if (associated(W%keep_q))  then
   ABI_DEALLOCATE(W%keep_q)
 end if
 !
 ! types ------------------------------------------
 !
 ! Here be careful with dangling pointers.
 ! First Fgg_qbz that might point to one of the %Fgg then %Fgg.
 SELECT CASE (W%fgg_qbz_stat)
 CASE (FGG_QBZ_ISALLOCATED)
   call fgg_free_0D(W%Fgg_qbz)
   W%fgg_qbz_stat=FGG_QBZ_ISPOINTER
   !
 CASE (FGG_QBZ_ISPOINTER)
   nullify(W%Fgg_qbz)
   W%fgg_qbz_stat=FGG_QBZ_ISPOINTER
   !
 CASE DEFAULT 
   continue
 END SELECT
 !
 ! Free the Fgg matrices.
 if (associated(W%Fgg)) then 
   call fgg_free(W%Fgg)
   ABI_DEALLOCATE(W%Fgg)
 end if
 !
 ! Free the plasmon pole tables.
 call ppm_free(W%PPm)

end subroutine screen_free
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_init
!! NAME
!!  screen_init
!!
!! FUNCTION
!!  Initialize basic dimensions and the important (small) arrays in an Epsilonm1_results data type
!!  starting from a file containing either epsilon^{-1} (_SCR) or chi0 (_SUSC).
!!
!! INPUTS
!!  W_Info<screen_info_t>=The list of parameters used to construct the screen function.
!!  Dtset<dataset_type>=Used to call rhohxc when kernel corrections are wanted
!!  Cryst<crystal_structure>=Info on the unit cell.
!!  Qmesh<bz_mesh_type>=Info on the Q-mesh.
!!  Gsph<gvectors_type>=Info on the plane-wave basis set used for the two-point function.
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction.
!!  ifname=The name of the external file used to read the matrix.
!!  id_required=Identifier used to specify the type of two-point function that is wanted.
!!  accesswff=Option defining the file format of the external file.
!!  mqmem=0 for out-of-core solution, /=0 if entire matrix has to be stored in memory. 
!!  npw_asked=Number of G-vector to be used in the calculation, if <=0 use Max allowed number.
!!  ngfftf(18)=Info on the (fine) mesh used for the density.
!!  nfftf_tot=Total number of point in the FFT mesh for ae_rhor 
!!  nsppol=Numer of independent spin polarizations.
!!  nspden=Number of spin density components in ae_rhor
!!  ae_rhor(nfftf_tot,nspden)
!!  prtvol=Verbosity level.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  W<screen_t>=The structure initialized with basic dimensions and arrays.
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine screen_init(W,W_Info,Dtset,Cryst,Qmesh,Gsph,Vcp,ifname,mqmem,npw_asked,&
&  accesswff,ngfftf,nfftf_tot,nsppol,nspden,ae_rhor,prtvol,comm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_init'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mqmem,accesswff,npw_asked,comm,prtvol,nsppol
 integer,intent(in) :: nfftf_tot,nspden
 character(len=fnlen),intent(in) :: ifname
 type(dataset_type),intent(in) :: Dtset
 type(crystal_structure),intent(in) :: Cryst
 type(gvectors_type),intent(in) :: Gsph
 type(vcoul_t),intent(in) :: Vcp
 type(bz_mesh_type),intent(in) :: Qmesh
 type(screen_info_t),intent(in) :: W_Info
 type(screen_t),intent(out) :: W
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: ae_rhor(nfftf_tot,nspden)

!Local variables-------------------------------
!scalars
 integer :: option_test,approx_type,ixc_required,id_required !,nkxc
 integer :: ios,fform,rdwr,my_rank,master,unt,mat_type_read
 integer :: nqibz,nomega,iq_ibz,npw,nqlwl !,istat
 integer :: nI,nJ,comm_self,iq_bz,mdf_type,ppmodel,ierr
 integer :: iw,qsort,ii
 real(dp) :: eps_inf,drude_plsmf
 logical :: deallocate_Fgg,found,has_file,is_qeq0,remove_dgg !only_one_kpt,
 character(len=500) :: msg                   
 character(len=fnlen) :: sus_fname,scr_fname
 type(ScrHdr_type) :: Hscr
!arrays
 integer :: g0(3),iperm(Qmesh%nibz)
 real(dp) :: wt_list(Qmesh%nibz)
 complex(gwpc),pointer :: em1_ggw(:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 ABI_UNUSED(nsppol)

 !@screen_t
 master    = 0
 my_rank   = xcomm_rank(comm)
 comm_self = xmpi_self

 call screen_nullify(W)
 !
 ! Initialize basic parameters ----------------------------------------
 !
 !@screen_info_t  
 W%Info = W_info ! Copy basic parameters.
 call screen_info_print(W%Info,header="W info",unit=std_out)

 id_required  = W_Info%mat_type
 approx_type  = W_Info%vtx_family
 option_test  = W_Info%vtx_test
 ixc_required = W_Info%ixc 

 !if (ALL(id_required /= (/MAT_W_M1, MAT_W/)) ) then
 if (ALL(id_required /= (/MAT_INV_EPSILON/)) ) then
   write(msg,'(a,i0,a)')"id_required ",id_required," not available"
   MSG_ERROR(msg)
 end if

 ! This part must be rationalized.
 remove_dgg = (id_required==MAT_W_M1) 

 if (W%Info%use_mdf == MDL_NONE) W%fname = ifname
 W%nI = 1
 W%nJ = 1
 !
 ! The Q-point sampling is inited from Qmesh.
 W%nqibz=Qmesh%nibz
 ABI_ALLOCATE(W%qibz,(3,W%nqibz))
 W%qibz= Qmesh%ibz

 W%mqmem=mqmem; if (W%mqmem/=0) W%mqmem=W%nqibz
                                                       
 ABI_ALLOCATE(W%keep_q,(W%nqibz))
 W%keep_q=.TRUE.
 !
 if (W%mqmem<W%nqibz) then  ! Keep in memory the most representative q-points.
   W%keep_q=.FALSE.
   wt_list = Qmesh%wt; iperm = (/(ii,ii=1,Qmesh%nibz)/)
   call sort_dp(Qmesh%nibz,wt_list,iperm,tol12)
   do qsort=Qmesh%nibz,Qmesh%nibz-mqmem+1,1
     iq_ibz = iperm(qsort)
     W%keep_q(iq_ibz) = .TRUE.
   end do
 end if

 W%accesswff = accesswff
 W%prtvol    = prtvol

 W%has_ppmodel=0; if (W%Info%use_ppm/=PPM_NONE) W%has_ppmodel=1
 !
 ! Copy the AE density for the model dielectric function or for the vertex corrections.
 W%nspden     = nspden
 W%ngfftf     = ngfftf
 W%nfftf_tot  = nfftf_tot

 ABI_ALLOCATE(W%ae_rhor,(nfftf_tot,nspden))
 W%ae_rhor = ae_rhor
 
 deallocate_Fgg=.FALSE.
                                                                             
 W%has_fgg=0; if ( ANY(W%Info%wint_method == (/WINT_CONTOUR, WINT_AC/)) ) W%has_fgg=1
                                                                             
 if (W%has_fgg>0 .and. W%has_ppmodel>0) then
   MSG_WARNING("Both PPmodel tables and F_(GG')(q,w) are stored in memory") 
 end if

 ! G-sphere for W and Sigma_c is initialized from the SCR file.
 !%% only_one_kpt=.FALSE.
 !%% call init_gsphere(Gsph,only_one_kpt,Cryst,Sigp%npwc,gvec=Er%gvec)
 !%% deallocate(gvec_p)

 ! Default values used if external file is not read.
 nqlwl  = 0
 nomega = 1
 W%npw  = npw_asked

 ! Model dielectric function does not require any external file.
 has_file = (W%Info%use_mdf==MDL_NONE)

 if (has_file) then ! Open file and check its content.
   if (my_rank==master) then
     call wrtout(std_out,' Testing file: '//TRIM(W%fname),'COLL')
     unt=get_unit()
     open(unit=unt,file=W%fname,status='old',form='unformatted',iostat=ios)
     ABI_CHECK(ios==0,'Opening old file: '//TRIM(W%fname))

     rdwr=5
     call scr_hdr_io(fform,rdwr,unt,comm,master,accesswff,Hscr)
     close(unt)
     if (W%prtvol>0) then ! Echo of the header
       call scr_hdr_io(fform,4,std_out,xmpi_self,master,accesswff,Hscr)
     end if
   end if
   !
   ! Communicate the header and copy basic parameters.
   call scrhdr_comm(Hscr,master,my_rank,comm)
   mat_type_read = Hscr%id
   nqlwl         = Hscr%nqlwl
   nomega        = Hscr%nomega
 end if

 W%nqlwl  = nqlwl
 W%nomega = nomega

 ABI_ALLOCATE(W%qlwl,(3,W%nqlwl))
 ABI_ALLOCATE(W%omega,(W%nomega))

 if (has_file) then
   W%qlwl  = Hscr%qlwl
   W%omega = Hscr%omega
   !
   ! G-vectors.
   W%npw=Hscr%npwe
   if (npw_asked>0) then
     if (npw_asked>Hscr%npwe) then
       write(msg,'(a,i8,2a,i8)')&
&       ' The number of G-vectors saved on file is less than the value required = ',npw_asked,ch10,&
&       ' Calculation will proceed with the Max available npw = ',Hscr%npwe
       MSG_WARNING(msg)
     else  
       W%npw=npw_asked ! Redefine the no. of G"s for W.
       write(msg,'(a,i8,2a,i8)')&
&       ' The Number of G-vectors saved on file is larger than the value required = ',npw_asked,ch10,&
&       ' Calculation will proceed with npw = ',W%npw
       MSG_COMMENT(msg)
     end if
   end if
   !
   ! consistency check on G-vectors and q-points.
   if ( ANY(Hscr%gvec(:,1:W%npw) /= Gsph%gvec(:,1:W%npw)) ) then
     !write(std_out) W%gvec, Gsph%gvec
     MSG_ERROR("Hscr%gvec /= Gsph%gvec(1:W%npw)")
   end if
   ABI_CHECK(Hscr%nqibz==Qmesh%nibz,"Mismatch in the number of q-points")
   ierr=0
   do iq_ibz=1,Hscr%nqibz
     if (ANY( ABS(Qmesh%ibz(:,iq_ibz) - Hscr%qibz(:,iq_ibz)) > tol6) ) then
       ierr=ierr+1
       write(std_out,'(i0,2(3f7.3,1x))')iq_ibz,Qmesh%ibz(:,iq_ibz),Hscr%qibz(:,iq_ibz)
     end if
   end do
   ABI_CHECK(ierr==0,"Wrong ordering in q-point list, aborting now")
   !
 end if
 !
 ABI_ALLOCATE(W%gvec,(3,W%npw))
 W%gvec = Gsph%gvec(:,1:W%npw)
 !
 ! Frequency mesh.
 W%nomega_r=1 
 W%nomega_i=0
 if (W%nomega==2) then
   W%nomega_r=1 
   W%nomega_i=1
 else 
   ! Real frequencies are packed in the first locations.
   W%nomega_r=1
   do iw=1,W%nomega
     if (DBLE(W%omega(iw))>0.001*Ha_eV) W%nomega_r=iw
   end do
   W%nomega_i=W%nomega-W%nomega_r
 end if     
 !
 ! ------------------------------ Initialization completed -------------------------------- 
 !
 ! Just to keep the code readable.
 npw    = W%npw
 nqibz  = W%nqibz
 nomega = W%nomega
 nI     = W%ni
 nJ     = W%nj
 ABI_ALLOCATE(W%Fgg,(nqibz))

 if (has_file) then  ! Read the ab-initio em1 from file.

   select case (mat_type_read)
   case (MAT_INV_EPSILON)
     call wrtout(std_out," Em1 will be initialized from SCR file: "//W%fname,"COLL")
     !
   case  (MAT_CHI0)  ! Write new SCR file.
     MSG_ERROR("Not coded yet")
     sus_fname = W%fname; scr_fname="TESTING_SUS2SCR"

     call sus2scr(sus_fname,scr_fname,W_info,Dtset,Cryst,Vcp,nspden,nfftf_tot,ngfftf,ae_rhor,prtvol,comm)
     W%fname = scr_fname  ! Change the name of the file associated to W.
     !
   case default
     write(msg,'(a,i0)')" Unsupported conversion from mat_type ",mat_type_read
     MSG_ERROR(msg)
   end select
   !
   !
   do iq_ibz=1,nqibz
     if (.not.W%keep_q(iq_ibz)) CYCLE 
     !
     nqlwl=0; is_qeq0=(normv(W%qibz(:,iq_ibz),Cryst%gmet,'G')<GW_TOLQ0)
     if (is_qeq0) nqlwl=W%nqlwl
     !
     ! Allocate F_{GG'}(w)
     call fgg_init(W%Fgg(iq_ibz),npw,nomega,nqlwl)
     !
     ! Read data from file.
     call read_screening(W%fname,npw,1,nomega,W%Fgg(iq_ibz)%mat,accesswff,comm,iqiA=iq_ibz)
     W%Fgg(iq_ibz)%has_mat = MAT_STORED
     !
     ! Copy heads and wings if q-->0.
     if (nqlwl>0) then 
       W%Fgg(iq_ibz)%lwing=Hscr%lwing(1:npw,:,:)
       W%Fgg(iq_ibz)%uwing=Hscr%uwing(1:npw,:,:)
       W%Fgg(iq_ibz)%head=HUGE(one)    ! TODO
       W%Fgg(iq_ibz)%has_wings = MAT_STORED
     end if
     ! W contains Em1 and is ready to use.
   end do

 else  
   !
   ! Model dielectric function. Only epsm-1 is supported here.
   call wrtout(std_out," Calculating model dielectric function... ","COLL")
   ABI_CHECK(W%nomega==1,"Cannot use nomega>1 in model dielectric function")
   !
   do iq_ibz=1,nqibz
     !
     if (.not.W%keep_q(iq_ibz)) CYCLE
     !
     ! The wings are not used here.
     nqlwl=0; is_qeq0= (normv(W%qibz(:,iq_ibz),Cryst%gmet,'G')<GW_TOLQ0)
     !
     ! Calculate the model. Note that mdielf awaits an index in the BZ.
     found = has_bz_item(Qmesh,Qmesh%ibz(:,iq_ibz),iq_bz,g0)
     if (.not.found.or.ANY(g0/=0)) then
       MSG_ERROR("Problem in retrieving ibz point")
     end if
     !
     ! Allocate F_{GG'}(w).
     call fgg_init(W%Fgg(iq_ibz),npw,nomega,nqlwl)

     eps_inf  =  W%Info%eps_inf
     mdf_type =  W%Info%use_mdf
     em1_ggw  => W%Fgg(iq_ibz)%mat

     ! Construct W TODO check the new implementation.
     call screen_mdielf(iq_bz,npw,nomega,mdf_type,eps_inf,Cryst,Qmesh,Vcp,Gsph,&
&      nspden,nfftf_tot,ngfftf,ae_rhor,"EM1",em1_ggw,comm)

     W%Fgg(iq_ibz)%has_mat = MAT_STORED

     do iw=1,nomega
       write(msg,'(a,i3,a,i4,a)')'  Model symmetrical e^{-1} (q=',iq_ibz,', omega=',iw,', G,G'')'
       call wrtout(std_out,msg,'COLL')
       call print_arr(em1_ggw(:,:,iw))
     end do
     !
   end do ! iq_ibz
   !
 end if
 !
 ! Init plasmon-pole parameters.
 if (W%has_ppmodel>0) then
   MSG_WARNING("Calculating PPmodel parameters")
   ppmodel     = W%Info%use_ppm
   drude_plsmf = W%Info%drude_plsmf 
   call ppm_init(W%PPm,W%mqmem,W%nqibz,W%npw,ppmodel,drude_plsmf)
   !
   do iq_ibz=1,nqibz
     if (screen_ihave_fgg(W,iq_ibz,how="Stored")) then
       em1_ggw => W%Fgg(iq_ibz)%mat
       call new_setup_ppmodel(W%PPm,iq_ibz,Cryst,Qmesh,npw,nomega,W%omega,em1_ggw,nfftf_tot,Gsph%gvec,ngfftf,W%ae_rhor(:,1))
     end if
   end do
   !
 end if
 !
 ! Deallocate Fgg if the matrices are not needed anymore.
 if (deallocate_Fgg) then
   call screen_fgg_qbz_set(W,0,0,"Pointer") ! Avoid dangling pointer.
   call fgg_free(W%Fgg,keep_q=W%keep_q)
 end if

 if (has_file) call free_scrhdr(Hscr)

 DBG_EXIT("COLL")

end subroutine screen_init
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_symmetrizer
!! NAME
!!  screen_symmetrizer
!!
!! FUNCTION
!!  Modify the status of the object so that the symmetrized component F(q_bz)_GG' is calculated
!!  (if needed) and is made available in the internal buffer. This routine must be called before 
!!  performing any operation that involves the symmetrized component of the two-point function.
!!
!! INPUTS
!!  W<screen_t>=Data structure used to represent the two-point function
!!  Cryst<crystal_structure>=Info on the crystal structure.
!!  Gsph<Gvectors_type>=data related to the G-sphere
!!  Qmesh<BZ_mesh_type>=Structure defining the q-mesh used for sample the BZ.
!!  iq_bz=Index of the q-point in the BZ where F(q_bz)_GG' is wanted. 
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   W%PPm
!!   W%Fgg_qbz
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the 
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere 
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required 
!!  to reconstruct the BZ.
!! 
!! PARENTS
!!      exc_build_block
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine screen_symmetrizer(W,iq_bz,Cryst,Gsph,Qmesh) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_symmetrizer'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz
 type(screen_t),intent(inout) :: W
 type(crystal_structure),intent(in) :: Cryst
 type(gvectors_type),intent(in) :: Gsph
 type(bz_mesh_type),intent(in) :: Qmesh
!arrays

!Local variables-------------------------------
!scalars
 integer,parameter :: nqlwl0=0
 integer :: iq_ibz,isym_q,itim_q,npw,nomega,accesswff,nqibz !istat
 integer :: mdf_type
 real(dp) :: eps_inf
 logical :: q_isirred
 !character(len=500) :: msg
!arrays
 !complex(gwpc),intent(out) :: epsm1_qbz(npw_r,npw_r,nomega)  
 real(dp) :: qbz(3)
 complex(gwpc),pointer :: em1_qbz(:,:,:),em1_qibz(:,:,:)
 !complex(gwpc),allocatable :: cbuff(npweA,npweA,nomegaA,nqibzA)

! *********************************************************************

 ! Store the index of the q-point in the BZ for checking purpose.
 W%fgg_qbz_idx = iq_bz

 npw       = W%npw
 nqibz     = W%nqibz
 nomega    = W%nomega
 accesswff = W%accesswff

 call get_bz_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q,isirred=q_isirred)
 !
 ! ========================================================
 ! ==== Branching for in-core or out-of-core solutions ====
 ! ========================================================
 if (screen_ihave_fgg(W,iq_ibz,how="Stored")) then
   !
   if (q_isirred) then ! Symmetrization is not needed. Point the data in memory.
     call screen_fgg_qbz_set(W,iq_bz,nqlwl0,"Pointer") 
     W%Fgg_qbz => W%Fgg(iq_ibz)  
   else 
     ! Allocate space. ! TODO Wings are not symmetrized but oh well
     call screen_fgg_qbz_set(W,iq_bz,nqlwl0,"Allocate")  ! Dimensions should not be changed.
     
     em1_qibz => W%Fgg(iq_ibz)%mat ! out-of-place symmetrization.
     em1_qbz  => W%Fgg_qbz%mat  

     call em1_symmetrize_op(iq_bz,npw,nomega,Gsph,Qmesh,em1_qibz,em1_qbz) 
   end if
   !
   if (W%has_ppmodel>0) then ! Symmetrize the ppmodel tables.
     em1_qibz => W%Fgg(iq_ibz)%mat
     call ppm_symmetrizer(W%PPm,iq_bz,Cryst,Qmesh,Gsph,npw,nomega,W%omega,em1_qibz,W%nfftf_tot,W%ngfftf,W%ae_rhor(:,1))
   end if
   !
 else if (screen_ihave_fgg(W,iq_ibz,how="Allocated")) then
   MSG_ERROR("Fgg_iqibz is allocated but not initialized!")
   !
 else  
   MSG_COMMENT("out of core with file "//TRIM(W%fname))
   MSG_ERROR("Not implemented error")
   !
   ! Ppmodel calculations with ppm tables in memory.
   ! TODO treat the case in which IBZ tables are stored in memory.
   if (W%has_ppmodel>0) then
     !%call ppm_symmetrizer(W%PPm,iq_bz,Cryst,Qmesh,Gsph,npw,nomega,W%omega,em1_qibz,W%nfftf_tot,W%ngfftf,W%ae_rhor(:,1))
   end if
   !
   ! Allocate the BZ buffer.
   call screen_fgg_qbz_set(W,iq_bz,nqlwl0,"Allocate") 
   !
   ! Read data from file.
   em1_qbz => W%Fgg_qbz%mat 

   call read_screening(W%fname,npw,1,nomega,em1_qbz,accesswff,xmpi_self,iqiA=iq_ibz)
   W%Fgg_qbz%has_mat = MAT_STORED

   ! Symmetrize the ppmodel using em1_qibz.
   if (W%has_ppmodel>0) then
     call ppm_symmetrizer(W%PPm,iq_bz,Cryst,Qmesh,Gsph,npw,nomega,W%omega,em1_qbz,W%nfftf_tot,W%ngfftf,W%ae_rhor(:,1))
   end if
   !
   ! In-place symmetrization.
   call em1_symmetrize_ip(iq_bz,npw,nomega,Gsph,Qmesh,em1_qbz) 
 end if

 ! Calculate model dielectric function for this q-point in the BZ.
 eps_inf  =  W%Info%eps_inf
 mdf_type =  W%Info%use_mdf

!%  call screen_mdielf(iq_bz,npw,nomega,mdf_type,eps_inf,Cryst,Qmesh,Vcp,Gsph,&
!% &  W%nspden,W%nfftf_tot,W%ngfftf,W%ae_rhor,"EM1",em1_qbz,xmpi_self)

end subroutine screen_symmetrizer
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_w0gemv
!! NAME
!! screen_w0gemv
!!
!! FUNCTION
!!
!! INPUTS
!!  W<screen_t>=
!!  in_npw=Number of G vectors in in_ket
!!  nspinor=Number of spinorial components.
!!  in_ket(in_npw)= |\phi> in reciprocal space. 
!!  trans= On entry, TRANS specifies the operation to be performed as follows:
!!  TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!! 
!!  TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!! 
!!  TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.
!!
!! OUTPUT
!!   out_ket(in_npw)= W |\phi\> in reciprocal space.
!!   ZGEMV  performs one of the matrix-vector operations
!!   *
!!   *     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
!!   *
!!   *     y := alpha*A**H*x + beta*y,
!!   *
!!   *  where alpha and beta are scalars, x and y are vectors and A is an m by n matrix.
!!
!! PARENTS
!!      exc_build_block
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine screen_w0gemv(W,trans,in_npw,nspinor,only_diago,alpha,beta,in_ket,out_ket)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_w0gemv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: in_npw,nspinor
 complex(gwpc),intent(in) :: alpha,beta
 logical,intent(in) :: only_diago
 character(len=*),intent(in) ::  trans
 type(screen_t),intent(in) :: W
!arrays
 complex(gwpc),intent(in) :: in_ket(in_npw*nspinor)
 complex(gwpc),intent(out) :: out_ket(in_npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: ig,lda
 !character(len=500) :: msg
!arrays
 complex(gwpc),pointer :: em1_qbz(:,:)

! *************************************************************************

 lda = W%npw
 em1_qbz => W%Fgg_qbz%mat(:,:,1)

 if (.not.only_diago) then
   !out_ket = MATMUL(TRANSPOSE(CONJG(em1_qbz(1:in_npw,1:in_npw)),in_ket)
   call xgemv(trans,in_npw,in_npw,alpha,em1_qbz,lda,in_ket,1,beta,out_ket,1)
   !
 else 
   !
   if (beta /= czero_gw) then
     !
     if ( starts_with(trans,(/"C"/)) ) then
       do ig=1,in_npw
         out_ket(ig) = alpha * CONJG(em1_qbz(ig,ig)) * in_ket(ig) + beta * out_ket(ig)
       end do
     else if ( starts_with(trans,(/"N","T"/)) ) then
       do ig=1,in_npw
         out_ket(ig) = alpha * em1_qbz(ig,ig) * in_ket(ig) + beta * out_ket(ig)
       end do
     else 
       MSG_ERROR("Wrong trans: "//trans)
     end if
     !
   else  ! beta==0
     !
     if ( starts_with(trans,(/"C"/)) ) then
       do ig=1,in_npw
         out_ket(ig) = alpha * CONJG(em1_qbz(ig,ig)) * in_ket(ig) 
       end do
     else if ( starts_with(trans,(/"N","T"/)) ) then
       do ig=1,in_npw
         out_ket(ig) = alpha * em1_qbz(ig,ig) * in_ket(ig) 
       end do
     else 
       MSG_ERROR("Wrong trans: "//trans)
     end if
     !
   end if
   !
 end if

end subroutine screen_w0gemv
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_times_ket
!! NAME
!! screen_times_ket
!!
!! FUNCTION
!!
!! INPUTS
!!  W<screen_t>=Data structure describing the two-point function in reciprocal space. See also SIDE EFFECTS.
!!  nomega=Total number of frequencies where $\Sigma_c$ matrix elements are evaluated.
!!  npwc=Number of G vectors for the correlation part.
!!  npwx=Number of G vectors in rhotwgp for each spinorial component.
!!  nspinor=Number of spinorial components.
!!  theta_mu_minus_e0i=1 if e0i is occupied, 0 otherwise. Fractional occupancy in case of metals. 
!!  omegame0i(nomega)=Contains $\omega-\epsilon_{k-q,b1,\sigma}$
!!  zcut=Small imaginary part to avoid the divergence in the ppmodel. (see related input variable)
!!  rhotwgp(npwx*nspinor)=Matrix elements: $<k-q,b1,\sigma|e^{-i(q+G)r} |k,b2,\sigma>*vc_sqrt$
!!
!! OUTPUT
!! ket(npwc,nomega)=Contains \Sigma_c(\omega)|\phi> in reciprocal space. 
!!
!! SIDE EFFECTS
!! npoles_missing=Incremented with the number of poles whose contribution has not been taken into account due to
!!  limited frequency mesh used for W.
!!
!! OUTPUT
!!   ket(npwc,nomega)=Contains \Sigma_c(\omega)|\phi> in reciprocal space. 
!!   sigcme(nomega) (to be described), only relevant if ppm3 or ppm4
!!
!! PARENTS
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine screen_times_ket(W,npwc,npwx,nspinor,nomega,rhotwgp,omegame0i,zcut,theta_mu_minus_e0i,ket,sigcme,npoles_missing)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_times_ket'
 use interfaces_70_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,npwx,nspinor
 integer,intent(inout) :: npoles_missing
 real(dp),intent(in) :: theta_mu_minus_e0i,zcut
 type(screen_t),intent(in) :: W
!arrays
 real(dp),intent(in) :: omegame0i(nomega)
 !complex(gwpc),intent(in) :: epsm1q(npwc,npwc,nomegae) 
 complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 complex(gwpc),intent(inout) :: ket(nspinor*npwc,nomega)
 complex(gwpc),intent(out) :: sigcme(nomega)

!Local variables-------------------------------
!scalars
 integer :: nomegae,nomegaei,nomegaer
 character(len=500) :: msg
!arrays
 complex(dpc),pointer :: w_omega(:)
 complex(gwpc),pointer :: epsm1_qbz(:,:,:)
! complex(gwpc),pointer :: botsq(:,:),eig(:,:),otq(:,:)
 !botsq(npwc,PPm%dm2_botsq),eig(PPm%dm_eig,PPm%dm_eig),otq(npwc,PPm%dm2_otq)

! *************************************************************************

 ABI_UNUSED(zcut)

 ! Mesh for W
 nomegae  = W%nomega
 nomegaei = W%nomega_i
 nomegaer = W%nomega_r
 w_omega => W%omega

 SELECT CASE (W%Info%wint_method) 
 CASE (WINT_CONTOUR) ! Contour deformation method.
   !
   ! Pass the symmetrized matrix.
   epsm1_qbz => W%Fgg_qbz%mat(1:npwc,1:npwc,1:nomegae) 
   !
   call calc_sigc_cd(npwc,npwx,nspinor,nomega,nomegae,nomegaer,nomegaei,rhotwgp,&
&     w_omega,epsm1_qbz,omegame0i,theta_mu_minus_e0i,ket,npoles_missing)

 CASE (WINT_AC)
   MSG_ERROR("Not implemented error")

 CASE (WINT_PPMODEL)
   MSG_ERROR("Not implemented error")
!   botsq => W%PPm%bigomegatwsq_qbz%value
!   otq   => W%PPm%omegatw_qbz%value
!   eig   => W%PPm%eigpot_qbz%value
!
!   call calc_sig_ppm(W%PPm,nspinor,npwc,nomega,rhotwgp,botsq,otq,&
!&     omegame0i,zcut,theta_mu_minus_e0i,eig,npwx,ket,sigcme)

!%    call ppm_times_ket(W%PPm,nspinor,npwc,nomega,rhotwgp,omegame0i,zcut,theta_mu_minus_e0i,npwx,ket,sigcme)

 CASE DEFAULT
   write(msg,'(a,i0)')" unknown value for wint_method ",W%Info%wint_method
   MSG_ERROR(msg)
 END SELECT 

end subroutine screen_times_ket
!!***

!----------------------------------------------------------------------

!!****f* m_screen/sus2scr
!! NAME
!! sus2scr
!!
!! FUNCTION
!!  This routine calculate the inverse of the symmetrical dielectric matrix starting
!!  from the irreducible polarizability stored on the external file sus_fname.
!!
!! INPUTS
!!  sus_fname=Name of the suscetibility file.
!!  scr_fname=Name of the output file that will containe the inverse of the symmetrized dielectric matrix.
!!  Dtset
!!  W_Info<screen_info_t>=Container storing the parameters defining the kind of symmetrized dielectric matrix that is wanted.
!!  Cryst<crystal_structure>=Info on the crystal structure.
!!  nspden=Number of spin Density components.
!!  Vcp<vcoul_t>=Structure gathering data on the Coulomb interaction including a possible cutoff.
!!    %nqibz=Number of q-points.
!!    %qibz(3,nqibz)=q-points in the IBZ.  
!!  nfftf_tot=Total number of points in the FFT mesh.
!!  ngfftf(18)=Info on the FFT mesh.
!!  ae_rhor(nfftf_tot,nspden)= Density on the ful real space FFT box used to construct the TDDFT kernel.
!!    NOTE that ae_rhor is given on the dense mesh since if has the PAW onsite contributions
!!  prtvol=Verbosity level.
!!  comm=MPI communicator.
!!
!! SIDE EFFECTS
!!  Output is written on file.
!!
!! PARENTS
!!      m_screen
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine sus2scr(sus_fname,scr_fname,W_info,Dtset,Cryst,Vcp,nspden,nfftf_tot,ngfftf,ae_rhor,prtvol,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sus2scr'
 use interfaces_14_hidewrite
 use interfaces_56_xc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: nspden,nfftf_tot,comm,prtvol
 character(len=fnlen),intent(in) :: sus_fname,scr_fname
 type(screen_info_t),intent(in) :: W_info
 type(vcoul_t),intent(in) :: Vcp
 type(dataset_type),intent(in) :: Dtset
 type(crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: ae_rhor(nfftf_tot,nspden)

!Local variables-------------------------------
!scalars
 integer :: approx_type,option_test,ixc_required
 integer :: nomega,npw,dim_wing,nqibz,iq_ibz,accesswff
 integer :: sus_unt,scr_unt,ios,rdwr,fform,nkxc,nI,nJ
 integer :: istat,master,my_rank,nprocs,comm_self !ierr,
 logical :: is_qeq0
 character(len=500) :: msg
 type(ScrHdr_type) :: Hscr,Hscr_cp
 type(spectra_type) :: Spectra
!arrays
 complex(dpc),allocatable :: chi0_lwing(:,:,:),chi0_uwing(:,:,:),chi0_head(:,:,:)
 complex(gwpc),allocatable :: chi0(:,:,:)
 complex(gwpc),allocatable :: fxc_ADA_qibz(:,:)
 !complex(gwpc),intent(in),optional :: fxc_ADA(inpw*nI,inpw*nJ,)
 complex(gwpc),allocatable :: kxcg(:,:)
 !complex(gwpc),intent(in) :: kxcg(nfftf_tot,nkxc)

! *************************************************************************

 DBG_ENTER("COLL")

 nprocs  = xcomm_size(comm)
 my_rank = xcomm_rank(comm)
 master=0
 ABI_CHECK(nprocs==1,"MPI not coded yet")

 accesswff = IO_MODE_FORTRAN

 ! MG TODO We use comm_self for the inversion as the single precision version is not yet available
 comm_self = xmpi_self

 approx_type  = W_info%vtx_family
 option_test  = W_info%vtx_test
 ixc_required = W_info%ixc
 !
 ! === Open the SUS file ===
 !if (my_rank==master) then
   call wrtout(std_out,' Testing SUS file: '//TRIM(sus_fname),'COLL')
   sus_unt=get_unit()
   open(unit=sus_unt,file=sus_fname,status='old',form='unformatted',iostat=ios)
   ABI_CHECK(ios==0,' Opening old file '//TRIM(sus_fname))
 !end if
                                                                        
 rdwr=5
 call scr_hdr_io(fform,rdwr,sus_unt,comm,master,accesswff,Hscr)
 close(sus_unt)

 if (my_rank==master.and.prtvol>0) then ! Echo of the header
   call scr_hdr_io(fform,4,std_out,comm,master,accesswff,Hscr)
 end if
 !
 if (Hscr%id/=MAT_CHI0) then ! Make sure this is a SUSC file.
   write(msg,'(a,i0)')" Expecting a SUS file but found Hdr%id= ",Hscr%id
   MSG_ERROR(msg)
 end if
 !
 ! Get dimension and basic parameters from the input header.
 npw    = Hscr%npwe  
 nqibz  = Hscr%nqibz
 nomega = Hscr%nomega
 nI     = Hscr%nI
 nJ     = Hscr%nJ

 if (nI/=1.or.nJ/=1) then 
   MSG_ERROR("nI or nJ=/1 not implemented yet")
 end if

 if ( Vcp%nqibz/=Hscr%nqibz) then
   write(msg,'(a,2i0)')" Vcp%nqibz/=Hscr%nqibz wrong dimensions passed by the caller ",Vcp%nqibz, Hscr%nqibz
   MSG_ERROR(msg)
 end if

 if ( ANY(ABS(Vcp%qibz-Hscr%qibz)>tol6) ) then
   write(std_out,*)Vcp%qibz 
   write(std_out,*)Hscr%qibz
   MSG_ERROR("Wrong list of q-points")
 end if
 !
 ! Open the output file.
 call wrtout(std_out," Writing epsilon^-1 on file "//TRIM(scr_fname),'COLL')
 scr_unt=get_unit()
 open(unit=scr_unt,file=scr_fname,status='new',form='unformatted',iostat=ios)
 ABI_CHECK(ios==0," Opening new file "//TRIM(scr_fname))
 !
 ! Update the entries in the new header (basic dimension are unchanged).
 ! TODO, write function to return title, just for info
 ! TODO nullification of pointers in fortran structures is needed to avoid problems during the copy.
 call copy_scrhdr(Hscr,Hscr_cp)

 Hscr_cp%id        = MAT_INV_EPSILON
 Hscr_cp%ikxc      = ixc_required
 Hscr_cp%test_type = option_test 
 Hscr_cp%title(1)  = 'SCR file: epsilon^-1'
 Hscr_cp%title(2)  = 'Unknown'
 if (option_test==VTX_TEST_PARTICLE) Hscr_cp%title(2)  = 'TESTPARTICLE'
 if (option_test==VTX_TEST_CHARGE  ) Hscr_cp%title(2)  = 'TESTCHARGE'
 !
 rdwr=2; fform=Hscr_cp%fform
 call scr_hdr_io(fform,rdwr,scr_unt,comm_self,master,accesswff,Hscr_cp)
 call free_scrhdr(Hscr_cp)
 !
 ! Allocate the matrix for a single q-point.
 ! TODO: Rewrite make_epsm1_driver so that the inversion is performed for a single frequency.
 ABI_ALLOCATE(chi0,(npw*nI,npw*nJ,nomega))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in chi0")
 !
 ! Build the symmetrized e^{-1}(q) then writing the results on file.
 do iq_ibz=1,nqibz
   !
   dim_wing=0; is_qeq0 = (normv(Hscr%qibz(:,iq_ibz),Cryst%gmet,'G')<GW_TOLQ0)
   if (is_qeq0) dim_wing=Hscr%nqlwl
   !if (is_qeq0) dim_wing=3 
   !
   ! Read chi0_q from file.
   call read_screening(sus_fname,npw,1,nomega,chi0,accesswff,comm,iqiA=iq_ibz)
   !
   ! Copy heads and wings.
   ABI_ALLOCATE(chi0_lwing,(npw*nI,nomega,dim_wing))
   ABI_ALLOCATE(chi0_uwing,(npw*nJ,nomega,dim_wing))
   ABI_ALLOCATE(chi0_head,(dim_wing,dim_wing,nomega))

   if (dim_wing>0) then 
     chi0_lwing=czero
     chi0_uwing=czero
     chi0_head =czero
     ! TODO
     !chi0_lwing=Hscr%lwing(1:npw,:,:)
     !chi0_uwing=Hscr%uwing(1:npw,:,:)
     !chi0_head =Hscr%uwing(1:npw,:,:)
     MSG_WARNING("Treatment of wings is not coded yet")
   end if

   select case (approx_type)
   case (VTX_FAMILY_NONE, VTX_FAMILY_TDDFT)
     MSG_COMMENT('Entering Kxc branch')
     !
     if (approx_type==VTX_FAMILY_TDDFT) then
       nkxc=1 ! LDA only
       ABI_ALLOCATE(kxcg,(nfftf_tot,nkxc))
       call xc_kernel(Dtset,Cryst,ixc_required,ngfftf,nfftf_tot,nspden,ae_rhor,npw,nkxc,kxcg,Hscr%gvec,comm)
       MSG_ERROR("Not tested")
     else 
       nkxc=0
       ABI_ALLOCATE(kxcg,(nfftf_tot,nkxc))
       istat = ABI_ALLOC_STAT
     end if

     call make_epsm1_driver(iq_ibz,dim_wing,npw,nI,nJ,nomega,Hscr%omega,&
&              approx_type,option_test,Vcp,nfftf_tot,ngfftf,nkxc,kxcg,Hscr%gvec,chi0_head,&
&              chi0_lwing,chi0_uwing,chi0,Spectra,comm)

     ABI_DEALLOCATE(kxcg)
  
   case (VTX_FAMILY_ADA)
     MSG_WARNING('Entering in-core fxc_ADA branch')
     ABI_ALLOCATE(fxc_ADA_qibz,(npw*nI,npw*nJ))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,"out of memory in fxc_ADA_qibz")
     MSG_ERROR("Not coded")

     !call xc_kernel_ADA(Dtset,ixc,ngfftf,nfft_tot,nspden,ae_rhor,Cryst%rprimd,&
     !&  npw,nqibz,qibz,fxc_ADA,gvec,Cryst%gprimd,comm,kappa_init=kappa_init)

     call make_epsm1_driver(iq_ibz,dim_wing,npw,nI,nJ,nomega,Hscr%omega,&
&              approx_type,option_test,Vcp,nfftf_tot,ngfftf,nkxc,kxcg,Hscr%gvec,chi0_head,&
&              chi0_lwing,chi0_uwing,chi0,Spectra,comm,fxc_ADA=fxc_ADA_qibz)

     ABI_DEALLOCATE(fxc_ADA_qibz)

   case default
     write(msg,'(a,i0)')" Unknown value for approx_type: ",approx_type
     MSG_ERROR(msg)
   end select

   ABI_DEALLOCATE(chi0_lwing)
   ABI_DEALLOCATE(chi0_uwing)
   ABI_DEALLOCATE(chi0_head)

   if (is_qeq0) then
     !% call repr_dielconst(Spectra,msg)
     !% call wrtout(std_out,msg,'PERS') 
     !% call wrtout(ab_out,msg,'PERS')
     !% call dump_spectra(Spectra,write_bits,spectra_fname)
   end if
   call destroy_spectra(Spectra)

   call write_screening(scr_unt,accesswff,npw,nomega,chi0)

 end do ! iq_ibz

 ABI_DEALLOCATE(chi0)
 call free_scrhdr(Hscr)

 close(scr_unt)

 call xbarrier_mpi(comm)

 DBG_EXIT("COLL")

end subroutine sus2scr
!!***

!----------------------------------------------------------------------

!!****f* m_screen/em1_symmetrize_ip
!! NAME
!!  em1_symmetrize_ip
!!
!! FUNCTION
!!  Symmetrizes the two-point function in G-space. Symmetrization is done 
!!  inplace thorugh an auxiliary work array of dimension (npwc,npwc) 
!!
!! INPUTS
!!  nomega=All frequencies from 1 up to nomega are symmetrized.
!!  npwc=Number of G vectors in the symmetrized matrix.
!!  Gsph<Gvectors_type>=data related to the G-sphere
!!  Qmesh<BZ_mesh_type>=Structure defining the q-mesh used for Er.
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required. 
!!
!! SIDE EFFECTS 
!!  epsm1(npwc,npwc,nomega) 
!!   input:  filled with the matrix at the q-point that has to be symmetrized.
!!   output: symmetrised matrix.
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the 
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere 
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required 
!!  to reconstruct the BZ.
!! 
!!  * Remember the symmetry properties of E
!!    If q_bz=Sq_ibz+G0:
!! 
!!    $ E_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau} E_{G1,G2)}(q)
!!
!!    The invariance under exchange of the real space position E(1,2) = E(2,1) leads to:
!!    $ E_{-G2,-G1}(-q) = E_{G1,G2)
!!
!! PARENTS
!!      m_screen
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine em1_symmetrize_ip(iq_bz,npwc,nomega,Gsph,Qmesh,epsm1) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'em1_symmetrize_ip'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npwc
 type(Gvectors_type),intent(in) :: Gsph
 type(BZ_mesh_type),intent(in) :: Qmesh
!arrays
 complex(gwpc),intent(inout) :: epsm1(npwc,npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: iw,g1,g2,isg1,isg2,iq_ibz,itim_q,isym_q,istat
 logical :: q_isirred
 character(len=500) :: msg
!arrays
 integer,pointer :: grottb(:)
 real(dp) :: qbz(3)
 complex(gwpc),pointer :: phmSgt(:)
 complex(gwpc),allocatable :: work(:,:)

! *********************************************************************

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q,isirred=q_isirred)

 if (q_isirred) RETURN ! Nothing to do

 grottb => Gsph%rottb (:,itim_q,isym_q)
 phmSgt => Gsph%phmSGt(:,isym_q) 

 ABI_ALLOCATE(work,(npwc,npwc))
 istat = ABI_ALLOC_STAT
 if (istat/=0) then
   write(msg,'(a,f8.2,a)')" out of memory in work , requiring ",npwc**2*gwpc*b2Mb," Mb"
   MSG_ERROR(msg)
 end if

 do iw=1,nomega
   do g2=1,npwc
     isg2 = grottb(g2)
     do g1=1,npwc
       isg1 = grottb(g1)
       work(isg1,isg2) = epsm1(g1,g2,iw) * phmSgt(g1) * CONJG(phmSgt(g2))
     end do
   end do
   epsm1(:,:,iw) = work
 end do

 ABI_DEALLOCATE(work)
 !
 ! Account for time-reversal ----------------------
 if (itim_q==2) then
   do iw=1,nomega
     epsm1(:,:,iw)=TRANSPOSE(epsm1(:,:,iw))
   end do
 end if

end subroutine em1_symmetrize_ip
!!***

!----------------------------------------------------------------------

!!****f* m_screen/em1_symmetrize_op
!! NAME
!!  em1_symmetrize_op
!!
!! FUNCTION
!!  Symmetrizes the two-point function in G-space. Symmetrization is done outofplace.
!!
!! INPUTS
!!  nomega=All frequencies from 1 up to nomega are symmetrized.
!!  npwc=Number of G vectors in the symmetrized matrix.
!!  Gsph<Gvectors_type>=data related to the G-sphere
!!  Qmesh<BZ_mesh_type>=Structure defining the q-mesh used for Er.
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required. 
!!  in_epsm1(npwc,npwc,nomega) 
!!
!! OUTPUT
!!  out_epsm1(npwc,npwc,nomega) 
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the 
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere 
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required 
!!  to reconstruct the BZ.
!! 
!!  * Remember the symmetry properties of E
!!    If q_bz=Sq_ibz+G0:
!! 
!!    $ E_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau} E_{G1,G2)}(q)
!!
!!    The invariance under exchange of the real space position E(1,2) = E(2,1) leads to:
!!    $ E_{-G2,-G1}(-q) = E_{G1,G2)
!!
!! PARENTS
!!      m_screen
!!
!! CHILDREN
!!      get_bz_item
!!
!! SOURCE

subroutine em1_symmetrize_op(iq_bz,npwc,nomega,Gsph,Qmesh,in_epsm1,out_epsm1) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'em1_symmetrize_op'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npwc
 type(Gvectors_type),intent(in) :: Gsph
 type(BZ_mesh_type),intent(in) :: Qmesh
!arrays
 complex(gwpc),intent(in) :: in_epsm1(npwc,npwc,nomega)
 complex(gwpc),intent(out) :: out_epsm1(npwc,npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: iw,g1,g2,isg1,isg2,iq_ibz,itim_q,isym_q !,istat
 logical :: q_isirred
!arrays
 integer,pointer :: grottb(:)
 real(dp) :: qbz(3)
 complex(gwpc),pointer :: phmSgt(:)

! *********************************************************************

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q,isirred=q_isirred)

 if (q_isirred) then
   out_epsm1 = in_epsm1
   RETURN 
 end if

 grottb => Gsph%rottb (:,itim_q,isym_q)
 phmSgt => Gsph%phmSGt(:,isym_q) 

 do iw=1,nomega
   do g2=1,npwc
     isg2 = grottb(g2)
     do g1=1,npwc
       isg1 = grottb(g1)
       out_epsm1(isg1,isg2,iw) = in_epsm1(g1,g2,iw) * phmSgt(g1) * CONJG(phmSgt(g2))
     end do
   end do
 end do
 !
 ! Account for time-reversal ----------------------
 if (itim_q==2) then
   do iw=1,nomega
     out_epsm1(:,:,iw)=TRANSPOSE(out_epsm1(:,:,iw))
   end do
 end if

end subroutine em1_symmetrize_op
!!***

!----------------------------------------------------------------------

END MODULE m_screen
!!***
