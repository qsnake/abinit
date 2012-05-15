!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_bz_mesh
!! NAME
!!  m_bz_mesh
!!
!! FUNCTION
!!  This module provides the definition of the BZ_mesh_type structure gathering information 
!!  on the sampling of the Brillouin zone. It also contains useful tools to operate on k-points.
!!  and the definition of the Little_group data type. The Little_group structure is used 
!!  to store tables and useful info on the set of k-points belonging
!!  to the irreducible wedge defined by the symmetry properties 
!!  of the point group that preserve the external q-point.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG, GMR, VO, LR, RWG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  One has to use a fixed ordering of the loops over nsym and time-reversal 
!!  when the full zone is reconstructed by symmetry starting from the IBZ. 
!!  This is especially important in systems with both time-reversal and 
!!  spatial inversion as the effect of the two operation in reciprocal 
!!  space is very similar the only difference being the possibly non-zero
!!  fractional translation associated to the spatial inversion.
!!  In the present implementation, the spatial inversion has the precedence 
!!  wrt time-reversal (i.e., do itim; do isym). 
!!  Note that this particular ordering should be used in any routine used to 
!!  symmetrize k-dependent quantities in the full BZ zone to avoid possible errors.
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

MODULE m_bz_mesh

 use m_profiling

 use defs_basis
 use m_errors 

 use m_io_tools,       only : flush_unit
 use m_numeric_tools,  only : is_zero, isinteger, imin_loc, imax_loc, bisect
 use m_geometry,       only : normv
 use m_crystal,        only : crystal_structure

 use m_tetrahedron

 implicit none

 private  

 real(dp),parameter :: TOL_KDIFF=0.0001_dp
 ! Tolerance below which two points are considered equal within a RL vector:
 ! for each reduced direction the absolute difference between the coordinates must be less that TOL_KDIFF

 integer,parameter :: NONE_KPTRLATT(3,3)=RESHAPE((/0,0,0,0,0,0,0,0,0/),(/3,3/))
!!***

!!****t* m_bz_mesh/bz_mesh_type
!! NAME
!! bz_mesh_type
!!
!! FUNCTION
!! The BZ_mesh structured datatype contains different information on the grid used to sample the BZ : 
!! the k-points in the full Brillouin zone BZ, the irreducible wedge IBZ as well as tables describing 
!! the symmetry relationship between the points.
!!
!! SOURCE

 type,public :: bz_mesh_type

  !scalars
  integer :: nshift=0

  integer :: nbz            
  ! Number of points in the BZ.

  integer :: nibz           
  ! Number of points in the IBZ.

  integer :: nsym           
  ! Number of symmetry operations.

  integer :: kptopt
  ! kptopt=option for the generation of k points (see input variable description)
  ! 
  ! 1  if both time-reversal and point group symmetries are used.
  ! 2  if only time-reversal symmetry is used.
  ! 3  do not take into account any symmetry (except the identity).
  ! 4  if time-reversal is not used (spin-orbit coupling).
  ! <0 number of segments used to construct the k-path for NSCF calculation.

  integer :: timrev         
  ! 2 if time reversal symmetry can be used, 1 otherwise.

 !arrays
  integer :: kptrlatt(3,3) = NONE_KPTRLATT
   ! Coordinates of three vectors in real space, expressed in reduced coordinates.
   ! They define a super-lattice in real space. The k point lattice is the reciprocal of
   ! this super-lattice, eventually shifted by shift.
   ! Not available if the structure is initialized from the points in the IBZ.

  integer,pointer :: rottb(:,:,:)    SET2NULL      
  ! rottb(nbz,timrev,nsym),
  ! Index of (IS)k in the BZ array where S is a sym operation in reciprocal space,
  ! I is the identity or the inversion operator (1,2 resp)

  integer,pointer :: rottbm1(:,:,:)  SET2NULL
  ! rottbm1(nbz,timrev,nsym)
  ! Index of IS^{-1} k in the BZ array.

  integer,pointer :: tab(:)          SET2NULL
  ! tab(nbz)
  ! For each point in the BZ, it gives the index of the symmetric irreducible point in the array ibz.

  integer,pointer :: tabi(:)         SET2NULL
  ! tabi(nbz)
  ! For each point in the BZ, tabi tells whether time-reversal has to be
  ! used to obtain k_BZ starting from the corresponding point in the IBZ  (1=>no, -1=>yes)

  integer,pointer :: tabo(:)         SET2NULL
  ! tabo(nbz)
  ! For each point in the BZ, it gives the index in the array symrec of the
  ! symmetry operation in reciprocal space which rotates k_IBZ onto \pm k_BZ (depending on tabi)

  integer,pointer :: umklp(:,:)   SET2NULL
   ! umklp(3,nbz)
   ! The Umklapp G0-vector such as kbz + G0 = (IS) k_ibz, where kbz is in the first BZ.

  real(dp) :: gmet(3,3)
  ! Reciprocal space metric ($\textrm{bohr}^{-2}$).

  real(dp) :: gprimd(3,3)
  ! Dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)

  real(dp),pointer :: bz(:,:)     SET2NULL
  ! bz(3,nbz)
  ! Points in the BZ in reduced coordinates.
  ! TODO should packed in shells.

  real(dp),pointer :: ibz(:,:)    SET2NULL
  ! ibz(3,nibz)
  ! Points in the IBZ in reduced coordinates.

  real(dp),pointer :: shift(:,:)  SET2NULL
  !  shift(3,nshift)
  !  Shift for k-points, not available is nshift=0. Usually nshift=1

  real(dp),pointer :: wt(:)   SET2NULL
  ! wt(nbz)
  ! Weights for each point in the IBZ.

  !%real(dp),pointer :: vbox(:)   SET2NULL
  ! vbox(nkbz)
  ! Volume of the small box centered on the k-point in the full BZ.
  ! Mainly used for inhomogeneous meshes.
  
  complex(dpc),pointer :: tabp(:)   SET2NULL   
  ! tabp(nkbz) 
  ! For each point in the BZ, this table gives the phase factors associated
  ! to non-symmorphic operations, i.e., e^{-i2\pi k_IBZ.R{^-1}t}=e^{-i2\pi k_BZ cdot t}
  ! where \transpose R{-1}=S and  (S k_IBZ)=\pm k_BZ (depending on tabi)

  type(t_tetrahedron), pointer :: tetrahedra

 end type BZ_mesh_type
!!***

!----------------------------------------------------------------------

 public :: init_kmesh            ! Main creation method.
 public :: nullify_BZ_mesh       ! Nullify all pointers in BZ_mesh_type.
 public :: destroy_bz_mesh_type  ! Free the structure.
 public :: print_bz_mesh         ! Printout of basic info on the object.
 public :: get_bz_item           ! Get point in the  BZ as well as useful quantities.
 public :: get_IBZ_item          ! Get point in the IBZ as well as useful quantities. 
 public :: get_BZ_diff           ! Get the difference k1-k2 in the BZ (if any).
 public :: isamek                ! Check whether two points are equal within an umklapp G0.
 public :: isequalk              ! Check whether two points are equal within an umklapp G0 (does not report G0)
 public :: has_BZ_item           ! Check if a point belongs to the BZ mesh.
 public :: has_IBZ_item          ! Check if a point is in the IBZ
 public :: bz_mesh_isirred       ! TRUE. if ikbz is in the IBZ (a non-zero umklapp is not allowed)
 public :: make_mesh             ! Initialize the mesh starting from kptrlatt and shift.
 public :: init_tetra_tab        ! Initialize tables for tetrahedron integration scheme. 
 public :: identk                ! Find the BZ starting from the irreducible k-points.
 public :: get_ng0sh             ! Calculate the smallest box in RSpace able to treat all possible umlapp processes.
 public :: make_path             ! Construct a normalized path.
 public :: find_qmesh            ! Find the Q-mesh defined as the set of all possible k1-k2 differences.
 public :: findnq                ! Helper routine returning the number of q-points. 
 public :: findq                 ! Helper routine returning the list of q-points. 
 public :: findqg0               ! Identify q + G0 = k1-k2.
 public :: box_len               ! Return the length of the vector connecting the origin with one the faces of the unit cell.

!----------------------------------------------------------------------

!!****t* m_bz_mesh/little_group
!! NAME
!! little_group
!!
!! FUNCTION
!! For the GW part of ABINIT. The little_group structured datatype gather information on
!! the little group associated to an external vector q. The little group associated to q
!! is defined as the subset of the space group that preserves q, modulo a G0 vector
!! (also called umklapp vector). Namely
!!
!!  Sq = q +G0,  where S is an operation in reciprocal space.
!!
!! If time reversal symmetry holds true, it is possible to enlarge the little group by
!! including the operations such as
!!  -Sq = q+ G0.
!!
!! The operations belongin to the little group define an irriducible wedge in the Brillouin zone
!! that is, usually, larger than the irredubile zone defined by the space group.
!! The two zone coincide when q=0
!!
!! TODO
!! Rationalize most of the arrays, in particular the tables
!! This structre shoud be rewritten almost from scratch, thus avoid using it
!! for your developments.
!!
!! SOURCE

 type,public :: little_group

  integer :: npw             ! No. of planewaves used to describe the wavefuntion, used to dimension igmG0
  integer :: nsym_sg         ! No. of operations in the space group (*NOT* the little group)
  integer :: nsym_ltg        ! No. of symmetry operations in the little group (time-reversal is included, if can be used)
  integer :: timrev          ! 2 if time-reversal is considered, 1 otherwise
  integer :: nbz             ! No. of kpoints in the full BZ
  integer :: nibz_ltg        ! No. of points in the irreducible wedge defined by the little group
  !integer :: use_umklp      ! 1 if umklapp processes are included

  real(dp) :: max_kin_gmG0
  ! Max kinetic energy of G-G0 in case of umklapp.

  integer,pointer :: G0(:,:,:)    SET2NULL
  ! g0(3,timrev,nsym_sg)
  ! Reduced coordinates of the umklapp G0 vector.

  integer,pointer :: ibzq(:)      SET2NULL
  ! ibzq(nbz)
  ! 1 if the point belongs to the IBZ_q defined by ext_pt, 0 otherwise.

  integer,pointer :: bz2ibz(:)    SET2NULL
  ! bz2ibz(nbz)
  ! Index of the point in the irreducible wedge defined by the little group, 0 otherwise.

  integer,pointer :: ibz2bz(:)    SET2NULL
  ! ibz2bz(nibz_ltg)
  ! The correspondind index in the BZ array

  integer,pointer :: igmG0(:,:,:)  SET2NULL
  ! iumklp(npw,timrev,nsym_sg)
  ! Index of G-G0 in the FFT array for each operations IS (I=\pm 1).

  integer,pointer :: flag_umklp(:,:)  SET2NULL
  ! flag_umklp(timrev,nsym_sg)
  ! 1 if the operation IS requires a non null G0 vector to preserve q, 0 otherwise.

  integer,pointer :: preserve(:,:)    SET2NULL
  ! preserve(timrev,nsym_sg)
  ! preserve(1,S) is 1 if the operation S in rec space preserves the external q-point i.e Sq=q+G0
  ! preserve(2,S) is 1 if -Sq=q+G0. G0 is a reciprocal lattice vector also called "umklapp vector".

  integer,pointer :: tab(:)   SET2NULL
  ! tab(nbz)
  ! For each point in BZ, the index of the irreducible point (kIBZ_q) in the irreducible
  ! wedge defined by the little group of q. kBZ= (IS) kIBZ where I is the inversion or the identity.

  integer,pointer :: tabo(:)  SET2NULL
  ! tabo(nbz)
  ! The index of the operation S in the little group that rotates kIBZ_q into \pm kBZ.

  integer,pointer :: tabi(:)  SET2NULL
  ! tabi(nbz)
  ! for each k-point in the BZ defines whether inversion has to be
  ! considered in the relation kBZ= IS kIBZ_q (1 => only S; -1 => -S).

  integer,pointer :: wtksym(:,:,:)  SET2NULL
  ! wtksym(timrev,nsym_sg,kbz)
  ! 1 if IS belongs to the little group, 0 otherwise !(should invert firt dimensions.

  real(dp) :: ext_pt(3)
  ! The external point defining the little group.

 end type little_group
!!***

 public :: setup_little_group  
 public :: nullify_little_group
 public :: destroy_little_group
 public :: print_little_group

 interface nullify_little_group
   module procedure nullify_little_group_0D
   module procedure nullify_little_group_1D
 end interface nullify_little_group

 interface destroy_little_group
   module procedure destroy_little_group_0D
   module procedure destroy_little_group_1D
 end interface destroy_little_group

CONTAINS  !=============================================================================
!!***

!!****f* m_bz_mesh/init_kmesh
!! NAME
!! init_kmesh
!!
!! FUNCTION
!!  Initialize and construct a bz_mesh_type datatype
!!  gathering information on the mesh in the Brilloin zone. 
!!
!! INPUTS
!!  nkibz=Number of irreducible k-points.
!!  kibz(3,nkibz)=Irreducible k-points in reduced coordinates.
!!  Cryst<Crystal_structure> = Info on unit cell and its symmetries
!!     %nsym=number of symmetry operations
!!     %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!     %tnons(3,nsym)=fractional translations
!!  kptopt=option for the generation of k points (see input variable description)
!!  [wrap_1zone]=If .TRUE., the points are wrapped in in the first BZ. Defaults to .FALSE. to preserve GW implementation.
!!  [ref_bz(:,:)]= Reference set of points in the full Brillouin zone
!!   used to prune k-points.
!!
!! OUTPUT
!!  Kmesh<bz_mesh_type>=Datatype gathering information on the k point sampling.
!!
!! PARENTS
!!      cchi0q0_intraband,m_bz_mesh,m_ebands,m_gwannier,m_screening,mlwfovlp_qp
!!      mrgscr,setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine init_kmesh(Kmesh,Cryst,nkibz,kibz,kptopt,wrap_1zone,ref_bz,break_symmetry)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_kmesh'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkibz,kptopt
 logical,optional,intent(in) :: wrap_1zone,break_symmetry
 type(BZ_mesh_type),intent(inout) :: Kmesh
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 real(dp),intent(in) :: kibz(3,nkibz)
 real(dp),optional,intent(in) :: ref_bz(:,:)

!Local variables-------------------------------
!scalars
 integer :: ik_bz,ik_ibz,isym,nkbz,nkbzX,nsym,timrev,itim
 real(dp) :: shift1,shift2,shift3
 logical :: ltest,do_wrap,do_hack
 character(len=500) :: msg
!arrays
 integer,allocatable :: ktab(:),ktabi(:),ktabo(:)
 integer,pointer :: symafm(:),symrec(:,:,:)
 real(dp) :: rm1t(3),kbz_wrap(3)
 real(dp),allocatable :: kbz(:,:),wtk(:)
 real(dp),pointer :: tnons(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 !@BZ_mesh_type
 call nullify_BZ_mesh(Kmesh)

 ! === Initial tests on input arguments ===
 ltest=(Cryst%timrev==1.or.Cryst%timrev==2)
 write(msg,'(a,i0)')'Wrong value for timrev= ',Cryst%timrev
 ABI_CHECK(ltest,msg)

 if (ALL(kptopt/=(/1,3/))) then
   write(msg,'(a,i0)')" Not allowed value for kptopt: ",kptopt
   MSG_BUG(msg)
 end if

 Kmesh%kptopt = kptopt

 nsym   = Cryst%nsym
 timrev = Cryst%timrev

 symrec => Cryst%symrec
 symafm => Cryst%symafm
 tnons  => Cryst%tnons
 !
 ! === Find BZ from IBZ and fill tables ===
 nkbzX=nkibz*nsym*timrev ! Maximum possible number
 ABI_ALLOCATE(kbz,(3,nkbzX))
 ABI_ALLOCATE(wtk,(nkibz))
 ABI_ALLOCATE(ktab,(nkbzX))
 ABI_ALLOCATE(ktabi,(nkbzX))
 ABI_ALLOCATE(ktabo,(nkbzX))

 if (PRESENT(ref_bz)) then
   call identk(kibz,nkibz,nkbzX,nsym,timrev,symrec,symafm,kbz,ktab,ktabi,ktabo,nkbz,wtk,ref_bz=ref_bz)
 else
   call identk(kibz,nkibz,nkbzX,nsym,timrev,symrec,symafm,kbz,ktab,ktabi,ktabo,nkbz,wtk)
 end if

 ! TODO: Force the k-points to be in the first Brillouin zone.
 !  Now the GW tests seem to be OK, additional tests have to be done though.
 do_wrap=.FALSE.; if (PRESENT(wrap_1zone)) do_wrap=wrap_1zone
 !do_wrap=.TRUE.

 if (do_wrap) then ! Wrap the BZ points in the interval ]-1/2,1/2]
   do ik_bz=1,nkbz
     call wrap2_pmhalf(kbz(1,ik_bz),kbz_wrap(1),shift1)
     call wrap2_pmhalf(kbz(2,ik_bz),kbz_wrap(2),shift2)
     call wrap2_pmhalf(kbz(3,ik_bz),kbz_wrap(3),shift3)
     kbz(:,ik_bz) = kbz_wrap
   end do
 end if
 !
 ! ================================================================
 ! ==== Create data structure to store information on k-points ====
 ! ================================================================
 !
 ! * Dimensions.
 Kmesh%nbz   =nkbz      ! Number of points in the full BZ
 Kmesh%nibz  =nkibz     ! Number of points in the IBZ
 Kmesh%nsym  =nsym      ! Number of operations
 Kmesh%timrev=timrev    ! 2 if time-reversal is used, 1 otherwise

 !
 ! * Arrays.
 Kmesh%gmet   = Cryst%gmet
 Kmesh%gprimd = Cryst%gprimd

 ABI_ALLOCATE(Kmesh%bz ,(3,nkbz))
 Kmesh%bz   =  kbz(:,1:nkbz )  ! Red. coordinates of points in full BZ.
 ABI_ALLOCATE(Kmesh%ibz,(3,nkibz))
 Kmesh%ibz  = kibz(:,1:nkibz)  ! Red. coordinates of points in IBZ.

 ABI_ALLOCATE(Kmesh%tab ,(nkbz))
 Kmesh%tab  = ktab (1:nkbz)    ! Index of the irred. point in the array IBZ.
 ABI_ALLOCATE(Kmesh%tabi,(nkbz))
 Kmesh%tabi = ktabi(1:nkbz)    !-1 if time reversal must be used to obtain this point,
 ABI_ALLOCATE(Kmesh%tabo,(nkbz))
 Kmesh%tabo = ktabo(1:nkbz)    ! Symm. operation that rotates k_IBZ onto \pm k_BZ
                                                             ! (depending on tabi)
 ABI_ALLOCATE(Kmesh%wt,(nkibz))
 Kmesh%wt(:)= wtk(1:nkibz)     ! Weight for each k_IBZ

 ABI_ALLOCATE(Kmesh%rottbm1,(nkbz,timrev,nsym))
 ABI_ALLOCATE(Kmesh%rottb  ,(nkbz,timrev,nsym))
                                                              ! identity or time-reversal

 do_hack = .FALSE.
 if (PRESENT(ref_bz) .and. PRESENT(break_symmetry)) then
   do_hack = break_symmetry
 end if

 if (do_hack) then
   MSG_WARNING("Hacking the rottb tables!")
   do ik_bz=1,nkbz
     Kmesh%rottbm1(ik_bz,:,:) = ik_bz
     Kmesh%rottb  (ik_bz,:,:) = ik_bz
   end do
 else 
   call setup_k_rotation(nsym,timrev,symrec,nkbz,Kmesh%bz,Cryst%gmet,Kmesh%rottb,Kmesh%rottbm1)
 end if

 ! TODO umklp can be calculated inside setup_k_rotation.
 ABI_ALLOCATE(Kmesh%umklp,(3,nkbz))
 do ik_bz=1,nkbz
   ik_ibz= Kmesh%tab (ik_bz)
   isym  = Kmesh%tabo(ik_bz) 
   itim  = (3-Kmesh%tabi(ik_bz))/2
   Kmesh%umklp(:,ik_bz) = NINT( -Kmesh%bz(:,ik_bz) + (3-2*itim)*MATMUL(symrec(:,:,isym),Kmesh%ibz(:,ik_ibz)) )
 end do

 ABI_ALLOCATE(Kmesh%tabp,(nkbz))
 do ik_bz=1,nkbz
   isym  =Kmesh%tabo(ik_bz) 
   ik_ibz=Kmesh%tab (ik_bz)
   rm1t=MATMUL(TRANSPOSE(symrec(:,:,isym)),tnons(:,isym))
   Kmesh%tabp(ik_bz)=EXP(-(0.,1.)*two_pi*DOT_PRODUCT(kibz(:,ik_ibz),rm1t))
 end do

 ABI_DEALLOCATE(kbz)
 ABI_DEALLOCATE(wtk)
 ABI_DEALLOCATE(ktab)
 ABI_DEALLOCATE(ktabi)
 ABI_DEALLOCATE(ktabo)

 ! allocate structure pointer for tetrahedra (subarrays are not allocated yet!)
 ABI_ALLOCATE(Kmesh%tetrahedra,)
 call nullify_tetra(Kmesh%tetrahedra)

 DBG_EXIT("COLL")

end subroutine init_kmesh
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/my_reset_bz_mesh
!! NAME
!! my_reset_bz_mesh
!!
!! FUNCTION
!! Reset dimensions in a BZ_mesh_type. [PRIVATE]
!!
!! INPUTS
!! Kmesh<BZ_mesh_type>= The datatype whose dimensions have to be set to zero
!!
!! SIDE EFFECTS
!! All dimension reset to 0.
!!
!! PARENTS
!!      m_bz_mesh
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine my_reset_bz_mesh(Kmesh)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'my_reset_bz_mesh'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(BZ_mesh_type),intent(inout) :: Kmesh

! *********************************************************************

 !@BZ_mesh_type

!dimensions
 Kmesh%nshift     = 0 
 Kmesh%nbz        = 0 
 Kmesh%nibz       = 0 

end subroutine my_reset_bz_mesh
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/nullify_BZ_mesh
!! NAME
!! nullify_BZ_mesh
!!
!! FUNCTION
!! Nullify the pointers in a BZ_mesh_type.
!!
!! INPUTS
!! Kmesh<BZ_mesh_type>= The datatype whose pointers have to be nullified.
!!
!! SIDE EFFECTS
!! All pointers set to null().
!!
!! PARENTS
!!      m_bz_mesh,m_gwannier,m_shexc,mlwfovlp_qp
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine nullify_BZ_mesh(Kmesh)

 use defs_basis
 use m_tetrahedron

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_BZ_mesh'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(BZ_mesh_type),intent(inout) :: Kmesh

! *********************************************************************

 !@BZ_mesh_type

! zero dimensions
 call my_reset_bz_mesh(Kmesh)

!integer
 nullify(Kmesh%rottb     )
 nullify(Kmesh%rottbm1   )
 nullify(Kmesh%tab       )
 nullify(Kmesh%tabi      )
 nullify(Kmesh%tabo      )
 nullify(Kmesh%umklp     )

!real
 nullify(Kmesh%bz   )
 nullify(Kmesh%ibz  )
 nullify(Kmesh%shift)
 nullify(Kmesh%wt   )

!complex
 nullify(Kmesh%tabp)

 nullify(Kmesh%tetrahedra)

end subroutine nullify_BZ_mesh
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/destroy_bz_mesh_type
!! NAME
!! destroy_bz_mesh_type
!!
!! FUNCTION
!! Deallocate all dynamics entities present in a BZ_mesh_type structure.
!!
!! INPUTS
!! Kmesh<BZ_mesh_type>=The datatype to be freed.
!!
!! SIDE EFFECTS
!! All allocated memory is released. 
!!
!! PARENTS
!!      bethe_salpeter,cchi0q0_intraband,exc_interp_ham,m_ebands,m_gwannier
!!      m_screening,m_shexc,mlwfovlp_qp,mrgscr,screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_bz_mesh_type(Kmesh)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_bz_mesh_type'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(BZ_mesh_type),intent(inout) :: Kmesh

! *********************************************************************

 DBG_ENTER("COLL")

 !@BZ_mesh_type

!integer
 if (associated(Kmesh%rottb     ))  then
   ABI_DEALLOCATE(Kmesh%rottb)
 end if
 if (associated(Kmesh%rottbm1   ))  then
   ABI_DEALLOCATE(Kmesh%rottbm1)
 end if
 if (associated(Kmesh%tab       ))  then
   ABI_DEALLOCATE(Kmesh%tab)
 end if
 if (associated(Kmesh%tabi      ))  then
   ABI_DEALLOCATE(Kmesh%tabi)
 end if
 if (associated(Kmesh%tabo      ))  then
   ABI_DEALLOCATE(Kmesh%tabo)
 end if
 if (associated(Kmesh%umklp     ))  then
   ABI_DEALLOCATE(Kmesh%umklp)
 end if

!real
 if (associated(Kmesh%ibz  ))  then
   ABI_DEALLOCATE(Kmesh%ibz)
 end if
 if (associated(Kmesh%bz   ))  then
   ABI_DEALLOCATE(Kmesh%bz)
 end if
 if (associated(Kmesh%shift))  then
   ABI_DEALLOCATE(Kmesh%shift)
 end if
 if (associated(Kmesh%wt   ))  then
   ABI_DEALLOCATE(Kmesh%wt)
 end if

!complex
 if (associated(Kmesh%tabp))  then
   ABI_DEALLOCATE(Kmesh%tabp)
 end if

 call destroy_tetra(Kmesh%tetrahedra)
 if (associated(Kmesh%tetrahedra))  then
   ABI_DEALLOCATE(Kmesh%tetrahedra)
 end if

!set dimensions to zero
 call my_reset_bz_mesh(Kmesh)

 DBG_EXIT("COLL")

end subroutine destroy_bz_mesh_type
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/print_bz_mesh
!! NAME
!! print_bz_mesh
!!
!! FUNCTION
!! Print the content of a bz_mesh_type datatype
!!
!! INPUTS
!! Kmesh<bz_mesh_type>=the datatype to be printed
!! [unit]=the unit number for output 
!! [prtvol]=verbosity level
!! [mode_paral]=either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only printing.
!!
!! PARENTS
!!      m_screening,m_shexc,mrgscr,setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_BZ_mesh(Kmesh,header,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_BZ_mesh'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(BZ_mesh_type),intent(in) :: Kmesh

!Local variables-------------------------------
!scalars
 integer,parameter :: nmaxk=50
 integer :: ii,ik,my_unt,my_prtvol
 character(len=100) :: fmt
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0    ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the Kmesh% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a,i5,3a)')&
&  ' Number of points in the irreducible wedge : ',Kmesh%nibz,ch10,&
&  ' Reduced coordinates and weights : ',ch10
 call wrtout(my_unt,msg,my_mode)

 write(fmt,*)'(1x,i5,a,2x,3es16.8,3x,f11.5)'
 do ik=1,Kmesh%nibz ! Add tol8 for portability reasons.
   write(msg,fmt) ik,') ',(Kmesh%ibz(ii,ik),ii=1,3),Kmesh%wt(ik)+tol8
   call wrtout(my_unt,msg,my_mode)
 end do

 SELECT CASE (Kmesh%timrev)

 CASE (1)
   write(msg,'(2a,i2,3a,i5,a)')ch10,&
&    ' Together with ',Kmesh%nsym,' symmetry operations (time-reversal symmetry not used) ',ch10,&
&    ' yields ',Kmesh%nbz,' points in the full Brillouin Zone.'

 CASE (2) 
   write(msg,'(2a,i2,3a,i5,a)')ch10,&
&    ' Together with ',Kmesh%nsym,' symmetry operations and time-reversal symmetry ',ch10,&
&    ' yields ',Kmesh%nbz,' points in the full Brillouin Zone.'

 CASE DEFAULT
   write(msg,'(a,i3)')' Wrong value for timrev = ',Kmesh%timrev
   MSG_BUG(msg)
 END SELECT

 call wrtout(my_unt,msg,my_mode)

 if (my_prtvol>0) then
   write(fmt,*)'(1x,i5,a,2x,3es16.8)'
   do ik=1,Kmesh%nbz
     if (my_prtvol==1 .and. ik>nmaxk) then
       write(msg,'(a)')' prtvol=1, do not print more points.'
       call wrtout(my_unt,msg,my_mode) ; EXIT
     end if
     write(msg,fmt)ik,') ',(Kmesh%bz(ii,ik),ii=1,3)
     call wrtout(my_unt,msg,my_mode)
   end do
 end if
 !
 ! === Additional printing ===
 if (my_prtvol>=10) then
   write(msg,'(2a)')ch10,&
&   '                  Full point  ------->    Irred point -->            through:  Symrec  Time-Rev (1=No,-1=Yes) G0(1:3) '
   call wrtout(my_unt,msg,my_mode)
   write(fmt,*)'(2x,i5,2x,2(3(f7.4,2x)),i3,2x,i2,3(i3))'
   do ik=1,Kmesh%nbz
     write(msg,fmt)ik,Kmesh%bz(:,ik),Kmesh%ibz(:,Kmesh%tab(ik)),Kmesh%tabo(ik),Kmesh%tabi(ik),Kmesh%umklp(:,ik)
     call wrtout(my_unt,msg,my_mode)
   end do
 end if

 write(msg,'(a)')ch10
 call wrtout(my_unt,msg,my_mode)

end subroutine print_BZ_mesh
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/setup_k_rotation
!! NAME
!! setup_k_rotation
!!
!! FUNCTION
!! Set up tables giving the correspondence btw a k-point and its rotated image. 
!!
!! INPUTS
!! timrev=2 if time-reversal can be used, 1 otherwise.
!! nsym=Number of symmetry operations
!! symrec(3,3,nsym)=Symmetry operations in reciprocal space in reduced coordinates.
!! nbz=Number of k-points
!! kbz(3,nbz)=k-points in reduced coordinates.
!! gmet(3,3)=Metric in reciprocal space.
!!
!! OUTPUT
!! krottb(k,I,S)=Index of (IS) k in the array bz
!! krottbm1(k,I,S)=Index of IS^{-1} k
!!
!! PARENTS
!!      m_bz_mesh
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine setup_k_rotation(nsym,timrev,symrec,nbz,kbz,gmet,krottb,krottbm1)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_k_rotation'
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbz,nsym,timrev
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(out) :: krottb(nbz,timrev,nsym),krottbm1(nbz,timrev,nsym)
 real(dp),intent(in) :: kbz(3,nbz),gmet(3,3)

!Local variables ------------------------------
!scalars
 integer :: ik,ikp,isym,itim,nsh,ik_st !,sh_start
 real(dp),parameter :: KTOL=tol6
 real(dp) :: norm_old,norm !,norm_rot !norm_base
 logical :: found,isok
 character(len=500) :: msg
!arrays
 integer :: g0(3)
 integer :: iperm(nbz),shlim(nbz+1)
 real(dp) :: shift(3),kbase(3),krot(3),knorm(nbz),kwrap(3),shlen(nbz+1)

!************************************************************************

 DBG_ENTER("COLL")

 ! Sort the k-points according to their norm to speed up the search below.
 do ik=1,nbz
   call wrap2_pmhalf(kbz(1,ik),kwrap(1),shift(1))
   call wrap2_pmhalf(kbz(2,ik),kwrap(2),shift(2))
   call wrap2_pmhalf(kbz(3,ik),kwrap(3),shift(3))
   knorm(ik) = normv(kwrap,gmet,"G")
   iperm(ik)= ik
 end do

 call sort_dp(nbz,knorm,iperm,KTOL)
 !
 ! The index of the initial (sorted) k-point in each shell
 nsh=1; norm_old=knorm(1)

 shlim(1)=1; shlen(1) = norm_old
 do ik_st=2,nbz
   norm = knorm(ik_st) 
   if (ABS(norm-norm_old) > KTOL) then 
     norm_old   = norm
     nsh        = nsh+1
     shlim(nsh) = ik_st
     shlen(nsh) = norm 
   end if
 end do
 shlim(nsh+1)=nbz+1
 shlen(nsh+1)=HUGE(one)
 !
 ! === Set up k-rotation tables ===
 ! * Use spatial inversion instead of time reversal whenever possible.
 !call wrtout(std_out," Begin sorting ","COLL")       
 !call flush_unit(std_out)

 isok=.TRUE.
 do ik=1,nbz
   kbase(:)=kbz(:,ik)
   !
   do itim=1,timrev
     do isym=1,nsym
       krot(:)=(3-2*itim)*MATMUL(symrec(:,:,isym),kbase)

       found=.FALSE.
#if 1
      ! Old code
       do ikp=1,nbz
         if (isamek(krot,kbz(:,ikp),g0)) then
           found=.TRUE.
           krottb  (ik ,itim,isym)=ikp
           krottbm1(ikp,itim,isym)=ik
           EXIT
         end if
       end do
#else
       !
       ! Locate the shell index with bisection.
       call wrap2_pmhalf(krot(1),kwrap(1),shift(1))
       call wrap2_pmhalf(krot(2),kwrap(2),shift(2))
       call wrap2_pmhalf(krot(3),kwrap(3),shift(3))
       norm_rot = normv(kwrap,gmet,"G")
       sh_start = bisect(shlen(1:nsh+1),norm_rot) 

       do ik_st=shlim(sh_start),nbz
         ikp = iperm(ik_st)
         if (isamek(krot,kbz(:,ikp),g0)) then
           found=.TRUE.
           krottb  (ik ,itim,isym)=ikp
           krottbm1(ikp,itim,isym)=ik
           !write(std_out,*)ik_st,shlim(sh_start),nbz
           EXIT
         end if
       end do
#endif
       if (.not.found) then
         isok=.FALSE.
         !write(std_out,*)" norm_base,norm_rot ",norm_base,norm_rot
         !write(std_out,*)normv(kbase,gmet,"G"),normv(krot,gmet,"G")
         write(msg,'(2(a,i4),2x,2(3f12.6,2a),i3,a,i2)')&
&          ' Initial k-point ',ik,'/',nbz,kbase(:),ch10,&
&          ' Rotated k-point (not found) ',krot(:),ch10,&
&          ' Through symmetry operation ',isym,' and itim ',itim
         !MSG_WARNING(msg)
         MSG_ERROR(msg)
       end if

     end do
   end do
 end do

 if (.not.isok) then
   MSG_ERROR('k-mesh not closed')
 end if

 DBG_EXIT("COLL")

end subroutine setup_k_rotation
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/get_bz_item
!! NAME
!! get_bz_item
!!
!! FUNCTION
!! Given the index of a point in the full BZ, this routine returns the index of the 
!! symmetric image in the IBZ, the index of the symmetry operation symrec needed, 
!! whether time-reversal has to be used.
!! Optionally the non-symmorphic phase and the umklapp vector is returned.
!!
!! INPUTS
!! ikbz=The index of the required point in the BZ
!! Kmesh<BZ_mesh_type>=Datatype gathering information on the k point sampling. 
!!
!! OUTPUT
!! kbz(3)=The k-point in the first BZ in reduced coordinated.
!! isym=Index of the symrec symmetry required to rotate ik_ibz onto ik_bz.
!! itim=2 is time-reversal has to be used, 1 otherwise
!! ik_ibz=The index of the corresponding symmetric point in the IBZ.
!! [ph_mkbzt]=The phase factor for non-symmorphic operations  e^{-i 2 \pi k_IBZ \cdot R{^-1}t}=e{-i 2\pi k_BZ cdot t}
!! [umklp(3)]=The umklapp G0 vector such as kbz + G0 = (IS) k_ibz, where kbz is in the BZ.
!! [isirred]=.TRUE. if the k-point belongs to IBZ.
!!
!! PARENTS
!!      calc_optical_mels,calc_sigc_me,calc_sigx_me,cchi0,cchi0q0
!!      cchi0q0_intraband,check_completeness,cohsex_me,debug_tools
!!      exc_build_block,exc_build_ham,m_dyson_solver,m_ebands,m_oscillators
!!      m_ppmodel,m_screen,m_screening,m_shexc,m_shirley,m_vcoul,m_wfs
!!      paw_symcprj,setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine get_bz_item(Kmesh,ik_bz,kbz,ik_ibz,isym,itim,ph_mkbzt,umklp,isirred)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_bz_item'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz
 integer,intent(out) :: ik_ibz,isym,itim
 complex(dpc),optional,intent(out) :: ph_mkbzt
 logical,optional,intent(out) :: isirred
 type(BZ_mesh_type),intent(in) :: Kmesh
!arrays
 integer,optional,intent(out) :: umklp(3)
 real(dp),intent(out) :: kbz(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *********************************************************************

 if (ik_bz>Kmesh%nbz.or.ik_bz<=0) then
   write(msg,'(a,i3)')' Wrong value for ik_bz: ',ik_bz
   MSG_BUG(msg)
 end if

 kbz    = Kmesh%bz(:,ik_bz)
 ik_ibz = Kmesh%tab(ik_bz)
 isym   = Kmesh%tabo(ik_bz)
 itim   = (3-Kmesh%tabi(ik_bz))/2

 if (PRESENT(ph_mkbzt)) ph_mkbzt=Kmesh%tabp(ik_bz)
 if (PRESENT(umklp))    umklp   =Kmesh%umklp(:,ik_bz)
 ! Be careful here as we assume a particular ordering of symmetries.
 if (PRESENT(isirred))  isirred = ( isym==1.and.itim==1.and.ALL(Kmesh%umklp(:,ik_bz)==(/0,0,0/)) ) 

end subroutine get_bz_item
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/get_IBZ_item
!! NAME
!! get_IBZ_item
!!
!! FUNCTION
!! Report useful information on a k-point in the IBZ starting from its sequential index in %ibz. 
!!
!! INPUTS
!! ik_ibz=The index of the required point in the IBZ
!! Kmesh<bz_mesh_type>=datatype gathering information on the k point sampling. 
!!
!! OUTPUT
!! kibz(3)=the k-point in reduced coordinated
!! wtk=the weight
!!
!! TODO
!!  Add mapping ibz2bz, ibz2star
!!
!! PARENTS
!!      paw_symcprj
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine get_IBZ_item(Kmesh,ik_ibz,kibz,wtk)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_IBZ_item'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz
 real(dp),intent(out) :: wtk
 type(bz_mesh_type),intent(in) :: Kmesh
!arrays
 real(dp),intent(out) :: kibz(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *********************************************************************

 if (ik_ibz>Kmesh%nibz.or.ik_ibz<=0) then
   write(msg,'(a,i3)')' wrong value for ik_ibz: ',ik_ibz
   MSG_BUG(msg)
 end if

 kibz=Kmesh%ibz(:,ik_ibz)
 wtk =Kmesh%wt(ik_ibz)

end subroutine get_IBZ_item
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/get_BZ_diff
!! NAME
!! get_BZ_diff
!!
!! FUNCTION
!! Given two points k1 and k2 where k1 belongs to the BZ, check if the difference
!! k1-k2 still belongs to the BZ reporting useful quantities
!!
!! INPUTS
!!  Kmesh<bz_mesh_type>=datatype gathering information on the k-mesh
!!  k1(3)=the first k-points (supposed to be in the BZ)
!!  k2(3)=the second point
!!
!! OUTPUT
!!  idiff_bz=the idex of k1-k2 in the BZ
!!  G0(3)=the umklapp G0 vector required to bring k1-k2 back to the BZ
!!  nfound= the number of points in the BZ that are equal to k1-k2 (should be 1 if everything is OK)
!!
!! PARENTS
!!      cchi0,check_completeness,m_ebands,m_oscillators
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine get_BZ_diff(Kmesh,k1,k2,idiff_bz,g0,nfound)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_BZ_diff'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: idiff_bz,nfound
 type(bz_mesh_type),intent(in) :: Kmesh
!arrays
 integer,intent(out) :: g0(3)
 real(dp),intent(in) :: k1(3),k2(3)

!Local variables-------------------------------
!scalars
 integer :: ikp
 character(len=500) :: msg
!arrays
 integer :: umklp(3)
 real(dp) :: kdiff(3),ktrial(3)

! *********************************************************************

 if (.not.has_BZ_item(Kmesh,k1,ikp,umklp)) then
   write(msg,'(a,3f12.6)')' first point must be in BZ: ',k1
   MSG_ERROR(msg)
 end if

 kdiff   = k1-k2
 nfound  = 0 
 idiff_bz= 0
 !
 ! === Find p such k1-k2=p+g0 where p in the BZ ===
 do ikp=1,Kmesh%nbz
   ktrial=Kmesh%bz(:,ikp)
   if (isamek(kdiff,ktrial,umklp)) then
     idiff_bz=ikp
     g0=umklp
     nfound=nfound+1
   end if
 end do
 !
 ! === Check if p has not found of found more than once ===
 ! * For extremely dense meshes, tol1q in defs_basis might be too large!
 if (nfound/=1) then
   if (nfound==0) then
     MSG_WARNING(" k1-k2-G0 not found in BZ")
   else
     write(msg,'(a,i3)')' Multiple k1-k2-G0 found in BZ, nfound= ',nfound
     MSG_WARNING(msg)
   end if
   write(msg,'(4a,3(a,3f12.6,a))')&
&    ' k1    = ',k1   ,ch10,&
&    ' k2    = ',k2   ,ch10,&
&    ' k1-k2 = ',kdiff,ch10
   MSG_WARNING(msg)
 end if

end subroutine get_BZ_diff
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/isamek
!! NAME
!! isamek
!!
!! FUNCTION
!! Test two k-points for equality. 
!! Return .TRUE. is they are equal within a reciprocal lattice vector G0.
!!
!! INPUTS
!!  k1(3),k2(3)=The two k points to be compared.
!!
!! OUTPUT
!! Return .TRUE. if they are the same within a RL vector,
!!        .FALSE. if they are different.
!! G0(3)=if .TRUE. G0(3) is the reciprocal lattice vector such as k1=k2+G0
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

logical function isamek(k1,k2,g0)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isamek'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(out) :: g0(3)
 real(dp),intent(in) :: k1(3),k2(3)

! *************************************************************************

 isamek = isinteger(k1-k2,TOL_KDIFF)

 if (isamek) then
   g0=NINT(k1-k2)
 else 
   g0=HUGE(1)
 end if

end function isamek
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/isequalk
!! NAME
!! is_equalk
!!
!! FUNCTION
!! Return .TRUE. if two points are equal within a reciprocal lattice vector.
!!
!! INPUTS
!!  q1(3),q2(3)=The two points to be compared for equivalence.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function isequalk(q1,q2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isequalk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 logical :: isequalk
 real(dp),intent(in) :: q1(3),q2(3)

!Local variables-------------------------------
!arrays
 integer :: g0(3)
! *************************************************************************

 isequalk = isamek(q1,q2,g0)

end function isequalk
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/has_BZ_item
!! NAME
!! has_BZ_item
!!
!! FUNCTION
!!  check if item belongs to the BZ  within a reciprocal lattice vector
!!  and return the index number and the reciprocal vector g0.
!!
!! INPUTS
!!  Kmesh<bz_mesh_type>=datatype gathering information on the k-mesh
!!  item(3)=the k-point to be checked
!!
!! OUTPUT
!!  .TRUE. if item is the BZ within a RL vector
!!  ikbz=Index of the k-point in the Kmesh%bz array
!!  g0(3)=Umklapp vector.
!!
!! PARENTS
!!
!! FIXME
!!  Switch to routine version. Due to side-effects the present implementation
!!  might be source of bugs in logical statements
!!
!! CHILDREN
!!
!! SOURCE

logical function has_BZ_item(Kmesh,item,ikbz,g0)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'has_BZ_item'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ikbz
 type(BZ_mesh_type),intent(in) :: Kmesh
!arrays
 integer,intent(out) :: g0(3)
 real(dp),intent(in) :: item(3)

!Local variables-------------------------------
!scalars
 integer :: ik_bz,yetfound
!arrays
 integer :: g0_tmp(3)

! *************************************************************************

 has_BZ_item=.FALSE.; ikbz=0; g0=0; yetfound=0
 do ik_bz=1,Kmesh%nbz
   if (isamek(item,Kmesh%bz(:,ik_bz),g0_tmp)) then
     has_BZ_item=.TRUE.
     ikbz=ik_bz
     g0 = g0_tmp
     yetfound=yetfound+1
     !EXIT
   end if
 end do

 if (yetfound/=0.and.yetfound/=1) then
   MSG_ERROR('Multiple k-points found')
 end if

end function has_BZ_item
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/has_IBZ_item
!! NAME
!! has_IBZ_item
!!
!! FUNCTION
!!  Check if item belongs to the IBZ within a reciprocal lattice vector
!!
!! INPUTS
!!  Kmesh<bz_mesh_type>=Datatype gathering information on the mesh in the BZ.
!!  item(3)=the k-point to be checked
!!
!! OUTPUT
!!  Return .TRUE. if item is the IBZ within a RL vector
!!  ikibz=The index of the k-point in the IBZ.
!!  g0(3)=The reciprocal lattice vector.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

logical function has_IBZ_item(Kmesh,item,ikibz,g0)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'has_IBZ_item'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ikibz
 type(BZ_mesh_type),intent(in) :: Kmesh
!arrays
 integer,intent(out) :: g0(3)
 real(dp),intent(in) :: item(3)

!Local variables-------------------------------
!scalars
 integer :: ik_ibz,yetfound
 !character(len=500) :: msg
!arrays
 integer :: g0_tmp(3)

! *************************************************************************

 has_IBZ_item=.FALSE.; ikibz=0; g0=0; yetfound=0
 do ik_ibz=1,Kmesh%nibz
   if (isamek(item,Kmesh%ibz(:,ik_ibz),g0_tmp)) then
     has_IBZ_item=.TRUE.
     ikibz=ik_ibz
     g0 = g0_tmp
     yetfound=yetfound+1
     !EXIT
   end if
 end do

 if (yetfound/=0.and.yetfound/=1) then
   MSG_BUG("multiple k-points found")
 end if

end function has_IBZ_item
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/bz_mesh_isirred
!! NAME
!! bz_mesh_isirred
!!
!! FUNCTION
!!  bz_mesh_isirred
!!
!! INPUTS
!!  ik_bz=Index of the k-point in the BZ.
!!
!! OUTPUT
!! bz_mesh_isirred=.TRUE. if the k-point is in the IBZ (a non-zero umklapp is not allowed)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function bz_mesh_isirred(Kmesh,ik_bz)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bz_mesh_isirred'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz
 logical :: bz_mesh_isirred
 type(BZ_mesh_type),intent(in) :: Kmesh
!arrays

!Local variables-------------------------------
!scalars
 integer :: isym,itim
 character(len=500) :: msg

! *********************************************************************

 if (ik_bz>Kmesh%nbz.or.ik_bz<=0) then
   write(msg,'(a,i3)')' Wrong value for ik_bz: ',ik_bz
   MSG_BUG(msg)
 end if

 isym   = Kmesh%tabo(ik_bz)
 itim   = (3-Kmesh%tabi(ik_bz))/2

 ! Be careful here as we assume a particular ordering of symmetries.
 bz_mesh_isirred = ( isym==1.and.itim==1.and.ALL(Kmesh%umklp(:,ik_bz)==(/0,0,0/)) ) 

end function bz_mesh_isirred
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/make_mesh
!! NAME
!! make_mesh
!!
!! FUNCTION
!! Initialize the BZ_mesh_type starting from kptrlatt and shiftk
!!
!! INPUTS
!! Cryst<Crystal_structure>=Info on the crystalline structure.
!! nshiftk=Number of shifts for the mesh.
!! kptrlatt(3,3)= Coordinates of three vectors in real space, expressed in reduced coordinates.
!!  They define a super-lattice in real space. The k point lattice is the reciprocal of
!!  this super-lattice, eventually shifted by shift.
!! shiftk(3,nshiftk)=Shifts for the k-mesh.
!! [want_tetra]= 1 if tables used for integration with tetrahedron methods are wanted.
!! [vacuum(3)]=For each direction, 0 if no vacuum, 1 if vacuum
!!
!! OUTPUT
!! Kmesh<BZ_mesh_type>=Object gathering info on the sampling of the Brillouin zone.
!!
!! PARENTS
!!      exc_interp_ham,m_gwannier,m_shexc,setup_bse
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine make_mesh(Kmesh,Cryst,kptopt,kptrlatt,nshiftk,shiftk,want_tetra,vacuum,break_symmetry)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_mesh'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nshiftk,kptopt
 integer,optional,intent(in) :: want_tetra
 logical,optional,intent(in) :: break_symmetry
 type(BZ_mesh_type),intent(inout) :: Kmesh
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(inout) :: kptrlatt(3,3) 
 integer,optional,intent(in) :: vacuum(3)
 real(dp),intent(in) :: shiftk(3,nshiftk)

!Local variables -------------------------
!scalars
 integer :: timrev,iscf
 integer :: nkbz,nkibz,nkpt_computed,my_nshiftk
 real(dp) :: kptrlen
 character(len=500) :: msg
 logical :: my_break_symmetry
!arrays
 integer :: my_vacuum(3)
 real(dp),allocatable :: kibz(:,:),wtk(:)
 real(dp),allocatable :: my_shiftk(:,:)
 real(dp),pointer :: ref_kbz(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 !@BZ_mesh_type
 call nullify_BZ_mesh(Kmesh)

 !if (ALL(kptopt/=(/1,2,3,4/))) then
 if (ALL(kptopt/=(/1,3/))) then
   write(msg,'(a,i0)')" Not allowed value for kptopt: ",kptopt
   MSG_BUG(msg)
 end if

 timrev=0; if (Cryst%timrev==2) timrev=1  !FIXME there an incompatibly between Cryst%timrev and symkpt
 !
 ! ======================================================================
 ! ==== First call to getkgrid to obtain nkibz as well as the BZ set ====
 ! ======================================================================
 iscf=7  ! use for the Weights in NSCF calculation. check it more carefully.
 nkibz=0 ! Compute number of k-points in the BZ and IBZ

 my_vacuum = (/0,0,0/); if (PRESENT(vacuum)) my_vacuum=vacuum
 
 my_nshiftk = nshiftk
 ABI_CHECK(my_nshiftk>0.and.my_nshiftk<=8,"Wrong nshiftk")
 ABI_ALLOCATE(my_shiftk,(3,8))
 my_shiftk=zero; my_shiftk(:,1:nshiftk) = shiftk(:,:)

 nullify(ref_kbz)

 write(std_out,*)" In make_mesh"
 write(std_out,*)" kptopt   ",kptopt
 write(std_out,*)" kptrlatt ",kptrlatt
 write(std_out,*)" nshiftk  ",nshiftk
 write(std_out,*)" shiftk   ",shiftk

 call getkgrid(0,0,iscf,kibz,kptopt,kptrlatt,kptrlen,Cryst%nsym,0,nkibz,my_nshiftk,&
&  Cryst%nsym,Cryst%rprimd,my_shiftk,Cryst%symafm,Cryst%symrel,my_vacuum,wtk,kbz_p=ref_kbz)

 nkbz = SIZE(ref_kbz,DIM=2)

 write(std_out,*)" after getkgrid1: nkbz = ",nkbz," nkibz=",nkibz
 write(std_out,*)" ref_kbz = ",ref_kbz

 ! ==============================================================
 ! ==== Recall getkgrid to get kibz(3,nkibz) and wtk(nkibz) =====
 ! ==============================================================
 ABI_ALLOCATE(kibz,(3,nkibz))
 ABI_ALLOCATE(wtk,(nkibz))

 call getkgrid(0,0,iscf,kibz,kptopt,kptrlatt,kptrlen,Cryst%nsym,nkibz,nkpt_computed,my_nshiftk,&
&  Cryst%nsym,Cryst%rprimd,my_shiftk,Cryst%symafm,Cryst%symrel,my_vacuum,wtk)

 ! Store quantities that cannot be easily (and safely) calculated if we only know the IBZ.
 Kmesh%nshift   = my_nshiftk
 Kmesh%kptrlatt = kptrlatt

 ! Call the main creation method to get the tables tabo, tabi, tabp, umklp...
 ! init_kmesh will reconstruct the BZ from kibz but pruning the k-points not in ref_bz
 ! TODO: solve problem with timrev
 my_break_symmetry=.FALSE.; if (PRESENT(break_symmetry)) my_break_symmetry=break_symmetry
 call init_kmesh(Kmesh,Cryst,nkibz,kibz,kptopt,ref_bz=ref_kbz,break_symmetry=my_break_symmetry)

 ABI_DEALLOCATE(ref_kbz)
 ABI_DEALLOCATE(kibz)
 ABI_DEALLOCATE(wtk)

 ABI_ALLOCATE(Kmesh%shift,(3,my_nshiftk))
 Kmesh%shift=my_shiftk(:,1:my_nshiftk)
                                      
 ABI_DEALLOCATE(my_shiftk)

 if (PRESENT(want_tetra)) then
   if (want_tetra==1) call init_tetra_tab(Kmesh,Cryst)
 end if

 DBG_EXIT("COLL")

end subroutine make_mesh
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/init_tetra_tab
!! NAME
!! init_tetra_tab
!!
!! FUNCTION
!!  Initialize tables for tetrahedron integration scheme. 
!!
!! INPUTS
!!  Cryst<Crystal_structure> = Info on unit cell and its symmetries
!!
!! SIDE EFFECTS
!!  Kmesh<bz_mesh_type>=Datatype gathering information on the k point sampling.
!!   In output tables for tetrahedron integrations scheme are allocated and initialized
!!
!! PARENTS
!!      m_bz_mesh
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine init_tetra_tab(Kmesh,Cryst)

 use defs_basis
 use m_tetrahedron

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_tetra_tab'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(inout) :: Kmesh

!Local variables-------------------------------
!scalars
 logical :: ltest
 integer :: ierr
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: my_fullbz(:,:)
 integer,allocatable :: bz2ibz(:)
 real(dp) :: qlatt(3,3),rlatt(3,3)
 character(len=80) :: errstr

! *************************************************************************

 ! ============================
 ! ==== Tetrahedron method ====
 ! ============================

 if (associated(Kmesh%tetrahedra)) then
   if (Kmesh%tetrahedra%ntetra/=0) then 
     msg = "Tetrahedrons already calculated, returning"
     MSG_WARNING(msg)
     RETURN
   end if
 end if

 if (ALL(Kmesh%kptrlatt==NONE_KPTRLATT))  then
   MSG_WARNING("%kptrlatt is not available, returning")
   RETURN
 end if
 !
 ! * convert kptrlatt to double and invert. 
 rlatt=DBLE(Kmesh%kptrlatt)
 call matr3inv(rlatt,qlatt)  ! qlatt refers to the shortest qpt vectors.
 !  
 ! === Make full kpoint grid and get equivalence to irred kpoints ===
 ! * Note: This routines scales badly wrt Kmesh%nbz 
 ! TODO should be rewritten, pass timrev and speed up by looping on shells.
 write(msg,'(a,i0)')" calling get_full_kgrid with nbz= ",Kmesh%nbz
 call wrtout(std_out,msg,"COLL")
 call pclock(0)

 ABI_CHECK(Kmesh%kptopt==1,"get_full_kgrid assumes kptopt==1")

 ! one can reuse %tabo instead of bz2ibz but have to make sure that ordering does not matter.
 ABI_ALLOCATE(bz2ibz,(Kmesh%nbz))
 ABI_ALLOCATE(my_fullbz,(3,Kmesh%nbz))

 call get_full_kgrid(bz2ibz,qlatt,Kmesh%ibz,my_fullbz,Kmesh%kptrlatt,Kmesh%nibz,Kmesh%nbz,&
&  Kmesh%nshift,Cryst%nsym,Kmesh%shift,Cryst%symrel)

 write(std_out,*) 'after get_full_kgrid'
 call pclock(9999)

 ltest = ALL(ABS(my_fullbz-Kmesh%bz)<tol6)
 ABI_CHECK(ltest,"Mismatch between my_fullbz and Kmesh%bz")
 ABI_DEALLOCATE(my_fullbz)
 
 call init_tetra(bz2ibz, Cryst%gprimd, qlatt, Kmesh%bz, Kmesh%nbz,&
&                Kmesh%tetrahedra, ierr, errstr)
 ABI_CHECK(ierr==0,errstr)

 write(msg,'(a,i0)')' Number of irreducible tetrahedrons= ',Kmesh%tetrahedra%ntetra
 call wrtout(std_out,msg,"COLL")
 ABI_DEALLOCATE(bz2ibz)

end subroutine init_tetra_tab
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/identk
!! NAME
!! identk
!!
!! FUNCTION
!! Identify k-points in the whole BZ starting from the IBZ. 
!! Generate also symmetry tables relating the BZ to the IBZ.
!!
!! INPUTS
!!  kibz(3,nkibz)=Coordinates of k-points in the IBZ.
!!  nkibz=Number of k points in IBZ.
!!  nkbzmx=Maximum number of k points in BZ.
!!  nsym=Number of symmetry operations.
!!  timrev=2 if time reversal symmetry can be used; 1 otherwise.
!!  symrec(3,3,nsym)=Symmetry operation matrices in reciprocal space.
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations.
!!  [ref_bz(:,:)]= Reference set of points in the full Brillouin zone.
!!
!! OUTPUT
!!  kbz(3,nkbzmx)= k-points in whole BZ
!!  ktab(nkbzmx)= table giving for each k-point in the BZ (array kbz), 
!!   the corresponding irreducible point in the array (kibz)
!!   k_BZ= (IS) kIBZ where S is one of the symrec operations and I is the inversion or the identity
!!    where k_BZ = (IS) k_IBZ and S = \transpose R^{-1} 
!!  ktabi(nkbzmx)= for each k-point in the BZ defines whether inversion has to be 
!!   considered in the relation k_BZ=(IS) k_IBZ (1 => only S; -1 => -S)  
!!  ktabo(nkbzmx)= the symmetry operation S that takes k_IBZ to each k_BZ
!!  nkbz= no. of k-points in the whole BZ
!!  wtk(nkibz)= weight for each k-point in IBZ for symmetric quantities:
!!              no. of distinct ks in whole BZ/(timrev*nsym)
!!
!! NOTES
!!  The logic of the routine relies on the assumption that kibz really represent an irreducible set. 
!!  If symmetrical points are present in the input list, indeed, some the output weights will turn out to be zero.
!!  An initial check is done at the beginning of the routine to trap this possible error.
!!
!! PARENTS
!!      m_bz_mesh,rdm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine identk(kibz,nkibz,nkbzmx,nsym,timrev,symrec,symafm,kbz,ktab,ktabi,ktabo,nkbz,wtk,ref_bz)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'identk'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkbzmx,nkibz,nsym,timrev
 integer,intent(out) :: nkbz
!arrays
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym)
 integer,intent(out) :: ktab(nkbzmx),ktabi(nkbzmx),ktabo(nkbzmx)
 real(dp),intent(in) :: kibz(3,nkibz)
 real(dp),intent(out) :: kbz(3,nkbzmx),wtk(nkibz)
 real(dp),optional,intent(in) :: ref_bz(:,:)

!Local variables ------------------------------
!scalars
 integer :: ik1,ik2,ikbz,ikibz,iold,isym,itim
 integer :: ikref,nkref,isym_swp,itim_swp,ikibz_swp
 logical :: is_irred_set
 logical :: found,ltest
 character(len=500) :: msg
!arrays
 integer :: g0(3)
 real(dp) :: knew(3),k1(3),k2(3)
 real(dp) :: kref(3),kbz_swp(3)

! *************************************************************************

 DBG_ENTER("COLL")
 !
 ! === Check whether kibz really forms an irreducible set ===
 is_irred_set=.TRUE.
 do ik1=1,nkibz-1
   k1=kibz(:,ik1)
   do ik2=ik1+1,nkibz
     k2=kibz(:,ik2)

     do itim=1,timrev
       do isym=1,nsym
         if (symafm(isym)==-1) CYCLE
         knew = (3-2*itim) * MATMUL(symrec(:,:,isym),k2)
         if (isamek(k1,knew,g0)) then
           is_irred_set=.FALSE.
           write(msg,'(2(a,3f8.4),2(a,i2))')&
&            ' k1 = ',k1,' is symmetrical of k2 = ',k2,' through sym = ',isym,' itim = ',itim
           MSG_WARNING(msg)
         end if
       end do
     end do

   end do
 end do

 !call klist_isirred(nkibz,kibz,Cryst,nimg)

 if (.not.is_irred_set) then
   msg = "Input array kibz does not constitute an irreducible set."
   MSG_WARNING(msg) 
 end if
 !
 ! === Loop over k-points in IBZ ===
 ! * Start with zero no. of k-points found.
 nkbz=0 
 do ikibz=1,nkibz
   wtk(ikibz) = zero

   ! === Loop over time-reversal I and symmetry operations S  ===
   ! * Use spatial inversion instead of time reversal whenever possible.
   do itim=1,timrev
     do isym=1,nsym
      if (symafm(isym)==-1) CYCLE
      !
      ! * Form IS k
      knew=(3-2*itim)*MATMUL(symrec(:,:,isym),kibz(:,ikibz))
      !
      ! * Check whether it has already been found (to within a RL vector).
      iold=0
      do ikbz=1,nkbz
        if (isamek(knew,kbz(:,ikbz),g0)) iold=iold+1
      end do
      !
      ! * If not yet found add to kbz and increase the weight.
      if (iold==0) then
        nkbz=nkbz+1
        wtk(ikibz)=wtk(ikibz)+one
        if (nkbz>nkbzmx) then
          write(msg,'(a,i0,a)')' nkbzmx too small, nkbzmx = ',nkbzmx,', increase nkbzmx !'
          MSG_BUG(msg)
        end if
        kbz(:,nkbz) = knew(:)
        ktab (nkbz) = ikibz
        ktabo(nkbz) = isym
        ktabi(nkbz) = 3-2*itim
      end if
     end do 
   end do 

 end do !ikibz

 if (PRESENT(ref_bz)) then 
   call wrtout(std_out," Pruning the k-points not in ref_bz then reordering tables","COLL")

   nkref = SIZE(ref_bz,DIM=2)
   ltest = (nkref>=nkbz.and.nkref<=nkbzmx)
   if (.not.ltest) then
     write(msg,'(3(a,i0))')" Wrong value for nkref: nkref= ",nkref," nkbz= ",nkbz," nkbzmx =",nkbzmx
     MSG_WARNING(msg)
   end if

   do ikref=1,nkref
     kref = ref_bz(:,ikref)
     found=.FALSE.

     do ikbz=1,nkbz ! Loop on the set of BZ points found above.
       if (isequalk(kref,kbz(:,ikbz))) then ! Swap indeces.
         kbz_swp   = kbz(:,ikref) 
         ikibz_swp = ktab (ikref) 
         isym_swp  = ktabo(ikref)
         itim_swp  = ktabi(ikref) 

         kbz(:,ikref) = kref
         ktab (ikref) = ktab (ikbz) 
         ktabo(ikref) = ktabo(ikbz)
         ktabi(ikref) = ktabi(ikbz) 

         kbz(:,ikbz) = kbz_swp 
         ktab (ikbz) = ikibz_swp
         ktabo(ikbz) = isym_swp
         ktabi(ikbz) = itim_swp

         found=.TRUE.; EXIT
       end if
     end do

     if (.not.found) then 
       write(msg,'(a,3es16.8)')" One of the k-point in ref_bz is not a symmetrical image of the IBZ: ",kref
       MSG_ERROR(msg)
     end if
   end do
   !
   ! Change nkbz to nkref, then get the new weights.
   nkbz=nkref; wtk=zero
   do ikref=1,nkref
     ikibz = ktab(ikref)
     wtk(ikibz) = wtk(ikibz) + 1
   end do
 end if ! PRESENT(ref_bz)
 !
 ! * Weights are normalized to 1.
 wtk = wtk/SUM(wtk)

 DBG_EXIT("COLL")

end subroutine identk
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/get_ng0sh
!! NAME
!! get_ng0sh
!!
!! FUNCTION
!!  Given two lists of k-points, kbz1 and kbz2, calculate any possible difference k1-k2.
!!  For each difference, find the umklapp g0 vector and the point k3 in the array kfold 
!!  such as k1-k2 = k3 + G0.
!!  The optimal value of G0 shells is returned, namely the smallest box around Gamma 
!!  which suffices to treat all possible umklapp processes.
!!
!! INPUTS
!!  nk1, nk2=Number of points in the arrays kbz1, kbz2.
!!  kbz1(3,nk1)=Reduced coordinates of the first set of points.
!!  kbz2(3,nk2)=Reduced coordinates of the second set of points.
!!  nkfold=Number of points in the array kfold.
!!  kfold(3,nkfol)=Reduced coordinated of the points in the BZ.
!!  gmet(3,3)=Reciprocal space metric ($\textrm{bohr}^{-2}$).
!!  tolq0=Tolerance below which a q-point is treated as zero.
!!  mg0sh=Integer defining the Max number of shells to be considered.
!!
!! OUTPUT
!!  opt_ng0(3)=Minimal reduced components of the G0 vectors to account for umklapps.
!!
!! PARENTS
!!      setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine get_ng0sh(nk1,kbz1,nk2,kbz2,nkfold,kfold,gmet,tolq0,mg0sh,opt_ng0)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_ng0sh'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mg0sh,nk1,nk2,nkfold
 real(dp),intent(in) :: tolq0
!arrays
 integer,intent(out) :: opt_ng0(3)
 real(dp),intent(in) :: gmet(3,3),kbz1(3,nk1),kbz2(3,nk2),kfold(3,nkfold)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ig1,ig2,ig3,ii,ikf,itr,ng0
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gtrial(3)
 integer,allocatable :: iperm(:),g0(:,:) 
 real(dp) :: k1mk2(3)
 real(dp),allocatable :: norm(:)

!************************************************************************
 
 ng0=(2*mg0sh+1)**3
 ABI_ALLOCATE(g0,(3,ng0))
 ABI_ALLOCATE(norm,(ng0))
 ABI_ALLOCATE(iperm,(ng0))

 ii=1
 do ig3=-mg0sh,mg0sh
   do ig2=-mg0sh,mg0sh
     do ig1=-mg0sh,mg0sh
       g0(1,ii)=ig1
       g0(2,ii)=ig2
       g0(3,ii)=ig3
       norm(ii)=normv(g0(:,ii),gmet,'G')
       iperm(ii)=ii
       ii=ii+1
     end do
   end do
 end do
 !
 ! === Order g0-vectors in order of increasing module ===
 ! * Should speed up the search in the while statement below.
 call sort_dp(ng0,norm,iperm,tol14)

 opt_ng0(:)=0 
 do i2=1,nk2
   ! This is used in case of screening calculation.
   ! If q is small treat it as zero. In this case, indeed, 
   ! we use q=0 to calculate the oscillator matrix elements. 
   if (is_zero(kbz2(:,i2),tolq0)) CYCLE
   do i1=1,nk1 
     k1mk2(:)=kbz1(:,i1)-kbz2(:,i2)
     itr=1; found=.FALSE.
     do while (itr<=ng0.and..not.found)
       do ikf=1,nkfold
         if (ALL(ABS(k1mk2(:)-kfold(:,ikf)-g0(:,iperm(itr))) < TOL_KDIFF)) then 
           found=.TRUE. 
           gtrial(:)=ABS(g0(:,iperm(itr)))
           opt_ng0(1)=MAX(opt_ng0(1),gtrial(1))
           opt_ng0(2)=MAX(opt_ng0(2),gtrial(2))
           opt_ng0(3)=MAX(opt_ng0(3),gtrial(3)); EXIT
         end if
       end do
       itr=itr+1
     end do
     if (.not.found) then 
       write(msg,'(a,i5,3a,i2,2(2a,i4,3es16.8),a)')&
&        ' Not able to found umklapp G0 vector among ',ng0,' vectors',ch10,&
&        ' Increase mg0sh such as k1-k2 = kf+G0, present value = ',mg0sh,ch10,&
&        ' point1 = ',i1,kbz1(:,i1),ch10,&
&        ' point2 = ',i2,kbz2(:,i2),ch10
       MSG_ERROR(msg)
     end if
   end do
 end do 

 write(msg,'(a,3i2)')' optimal value for ng0sh = ',opt_ng0(:)
 call wrtout(std_out,msg,'COLL')

 ABI_DEALLOCATE(g0)
 ABI_DEALLOCATE(norm)
 ABI_DEALLOCATE(iperm)

end subroutine get_ng0sh
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/make_path
!! NAME
!! make_path
!!
!! FUNCTION
!!  Create a normalized path given the extrema.
!!
!! INPUTS
!!  nbounds=Number of extrema defining the path.
!!  bounds(3,nbounds)=The points defining the path in reduced coordinates.
!!  met(3,3)=Metric matrix. 
!!  space='R' for real space, G for reciprocal space.
!!  ndiv_small=Number of divisions to be used for the smallest segment.
!!
!! OUTPUT
!!  npt_tot=Total number of points in the normalized circuit.
!!  ndiv(nbounds-1)=Number of division for each segment
!!  See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  In input path(:,:) is a NULL pointer.
!!  It is allocated inside the routine. When the subroutine returns, path(3,npt_tot) will
!!  contain the path in reduced coordinates.
!! 
!! PARENTS
!!      exc_interp_ham,m_gamma,m_gwannier,mkph_linwid,mkphbs,outnesting
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine make_path(nbounds,bounds,met,space,ndiv_small,ndiv,npt_tot,path)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_path'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbounds,ndiv_small
 integer,intent(out) :: npt_tot
 character(len=1),intent(in) :: space
!arrays
 integer,intent(out) :: ndiv(nbounds-1)
 real(dp),intent(in) :: bounds(3,nbounds),met(3,3)
 real(dp),pointer :: path(:,:)

!Local variables-------------------------------
!scalars
 integer :: idx,ii,jp
 real(dp) :: nfact
 character(len=500) :: msg
!arrays
 real(dp) :: diff(3),lng(nbounds-1)

! *************************************************************************

 ABI_CHECK(ndiv_small>0,'ndiv_small <=0')

 do ii=1,nbounds-1 
   diff(:)=bounds(:,ii+1)-bounds(:,ii)
   lng(ii) = normv(diff,met,space)
 end do

 ! Avoid division by zero if any k(:,i+1)=k(:,i).
 nfact=MINVAL(lng)
 if (ABS(nfact)<tol6) then 
   write(msg,'(3a)')&
&    ' Found two equal consecutive points in the path ',ch10,&
&    ' This is not allowed, please modify the path in your input file'
   MSG_ERROR(msg)
 end if

 nfact=nfact/ndiv_small
 ndiv(:)=NINT(lng(:)/nfact) 
 npt_tot=SUM(ndiv)+1 !1 for the first point

 write(msg,'(2a,i6,2a)')ch10,&
&  ' Total number of points in the path : ',npt_tot,ch10,&
&  ' Number of divisions for each segment of the normalized path : '
 call wrtout(std_out,msg,'COLL') 

 do ii=1,nbounds-1
   write(msg,'(2(3f8.5,a),i5,a)')bounds(:,ii),' ==> ',bounds(:,ii+1),' ( ndiv : ',ndiv(ii),' )' 
   call wrtout(std_out,msg,'COLL')
 end do 
 write(msg,'(a)')ch10
 call wrtout(std_out,msg,'COLL') 

 ! Allocate and construct the path.
 ABI_ALLOCATE(path,(3,npt_tot))

 call wrtout(std_out,' Normalized Path: ','COLL')
 idx=0
 do ii=1,nbounds-1
   do jp=1,ndiv(ii)
     idx=idx+1
     path(:,idx)=bounds(:,ii)+(jp-1)*(bounds(:,ii+1)-bounds(:,ii))/ndiv(ii)
     write(msg,'(i4,4x,3(f8.5,1x))')idx,path(:,idx)
     call wrtout(std_out,msg,'COLL')
   end do 
 end do 

 path(:,npt_tot)=bounds(:,nbounds)
 write(msg,'(i4,4x,3(f8.5,1x))')npt_tot,path(:,npt_tot)
 call wrtout(std_out,msg,'COLL')

end subroutine make_path
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/find_qmesh
!! NAME
!! find_qmesh
!!
!! FUNCTION
!!  Find the q-mesh defined as all the possible differences between k-points
!!  Find the irreducible q-points using a special treatment for the Gamma point.
!!  Then call setup_kmesh to initialize the Qmesh datatype
!!
!! INPUTS
!!  Cryst<crystal_structure>=datatype gathering info on the unit cell and symmetries
!!    %nsym=number of symmetry operations
!!    %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  Kmesh<bz_mesh_type>=datatype gathering information on the k-mesh
!!
!! OUTPUT
!!  Qmesh<bz_mesh_type>=datatype gathering information on the q point sampling. 
!!
!! PARENTS
!!      m_ebands,m_screening,m_shexc,mrgscr,setup_bse,setup_screening
!!      setup_sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine find_qmesh(Qmesh,Cryst,Kmesh)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_qmesh'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(BZ_mesh_type),intent(inout) :: Qmesh
 type(Crystal_structure),intent(in) :: Cryst

!Local variables-------------------------------
!scalars
 integer :: nqibz,kptopt
!arrays
 real(dp),allocatable :: qibz(:,:)

! *************************************************************************

 !@BZ_mesh_type
 call nullify_BZ_mesh(Qmesh)

 ! * Find the number of q-points such that q = k1-k2.
 call findnq(Kmesh%nbz,Kmesh%bz,Cryst%nsym,Cryst%symrec,Cryst%symafm,nqibz,Cryst%timrev) 
 !
 ! * Find the coordinates of the q-points in the IBZ.
 ABI_ALLOCATE(qibz,(3,nqibz))
 call findq(Kmesh%nbz,Kmesh%bz,Cryst%nsym,Cryst%symrec,Cryst%symafm,Cryst%gprimd,nqibz,qibz,Cryst%timrev)
 !
 ! * Create the Qmesh object starting from the IBZ.
 kptopt = Kmesh%kptopt
 call init_kmesh(Qmesh,Cryst,nqibz,qibz,kptopt)

 ABI_DEALLOCATE(qibz)

end subroutine find_qmesh
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/findnq
!! NAME
!! findnq
!!
!! FUNCTION
!! Identify the number of q-points in the IBZ by which the k-points in BZ differ
!! (count the q points in the k-point difference set)
!!
!! INPUTS
!!  nkbz=number of k points in Brillouin zone
!!  kbz(3,nkbz)=coordinates of k points in BZ
!!  nsym=number of symmetry operations
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  symafm(nsym)=-1 if AFM symmetry.
!!  timrev=2 if time-reversal symmetry is used, 1 otherwise
!!
!! OUTPUT
!!  nqibz=number of q points
!!
!! PARENTS
!!      m_bz_mesh,rdm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine findnq(nkbz,kbz,nsym,symrec,symafm,nqibz,timrev)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'findnq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: timrev,nkbz,nsym
 integer,intent(out) :: nqibz
!arrays
 integer,intent(in) :: symrec(3,3,nsym),symafm(nsym)
 real(dp),intent(in) :: kbz(3,nkbz)

!Local variables ------------------------------
!scalars
 integer :: ifound,ik,isym,iq,istat,memory_exhausted,nqall,nqallm,itim
!arrays
 integer :: g0(3)
 real(dp) :: qposs(3),qrot(3)
 real(dp),allocatable :: qall(:,:)
!************************************************************************

 ! === Infinite do-loop to be able to allocate sufficient memory ===
 nqallm=1000
 do 
   memory_exhausted=0
   ABI_ALLOCATE(qall,(3,nqallm))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,'out-of-memory qall')
   nqall=0
   ! === Loop over all k-points in BZ, forming k-k1 ===
   do ik=1,nkbz
     qposs(:)=kbz(:,ik)-kbz(:,1)
     !
     ! * Check whether this q (or its equivalent) has already been found within a reciprocal lattice vector.
     ! * Use spatial inversion instead of time reversal whenever possible.
     ifound=0
     do iq=1,nqall
       do itim=1,timrev
         do isym=1,nsym
           if (symafm(isym)==-1) CYCLE
           qrot = (3-2*itim) * MATMUL(symrec(:,:,isym),qall(:,iq))
           if (isamek(qrot,qposs,g0)) ifound=ifound+1
         end do
       end do
     end do

     if (ifound==0) then
       nqall=nqall+1
       !
       ! === If not yet found, check that the allocation is big enough ===
       if (nqall>nqallm) then
         memory_exhausted=1 
         ABI_DEALLOCATE(qall)
         nqallm=nqallm*2    ; EXIT ! Exit the do ik=1 loop
       end if
       ! * Add it to the list.
       qall(:,nqall)=qposs(:)
     end if
   end do

   if (memory_exhausted==0) EXIT
 end do !infinite loop

 ABI_DEALLOCATE(qall)
 nqibz=nqall

end subroutine findnq
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/findq
!! NAME
!! findq
!!
!! FUNCTION
!! Identify the q-points by which the k-points in BZ differ
!!
!! INPUTS
!!  nkbz=number of k points in Brillouin zone
!!  kbz(3,nkbz)=coordinates of k points in BZ
!!  nsym=number of symmetry operations
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  symafm(nsym)=-1 is symmetry is AFM, +1 otherwise.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  nqibz=number of q points in the IBZ by which k points differ (computed in findnq)
!!  timrev=2 if time-reversal symmetry is used, 1 otherwise
!!
!! OUTPUT
!!  qibz(3,nqibz)=coordinates of q points by which k points differ
!!
!! PARENTS
!!      m_bz_mesh,rdm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE


subroutine findq(nkbz,kbz,nsym,symrec,symafm,gprimd,nqibz,qibz,timrev)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'findq'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkbz,nqibz,nsym,timrev
!arrays
 integer,intent(in) :: symrec(3,3,nsym),symafm(nsym)
 real(dp),intent(in) :: gprimd(3,3),kbz(3,nkbz)
 real(dp),intent(out) :: qibz(3,nqibz)

!Local variables ------------------------------
!scalars
 integer :: ii,ik,iq,iqp,isym,itim
 real(dp) :: shift1
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: g0(3)
 real(dp) :: gmet(3,3),qposs(3),qrot(3)

!************************************************************************

 ! Compute reciprocal space metrics
 do ii=1,3
   gmet(ii,:)=gprimd(1,ii)*gprimd(1,:)+&
&             gprimd(2,ii)*gprimd(2,:)+&
&             gprimd(3,ii)*gprimd(3,:)
 end do
 !
 ! === Loop over k-points in BZ, form k-k1 and translate in first BZ ===
 ! iq is the no. of q-points found, zero at the beginning
 iq=0
 do ik=1,nkbz
   qposs(:)=kbz(:,ik)-kbz(:,1)
   ! === Check whether this q (or its equivalent) has already been found ===
   ! * Use spatial inversion instead of time reversal whenever possible.
   found=.FALSE.
   do iqp=1,iq
     do itim=1,timrev
       do isym=1,nsym
        if (symafm(isym)==-1) CYCLE
        qrot = (3-2*itim) * MATMUL(symrec(:,:,isym),qibz(:,iqp))
        if (isamek(qrot,qposs,g0)) found=.TRUE.
       end do
     end do
   end do
   if (.not.found) then
     iq=iq+1
     if (iq>nqibz) then 
       write(msg,'(a,i5)')' iq > nqibz= ',nqibz
       MSG_BUG(msg)
     end if 
     qibz(:,iq)=qposs(:)
   end if
 end do

 if (iq/=nqibz) then 
   write(msg,'(2(a,i5))')' iq= ',iq,'/= nqibz= ',nqibz
   MSG_BUG(msg)
 end if 
 !
 ! * Translate q-points to 1st BZ in the interval [-1/2,1/2[
 do iq=1,nqibz
   do ii=1,3
     call wrap2_pmhalf(qibz(ii,iq),qibz(ii,iq),shift1)
   end do
 end do 

end subroutine findq
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/findqg0
!! NAME
!! findqg0
!!
!! FUNCTION
!! Identify q + g0 = k - kp
!!
!! INPUTS
!!  kmkp(3)= k - kp input vector
!!  nqbz=number of q points
!!  qbz(3,nqbz)=coordinates of q points
!!  mG0(3)= For each reduced direction gives the max G0 component to account for umklapp processes
!!
!! OUTPUT
!!  iq=index of q vector in BZ table
!!  g0(3)=reciprocal space vector, to be used in igfft
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cohsex_me,exc_build_block,gw_tools,rdm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE


subroutine findqg0(iq,g0,kmkp,nqbz,qbz,mG0)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'findqg0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqbz
 integer,intent(out) :: iq
!arrays
 integer,intent(in) :: mG0(3)
 integer,intent(out) :: g0(3)
 real(dp),intent(in) :: kmkp(3),qbz(3,nqbz)

!Local variables-------------------------------
!FIXME if I use 1.0d-4 the jobs crash, should understand why
!scalars
 integer :: iqbz,jg01,jg02,jg03,nq0found
 real(dp) :: tolq0=1.0D-3
 character(len=500) :: msg
!arrays
 real(dp) :: qpg0(3)

! *************************************************************************

 iq=0
 if (ALL(ABS(kmkp(:))<EPSILON(one))) then ! Find q close to 0 ===
   nq0found=0
   do iqbz=1,nqbz
     if (ALL(ABS(qbz(:,iqbz))<tolq0)) then 
       iq=iqbz
       nq0found=nq0found+1
     end if 
   end do
   if (iq==0) then 
     msg= 'Wrong list of q-points: q=0 not present.'
     MSG_BUG(msg)
   end if
   if (nq0found/=1) then 
     write(msg,'(3a)')&
&      ' The input q-list contains more than one small q-point',ch10,&
&      ' ACTION: reduce the value of tolq0 '
     MSG_BUG(msg)
   end if 
   g0(:)=0; RETURN

 else ! q is not zero, find q such as k-kp=q+G0.

   do iqbz=1,nqbz
     do jg01=-mG0(1),mG0(1)
       do jg02=-mG0(2),mG0(2)
         do jg03=-mG0(3),mG0(3)
           ! * Form q-G0 and check if it is the one.
           qpg0(1)= qbz(1,iqbz)+jg01
           qpg0(2)= qbz(2,iqbz)+jg02
           qpg0(3)= qbz(3,iqbz)+jg03
           if (ALL(ABS(qpg0(:)-kmkp(:))<TOL_KDIFF)) then
             if (iq/=0) then
               write(std_out,*)iqbz,qbz(:,iqbz),jg01,jg02,jg03
               write(std_out,*)iq,qbz(:,iq),g0
               msg='Found duplicated q+g0.'
               MSG_ERROR(msg)
             end if
             iq=iqbz
             g0(1)=jg01
             g0(2)=jg02
             g0(3)=jg03
           end if
         end do
       end do
     end do
   end do
   if (iq==0) then
     write(msg,'(a,3f9.5)')' q = k-kp+G0 not found. kmkp = ',kmkp
     MSG_ERROR(msg)
   end if
 end if

end subroutine findqg0
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/setup_little_group
!! NAME
!! setup_little_group
!!
!! FUNCTION
!!  Finds symmetry operations belonging to the little group associated to an external 
!!  point ext_pt and fills symmetry tables.
!!
!! INPUTS
!! ext_pt(3)= External point in the Brillouin zone in reduce coordinated 
!! Kmesh<BZ_mesh_type>
!!   %nbz=number of points in the full BZ
!!   %kbz(3,nbz)=points in the full Brillouin Zon
!! Cryst<Crystal_structure>= Info on symmetries and unit cell
!!   %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!   %nsym=number of symmetry operations in the space group
!!   %timrev=if 2 time-reversal can be used; 1 otherwise
!!   %gmet(3,3)=reciprocal space metric (bohr**-2).
!! npwvec=number of G vectors
!! gvec(3,npwvec) coordinates of G vectors
!! npwe=If greater than 0, the index of G-Go in the gvec(:,1:npwvec) array will be calculated
!!  and stored in %igmG0 for each symmetry preserving the external q. Note that G is one of the npwe vectors. 
!! use_umklp=flag to include umklapp G0 vectors in the definition of the little group (0:n0,1:yes)
!!
!! OUTPUT
!! Ltg% <little_group_datatype>.
!!  %ibzq(nbz)= 1 if the kpoint belongs to the IBZ defined by ext_pt, 0 otherwise
!!  %bz2ibz(nbz)= sequential index of the point in the IBZ defined by ext_pt
!!  %ibz2bz(nibz_Ltg) For each nibz_Ltg the correspondind index in the BZ array
!!  %igmG0(npwepstimrev,nsym)= index of the uklapp vector G_o in the FFT array
!!  %flag_umklp(timrev,nsym)= flag for umklapp processes 
!!    1 if operation (IS) requires a G_o to preserve ext_pt, 0 otherwise 
!!  %tab(nbz)=table giving, for each k-point in the BZ (kBZ), the corresponding 
!!   irreducible point (kIBZ) in the irreducible zone defined by the little group of ext_pt,
!!   i.e kBZ= (IS) kIBZ where I is either the inversion or the identity and S is an
!!   operation in the little group defined by ext_pt
!!  %tabo(nbz)=the symmetry operation S in the little group that takes kIBZ to each kBZ
!!  %tabi(nbz)= defines whether inversion has to be considered in the 
!!   relation kBZ=(IS) kIBZ (1 => only S; -1 => -S)  
!!  %preserve(2,nsym)= 1 if ISq=q, 0 otherwise, the first index is for the identity or the time reversal symmetry, 
!!  %wtksym(2,nsym,nbz)= for each kpoint is equal to 1 if the symmetry operation (with or without time reversal)  
!!   must be considered in the calculation of \chi_o, 0 otherwise  
!!
!! PARENTS
!!      cchi0q0_intraband,setup_screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine setup_little_group(ext_pt,Kmesh,Cryst,use_umklp,Ltg,npwe,gvec)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_little_group'
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwe,use_umklp
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Little_group),intent(inout) :: Ltg
!arrays
 !integer,optional,intent(in) :: gvec(3,npwvec)
 integer,optional,intent(in) :: gvec(:,:)
 real(dp),intent(in) :: ext_pt(3)

!Local variables-------------------------------
!scalars
 integer :: dummy_timrev,enough,idx,ige,igpw,ik,ind,iold,iout,isym,itest,itim
 integer :: nbz,nkibzq,nsym,nsym_Ltg,ntest,timrev,ierr,npwvec
 real(dp) :: G0len,kin,mG0len,max_kin
 logical :: found,found_identity,ltest,use_antiferro
 character(len=500) :: msg
!arrays
 integer :: g0(3),gg(3),gmG0(3),identity(3,3),nop(Cryst%timrev),nopg0(2)
 integer :: symxpt(4,2,Cryst%nsym)
 integer,allocatable :: indkpt1(:),symafm_ltg(:),symrec_Ltg(:,:,:)
 integer,pointer :: symafm(:),symrec(:,:,:)
 real(dp) :: knew(3)
 real(dp),allocatable :: ktest(:,:),wtk(:),wtk_folded(:)
 real(dp),pointer :: kbz(:,:)

!************************************************************************

 DBG_ENTER("COLL")
 !
 ! !@little_group
 ! === Initial check ====
 ltest=(Cryst%timrev==1.or.Cryst%timrev==2)
 ABI_CHECK(ltest,'Wrong value for timrev')
 !
 ! === Get useful data ===
 nsym          =  Cryst%nsym
 !timrev       = 1
 timrev        =  Cryst%timrev
 symrec        => Cryst%symrec
 symafm        => Cryst%symafm
 use_antiferro =  Cryst%use_antiferro

 nbz =  Kmesh%nbz 
 kbz => Kmesh%bz(1:3,1:nbz)
 !
 ! === Destroy structure if it already exists ===
 call destroy_little_group(Ltg)
 !
 ! === Store dimensions and useful info ===
 Ltg%nsym_sg  =nsym
 Ltg%timrev   =timrev
 Ltg%nbz      =nbz
 !Ltg%use_umklp=use_umklp ! 1 if umklapp processes are used
 Ltg%ext_pt(:)=ext_pt(:)

 ABI_ALLOCATE(Ltg%G0,(3,timrev,nsym))
 ABI_ALLOCATE(Ltg%ibzq,(nbz))
 ABI_ALLOCATE(Ltg%bz2ibz,(nbz))
 ABI_ALLOCATE(Ltg%preserve,(timrev,nsym))
 ABI_ALLOCATE(Ltg%wtksym,(timrev,nsym,nbz))
 ABI_ALLOCATE(Ltg%tab,(nbz))
 ABI_ALLOCATE(Ltg%tabi,(nbz))
 ABI_ALLOCATE(Ltg%tabo,(nbz))
 ABI_ALLOCATE(Ltg%flag_umklp,(timrev,nsym))
 !
 ! In the old GW implementation we were removing symmetries related by time-reversal and 
 ! sometimes it happened that only the inversion was reported in the KSS file (see outkss.F90).
 identity(:,:)=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/)) ; found_identity=.FALSE. 
 do isym=1,nsym
   if (ALL(symrec(:,:,isym)==identity)) then 
    found_identity=.TRUE. ; EXIT
   end if
 end do 
 if (.not.found_identity) then 
   write(msg,'(5a)')&
&    ' Only the inversion was found in the set of symmetries read from the KSS file ',ch10,&
&    ' Likely you are using a KSS file generated with an old version of Abinit, ',ch10,&
&    ' To run a GW calculation with an old KSS file, use version < 5.5 '
   MSG_ERROR(msg)
 end if 
 !
 ! === Find operations in the little group as well as umklapp vectors G0 === 
 call symq3(nsym,ext_pt,symxpt,symrec,dummy_timrev,0) 

 Ltg%preserve(:,:)=0; Ltg%g0(:,:,:)=0; Ltg%flag_umklp(:,:)=0; mG0len=zero

 do itim=1,timrev
   do isym=1,nsym
     if (symafm(isym)==-1) CYCLE

     if (symxpt(4,itim,isym)==1) then  !\pm Sq = q+g0
       if (ANY(symxpt(1:3,itim,isym)/=0).and.use_umklp==0) CYCLE ! Exclude non zero G0 vectors
       Ltg%preserve(itim,isym)=1
       g0(:)=symxpt(1:3,itim,isym); Ltg%g0(:,itim,isym)=g0(:)
       if (ANY(Ltg%g0(:,itim,isym)/=0)) Ltg%flag_umklp(itim,isym)=1
       ! Max radius to be considered to include all G0s
       G0len=normv(g0,Cryst%gmet,'G') 
       mG0len=MAX(mG0len,G0len) 
     end if
   end do
 end do 

 nop(:)=0; nopg0(:)=0
 do itim=1,timrev
   nop  (itim)=SUM(Ltg%preserve  (itim,:))
   nopg0(itim)=SUM(Ltg%flag_umklp(itim,:))
 end do
 nsym_Ltg=SUM(nop(:))
 !
 ! === Store little group operations, include time-reversal if present ===
 Ltg%nsym_Ltg=nsym_Ltg
 ABI_ALLOCATE(symrec_Ltg,(3,3,Ltg%nsym_Ltg))

 ind=1
 do itim=1,timrev
   do isym=1,nsym
     if (Ltg%preserve(itim,isym)==1) then  
      if (itim==1) symrec_Ltg(:,:,ind)= symrec(:,:,isym)
      if (itim==2) symrec_Ltg(:,:,ind)=-symrec(:,:,isym)
      ind=ind+1
     end if
   end do 
 end do
 !
 ! === Check the closure of the (ferromagnetic) little group ===
 ABI_ALLOCATE(symafm_ltg,(Ltg%nsym_Ltg))
 symafm_ltg(:)=1
 call chkgrp(Ltg%nsym_Ltg,symafm_ltg,symrec_Ltg,ierr)
 ABI_CHECK(ierr==0,"Error in group closure")

 ABI_DEALLOCATE(symafm_ltg)
 !
 ! === Find the irreducible zone associated to ext_pt ===
 ! * Do not use time-reversal since it has been manually introduced previously
 ABI_ALLOCATE(indkpt1,(nbz))
 ABI_ALLOCATE(wtk_folded,(nbz))
 ABI_ALLOCATE(wtk,(nbz))
 wtk=one; iout=0; dummy_timrev=0

 call symkpt(0,Cryst%gmet,indkpt1,iout,kbz,nbz,nkibzq,Ltg%nsym_Ltg,symrec_Ltg,dummy_timrev,wtk,wtk_folded)

 ABI_DEALLOCATE(indkpt1)
 ABI_DEALLOCATE(wtk)

 Ltg%nibz_Ltg=nkibzq
 !
 ! === Set up table in the BZ === 
 ! * 0 if the point does not belong to IBZ_xpt, 1 otherwise
 ABI_ALLOCATE(Ltg%ibz2bz,(nkibzq))
 Ltg%ibzq(:)=0; Ltg%bz2ibz(:)=0; Ltg%ibz2bz(:)=0

 ind=0; enough=0
 do ik=1,nbz
   if (wtk_folded(ik)>tol8) then
     ind=ind+1
     Ltg%ibzq(ik)=1
     Ltg%bz2ibz(ik) =ind
     Ltg%ibz2bz(ind)=ik
   end if
 end do
 ABI_CHECK(ind==Ltg%nibz_Ltg," BUG ind/=Ltg%nibz_Ltg")
 !
 ! === Reconstruct full BZ starting from IBZ_q  ===
 ! Calculate appropriate weight for each item (point,symmetry operation,time-reversal)
 Ltg%tab=0; Ltg%tabo=0; Ltg%tabi=0; Ltg%wtksym(:,:,:)=0

 ! === Zero no. of k-points found ===
 ntest=0
 ABI_ALLOCATE(ktest,(3,nbz))
 ktest=zero  

 do ik=1,nbz
   if (Ltg%ibzq(ik)/=1) CYCLE
   ! * Loop over symmetry operations S and time-reversal.
   ! * Use spatial inversion instead of time reversal whenever possible.
   do itim=1,timrev
     do isym=1,nsym

      ! Form IS k only for (IS) pairs in the (ferromagnetic) little group.
      if (symafm(isym)==-1) CYCLE
      if (Ltg%preserve(itim,isym)==0) CYCLE
      knew(:)=(3-2*itim)*MATMUL(symrec(:,:,isym),kbz(:,ik))
      !
      ! === Check whether it has already been found (to within a RL vector) ===
      iold=0
      do itest=1,ntest
        if (isamek(knew(:),ktest(:,itest),gg)) iold=iold+1
      end do

      if (iold==0) then 
        ! == Found new BZ point ===
        ! For this point the operation (isym,itim) must be considered to reconstruct the full BZ
        Ltg%wtksym(itim,isym,ik)=1
        ntest=ntest+1
        ktest(:,ntest)=knew(:)
        ! 
        ! === Now find knew in the BZ array ===
        found=.FALSE.
        do idx=1,nbz
          if (isamek(knew(:),kbz(:,idx),gg)) then ! They are the same within a RL vector
            Ltg%tab (idx)=ik
            Ltg%tabo(idx)=isym
            Ltg%tabi(idx)=3-2*itim
            found=.TRUE.; EXIT
          end if 
        end do
        if (.not.found) then 
          write(msg,'(a,3f12.6,a)')'Not able to find the ',knew(:),' in the array BZ '
          MSG_ERROR(msg)
        end if
      end if

     end do !isym
   end do !itim 

 end do !nbz

 ABI_DEALLOCATE(ktest)

 if (ntest/=nbz) then 
   write(msg,'(a,i5)')' ntest-nbz = ',ntest-nbz
   MSG_BUG(msg)
 end if 

 if (SUM(Ltg%wtksym)/=nbz) then 
   write(msg,'(a,i5)')' sum(Ltg%wtksym)-nbz = ',SUM(Ltg%wtksym)-nbz
   MSG_BUG(msg)
 end if 

 nullify(Ltg%igmG0); Ltg%max_kin_gmG0=zero

 if (npwe>0.and.PRESENT(gvec)) then 
   npwvec = SIZE(gvec,DIM=2)
   ! This correspond to the case in which we need to know the index of G-Go in the gvec array
   ! where G is one of the npwe vectors. This is required in screening but not in sigma.
   ! The drawback is that the effective G sphere used to calculate the oscillators must be smaller 
   ! that gvec if we want to avoid possible aliasing effects. Lifting this constraint would require 
   ! a lot of boring coding. (no need to do this if ext_pt=zero, but oh well)
   ABI_ALLOCATE(Ltg%igmG0,(npwe,timrev,nsym))
   Ltg%igmG0(:,:,:)=0
   max_kin=zero
   !
   ! === Loop over symmetry operations S and time-reversal ===
   do itim=1,timrev
     do isym=1,nsym
       ! * Form IS k only for (IS) pairs in the little group
       if (symafm(isym)==-1) CYCLE
       if (Ltg%preserve(itim,isym)/=0) then 
         g0(:)=Ltg%g0(:,itim,isym)
         do ige=1,npwe
           gmG0(:)=gvec(:,ige)-g0(:)
           kin=half*normv(gmG0,Cryst%gmet,'G')**2 
           max_kin=MAX(max_kin,kin)

           found=.FALSE.
           do igpw=1,npwvec
             if (ALL(gvec(:,igpw)-gmG0(:)==0)) then 
               Ltg%igmG0(ige,itim,isym)=igpw
               found=.TRUE.; EXIT
             end if
           end do
           if (.not.found) then 
             write(msg,'(5a,f8.3,2a,3i5)')&
&              ' Not able to found G-G0 in the largest G-spere ',ch10,&
&              ' Decrease the size of epsilon or, if possible, increase ecutwfn (>ecuteps) ',ch10,&
&              ' Minimum required cutoff energy for G-G0 sphere= ',kin,ch10,&
&              ' G0 = ',g0(:)
             MSG_ERROR(msg)
           end if 
         end do 
       end if
     end do 
   end do 
   Ltg%max_kin_gmG0=max_kin
 end if
 ABI_DEALLOCATE(symrec_Ltg)

#ifdef DEBUG_MODE
 do ik=1,nbz
   if (ABS(SUM(Ltg%wtksym(1,:,ik)+Ltg%wtksym(2,:,ik))-wtk_folded(ik))>tol6) then 
     write(std_out,*)' sum(Ltg%wtksym,ik)-wtk_folded(ik) = ',sum(Ltg%wtksym(1,:,ik)+Ltg%wtksym(2,:,ik))-wtk_folded(ik)
     write(std_out,*)Ltg%wtksym(1,:,ik),Ltg%wtksym(2,:,ik),wtk_folded(ik)
     write(std_out,*)ik,kbz(:,ik) 
     MSG_BUG("Wrong weight")
   end if 
 end do
 do ik=1,nbz 
   knew = Ltg%tabi(ik) * MATMUL(symrec(:,:,Ltg%tabo(ik)),kbz(:,Ltg%tab(ik)))
   if (.not.isamek(knew,kbz(:,ik),gg)) then 
     write(std_out,*)knew,kbz(:,ik) 
     write(std_out,*)Ltg%tabo(ik),Ltg%tabi(ik),Ltg%tab(ik)
     MSG_BUG("Wrong tables")
   end if 
 end do
#endif

 ABI_DEALLOCATE(wtk_folded)

 DBG_EXIT("COLL")

end subroutine setup_little_group
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/nullify_little_group_0D
!! NAME
!! nullify_little_group_0D
!!
!! FUNCTION
!!  Set all pointer in the little_group datatype to NULL.
!!
!! INPUTS
!! Ltg= The datatype whose pointers have to be nullified.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_bz_mesh
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine nullify_little_group_0D(Ltg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_little_group_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Little_group),intent(inout) :: Ltg

! *********************************************************************
 
 !@little_group
 nullify(Ltg%g0)
 nullify(Ltg%ibzq)
 nullify(Ltg%bz2ibz)
 nullify(Ltg%ibz2bz)
 nullify(Ltg%igmG0)
 nullify(Ltg%flag_umklp)
 nullify(Ltg%preserve)          
 nullify(Ltg%tab)
 nullify(Ltg%tabo)
 nullify(Ltg%tabi)
 nullify(Ltg%wtksym)

end subroutine nullify_little_group_0D
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/nullify_little_group_1D
!! NAME
!! nullify_little_group_1D
!!
!! FUNCTION
!!  Set all pointer in the little_group datatype to NULL.
!!
!! INPUTS
!! Ltg= The datatype whose pointers have to be nullified.
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine nullify_little_group_1D(Ltg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_little_group_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Little_group),intent(inout) :: Ltg(:)

!Local variables-------------------------------
!scalars
 integer :: ii

! *********************************************************************
 
 do ii=1,SIZE(Ltg)
   call nullify_little_group_0D(Ltg(ii))
 end do

end subroutine nullify_little_group_1D
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/destroy_little_group_0D
!! NAME
!! destroy_little_group_0D
!!
!! FUNCTION
!!  Deallocate all associated pointers defined in the Little_group data type.
!!
!! INPUTS
!!   Ltg=datatype to be freed
!!
!! OUTPUT
!!
!! PARENTS
!!      m_bz_mesh
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_little_group_0D(Ltg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_little_group_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Little_group),intent(inout) :: Ltg

! *********************************************************************

 DBG_ENTER("COLL")

 !@little_group
 if (associated(Ltg%g0      ))     then
   ABI_DEALLOCATE(Ltg%g0)
 end if
 if (associated(Ltg%ibzq    ))     then
   ABI_DEALLOCATE(Ltg%ibzq)
 end if
 if (associated(Ltg%bz2ibz  ))     then
   ABI_DEALLOCATE(Ltg%bz2ibz)
 end if
 if (associated(Ltg%ibz2bz  ))     then
   ABI_DEALLOCATE(Ltg%ibz2bz)
 end if
 if (associated(Ltg%igmG0    ))    then
   ABI_DEALLOCATE(Ltg%igmG0)
 end if
 if (associated(Ltg%flag_umklp))   then
   ABI_DEALLOCATE(Ltg%flag_umklp)
 end if
 if (associated(Ltg%preserve))     then
   ABI_DEALLOCATE(Ltg%preserve)
 end if
 if (associated(Ltg%tab     ))     then
   ABI_DEALLOCATE(Ltg%tab)
 end if
 if (associated(Ltg%tabo    ))     then
   ABI_DEALLOCATE(Ltg%tabo)
 end if
 if (associated(Ltg%tabi    ))     then
   ABI_DEALLOCATE(Ltg%tabi)
 end if
 if (associated(Ltg%wtksym  ))     then
   ABI_DEALLOCATE(Ltg%wtksym)
 end if

 DBG_EXIT("COLL")

end subroutine destroy_little_group_0D
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/destroy_little_group_1D
!! NAME
!! destroy_little_group_1D
!!
!! FUNCTION
!!  Deallocate all the associated pointers defined in the Little_group data type.
!!
!! INPUTS
!!   Ltg=datatype to be freed
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_little_group_1D(Ltg)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_little_group_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Little_group),intent(inout) :: Ltg(:)

!Local variables-------------------------------
 integer :: ipt

! *********************************************************************

 do ipt=1,SIZE(Ltg)
   call destroy_little_group_0D(Ltg(ipt))
 end do

end subroutine destroy_little_group_1D
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/print_little_group
!! NAME
!! print_little_group
!!
!! FUNCTION
!!  Print info on the little_group data type.
!!
!! INPUTS
!!  Ltg=the datatype to be printed
!!  [unit]=the unit number for output 
!!  [prtvol]=verbosity level
!!  [mode_paral]=either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,cchi0q0_intraband
!!      check_completeness,cohsex_me
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_little_group(Ltg,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_little_group'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 type(Little_group),intent(in) :: Ltg

!Local variables-------------------------------
!scalars
 integer :: itim,my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg
!arrays
 integer :: nop(Ltg%timrev),nopg0(Ltg%timrev)

! *********************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 write(msg,'(4a,3es16.8,2a,i5,a,i5,2a,i3,a,i3)')ch10,&
&  ' ==== Little Group Info ==== ',ch10,&
&  '  External point ',Ltg%ext_pt,ch10,&
&  '  Number of points in the IBZ defined by little group  ',Ltg%nibz_Ltg,'/',Ltg%nbz,ch10,&
&  '  Number of operations in the little group : ',Ltg%nsym_Ltg,'/',Ltg%nsym_sg
 call wrtout(my_unt,msg,my_mode)

 nop=0 ; nopg0=0
 do itim=1,Ltg%timrev
   nop  (itim)=SUM(Ltg%preserve  (itim,:))
   nopg0(itim)=SUM(Ltg%flag_umklp(itim,:))
 end do

 do itim=1,Ltg%timrev
   if (itim==1) then 
     write(msg,'(a,2(a,i2,a))')ch10,&
&      '  No time-reversal symmetry with zero umklapp: ',nop(1)-nopg0(1),ch10,&
&      '  No time-reversal symmetry with non-zero umklapp: ',nopg0(1),ch10
     call wrtout(my_unt,msg,my_mode)
   else if (itim==2) then 
     write(msg,'(a,2(a,i2,a))')ch10,&
&      '  time-reversal symmetry with zero umklapp: ',nop(2)-nopg0(2),ch10,&
&      '  time-reversal symmetry with non-zero umklapp: ',nopg0(2),ch10
     call wrtout(my_unt,msg,my_mode)
   end if
 end do

end subroutine print_little_group 
!!***

!----------------------------------------------------------------------

!!****f* m_bz_mesh/box_len
!! NAME
!!  box_len
!!
!! FUNCTION
!!   Given a direction in q-space defined by the q-point qpt, this function returns
!!   the length of the vector connecting the origin with one the faces of the cell 
!!   defined by the lattice vectors gprimd.
!!
!! INPUTS
!!   qpt(3)=The reduced coordinates of the q-point defining the direction. Normalization is not mandatory.
!!   gprimd(3,3)=Cartesian coordinates of the vectors defining the lattice.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function box_len(qpt,gprimd)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'box_len'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: box_len
!arrays
 real(dp),intent(in) :: qpt(3),gprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: idir,iplane
 real(dp) :: x1,x2,x3
!arrays
 real(dp) :: my_qpt(3),gmet(3,3),q0box(3)

! *************************************************************************

!Compute reciprocal space metric
 gmet = MATMUL(TRANSPOSE(gprimd),gprimd)

!Rotate the input q-point such that it is always in the first octant then normalize it.
!Bravais lattices are invariant under inversion of any of the basis vectors.
 my_qpt = ABS(qpt)/normv(qpt,gmet,"G")

!Check whether the q is along one of the reduced directions.
 idir=0; if (COUNT(ABS(qpt)<tol16) == 2) idir = imax_loc(ABS(qpt))

 if (idir/=0) then ! easy as q is along vector idir.
   box_len =  normv(gprimd(:,idir),gmet,"G")
   RETURN

 else
   !iplane/=0 means that q is placed on one the planes defined by two reciprocal lattice vectors.
   iplane=0; if (COUNT(ABS(qpt)<tol16) == 1) iplane = imin_loc(ABS(qpt))
   q0box = (/-1,-1,-1/)

   if (iplane/=1) then
     x1=one
     x2=my_qpt(2)/my_qpt(1)
     x3=my_qpt(3)/my_qpt(1)
     if (x2<=one+tol16 .and. x3<=one+tol16) q0box=(/x1,x2,x3/)
   end if
   if (iplane/=2) then
     x1=my_qpt(1)/my_qpt(2)
     x2=one
     x3=my_qpt(3)/my_qpt(2)
     if (x1<=one+tol16 .and. x3<=one+tol16) q0box=(/x1,x2,x3/)
   end if
   if (iplane/=3) then
     x1=my_qpt(1)/my_qpt(3)
     x2=my_qpt(2)/my_qpt(3)
     x3=one
     if (x1<=one+tol16 .and. x2<=one+tol16) q0box=(/x1,x2,x3/)
   end if

   if (ALL(q0box==(/-1,-1,-1/) )) then 
     MSG_BUG("Cannot found q0box")
   end if

   box_len = normv(q0box,gmet,"G")
   RETURN
 end if

end function box_len
!!***

!----------------------------------------------------------------------

END MODULE m_bz_mesh
!!***
