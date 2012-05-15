!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_bands_sym
!! NAME
!! m_bands_sym
!!
!! FUNCTION
!! This module defines structures and provides procedures used to find 
!! the irreducible representations associated to electronic eigenstates.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_bands_sym

 use m_profiling
 
 use defs_basis
 use m_errors 

 use defs_datatypes,   only : coeffi1_type
 use m_io_tools,       only : file_exist, flush_unit
 use m_numeric_tools,  only : print_arr, set2unit, get_trace
 use m_abilasi,        only : xgeev, xginv
 use m_crystal,        only : crystal_structure, idx_spatial_inversion
 use m_defs_ptgroups,  only : point_group_t, irrep_t
 use m_ptgroups,       only : get_classes, init_point_group, destroy_irrep,&
&                             copy_irrep, init_irrep, mult_table, sum_irreps !, polish_irreps

 implicit none

 private 
!!***

 integer,private,parameter :: BSYM_NOERROR             = 0
 integer,private,parameter :: BSYM_ACCDEG_ERROR        = 10 
 integer,private,parameter :: BSYM_CLASSIFICATION_ERROR= 11
 integer,private,parameter :: BSYM_ORTHO_ERROR         = 12
 integer,private,parameter :: BSYM_UNITARY_ERROR       = 13
 integer,private,parameter :: BSYM_PTG_WRONG_MAPPING   = 20
 integer,private,parameter :: BSYM_HERRING_WRONG_TEST  = 30
 integer,private,parameter :: BSYM_HEUR_WRONG_NCLASSES = 40
 integer,private,parameter :: BSYM_HEUR_WRONG_DIMS     = 41

!----------------------------------------------------------------------

!!****t* m_bands_sym/bands_symmetries
!! NAME
!!  bands_symmetries
!!
!! FUNCTION
!!  Dataype gathering data and tables needed to analize the symmetries 
!!  of electronic states at a given k-point via Group Theory.
!!
!! SOURCE

 type,public :: bands_symmetries

  integer :: nspinor                      
  ! Number of spinorial components.

  integer :: first_ib
  ! Index of the first treated band.

  integer :: nbnds
  ! Number of bands for this k-point and spin.

  integer :: nclass
  ! The number of classes in the group of k.

  integer :: nsym_gk
  ! Number of symmetries in the group of k. Namely that the set of symmetries such that Sk = k +G0.

  integer :: nsym_trgk
  ! Number of symmetries in the extended group of k. Namely that the set of symmetries such that -Sk = k + G0.

  integer :: err_status = BSYM_NOERROR
  ! Flag signaling if the classification algorithm succeed or not.

  real(dp) :: tol_deg
  ! Energy tolerance below which two states are considered degenerate.

  logical :: can_use_tr
  ! .TRUE. if time-reversal can be used

  logical :: only_trace                   
  ! if .TRUE. only the trace of a single matrix per class is calculated
  ! this is the standard way used to analyze bands symmetries. If .FALSE.
  ! the full matrices of the irreducible representations are calculated and stored

  logical :: has_spatial_inv
  ! .TRUE. if the inversion belongs to the space group

  logical :: nonsymmorphic_at_zoneborder
  ! if .TRUE. analysis cannot be performed since kpt is
  ! at border zone and non-zero fractional translations are present in the space group

  logical :: has_chtabs
  ! True if Ref_irreps and character tables are available (tables are initialized either 
  ! from point group irreps or from an external database downloaded from the Bilbao server)

  real(dp) :: kpt(3)
  ! The crystalline momentum of the wavefunctions in reduced coordinates.

  character(len=500) :: err_msg="None"

  integer,pointer :: g0(:,:)  SET2NULL
  ! g0(3,nsym_gk)
  ! The umklapp g0 vector associated to each little group operation.

  integer,pointer :: tr_g0(:,:)  SET2NULL
  ! tr_g0(3,nsym_trgk)
  ! The umklapp g0 vector associated to each little group operation.

  integer :: ndegs
  ! Number of degenerate states.

  integer,pointer :: nelements(:)  SET2NULL
  ! nelements(nclass)
  ! Number of symmetry operations in each class.

  integer,pointer :: sgk2symrec(:)  SET2NULL
  ! sgk2symrec(nsym_gk)
  ! Mapping between the symmetries of the group of k and the symrec(l) array. 
  ! The symmetries of the little group are always packed in classes to facilitate
  ! the calculation of the character of the irrep. Abinit symmetries are randomly ordered. 

  integer,pointer :: tr_sgk2symrec(:)  SET2NULL
  ! trsgk2symrec(nsym_trgk)
  ! Mapping between the symmetries of the group of k and the symrec(l) array. 
  ! The symmetries of the little group are always packed in classes to facilitate
  ! the calculation of the character of the irrep. Abinit symmetries are randomly ordered. 

  integer,pointer :: herring_test(:) SET2NULL
  ! herring_test(nclass)
  ! The result of Herring test for each irreducible representantion of the group of k.
  ! Possible values are:
  ! +1 
  ! 0
  ! -1

  integer,pointer :: b2irrep(:)  SET2NULL
  ! b2irrep(nbnds)
  ! For each band, it gives the index of the irreducible representation in Ref_irreps.

  type(coeffi1_type),pointer :: irrep2b(:)  SET2NULL
  ! irrep2b(0:nclass)%value(:)
  ! Ragged arrays with the mapping between the set of irreducible representation and the band indices. 
  ! irrep2b(irp)%value(:) gives the indeces of the states belonging to irrep irp, irp=1,nclass
  ! irrep2b(0)%value(:) stores the indeces of the states that have not been classified due to 
  !   the presence of an accidental degeneracy.

  integer,pointer :: degs_bounds(:,:)  SET2NULL
  ! degs_bounds(2,ndegs)
  !   degs_bounds(1,idg)= first band index of the degenerate set idg=1,ndegs
  !   degs_bounds(2,idg)= final band index of the degenerate set idg=1,ndegs

  integer,pointer :: degs_dim(:)  SET2NULL
  ! degs_dim(ndegs) 
  ! Number of states in each degenerate subspace. Cannot be larger that nclass provided
  ! that no accidental degeneracy occurs.

  !% integer,pointer :: class_ids(:,:)   SET2NULL 
  ! class_ids(2,nclass)
  ! (1,icl) = index of the first symmetry of class icl
  ! (2,icl) = index of the last symmetry of class icl
  ! Note that symmetries in sym are packed in classes.

  type(irrep_t),pointer :: Calc_irreps(:)   SET2NULL
  ! Calc_irreps(ndegs)
  !  The representations of the little group of k calculated from the wavefunctions. <\phi_nk|R_t|\phi_mk>
  !  where R_t belong to the little group of k.
  !  They represent an unitary irreducible representation provided that no accidental degeneracy occurs.

  type(irrep_t),pointer :: trCalc_irreps(:)   SET2NULL
  ! trCalc_irreps(ndegs)
  !  The representations of the little group of k calculated from the wavefunctions. <\phi_nk|R_t|\phi_mk>
  !  where R_t belong to the little group of k.
  !  They represent an unitary irreducible representation provided that no accidental degeneracy occurs.

  type(irrep_t),pointer :: Ref_irreps(:)   SET2NULL
  ! Irreps(nclass)
  !   Reference irreducible representations of the group of k derived from the point group
  !   or from the external database downloaded from the Bilbao web site.

 end type bands_symmetries
!!***

!----------------------------------------------------------------------

 public :: init_bands_symmetries
 public :: nullify_bands_symmetries
 public :: print_bands_symmetries
 public :: destroy_bands_symmetries
 public :: finalize_bands_sym
 public :: symmetrize_me
 public :: bsym_failed

 !public :: polish_irreps ! TODO method of Irreps_t, therefore should be moved to m_ptgroups. 
                          ! but first one has to solve the dependency on m_abilasi and scalapack

 interface destroy_bands_symmetries
   module procedure destroy_bands_symmetries_0D
   module procedure destroy_bands_symmetries_2D
 end interface destroy_bands_symmetries

 interface nullify_bands_symmetries
   module procedure nullify_bands_symmetries_0D
   module procedure nullify_bands_symmetries_2D
 end interface nullify_bands_symmetries

CONTAINS  !==============================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_bands_sym/init_bands_symmetries
!! NAME
!! init_bands_symmetries
!!
!! FUNCTION
!!  Initialize a bands_symmetries datatype containing data and parameters 
!!  needed to analyze the irreducible representations at a particular k-point....
!!
!! INPUTS
!!  kpt_in(3)=The k-point where the classification of bands is required.
!!  Cryst<Crystal_structure>=Datatype describing the unit cell and its symmetries.
!!  nspinor=number of spinorial components
!!  nsppol=number of independent polarizations
!!  first_ib=Index of the first band.
!!  nbnds=Number of bands for this k-point.
!!  ene_k(nbnds)=energies for this k-point. ene_k(1) corresponds to band first_ib.
!!  EDIFF_TOL=tolerance below which two states are considered to belong to the same irreducible representation 
!!
!! OUTPUT
!!  Bsym<bands_symmetries>= Initialized data type gathering information of the small group
!!     of the k-point as well as the irreducible representations.
!!
!! NOTES
!!   The present implementation does NOT work at zone border if the little group of
!!   kpt_in is non-symmorphic namely thers is at lest a symmetry operation with non-zero tnons.
!!
!! PARENTS
!!      classify_bands
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine init_bands_symmetries(Bsym,kpt_in,Cryst,only_trace,nspinor,first_ib,nbnds,EDIFF_TOL,ene_k,tolsym)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_bands_symmetries'
 use interfaces_14_hidewrite
 use interfaces_27_toolbox_oop
 use interfaces_32_util
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,nspinor,first_ib
 real(dp),intent(in) :: EDIFF_TOL,tolsym
 logical,intent(in) :: only_trace
 type(crystal_structure),intent(in) :: Cryst
 type(bands_symmetries),intent(out) :: BSym
!arrays
 real(dp),intent(in) :: ene_k(nbnds),kpt_in(3)

!Local variables-------------------------------
!scalars
 integer :: dim_degs,iband,idg,irp,nacc_deg,isym_gk,grp_ierr
 integer :: nsym_afm,idx_fm,idx_gk,idx_trgk,isym,jsym,dummy_timrev !,iholohedry
 integer :: iel,icls,msym,iord !isym1,!iprod,dim_irrep,icls2, isym2,isym_tr,
 integer :: spgroup,chkprim !,ptgroupma
 real(dp) :: mkt
 !complex(dpc) :: phase_k 
 character(len=5) :: ptgroup,ptgroup_name
 character(len=10) :: spgroup_str
 character(len=1000) :: msg
 character(len=fnlen) :: lgroup_fname
!arrays
 integer :: inversion(3,3) 
 integer,allocatable :: degs_bounds(:,:),dim_irreps(:)
 integer :: bravais(11),sym_axis(3)
 real(dp) :: pmat1(3,3),pmat2(3,3),pmat3(3,3),pmat4(3,3),pmat5(3,3),pmat6(3,3)
 !real(dp) :: genafm(3)
 !integer :: rot2(3,3)
 !integer,allocatable :: mtab(:,:)
 integer,allocatable :: elements_idx(:,:),tmp_nelements(:)
 integer,allocatable :: found(:),symrec_fm(:,:,:),fm2symrec(:)
 integer,allocatable :: ksym_table(:,:,:),sgk(:,:,:),tr_sgk(:,:,:),dum_symafm(:) 
 integer,allocatable :: new_idx(:),new_g0(:,:),tmp_symrec(:,:,:),conv_symrec(:,:,:) !,tr_conv_symrec(:,:,:)
 real(dp) :: conv_gprimd(3,3),axes(3,3) !,tau2(3) 
 !complex(dpc),allocatable :: her_test(:) !,mat_test(:,:)
 complex(dpc),allocatable :: phase_mkt(:)
 type(point_group_t) :: Ptg

! *************************************************************************

 DBG_ENTER("COLL")

 !@bands_symmetries
 call nullify_Bands_Symmetries(BSym)

 BSym%err_status= BSYM_NOERROR
 inversion=RESHAPE((/-1,0,0,0,-1,0,0,0,-1/),(/3,3/)) 
 !
 ! ====================================
 ! ==== Initialize basic variables ====
 ! ====================================
 BSym%nspinor        = nspinor
 BSym%first_ib       = first_ib
 BSym%nbnds          = nbnds
 BSym%only_trace     = only_trace
 BSym%tol_deg        = EDIFF_TOL
 BSym%has_spatial_inv= (idx_spatial_inversion(Cryst) /= 0)
 BSym%can_use_tr     = .TRUE. !TODO this should be input 
 BSym%has_chtabs     = .FALSE.
 BSym%kpt            = kpt_in(:)
 Bsym%nonsymmorphic_at_zoneborder=.FALSE.
 !
 ! ===============================
 ! === Locate degenerate_bands ===
 ! ===============================
 BSym%ndegs=1 

 ABI_ALLOCATE(degs_bounds,(2,nbnds))
 degs_bounds=0; degs_bounds(1,1)=1

 do iband=2,nbnds
   if (ABS(ene_k(iband)-ene_k(iband-1))>EDIFF_TOL) then
     degs_bounds(2,BSym%ndegs) = iband-1 + (first_ib-1)
     BSym%ndegs=BSym%ndegs+1
     degs_bounds(1,BSym%ndegs) = iband + (first_ib-1)
   end if
 end do
 degs_bounds(2,BSym%ndegs)=nbnds + (first_ib-1)

 ABI_ALLOCATE(BSym%degs_bounds,(2,BSym%ndegs))
 BSym%degs_bounds = degs_bounds(:,1:BSym%ndegs)
 ABI_DEALLOCATE(degs_bounds)
 !
 ! Each band is initialized as "Unknown".
 ABI_ALLOCATE(BSym%b2irrep,(BSym%nbnds))
 BSym%b2irrep = 0
 !
 ! ==================================
 ! ==== Find the group of kpt_in ====
 ! ==================================
 ! * The small point group is the subset of symrec such that $ S q = q + g0 $
 ! * Symmetries are packed in classes.
 ! * For the time being, AFM symmetries are not treated.

 write(msg,'(a,3(1x,f7.4))')" Finding the little group of k-point: ",Bsym%kpt
 call wrtout(std_out,msg,"COLL")
 !
 ! * Only FM symmetries are used.
 nsym_afm = COUNT(Cryst%symafm==1)

 if (nsym_afm/=Cryst%nsym) then
   write(msg,'(4a)')ch10,&
&    " Band classification in terms of magnetic space groups not coded! ",ch10,&
&    " Only the ferromagnetic subgroup will be used "
   MSG_COMMENT(msg)
 end if

 ABI_ALLOCATE(symrec_fm,(3,3,nsym_afm))
 ABI_ALLOCATE(fm2symrec,(nsym_afm))

 idx_fm = 0
 do isym=1,Cryst%nsym
   if (Cryst%symafm(isym) == 1) then 
     idx_fm = idx_fm+1
     symrec_fm(:,:,idx_fm) = Cryst%symrec(:,:,isym)
     fm2symrec(idx_fm) = isym
   end if
 end do
 !
 ! * Find symmetries that preserve k.   
 ABI_ALLOCATE(ksym_table,(4,2,nsym_afm))

 call symq3(nsym_afm,Bsym%kpt,ksym_table,symrec_fm,dummy_timrev,prtvol=0)
                                                                                                 
 Bsym%nsym_gk  =COUNT(ksym_table(4,1,:)==1)  ! # S such that  S k = k +G0

 Bsym%nsym_trgk=0
 if (Bsym%can_use_tr) Bsym%nsym_trgk=COUNT(ksym_table(4,2,:)==1)  ! # S such that -S k = k +G0

 ! Allocate workspace.
 ABI_ALLOCATE(sgk,(3,3,Bsym%nsym_gk))
 ABI_ALLOCATE(tr_sgk,(3,3,Bsym%nsym_trgk))

 ! Allocate mapping little-group --> symrec and table for umklapps.
 ABI_ALLOCATE(Bsym%sgk2symrec,(Bsym%nsym_gk))
 ABI_ALLOCATE(Bsym%g0,(3,Bsym%nsym_gk))
 ABI_ALLOCATE(Bsym%tr_sgk2symrec,(Bsym%nsym_trgk))
 ABI_ALLOCATE(Bsym%tr_g0,(3,Bsym%nsym_trgk))
 !
 ! Important NOTE: 
 ! If nonsymmorphic_at_zoneborder symmetry analysis cannot be performed unless
 ! an external database retrieved from the bilbao server (REPRES) is found.
 !
 idx_gk=0; idx_trgk=0
 Bsym%sgk2symrec=-999; Bsym%tr_sgk2symrec=-999
 do isym=1,nsym_afm

   if (ksym_table(4,1,isym)==1) then ! S k = k +G0
     idx_gk=idx_gk+1
     sgk(:,:,idx_gk)=symrec_fm(:,:,isym)
     Bsym%g0(:,idx_gk)=ksym_table(1:3,1,isym)
     Bsym%sgk2symrec(idx_gk)=fm2symrec(isym)
     if (ANY(ksym_table(1:3,1,isym)/=0).and.(ANY(ABS(Cryst%tnons(:,fm2symrec(isym)))>tol6))) then 
        Bsym%nonsymmorphic_at_zoneborder=.TRUE.
     end if
   end if

   if (Bsym%can_use_tr.and.ksym_table(4,2,isym)==1) then ! -S k = k +G0
     idx_trgk=idx_trgk+1
     tr_sgk(:,:,idx_trgk)=symrec_fm(:,:,isym)
     Bsym%tr_g0(:,idx_trgk)=ksym_table(1:3,2,isym)
     Bsym%tr_sgk2symrec(idx_trgk)=fm2symrec(isym)
   end if
 end do

 ABI_DEALLOCATE(ksym_table)
 ABI_DEALLOCATE(symrec_fm)
 ABI_DEALLOCATE(fm2symrec)

! ==========================================
! ==== Divide the operations in classes ====
! ==========================================
 ABI_ALLOCATE(dum_symafm,(Bsym%nsym_gk))
 dum_symafm=1

 call chkgrp(Bsym%nsym_gk,dum_symafm,sgk,grp_ierr)

 ABI_CHECK(grp_ierr==0,"chkgrp failed")
 ABI_DEALLOCATE(dum_symafm)

 ABI_ALLOCATE(tmp_nelements,(Bsym%nsym_gk))
 ABI_ALLOCATE(elements_idx,(Bsym%nsym_gk,Bsym%nsym_gk))

 call get_classes(Bsym%nsym_gk,sgk,Bsym%nclass,tmp_nelements,elements_idx)

 ABI_ALLOCATE(Bsym%nelements,(Bsym%nclass))
 Bsym%nelements = tmp_nelements(1:Bsym%nclass)
 ABI_DEALLOCATE(tmp_nelements)
 !
 ! From the list of symmetry operations and the lattice vectors, determine the 
 ! Bravais information including the holohedry, the centering, the coordinate of 
 ! the primitive vectors in the conventional vectors, as well as the point group,
 !
 !FIXME: this should be generalized - for supercells and so on nsym can be arbitrarily large.
 msym=192
 ABI_ALLOCATE(tmp_symrec,(3,3,msym))
 tmp_symrec(:,:,1:Bsym%nsym_gk)=sgk
                                                                                                     
 call symbrav(bravais,msym,Bsym%nsym_gk,ptgroup,Cryst%gprimd,tmp_symrec,tolsym,axis=sym_axis)
                                                                                                     
 ABI_DEALLOCATE(tmp_symrec)

 write(std_out,'(a)')" symptgroup returned point group: "//TRIM(ptgroup)
 write(std_out,'(a,i2)')" iholohedry ",bravais(1)
 write(std_out,'(a,i2)')" center     ",bravais(2)
 write(std_out,'(a,9i3)')" gprimd in the axes of the conventional bravais lattice (*2 if center/=0)",bravais(3:11)
 write(std_out,'(a,3i3)')" sym_axis ",sym_axis

 ! Branching:
 ! 1) If the little group is not symmorphic_at_zoneborder we can 
 !    classify the states using the irreducible representation of the point group.
 !
 ! 2) If the little group is symmorphic_at_zoneborder, we have to rely on 
 !    an external database retrieved from the Bilbao server in order to classify the states.
 !    If the file is not available, we only know the number of classes but neither their 
 !    character nor the dimension of the irreducible representation.
 ! 
 if (Bsym%nonsymmorphic_at_zoneborder) then 

   spgroup=0
   chkprim=1 ! Cell must be primitive.
   !call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
   !call symspgr(bravais,Cryst%nsym,spgroup,Cryst%symrel,Cryst%tnons,tolsym)

   !call symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)

   call int2char(spgroup,spgroup_str)
   lgroup_fname = "lgroup_"//TRIM(spgroup_str)

   if (file_exist(lgroup_fname)) then
     MSG_ERROR("Not coded")

     ! Read little groups from the external database.
     !% call init_groupk_from_file(Lgrp,spgroup,lgroup_fname,ierr)

     ! Save the irreducible representations in BSym.
     ! Reorder symmetries such that they correspond to the Bilbao database.
     !% allocate(BSym%Ref_irreps(BSym%nclass)) 
     !% call copy_irrep(Irreps, BSym%Ref_irreps)

   else 
     write(msg,'(7a)')&
&      " Non-symmorphic small group and zone border. ",ch10,&
&      " External file: ",TRIM(lgroup_fname)," containing Bilbao tables not found ",ch10,&
&      " Character analysis cannot be performed. Accidental degeneracies cannot be detected. "
     MSG_WARNING(msg)

     Bsym%has_chtabs = .FALSE.

     ! Reorder indeces such that symmetries are packed in classes.
     ABI_ALLOCATE(new_idx,(Bsym%nsym_gk))
     ABI_ALLOCATE(new_g0,(3,Bsym%nsym_gk))
     new_g0=0; iord = 0
     do icls=1,Bsym%nclass
       do iel=1,Bsym%nelements(icls)
         iord = iord+1
         jsym = elements_idx(iel,icls)
         new_idx(iord)  = Bsym%sgk2symrec(jsym)
         new_g0(:,iord) = Bsym%g0(:,jsym)
       end do
     end do
                                                                                               
     Bsym%sgk2symrec = new_idx
     Bsym%g0 = new_g0
                                                                                               
     ABI_DEALLOCATE(new_idx)
     ABI_DEALLOCATE(new_g0)
   end if ! file exists 

 else
   !
   ! **** This part is still under development. It might not work for particular **** 
   ! **** orientations of the unit cell or particular lattices.                  ****
   !
   ! The symmetries in the Bilbao database refer to the conventional unit cells.
   ! Therefore we have to map the abinit symmetries (in reduced coordinates)
   ! onto the Bilbao dataset. Bilbao standard settings are:
   !
   ! * unique axis b (cell choice 1) for space groups withing the monoclinic system
   ! * obverse triple hexagonal unit cell R space groups.
   ! * origin choice two - inversion center at (0, 0, 0) - for the centrosymmetric
   !   space groups for which there are two origins choices, within the 
   !   orthorombic, tetragonal and cubic system.

   ! 1) Retrieve the rotation matrices and the irreducible representations (Bilbao setting).
   call init_point_group(Ptg,ptgroup)

   BSym%has_chtabs = .TRUE.

   ABI_CHECK(BSym%nclass==Ptg%nclass,"BSym%nclass/=Ptg%nclass!")

   do icls=1,BSym%nclass ! FIXME this is awful, should be done in a cleaner way.
     BSym%nelements(icls)=Ptg%class_ids(2,icls) - Ptg%class_ids(1,icls) + 1
   end do

   ! 2) Generate the symmetry operations in the conventional vector coordinates.
   conv_gprimd(:,1)=bravais(3:5)
   conv_gprimd(:,2)=bravais(6:8)
   conv_gprimd(:,3)=bravais(9:11)

   axes = conv_gprimd
   call matr3inv(conv_gprimd,axes) !; axes=TRANSPOSE(axes)

   conv_gprimd=MATMUL(Cryst%gprimd,TRANSPOSE(axes))
   !conv_gprimd=MATMUL(axes,Cryst%gprimd)
   !conv_gprimd=MATMUL(TRANSPOSE(axes),Cryst%gprimd)
   !write(std_out,*)"conv_gprimd:", conv_gprimd

   ptgroup_name = ADJUSTL(ptgroup)

   select case (ptgroup_name)

   case ("3m","-3m") 
     call wrtout(std_out," Changing the conventional cell: rhombohedral --> triple hexagonal","COLL")
     ! Transformation matrices: primitive rhombohedral --> triple hexagonal cell obverse setting. Table 5.1.3.1 ITA page 81.
     pmat1 = RESHAPE( (/ 1,-1, 0, 0, 1,-1, 1, 1, 1/), (/3,3/) ) ! R1
     pmat2 = RESHAPE( (/ 0, 1,-1,-1, 0, 1, 1, 1, 1/), (/3,3/) ) ! R2
     pmat3 = RESHAPE( (/-1, 0, 1, 1,-1, 0, 1, 1, 1/), (/3,3/) ) ! R3
     pmat4 = RESHAPE( (/-1, 1, 0, 0,-1, 1, 1, 1, 1/), (/3,3/) ) ! R1 reverse setting.
     pmat5 = RESHAPE( (/ 0,-1, 1, 1, 0,-1, 1, 1, 1/), (/3,3/) ) ! R2 reverse setting.
     pmat6 = RESHAPE( (/ 1, 0,-1,-1, 1, 0, 1, 1, 1/), (/3,3/) ) ! R3 reverse setting.
     conv_gprimd = MATMUL(conv_gprimd,pmat1)
     !conv_gprimd = MATMUL(conv_gprimd,pmat2)
     !conv_gprimd = MATMUL(conv_gprimd,pmat3)
     !conv_gprimd = MATMUL(conv_gprimd,pmat4)
     !conv_gprimd = MATMUL(conv_gprimd,pmat5)
     !conv_gprimd = MATMUL(conv_gprimd,pmat6)
     !write(std_out,*)" New conv_gprimd:", conv_gprimd

   case ("mm2")
     call wrtout(std_out," Changing the conventional cell: unconventional orthorhombic setting --> conventional","COLL")
     ! Transformation matrices: unconvential orthorhombic --> conventional orthorhombic. Table 5.1.3.1 ITA page 81.
     pmat1 = RESHAPE( (/ 0, 1, 0, 1, 0, 0, 0, 0,-1/), (/3,3/) )  ! ( b, a,-c) --> (a,b,c)
     pmat2 = RESHAPE( (/ 0, 1, 0, 0, 0, 1, 1, 0, 0/), (/3,3/) )  ! ( c, a, b) --> (a,b,c)
     pmat3 = RESHAPE( (/ 0, 0, 1, 0, 1, 0,-1, 0, 0/), (/3,3/) )  ! (-c, b, a) --> (a,b,c)
     pmat4 = RESHAPE( (/ 0, 0, 1, 1, 0, 0, 0, 1, 0/), (/3,3/) )  ! ( b, c, a) --> (a,b,c)
     pmat5 = RESHAPE( (/ 1, 0, 0, 0, 0, 1, 0,-1, 0/), (/3,3/) )  ! ( a,-c, b) --> (a,b,c)
     conv_gprimd = MATMUL(conv_gprimd,pmat2)
     !write(std_out,*)" New conv_gprimd:", conv_gprimd
   case default
     continue
   end select

   ABI_ALLOCATE(conv_symrec,(3,3,BSym%nsym_gk))
   conv_symrec = sgk

   !axes=zero; axes(1,1)=one ; axes(2,2)=one ; axes(3,3)=one
   !call symrelrot(BSym%nsym_gk,conv_gprimd,axes,conv_symrec,tolsym)
   call symrelrot(BSym%nsym_gk,Cryst%gprimd,conv_gprimd,conv_symrec,tolsym)

   ! 3) Reorder indeces such that symmetries are packed in classes.
   ABI_ALLOCATE(found,(BSym%nsym_gk))
   ABI_ALLOCATE(new_idx,(BSym%nsym_gk))
   ABI_ALLOCATE(new_g0,(3,BSym%nsym_gk))
   new_g0=0; found=0

   do isym=1,BSym%nsym_gk
     do jsym=1,BSym%nsym_gk
       if (ALL(Ptg%sym(:,:,isym) == conv_symrec(:,:,jsym) ))  then
         found(isym)    = found(isym) + 1
         new_idx(isym)  = BSym%sgk2symrec(jsym)
         new_g0(:,isym) = BSym%g0(:,jsym)
         !EXIT
       end if
     end do
   end do
   !
   ! DEBUGGING SECTION
   !do isym=1,BSym%nsym_gk
   !  jsym=BSym%sgk2symrec(isym)
   !  call print_symmetries(1,Cryst%symrec(:,:,jsym),Cryst%tnons(:,jsym),Cryst%symafm(jsym))
   !  write(std_out,*)BSym%g0(:,isym)
   !end do
     
   if ( Ptg%nsym/=BSym%nsym_gk .or. ANY(found/=1) ) then 
     !write(std_out,*)Ptg%nsym, BSym%nsym_gk
     !write(std_out,'(a,(i2))')" found = ",found
     write(std_out,*)" Ptg%sym list, conv_symrec list,  found Ptg% "
     do isym=1,Ptg%nsym
       write(std_out,'(a,i2,a,9i2,4x,a,9i2)')" found ",found(isym)," Ptg ",Ptg%sym(:,:,isym),"conv_symrec ",conv_symrec(:,:,isym)
     end do
     msg = " sgk and BSym%Ptg are inconsistent. Check tables or source"
     MSG_WARNING(msg)
     BSym%err_msg = msg
     BSym%err_status = BSYM_PTG_WRONG_MAPPING
     BSym%has_chtabs = .FALSE.

   else ! Reorder symmetries. 
     BSym%sgk2symrec = new_idx
     BSym%g0 = new_g0
   end if
   
   ABI_DEALLOCATE(new_idx)
   ABI_DEALLOCATE(new_g0)
   ABI_DEALLOCATE(found)
   ABI_DEALLOCATE(conv_symrec)

   if (BSym%has_chtabs) then
     ! Multiply the point group irreps by e^{-ik.\tau} to have the irreps of the little group.
     ! Store the results in BSym%Ref_irreps so that one can classify the states afterwards.
     ABI_ALLOCATE(BSym%Ref_irreps,(BSym%nclass))
     ABI_ALLOCATE(phase_mkt,(BSym%nsym_gk))

     do isym_gk=1,BSym%nsym_gk
       isym =  BSym%sgk2symrec(isym_gk)
       mkt = -two_pi * DOT_PRODUCT(BSym%kpt, Cryst%tnons(:,isym))
       phase_mkt(isym_gk) = CMPLX(DCOS(mkt), DSIN(mkt))
     end do

     call copy_irrep(Ptg%Irreps,BSym%Ref_irreps,phase_mkt)
     ABI_DEALLOCATE(phase_mkt)
   end if

#if 0
   ! Herring test requires the evaluation of the expression:
   !
   !   sum_{S,\tau} \chi^{k,\alpha} ({S|\tau}^2)
   !
   ! where Sk = -k + g0, and \chi is the trace of the \alpha-th 
   ! irreducible representation of the little group of k.
   ! \chi^{k,\alpha} = e^{-ik.\tau} \chi(\alpha) provided that
   ! we are not at zone border with a non-symmorphic operation.
   ! The expression is always real and it can only be equal to \pm Ptg%nsym or zero.
   ! FIXME this part has to be rewritten from scratch.
   !if (BSym%err_status/=BSYM_NOERROR) then
   !  write(std_out,*)" Skipping Herring test"
   !  goto 110
   !end if

   if (Bsym%can_use_tr) then
     ABI_ALLOCATE(her_test,(BSym%nclass))

     ABI_ALLOCATE(tr_conv_symrec,(3,3,BSym%nsym_trgk))
     do isym_tr=1,BSym%nsym_trgk
       isym = BSym%tr_sgk2symrec(isym_tr)
       tr_conv_symrec(:,:,isym_tr)=Cryst%symrec(:,:,isym)
     end do

     call symrelrot(BSym%nsym_trgk,Cryst%gprimd,conv_gprimd,tr_conv_symrec_tr,tolsym)

     do isym_tr=1,BSym%nsym_trgk
       isym = BSym%tr_sgk2symrec(isym_tr)
       !rot2 = MATMUL(tr_sgk(:,:,isym),tr_sgk(:,:,isym))
       !tau2 = MATMUL(tr_sgk(:,:,isym),Cryst%tnons(:,isym)) + Cryst%tnons(:,isym)

       rot2 = MATMUL(tr_conv_symrec(:,:,isym_tr),tr_conv_symrec(:,:,isym_tr))
       tau2 = MATMUL(tr_conv_symrec(:,:,isym_tr),Cryst%tnons(:,isym)) + Cryst%tnons(:,isym)

       phase_k = EXP(-j_dpc*two_pi*DOT_PRODUCT(kpoint,tau2))
       call locate_sym(Ptg,rot2,isym2,icls2)

       do irp=1,BSym%nclass
         her_test(irp) = her_test(irp) + phase_k * Ptg%Irreps(irp)%trace(icls2)
       end do
     end do

     ABI_DEALLOCATE(tr_conv_symrec)

     ! FIXME
     ABI_ALLOCATE(BSym%herring_test,(BSym%nclass))

     do irp=1,BSym%nclass
       if ( ABS(her_test(irp) - Ptg%nsym) < tol6 ) then
         BSym%herring_test(irp) = +1
       else if ( ABS(her_test(irp)) < tol6 ) then
         BSym%herring_test(irp) =  0
       else if ( ABS(her_test(irp) + Ptg%nsym) < tol6 ) then
         BSym%herring_test(irp) = -1
       else
         write(msg,'(a,i2,2a,i0,a,i2)')&
&          " Herring test for the irreducible representation number ",irp,ch10,&
&          " gave ",BSym%herring_test(irp),", while it should be 0 or +- ",Ptg%nsym 
          MSG_WARNING(msg)
          BSym%err_msg   =msg
          BSym%err_status=BSYM_HERRING_WRONG_TEST
       end if
     end do

     ABI_DEALLOCATE(her_test)
   end if ! can_use_tr
#endif
   !
   ! Final check
   !allocate(mtab(BSym%nsym_gk,BSym%nsym_gk))
   !call mult_table(BSym%nsym_gk,Ptg%sym,mtab)

   !do isym=1,BSym%nsym_gk
   !  isym1 = BSym%sgk2symrec(isym)  
   !  do jsym=1,BSym%nsym_gk
   !    isym2 = BSym%sgk2symrec(jsym)  
   !    rot2 = MATMUL(Cryst%symrec(:,:,isym1),Cryst%symrec(:,:,isym2))

   !    iprod = mtab(isym,jsym)

   !    do irp=1,BSym%nclass
   !       dim_irrep = Ptg%Irreps(irp)%dim 
   !       allocate(mat_test(dim_irrep,dim_irrep))
   !       mat_test = Ptg%Irreps(irp)%mat(:,:,isym) * Ptg%Irreps(irp)%mat(:,:,jsym)
   !       !call locate_sym(Ptg,rot2,isym2,icls2)
   !       write(std_out,*)mat_test - Ptg%Irreps(irp)%mat(:,:,iprod)
   !       deallocate(mat_test)
   !    end do
   !                                                                                  
   !  end do
   !end do
   !                                                                                     
   !deallocate(mtab)
 end if

 ABI_DEALLOCATE(sgk)
 ABI_DEALLOCATE(tr_sgk)
 ABI_DEALLOCATE(elements_idx)

 !% allocate(BSym%irrep2b(0:BSym%nclass))
 !% call nullify_coeff(BSym%irrep2b)
 !
 ! 1) Allocate space for the irreducible representations.

 ! 2) Try to determine if we are in presence of an accidental degeneracy. Sufficient condition: 
 !    There exists a set of degenerate states whose dimension is greater than the dimension 
 !    of the irreducible representations of the point group. The check can be done only 
 !    if Character tables are available.

 if (BSym%has_chtabs) then
   ABI_ALLOCATE(dim_irreps,(BSym%nclass))
   dim_irreps = (/(BSym%Ref_irreps(irp)%dim, irp=1,BSym%nclass)/)
 end if

 nacc_deg=0
 ABI_ALLOCATE(BSym%degs_dim,(BSym%ndegs))
 ABI_ALLOCATE(BSym%Calc_irreps,(BSym%ndegs))

 if (BSym%can_use_tr)  then
   ABI_ALLOCATE(BSym%trCalc_irreps,(BSym%ndegs))
 end if

 do idg=1,Bsym%ndegs
   dim_degs=Bsym%degs_bounds(2,idg)-Bsym%degs_bounds(1,idg)+1

   if (Bsym%has_chtabs) then
     if (ALL(dim_degs /= dim_irreps)) then ! An accidental degeneracy is present.
       nacc_deg=nacc_deg+1 
     end if
   end if

   BSym%degs_dim(idg) = dim_degs

   call init_irrep(BSym%Calc_irreps(idg),BSym%nsym_gk,dim_degs)
   if (BSym%can_use_tr) call init_irrep(BSym%trCalc_irreps(idg),BSym%nsym_trgk,dim_degs)
 end do ! idg

 if (BSym%has_chtabs) then
   ABI_DEALLOCATE(dim_irreps)
   if (nacc_deg/=0) then
     write(msg,'(a,i0,a)')" Detected ",nacc_deg," accidental degeneracies."
     MSG_WARNING(msg)
     BSym%err_status=BSYM_ACCDEG_ERROR  ! TODO this should signal to the caller that we have to decompose the calculated representation.
     BSym%err_msg   =msg
   end if
 end if

 DBG_EXIT("COLL")

end subroutine init_bands_symmetries
!!***

!----------------------------------------------------------------------

!!****f* m_bands_sym/print_bands_symmetries
!! NAME
!! print_bands_symmetries
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      classify_bands
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine print_bands_symmetries(Bsym,unit,mode_paral,prtvol)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_bands_symmetries'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 type(bands_symmetries),intent(in) :: BSym

!Local variables-------------------------------
!scalars
 integer :: icl,idg,my_unt,my_prtvol
 integer :: irr_idx,nstates,nunknown,istart,istop,ii
 character(len=4) :: my_mode
 character(len=1000) :: fmt,msg,msg0 
!arrays

! *********************************************************************

 DBG_ENTER("COLL")

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 write(fmt,*)'(2a,3f8.4,3a,i4,2a,i3,2a,i2,2a,i2,a,',Bsym%nclass,'i2,a)'
 write(msg,fmt)ch10,&
&  ' ===== Character of bands at k-point: ',BSym%kpt,' ===== ',ch10,&
&  '   Total number of bands analyzed .................. ',BSym%nbnds,ch10,&
&  '   Number of degenerate sets detected .............. ',Bsym%ndegs,ch10,&
&  '   Number of operations in the little group of k ... ',Bsym%nsym_gk,ch10,&
&  '   Number of classes (irreps) in the group of k .... ',Bsym%nclass,' (',(BSym%nelements(icl),icl=1,BSym%nclass),' )' 
 call wrtout(my_unt,msg,my_mode)

 if (Bsym%nonsymmorphic_at_zoneborder) then 
   call wrtout(my_unt," Non-symmorphic small group at zone border. Character analysis not available ",my_mode)
 end if

 if (bsym_failed(Bsym)) then
   write(std_out,'(3a)')" WARNING: Band classification algorithm failed with the error:",ch10,&
&    TRIM(Bsym%err_msg)
   write(msg,'(3a)')" WARNING: Band classification algorithm failed with the error:",ch10,&
&    TRIM(Bsym%err_msg)
   call wrtout(my_unt,msg,my_mode)
 end if

 !nunknown=0
 !do iband=1,BSym%nbnds
 !  irr_idx = BSym%b2irrep(iband)
 !  if (irr_idx /= 0) then
 !    if (     BSym%has_chtabs) irr_name = BSym%Ref_Irreps(irr_idx)%name
 !    if (.not.BSym%has_chtabs) write(irr_name,'(i0)')irr_idx ! use the index instead of the name.
 !  else
 !    irr_name = "???"
 !    nunknown = nunknown +1
 !  end if 
 !  write(msg,'(a,i3,2a)')' Band ',iband,' belongs to irrep ',TRIM(irr_name)
 !  call wrtout(my_unt,msg,my_mode)
 !end do

 do irr_idx=1,Bsym%nclass
   nstates = BSym%irrep2b(irr_idx)%size
   if (BSym%has_chtabs) then
     write(msg0,'(a,i0,3a)')"  Found ",nstates," states with character ",TRIM(BSym%Ref_irreps(irr_idx)%name),": "
   else
     write(msg0,'(2(a,i0),a)')"   Found ",nstates," states with character index ",irr_idx,": "
   end if
   do istart=1,nstates,20
     istop=istart+11; if (istop>nstates) istop=nstates
     write(msg,'(20(1x,i0))')(BSym%irrep2b(irr_idx)%value(ii), ii=istart,istop)
     if (istart==1) msg = TRIM(msg0)//TRIM(msg)
     if (istart/=1) msg = "   "//TRIM(msg)
     call wrtout(my_unt,msg,my_mode)
   end do
 end do

 nunknown = BSym%irrep2b(0)%size
 if (nunknown > 0) then
   write(msg0,'(a,i0,a)')" WARNING: ",nunknown," states have not been classified:"
   do istart=1,nunknown,20
     istop=istart+11; if (istop>nunknown) istop=nunknown
     write(msg,'(20(1x,i0))')(BSym%irrep2b(0)%value(ii), ii=istart,istop)
     if (istart==1) msg = TRIM(msg0)//TRIM(msg)
     if (istart/=1) msg = "   "//TRIM(msg)
     call wrtout(my_unt,msg,my_mode)
   end do
 end if

 if (my_prtvol>0 .or. nunknown>0 .or. .not.BSym%has_chtabs) then ! print the calculated character table.
   call wrtout(my_unt,ch10//" Calculated character table ",my_mode)
   !write(fmt,*)'(i2,a,i2,1x,',Bsym%nclass,'(a,2f6.3),a)'
   write(fmt,*)'(i2,a,i2,1x,',Bsym%nclass,'(a,2f5.2),a)'
   do idg=1,Bsym%ndegs
     write(msg,fmt)&
&      BSym%degs_bounds(1,idg),'-',BSym%degs_bounds(2,idg),&
&      ('|',BSym%Calc_irreps(idg)%trace(BSym%nelements(icl)), icl=1,BSym%nclass),'|'
     call wrtout(my_unt,msg,my_mode)
   end do
 end if

 call flush_unit(my_unt)

 DBG_EXIT("COLL")

end subroutine print_bands_symmetries
!!***

!----------------------------------------------------------------------

!!****f* m_bands_sym/destroy_bands_symmetries_0d
!! NAME
!! destroy_bands_symmetries_0d
!!
!! FUNCTION
!!  Deallocate the memory allocated in the bands_symmetries datatype (scalar version)
!!
!! PARENTS
!!      m_bands_sym
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine destroy_bands_symmetries_0D(Bsym)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_bands_symmetries_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(bands_symmetries),intent(inout) :: BSym

!Local variables ------------------------------
!scalars
 integer :: ii
! *************************************************************************

 !@bands_symmetries
 if (associated(BSym%g0           ))  then
   ABI_DEALLOCATE(BSym%g0)
 end if
 if (associated(BSym%tr_g0        ))  then
   ABI_DEALLOCATE(BSym%tr_g0)
 end if
 if (associated(Bsym%nelements    ))  then
   ABI_DEALLOCATE(Bsym%nelements)
 end if
 if (associated(Bsym%sgk2symrec   ))  then
   ABI_DEALLOCATE(Bsym%sgk2symrec)
 end if
 if (associated(Bsym%tr_sgk2symrec))  then
   ABI_DEALLOCATE(Bsym%tr_sgk2symrec)
 end if
 if (associated(Bsym%herring_test ))  then
   ABI_DEALLOCATE(Bsym%herring_test)
 end if
 if (associated(BSym%b2irrep      ))  then
   ABI_DEALLOCATE(BSym%b2irrep)
 end if
 if (associated(BSym%degs_bounds  ))  then
   ABI_DEALLOCATE(BSym%degs_bounds)
 end if
 if (associated(BSym%degs_dim     ))  then
   ABI_DEALLOCATE(BSym%degs_dim)
 end if

 if (associated(BSym%irrep2b)) then
   do ii=LBOUND(BSym%irrep2b,DIM=1),UBOUND(BSym%irrep2b,DIM=1)
     ABI_DEALLOCATE(BSym%irrep2b(ii)%value)
   end do
   ABI_DEALLOCATE(BSym%irrep2b)
 end if

 if (associated(BSym%Calc_irreps)) then
   call destroy_irrep(BSym%Calc_irreps)
 end if

 if (associated(BSym%trCalc_irreps)) then
   call destroy_irrep(BSym%trCalc_irreps)
 end if

 if (associated(BSym%Ref_irreps)) then
   call destroy_irrep(BSym%Ref_irreps)
 end if

end subroutine destroy_bands_symmetries_0D
!!***

!----------------------------------------------------------------------

!!****f* m_bands_sym/destroy_bands_symmetries_2d
!! NAME
!! destroy_bands_symmetries_2d
!!
!! FUNCTION
!!  Deallocate the memory allocated in the bands_symmetries datatype (2D version)
!!
!! PARENTS
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine destroy_bands_symmetries_2D(Bsym)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_bands_symmetries_2D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(bands_symmetries),intent(inout) :: BSym(:,:)

!Local variables ------------------------------
 integer :: id1,id2
! *************************************************************************

 do id2=1,SIZE(BSym,DIM=2)
   do id1=1,SIZE(BSym,DIM=1)
     call destroy_bands_symmetries_0D(BSym(id1,id2))
   end do
 end do

end subroutine destroy_bands_symmetries_2D
!!***

!----------------------------------------------------------------------

!!****f* m_bands_sym/nullify_Bands_Symmetries_0D
!! NAME
!! nullify_Bands_Symmetries_0D
!!
!! FUNCTION
!!  Nullify all the pointers defined in the bands_symmetries datatype.
!!
!! PARENTS
!!      m_bands_sym
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine nullify_bands_symmetries_0D(BSym)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_bands_symmetries_0D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(bands_symmetries),intent(inout) :: Bsym
! *************************************************************************

 !@bands_symmetries
 nullify(BSym%g0           )
 nullify(BSym%tr_g0        )
 nullify(BSym%nelements    )
 nullify(BSym%sgk2symrec   )
 nullify(BSym%tr_sgk2symrec)
 nullify(BSym%herring_test )
 nullify(BSym%b2irrep      )
 nullify(BSym%degs_bounds  )
 nullify(BSym%degs_dim     )

 nullify(BSym%irrep2b)

! types
 nullify(BSym%Calc_irreps  )
 nullify(BSym%trCalc_irreps)
 nullify(BSym%Ref_irreps   )

end subroutine nullify_bands_symmetries_0D
!!***

!----------------------------------------------------------------------

!!****f* m_bands_sym/nullify_Bands_Symmetries_2D
!! NAME
!! nullify_Bands_Symmetries_2D
!!
!! FUNCTION
!!  Nullify all the pointers defined in the bands_symmetries datatype.
!!
!! PARENTS
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine nullify_bands_symmetries_2D(BSym)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_bands_symmetries_2D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(bands_symmetries),intent(inout) :: Bsym(:,:)

!Local variables ------------------------------
 integer :: ik,is
! *************************************************************************

 do is=1,SIZE(Bsym,DIM=2)
   do ik=1,SIZE(Bsym,DIM=1)
     call nullify_bands_symmetries_0D(Bsym(ik,is))
   end do
 end do

end subroutine nullify_bands_symmetries_2D
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/finalize_bands_sym
!! NAME
!! finalize_bands_sym
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      classify_bands
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine finalize_bands_sym(Bsym,prtvol)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'finalize_bands_sym'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: prtvol
 type(bands_symmetries),intent(inout) :: BSym

!Local variables-------------------------------
 integer :: idg,ib1,ib2,idx,nunknown,dg_dim
 integer :: try,irep,nitems,istat,nseen,isn 
 integer :: isym,idg1,idg2,dim_mat,irr_idx2,irr_idx1
 real(dp),parameter :: TOL_TRACE=0.1_dp,TOL_ORTHO=0.1_dp,TOL_UNITARY=0.1_dp ! Large tolerance is needed to avoid problems. 
 !real(dp),parameter :: TOL_TRACE=0.01_dp,TOL_ORTHO=0.01_dp,TOL_UNITARY=0.01_dp ! Large tolerance is needed to avoid problems. 
 !real(dp),parameter :: TOL_TRACE=tol3,TOL_ORTHO=tol3,TOL_UNITARY=tol3 ! Large tolerance is needed to avoid problems. 
 real(dp) :: uerr,max_err
 complex(dpc) :: ctest
 logical :: isnew
 character(len=500) :: msg
!arrays
 integer,allocatable :: dims_seen(:)
 complex(dpc),allocatable :: traces_seen(:,:)
 complex(dpc),pointer :: trace(:)
 complex(dpc),pointer :: calc_mat(:,:),trace1(:),trace2(:)
 complex(dpc),allocatable :: cidentity(:,:)

! *************************************************************************

 !@bands_symmetries

 ! Each band is initialized as "Unknown".
 BSym%b2irrep = 0

 ! Force the matrices to be unitary.
 call polish_irreps(BSym%Calc_irreps)

 if (.not.BSym%has_chtabs) then

   write(msg,'(5a)')&
&    " Reference character table not available. ",ch10,&
&    " Symmetry analysis not available. Using heuristic method to classify the states.",ch10,&
&    " It might not work, especially if accidental degeneracies are present."
   MSG_WARNING(msg)
   !
   ! The simplest thing we can do here is using the calculated matrices to get the 
   ! character and comparing the results hoping everything is OK.
   ABI_ALLOCATE(traces_seen,(BSym%nsym_gk,BSym%ndegs))
   ABI_ALLOCATE(dims_seen,(BSym%ndegs))

   traces_seen=czero; nseen=1
   traces_seen(:,1) = BSym%Calc_irreps(1)%trace
   dims_seen(1)     = BSym%Calc_irreps(1)%dim

   do idg=2,BSym%ndegs
     dg_dim = BSym%Calc_irreps(idg)%dim
     trace => BSym%Calc_irreps(idg)%trace
     isnew=.TRUE.
     do isn=1,nseen
       if (ALL (ABS(trace - traces_seen(:,isn)) < TOL_TRACE) ) then
         isnew=.FALSE.; EXIT
       end if
     end do

     if (isnew) then
       nseen = nseen+1
       traces_seen(:,nseen) = trace
       dims_seen(nseen) = dg_dim
     end if
   end do

   if (nseen>BSym%nclass) then
     write(msg,'(3a)')&
&      " The number of different calculated traces is found to be greater than nclasses!",ch10,&
&      " Heuristic method clearly failed. Symmetry analysis cannot be performed."
     MSG_WARNING(msg)
     BSym%err_status = BSYM_HEUR_WRONG_NCLASSES
     BSym%err_msg    = msg

     do isn=1,nseen
       write(msg,'(a,i0)')" Representation: ",isn
       call wrtout(std_out,msg,"COLL")
       call print_arr(traces_seen(:,isn),max_r=BSym%nsym_gk,unit=std_out,mode_paral="COLL")
     end do

   else  ! It seems that the Heuristic method succeeded.
     do idg=1,BSym%ndegs
       ib1=BSym%degs_bounds(1,idg)
       ib2=BSym%degs_bounds(2,idg)
       trace => BSym%Calc_irreps(idg)%trace
       do isn=1,nseen
         if (ALL (ABS(trace - traces_seen(:,isn)) < TOL_TRACE) ) then
           BSym%b2irrep(ib1:ib2)=isn
           if (BSym%Calc_irreps(idg)%dim /= dims_seen(isn)) then
             write(msg,'(3a)')&
&              " Found two set of degenerate states with same character but different dimension!",ch10,&
&              " heuristic method clearly failed. Symmetry analysis cannot be performed."
             MSG_ERROR(msg)
             BSym%err_status = BSYM_HEUR_WRONG_DIMS
             BSym%err_msg    = msg
           end if
           EXIT
         end if
       end do
     end do
   end if

   ABI_DEALLOCATE(traces_seen)
   ABI_DEALLOCATE(dims_seen)

 else 
   !
   ! * Search in the lookup table definining the irreducible representation 
   nunknown = 0
   do idg=1,BSym%ndegs

     ib1=BSym%degs_bounds(1,idg)
     ib2=BSym%degs_bounds(2,idg)
     trace => BSym%Calc_irreps(idg)%trace

     try = which_irrep(BSym, trace, tol3)
     if (try==0) try = which_irrep(BSym, trace, 0.1_dp) ! try again with increased tolerance.
     if (try/=0) then
       BSym%b2irrep(ib1:ib2)=try
     else
       BSym%b2irrep(ib1:ib2)=0
       nunknown = nunknown + (ib2-ib1+1)
     end if
   end do
 end if
 !
 ! %irrep2b(0)) gives the indeces of the states that have not been classified.
 ABI_ALLOCATE(BSym%irrep2b,(0:BSym%nclass))

 !write(std_out,*)"b2irrep",BSym%b2irrep

 do irep=0,BSym%nclass
   nitems = COUNT(BSym%b2irrep==irep)
   ABI_ALLOCATE(BSym%irrep2b(irep)%value,(nitems))
   istat = ABI_ALLOC_STAT
   BSym%irrep2b(irep)%size=nitems
   idx=0
   do ib1=1,BSym%nbnds
     if (BSym%b2irrep(ib1) == irep) then
       idx = idx + 1
       BSym%irrep2b(irep)%value(idx) = ib1
     end if
   end do
 end do

 if (BSym%irrep2b(0)%size /= 0) then
   write(msg,'(a,i0,a)')" Band classification algorithm was not able to classify ",BSym%irrep2b(0)%size," states."
   MSG_WARNING(msg)
   BSym%err_status = BSYM_CLASSIFICATION_ERROR
   BSym%err_msg    = msg
 end if
 !
 ! ==============================================================
 ! ==== Test basic properties of irreducible representations ====
 ! ==============================================================

 if (.not.bsym_failed(Bsym)) then 
   !
   ! 1) \sum_R \chi^*_a(R)\chi_b(R)= N_R \delta_{ab} 
   !
   !call wrtout(std_out," \sum_R \chi^*_a(R)\chi_b(R) = N_R \delta_{ab} ","COLL")
   max_err=zero
   do idg2=1,Bsym%ndegs
     trace2 => BSym%Calc_irreps(idg2)%trace(1:BSym%nsym_gk)
     ib2 = BSym%degs_bounds(1,idg2)
     irr_idx2 = BSym%b2irrep(ib2)
     if (irr_idx2 == 0) CYCLE

     do idg1=1,idg2
       trace1 => BSym%Calc_irreps(idg1)%trace(1:BSym%nsym_gk)
       ib1 = BSym%degs_bounds(1,idg1)
       irr_idx1 = BSym%b2irrep(ib1)
       if (irr_idx1 == 0) CYCLE
       ctest=DOT_PRODUCT(trace1,trace2)/Bsym%nsym_gk
       if (irr_idx1==irr_idx2) ctest=ctest-one
       max_err = MAX(max_err,ABS(ctest))
       if (.FALSE..and.ABS(ctest)>tol3) then
         write(msg,'(a,4i3,2es16.8)')&
&          ' WARNING: should be delta_ij: cx1 cx2, irr1, irr2, ctest: ',idg1,idg2,irr_idx1,irr_idx2,ctest
         call wrtout(std_out,msg,"COLL")
       end if
     end do
   end do

   if (max_err>TOL_ORTHO) then
     write(msg,'(a,es10.2)')" Too large maximum error on \sum_R \chi^*_a(R)\chi_b(R) = N_R \delta_{ab}: ",max_err
     MSG_WARNING(msg)
     Bsym%err_status =  BSYM_ORTHO_ERROR
     Bsym%err_msg    =  msg
   else
     write(msg,'(a,es10.2)')" maximum error on \sum_R \chi^*_a(R)\chi_b(R) = N_R \delta_{ab}: ",max_err
     call wrtout(std_out,msg,"COLL")
   end if

   if (.not.Bsym%only_trace) then 
     !call wrtout(std_out," **** Testing the unitary of the calculated irreps ****",my_mode)
     max_err=zero
     do idg1=1,Bsym%ndegs
       ib1 = BSym%degs_bounds(1,idg1)
       irr_idx1 = BSym%b2irrep(ib1)
       if (irr_idx1 == 0) CYCLE

       do isym=1,Bsym%nsym_gk
         calc_mat => BSym%Calc_irreps(idg1)%mat(:,:,isym)
         dim_mat  =  BSym%Calc_irreps(idg1)%dim
         ABI_ALLOCATE(cidentity,(dim_mat,dim_mat))
         call set2unit(cidentity) 
         uerr = MAXVAL( ABS(MATMUL(calc_mat,TRANSPOSE(DCONJG(calc_mat))) - cidentity) )
         max_err = MAX(max_err,uerr)
         ABI_DEALLOCATE(cidentity)
         if (.FALSE..and.prtvol>=10) then 
           write(std_out,'(a,i3,a,i2,a,es16.8,a)')&
&          " === idg: ",idg1,", isym: ",isym,", Error on U^* U = 1: ",uerr," ==="
           call print_arr(calc_mat,dim_mat,dim_mat,unit=std_out,mode_paral="COLL")
         end if
       end do
     end do

     if (max_err>TOL_UNITARY) then
       write(msg,'(a,es10.2)')" Too large maximum error on the unitary of representions matrices: ",max_err
       MSG_WARNING(msg)
       Bsym%err_msg    = msg
       Bsym%err_status = BSYM_UNITARY_ERROR
     else 
       write(msg,'(a,es10.2)')" maximum error on the unitary of representions matrices: ",max_err
       call wrtout(std_out,msg,"COLL")
     end if

   end if

 end if 

end subroutine finalize_bands_sym
!!***

!----------------------------------------------------------------------

!!****f* m_bands_sym/which_irrep
!! NAME
!!  m_bands_sym
!!
!! FUNCTION
!!  Return the index of the irreducible representation with character charact. 0 if not found.
!!
!! INPUTS
!!  BSym<bands_symmetries> 
!!  trace(%nsym_gk)=The trace of the representation to be compared with the internal database (if present).
!!  tolerr=Absolute error on the character.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function which_irrep(BSym,trace,tolerr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'which_irrep'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: which_irrep
 real(dp),intent(in) :: tolerr
 type(bands_symmetries),intent(in) :: BSym
!arrays
 complex(dpc),intent(in) :: trace(BSym%nsym_gk)

!Local variables-------------------------------
!scalars
 integer :: irp

! *********************************************************************

 which_irrep = 0
 if (BSym%has_chtabs) then ! Symmetry analysis can be performed.
   do irp=1,BSym%nclass
     if ( ALL( ABS(BSym%Ref_irreps(irp)%trace(:) - trace(:)) < tolerr)) then
       which_irrep = irp; EXIT
     end if
   end do
 end if

end function which_irrep
!!***

!----------------------------------------------------------------------


!!****f* m_bands_sym/symmetrize_me
!! NAME
!!  symmetrize_me       
!!
!! FUNCTION
!!
!! INPUTS
!!  BSym<bands_symmetries> 
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cohsex_me
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE

subroutine symmetrize_me(BSym,lbnd,ubnd,in_me,out_me)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symmetrize_me'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: lbnd,ubnd
 type(bands_symmetries),intent(in) :: BSym
!arrays
 complex(dpc),intent(in) :: in_me(2,lbnd:ubnd,lbnd:ubnd)
 complex(dpc),intent(out) :: out_me(lbnd:ubnd,lbnd:ubnd)

!Local variables-------------------------------
!scalars
 integer :: idg1,b1_start,b1_stop,irp1
 integer :: idg2,b2_start,b2_stop,irp2
 integer :: ii,jj,ib,jb,kk,kb,lb,ll
 complex(dpc) :: tr_ofd,ofd,dsd,tr_dsd
 type(irrep_t),pointer :: Irrep1,Irrep2
 type(irrep_t),pointer :: tr_Irrep1,tr_Irrep2

! *********************************************************************

 if (bsym_failed(Bsym)) then
   MSG_ERROR("Symmetrization cannot be performed. You should not be here!")
 end if

 do idg1=1,Bsym%ndegs  ! First loop over set of degenerate states.
   b1_start = Bsym%degs_bounds(1,idg1)  
   b1_stop  = Bsym%degs_bounds(2,idg1) 

   !if (b1_stop<lbnd .or. b2_start >ubnd) then
   !  MSG_ERROR("Wrong band indeces, check Bsym initialization")
   !end if

   Irrep1 => Bsym%Calc_irreps(idg1)
   if (Bsym%can_use_tr) tr_Irrep1 => Bsym%trCalc_irreps(idg1)
   irp1 = Bsym%b2irrep(b1_start)

   do idg2=1,Bsym%ndegs ! Second loop over set of degenerate states.
     !write(std_out,*)" ==> Symmetrizing degenerate set ",idg1,idg2
     b2_start = Bsym%degs_bounds(1,idg2)
     b2_stop  = Bsym%degs_bounds(2,idg2)
     irp2 = Bsym%b2irrep(b2_start)

     if (irp1/=irp2 .or. idg1==idg2) CYCLE  ! Skip diago elements or elements belonging to different irreps.

     Irrep2 => Bsym%Calc_irreps(idg2)
     if (Bsym%can_use_tr) tr_Irrep2 => Bsym%trCalc_irreps(idg2)
     !
     ! Symmetrize the off-diagonal matrix elements.
     ! summing over kk and ll. ii and jj are the indeces of the bands that are symmetrized 
     do ii=1,b1_stop-b1_start+1
       ib= ii+b1_start-1
       do jj=1,b2_stop-b2_start+1
         jb= jj+b2_start-1
         !write(std_out,*)" ====> Symmetrizing ",ib,jb

         ofd= czero; tr_ofd=czero
         do kk=1,b1_stop-b1_start+1
           kb= kk+b1_start-1
           do ll=1,b2_stop-b2_start+1
             lb= ll+b2_start-1
             dsd = sum_irreps(Irrep1,Irrep2,kk,ii,ll,jj)
             ofd = ofd + dsd * in_me(1,kb,lb) 
             if (Bsym%can_use_tr) then
               tr_dsd = sum_irreps(tr_Irrep1,tr_Irrep2,kk,jj,ll,ii) ! Exchange of band indeces.
               tr_ofd = tr_ofd + tr_dsd * in_me(2,kb,lb)            ! Contribution obtained from TR.
             end if
           end do
         end do

         out_me(ib,jb)= ofd/Bsym%nsym_gk
         if (Bsym%can_use_tr .and. Bsym%nsym_trgk>0) out_me(ib,jb)= out_me(ib,jb) + tr_ofd/Bsym%nsym_trgk
       end do
     end do
   end do
 end do

end subroutine symmetrize_me
!!***

!----------------------------------------------------------------------

!!****f* m_bands_sym/bsym_failed
!! NAME
!!  bsym_failed
!!
!! FUNCTION
!!
!! INPUTS
!!  BSym<bands_symmetries> 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function bsym_failed(BSym) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bsym_failed'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical :: bsym_failed
 type(bands_symmetries),intent(in) :: BSym

! *********************************************************************

 bsym_failed = (BSym%err_status /= BSYM_NOERROR)

end function bsym_failed
!!***

!----------------------------------------------------------------------

!!****f* m_bands_sym/polish_irreps
!! NAME
!!  polish_irreps
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_bands_sym
!!
!! CHILDREN
!!      xgeev,xginv,zpotrf,ztrsm
!!
!! SOURCE
  
subroutine polish_irreps(Irreps)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polish_irreps'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(irrep_t),intent(inout) :: Irreps(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: ldvl1=1,ldvr1=1
 integer :: irp,sym,dim,ldvr,ii,ivec,jvec,info
 character(len=500) :: msg
!arrays 
 complex(dpc),allocatable :: vl(:,:),vr(:,:),vrm1(:,:),overlap(:,:)
 complex(dpc),allocatable :: cmat(:,:),eigval(:)

! *********************************************************************

 ! Eigen decomposition: A = V D V^{-1}.
 do irp=1,SIZE(Irreps)
   dim = Irreps(irp)%dim
   ABI_ALLOCATE(cmat,(dim,dim))
   ABI_ALLOCATE(eigval,(dim))
   ldvr=dim
   ABI_ALLOCATE(vl,(ldvl1,dim))
   ABI_ALLOCATE(vr,(ldvr,dim))
   ABI_ALLOCATE(vrm1,(dim,dim))
   ABI_ALLOCATE(overlap,(dim,dim))
   do sym=1,Irreps(irp)%nsym
     cmat = Irreps(irp)%mat(:,:,sym)
     call xgeev("No vectors","Vectors",dim,cmat,dim,eigval,vl,ldvl1,vr,ldvr)
     !
     ! Orthogonalize the eigenvectors using Cholesky orthogonalization.
     do jvec=1,dim
       do ivec=1,jvec
         overlap(ivec,jvec) = DOT_PRODUCT(vr(:,ivec),vr(:,jvec))  
       end do
     end do
     !
     ! 2) Cholesky factorization: overlap = U^H U with U upper triangle matrix.
     call ZPOTRF('U',dim,overlap,dim,info)
     if (info/=0)  then
       write(msg,'(a,i3)')' ZPOTRF returned info= ',info
       MSG_ERROR(msg)
     end if 
     !
     ! 3) Solve X U = Vr, on exit the Vr treated by this node is orthonormalized.
     call ZTRSM('R','U','N','N',dim,dim,cone,overlap,dim,vr,dim)

     !write(std_out,*)"After ortho",MATMUL(TRANSPOSE(CONJG(vr)),vr)

     vrm1 = vr
     call xginv(vrm1,dim)
     do ii=1,dim
       eigval(ii) = eigval(ii)/ABS(eigval(ii)) ! Rescale the eigevalues.
       vrm1(ii,:) =  eigval(ii) * vrm1(ii,:)
     end do
     Irreps(irp)%mat(:,:,sym) = MATMUL(vr,vrm1)
     Irreps(irp)%trace(sym) = get_trace(Irreps(irp)%mat(:,:,sym))
   end do
   ABI_DEALLOCATE(cmat)
   ABI_DEALLOCATE(eigval)
   ABI_DEALLOCATE(vl)
   ABI_DEALLOCATE(vr)
   ABI_DEALLOCATE(vrm1)
   ABI_DEALLOCATE(overlap)
 end do

end subroutine polish_irreps
!!***

!----------------------------------------------------------------------

END MODULE m_bands_sym

