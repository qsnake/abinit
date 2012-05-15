!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gsphere
!! NAME
!!  m_gsphere
!!
!! FUNCTION
!!  This module defines two objects:
!!
!!   1) The Gsphere data type defining the set of G-vectors
!!      centered on Gamma used to describe (chi0|epsilon|W) in the GW code.
!!      Note that, unlike the kg_k arrays used for wavefunctions, here the
!!      G-vectors are ordered in shells (increasing length). Moreover
!!      the sphere can be enlarged to take into account umklapps for which
!!      one need the knowledge of several quantities at G-G0.
!!
!!   2) The Gpairs_q object used to define the set of independent G1-G2 Fourier
!!      components of a two-point function which is invariant under the symmeties
!!      of the space group
!!
!!  Methods used to initialize and destroy these two objects are defined here.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MG, GMR, VO, LR, RWG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!! Change name (gvectors >>> Gsphere)
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

MODULE m_gsphere

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

 use m_gwdefs,        only : GW_TOLQ
 use defs_abitypes,   only : MPI_type
 use m_numeric_tools, only : bisect
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_structure

 implicit none

 private

! Low-level procedures.
 public :: merge_and_sort_kg   ! Merges a set of k-centered G-spheres of cutoff energy ecut and
                               ! returns a Gamma-centered G-spheres.

 public :: table_gbig2kg       ! Associate the kg_k set of G-vectors with Gamma-centered G-sphere.

 public :: get_kg              ! Helper function to calculate the set of G-vectors at a given kpoint.
!!***

!----------------------------------------------------------------------

!!****t* m_gsphere/gvectors_type
!! NAME
!! gvectors_type
!!
!! FUNCTION
!! The gvectors_type data type contains information related to the set of G vectors
!! used during a screening or a GW calculation, as well as symmetry tables relating
!! these vectors. Presently the following quantities are stored
!!
!! 1) The reduced coordinates of the G vectors (arrays gvec)
!! 2) Tables giving the correspondence between a G-vector and its rotated image
!!    through a symmetry operation in reciprocal space.
!! 3) List of the irreducible G pairs
!! 4) Tables giving, for each pair in the full reciprocal space, the corresponding
!!    irreducible pair as well as the symmetry operation in reciprocal space
!!
!! Note that, unlike the GS part, the basis set does not depend on the k-point.
!!
!! NOTES
!! To indicate the indices in the arrays grottb, grottbm1 we use the following notation :
!!
!!  g defines the index of the reciprocal lattice vector in the array gvec
!!  s  indicates the index of the symmetry operation in reciprocal space
!!  i  can  be one or two. 1 is used to indicate the identity operator
!!
!! SOURCE

 type,public :: gvectors_type

  integer :: ng
  ! Total number of G vectors in the sphere taking into account umklapps
  ! it defines the size of the array gvec and it accounts for possible umklapps for which
  ! one has to shift the sphere.

  ! TODO: The sphere should be enlarged in init_Gpairs_type using mg0 and (ecut|ng) as input.
  ! For the time being we keep the old implementation.
  !%integer :: ng_eff
  ! Effective number of G vectors, i.e. the number of G in the smaller sphere without umklapps.
  ! ng_eff<=ng and should be used to loop over the elements of (chi0|epsilon|W).
  !
  ! TODO: Add info on FFT including zero-padding algorithm.
  ! table must be recalculated for each G0 and rho_tw_g should accept gpshere in input.

  integer :: nsh                   ! Number of shells
  integer :: nsym                  ! The number of symmetry operations
  integer :: timrev                ! 2 if time-reversal is used, 1 otherwise
  integer :: istwfk=1              ! Storage mode. At present time-reversal is not used.

  !integer :: mg0(3)=0
  ! For each reduced direction gives the max G0 component to account for umklapp processes.

  real(dp) :: ecut
  ! Cutoff energy of the sphere.

  real(dp) :: gmet(3,3)
  ! Reciprocal space metric ($\textrm{bohr}^{-2}$).

  real(dp) :: gprimd(3,3)
  ! Dimensional reciprocal space primitive translations (bohr^{-1})

  integer,pointer :: g2sh(:)    SET2NULL
  ! g2sh(ng)
  ! For each G, it gives the index of the shell to which it belongs.

  integer,pointer :: gvec(:,:)  SET2NULL
  ! gvec(3,ng)
  ! Reduced coordinates of G vectors.

  integer,pointer :: rottb(:,:,:)   SET2NULL
  ! rottb(ng,timrev,nsym)
  ! rottb(G,I,S) is the index of (SI) G in the array gvec
  ! where I is either the identity or the inversion.

  integer,pointer :: rottbm1(:,:,:)  SET2NULL
  ! rottb(ng,timrev,nsym)
  ! rottbm1(G,I,S) is the index of IS{^-1} G in the array gvec

  integer,pointer :: shlim(:)  SET2NULL
  ! shlim(nsh+1)
  ! Index of the first G vector in each shell, =ng+1 for nsh+1

  real(dp),pointer :: shlen(:)   SET2NULL
  ! shlen(nsh)
  ! Radius of each shell.

  !TODO switch to dpc
  complex(gwpc),pointer :: phmGt(:,:)  SET2NULL
  ! phmGt(ng,nsym)
  ! Phase factor e^{-i2\pi(G.\tau)} where $\tau$ is the fractional translation associated to isym.

  complex(gwpc),pointer :: phmSGt(:,:)  SET2NULL
  ! phmSGt(ng,nsym)
  ! Phase factor e^{-i2\pi(SG.\tau)} where S is one of the symmetry properties in reciprocal space.

 end type gvectors_type
!!***

 public :: init_gsphere       ! Initialize the G-sphere.
 public :: gsph_fft_tabs      ! Returns useful tables for FFT (with or without padding).
 public :: gsph_in_fftbox     ! Initialize the largest Gsphere contained in the FFT box.
 public :: nullify_gsphere    ! Nullify all pointers.
 public :: print_gsphere      ! Printout of basic dimensions.
 public :: destroy_gsphere    ! Free dynamics memory allocated in the object.
 public :: gsph_g_idx         ! Returns the index of G from its reduced coordinates.
 public :: gsph_gmg_idx       ! Returns the index of G1-G2 from their indeces
 public :: gsph_gmg_fftidx    ! Returns the index of G1-G2 in the FFT mesh defined by ngfft.

!----------------------------------------------------------------------

!!****t* m_gsphere/Gpairs_type
!! NAME
!! Gpairs_type
!!
!! FUNCTION
!! The Gpairs_type data type contains information useful to symmetrize in reciprocal
!! space any two-point function, f_q(G1,G2), which has the same symmetry of the crystal.
!! In particular the structure contains:
!!
!! 1) The List of the irreducible G pairs
!! 2) Tables giving, for each G1 G2 pair in reciprocal space, the corresponding
!!    irreducible pair and the reqired symmetry.
!!
!! SOURCE

 type,public :: Gpairs_type

  integer :: ng                      ! The number of G vectors
  integer :: ngpairs                 ! Total number of G pairs ie ng**2, used to dimension gptab and gptabo
  integer :: niggp                   ! Number of irreducible G-Gp pairs for this q-point, used to dimension gptab and gptabo
  integer :: nsym                    ! Number of operations
  integer :: timrev                  ! If time-reversal has been considered to generate the G-sphere
  logical :: can_use_timrev          ! .TRUE. If time-reversal can be used (actually only at Gamma)

  integer,pointer :: fp2ip(:,:,:)   SET2NULL
   ! fp2ip(2,T1,T2) gives the sequential index, in the array gvec, of the
   ! indendent pair (G1,G2) such as (T1,T2)= S (G1,G2)

  integer,pointer :: fptabo(:,:)    SET2NULL
   ! fptabo(ng,ng)=index of the symmetry operation in the array symrec
   ! such as (T1,T2)= S (G1,G2). At Gamma, if time-reversal is used, the index is negative

  integer,pointer :: ip2fp(:,:)     SET2NULL
   ! ip2fp(2,niggp)= index of G1 and G2 in the array gvec for each niggp independent pair.

  real(dp) :: qpt(3)
  ! The point in the BZ where the two-point function has to be symmetrized

 end type Gpairs_type
!!***

!----------------------------------------------------------------------

 public :: init_Gpairs_type       ! Init the G-pairs object.
 public :: destroy_Gpairs_type    ! Free the memory allocate in the G-pairs object.
 public :: nullify_Gpairs_type
 public :: prune_g1mg2            ! Evalute all possible differences G1-G2, remove duplicated differences and
                                  ! return the list of different differences.


CONTAINS  !=========================================================================================================================
!!***

!!****f* m_gsphere/setup_G_rotation
!! NAME
!! setup_G_rotation
!!
!! FUNCTION
!! Set up tables indicating rotation of G-vectors.
!!
!! INPUTS
!! only_one_kpt=TRUE if only the Gamma point is sampled [FIXME to be removed]
!! nsym=Number of symmetry operations.
!! symrec(3,3,nsym)=Symmetry operations in reciprocal space.
!! timrev=2 if time reversal can be used, 1 otherwise.
!! npw=Number of planewaves in the sphere.
!! gvec(3,npw)=Coordinates of plane waves, supposed to be ordered in increasing modulus
!! g2sh(npw)=For each G, it gives the index of the shell to which it belongs.
!! nsh=Number of shells
!! shlim(nsh+1)=Index of the first G vector in each shell, =npw+1 for nsh+1
!!
!! OUTPUT
!!  grottb  (npw,2,nsym)= grottb(G,I,S) is the index of (SI) G in the array gvec.
!!  grottbm1(npw,2,nsym)= index of IS^{-1} G.
!!
!! NOTES:
!!  I is either the identity or the inversion (time reversal in reciprocal space).
!!  S is one of the symmetry operation in reciprocal space belonging to the Space group.
!!
!! PARENTS
!!      m_gsphere
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine setup_G_rotation(only_one_kpt,nsym,symrec,timrev,npw,gvec,g2sh,nsh,shlim,grottb,grottbm1)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_G_rotation'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nsh,nsym,timrev
 logical,intent(in) :: only_one_kpt
!arrays
 integer,intent(in) :: g2sh(npw),gvec(3,npw),shlim(nsh+1),symrec(3,3,nsym)
 integer,intent(inout) :: grottb  (npw,timrev,nsym)
 integer,intent(inout) :: grottbm1(npw,timrev,nsym)

!Local variables ------------------------------
!scalars
 integer :: ee,ig1,ig2,ish1,isym,itim,ss
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gbase(3),grot(3)

!************************************************************************
 !
 ! === Set up G-rotation table ===
 ! * This loop might be CPU consuming in isolated systems.
 !   Therefore we skip it in case of one single k-point.
 if (only_one_kpt) then
   ! As only_one_kpt is true, the only symmetry needed is identity
   do ig1=1,npw
     grottb  (ig1,1,1)=ig1
     grottbm1(ig1,1,1)=ig1
     !TODO check if inversion can enter somewhere!!!
   end do

 else  ! Several k-points.
  do ig1=1,npw
    ish1=g2sh(ig1) ; ss=shlim(ish1) ; ee=shlim(ish1+1)-1
    gbase(:)=gvec(:,ig1)

    do itim=1,timrev
      do isym=1,nsym
        grot=(3-2*itim)*MATMUL(symrec(:,:,isym),gbase)
        found=.FALSE.
        ! * Loop on the shell of ig1 to speed up the search.
        do ig2=ss,ee
          if (ALL(ABS(grot(:)-gvec(:,ig2))==0)) then
            found=.TRUE.
            grottb  (ig1,itim,isym)=ig2
            grottbm1(ig2,itim,isym)=ig1
          end if
        end do
        if (.not.found) then
          write(msg,'(3a,i5,a,i5,1x,2(3i5,a),a,i3,a,i3)')&
&          ' G-shell not closed',ch10,&
&          '  Initial G vector ',ig1,'/',npw,gbase(:),' Rotated G vector ',grot(:),ch10,&
&          '  Through sym ',isym,' and itim ',itim
          MSG_ERROR(msg)
        end if
      end do
    end do

  end do !ig1
 end if !only_one_kpt

end subroutine setup_G_rotation
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/init_gsphere
!! NAME
!! init_gsphere
!!
!! FUNCTION
!!  Main creation method for the Gvectors data type
!!
!! INPUTS
!!  only_one_kpt=TRUE if only the Gamma point is sampled [FIXME to be removed]
!!  Cryst<Crystal_structure> = Info on unit cell and its symmetries
!!     %nsym=number of symmetry operations
!!     %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!     %tnons(3,nsym)=fractional translations
!!     %gmet(3,3)=reciprocal space metric (bohr**-2).
!!     %gprimd(3,3)=dimensional reciprocal space primitive translations
!!  ng=number of G vectors, needed only if gvec is passed.
!!  [gvec(3,ng)]=coordinates of G vectors
!!  [ecut]=Cutoff energy for G-sphere. gvec and ecut are mutually exclusive.
!!
!! OUTPUT
!!  Gsph<gvectors_type>=Data type containing information related to the set of G vectors
!!   completetly initialized in output.
!!
!! NOTES
!!  gvec are supposed to be ordered with increasing norm.
!!
!! PARENTS
!!      m_gsphere,m_screening,mrgscr,setup_bse,setup_screening,setup_sigma
!!      sigma
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine init_gsphere(Gsph,only_one_kpt,Cryst,ng,gvec,ecut)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_gsphere'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng
 logical,intent(in) :: only_one_kpt
 real(dp),optional,intent(in) :: ecut
 type(Crystal_structure),intent(in) :: Cryst
 type(gvectors_type),intent(out) :: Gsph
!arrays
 integer,optional,intent(in) :: gvec(3,ng)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt1=1
 integer :: ig,istat,isym,nsh,nsym,timrev,pinv
 integer :: g1,g2,g3
 real(dp) :: eps,norm,norm_old,max_ecut,gsq
 logical :: ltest
!arrays
 real(dp),parameter :: k_gamma(3)=(/zero,zero,zero/)
 integer :: sg(3)
 integer,allocatable :: shlim(:)
 integer,pointer :: symrec(:,:,:)
 real(dp) :: kptns1(3,nkpt1)
 real(dp),allocatable :: shlen(:)
 real(dp),pointer :: tnons(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 !@gvectors_type
 call nullify_gsphere(Gsph)

 ! === Copy info on symmetries ===
 nsym   =  Cryst%nsym
 timrev =  Cryst%timrev
 symrec => Cryst%symrec
 tnons  => Cryst%tnons
 !
 ! === Initialize the object ===
 Gsph%istwfk = 1           ! Time reversal is not used here.
 Gsph%nsym   = nsym
 Gsph%timrev = timrev

 Gsph%gmet   = Cryst%gmet
 Gsph%gprimd = Cryst%gprimd

 if (PRESENT(gvec)) then
   if (PRESENT(ecut)) then
     MSG_BUG("ecut cannot be present when gvec is used")
   end if
   Gsph%ng= ng
   ABI_ALLOCATE(Gsph%gvec,(3,ng))
   Gsph%gvec=gvec
   !
   ! Calculate cutoff energy.of the sphere.
   max_ecut=-one
   do ig=1,ng
     g1=gvec(1,ig)
     g2=gvec(2,ig)
     g3=gvec(3,ig)
     gsq=       Cryst%gmet(1,1)*g1**2+Cryst%gmet(2,2)*g2**2+Cryst%gmet(3,3)*g3**2+ &
&          two*(Cryst%gmet(1,2)*g1*g2+Cryst%gmet(1,3)*g1*g3+Cryst%gmet(2,3)*g2*g3)
     max_ecut=MAX(max_ecut,gsq)
   end do
   max_ecut=two*max_ecut*pi**2
   Gsph%ecut= max_ecut

 else
   ! To be consistent with the previous implementation.
   MSG_WARNING("Init from ecut has to be tested")
   Gsph%ecut = ecut
   pinv=+1; kptns1(:,1)=k_gamma
   nullify(Gsph%gvec)
   call merge_and_sort_kg(nkpt1,kptns1,ecut,Cryst%nsym,pinv,Cryst%symrel,Cryst%gprimd,Gsph%gvec,0)
   Gsph%ng = SIZE(Gsph%gvec,DIM=2)
 end if
 !
 ! === Calculate phase exp{-i2\pi G.\tau} ===
 ABI_ALLOCATE(Gsph%phmGt,(Gsph%ng,nsym))
 do ig=1,Gsph%ng
   do isym=1,nsym
    Gsph%phmGt(ig,isym)=EXP(-j_dpc*two_pi*DOT_PRODUCT(Gsph%gvec(:,ig),tnons(:,isym)))
   end do
 end do
 !
 ! === Calculate phase phsgt= exp{-i2\pi SG\cdot t} ===
 ! TODO Here we can store only one of this arrays but I have to rewrite screeening!
 ABI_ALLOCATE(Gsph%phmSGt,(Gsph%ng,nsym))
 do ig=1,Gsph%ng
   do isym=1,nsym
     sg=MATMUL(symrec(:,:,isym),Gsph%gvec(:,ig))
     Gsph%phmSGt(ig,isym)=EXP(-j_dpc*two_pi*DOT_PRODUCT(sg,tnons(:,isym)))
   end do
 end do
 !
 ! === Calculate number of shells and corresponding starting index ===
 ! * Shells are useful to speed up search algorithms see e.g setup_G_rotation.
 ! * The last shell ends at ng+1, thus gvec is supposed to be closed.
 ltest=ALL(Gsph%gvec(:,1)==0)
 ABI_CHECK(ltest,'First G must be 0')

 ABI_ALLOCATE(Gsph%g2sh,(Gsph%ng))
 Gsph%g2sh(1)=1 ! This table is useful if we dont loop over shell

 ! For each shell, gives the index of the initial G-vector.
 ABI_ALLOCATE(shlim,(Gsph%ng+1))
 shlim(1)=1

 ! For each shell, gives the radius of the shell.
 ABI_ALLOCATE(shlen,(Gsph%ng))
 shlen(1)=zero

 nsh=1; norm_old=zero
 do ig=2,Gsph%ng
   norm=two_pi*SQRT(DOT_PRODUCT(Gsph%gvec(:,ig),MATMUL(Cryst%gmet,Gsph%gvec(:,ig))))
   eps=norm*tol8
   if (ABS(norm-norm_old)>eps) then
     norm_old=norm
     nsh=nsh+1
     shlim(nsh)=ig
     shlen(nsh)=norm
   end if
   Gsph%g2sh(ig)=nsh
 end do
 shlim(nsh+1)=Gsph%ng+1

 ! === Save info on the shells ===
 Gsph%nsh=nsh
 ABI_ALLOCATE(Gsph%shlim,(nsh+1))
 Gsph%shlim=shlim(1:nsh+1)
 ABI_ALLOCATE(Gsph%shlen,(nsh  ))
 Gsph%shlen=shlen(1:nsh)
 ABI_DEALLOCATE(shlim)
 ABI_DEALLOCATE(shlen)
 !
 ! === Calculate tables for rotated G"s ===
 ABI_ALLOCATE(Gsph%rottb  ,(Gsph%ng,timrev,nsym))
 istat = ABI_ALLOC_STAT
 ABI_ALLOCATE(Gsph%rottbm1,(Gsph%ng,timrev,nsym))
 istat = ABI_ALLOC_STAT

 call setup_G_rotation(only_one_kpt,nsym,symrec,timrev,Gsph%ng,Gsph%gvec,&
&  Gsph%g2sh,Gsph%nsh,Gsph%shlim,Gsph%rottb,Gsph%rottbm1)

 !call print_gsphere(Gsph,unit=std_out,prtvol=1)

 DBG_EXIT("COLL")

end subroutine init_gsphere
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_fft_tabs
!! NAME
!! gsph_fft_tabs
!!
!! FUNCTION
!!
!! INPUTS
!!  Gsph<gvectors_type>=Info on the G-sphere
!!  g0(3)
!!  mgfft=MAXVAL(ngfft(1:3))
!!  ngfftf(18)=Info on the FFT mesh.
!!
!! OUTPUT
!!  use_padfft=1 if padded FFT can be used, 0 otherwise.
!!  gmg0_gbound(2*mgfft+8,2)=Tables for improved zero-padded FFTS. Calculated only if use_padfft==1
!!  gmg0_ifft(Gsph%ng)=Index of G-G0 in the FFT mesh defined by ngfft.
!!
!! NOTES
!!  The routine will stop if any G-G0 happens to be outside the FFT box.
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,cchi0q0_intraband
!!      check_completeness,cohsex_me,exc_build_block,exc_build_ham
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine gsph_fft_tabs(Gsph,g0,mgfft,ngfft,use_padfft,gmg0_gbound,gmg0_ifft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsph_fft_tabs'
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft
 integer,intent(out) :: use_padfft
 type(gvectors_type),intent(in) :: Gsph
!arrays
 integer,intent(in) :: g0(3),ngfft(18)
 integer,intent(out) :: gmg0_gbound(2*mgfft+8,2),gmg0_ifft(Gsph%ng)

!Local variables-------------------------------
!scalars
 integer :: ig,ng,ierr
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer,allocatable :: gmg0(:,:)
 logical,allocatable :: kg_mask(:)

! *************************************************************************

 if (mgfft/=MAXVAL(ngfft(1:3))) then
   MSG_ERROR("mgfft/-MAXVAL(ngfft(1:3)")
 end if

 ng = Gsph%ng

 ierr=0; use_padfft=0
 ABI_ALLOCATE(gmg0,(3,ng))
 do ig=1,ng
   gmg0(:,ig) = Gsph%gvec(:,ig)-g0
   ! Consider possible wrap around errors.
   if ( ANY(gmg0(:,ig)>ngfft(1:3)/2) .or. ANY(gmg0(:,ig)<-(ngfft(1:3)-1)/2) ) then
    !gmg0_ifft(ig,ig01+mg0(1)+1,ig02+mg0(2)+1,ig03+mg0(3)+1) = 0
    write(std_out,*)" outside FFT box ",gmg0(:,ig)
    ierr=ierr+1
   end if
   if (ALL(gmg0(:,ig) == 0)) use_padfft=1
 end do

 if (ierr/=0) then
   write(msg,'(a,i0,a)')' Found ',ierr,' G-G0 vectors falling outside the FFT box. This is not allowed '
   MSG_ERROR(msg)
 end if
 !
 ! Evaluate the tables needed for the padded FFT performed in rhotwg. Note that we have
 ! to pass G-G0 to sphereboundary instead of G as we need FFT results on the shifted G-sphere,
 ! If Gamma is not inside G-G0 one has to disable FFT padding as sphereboundary will give wrong tables.
 if (use_padfft==1) call sphereboundary(gmg0_gbound,1,gmg0,mgfft,ng)

 call initmpi_seq(MPI_enreg_seq) ! No FFT parallelism.

 ABI_ALLOCATE(kg_mask,(ng))
 call kgindex(gmg0_ifft,gmg0,kg_mask,MPI_enreg_seq,ngfft,ng)
 ABI_CHECK(ALL(kg_mask),"FFT para not yet implemented")
 ABI_DEALLOCATE(kg_mask)

 ABI_DEALLOCATE(gmg0)

end subroutine gsph_fft_tabs
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_in_fftbox
!! NAME
!! gsph_fftbox
!!
!! FUNCTION
!!  Initialize the largest Gsphere contained in the FFT box.
!!
!! INPUTS
!!  Cryst<Crystal_structure> = Info on unit cell and its symmetries.
!!  ngfft(18)=Info on the FFT box.
!!
!! OUTPUT
!!  Gsph<gvectors_type>=Data type containing information related to the set of G vectors
!!   completetly initialized in output.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gsph_in_fftbox(Gsph,Cryst,ngfft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsph_in_fftbox'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(gvectors_type),intent(out) :: Gsph
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: dir1,dir2,dir3,npw,ig,ist
 real(dp) :: ecut,trial_ene
 logical :: only_one_kpt=.FALSE.
!arrays
 integer :: n1_max(3),n2_max(3),n3_max(3),vec(3)
 integer,allocatable :: gvec(:,:)

!************************************************************************

 ! Find ecut for the largest G-sphere contained in the FFT box.
 n1_max(1) = -(ngfft(1)-1)/2
 n2_max(1) = -(ngfft(2)-1)/2
 n3_max(1) = -(ngfft(3)-1)/2

 n1_max(2) = 0
 n2_max(2) = 0
 n3_max(2) = 0

 n1_max(3) = ngfft(1)/2
 n2_max(3) = ngfft(2)/2
 n3_max(3) = ngfft(3)/2

 ecut = HUGE(one)
 do dir3=1,3
   vec(3) = n1_max(dir3)
   do dir2=1,3
     vec(2) = n2_max(dir2)
     do dir1=1,3
       vec(1) = n1_max(dir1)
       if (ANY(vec/=0)) then
         trial_ene = half * normv(vec,Cryst%gmet,"G")**2
         ecut = MIN(ecut,trial_ene)
         !write(std_out,*)vec(:),trial_ene
       end if
     end do
   end do
 end do
 !
 ! Init sphere from ecut.
 call init_gsphere(Gsph,only_one_kpt,Cryst,0,ecut=ecut)
 !
 ! Make sure that Gsph does not contain G vectors outside the FFT box.
 ! kpgsph might return G whose energy is larger than the input ecut.
 npw = Gsph%ng
 star_loop: do ist=1,Gsph%nsh-1
   do ig=Gsph%shlim(ist),Gsph%shlim(ist+1)
     if ( ANY(Gsph%gvec(:,ig)>ngfft(1:3)/2) .or. ANY(Gsph%gvec(:,ig)<-(ngfft(1:3)-1)/2) ) then
       npw = Gsph%shlim(ist)-1  ! Gsph exceeds the FFT box. Only the shells up to npw will be used.
       EXIT star_loop
     end if
   end do
 end do star_loop

 if (npw<Gsph%ng) then
   MSG_COMMENT("Have to reinit Gpshere")
   ABI_ALLOCATE(gvec,(3,npw))
   gvec =Gsph%gvec(:,1:npw)
   call destroy_gsphere(Gsph)
   call init_gsphere(Gsph,only_one_kpt,Cryst,npw,gvec=gvec)
   ABI_DEALLOCATE(gvec)
 end if

end subroutine gsph_in_fftbox
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/print_gsphere
!! NAME
!! print_gsphere
!!
!! FUNCTION
!!  Print the content of a gvectors data type
!!
!! INPUTS
!!  Gsph<gvectors_type>=Info on the G-sphere
!!  unit=the unit number for output
!!  prtvol = verbosity level
!!  mode_paral =either "COLL" or "PERS"
!!
!! OUTPUT
!!
!! PARENTS
!!      cchi0q0,setup_bse
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine print_gsphere(Gsph,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_gsphere'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(gvectors_type),intent(in) :: Gsph

!Local variables-------------------------------
!scalars
 integer :: ish,nsc,my_unt,my_prtvol
 real(dp) :: fact,kin
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt    =std_out; if (PRESENT(unit      )) my_unt    =unit
 my_prtvol=0       ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode   ='COLL' ; if (PRESENT(mode_paral)) my_mode   =mode_paral

 write(msg,'(3a,2(a,i8,a))')ch10,&
& ' ==== Info on the G-sphere ==== ',ch10,&
& '  Number of G vectors ... ',Gsph%ng,ch10,&
& '  Number of shells ...... ',Gsph%nsh,ch10
 call wrtout(my_unt,msg,my_mode)

 SELECT CASE (Gsph%timrev)
 CASE (1)
   call wrtout(my_unt,' Time reversal symmetry cannot be used',my_mode)
 CASE (2)
   call wrtout(my_unt,' Time reversal symmetry is used',my_mode)
 CASE DEFAULT
   MSG_BUG("Wrong timrev")
 END SELECT

 if (my_prtvol/=0) then
   fact=half*two_pi**2
   write(msg,'(a)')
   call wrtout(my_unt,' Shell   Tot no. of Gs   Cutoff [Ha]',my_mode)
   do ish=1,Gsph%nsh
     nsc=Gsph%shlim(ish+1)-1
     kin=half*Gsph%shlen(ish)**2
     write(msg,'(2x,i4,10x,i6,5x,f8.3)')ish,nsc,kin
     call wrtout(my_unt,msg,'COLL')
   end do
   call wrtout(my_unt,ch10,my_mode)
 end if

end subroutine print_gsphere
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/destroy_gsphere
!! NAME
!! destroy_gsphere
!!
!! FUNCTION
!!  Deallocate all associated pointers defined in the structure.
!!
!! INPUTS
!!   Gsph = datatype to be freed
!!
!! OUTPUT
!!
!! PARENTS
!!      bethe_salpeter,cchi0,cchi0q0,check_completeness,m_gsphere,m_screening
!!      mrgscr,screening,setup_screening,sigma
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine destroy_gsphere(Gsph)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_gsphere'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(gvectors_type),intent(inout) :: Gsph

! *************************************************************************

 DBG_ENTER("COLL")

 !@gvectors_type

! integer pointers.
 if (associated(Gsph%g2sh   ))   then
   ABI_DEALLOCATE(Gsph%g2sh)
 end if
 if (associated(Gsph%gvec   ))   then
   ABI_DEALLOCATE(Gsph%gvec)
 end if
 if (associated(Gsph%rottb  ))   then
   ABI_DEALLOCATE(Gsph%rottb)
 end if
 if (associated(Gsph%rottbm1))   then
   ABI_DEALLOCATE(Gsph%rottbm1)
 end if
 if (associated(Gsph%shlim  ))   then
   ABI_DEALLOCATE(Gsph%shlim)
 end if

! integer pointers.
 if (associated(Gsph%shlen  ))   then
   ABI_DEALLOCATE(Gsph%shlen)
 end if

! complex pointers
 if (associated(Gsph%phmGt  ))   then
   ABI_DEALLOCATE(Gsph%phmGt)
 end if
 if (associated(Gsph%phmSGt ))   then
   ABI_DEALLOCATE(Gsph%phmSGt)
 end if

 call nullify_gsphere(Gsph)

 DBG_EXIT("COLL")

end subroutine destroy_gsphere
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/nullify_gsphere
!! NAME
!! nullify_gsphere
!!
!! FUNCTION
!!  Initialize all pointers of the structure to NULL.
!!
!! INPUTS
!!   Gsph = datatype whose pointers have to be nullified
!!
!! OUTPUT
!!
!! PARENTS
!!      cchi0,cchi0q0,check_completeness,m_gsphere
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine nullify_gsphere(Gsph)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_gsphere'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(gvectors_type),intent(inout) :: Gsph

! *************************************************************************

 !@gvectors_type

! integer pointers.
 nullify(Gsph%g2sh   )
 nullify(Gsph%gvec   )
 nullify(Gsph%rottb  )
 nullify(Gsph%rottbm1)
 nullify(Gsph%shlim  )

! real pointers.
 nullify(Gsph%shlen  )

! complex pointers.
 nullify(Gsph%phmGt  )
 nullify(Gsph%phmSGt )

end subroutine nullify_gsphere
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_g_idx
!! NAME
!! gsph_g_idx
!!
!! FUNCTION
!! Return the index of G in the sphere. zero if not in the sphere
!!
!! INPUTS
!!  Gsph<gvectors_type>=Info on the G-sphere
!!  gg(3)=Reduced coordinates of the G-vector.
!!
!! NOTES
!!  The function assumes that the G-vectors are ordered with increasing lenght.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function gsph_g_idx(Gsph,gg) result(g_idx)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsph_g_idx'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(gvectors_type),intent(in) :: Gsph
 integer :: g_idx
!arrays
 integer,intent(in) :: gg(3)

!Local variables-------------------------------
!scalars
 integer :: ishbsc,igs,ige
 real(dp) :: glen
 logical :: found

! *************************************************************************

 ! * Use shells and bisect to find the star and stop index thus avoiding the storage of a table (ig1,ig2)
 glen = two_pi*SQRT(DOT_PRODUCT(gg,MATMUL(Gsph%gmet,gg)))

 ishbsc = bisect(Gsph%shlen,glen)
 if ( ANY(ishbsc==(/0,Gsph%nsh/)) ) then ! glen out of range.
   g_idx=0; RETURN
 end if

 igs = Gsph%shlim(ishbsc)
 ige = Gsph%shlim(MIN(ishbsc+2,Gsph%nsh+1))-1

 g_idx=igs-1; found=.FALSE.
 do while (.not.found .and. g_idx<ige)
   g_idx=g_idx+1
   found=(ALL(Gsph%gvec(:,g_idx)==gg(:)))
 end do
 if (.not.found) g_idx=0

end function gsph_g_idx
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_gmg_idx
!! NAME
!! gsph_gmg_idx
!!
!! FUNCTION
!! Return the index of G1-G2 in the sphere. zero if not in the sphere
!!
!! INPUTS
!!  Gsph<gvectors_type>=Info on the G-sphere
!!  ig1,ig2 index of g1 and g2 in the G-sphere.
!!
!! NOTES
!!  The function assumes that the G-vectors are ordered with increasing lenght.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function gsph_gmg_idx(Gsph,ig1,ig2) result(ig1mg2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsph_gmg_idx'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(gvectors_type),intent(in) :: Gsph
 integer,intent(in) :: ig1,ig2
 integer :: ig1mg2

!Local variables-------------------------------
!scalars
 integer :: ishbsc,igs,ige
 real(dp) :: difflen
 logical :: found
!arrays
 integer :: g1mg2(3)

! *************************************************************************

 g1mg2 = Gsph%gvec(:,ig1)-Gsph%gvec(:,ig2)

 ! * Use shells and bisect to find the star and stop index thus avoiding the storage of a table (ig1,ig2)
 difflen = two_pi*SQRT(DOT_PRODUCT(g1mg2,MATMUL(Gsph%gmet,g1mg2)))

 ! FIXME It seems bisect is not portable, on my laptop test v5/t72 the number of skipped G-vectors is > 0
 ishbsc = bisect(Gsph%shlen,difflen)
 if ( ANY(ishbsc==(/0,Gsph%nsh/)) ) then ! difflen out of range.
   ig1mg2=0; RETURN
 end if

 igs = Gsph%shlim(ishbsc)
 ige = Gsph%shlim(MIN(ishbsc+2,Gsph%nsh+1))-1

 ig1mg2=igs-1; found=.FALSE.
 do while (.not.found .and. ig1mg2<ige)
   ig1mg2=ig1mg2+1
   found=(ALL(Gsph%gvec(:,ig1mg2)==g1mg2(:)))
 end do
 if (.not.found) ig1mg2=0

end function gsph_gmg_idx
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_gmg_fftidx
!! NAME
!! gsph_gmg_fftidx
!!
!! FUNCTION
!! Return the index of G1-G2 in the FFT mesh defined by ngfft. zero if not found.
!!
!! INPUTS
!!  Gsph<gvectors_type>=Info on the G-sphere
!!  ig1,ig2 index of g1 and g2 in the G-sphere.
!!  ngfft(18)=Info on the FFT mesh.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function gsph_gmg_fftidx(Gsph,ig1,ig2,ngfft) result(fft_idx)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsph_gmg_fftidx'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(gvectors_type),intent(in) :: Gsph
 integer,intent(in) :: ig1,ig2
 integer :: fft_idx
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: id1,id2,id3
!arrays
 integer :: g1mg2(3)

! *************************************************************************

 g1mg2(:)=Gsph%gvec(:,ig1)-Gsph%gvec(:,ig2)

 ! Make sure G1-G2 is still in the FFT mesh.
 ! MODULO wraps G1-G2 in the FFT box but the Fourier components are not periodic!
 if (ANY(g1mg2(:)>ngfft(1:3)/2) .or. ANY(g1mg2(:)<-(ngfft(1:3)-1)/2)) then
   fft_idx=0; RETURN
 end if

 id1=MODULO(g1mg2(1),ngfft(1))
 id2=MODULO(g1mg2(2),ngfft(2))
 id3=MODULO(g1mg2(3),ngfft(3))
 fft_idx= 1 + id1 + id2*ngfft(1) + id3*ngfft(1)*ngfft(2)

end function gsph_gmg_fftidx
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/findggp_
!! NAME
!! findggp_
!!
!! FUNCTION
!!  Fing the independent (G1,G2) pairs that are sufficient to reconstruct using
!!  symmetry properties the Fourier components f_q(G1,G2) of a two-point function
!!  which has the same symmetry of the crystal
!!
!! INPUTS
!!  nsym=number of symmetry operations
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  Gsphere<gvectors_type>=Info on the G-sphere
!!  qpt(3)=q-point in reciprocal coordinated
!!   %ng=number of G vectors in the f matrix
!!   %gmet(3,3)=metric tensor in reciprocal space
!!   %gvec(3,ng)=G vectors in reduced coordinates
!!
!! OUTPUT
!!  Gpairs_q<Gpairs_type>= Structure containing information on the irreducible pairs
!!   %niggp=nuber of independent (G1,G1) pairs
!!   %fp2ip(2,ng,ng)= for given T1 and T1 reports the sequential index of
!!     the indendent pair (G1,G2) such as (T1,T2)= S (G1,G2)
!!   %fptabo(ng,ng)= index of the symmetry operation S in the array symrec such as (T1,T2)= S (G1,G2)
!!   %ip2fp(2,niggp)= index of G1 and G2 in the array gvec for each niggp independent pair.
!!
!! PARENTS
!!      m_gsphere
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine findggp_(nsym,symrec,Gsphere,qpt,Gpairs_q)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'findggp_'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 type(gvectors_type),intent(in) :: Gsphere
 type(Gpairs_type),intent(inout) :: Gpairs_q
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 real(dp),intent(in) :: qpt(3)

!Local variables-------------------------------
!scalars
 integer :: igpi1,igpi2,isym,itim,ig,ig1,ig2,igpi,ng,niggp
 integer :: dummy,timrev_,iid,ic,istat
 real(dp) :: diff2,norm1,norm2,eps1,eps2
 real(dp) :: a1,a2,a3,d1,d2,d3,g12,g23,g31
 logical :: found,found_identity
 character(len=500) :: msg
!arrays
 integer :: identity(3,3),symqpt(4,2,nsym),ltg(nsym)
 integer :: gposs1(3),gposs2(3),gxx1(3),gxx2(3)
 integer, allocatable :: iperm(:)
 real(dp) :: gmet(3,3)
 real(dp),allocatable :: gnorm(:),gdiff2(:)
 integer,allocatable :: ip2fp(:,:)
 integer,pointer :: gvec(:,:)
 character(len=50),parameter :: my_name='findggp_'

!************************************************************************

 DBG_ENTER("COLL")
 !
 ! === Check the presence of identity and save its index ===
 identity(:,:)=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
 found_identity=.FALSE.
 do isym=1,nsym
   if (ALL(symrec(:,:,isym)==identity)) then
     iid=isym
     found_identity=.TRUE.; EXIT
   end if
 end do

 if (.not.found_identity) then
   write(msg,'(3a)')&
&   ' Only the inversion was found in the set of symmetries read from the KSS file ',ch10,&
&   ' Likely you are using an old version of the KSS file! '
   MSG_ERROR(msg)
 else
   write(msg,'(a,i2)')' findggp_ : found identity with index: ',iid
   call wrtout(std_out,msg,'COLL')
 end if
 !
 ! === Only at Gamma time-reversal can be included ===
 timrev_=1; if (ALL(ABS(qpt)<GW_TOLQ).and.Gsphere%timrev==2) timrev_=2

 write(msg,'(2a,3f8.5)')ch10,' Analyzing symmetries at point : ',qpt(:)
 if (timrev_==2) msg=TRIM(msg)//' (including time-reversal) '
 call wrtout(std_out,msg,'COLL')
 !
 ! === Find operations in the little group ===
 ! ltg is 1 if the symmetry  preserves q with a zero umklapp vector
 call symq3(nsym,qpt,symqpt,symrec,dummy,prtvol=0)

 ltg(:)=0
 do itim=1,timrev_
   do isym=1,nsym
     if (symqpt(4,itim,isym)==1.and.ALL(symqpt(1:3,itim,isym)==0)) ltg(isym)=1
   end do
 end do

 if (SUM(ltg)==1.and.timrev_==1) then ! In this case do not allocate anything, set niggp=ng**2 and exit ===ABI_ALLOCATE(ltg,)
  call wrtout(std_out,'findggp_: not enough symmetries to reduce the number of G-pairs','COLL')
  Gpairs_q%niggp=Gsphere%ng**2
  RETURN
 end if
 !
 ! === We can use symmetries to reduce the number of pairs ===
 ng=Gsphere%ng
 gmet =  Gsphere%gmet(:,:)
 gvec => Gsphere%gvec(1:3,1:ng)

 ABI_ALLOCATE(Gpairs_q%fp2ip,(2,ng,ng))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out-of-memory fp2ip')

 ABI_ALLOCATE(Gpairs_q%fptabo,(ng,ng))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out-of-memory fptabo')

 ABI_ALLOCATE(ip2fp,(2,ng**2))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out-of-memory ip2fp')
 !
 ! === Precalculate G norm to speed the loop over G-pairs ===
 ! Use the scalar-valued function to work around a bug in sunstudio12.
 ABI_ALLOCATE(gnorm,(ng))
 do ig=1,ng
   gnorm(ig)=normv(Gsphere%gvec(:,ig),gmet,'G')
   !gnorm(ig)=normv(DBLE(Gsphere%gvec(:,ig)),gmet,'G')
 end do
 !
 ! === Precalculate |G1-G2|^2 square to speed loop over G pairs ===
 ABI_ALLOCATE(gdiff2,(ng**2))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out-of-memory gdiff2')

 ABI_ALLOCATE(iperm,(ng**2))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out-of-memory iperm')

 ic=1
 g12=two*gmet(1,2) ; g23=two*gmet(2,3) ; g31=two*gmet(3,1)
 do ig1=1,ng
   a1=gvec(1,ig1) ; a2=gvec(2,ig1) ; a3=gvec(3,ig1)
   do ig2=1,ng
     d1=gvec(1,ig2)-a1 ; d2=gvec(2,ig2)-a2 ; d3=gvec(3,ig2)-a3
     gdiff2(ic)= d1*(gmet(1,1)*d1+g12*d2) &
&               +d2*(gmet(2,2)*d2+g23*d3) &
&               +d3*(gmet(3,3)*d3+g31*d1)
     iperm(ic)=ic; ic=ic+1
   end do
 end do

 ! === Sort differences ===
 call sort_dp(ng**2,gdiff2,iperm,tol10)
 !
 ! === Loop over all all possible pairs (G1,G2), finding the irreducible ones ===
 ! * Note that the pairs are addressed by ascending order of the norm of G1
 niggp=0 ! no. of pairs found
 do ic=1,ng**2

  ig1=(iperm(ic)-1)/ng+1
  ig2=iperm(ic)-(ig1-1)*ng
  diff2=gdiff2(ic)

  norm1=gnorm(ig1) ; eps1=tol8*norm1
  norm2=gnorm(ig2) ; eps2=tol8*norm2
  gposs1(:)=gvec(:,ig1)
  gposs2(:)=gvec(:,ig2)
  !
  ! === Check if this pair is the image through the same operation of a pair already found ===
  ! * Consider only vectors with the same length
  found=.FALSE.
  if (niggp>0) then
   ip : do igpi=niggp,1,-1

    if (diff2-gdiff2(igpi)>eps1+eps2+48*tol10) EXIT  ! This makes the algorithm scale as N**2
    igpi1=ip2fp(1,igpi) ; if (ABS(norm1-gnorm(igpi1))>eps1) CYCLE
    igpi2=ip2fp(2,igpi) ; if (ABS(norm2-gnorm(igpi2))>eps2) CYCLE

    do itim=1,timrev_
     do isym=1,nsym
      ! === Calculate IS G1 and IS G2 ===
      ! * Do this only for operations in the little group such as Sq=q
      ! TODO Calculate SG outside the loop, requires more memory but should be faster!
      !       avoid check on ltg, could make a table full==> little_group
      if (ltg(isym)==1) then
       gxx1=(3-2*itim)*MATMUL(symrec(:,:,isym),gvec(:,igpi1))
       if (ALL(ABS(gxx1-gposs1)==0)) then
        gxx2=(3-2*itim)*MATMUL(symrec(:,:,isym),gvec(:,igpi2))
        if (ALL(ABS(gxx2-gposs2)==0)) then
         found=.TRUE.
         Gpairs_q%fp2ip(1,ig1,ig2)=igpi1
         Gpairs_q%fp2ip(2,ig1,ig2)=igpi2
         Gpairs_q%fptabo(ig1,ig2)=isym*(3-2*itim) ! Minus is time-reversal is considered (Only at Gamma)
         EXIT ip
        end if
       end if
      end if
     end do !isym
    end do !itim
   end do ip
  end if

  if (.not.found) then ! Increment counter and fill tables ===
    niggp=niggp+1
    gdiff2(niggp)=diff2
    ip2fp(1,niggp)=ig1
    ip2fp(2,niggp)=ig2
    Gpairs_q%fp2ip(1,ig1,ig2)=ig1
    Gpairs_q%fp2ip(2,ig1,ig2)=ig2
    Gpairs_q%fptabo(ig1,ig2)=iid  ! This irreducible pair comes from the identity!
  end if
 end do !ic

 ABI_DEALLOCATE(gnorm)
 ABI_DEALLOCATE(gdiff2)
 ABI_DEALLOCATE(iperm)

 if (niggp>ng**2) then
   MSG_BUG('niggp>ng**2')
 end if

 write(msg,'(2a,i8,a,i8,a)')ch10,&
&  ' findggp_ : number of independent (G,G'') pairs found = ',niggp,' / ',ng**2,ch10
 call wrtout(std_out,msg,'COLL')
 !
 ! === Save final values ===
 Gpairs_q%niggp=niggp
 ABI_ALLOCATE(Gpairs_q%ip2fp,(2,niggp))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out-of-memory ip2fp')

 Gpairs_q%ip2fp=ip2fp(1:2,1:niggp)
 ABI_DEALLOCATE(ip2fp)

 call DEBUG_Gpairs__(Gpairs_q,Gsphere,nsym,symrec)

 DBG_EXIT("COLL")

end subroutine findggp_
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/init_Gpairs_type
!! NAME
!! init_Gpairs_type
!!
!! FUNCTION
!!  Main creation method for the Gpairs_type data type.
!!
!! INPUTS
!!  Gpairs_q<Gpairs_type>= Structure containing information on the irreducible pairs
!!   %niggp=nuber of independent (G1,G1) pairs
!!   %fp2ip(2,ng,ng)= for given T1 and T1 reports the sequential index of
!!     the indendent pair (G1,G2) such as (T1,T2)= S (G1,G2)
!!   %fptabo(ng,ng)= index of the symmetry operation S in the array symrec such as (T1,T2)= S (G1,G2)
!!   %ip2fp(2,niggp)= index of G1 and G2 in the array gvec for each niggp independent pair.
!!  qpt(3)=q-point in reciprocal coordinated
!!  Gsphere<gvectors_type>=Info on the G-sphere
!!   %ng=number of G vectors in the f matrix
!!   %gmet(3,3)=metric tensor in reciprocal space
!!   %gvec(3,ng)=G vectors in reduced coordinates
!!  Cryst<Crystal_structure> = Info on unit cell and its symmetries
!!     %nsym=number of symmetry operations
!!     %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!     %tnons(3,nsym)=fractional translations
!!
!! OUTPUT
!!
!! PARENTS
!!      mrgscr,screening
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine init_Gpairs_type(Gpairs_q,qpt,Gsphere,Cryst)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_Gpairs_type'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(gvectors_type),intent(in) :: Gsphere
 type(Crystal_structure),intent(in) :: Cryst
 type(Gpairs_type),intent(out) :: Gpairs_q
!arrays
 real(dp),intent(in) :: qpt(3)

!Local variables ------------------------------
!scalars
 integer :: ng,nsym
 integer,pointer :: symrec(:,:,:)
!************************************************************************

 !@gpairs_type
 call destroy_Gpairs_type(Gpairs_q)

 nsym   =  Cryst%nsym
 symrec => Cryst%symrec

 ng=Gsphere%ng
 !
 ! === Dimensions ===
 Gpairs_q%ng            = ng
 Gpairs_q%ngpairs       = ng**2
 Gpairs_q%nsym          = Gsphere%nsym
 Gpairs_q%timrev        = Gsphere%timrev
 Gpairs_q%can_use_timrev= .FALSE.

 ! === The point under investigation ===
 Gpairs_q%qpt=qpt
 if (ALL(ABS(qpt)<GW_TOLQ).and.Gsphere%timrev==2) Gpairs_q%can_use_timrev=.TRUE.

 call findggp_(nsym,symrec,Gsphere,qpt,Gpairs_q)

end subroutine init_Gpairs_type
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/destroy_Gpairs_type
!! NAME
!! destroy_Gpairs_type
!!
!! FUNCTION
!!  Deallocate all associated pointers defined in the Gpairs_type structure.
!!
!! INPUTS
!!  Gpairs_q<Gpairs_type>= Structure containing information on the irreducible pairs
!!   %niggp=nuber of independent (G1,G1) pairs
!!   %fp2ip(2,ng,ng)= for given T1 and T1 reports the sequential index of
!!     the indendent pair (G1,G2) such as (T1,T2)= S (G1,G2)
!!   %fptabo(ng,ng)= index of the symmetry operation S in the array symrec such as (T1,T2)= S (G1,G2)
!!   %ip2fp(2,niggp)= index of G1 and G2 in the array gvec for each niggp independent pair.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gsphere,mrgscr,screening
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine destroy_Gpairs_type(Gpairs_q)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_Gpairs_type'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Gpairs_type),intent(inout) :: Gpairs_q

!************************************************************************

 DBG_ENTER("COLL")

 !@gpairs_type
 if (associated(Gpairs_q%fp2ip ))  then
   ABI_DEALLOCATE(Gpairs_q%fp2ip)
 end if
 if (associated(Gpairs_q%fptabo))  then
   ABI_DEALLOCATE(Gpairs_q%fptabo)
 end if
 if (associated(Gpairs_q%ip2fp ))  then
   ABI_DEALLOCATE(Gpairs_q%ip2fp)
 end if

 DBG_EXIT("COLL")

end subroutine destroy_Gpairs_type
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/nullify_Gpairs_type
!! NAME
!! nullify_Gpairs_type
!!
!! FUNCTION
!!  Initialize all pointers to NULL.
!!
!! INPUTS
!!  Gpairs_q<Gpairs_type>= Structure containing information on the irreducible pairs
!!   %niggp=nuber of independent (G1,G1) pairs
!!   %fp2ip(2,ng,ng)= for given T1 and T1 reports the sequential index of
!!     the indendent pair (G1,G2) such as (T1,T2)= S (G1,G2)
!!   %fptabo(ng,ng)= index of the symmetry operation S in the array symrec such as (T1,T2)= S (G1,G2)
!!   %ip2fp(2,niggp)= index of G1 and G2 in the array gvec for each niggp independent pair.
!!
!! OUTPUT
!!
!! PARENTS
!!      mrgscr,screening
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine nullify_Gpairs_type(Gpairs_q)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_Gpairs_type'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Gpairs_type),intent(inout) :: Gpairs_q

!************************************************************************

 !@gpairs_type
 nullify(Gpairs_q%fp2ip )
 nullify(Gpairs_q%fptabo)
 nullify(Gpairs_q%ip2fp )

end subroutine nullify_Gpairs_type
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/DEBUG_Gpairs__
!! NAME
!! DEBUG_Gpairs__
!!
!! FUNCTION
!!  Mainly used for debugging purpose
!!
!! INPUTS
!!  Gpairs_q<Gpairs_type>= Structure containing information on the irreducible pairs
!!   %niggp=nuber of independent (G1,G1) pairs
!!   %fp2ip(2,ng,ng)= for given T1 and T1 reports the sequential index of
!!     the indendent pair (G1,G2) such as (T1,T2)= S (G1,G2)
!!   %fptabo(ng,ng)= index of the symmetry operation S in the array symrec such as (T1,T2)= S (G1,G2)
!!   %ip2fp(2,niggp)= index of G1 and G2 in the array gvec for each niggp independent pair.
!!  Gsphere<gvectors_type>=Info on the G-sphere
!!   %ng=number of G vectors in the f matrix
!!   %gmet(3,3)=metric tensor in reciprocal space
!!   %gvec(3,ng)=G vectors in reduced coordinates
!!  nsym=number of symmetry operations
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gsphere
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine DEBUG_Gpairs__(Gpairs_q,Gsphere,nsym,symrec)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'DEBUG_Gpairs__'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 type(Gpairs_type),intent(in) :: Gpairs_q
 type(gvectors_type),intent(in) :: Gsphere
!arrays
 integer,intent(in) :: symrec(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: ig1,ig2,ng,ir1,ir2,isym,itim
 character(len=500) :: msg
!arrays
 integer :: girr1(3),girr2(3),gxx1(3),gxx2(3),gtest1(3),gtest2(3)
 integer,pointer :: gvec(:,:)

!************************************************************************

 ng   =  Gpairs_q%ng
 gvec => Gsphere%gvec

 do ig1=1,ng
   gtest1=gvec(:,ig1)
   do ig2=1,ng
     gtest2=gvec(:,ig2)
     ir1=Gpairs_q%fp2ip(1,ig1,ig2)
     ir2=Gpairs_q%fp2ip(2,ig1,ig2)
     girr1=gvec(:,ir1)
     girr2=gvec(:,ir2)
     isym=Gpairs_q%fptabo(ig1,ig2) ; itim=1
     if (isym<0) then
       isym=-isym
       itim=2
     end if
     gxx1=(3-2*itim)*MATMUL(symrec(:,:,isym),girr1)
     gxx2=(3-2*itim)*MATMUL(symrec(:,:,isym),girr2)
     if (ANY((gtest1-gxx1)/=0).or.ANY((gtest2-gxx2)/=0)) then
       write(msg,'(a,2i8,a)')' G1-G2 pair',ig1,ig2,' has no corresponding irreducible pair '
       write(std_out,*)' G1, G2 ',gtest1,gtest2
       write(std_out,*)' independent pair? ',girr1,girr2
       write(std_out,*)' operation ',isym,symrec(:,:,isym),' with itim ',itim
       MSG_BUG(msg)
     end if
   end do
 end do

end subroutine DEBUG_Gpairs__
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/prune_g1mg2
!! NAME
!! prune_g1mg2
!!
!! FUNCTION
!! Given a list of G-vectors, evalute any possible difference G1-G2
!! remove duplicated differences and report the list of inequivalent G-vectors.
!!
!! INPUTS
!!  npw=Number of plane waves
!!  gvec(3,npw)= the reciprocal lattice vectors of the PW
!!
!! OUTPUT
!!  ngdiff=Number of inequivalent differences G1-G2
!!  g1mg2(3,ngdiff)=The set of inequivalent G1-G2 vectors.
!!
!! TODO
!!  Loop by shells to have better scaling.
!!
!! PARENTS
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine prune_g1mg2(npw,gvec,ngdiff,g1mg2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prune_g1mg2'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw
 integer,intent(out) :: ngdiff
!arrays
 integer,intent(in) :: gvec(3,npw)
 integer,pointer :: g1mg2(:,:)

!Local variables ------------------------------
!scalars
 integer :: ii,ig1,ig2
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gdiff(3)
 integer,allocatable :: g1mg2_tmp(:,:)

!************************************************************************

 ABI_ALLOCATE(g1mg2_tmp,(3,9*npw))

 ngdiff=0
 do ig2=1,npw
   do ig1=1,npw
    gdiff = gvec(:,ig1) - gvec(:,ig2)

    found=.FALSE. ; ii=0
    do while (.not.found .and. ii<ngdiff)
      ii = ii+1
      found = ALL(gdiff == g1mg2_tmp(:,ii) )
    end do

    if (.not.found) then
      ngdiff = ngdiff + 1
      if (ngdiff > 9*npw) GOTO 100
      g1mg2_tmp(:,ngdiff) = gdiff
    end if

   end do
 end do

 ! * Save results
 ABI_ALLOCATE(g1mg2,(3,ngdiff))
 g1mg2 = g1mg2_tmp(:,1:ngdiff)
 ABI_DEALLOCATE(g1mg2_tmp)

 RETURN

100 continue
 write(msg,'(2(a,i6))')' ngdiff = ',ngdiff,' > 9*npw = ',9*npw
 MSG_BUG(msg)

end subroutine prune_g1mg2
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/merge_and_sort_kg
!! NAME
!!  merge_and_sort_kg
!!
!! FUNCTION
!!  This routine merges a set of k-centered G-spheres of cutoff energy ecut and
!!  returns a Gamma-centered G-spheres. The elements in the final G-spheres are packed with increasing module.
!!
!! INPUTS
!!  nkpt=Number of k-points
!!  kptns(3,nkpt)=The k-points in reduced coordinates defining the k-centered G-spheres.
!!  ecut=Cutoff energy for the k-centered G-spheres.
!!  nsym2=Number of symmetry operations.
!!  pinv=-1 if time-reversal can be used, 1 otherwise
!!  symrel2(3,3,nsym2)=symmetry operations in real space.
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  prtvol=Flag defining the verbosity level.
!!
!! SIDE EFFECTS
!!  gbig(:,:)
!!    in input : pointer to NULL
!!    in output: gbig(3,1:npw) contains the set of G-vectors ordered by shell obtained by
!!               merging the k-centered sphere.
!!  shlim_p(:)
!!    in input : pointer to NULL
!!    in output: shlim_p(nbase)=Cumulative number of G-vectors for each shell.
!!               where nbase is the number of irreducible G"s found.
!!
!! PARENTS
!!      m_gsphere,outkss,setup_screening,setup_sigma
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine merge_and_sort_kg(nkpt,kptns,ecut,nsym2,pinv,symrel2,gprimd,gbig,prtvol,shlim_p)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'merge_and_sort_kg'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsym2,pinv,prtvol
 real(dp),intent(in) :: ecut
!arrays
 integer,intent(in) :: symrel2(3,3,nsym2)
 real(dp),intent(in) :: kptns(3,nkpt),gprimd(3,3)
 integer,pointer :: gbig(:,:)
 integer,optional,pointer :: shlim_p(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: mkmem_=1
 integer :: ikg,ig,ikpt,nbase,sizepw,in,maxpw,is,iinv,ish,ilim,mpw
 integer :: exchn2n3d,istwf_k,onpw_k,ierr,npw_k,ii,isym
 logical :: found
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: gcur(3),geq(3)
 integer :: symrec2t(3,3,nsym2)
 integer :: dum_kg(3,0)
 integer,allocatable :: gbase(:,:),gbasek(:,:,:)
 integer,allocatable :: gcurr(:,:),gshell(:,:),insort(:),gtmp(:,:)
 integer,allocatable :: nbasek(:),nshell(:),shlim(:)
 integer,allocatable :: npwarr(:)
 real(dp) :: kpoint(3),gmet(3,3)
 real(dp),allocatable :: cnorm(:),cnormk(:,:),ctmp(:)

! *********************************************************************

 ! * Fake MPI_type for the sequential part.
 ! This routine should not be parallelized as communicating gbig and other
 ! tables takes more time than recalculating them in sequential.
 call initmpi_seq(MPI_enreg_seq)

!Compute reciprocal space metrics
 do ii=1,3
   gmet(ii,:)=gprimd(1,ii)*gprimd(1,:)+&
&   gprimd(2,ii)*gprimd(2,:)+&
&   gprimd(3,ii)*gprimd(3,:)
 end do

!* Here we use TRANSPOSE(symrel2) instead of the more intuitive symrel2^{-1t} for historical reasons
!It does not affect the results since in the code below we only check the module of G
 do isym=1,nsym2
   symrec2t(:,:,isym)=TRANSPOSE(symrel2(:,:,isym))
 end do
 !
 ! ==============================================
 ! ==== Find irreducible G-vectors at each k ====
 ! ==============================================

 ABI_ALLOCATE(npwarr,(nkpt))
 exchn2n3d=0; ikg=0
 do ikpt=1,nkpt
   kpoint=kptns(:,ikpt); istwf_k=1
   call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,dum_kg,kpoint,0,MPI_enreg_seq,0,npwarr(ikpt))
 end do
 mpw = MAXVAL(npwarr)

 ABI_ALLOCATE(nbasek,(nkpt))
 ABI_ALLOCATE(gbasek,(3,mpw,nkpt))
 ABI_ALLOCATE(cnormk,(mpw,nkpt))
 nbasek=0     ! # of irreducible G at each k.
 cnormk=zero  ! Norm of each irreducible G.
 gbasek=0     ! The set of irreducible G"s at each k.

 do ikpt=1,nkpt

   kpoint = kptns(:,ikpt)
   npw_k  = npwarr(ikpt)

   exchn2n3d=0; ikg=0; istwf_k=1
   ABI_ALLOCATE(gcurr,(3,npw_k))
   call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,gcurr,kpoint,mkmem_,MPI_enreg_seq,npw_k,onpw_k)
   if (ANY(gcurr(:,1)/=0)) stop 'bug gcurr in outkss'
   !
   ! * Search for the G"s generating the others by symmetry.
   !  NB: Here we use symrec2t=TRANSPOSE(symrel2) for historical reasons, see note above
   call get_irredg(npw_k,nsym2,pinv,gprimd,symrec2t,gcurr,nbasek(ikpt),gbasek(:,:,ikpt),cnormk(:,ikpt))

   ABI_DEALLOCATE(gcurr)
 end do
 !
 ! === Reduce info over k-points ===
 ! * Here symrec2t=TRANSPOSE(symrel2) for historical reasons, see note above
 sizepw=2*mpw
 ABI_ALLOCATE(gbase,(3,sizepw))
 ABI_ALLOCATE(cnorm,(sizepw))
 nbase=0                    ! # of irred G found.

 call merge_kgirr(nsym2,pinv,nkpt,mpw,sizepw,symrec2t,nbasek,cnormk,gbasek,nbase,gbase,cnorm,ierr)
 if (ierr/=0) then
   MSG_ERROR(' merge_kgirr returned a non-zero status error')
 end if

 ABI_DEALLOCATE(nbasek)
 ABI_DEALLOCATE(cnormk)
 ABI_DEALLOCATE(gbasek)
 !
 !=== Reorder base G-vectors in order of increasing module ===
 !
 !Generate all shells of G-vectors: star of a g==set of all symetrics of this g
 ABI_ALLOCATE(gshell,(3,2*nsym2))
 ABI_ALLOCATE(shlim,(nbase))

 ABI_ALLOCATE(gbig,(3,sizepw))
!
!TODO
#if 0
!* Here symrec2t=TRANSPOSE(symrel2) for historical reasons, see note above
 call getfullg(nbase,nsym2,pinv,sizepw,gbase,symrec2t,cnorm,maxpw,gbig,shlim,ierr)
 if (ierr/0) RETURN

#else
 ABI_ALLOCATE(insort,(nbase))
 ABI_ALLOCATE(nshell,(nbase))
 do in=1,nbase
   insort(in)=in
 end do
 call sort_dp(nbase,cnorm,insort,tol14)
!
!Loop over all different modules of g''s (=shells):
 maxpw=0
 do in=1,nbase
   nshell(in)=0
   gcur(:)=gbase(:,insort(in))

   do is=1,nsym2 !  Loop over all symetries:
     do iinv=pinv,1,2
       geq(:)=iinv*(symrel2(1,:,is)*gcur(1)+symrel2(2,:,is)*gcur(2)+symrel2(3,:,is)*gcur(3))

       found=.FALSE.; ish=1
       do while ((.not.found) .and. (ish<=nshell(in))) ! Search for symetric of g and eventually add it:
         found=ALL(geq(:)==gshell(:,ish))
         ish=ish+1
       end do
       if (.not.found) then
         nshell(in)=nshell(in)+1
         gshell(:,nshell(in))=geq(:)
       end if
     end do
   end do

   if ((maxpw+nshell(in)) > sizepw) then
     ! We need to increase the size of the gbase, gbig and cnorm arrays while still keeping their content.
     ! This is done using two temporary arrays gtmp and ctmp
     MSG_WARNING("Had to reallocate gbase, gbig, cnorm")
     ABI_ALLOCATE(ctmp,(sizepw))
     ABI_ALLOCATE(gtmp,(3,sizepw))
     sizepw=maxpw+nshell(in)

     ctmp(:)=cnorm(:)
     gtmp(:,:)=gbase(:,:)

     ABI_DEALLOCATE(cnorm)
     ABI_ALLOCATE(cnorm,(sizepw))
     cnorm(:)=ctmp(:)
     ABI_DEALLOCATE(ctmp)

!    MG why this? gbase should not be changed!
     ABI_DEALLOCATE(gbase)
     ABI_ALLOCATE(gbase,(3,sizepw))
     gbase(:,:)=gtmp(:,:)
     gtmp(:,:)=gbig(:,:)

     ABI_DEALLOCATE(gbig)
     ABI_ALLOCATE(gbig,(3,sizepw))
     gbig(:,:)=gtmp(:,:)
     ABI_DEALLOCATE(gtmp)
   end if
   !
   ! Store this shell of g''s in a big array of g (gbig):
   do ig=1,nshell(in)
     gbig(:,ig+maxpw)=gshell(:,ig)
   end do
   maxpw=maxpw+nshell(in)

 end do ! End loop over shells
 !
 ! * Compute shell limits
 ilim=0
 do in=1,nbase
   ilim=ilim+nshell(in)
   shlim(in)=ilim
 end do

 if (PRESENT(shlim_p)) then ! Return shlim_p
  ABI_ALLOCATE(shlim_p,(nbase))
  shlim_p = shlim
 end if

 ! Re-allocate gbig with correct sizes so that caller can inquire the size
 ABI_ALLOCATE(gtmp,(3,ilim))
 gtmp = gbig(:,1:ilim)
 ABI_DEALLOCATE(gbig)
 ABI_ALLOCATE(gbig,(3,ilim))
 gbig=gtmp
 ABI_DEALLOCATE(gtmp)

 if (prtvol>10) then ! Print out shell limits
   write(msg,'(3a)')&
&    ' Shells found:',ch10,&
&    ' number of shell    number of G vectors      cut-off energy [Ha} '
   call wrtout(std_out,msg,'COLL')
   do in=1,nbase
     write(msg,'(12x,i4,17x,i6,12x,f8.3)')in,shlim(in),2*pi**2*cnorm(in)
     call wrtout(std_out,msg,'COLL')
   end do
   call wrtout(std_out,ch10,'COLL')
 end if

 ABI_DEALLOCATE(gshell)
 ABI_DEALLOCATE(insort)
 ABI_DEALLOCATE(nshell)
#endif

 ABI_DEALLOCATE(gbase)
 ABI_DEALLOCATE(cnorm)
 ABI_DEALLOCATE(npwarr)

end subroutine merge_and_sort_kg
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/table_gbig2kg
!! NAME
!!  table_gbig2kg
!!
!! FUNCTION
!!  Associate the kg_k set of g-vectors with the big array of gbig
!!  The array gbig(3,maxpw) contains all g-vectors used for all k-points, in order of
!!  increasing shells. For a each k-point, the wave-functions are defined only on a particular set
!!  of g-vectors kg_k (included in gbig). This set is defined by array gamma2k:
!!  The array gamma2k(ig=1,maxpw) translates the index of the gbig (from 1 to maxpw) into the corresponding
!!  index in array kg_k. If gbig(ig) does not exist in kg_k, gamma2k(ig) contains npw_k+1.
!!
!! INPUTS
!!  npw_k=Number of planewaves in the k-centered basis set
!!  kg_k(3,npw_k)=The k-centered basis set
!!  maxpw=Number of G in gbig
!!  gbig(3,maxpw)=The union of the G-spheres at different k-points.
!!
!! OUTPUT
!!  ierr=Status error. It gives the number of G of kg_k not contained in gbig.
!!  gamma2k(maxpw)=Mapping gbig -> kg_k
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine table_gbig2kg(npw_k,kg_k,maxpw,gbig,gamma2k,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'table_gbig2kg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,maxpw
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 integer,intent(in) :: gbig(3,maxpw)
 integer,intent(out) :: gamma2k(maxpw)

!Local variables-------------------------------
!scalars
 integer :: ig,igp
 logical :: found
!arrays
 integer :: gcur(3)

! *********************************************************************

 ierr=0
 gamma2k(:)=npw_k+1  ! Initialize array gamma2k

 do ig=1,npw_k       ! Loop over g-vectors, for this k point.
   gcur(:)=kg_k(:,ig)
   igp=0; found=.FALSE.
   do while ((.not.found) .and. igp<maxpw) ! Search selected vector in array gbig: TODO this part can be optimized
     igp=igp+1
     found=ALL(gcur(:)==gbig(:,igp))
   end do
   if (found) then ! Store it if found:
     gamma2k(igp)=ig
   else
     ierr=ierr+1
   end if
 end do

end subroutine table_gbig2kg
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/get_kg
!! NAME
!!  get_kg
!!
!! FUNCTION
!!  Helper function to calculate the set of G-vectors at a given kpoint.
!!  without taking advantage of FFT parallelism and G-vector distributions.
!!
!! INPUTS
!!  kpoint(3)=The k-point in reduced coordinates.
!!  ecut=Cutoff energy for planewave basis set.
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!  istwfk=Options defining if time-reversal is used to decrease the number of G"s.
!!
!! OUTPUT
!!  npw_k=Total number of G-vectors in the full G-sphere.
!!
!! SIDE EFFECTS
!!  kg_k(:,:):
!!   input : NULL pointer.
!!   output: kg_k(3,npw_k) contains the list of G-vectors.
!!
!! PARENTS
!!      fftprof,kss2wfk,m_fft_prof,m_io_kss,m_shirley,m_wfs,outkss
!!
!! CHILDREN
!!      initmpi_seq,kpgsph
!!
!! SOURCE

subroutine get_kg(kpoint,istwf_k,ecut,gmet,npw_k,kg_k)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_kg'
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k
 integer,intent(out) :: npw_k
 real(dp),intent(in) :: ecut
!arrays
 integer,pointer :: kg_k(:,:)
 real(dp),intent(in) :: gmet(3,3),kpoint(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: mkmem_=1
 integer :: exchn2n3d,ikg,npw_k_test
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: kg_dum(3,0)

! *********************************************************************

 call initmpi_seq(MPI_enreg_seq)
 !
 ! * Calculate the number of G-vectors for this k-point.
 exchn2n3d=0; ikg=0
 call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,kg_dum,kpoint,0,MPI_enreg_seq,0,npw_k)
 !
 ! * Allocate and calculate the set of G-vectors.
 ABI_ALLOCATE(kg_k,(3,npw_k))
 call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,kg_k,kpoint,mkmem_,MPI_enreg_seq,npw_k,npw_k_test)

end subroutine get_kg
!!***

!----------------------------------------------------------------------

END MODULE m_gsphere
!!***
