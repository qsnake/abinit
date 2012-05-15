!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hamiltonian
!! NAME
!! m_hamiltonian
!!
!! FUNCTION
!!  This module provides the definition of the gs_hamiltonian_type datastructure
!!  used in the "getghc" routine to apply the Hamiltonian on a wavefunction.
!!  methods to Initialize or destroy the object are defined here.
!!  It also defines the ddiago_ctl_type structures datatype used to control the algorithm
!!  used in ks_ddiago for performing the direct diagonalization of the KS Hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!  The definition of gs_hamiltonian_type should be done here instead of defs_datatypes.
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

module m_hamiltonian

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_fstrings,           only : toupper

 implicit none

 private

 public ::  nullify_hamiltonian
 public ::  init_hamiltonian
 public ::  destroy_hamiltonian
 public ::  finalize_hamiltonian
!!***

!!****t* m_hamiltonian/ddiago_ctl_type
!! NAME
!!  ddiago_ctl_type
!!
!! FUNCTION
!!  Structure storing the variables controlling the direct diagonalization of the Kohn-Sham Hamiltonian.
!!
!! SOURCE

 type, public :: ddiago_ctl_type

  integer :: isppol
   ! The spin component of the Hamiltonian (1 if nspinor==1 or nsppol==1).

  integer :: istwf_k
   ! Option defining whether time-reversal symmetry is used at particular k-points
   ! If 0, the code will automatically use TR symmetry if possible (depending on the k-point)

  integer :: nband_k
   ! Number of bands to be calculated.

  integer :: npw_k
  ! The number of planes waves for the wavefunctions taking into account time-reversal symmetry.

  integer :: npwtot
  ! The number of planes waves in the Hamiltonian without taking into account istwf_k

  integer :: nspinor
  ! Number of spinorial components.

  integer :: prtvol
   ! Flag controlling the verbosity.

  integer :: use_scalapack
  ! 0 if diagonalization is done in sequential on each node.
  ! 1 to use scalapack

  real(dp) :: abstol
   ! used fro RANGE="V","I", and "A" when do_full_diago=.FALSE.
   ! The absolute error tolerance for the eigenvalues. An approximate eigenvalue is accepted
   ! as converged when it is determined to lie in an interval [a,b] of width less than or equal to
   !
   !         ABSTOL + EPS *   max( |a|,|b| ) ,
   !
   ! where EPS is the machine precision.  If ABSTOL is less than or equal to zero, then  EPS*|T|  will be used in its place,
   ! where |T| is the 1-norm of the tridiagonal matrix obtained by reducing A to tridiagonal form.
   !
   ! Eigenvalues will be computed most accurately when ABSTOL is
   ! set to twice the underflow threshold 2*DLAMCH('S'), not zero.
   ! If this routine returns with INFO>0, indicating that some
   ! eigenvectors did not converge, try setting ABSTOL to 2*DLAMCH('S').

  real(dp) :: ecut
   ! The cutoff energy for the plane wave basis set.

  real(dp) :: ecutsm
   ! Smearing energy for plane wave kinetic energy (Ha)

  real(dp) :: effmass
   ! Effective mass for electrons (usually one).

  logical :: do_full_diago
  ! Specifies whether direct or partial diagonalization will be performed.
  ! Meaningful only if RANGE='A'.

  integer :: ilu(2)
   ! If RANGE='I', the indices (in ascending order) of the smallest and largest eigenvalues to be returned.
   ! il=ilu(1), iu=ilu(2) where
   ! 1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0. NOT used if RANGE = 'A' or 'V'.

  integer :: nloalg(5)

  real(dp) :: kpoint(3)
   ! The k-point in reduced coordinates at which the Hamiltonian is diagonalized.

  real(dp) :: vlu(2)
   ! If RANGE='V', the lower and upper bounds of the interval to
   ! be searched for eigenvalues. vl=vlu(1) and vu=vlu(2) with VL < VU.
   ! Not referenced if RANGE = 'A' or 'I'.

  character(len=1) :: jobz
   ! character defining whether wavefunctions are required (lapack option).
   ! "N":  Compute eigenvalues only;
   ! "V":  Compute eigenvalues and eigenvectors.

  character(len=1) :: range
   ! character defining the subset of eigenstates that will be calculated (lapack option).
   ! "A": all eigenvalues will be found.
   ! "V": all eigenvalues in the half-open interval (VL,VU] will be found.
   ! "I": the IL-th through IU-th eigenvalues will be found.

  !$character(len=fnlen) :: fname
  ! The name of the file storing the eigenvectors and eigenvalues (only if jobz="V")

 end type ddiago_ctl_type

 public :: init_ddiago_ctl

CONTAINS  !===========================================================
!!***

!!****f* m_hamiltonian/nullify_hamiltonian
!! NAME
!!  nullify_hamiltonian
!!
!! FUNCTION
!!  Set all pointers in a gs_hamiltonian_type structure to NULL
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=Structure will all pointers Initialized to NULL.
!!
!! PARENTS
!!      m_hamiltonian
!!
!! CHILDREN
!!      initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine nullify_hamiltonian(Ham)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_hamiltonian'
!End of the abilint section

 type(gs_hamiltonian_type),intent(inout) :: Ham

! *************************************************************************

 !@gs_hamiltonian_type

! Integer pointers.
 nullify(Ham%atindx )
 nullify(Ham%atindx1)
 nullify(Ham%gbound )
 nullify(Ham%indlmn )
 !nullify(Ham%indpw_k)
 !nullify(Ham%kg_k  )
 nullify(Ham%nattyp )
 nullify(Ham%pspso  )
 nullify(Ham%typat  )

! Real pointers
 nullify(Ham%ekb    )
 nullify(Ham%sij    )
 !nullify(Ham%ffnl  )
 !nullify(Ham%kinpw )
 nullify(Ham%phkxred)
 nullify(Ham%ph1d   )
 !nullify(Ham%ph3d  )
 !nullify(Ham%vlocal)
 nullify(Ham%xred   )
 !nullify(Ham%ylm   )

end subroutine nullify_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/destroy_hamiltonian
!! NAME
!!  destroy_hamiltonian
!!
!! FUNCTION
!!  Clean and destroy gs_hamiltonian_type datastructure
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=All dynamic memory defined in the structure is deallocated.
!!
!! PARENTS
!!      ks_ddiago,m_shirley,nstdy3,nstpaw3,vtorho,vtorho3
!!
!! CHILDREN
!!      initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine destroy_hamiltonian(Ham)

!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_hamiltonian'
!End of the abilint section

 type(gs_hamiltonian_type),intent(inout) :: Ham

! *************************************************************************

 !@gs_hamiltonian_type

! Integer pointers.
 if (associated(Ham%atindx ))  then
   ABI_DEALLOCATE(Ham%atindx)
 end if
 if (associated(Ham%atindx1))  then
   ABI_DEALLOCATE(Ham%atindx1)
 end if
 if (associated(Ham%gbound ))  then
   ABI_DEALLOCATE(Ham%gbound)
 end if
 if (associated(Ham%indlmn ))  then
   ABI_DEALLOCATE(Ham%indlmn)
 end if
 !if (associated(indpw_k   )) deallocate(indpw_k    )
 !if (associated(kg_k      )) deallocate(kg_k       )
 if (associated(Ham%nattyp ))  then
   ABI_DEALLOCATE(Ham%nattyp)
 end if
 if (associated(Ham%pspso  ))  then
   ABI_DEALLOCATE(Ham%pspso)
 end if
 if (associated(Ham%typat  ))  then
   ABI_DEALLOCATE(Ham%typat)
 end if

! Real pointers
 if (associated(Ham%ekb    ))   then
   ABI_DEALLOCATE(Ham%ekb)
 end if
 if (associated(Ham%sij    ))   then
   ABI_DEALLOCATE(Ham%sij)
 end if
 !if (associated(Ham%ffnl  ))  deallocated(Ham%ffnl  )
 !if (associated(Ham%kinpw ))  deallocated(Ham%kinpw )
 if (associated(Ham%phkxred))   then
   ABI_DEALLOCATE(Ham%phkxred)
 end if
 if (associated(Ham%ph1d   ))   then
   ABI_DEALLOCATE(Ham%ph1d)
 end if
 !if (associated(Ham%ph3d  ))  deallocated(Ham%ph3d  )
 !if (associated(Ham%vlocal))  deallocated(Ham%vlocal)
 if (associated(Ham%xred   ))   then
   ABI_DEALLOCATE(Ham%xred)
 end if
 !if (associated(Ham%ylm   ))  deallocated(Ham%ylm   )

end subroutine destroy_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/init_hamiltonian
!! NAME
!!  init_hamiltonian
!!
!! FUNCTION
!!  Creation method for the gs_hamiltonian_type structure.
!!  It allocates memory and initializes all quantities that do not depend on the k-point or spin.
!!
!! INPUTS
!!  nfft=Number of FFT grid points (for this processors)
!!  natom=Number of atoms in the unit cell.
!!  ntypat=Number of type of atoms.
!!  nspinor=Number of spinorial components
!!  nspden=Number of spin density components.
!!  mgfft=Maximum size for 1D FFTs i.e., MAXVAL(ngfft(1:3))
!!  psps<pseudopotential_type>=structure datatype gathering data on the pseudopotentials.
!!  [electronpositron<electronpositron_type>]=Structured datatype storing data for the
!!    electron-positron two-component DFT (optional).
!!  ngfft(18)=integer array with FFT box dimensions and other information on FFTs, for the FINE rectangular grid.
!!  nloalg(5)=governs the choice of the algorithm for non-local operator
!!  typat(natom)=Type of each atom.
!!  [ph1d(2,3*(2*mgfft+1)*natom)]= 1-dimensions phase arrays for structure factor (see getph.f).
!!    Recalculated inside the routine if not present in input.
!!  rprimd(3,3)=Direct lattice vectors in Bohr.
!!  xred(3,natom)=Reduced coordinates of the atoms.
!!  paw_ij(natom*psps%usepaw)<paw_ij_type>=Various arrays given on (i,j) (partial waves) channels.
!!  pawtab(ntypat*psps%usepaw)<pawtab_type>=PAW TABulated data initialized at start.
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=Structured datatype almost completely initialized:
!!   * Basic variables and dimensions are transfered to the structure.
!!   * All pointers are allocated with correct dimensions.
!!   * Quantities that do not depend on the k-point or spin are initialized.
!!
!! PARENTS
!!      ks_ddiago,m_shirley,nstdy3,nstpaw3,vtorho,vtorho3
!!
!! CHILDREN
!!      initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine init_hamiltonian(gs_hamk,Psps,paw_ij,pawtab,nspinor,nspden,natom,ntypat,typat,xred,&
&                           nfft,mgfft,ngfft,rprimd,nloalg,ph1d,electronpositron,use_gpu_cuda)

 use defs_basis
 use m_electronpositron,   only : electronpositron_type, electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_hamiltonian'
 use interfaces_42_geometry
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,natom,ntypat,nspinor,nspden,mgfft
 integer,optional,intent(in) :: use_gpu_cuda
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(pseudopotential_type),intent(in) :: psps
 type(electronpositron_type),optional,pointer :: electronpositron
!arrays
 integer,intent(in) :: ngfft(18),nloalg(5)
 integer,intent(in) :: typat(natom)
 real(dp),optional,intent(in) :: ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 type(paw_ij_type),intent(in) :: paw_ij(natom*psps%usepaw)
 type(pawtab_type),intent(in)  :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: itypat,iat,indx,ilmn,cplex_dij
 real(dp) :: ucvol
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)

! *************************************************************************

 !@gs_hamiltonian_type
 call nullify_hamiltonian(gs_hamk)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ABI_CHECK(mgfft==MAXVAL(ngfft(1:3)),"Wrong mgfft")

!Allocate the arrays of the Hamiltonian whose dimensions do not depend on k
 ABI_ALLOCATE(gs_hamk%atindx,(natom))
 ABI_ALLOCATE(gs_hamk%atindx1,(natom))
 ABI_ALLOCATE(gs_hamk%typat,(natom))
 gs_hamk%typat=typat(1:natom)
 ABI_ALLOCATE(gs_hamk%gbound,(2*mgfft+8,2))
 gs_hamk%gbound(:,:)=0
 ABI_ALLOCATE(gs_hamk%indlmn,(6,psps%lmnmax,ntypat))
 ABI_ALLOCATE(gs_hamk%nattyp,(ntypat))
 ABI_ALLOCATE(gs_hamk%phkxred,(2,natom))
 ABI_ALLOCATE(gs_hamk%ph1d,(2,3*(2*mgfft+1)*natom))
 ABI_ALLOCATE(gs_hamk%pspso,(ntypat))
 ABI_ALLOCATE(gs_hamk%xred,(3,natom))

!Initialize most of the Hamiltonian
 indx=1
 do itypat=1,ntypat
   gs_hamk%nattyp(itypat)=0
   do iat=1,natom
     if (typat(iat)==itypat) then
       gs_hamk%atindx (iat )=indx
       gs_hamk%atindx1(indx)=iat
       indx=indx+1
       gs_hamk%nattyp(itypat)=gs_hamk%nattyp(itypat)+1
     end if
   end do
 end do

 gs_hamk%gmet(:,:)  =gmet(:,:)
 gs_hamk%gprimd(:,:)=gprimd(:,:)
 gs_hamk%indlmn(:,:,:)=psps%indlmn(:,:,:)
 gs_hamk%lmnmax     =psps%lmnmax
 gs_hamk%mgfft      =mgfft
 gs_hamk%mproj      =psps%mproj
 gs_hamk%mpsang     =psps%mpsang
 gs_hamk%mpssoang   =psps%mpssoang
 gs_hamk%natom      =natom
 gs_hamk%nfft       =nfft
 gs_hamk%ngfft(:)   =ngfft(:)
 gs_hamk%nloalg(:)  =nloalg(:)
 !$gs_hamk%matblk=nloalg(4); if (nloalg(1)>0) gs_hamk%matblk=natom
 gs_hamk%nspinor    =nspinor
 gs_hamk%ntypat     =ntypat

 gs_hamk%nvloc=1; if(nspden==4)gs_hamk%nvloc=4
 gs_hamk%n4         =ngfft(4)
 gs_hamk%n5         =ngfft(5)
 gs_hamk%n6         =ngfft(6)
 gs_hamk%usepaw     =psps%usepaw
 if (PRESENT(ph1d)) then
   gs_hamk%ph1d(:,:)  =ph1d(:,:)
 else ! Recalculate structure factor phases
   call getph(gs_hamk%atindx,natom,ngfft(1),ngfft(2),ngfft(3),gs_hamk%ph1d,xred)
 end if
 gs_hamk%pspso(:)   =psps%pspso(1:ntypat)
 gs_hamk%ucvol      =ucvol
 gs_hamk%useylm     =psps%useylm
 gs_hamk%use_gpu_cuda=0
 if(PRESENT(use_gpu_cuda)) gs_hamk%use_gpu_cuda=use_gpu_cuda
 gs_hamk%xred(:,:)  =xred(:,:)

! ===========================
! ==== Non-local factors ====
! ===========================

 if (psps%usepaw==0) then ! Norm-conserving: use constant Kleimann-Bylander energies.
   gs_hamk%dimekb1=psps%dimekb
   gs_hamk%dimekb2=ntypat
   ABI_ALLOCATE(gs_hamk%ekb,(psps%dimekb,ntypat,nspinor**2))
   ABI_ALLOCATE(gs_hamk%sij,(0,0))
   gs_hamk%ekb(:,:,1)=psps%ekb(:,:)
   if (nspinor==2) then
     gs_hamk%ekb(:,:,2)=psps%ekb(:,:)
     gs_hamk%ekb(:,:,3:4)=zero
   end if
   if (PRESENT(electronpositron)) then
     if (electronpositron_calctype(electronpositron)==1) gs_hamk%ekb(:,:,:)=-gs_hamk%ekb(:,:,:)
   end if

 else ! PAW: store overlap coefficients and allocate memory for Dij coefficients (spin dependent)
   cplex_dij=paw_ij(1)%cplex_dij
   gs_hamk%dimekb1=psps%dimekb*cplex_dij
   gs_hamk%dimekb2=natom
   ABI_ALLOCATE(gs_hamk%ekb,(gs_hamk%dimekb1,gs_hamk%dimekb2,nspinor**2))
   ABI_ALLOCATE(gs_hamk%sij,(gs_hamk%dimekb1,ntypat))
   do itypat=1,ntypat
     if (cplex_dij==1) then
       gs_hamk%sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
     else
       do ilmn=1,pawtab(itypat)%lmn2_size
         gs_hamk%sij(2*ilmn-1,itypat)=pawtab(itypat)%sij(ilmn)
         gs_hamk%sij(2*ilmn  ,itypat)=zero
       end do
     end if
   end do
 end if

end subroutine init_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/finalize_hamiltonian
!! NAME
!!  finalize_hamiltonian
!!
!! FUNCTION
!!  Setup of the k-dependent part of the Hamiltonian.
!!
!! SIDE EFFECTS
!!  Ham<gs_hamiltonian_type>=
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine finalize_hamiltonian(gs_hamk,isppol,npw_k,istwfk,kpoint,paw_ij)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'finalize_hamiltonian'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,istwfk,isppol
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(paw_ij_type),intent(in) :: paw_ij(gs_hamk%natom*gs_hamk%usepaw)
!arrays
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer :: iat,iatom,ispden,isp,nspinor,natom,dimdij,ilmn
 real(dp) :: arg

! *************************************************************************

 !@gs_hamiltonian_type

 ! Setup of the k-dependent part of the Hamiltonian.
 gs_hamk%npw      = npw_k
 gs_hamk%istwf_k  = istwfk
 gs_hamk%kpoint(:)= kpoint(:)

 natom   = gs_hamk%natom
 nspinor = gs_hamk%nspinor

!  Allocate the arrays phkxred and ph3d, compute phkxred and eventually ph3d.
 do iat=1,natom
   iatom=gs_hamk%atindx(iat)
   arg=two_pi*DOT_PRODUCT(kpoint,gs_hamk%xred(:,iat))
   gs_hamk%phkxred(1,iatom)=DCOS(arg)
   gs_hamk%phkxred(2,iatom)=DSIN(arg)
 end do

 if (gs_hamk%usepaw==1) then ! Retrieve PAW Dij coefficients for this spin component

  do ispden=1,nspinor**2
    isp=isppol; if (nspinor==2) isp=ispden
    do iatom=1,natom
      dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
      do ilmn=1,dimdij
        gs_hamk%ekb(ilmn,iatom,ispden)=paw_ij(iatom)%dij(ilmn,isp)
      end do
      if(dimdij+1<=gs_hamk%dimekb1) gs_hamk%ekb(dimdij+1:gs_hamk%dimekb1,iatom,ispden)=zero
    end do
  end do

 end if

end subroutine finalize_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/init_ddiago_ctl
!! NAME
!!  init_ddiago_ctl
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_shirley,outkss
!!
!! CHILDREN
!!      initmpi_seq,kpgsph,wrtout
!!
!! SOURCE

subroutine init_ddiago_ctl(Dctl,jobz,isppol,nspinor,ecut,kpoint,nloalg,gmet,&
& nband_k,istwf_k,ecutsm,effmass,abstol,range,ilu,vlu,use_scalapack,prtvol)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_ddiago_ctl'
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol,nspinor
 integer,optional,intent(in) :: istwf_k,prtvol,use_scalapack,nband_k
 real(dp),intent(in) :: ecut
 real(dp),optional,intent(in) :: ecutsm,effmass
 real(dp),optional,intent(in) :: abstol
 character(len=*),intent(in) :: jobz
 character(len=*),optional,intent(in) :: range
 type(ddiago_ctl_type),intent(out) :: Dctl
!arrays
 integer,intent(in) :: nloalg(5)
 integer,optional,intent(in) :: ilu(2)
 real(dp),intent(in) :: kpoint(3)
 real(dp),optional,intent(in) :: vlu(2)
 real(dp),intent(in) :: gmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: npw_k
 logical :: ltest
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer,allocatable :: kg_k(:,:)

! *************************************************************************

 call initmpi_seq(MPI_enreg_seq) ! Fake MPI_type.

 Dctl%isppol  = isppol
 Dctl%nspinor = nspinor
 Dctl%kpoint  = kpoint

 if (PRESENT(istwf_k)) then
  Dctl%istwf_k = istwf_k
 else
  Dctl%istwf_k = set_istwfk(kpoint)
 end if

 ABI_CHECK(Dctl%istwf_k==1,"istwf_k/=1 not coded")

 Dctl%jobz   = toupper(jobz(1:1))
 Dctl%range  = "A"
 if (PRESENT(range)) then
  Dctl%range = toupper(range)
 end if

 Dctl%ecut = ecut

 Dctl%ecutsm = zero; if (PRESENT(ecutsm)) Dctl%ecutsm = ecutsm

 Dctl%effmass = one; if (PRESENT(effmass)) Dctl%effmass = effmass


 Dctl%nloalg  = nloalg

 Dctl%prtvol = 0; if (PRESENT(prtvol)) Dctl%prtvol = prtvol

 Dctl%abstol = -tol8; if (PRESENT(abstol)) Dctl%abstol = abstol

 ABI_ALLOCATE(kg_k,(3,0))

! * Total number of G-vectors for this k-point with istwf_k=1.
 call kpgsph(ecut,0,gmet,0,0,1,kg_k,kpoint,0,MPI_enreg_seq,0,Dctl%npwtot)

! * G-vectors taking into account time-reversal symmetry.
 call kpgsph(ecut,0,gmet,0,0,istwf_k,kg_k,kpoint,0,MPI_enreg_seq,0,npw_k)

 Dctl%npw_k = npw_k
 ABI_DEALLOCATE(kg_k)

 Dctl%do_full_diago = .FALSE.

 SELECT CASE (Dctl%range)
  CASE ("A")

  ! Check on the number of stored bands.
  Dctl%nband_k=-1
  if (PRESENT(nband_k)) Dctl%nband_k=nband_k

  if (Dctl%nband_k==-1.or.Dctl%nband_k>=npw_k*nspinor) then
    Dctl%nband_k=npw_k*nspinor
    write(msg,'(4a)')ch10,&
&    ' Since the number of bands to be computed was (-1) or',ch10,&
&    ' too large, it has been set to the max. value npw_k*nspinor. '
    if (Dctl%prtvol>0) call wrtout(std_out,msg,'COLL')
  else
    Dctl%nband_k=nband_k
  end if

  Dctl%do_full_diago = (Dctl%nband_k==npw_k*nspinor)

  if (Dctl%do_full_diago) then
    write(msg,'(6a)')ch10,&
&    ' Since the number of bands to be computed',ch10,&
&    ' is equal to the number of G-vectors found for this k-point,',ch10,&
&    ' the program will perform complete diagonalization.'
  else
    write(msg,'(6a)')ch10,&
&     ' Since the number of bands to be computed',ch10,&
&     ' is less than the number of G-vectors found,',ch10,&
&     ' the program will perform partial diagonalization.'
  end if
  if (Dctl%prtvol>0) call wrtout(std_out,msg,'COLL')

 CASE ("I")
  if (.not.PRESENT(ilu)) then
    MSG_ERROR(" ilu must be specified when range=I ")
  end if
  Dctl%ilu = ilu

  ltest = ( ( ilu(2)>=ilu(1) ) .and. ilu(1)>=1 .and. ilu(2)<=Dctl%npwtot )
  write(msg,'(a,2i0)')" Illegal value for ilu: ",ilu
  ABI_CHECK(ltest,msg)
  Dctl%nband_k= ilu(2)-ilu(1)+1

 CASE ("V")
  if (.not.PRESENT(vlu)) then
    MSG_ERROR(" vlu must be specified when range=V ")
  end if
  Dctl%vlu = vlu

  Dctl%nband_k=-1 !??

  ltest = (vlu(2)>vlu(1))
  write(msg,'(a,2f0.3)')" Illegal value for vlu: ",vlu
  ABI_CHECK(ltest,msg)

 CASE DEFAULT
   msg = " Unknown value for range: "//TRIM(Dctl%range)
   MSG_ERROR(msg)
 END SELECT

 ! Consider the case in which we asked for the entire set of eigenvectors
 ! but the number of bands is less that npw_k. Therefore have to prepare the call to ZHEEVX.
 ! TODO this has to be done in a cleaner way.
 if (Dctl%range=="A".and. (.not.Dctl%do_full_diago)) then
   Dctl%range="I"
   Dctl%ilu(1) = 1
   Dctl%ilu(2) = npw_k*nspinor
   Dctl%nband_k= npw_k*nspinor
 end if

 Dctl%use_scalapack=0
 if (PRESENT(use_scalapack)) then
   Dctl%use_scalapack=use_scalapack
 end if
 ABI_CHECK(Dctl%use_scalapack==0," scalapack mode not coded")

end subroutine init_ddiago_ctl
!!***

!----------------------------------------------------------------------

END MODULE m_hamiltonian
!!***
