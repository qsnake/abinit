!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gwdefs
!! NAME
!! m_gwdefs
!!
!! FUNCTION
!! This module contains definitions for a number of named constants used in the GW part of abinit
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG, FB, GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_gwdefs

 use m_profiling

 use defs_basis
 use m_errors

 use m_fstrings,      only : tolower

 implicit none

 private

! Unit number for formatted files produced by GW calculations.
! These files are not supposed to be read by abinit therefore
! their names and unit numbers are not defined in the Dtfil% structure.

 !logical,public,parameter :: gw_uses_wfk_file=.TRUE.
 logical,public,parameter :: gw_uses_wfk_file=.FALSE.

 integer,public,parameter :: unt_gw  = 21  ! GW corrections
 integer,public,parameter :: unt_sig = 22  ! Self-energy as a function of frequency
 integer,public,parameter :: unt_sgr = 23  ! Derivative wrt omega of the Self-energy
 integer,public,parameter :: unt_sgm = 20  ! Sigma on the Matsubara axis

 real(dp),public,parameter :: GW_TOLQ =0.0001_dp  
 ! Tolerance below which two BZ points are considered equal within a RL vector:
 ! for each red. direct. the abs value of the difference btw the two coord must be smaller that tolq.

 real(dp),public,parameter :: GW_TOLQ0=0.001_dp   
 ! Tolerance below which a q-point is treated as zero (long wavelength limit)

 real(dp),public,parameter :: GW_TOL_DOCC=0.01_dp 
 ! Tolerance on the difference between two occupation numbers.
 ! below this value, the contribution of the transition is neglected in the evaluation of chi0

 real(dp),public,parameter :: GW_TOL_W0=0.001_dp 
 ! Tolerance on the real part of the frequency appearing in the denominator of the 
 ! non-interacting Green function G0. Above this value, a small purely imaginary 
 ! complex shift is added to the denominator during the evaluation of chi0.

 complex(gwpc),public,parameter :: czero_gw=(0._gwp,0._gwp)
 complex(gwpc),public,parameter :: cone_gw =(1._gwp,0._gwp)
 complex(gwpc),public,parameter :: j_gw    =(0._gwp,1._gwp)

!arrays
 real(dp),public,parameter :: GW_Q0_DEFAULT(3) = (/0.00001_dp, 0.00002_dp, 0.00003_dp/)

! Flags for self-consistent GW calculations used in gw_driver and for parsing the input file.
 integer,public,parameter :: GWSC_one_shot      =1
 integer,public,parameter :: GWSC_only_W        =2
 integer,public,parameter :: GWSC_only_G        =3
 integer,public,parameter :: GWSC_both_G_and_W  =4


! Flags defining the approximation used for the self-energy (used in csigme).
 integer,public,parameter :: SIG_GW_PPM      =0  ! standard GW with PPM
 integer,public,parameter :: SIG_GW_AC       =1  ! standard GW without PPM (analytical continuation)
 integer,public,parameter :: SIG_GW_CD       =2  ! standard GW without PPM (contour deformation)
 integer,public,parameter :: SIG_HF          =5  ! Hartree-Fock calculation
 integer,public,parameter :: SIG_SEX         =6  ! Screened Exchange calculation
 integer,public,parameter :: SIG_COHSEX      =7  ! COHSEX calculation
 integer,public,parameter :: SIG_QPGW_PPM    =8  ! model GW with PPM
 integer,public,parameter :: SIG_QPGW_CD     =9  ! model GW without PPM

 public :: sigma_type_from_key
 public :: sigma_is_herm
 public :: sigma_needs_w
 public :: sigma_needs_ppm
 !public :: sigma_sc_on_wfs
 !public :: sigma_sc_on_ene
 public :: g0g0w

! Private variables
 integer,private,parameter :: STR_LEN=500
!!***

!----------------------------------------------------------------------

!!****t* m_gwdefs/epsilonm1_parameters
!! NAME
!! epsilonm1_parameters
!!
!! FUNCTION
!! For the GW part of ABINIT, the  epsilonm1_parameters structured datatype
!! gather different parameters used to calculate the inverse dielectric matrices.
!!
!! SOURCE

 type,public :: epsilonm1_parameters

!scalars
  integer :: awtr                   ! If 1 the Adler-Wiser expression for Chi_0 is evaluated
                                    !  taking advantage of time-reversal symmetry
  integer :: gwcalctyp              ! Calculation type (see input variable)
  integer :: gwcomp                 ! 1 if extrapolar technique is used. 0 otherwise.
  integer :: inclvkb                ! Integer flag related to the evaluation of the commutator for q-->0
  integer :: spmeth                 ! Method used to approximate the delta function in the expression for Im Chi_0
  integer :: nI                     ! Number of components (rows) in the chi0 matrix.
  integer :: nJ                     ! Number of components (columns) in the chi0 matrix.
  integer :: npwvec                 ! Max between npwe and npwwfn, used to pass the dimension of arrays e.g gvec
  integer :: npwwfn                 ! Number of planewaves for wavefunctions
  integer :: npwe                   ! Number of planewaves for $\tilde \epsilon$
  integer :: npwepG0                ! Number of planewaves in the enlarged sphere G-G0, to account for umklapp G0 vectors
  integer :: nbnds                  ! Number of bands used to evaluate $\tilde \epsilon$
  integer :: nkibz                  ! Number of k-points in the IBZ
  integer :: nsppol                 ! 1 for spin unpolarized, 2 for collinear spin polarized
  integer :: nqcalc                 ! Number of q-points that are calculated (subset of qibz)
  integer :: nqibz                  ! Number of q-points in the IBZ
  integer :: nqlwl                  ! Number of directions to analyze the non analytical behavior for q-->0
  integer :: nomega                 ! Number of frequencies where evaluate $\tilde \epsilon (\omega)$
  integer :: nomegaer,nomegaei      ! Number of real and imaginary frequencies, respectively
  integer :: nomegaec               ! Number of frequencies on a grid in the complex plane
                                    !    nomegaec = nomegaei*(nomegaer-1)
  integer :: npoles                 ! Number of poles for pole-fit screening
  integer :: ncoeff                 ! Number of coefficients, npoles*3+1 for phase
  integer :: nomegasf               ! Number of frequencies used for the spectral function
  integer :: symchi                 ! 0 ==> do not use symmetries to reduce the k-points summed over in chi0
                                    ! 1 ==> take advantage of point group symmetries as well as time-reversal

  real(dp) :: gwencomp              ! Extrapolar energy used if gwcomp==1.
  real(dp) :: omegaermin            ! Minimum real frequency used in the contour deformation method
  real(dp) :: omegaermax            ! Maximum real frequency used in the contour deformation method
  real(dp) :: soenergy              ! Scissor energy used in chi0
  real(dp) :: spsmear               ! Smearing of the delta in case of spmeth==2
  real(dp) :: zcut                  ! Small imaginary shift to avoid poles in chi0

  logical :: analytic_continuation  ! if true calculate chi0 only along the imaginary axis
  logical :: contour_deformation    ! if true calculated chi0 both along the real and the imaginary axis
  logical :: plasmon_pole_model     ! if true a plasmonpole model is used (only 1 or 2 frequencies are calculated)

!arrays
  integer :: mG0(3)
  ! For each reduced direction gives the max G0 component to account for umklapp processes

  real(dp),pointer :: qcalc(:,:)   SET2NULL
  ! qcalc(3,nqcalc)
  ! q-points that are explicitely calculated (subset of qibz).

  real(dp),pointer :: qibz(:,:)    SET2NULL
  ! qibz(3,nqibz)
  ! q-points in the IBZ.

  real(dp),pointer :: qlwl(:,:)    SET2NULL
  ! qlwl(3,nqlwl)
  ! q-points used for the long-wavelength limit.

  real(dp),pointer :: omegasf(:)   SET2NULL
  ! omegasf(nomegasf)
  ! real frequencies used to calculate the imaginary part of chi0.

  complex(dpc),pointer :: omega(:) SET2NULL
  ! omega(nomegasf)
  ! real and imaginary frequencies in chi0,epsilon and epsilonm1.

 end type epsilonm1_parameters
!!***

 public :: nullify_epsilonm1_parameters
 public :: destroy_epsilonm1_parameters 

 type,public :: sigij_col_t
   integer :: size1
   integer,pointer :: bidx(:)  SET2NULL
 end type sigij_col_t

 type,public :: sigijtab_t
   type(sigij_col_t),pointer :: col(:)   SET2NULL
 end type sigijtab_t

 public :: destroy_sigijtab

!----------------------------------------------------------------------

!!****t* m_gwdefs/sigma_parameters
!! NAME
!! sigma_parameters
!!
!! FUNCTION
!! For the GW part of ABINIT, the sigma_parameters structured datatype
!! gather different parameters that characterize the calculation of the matrix
!! elements of the self-energy operator.
!!
!! SOURCE

 type,public :: sigma_parameters

  integer :: gwcalctyp                   ! Calculation type
  integer :: minbdgw,maxbdgw             ! Minimum and maximum band index (considering the spin) defining
                                         ! The set of bands where GW corrections are evaluated
  integer :: gwgamma                     ! If 1 include vertex correction (GWGamma)
  integer :: gwcomp                      ! 1 if the extrapolar technique is used.

  integer :: npwvec                      ! Max betwenn npwe and npwwfn, used to pass the dimension of arrays e.g gvec
  integer :: npwwfn                      ! No. of planewaves for wavefunctions
  integer :: npwx                        ! No. of planewaves for $\Sigma_x$
  integer :: npwc                        ! No. of planewaves for $\Sigma_c$ and W
  integer :: nbnds                       ! No. of bands summed over.
  integer :: nomegasr                    ! No. of frequencies on the real axis to evaluate the spectral function
  integer :: nomegasrd                   ! No. of frequencies on the real axis to evaluate $\Sigma(E)$
  integer :: nomegasi                    ! No. of frequencies along the imaginary axis for Sigma in case of AC
  integer :: nsig_ab                     ! No. of components in the self-energy operator (1 if nspinor==1, 4 if nspinor==2)
  integer :: nspinor                     ! No. of spinorial components.
  integer :: nsppol                      ! 1 for unpolarized, 2 for spin-polarized calculation
  integer :: nkptgw                      ! No. of k-points where GW corrections have been calculated
  integer :: ppmodel                     ! Integer defining the plasmon pole model used, 0 for None.
  integer :: symsigma                    ! 0 ==> do not use symmetries to reduce the k-points summed over in sigma
                                         ! 1 ==> take advantage of space group symmetries as well as time-reversal
  integer :: use_sigxcore                ! 1 if core contribution to sigma is estimated by using Hartree-Fock 

  real(dp) :: soenergy                   ! Scissor energy used in G0
  real(dp) :: gwencomp                   ! Extrapolar energy used if gwcomp==1.

  integer :: mG0(3)                      ! For each reduced direction gives the max G0 component
                                         ! to account for umklapp processes

  real(dp) :: deltae                     ! Energy step used to evaluate numerically the derivative of the self energy
                                         ! $\frac{\partial \Re \Sigma(E)}{\partial E_o}$
  real(dp) :: minomega_r                 ! Minimum real frequency for the evaluation of the spectral function
  real(dp) :: maxomega_r                 ! Maximum real frequency for the evaluation of the spectral function
  real(dp) :: maxomega4sd                ! Maximum displacement around the KS energy where evaluate the diagonal
                                         ! Elements of $ \Sigma(E)$
  real(dp) :: omegasimax                 ! Max omega for Sigma along the imag axis in case of analytic continuation
  real(dp) :: omegasimin                 ! min omega for Sigma along the imag axis in case of analytic continuation
  real(dp) :: zcut                       ! Value of $\delta$ used to avoid the divergences (see related input variable)

  integer,pointer :: kptgw2bz(:)   SET2NULL
  ! kptgw2bz(nkptgw)
  ! For each k-point where GW corrections are calculated, the corresponding index in the BZ.

  integer,pointer :: minbnd(:,:)  SET2NULL
  integer,pointer :: maxbnd(:,:)  SET2NULL
  ! minbnd(nkptgw,nsppol), maxbnd(nkptgw,nsppol)
  ! For each k-point at which GW corrections are calculated, the min and Max band index considered
  ! (see also input variable dtset%bdgw).

  real(dp),pointer :: kptgw(:,:)  SET2NULL
  ! kptgw(3,nkptgw)
  ! k-points for the GW corrections in reduced coordinates.

  !TODO should be removed, everything should be in Sr%

  complex(dpc),pointer :: omegasi(:)  SET2NULL
  ! omegasi(nomegasi)
  ! Frequencies along the imaginary axis used for the analytical continuation.

  complex(dpc),pointer :: omega_r(:)  SET2NULL
  ! omega_r(nomegasr)
  ! Frequencies used to evaluate the spectral function.

  type(sigijtab_t),pointer :: Sigcij_tab(:,:)  SET2NULL
  ! Sigcij_tab(nkptgw,nsppol)%col(kb)%bidx(ii)  gived the index of the left wavefunction.
  ! in the <i,kgw,s|\Sigma_c|j,kgw,s> matrix elements that has to be calculated in cisgme. 
  ! in the case of self-consistent GW on wavefunctions.

  type(sigijtab_t),pointer :: Sigxij_tab(:,:)  SET2NULL
  ! Save as Sigcij_tab but for the Hermitian \Sigma_x where only the upper triangle is needed.

 end type sigma_parameters
!!***

 public :: nullify_sigma_parameters
 public :: destroy_sigma_parameters 

CONTAINS  !==============================================================================
!!***

!!****f* m_gwdefs/nullify_epsilonm1_parameters
!! NAME
!! nullify_epsilonm1_parameters
!!
!! FUNCTION
!!  Nullify the pointers in the epsilonm1_parameters structure.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!  All the pointer defined in the data types could be nullified 
!!  using the null() functions. Unfortunately null() has been 
!!  introduced in the F95 specifications and this might lead 
!!  to portability problems.
!!
!! TODO 
!!  write other methods to write the content of the data type.
!!
!! PARENTS
!!      setup_screening
!!
!! CHILDREN
!!      destroy_sigijtab
!!
!! SOURCE

subroutine nullify_epsilonm1_parameters(Ep)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_epsilonm1_parameters'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(epsilonm1_parameters),intent(inout) :: Ep

! *************************************************************************

 !@epsilonm1_parameters

!real
 nullify(Ep%qcalc  )
 nullify(Ep%qibz   )
 nullify(Ep%qlwl   )
 nullify(Ep%omegasf)

!complex
 nullify(Ep%omega) 

end subroutine nullify_epsilonm1_parameters
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/destroy_epsilonm1_parameters
!! NAME
!! destroy_epsilonm1_parameters
!!
!! FUNCTION
!!  Free all associated pointer in the structure.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      destroy_sigijtab
!!
!! SOURCE

subroutine destroy_epsilonm1_parameters(Ep)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_epsilonm1_parameters'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_parameters),intent(inout) :: Ep

! *************************************************************************

 !@epsilonm1_parameters

!real
 if (associated(Ep%qcalc  ))   then
   ABI_DEALLOCATE(Ep%qcalc)
 end if
 if (associated(Ep%qibz   ))   then
   ABI_DEALLOCATE(Ep%qibz)
 end if
 if (associated(Ep%qlwl   ))   then
   ABI_DEALLOCATE(Ep%qlwl)
 end if
 if (associated(Ep%omegasf))   then
   ABI_DEALLOCATE(Ep%omegasf)
 end if

!complex
 if (associated(Ep%omega  ))   then
   ABI_DEALLOCATE(Ep%omega)
 end if

end subroutine destroy_epsilonm1_parameters
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/destroy_sigijtab
!! NAME
!! destroy_sigijtab
!!
!! FUNCTION
!!   deallocate all memory in a sigijtab_t datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwdefs,setup_sigma
!!
!! CHILDREN
!!      destroy_sigijtab
!!
!! SOURCE

subroutine destroy_sigijtab(Sigijtab)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_sigijtab'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(sigijtab_t),intent(inout) :: Sigijtab(:,:)

!Local variables
 integer :: ii,jj,kk,ilow,iup
! *************************************************************************

 !@sigijtab_t
  do jj=1,SIZE(Sigijtab,DIM=2)
    do ii=1,SIZE(Sigijtab,DIM=1)
                                                                       
      ilow=LBOUND(Sigijtab(ii,jj)%col,DIM=1)
      iup =UBOUND(Sigijtab(ii,jj)%col,DIM=1)
      do kk=ilow,iup
        ABI_DEALLOCATE(Sigijtab(ii,jj)%col(kk)%bidx)
      end do
      ABI_DEALLOCATE(Sigijtab(ii,jj)%col)
                                                                       
    end do
  end do

end subroutine destroy_sigijtab
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/nullify_sigma_parameters
!! NAME
!! nullify_sigma_parameters
!!
!! FUNCTION
!!  Initialize all pointers to NULL.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      setup_sigma
!!
!! CHILDREN
!!      destroy_sigijtab
!!
!! SOURCE

subroutine nullify_sigma_parameters(Sigp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_sigma_parameters'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(sigma_parameters),intent(inout) :: Sigp

! *************************************************************************

 !@sigma_parameters

!integer
 nullify(Sigp%kptgw2bz)
 nullify(Sigp%minbnd  )
 nullify(Sigp%maxbnd  )

!real
 nullify(Sigp%kptgw) 

!complex
 nullify(Sigp%omegasi)
 nullify(Sigp%omega_r)

!type 
 nullify(Sigp%Sigcij_tab)
 nullify(Sigp%Sigxij_tab)

end subroutine nullify_sigma_parameters
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/destroy_sigma_parameters
!! NAME
!! destroy_sigma_parameters
!!
!! FUNCTION
!!  Destroy all associated pointers defined in the structure.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      destroy_sigijtab
!!
!! SOURCE

subroutine destroy_sigma_parameters(Sigp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_sigma_parameters'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(sigma_parameters),intent(inout) :: Sigp

! *************************************************************************

 !@sigma_parameters

!integer
 if (associated(Sigp%kptgw2bz))   then
   ABI_DEALLOCATE(Sigp%kptgw2bz)
 end if
 if (associated(Sigp%minbnd  ))   then
   ABI_DEALLOCATE(Sigp%minbnd)
 end if
 if (associated(Sigp%maxbnd  ))   then
   ABI_DEALLOCATE(Sigp%maxbnd)
 end if

!real
 if (associated(Sigp%kptgw))   then
   ABI_DEALLOCATE(Sigp%kptgw)
 end if

!complex
 if (associated(Sigp%omegasi))  then
   ABI_DEALLOCATE(Sigp%omegasi)
 end if
 if (associated(Sigp%omega_r))  then
   ABI_DEALLOCATE(Sigp%omega_r)
 end if

!logical 

!types
 if (associated(Sigp%Sigcij_tab)) then 
   call destroy_sigijtab(Sigp%Sigcij_tab)
   ABI_DEALLOCATE(Sigp%Sigcij_tab)
 end if

 if (associated(Sigp%Sigxij_tab)) then 
   call destroy_sigijtab(Sigp%Sigxij_tab)
   ABI_DEALLOCATE(Sigp%Sigxij_tab)
 end if

end subroutine destroy_sigma_parameters
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/sigma_type_from_key
!! NAME
!! sigma_type_from_key
!!
!! FUNCTION
!!  Return a string definining the particular approximation used for the self-energy.
!!  Stops if the key is not in the list of allowed possibilities.
!!
!! INPUTS
!!  key=Integer
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function sigma_type_from_key(key) result(sigma_type)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigma_type_from_key'
!End of the abilint section

 implicit none

 integer,intent(in) :: key
 character(len=STR_LEN) :: sigma_type

!Local variables ------------------------------
!scalars
 character(len=500) :: msg

!************************************************************************

 sigma_type = "None"
 if (key==SIG_GW_PPM  )  sigma_type = ' standard GW with PPM'
 if (key==SIG_GW_AC   )  sigma_type = ' standard GW without PPM (analytical continuation)'
 if (key==SIG_GW_CD   )  sigma_type = ' standard GW without PPM (contour deformation)'
 if (key==SIG_HF      )  sigma_type = ' Hartree-Fock calculation'
 if (key==SIG_SEX     )  sigma_type = ' Screened Exchange calculation'
 if (key==SIG_COHSEX  )  sigma_type = ' COHSEX calculation'
 if (key==SIG_QPGW_PPM)  sigma_type = ' model GW with PPM'
 if (key==SIG_QPGW_CD )  sigma_type = ' model GW without PPM'

 if (sigma_type == "None") then
   write(msg,'(a,i0)')" Unknown value for key= ",key
   MSG_ERROR(msg)
 end if

end function sigma_type_from_key
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/sigma_is_herm
!! NAME
!! sigma_is_herm
!!
!! FUNCTION
!!  Return .TRUE. if the approximated self-energy is hermitian.
!!
!! INPUTS
!!  Sigp<sigma_parameters>=datatype gathering data and info on the self-energy run.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function sigma_is_herm(Sigp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigma_is_herm'
!End of the abilint section

 implicit none

!Arguments ------------------------------
!scalars
 logical :: sigma_is_herm
 type(sigma_parameters),intent(in) :: Sigp

!Local variables ------------------------------
 integer :: mod10

!************************************************************************

 mod10=MOD(Sigp%gwcalctyp,10)
 sigma_is_herm = ANY(mod10 == (/SIG_HF, SIG_SEX, SIG_COHSEX/) )

end function sigma_is_herm  
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/sigma_needs_w
!! NAME
!! sigma_needs_w
!!
!! FUNCTION
!!  Return .TRUE. if self-energy requires the screened interaction W. 
!!  For example HF does not need the SCR file.
!!
!! INPUTS
!!  Sigp<sigma_parameters>=datatype gathering data and info on the self-energy run.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function sigma_needs_w(Sigp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigma_needs_w'
!End of the abilint section

 implicit none

!Arguments ------------------------------
!scalars
 logical :: sigma_needs_w
 type(sigma_parameters),intent(in) :: Sigp

!Local variables ------------------------------
 integer :: mod10

!************************************************************************

 mod10=MOD(Sigp%gwcalctyp,10)
 sigma_needs_w = (mod10/=SIG_HF)

end function sigma_needs_w
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/sigma_needs_ppm
!! NAME
!! sigma_needs_ppm
!!
!! FUNCTION
!!  Return .TRUE. if the self-energy run requires a plasmon-pole model.
!!
!! INPUTS
!!  Sigp<sigma_parameters>=datatype gathering data and info on the self-energy run.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function sigma_needs_ppm(Sigp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigma_needs_ppm'
!End of the abilint section

 implicit none

!Arguments ------------------------------
!scalars
 logical :: sigma_needs_ppm
 type(sigma_parameters),intent(in) :: Sigp

!Local variables ------------------------------
 integer :: mod10

!************************************************************************

 mod10=MOD(Sigp%gwcalctyp,10)
 sigma_needs_ppm = ( ANY(mod10 == (/SIG_GW_PPM, SIG_QPGW_PPM/)) .or. &
&                   Sigp%gwcomp==1                                   &
&                  )

end function sigma_needs_ppm
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/g0g0w
!! NAME
!! g0g0w
!!
!! FUNCTION
!!  Calculates the frequency-dependent part of the RPA polarizability G0G0.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function g0g0w(omega,numerator,delta_ene,zcut,TOL_W0,opt_poles)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'g0g0w'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in):: opt_poles
 real(dp),intent(in) :: TOL_W0,delta_ene,numerator,zcut
 complex(dpc) :: g0g0w 
 complex(dpc),intent(in) :: omega

!Local variables ------------------------------
!scalars
 real(dp) :: sgn
 character(len=500) :: msg

!************************************************************************

 if (delta_ene**2>tol14) then
   sgn=delta_ene/ABS(delta_ene)
   !
   if (opt_poles == 2) then ! Resonant and anti-resonant contributions.
     if (DABS(REAL(omega))>TOL_W0) then
       g0g0w =  numerator / (omega + delta_ene - j_dpc*sgn*zcut)&
&              -numerator / (omega - delta_ene + j_dpc*sgn*zcut)
     else
       g0g0w =  numerator / (omega + delta_ene)&
&              -numerator / (omega - delta_ene)
     end if
   
   else if (opt_poles == 1) then ! Only resonant contribution is included.
     if (DABS(REAL(omega))>TOL_W0) then
       g0g0w =  numerator / (omega + delta_ene - j_dpc*sgn*zcut)
     else
       g0g0w =  numerator / (omega + delta_ene)
     end if

   else
     write(msg,'(a,i0)')" Wrong value for opt_poles: ",opt_poles
     MSG_ERROR(msg)
   end if ! opt_poles

 else ! delta_ene**2<tol14
   g0g0w = czero
 end if

end function g0g0w
!!***

END MODULE m_gwdefs
!!***
