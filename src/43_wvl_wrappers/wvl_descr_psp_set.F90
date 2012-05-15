!!****f* defs_wvltypes/wvl_descr_psp_set
!!
!! NAME
!! wvl_descr_psp_set
!!
!! FUNCTION
!! Defines the part of the wvl%atoms%-datastructure which
!! depends on psps
!!
!! INPUTS
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!! wvl <type(wvl_internal_type)> = wavelet type
!!                 | psppar   = The covalence radii for each pseudo 
!!                 | pspcod   = the format -or code- of psp generation
!!                 | iasctype = semicore code (see defs_datatype)
!!                 | nzatom   = charge of the nucleus
!!                 | nelpsp   = the ionic pseudo-charge
!!                 | natsc    = number of atoms with semicore 
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      atomic_occupation_numbers,eleconf
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_descr_psp_set(nsppol, psps, wvl)

 use m_profiling

  use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only: eleconf, atomic_occupation_numbers, UNINITIALIZED
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_descr_psp_set'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in)                    :: nsppol
  type(wvl_internal_type), intent(inout) :: wvl
  type(pseudopotential_type), intent(in) :: psps
!arrays

!Local variables-------------------------------
!scalars
#if defined HAVE_DFT_BIGDFT
  integer, parameter :: nelecmax=32,nmax=6,lmax=4
  integer :: iat, ityp, mxpl, mxchg, nsccode
  real(dp) :: rcov,rprb,ehomo
  real(kind=8), dimension(nmax,0:lmax-1) :: neleconf
  character(len=2) :: symbol
#endif

! *********************************************************************

#if defined HAVE_DFT_BIGDFT
!We create the atoms_data structure, the part that is dependent from psp.
 wvl%atoms%psppar   => psps%gth_params%psppar
 wvl%atoms%npspcode => psps%pspcod
 wvl%atoms%ixcpsp   => psps%pspxc
 ABI_ALLOCATE(wvl%atoms%iasctype,(wvl%atoms%nat))
 ABI_ALLOCATE(wvl%atoms%nzatom,(wvl%atoms%ntypes))
 wvl%atoms%nzatom   =  psps%znucltypat
 ABI_ALLOCATE(wvl%atoms%nelpsp,(wvl%atoms%ntypes))
 wvl%atoms%nelpsp   =  psps%ziontypat
 wvl%atoms%natsc    =  0
 do iat = 1, wvl%atoms%nat, 1
    wvl%atoms%iasctype(iat) = psps%gth_params%semicore(wvl%atoms%iatype(iat))
    if (wvl%atoms%iasctype(iat) > 0) then
       wvl%atoms%natsc = wvl%atoms%natsc + 1
    end if
 end do
 ABI_ALLOCATE(wvl%atoms%aocc,(nelecmax,wvl%atoms%nat))
 ABI_ALLOCATE(wvl%atoms%amu,(wvl%atoms%ntypes))
 ABI_ALLOCATE(wvl%atoms%nlcc_ngc,(wvl%atoms%ntypes))
 ABI_ALLOCATE(wvl%atoms%nlcc_ngv,(wvl%atoms%ntypes))
 ABI_ALLOCATE(wvl%atoms%nlccpar,(0:4,1))
 do ityp = 1, wvl%atoms%ntypes, 1
    call eleconf(wvl%atoms%nzatom(ityp), wvl%atoms%nelpsp(ityp), &
         & symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,wvl%atoms%amu(ityp))
    write(wvl%atoms%atomnames(ityp), "(A)") symbol
    call atomic_occupation_numbers("input.occup",ityp,nsppol,wvl%atoms, &
         & nmax,lmax,nelecmax,neleconf,nsccode,mxpl,mxchg)
    wvl%atoms%nlcc_ngv(ityp) = UNINITIALIZED(1.0_dp)
    wvl%atoms%nlcc_ngc(ityp) = UNINITIALIZED(1.0_dp)
 end do
 ! Missing currently read radii_cf.
 wvl%atoms%donlcc = .false.
 ! TODO: add symmetry support
 wvl%atoms%symObj = -1
#endif  
end subroutine wvl_descr_psp_set
!!***
