!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_screening
!! NAME
!!  m_screening
!!
!! FUNCTION
!!  This module contains the definition of the object used to deal 
!!  with the inverse dielectric matrix as well as related methods.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG)
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

MODULE m_screening

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_copy
 use m_splines
 use m_lebedev

 use m_gwdefs,         only : GW_TOLQ0, czero_gw
 use m_fstrings,       only : toupper
 use m_io_tools,       only : get_unit
 use m_numeric_tools,  only : print_arr, hermitianize, bisect
 use m_blas,           only : xdotu
 use m_geometry,       only : normv
 use m_abilasi,        only : xginv
 use m_crystal,        only : crystal_structure , destroy_crystal       
 use m_crystal_io,     only : init_crystal_from_hdr
 use m_bz_mesh,        only : bz_mesh_type, init_kmesh, get_BZ_item, box_len, has_bz_item, print_bz_mesh,&
&                             destroy_bz_mesh_type,find_qmesh
 use m_fft_mesh,       only : g2ifft
 use m_gsphere,        only : gpairs_type, gvectors_type, init_gsphere, gsph_g_idx, gsph_gmg_idx, destroy_gsphere
 use m_vcoul,          only : vcoul_t
 use m_io_screening,   only : free_scrhdr, nullify_HScr, scr_hdr_io, read_screening, write_screening, &
&                             copy_scrhdr, HSCR_LATEST_HEADFORM, scrhdr_type
 use m_spectra,        only : spectra_type, init_spectra, destroy_spectra, repr_dielconst
 use m_special_funcs,  only : ylmc

 implicit none

 private 
!!***

!----------------------------------------------------------------------

!!****t* m_screening/epsilonm1_results
!! NAME
!! epsilonm1_results
!!
!! FUNCTION
!! For the GW part of ABINIT, the epsilonm1_results structured datatype
!! gather the results of screening : the inverse dielectric matrix,
!! and the omega matrices .
!!
!! SOURCE

 type,public :: Epsilonm1_results

! WARNING : if you modify this datatype, please check there there is no creation/destruction/copy routine,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: ID                          ! Matrix identifier: O if not yet defined, 1 for chi0,
                                         ! 2 for chi, 3 for epsilon, 4 for espilon^{-1}, 5 for W.
  integer :: ikxc                        ! Kxc kernel used, 0 for None (RPA), >0 for static TDDFT (=ixc), <0 for TDDFT
  integer :: fform                       ! File format: 1002 for SCR|SUSC files. 2002 for pole fit
  integer :: mqmem                       ! =0 for out-of-core solution, =nqibz if entire matrix is stored in memory.
  integer :: nI,nJ                       ! Number of components (rows,columns) in chi|eps^-1. (1,1) if collinear.
  integer :: nqibz                       ! Number of q-points in the IBZ used.
  integer :: nqlwl                       ! Number of point used for the treatment of the long wave-length limit.
  integer :: nomega                      ! Number of frequencies used.
  integer :: nomega_i                    ! Number of purely imaginary frequencies used.
  integer :: nomega_r                    ! Number of real frequencies used.
  integer :: nomega_c                    ! Number of frequencies in the complex plane.
  integer :: npoles                      ! Number of poles for pole-fit screening
  integer :: ncoeff                      ! Number of coefficients = npoles*3+1 for phase
  integer :: npwe                        ! Number of G vectors used.
  integer :: test_type                   ! 0 for None, 1 for TEST-PARTICLE, 2 for TEST-ELECTRON (only for TDDFT)
  integer :: Tordering                   ! 0 if not defined, 1 for Time-Ordered, 2 for Advanced, 3 for Retarded.

  character(len=fnlen) :: fname          ! Name of the file from which epsm1 is read.

!arrays
  integer,pointer  :: gvec(:,:)   SET2NULL
  ! gvec(3,npwe)
  ! G-vectors used to describe the two-point function (r.l.u.).

  real(dp),pointer :: qibz(:,:)  SET2NULL
  ! qibz(3,nqibz)
  ! q-points in reduced coordinates

  real(dp),pointer :: qlwl(:,:) SET2NULL
  ! qlwl(3,nqlwl)
  ! q-points used for the long wave-length limit treatment.

  real(gwp),pointer :: epsm1_pole(:,:,:,:)  SET2NULL
  ! epsm1(npwe,npwe,nomega,nqibz)
  ! Contains the two-point function $\epsilon_{G,Gp}(q,omega)$ in frequency and reciprocal space.

  complex(gwpc),pointer :: epsm1(:,:,:,:)  SET2NULL
  ! epsm1(npwe,npwe,nomega,nqibz)
  ! Contains the two-point function $\epsilon_{G,Gp}(q,omega)$ in frequency and reciprocal space.

  complex(dpc),pointer :: lwing(:,:,:)  SET2NULL
  ! lwing(npwe,nomega,nqlwl)
  ! Lower wings for the different q"s -->0

  complex(dpc),pointer :: omega(:)  SET2NULL
  ! omega(nomega)
  ! Frequencies used both along the real and the imaginary axis.

  complex(dpc),pointer :: uwing(:,:,:)  SET2NULL
  ! uwing(npwe,nomega,nqlwl)
  ! Upper wings for the different q"s -->0

  type(ScrHdr_type) :: Hscr
  ! The header reported in the _SCR of _SUSC file.
  ! This object contains information on the susceptibility or the inverse dielectric matrix
  ! as stored in the external file. These quantities do *NOT* correspond to the quantities
  ! used during the GW calculation since some parameters might differ, actually they might be smaller.
  ! For example, the number of G-vectors used can be smaller than the number of G"s stored on file.

 end type Epsilonm1_results
!!***

!----------------------------------------------------------------------

 public :: nullify_epsilonm1_results     ! Nullify all pointers before use.
 public :: destroy_epsilonm1_results     ! Free all associated pointers
 public :: print_epsilonm1_results       ! Print basic info 
 public :: Epsm1_symmetrizer             ! Symmetrize two-point function at a q-point in the BZ.
 public :: Epsm1_symmetrizer_inplace     ! In-place version of the above
 public :: Epsm1_pole_symmetrizer         ! Symmetrize pole fit screening
 public :: Epsm1_pole_symmetrizer_inplace ! In-place version of the above
 public :: init_Er_from_file             ! Initialize the object from file
 public :: mkdump_Er                     ! Dump the object to a file.
 public :: get_epsm1                    
 public :: get_pole_epsm1                    
 public :: decompose_epsm1
 public :: outeps
 public :: symf12
 public :: make_epsm1_driver
 public :: make_W                        ! Calculate W from the data stored in Er (in place, content of Er% is changed).
 public :: mkem1_q0
 public :: screen_mdielf
 public :: interpolate_w
 public :: recalculate_epsm1_freq_grid

CONTAINS  !========================================================================================
!!***

!!****f* m_screening/nullify_epsilonm1_results
!! NAME
!! nullify_epsilonm1_results
!!
!! FUNCTION
!! Initialize the pointer to null()
!!
!! INPUTS
!! Er<Epsilonm1_results>=The data structure.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screening,mrgscr,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_epsilonm1_results(Er)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_epsilonm1_results'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_results),intent(inout) :: Er
! *************************************************************************

 !@Epsilonm1_results
 !integer
 nullify(Er%gvec)

 !real
 nullify(Er%qibz)
 nullify(Er%qlwl) 
 nullify(Er%epsm1_pole) 

 !complex
 nullify(Er%epsm1) 
 nullify(Er%lwing)
 nullify(Er%omega)             
 nullify(Er%uwing)       

 ! datatypes
 call nullify_HScr(Er%HScr) ! All pointers are always allocates (except PAW rhoij) but nullification might be useful!

end subroutine nullify_epsilonm1_results
!!***

!----------------------------------------------------------------------

!!****f* m_screening/destroy_epsilonm1_results
!! NAME
!! destroy_epsilonm1_results
!!
!! FUNCTION
!! Deallocate all the pointers in Er that result to be associated.
!! Perform also a cleaning of the Header.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screening,mrgscr,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_epsilonm1_results(Er)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_epsilonm1_results'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_results),intent(inout) :: Er
! *************************************************************************

 !@Epsilonm1_results
 !integer
 if (associated(Er%gvec ))  then
   ABI_DEALLOCATE(Er%gvec)
 end if

 !real
 if (associated(Er%qibz ))  then
   ABI_DEALLOCATE(Er%qibz)
 end if
 if (associated(Er%qlwl ))  then
   ABI_DEALLOCATE(Er%qlwl)
 end if
 if (associated(Er%epsm1_pole ))  then
   ABI_DEALLOCATE(Er%epsm1_pole)
 end if

 !complex
 if (associated(Er%epsm1))  then
   ABI_DEALLOCATE(Er%epsm1)
 end if
 if (associated(Er%lwing))  then
   ABI_DEALLOCATE(Er%lwing)
 end if
 if (associated(Er%omega))  then
   ABI_DEALLOCATE(Er%omega)
 end if
 if (associated(Er%uwing))  then
   ABI_DEALLOCATE(Er%uwing)
 end if

 !datatypes
 call free_scrhdr(Er%Hscr)

end subroutine destroy_epsilonm1_results
!!***

!----------------------------------------------------------------------

!!****f* m_screening/print_epsilonm1_results
!! NAME
!!  print_epsilonm1_results
!!
!! FUNCTION
!! Print the basic dimensions and the most important 
!! quantities reported in the Epsilonm1_results data type.
!!
!! INPUTS
!!  Er<Epsilonm1_results>=The data type.
!!  unit[optional]=the unit number for output.
!!  prtvol[optional]=verbosity level.
!!  mode_paral[optional]=either COLL or PERS.
!!
!! OUTPUT
!!  Only printing. 
!!
!! PARENTS
!!      m_screening,mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_epsilonm1_results(Er,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_epsilonm1_results'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 type(Epsilonm1_results),intent(in) :: Er

!Local variables-------------------------------
 integer :: iomega,iqibz,iqlwl,unt,verbose,rdwr
 character(len=50) :: rfname,rforder,rfapprox,rftest,kxcname
 character(len=500) :: msg
 character(len=4) :: mode
! *************************************************************************

 unt    =std_out; if (PRESENT(unit      )) unt    =unit
 verbose=0      ; if (PRESENT(prtvol    )) verbose=prtvol
 mode   ='COLL' ; if (PRESENT(mode_paral)) mode   =mode_paral

 ! === chi0 or \epsilon^{-1} ? ===
 SELECT CASE (Er%ID)
 CASE (0)
   rfname='Undefined'
 CASE (1)
   rfname='Irreducible Polarizability'
 CASE (2)
   rfname='Polarizability'
 CASE (3)
   rfname='Symmetrical Dielectric Matrix'
 CASE (4)
   rfname='Symmetrical Inverse Dielectric Matrix'
 CASE DEFAULT 
   write(msg,'(a,i3)') 'Wrong value of Er%ID= ',Er%ID
   MSG_BUG(msg)
 END SELECT

 ! === For chi, \espilon or \epsilon^{-1}, define the approximation ===
 rfapprox='None'
 if (Er%ID>=2.or.Er%ID<=4) then
   if (Er%ikxc==0) then 
     rfapprox='RPA'
   else if (Er%ikxc>0) then 
     rfapprox='Static TDDFT'
   else  
     rfapprox='TDDFT'
   end if
 end if

 ! === If TDDFT and \epsilon^{-1}, define the type ===
 rftest='None'
! if (Er%ID==0) then
!  if (Er%test_type==0) then 
!   rftest='TEST-PARTICLE'
!  else if (Er%test_type==1) then 
!   rftest='TEST-ELECTRON'
!  else 
!   write(msg,'(4a,i3)')ch10,&
!&   ' print_epsilonm1_results : BUG - ',ch10,&
!&   ' Wrong value of Er%test_type = ',Er%test_type
!   MSG_ERROR(msg)
!  end if
! end if

 ! === Define time-ordering ===
 rforder='Undefined'
 if (Er%Tordering==1) then 
   rforder='Time-Ordered'
 else if (Er%Tordering==2) then 
   rforder='Advanced'
 else if (Er%Tordering==3) then 
   rforder='Retarded'
 else  
   write(msg,'(a,i3)')'Wrong value of Er%Tordering= ',Er%Tordering
   MSG_BUG(msg)
 end if

 kxcname='None'
 if (Er%ikxc/=0) then 
   !TODO Add function to retrieve kxc name
   STOP 'Add function to retrieve kxc name'
   kxcname='XXXXX'
 end if

 write(msg,'(6a,5(3a))')ch10,&
&  ' ==== Info on the Response Function ==== ',ch10,&
&  '  Associated File ................  ',TRIM(Er%fname),ch10,&
&  '  Response Function Type .......... ',TRIM(rfname),ch10,&
&  '  Type of Approximation ........... ',TRIM(rfapprox),ch10,&
&  '  XC kernel used .................. ',TRIM(kxcname),ch10,&
&  '  Type of probing particle ........ ',TRIM(rftest),ch10,&
&  '  Time-Ordering ................... ',TRIM(rforder),ch10
 call wrtout(unt,msg,mode) 
 write(msg,'(a,2i4,a,3(a,i4,a),a,3i4,2a,i4,a)')&
&  '  Number of components ............ ',Er%nI,Er%nJ,ch10,&
&  '  Number of q-points in the IBZ ... ',Er%nqibz,ch10,&
&  '  Number of q-points for q-->0 .... ',Er%nqlwl,ch10,&
&  '  Number of G-vectors ............. ',Er%npwe,ch10,&
&  '  Number of frequencies ........... ',Er%nomega,Er%nomega_r,Er%nomega_i,ch10,&
&  '  Value of mqmem .................. ',Er%mqmem,ch10
 call wrtout(unt,msg,mode) 
 if (Er%npoles>0) then ! We have a pole-fit
   write(msg,'(3a,i4,2a,i4,3a)')&
&    '  - THIS IS A POLE-FIT SCREENING -  ',ch10,&
&    '  Number of poles ................. ',Er%npoles,ch10,&
&    '  Number of coefficients .......... ',Er%ncoeff,ch10,&
&    '  (npoles*3 + 1 value for phase)... ',ch10
   call wrtout(unt,msg,mode) 
 end if

 if (Er%nqlwl/=0) then 
   write(msg,'(a,i3)')' q-points for long wavelength limit: ',Er%nqlwl
   call wrtout(unt,msg,mode) 
   do iqlwl=1,Er%nqlwl
     write(msg,'(1x,i5,a,3es16.8)')iqlwl,') ',Er%qlwl(:,iqlwl)
     call wrtout(unt,msg,mode)
   end do
 end if

 if (verbose>0) then ! Print out head and wings in the long-wavelenght limit.
   ! TODO add additional stuff.
   if (Er%nqlwl>0) then 
     write(msg,'(1x,2a)')' Heads and wings of chi0(G,G'')',ch10
     call wrtout(unt,msg,mode)
     do iqlwl=1,Er%nqlwl
       write(msg,'(1x,a,i2,a)')' chi0(qlwl =',iqlwl,')'
       call wrtout(unt,msg,mode)
       do iomega=1,Er%nomega
         write(msg,'(2x,a,i4,a,2f9.4,a)')&
&          ' Upper and lower wings at the ',iomega,' th omega',Er%omega(iomega)*Ha_eV,' [eV]'
         call wrtout(unt,msg,mode)
         call print_arr(Er%uwing(:,iomega,iqlwl),max_r=9,unit=unt)
         call print_arr(Er%lwing(:,iomega,iqlwl),max_r=9,unit=unt)
       end do
     end do
   end if

   write(msg,'(a,i4)')' Calculated Frequencies: ',Er%nomega
   call wrtout(unt,msg,mode) 
   do iomega=1,Er%nomega 
     write(msg,'(i4,es14.6)')iomega,Er%omega(iomega)*Ha_eV
     call wrtout(unt,msg,mode) 
   end do

   write(msg,'(a,i4)')' Calculated q-points: ',Er%nqibz
   call wrtout(unt,msg,mode) 
   do iqibz=1,Er%nqibz
     write(msg,'(1x,i4,a,3es16.8)')iqibz,') ',Er%qibz(:,iqibz)
     call wrtout(unt,msg,mode)
   end do

   rdwr=4
   !$call hdr_io_int(Er%fform,Er%Hscr%Hdr,rdwr,unt)
 end if ! verbose>0

end subroutine print_epsilonm1_results 
!!***

!----------------------------------------------------------------------

!!****f* m_screening/Epsm1_symmetrizer
!! NAME
!!  Epsm1_symmetrizer
!!
!! FUNCTION
!!  Symmetrize the inverse dielectric matrix, namely calculate epsilon^{-1} at a generic 
!!  q-point in the BZ starting from the knowledge of the matrix at a q-point in the IBZ.
!!  The procedure is quite generic and can be used for every two-point function which has 
!!  the same symmetry as the crystal. 
!!
!! INPUTS
!!  nomega=Number of frequencies required. All frequencies from 1 up to nomega are symmetrized.
!!  npwc=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  remove_exchange=If .TRUE., return e^{-1}-1 namely remove the exchange part.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<Gvectors_type>=data related to the G-sphere
!!    %grottb
!!    %phmSGt 
!!  Qmesh<BZ_mesh_type>=Structure defining the q-mesh used for Er.
!!    %nbz=Number of q-points in the BZ
!!    %tab(nbz)=Index of the symmetric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=The operation that rotates q_ibz onto \pm q_bz (depending on tabi) 
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required. 
!!
!! OUTPUT
!!  epsm1_qbz(npwc,npwc,nomega)=The inverse dielectric matrix at the q-point defined by iq_bz. 
!!   Exchange part can be subtracted out.
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the 
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere 
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required 
!!  to reconstruct the BZ.
!! 
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!! 
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine Epsm1_symmetrizer(iq_bz,nomega,npwc,Er,Gsph,Qmesh,remove_exchange,epsm1_qbz) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Epsm1_symmetrizer'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npwc
 logical,intent(in) :: remove_exchange
 type(Epsilonm1_results),intent(in) :: Er
 type(Gvectors_type),intent(in) :: Gsph
 type(BZ_mesh_type),intent(in) :: Qmesh
!arrays
 complex(gwpc),intent(out) :: epsm1_qbz(npwc,npwc,nomega)  

!Local variables-------------------------------
!scalars
 integer :: iomega,ii,jj,iq_ibz,itim_q,isym_q,iq_loc
!arrays
 integer,pointer :: grottb(:)
 real(dp) :: qbz(3)
 complex(gwpc),pointer :: phmSgt(:)

! *********************************************************************

 ABI_CHECK(Er%nomega>=nomega,'Too many frequencies required')
 ABI_CHECK(Er%npwe  >=npwc , 'Too many G-vectors required')

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

 grottb => Gsph%rottb (1:npwc,itim_q,isym_q)
 phmSgt => Gsph%phmSGt(1:npwc,isym_q) 

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled. 
 iq_loc=iq_ibz ; if (Er%mqmem==0) iq_loc=1 

 do iomega=1,nomega
   do jj=1,npwc
     do ii=1,npwc
       epsm1_qbz(grottb(ii),grottb(jj),iomega)=Er%epsm1(ii,jj,iomega,iq_loc)*phmSgt(ii)*CONJG(phmSgt(jj))
     end do
   end do
 end do
 !
 ! === Account for time-reversal ===
 if (itim_q==2) then
   do iomega=1,nomega
     epsm1_qbz(:,:,iomega)=TRANSPOSE(epsm1_qbz(:,:,iomega))
   end do
 end if

 if(remove_exchange)then
   ! === Subtract the exchange contribution ===
   ! If it's a pole screening, the exchange contribution is already removed
   do iomega=1,nomega
     do ii=1,npwc
       epsm1_qbz(ii,ii,iomega)=epsm1_qbz(ii,ii,iomega)-1.0_gwp
     end do
   end do
 endif

end subroutine Epsm1_symmetrizer
!!***

!----------------------------------------------------------------------

!!****f* m_screening/Epsm1_symmetrizer_inplace
!! NAME
!!  Epsm1_symmetrizer_inplace
!!
!! FUNCTION
!!  Same function as Epsm1_symmetrizer, ecept now the array Ep%epsm1 is modified inplace
!!  thorugh an auxiliary work array of dimension (npwc,npwc) 
!!
!! INPUTS
!!  nomega=Number of frequencies required. All frequencies from 1 up to nomega are symmetrized.
!!  npwc=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  remove_exchange=If .TRUE., return e^{-1}-1 namely remove the exchange part.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<Gvectors_type>=data related to the G-sphere
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<Gvectors_type>=data related to the G-sphere
!!    %grottb
!!    %phmSGt 
!!  Qmesh<BZ_mesh_type>=Structure defining the q-mesh used for Er.
!!    %nbz=Number of q-points in the BZ
!!    %tab(nbz)=Index of the symmetric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=The operation that rotates q_ibz onto \pm q_bz (depending on tabi) 
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required. 
!!
!! OUTPUT
!!  Er%epsm1(npwc,npwc,nomega,iq_loc) symmetrised
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the 
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere 
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required 
!!  to reconstruct the BZ.
!! 
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!! 
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine Epsm1_symmetrizer_inplace(iq_bz,nomega,npwc,Er,Gsph,Qmesh,remove_exchange) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Epsm1_symmetrizer_inplace'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npwc
 logical,intent(in) :: remove_exchange
 type(Epsilonm1_results) :: Er
 type(Gvectors_type),intent(in) :: Gsph
 type(BZ_mesh_type),intent(in) :: Qmesh

!Local variables-------------------------------
!scalars
 integer :: iomega,ii,jj,iq_ibz,itim_q,isym_q,iq_loc
!arrays
 integer,pointer :: grottb(:)
 real(dp) :: qbz(3)
 complex(gwpc),pointer :: phmSgt(:)
 complex(gwpc),allocatable :: work(:,:)

! *********************************************************************

 ABI_CHECK(Er%nomega>=nomega,'Too many frequencies required')
 ABI_CHECK(Er%npwe  >=npwc , 'Too many G-vectors required')

 ABI_ALLOCATE(work,(npwc,npwc))

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

 grottb => Gsph%rottb (1:npwc,itim_q,isym_q)
 phmSgt => Gsph%phmSGt(1:npwc,isym_q) 

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled. 
 iq_loc=iq_ibz ; if (Er%mqmem==0) iq_loc=1 

 do iomega=1,nomega
   do jj=1,npwc
     do ii=1,npwc
       work(grottb(ii),grottb(jj))=Er%epsm1(ii,jj,iomega,iq_loc)*phmSgt(ii)*CONJG(phmSgt(jj))
     end do
   end do
   Er%epsm1(:,:,iomega,iq_loc) = work(:,:)
 end do
 !
 ! === Account for time-reversal ===
 if (itim_q==2) then
   do iomega=1,nomega
     Er%epsm1(:,:,iomega,iq_loc)=TRANSPOSE(Er%epsm1(:,:,iomega,iq_loc))
   end do
 end if

 ! === Subtract the exchange contribution ===
 if(remove_exchange)then
   do iomega=1,nomega
     do ii=1,npwc
       Er%epsm1(ii,ii,iomega,iq_loc)=Er%epsm1(ii,ii,iomega,iq_loc)-1.0_gwp
     end do
   end do
 endif 

 ABI_DEALLOCATE(work)

end subroutine Epsm1_symmetrizer_inplace
!!***

!----------------------------------------------------------------------

!!****f* m_screening/Epsm1_pole_symmetrizer
!! NAME
!!  Epsm1_symmetrizer
!!
!! FUNCTION
!!  Symmetrize the inverse dielectric matrix, namely calculate epsilon^{-1} at a generic 
!!  q-point in the BZ starting from the knowledge of the matrix at a q-point in the IBZ.
!!  The procedure is quite generic and can be used for every two-point function which has 
!!  the same symmetry as the crystal. This routine is meant for pole-fit screening 
!!
!! INPUTS
!!  ncoeff=Number of coefficients in pole-fit + phase value
!!  npwc=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<Gvectors_type>=data related to the G-sphere
!!    %grottb
!!    %phmSGt 
!!  Qmesh<BZ_mesh_type>=Structure defining the q-mesh used for Er.
!!    %nbz=Number of q-points in the BZ
!!    %tab(nbz)=Index of the symmetric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=The operation that rotates q_ibz onto \pm q_bz (depending on tabi) 
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required. 
!!
!! OUTPUT
!!  epsm1_pole_qbz(npwc,npwc,nomega)=The inverse dielectric matrix at the q-point defined by iq_bz. 
!!   Exchange part can be subtracted out.
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the 
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere 
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required 
!!  to reconstruct the BZ.
!! 
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!! 
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! CHILDREN
!!      copy_scrhdr,deep_copy,destroy_bz_mesh_type,destroy_crystal
!!      destroy_gsphere,find_qmesh,free_scrhdr,get_bz_item
!!      init_crystal_from_hdr,init_gsphere,init_kmesh,print_bz_mesh
!!      read_screening,scr_hdr_io,sort_dp,wrap2_pmhalf,write_screening,wrtout
!!
!! SOURCE

subroutine Epsm1_pole_symmetrizer(iq_bz,ncoeff,npwc,Er,Gsph,Qmesh,epsm1_pole_qbz) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Epsm1_pole_symmetrizer'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,ncoeff,npwc
 type(Epsilonm1_results),intent(in) :: Er
 type(Gvectors_type),intent(in) :: Gsph
 type(BZ_mesh_type),intent(in) :: Qmesh
!arrays
 real(gwp),intent(out) :: epsm1_pole_qbz(npwc,npwc,ncoeff)  

!Local variables-------------------------------
!scalars
 integer :: ii,jj,iq_ibz,itim_q,isym_q,iq_loc,icoeff
 real(gwp) :: phase
 complex   :: cphase
!arrays
 integer,pointer :: grottb(:)
 real(dp) :: qbz(3)
 complex(gwpc),pointer :: phmSgt(:)

! *********************************************************************

 ABI_CHECK(Er%npwe  >=npwc , 'Too many G-vectors required')

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

 grottb => Gsph%rottb (1:npwc,itim_q,isym_q)
 phmSgt => Gsph%phmSGt(1:npwc,isym_q) 

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled. 
 iq_loc=iq_ibz ; if (Er%mqmem==0) iq_loc=1 

 do jj=1,npwc
   do ii=1,npwc
     do icoeff=1,ncoeff-1,3
       epsm1_pole_qbz(grottb(ii),grottb(jj),icoeff)=Er%epsm1_pole(ii,jj,icoeff,iq_loc)
       epsm1_pole_qbz(grottb(ii),grottb(jj),icoeff+1)=Er%epsm1_pole(ii,jj,icoeff+1,iq_loc)
       epsm1_pole_qbz(grottb(ii),grottb(jj),icoeff+2)=Er%epsm1_pole(ii,jj,icoeff+2,iq_loc)
     end do
     ! Add phase
     cphase = phmSgt(ii)*CONJG(phmSgt(jj))
     phase  = ATAN(AIMAG(cphase)/REAL(cphase))
     epsm1_pole_qbz(grottb(ii),grottb(jj),ncoeff)=Er%epsm1_pole(ii,jj,ncoeff,iq_loc)-phase
   end do
 end do
 !
 ! === Account for time-reversal ===
 if (itim_q==2) then
   do icoeff=1,ncoeff
     epsm1_pole_qbz(:,:,icoeff)=TRANSPOSE(epsm1_pole_qbz(:,:,icoeff))
   end do
 end if

end subroutine Epsm1_pole_symmetrizer
!!***

!----------------------------------------------------------------------

!!****f* m_screening/Epsm1_pole_symmetrizer_inplace
!! NAME
!!  Epsm1_symmetrizer
!!
!! FUNCTION
!!  Symmetrize the inverse dielectric matrix, namely calculate epsilon^{-1} at a generic 
!!  q-point in the BZ starting from the knowledge of the matrix at a q-point in the IBZ.
!!  The procedure is quite generic and can be used for every two-point function which has 
!!  the same symmetry as the crystal. This routine is meant for pole-fit screening 
!!
!! INPUTS
!!  ncoeff=Number of coefficients in pole-fit + phase value
!!  npwc=Number of G vectors in symmetrized matrix, has to be smaller than Er%npwe.
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix.
!!  Gsph<Gvectors_type>=data related to the G-sphere
!!    %grottb
!!    %phmSGt 
!!  Qmesh<BZ_mesh_type>=Structure defining the q-mesh used for Er.
!!    %nbz=Number of q-points in the BZ
!!    %tab(nbz)=Index of the symmetric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=The operation that rotates q_ibz onto \pm q_bz (depending on tabi) 
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required. 
!!
!! OUTPUT
!!  epsm1_pole_qbz(npwc,npwc,nomega)=The inverse dielectric matrix at the q-point defined by iq_bz. 
!!   Exchange part can be subtracted out.
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the 
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere 
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required 
!!  to reconstruct the BZ.
!! 
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!! 
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! CHILDREN
!!      copy_scrhdr,deep_copy,destroy_bz_mesh_type,destroy_crystal
!!      destroy_gsphere,find_qmesh,free_scrhdr,get_bz_item
!!      init_crystal_from_hdr,init_gsphere,init_kmesh,print_bz_mesh
!!      read_screening,scr_hdr_io,sort_dp,wrap2_pmhalf,write_screening,wrtout
!!
!! SOURCE

subroutine Epsm1_pole_symmetrizer_inplace(iq_bz,ncoeff,npwc,Er,Gsph,Qmesh) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Epsm1_pole_symmetrizer_inplace'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,ncoeff,npwc
 type(Epsilonm1_results) :: Er
 type(Gvectors_type),intent(in) :: Gsph
 type(BZ_mesh_type),intent(in) :: Qmesh
!arrays

!Local variables-------------------------------
!scalars
 integer :: ii,jj,iq_ibz,itim_q,isym_q,iq_loc,icoeff
 real(gwp) :: phase
 complex   :: cphase
!arrays
 integer,pointer :: grottb(:)
 real(dp) :: qbz(3)
 complex(gwpc),pointer :: phmSgt(:)
 real(gwp), allocatable :: work(:,:)

! *********************************************************************

 ABI_CHECK(Er%npwe  >=npwc , 'Too many G-vectors required')

 ! * Get iq_ibz, and symmetries from iq_ibz.
 call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)

 ABI_ALLOCATE(work,(npwc,npwc))

 grottb => Gsph%rottb (1:npwc,itim_q,isym_q)
 phmSgt => Gsph%phmSGt(1:npwc,isym_q) 

 ! If out-of-memory, only Er%espm1(:,:,:,1) has been allocated and filled. 
 iq_loc=iq_ibz ; if (Er%mqmem==0) iq_loc=1 

 do icoeff=1,ncoeff-1
   do jj=1,npwc
     do ii=1,npwc
         work(grottb(ii),grottb(jj))=Er%epsm1_pole(ii,jj,icoeff,iq_loc)
     end do
   end do
   Er%epsm1_pole(:,:,icoeff,iq_loc)=work(:,:)
 end do
 ! Add phase
 do jj=1,npwc
   do ii=1,npwc
     cphase = phmSgt(ii)*CONJG(phmSgt(jj))
     phase  = ATAN(AIMAG(cphase)/REAL(cphase))
     work(grottb(ii),grottb(jj))=Er%epsm1_pole(ii,jj,ncoeff,iq_loc)+phase
   end do
 end do
 Er%epsm1_pole(:,:,ncoeff,iq_loc) = work(:,:)
 !
 ! === Account for time-reversal ===
 if (itim_q==2) then
   do icoeff=1,ncoeff
     Er%epsm1_pole(:,:,icoeff,iq_loc)=TRANSPOSE(Er%epsm1_pole(:,:,icoeff,iq_loc))
   end do
 end if

 ABI_DEALLOCATE(work)

end subroutine Epsm1_pole_symmetrizer_inplace
!!***

!----------------------------------------------------------------------

!!****f* m_screening/init_Er_from_file
!! NAME
!!  init_Er_from_file
!!
!! FUNCTION
!!  Initialize basic dimensions and the important (small) arrays in an Epsilonm1_results data type
!!  starting from a file containing either epsilon^{-1} (_SCR) or chi0 (_SUSC).
!!
!! INPUTS
!!  accesswff=Option defining the file format of the external file.
!!  fname=The name of the external file used to read the matrix.
!!  mqmem=0 for out-of-core solution, /=0 if entire matrix has to be stored in memory. 
!!  npwe_asked=Number of G-vector to be used in the calculation, if <=0 use Max allowed number.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Er<Epsilonm1_results>=The structure initialized with basic dimensions and arrays.
!!
!! PARENTS
!!      m_screening,mrgscr,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_Er_from_file(Er,fname,mqmem,npwe_asked,accesswff,comm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_Er_from_file'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: mqmem,accesswff,npwe_asked,comm
 character(len=fnlen),intent(in) :: fname
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: iomega,ios,fform,rdwr,my_rank,master,unt
 character(len=500) :: msg                   

! *********************************************************************

 DBG_ENTER("COLL")

 !@Epsilonm1_results
 my_rank = xcomm_rank(comm)
 master=0

 call nullify_epsilonm1_results(Er)

 ! === Open file ===
 !if (my_rank==master.or.localrdwf==1) then
   write(msg,'(3a)')' init_Er_from_file : testing file ',TRIM(fname),ch10
   call wrtout(std_out,msg,'COLL')
   unt=get_unit()
   open(unit=unt,file=fname,status='old',form='unformatted',iostat=ios)
   msg = ' Opening file '//TRIM(fname)//' as old '
   ABI_CHECK(ios==0,msg)
 !end if

 rdwr=5
 call scr_hdr_io(fform,rdwr,unt,comm,master,accesswff,Er%Hscr)

 !if (my_rank==master.or.localrdwf==1) close(unt)
 close(unt)

 ! === Master echoes the header ===
 if (my_rank==master) then
   rdwr=4
   call scr_hdr_io(fform,rdwr,std_out,comm,master,accesswff,Er%Hscr)
 end if

 ! === Generic Info ===
 Er%ID         =0       ! Not yet initialized as epsm1 is calculated in mkdump_Er.F90
 Er%fname      =fname
 Er%fform      =fform
 if (fform==2002) then
   Er%npoles     =Er%Hscr%npoles
   Er%ncoeff     =Er%Hscr%ncoeff
 else
   Er%npoles     =0
   Er%ncoeff     =0
 end if

 Er%Tordering=Er%Hscr%Tordering

!TODO these quantitities should be checked and initiliazed in mkdump_Er
!BEGIN HARCODED
 Er%nI       = 1
 Er%nJ       = 1
 Er%ikxc     = 0
 Er%test_type=-1    

 Er%Hscr%headform=HSCR_LATEST_HEADFORM   ! XG20090912 
!END HARDCODED

 Er%nqibz=Er%Hscr%nqibz
 Er%mqmem=mqmem ; if (mqmem/=0) Er%mqmem=Er%nqibz
 ABI_ALLOCATE(Er%qibz,(3,Er%nqibz))
 Er%qibz(:,:)=Er%Hscr%qibz(:,:)

 Er%nqlwl=Er%Hscr%nqlwl
 ABI_ALLOCATE(Er%qlwl,(3,Er%nqlwl))
 Er%qlwl(:,:)=Er%Hscr%qlwl(:,:)

 Er%nomega=Er%Hscr%nomega
 ABI_ALLOCATE(Er%omega,(Er%nomega))
 Er%omega(:)=Er%Hscr%omega(:)

 if (Er%nomega==2) then
   Er%nomega_r=1 
   Er%nomega_i=1
 else ! Real frequencies are packed in the first locations.
   Er%nomega_r=1
   do iomega=1,Er%nomega
     if ((REAL(Er%omega(iomega))>0.001*Ha_eV).AND.&
&     (AIMAG(Er%omega(iomega))<0.001*Ha_eV)) Er%nomega_r=iomega
   end do
   Er%nomega_i=0
   do iomega=Er%nomega_r+1,Er%nomega
     if ((REAL(Er%omega(iomega))<0.001*Ha_eV).AND.&
&     (AIMAG(Er%omega(iomega))>0.001*Ha_eV)) Er%nomega_i=Er%nomega_i+1
   end do
   Er%nomega_c=Er%nomega-Er%nomega_r-Er%nomega_i
 end if     

 ! === Get G-vectors ===
 Er%npwe=Er%Hscr%npwe
 if (npwe_asked>0) then
   if (npwe_asked>Er%Hscr%npwe) then
     write(msg,'(a,i8,2a,i8)')&
&     '  Number of G-vectors saved on file is less than the value required = ',npwe_asked,ch10,&
&     '  Calculation will proceed with Max available npwe = ',Er%Hscr%npwe
     MSG_WARNING(msg)
   else  ! Redefine the no. of G"s for W.
     Er%npwe=npwe_asked
   end if
 end if

 ! pointer to Er%Hscr%gvec ?
 ABI_ALLOCATE(Er%gvec,(3,Er%npwe))
 Er%gvec=Er%Hscr%gvec(:,1:Er%npwe)

 DBG_EXIT("COLL")

end subroutine init_Er_from_file
!!***

!----------------------------------------------------------------------

!!****f* m_screening/mkdump_Er
!! NAME
!!  mkdump_Er
!!
!! FUNCTION
!!  Dump the content of an Epsilonm1_results data type on file.
!!
!! INPUTS
!!  id_required=Identifier of the matrix to be calculated
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  ngfft(18)=Info on the FFT mesh.
!!  nfftot=Total number of point on the FFT mesh.
!!  gvec(3,npwe)=Reduced coordinates of plane waves for the response functions
!!  npwe=Number of plane waves.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!      mrgscr,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkdump_Er(Er,Vcp,npwe,gvec,nkxc,kxcg,id_required,approx_type,&
&                    ikxc_required,option_test,fname_dump,accesswff,&
&                    nfftot,ngfft,comm,fxc_ADA,reconstruct_scr)

 use defs_basis
 use m_io_screening, only : read_pole_screening
 use m_model_screening, only : re_and_im_screening_with_phase

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkdump_Er'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: id_required,approx_type,option_test,ikxc_required,nkxc
 integer,intent(in) :: accesswff,nfftot,npwe,comm
 integer,intent(in),optional :: reconstruct_scr
 type(Epsilonm1_results),intent(inout) :: Er
 type(vcoul_t),intent(in) :: Vcp
 character(len=fnlen),intent(in) :: fname_dump
!arrays 
 integer,intent(in) :: ngfft(18),gvec(3,npwe)
 complex(gwpc),intent(in) :: kxcg(nfftot,nkxc) 
 complex(gwpc),intent(in), optional :: fxc_ADA(Er%npwe*Er%nI,Er%npwe*Er%nJ,Er%nqibz)

!Local variables-------------------------------
!scalars
 integer :: dim_wing,ios,istat,iqibz,is_qeq0,mqmem_,npwe_asked
 integer :: unt_dump,fform,rdwr,ig1,ig2 
 integer :: master,my_rank,comm_self
 real(dp) :: ucvol
 character(len=500) :: msg                   
 type(ScrHdr_type) :: Hscr_cp
 type(Spectra_type) :: Spectra
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 complex(gwpc),pointer :: epsm1(:,:,:)
 complex(dpc),allocatable :: dummy_lwing(:,:,:),dummy_uwing(:,:,:),dummy_head(:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")
 
 ABI_CHECK(id_required==4,'Value of id_required not coded')
 ABI_CHECK(npwe==Er%npwe,"mismatch in npwe")

 my_rank = xcomm_rank(comm)
 master=0
 comm_self = xmpi_self

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 !! if (Er%ID/=0) call reset_Epsilonm1(Er)
 Er%ID=id_required

 write(std_out,*) 'Er%ID:',Er%ID
 write(std_out,*) 'Er%Hscr%ID:',Er%Hscr%ID

 if (Er%ID==Er%Hscr%ID) then
   ! === The two-point function we are asking for is already stored on file ===
   ! * According to mqmem either read and store the entire matrix in memory or do nothing.
   !
   if (Er%mqmem>0) then ! In-core solution.
     ABI_ALLOCATE(Er%lwing,(Er%npwe,Er%nomega,Er%nqlwl))
     ABI_ALLOCATE(Er%uwing,(Er%npwe,Er%nomega,Er%nqlwl))
     if (Er%nqlwl>0) then 
       !pointer ?
       Er%lwing(:,:,:)=Er%Hscr%lwing(1:Er%npwe,1:Er%nomega,1:Er%nqlwl)
       Er%uwing(:,:,:)=Er%Hscr%uwing(1:Er%npwe,1:Er%nomega,1:Er%nqlwl)
     end if

     if (Er%fform/=2002) then ! Check for pole-fit scr
       ABI_ALLOCATE(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,Er%nqibz))
       istat = ABI_ALLOC_STAT
       if (istat/=0) then 
         MSG_ERROR('Out-of-memory in Er%epsm1 (in-core)')
       end if
       call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,accesswff,comm)
     else
       !write(std_out,*) Er%fform
       !write(std_out,*) Er%npwe,Er%npwe,Er%ncoeff,Er%nqibz
       ABI_ALLOCATE(Er%epsm1_pole,(Er%npwe,Er%npwe,Er%ncoeff,Er%nqibz))
       istat = ABI_ALLOC_STAT
       if (istat/=0) then 
         MSG_ERROR('Out-of-memory in Er%epsm1_pole (in-core)')
       end if
       call read_pole_screening(Er%fname,Er%npwe,Er%nqibz,Er%ncoeff,Er%epsm1_pole,accesswff,comm)

       if (present(reconstruct_scr)) then
         if (reconstruct_scr == 1) then
           ! We are supposed to reconstruct the screening file from the
           ! pole-fit for each omega as if it was all read from file.
           ABI_ALLOCATE(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,Er%nqibz))
           istat = ABI_ALLOC_STAT
           if (istat/=0) then 
             MSG_ERROR('Out-of-memory in Er%epsm1 (in-core, reconstruct screening)')
           end if
           do iqibz=1,Er%nqibz
             do ig2=1,Er%npwe
               do ig1=1,Er%npwe
                 call re_and_im_screening_with_phase(Er%omega,Er%epsm1(ig1,ig2,:,iqibz),Er%nomega, &
&                 Er%epsm1_pole(ig1,ig2,:,iqibz),Er%ncoeff)
                 if (ig1==ig2) Er%epsm1(ig1,ig2,:,iqibz) = &
&                 Er%epsm1(ig1,ig2,:,iqibz) + CMPLX(1.0_gwp,0.0_gwp)
               end do
             end do
           end do
           if (associated(Er%epsm1_pole))  then
             ABI_DEALLOCATE(Er%epsm1_pole)
           end if
           Er%fform = 1002
           Er%npoles=0
           Er%ncoeff=0
           msg = "--- NOTE: Reconstructed full screening from pole-fit."
           MSG_COMMENT(msg)
         end if
       end if

  end if
     ! 
   else  ! Out-of-core solution ===
     msg = " mqmem==0 => allocating a single q-slice of (W|chi0) (slower but less memory)."
     MSG_COMMENT(msg)
     continue 
   end if

   RETURN

 else 
  ! === The matrix stored on file do not correspond to the quantity required ===
  ! * Presently only the transformation chi0 => e^-1 is coded
  ! * According to Er%mqmem either calculate e^-1 dumping the result to a file 
  !   for a subsequent use or calculate e^-1 keeping everything in memory.

  if (Er%mqmem==0) then 
   ! === Open file and write the header for the SCR file ===
   ! * For the moment only master works.
   if (my_rank==master) then
    write(msg,'(3a)')ch10,&
&    ' mkdump_Er : calculating and writing epsilon^-1 matrix on file ',TRIM(fname_dump)
    call wrtout(std_out,msg,'COLL')
    unt_dump=get_unit()
    open(unit=unt_dump,file=fname_dump,status='unknown',form='unformatted',iostat=ios)
    if (ios/=0) then 
     write(msg,'(3a)')' Opening file ',TRIM(fname_dump),' as new-unformatted'
     MSG_ERROR(msg)
    end if  

    ! === Update the entries in the header that have been modified ===
    ! TODO, write function to return title, just for info
    ! TODO nullification of pointers in fortran structures is needed to avoid problems during the copy.
    call copy_scrhdr(Er%Hscr,Hscr_cp)
    Hscr_cp%ID        = id_required
    Hscr_cp%ikxc      = ikxc_required
    Hscr_cp%test_type = option_test 
    Hscr_cp%title(1)  = 'SCR file: epsilon^-1'
    Hscr_cp%title(2)  = 'TESTPARTICLE'

    ! Treat the case in which a smaller matrix is used.
    Hscr_cp%npwe      = Er%npwe  
    if (Hscr_cp%nqlwl>0) then ! Reallocate wing with correct dimension keeping previous values. ABI_ALLOCATE(Hscr_cp%nqlwl>0,)
      ABI_ALLOCATE(dummy_lwing,(Er%npwe*Er%nI,Er%nomega,Hscr_cp%nqlwl))
      dummy_lwing(:,:,:) = Hscr_cp%lwing(1:Er%npwe*Er%nI,:,:)
      ABI_DEALLOCATE(Hscr_cp%lwing)
      ABI_ALLOCATE(Hscr_cp%lwing,(Hscr_cp%npwe,Hscr_cp%nomega,Hscr_cp%nqlwl))
      Hscr_cp%lwing = dummy_lwing
      dummy_lwing(:,:,:) = Hscr_cp%uwing(1:Er%npwe*Er%nI,:,:)
      ABI_DEALLOCATE(Hscr_cp%uwing)
      ABI_ALLOCATE(Hscr_cp%uwing,(Hscr_cp%npwe,Hscr_cp%nomega,Hscr_cp%nqlwl))
      Hscr_cp%uwing = dummy_lwing
      ABI_DEALLOCATE(dummy_lwing)
    end if    

    rdwr=2 ; fform=Hscr_cp%fform
    call scr_hdr_io(fform,rdwr,unt_dump,comm_self,master,accesswff,Hscr_cp)
    call free_scrhdr(Hscr_cp)

    !allocate(Er%epsm1(Er%npwe,Er%npwe,Er%nomega,1))
    !epsm1 => Er%epsm1(:,:,:,1)
    ABI_ALLOCATE(epsm1,(Er%npwe,Er%npwe,Er%nomega))

    do iqibz=1,Er%nqibz
     is_qeq0=0
     if (normv(Er%qibz(:,iqibz),gmet,'G')<GW_TOLQ0) is_qeq0=1

     call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,epsm1,&
&                        accesswff,comm_self,iqiA=iqibz)

     dim_wing=0; if (is_qeq0==1) dim_wing=3
     ABI_ALLOCATE(dummy_lwing,(Er%npwe*Er%nI,Er%nomega,dim_wing))
     ABI_ALLOCATE(dummy_uwing,(Er%npwe*Er%nJ,Er%nomega,dim_wing))
     ABI_ALLOCATE(dummy_head,(dim_wing,dim_wing,Er%nomega))

     if (approx_type<2) then
       MSG_WARNING('Entering out-of core RPA or Kxc branch')
       call make_epsm1_driver(iqibz,dim_wing,Er%npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&                approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&                dummy_lwing,dummy_uwing,epsm1,Spectra,comm_self)
     else
       MSG_WARNING('Entering out-of core fxc_ADA branch')
       call make_epsm1_driver(iqibz,dim_wing,Er%npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&                approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&                dummy_lwing,dummy_uwing,epsm1,Spectra,comm_self,fxc_ADA(:,:,iqibz))
     end if

     ABI_DEALLOCATE(dummy_head)
     ABI_DEALLOCATE(dummy_uwing)
     ABI_DEALLOCATE(dummy_lwing)

     if (is_qeq0==1) then
      call repr_dielconst(Spectra,msg)
      call wrtout(std_out,msg,'COLL') 
      call wrtout(ab_out,msg,'COLL')
     end if

     call destroy_spectra(Spectra)

     call write_screening(unt_dump,accesswff,Er%npwe,Er%nomega,epsm1)
    end do

    close(unt_dump)
    ABI_DEALLOCATE(epsm1)
   end if !master
 
   ! A synchronization is required here, else the other procs start to read the
   ! file _SCR before it is written by the master
   call xbarrier_mpi(comm)

   ! Now Er% "belongs" to the file "fname_dump", thus 
   ! each proc has to destroy and re-initialize the object.
   call destroy_Epsilonm1_results(Er)

   mqmem_=Er%mqmem; npwe_asked=Er%npwe
   call init_Er_from_file(Er,fname_dump,mqmem_,npwe_asked,accesswff,comm)

   !Now Er% has been reinitialized and ready-to-use.
   Er%ID=id_required
   call print_epsilonm1_results(Er)

  else 
   ! ========================
   ! === In-core solution ===
   ! ========================
   ABI_ALLOCATE(Er%lwing,(Er%npwe,Er%nomega,Er%nqlwl))
   ABI_ALLOCATE(Er%uwing,(Er%npwe,Er%nomega,Er%nqlwl))
   ABI_ALLOCATE(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,Er%nqibz))
   istat = ABI_ALLOC_STAT
   if (istat/=0) then 
     MSG_ERROR(' Out-of-memory in Er%epsm1 (in-core)')
   end if

   call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,accesswff,comm)

   do iqibz=1,Er%nqibz
    is_qeq0=0
    if (normv(Er%qibz(:,iqibz),gmet,'G')<GW_TOLQ0) is_qeq0=1
    epsm1 => Er%epsm1(:,:,:,iqibz)
                                                                                                           
    dim_wing=0; if (is_qeq0==1) dim_wing=3 ! FIXME
    ABI_ALLOCATE(dummy_lwing,(Er%npwe*Er%nI,Er%nomega,dim_wing))
    ABI_ALLOCATE(dummy_uwing,(Er%npwe*Er%nJ,Er%nomega,dim_wing))
    ABI_ALLOCATE(dummy_head,(dim_wing,dim_wing,Er%nomega))

    if (approx_type<2) then
      MSG_WARNING('Entering in-core RPA and Kxc branch')
      call make_epsm1_driver(iqibz,dim_wing,Er%npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&               approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&               dummy_lwing,dummy_uwing,epsm1,Spectra,comm)
    else
      MSG_WARNING('Entering in-core fxc_ADA branch')
      call make_epsm1_driver(iqibz,dim_wing,Er%npwe,Er%nI,Er%nJ,Er%nomega,Er%omega,&
&               approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,dummy_head,&
&               dummy_lwing,dummy_uwing,epsm1,Spectra,comm,fxc_ADA=fxc_ADA(:,:,iqibz))
    end if

    ABI_DEALLOCATE(dummy_lwing)
    ABI_DEALLOCATE(dummy_uwing)
    ABI_DEALLOCATE(dummy_head)

    if (is_qeq0==1) then
      call repr_dielconst(Spectra,msg)
      call wrtout(std_out,msg,'COLL') 
      call wrtout(ab_out,msg,'COLL')
    end if

    call destroy_spectra(Spectra)
   end do

   Er%ID=id_required
   call print_epsilonm1_results(Er)

  end if

 end if

 DBG_EXIT("COLL")

end subroutine mkdump_Er
!!***

!----------------------------------------------------------------------

!!****f* m_screening/get_epsm1
!! NAME
!!  get_epsm1
!!
!! FUNCTION
!!  Work in progress but the main is idea is as follows:
!!
!!  Return the symmetrized inverse dielectric matrix.
!!  This method implements both in-core and the out-of-core solution 
!!  In the later, epsilon^-1 or chi0 are read from file.
!!  It is possible to specify options to retrieve (RPA |TDDDT, [TESTCHARGE|TESTPARTICLE]).
!!  All dimensions are already initialized in the Er% object, this method 
!!  should act as a wrapper around rdscr and make_epsm1_driver. A better 
!!  implementation will be done in the following once the coding of file handlers is completed.
!!
!! INPUTS
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  iqibzA[optional]=Index of the q-point to be read from file (only for out-of-memory solutions)
!!  accesswff=option definig the file format.
!!  option_test
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Er%epsm1
!!
!! TODO
!!  Remove this routine. Now everything should be done with mkdump_Er
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_epsm1(Er,Vcp,approx_type,option_test,accesswff,comm,iqibzA)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_epsm1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff,option_test,approx_type,comm
 integer,optional,intent(in) :: iqibzA
 type(vcoul_t),intent(in) :: Vcp
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: istat,my_approx_type,my_option_test,ng
 character(len=500) :: msg                   

! *********************************************************************

 DBG_ENTER("COLL")

 my_approx_type = approx_type
 my_option_test = option_test

 ! Vcp not yet used.
 ng = Vcp%ng        

 select case (Er%mqmem)

 case (0) !  Out-of-core solution
   if (associated(Er%lwing))  then
     ABI_DEALLOCATE(Er%lwing)
   end if
   if (associated(Er%uwing))  then
     ABI_DEALLOCATE(Er%uwing)
   end if
   if (associated(Er%epsm1))  then
     ABI_DEALLOCATE(Er%epsm1)
   end if
   ABI_ALLOCATE(Er%lwing,(Er%npwe,Er%nomega,Er%nqlwl))
   ABI_ALLOCATE(Er%uwing,(Er%npwe,Er%nomega,Er%nqlwl))
   ABI_ALLOCATE(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,1))
   istat = ABI_ALLOC_STAT
   if (istat/=0) then 
     msg = ' Out-of-memory in Er%epsm1 (out-of-core)'
     MSG_ERROR(msg)
   end if

   call read_screening(Er%fname,Er%npwe,Er%nqibz,Er%nomega,Er%epsm1,accesswff,comm,iqiA=iqibzA)

   if (Er%ID==4) then 
     ! === If q-slice of epsilon^-1 has been read then return === 
     !call print_epsilonm1_results(Er)
     RETURN 
   else 
     MSG_ERROR('Wrong Er%ID')
   end if

 case default
   ! ========================
   ! === In-core solution ===
   ! ========================
   MSG_ERROR("you should not be here")
 end select

 DBG_EXIT("COLL")

end subroutine get_epsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/get_pole_epsm1
!! NAME
!!  get_epsm1
!!
!! FUNCTION
!!  Work in progress but the main is idea is as follows:
!!
!!  Return the symmetrized inverse dielectric matrix.
!!  This method implements both in-core and the out-of-core solution 
!!  In the later, epsilon^-1 or chi0 are read from file.
!!  It is possible to specify options to retrieve (RPA |TDDDT, [TESTCHARGE|TESTPARTICLE]).
!!  All dimensions are already initialized in the Er% object, this method 
!!  should act as a wrapper around rdscr and make_epsm1_driver. A better 
!!  implementation will be done in the following once the coding of file handlers is completed.
!!
!! INPUTS
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  iqibzA[optional]=Index of the q-point to be read from file (only for out-of-memory solutions)
!!  accesswff=option definig the file format.
!!  option_test
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Er%epsm1
!!
!! TODO
!!  Remove this routine. Now everything should be done with mkdump_Er
!!
!! PARENTS
!!      calc_sigc_me,cohsex_me
!!
!! CHILDREN
!!      copy_scrhdr,deep_copy,destroy_bz_mesh_type,destroy_crystal
!!      destroy_gsphere,find_qmesh,free_scrhdr,get_bz_item
!!      init_crystal_from_hdr,init_gsphere,init_kmesh,print_bz_mesh
!!      read_screening,scr_hdr_io,sort_dp,wrap2_pmhalf,write_screening,wrtout
!!
!! SOURCE

subroutine get_pole_epsm1(Er,Vcp,approx_type,option_test,accesswff,comm,iqibzA,reconstruct_scr)

 use defs_basis
 use m_io_screening, only : read_pole_screening
 use m_model_screening, only : re_and_im_screening_with_phase

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_pole_epsm1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff,option_test,approx_type,comm
 integer,optional,intent(in) :: iqibzA,reconstruct_scr
 type(vcoul_t),intent(in) :: Vcp
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: istat,my_approx_type,my_option_test,ng,ig1,ig2
 character(len=500) :: msg                   

! *********************************************************************

 DBG_ENTER("COLL")

 my_approx_type = approx_type
 my_option_test = option_test

 ! Vcp not yet used.
 ng = Vcp%ng        

 select case (Er%mqmem)

 case (0) !  Out-of-core solution
   if (associated(Er%lwing))  then
     ABI_DEALLOCATE(Er%lwing)
   end if
   if (associated(Er%uwing))  then
     ABI_DEALLOCATE(Er%uwing)
   end if
   if (associated(Er%epsm1))  then
     ABI_DEALLOCATE(Er%epsm1)
   end if
   if (associated(Er%epsm1_pole))  then
     ABI_DEALLOCATE(Er%epsm1_pole)
   end if
   ABI_ALLOCATE(Er%lwing,(Er%npwe,Er%nomega,Er%nqlwl))
   ABI_ALLOCATE(Er%uwing,(Er%npwe,Er%nomega,Er%nqlwl))
   ABI_ALLOCATE(Er%epsm1_pole,(Er%npwe,Er%npwe,Er%ncoeff,1))
   istat = ABI_ALLOC_STAT
   if (istat/=0) then 
     msg = ' Out-of-memory in Er%epsm1_pole (out-of-core)'
     MSG_ERROR(msg)
   end if

   call read_pole_screening(Er%fname,Er%npwe,Er%nqibz,Er%ncoeff,Er%epsm1_pole,accesswff,comm,iqiA=iqibzA)

   if (present(reconstruct_scr)) then
     if (reconstruct_scr==1) then
       ABI_ALLOCATE(Er%epsm1,(Er%npwe,Er%npwe,Er%nomega,1))
       istat = ABI_ALLOC_STAT
       if (istat/=0) then 
         MSG_ERROR('Out-of-memory in Er%epsm1 (out-of-core, reconstruct screening)')
       end if
       do ig2=1,Er%npwe
         do ig1=1,Er%npwe
           call re_and_im_screening_with_phase(Er%omega,Er%epsm1(ig1,ig2,:,1),Er%nomega, &
&           Er%epsm1_pole(ig1,ig2,:,1),Er%ncoeff)
           if (ig1==ig2) Er%epsm1(ig1,ig2,:,1) = &
&            Er%epsm1(ig1,ig2,:,1) + CMPLX(1.0_gwp,0.0_gwp)
         end do
       end do
       if (associated(Er%epsm1_pole))  then
         ABI_DEALLOCATE(Er%epsm1_pole)
       end if
       msg = "--- NOTE: Reconstructed full screening from pole-fit. (out-of-core)"
       MSG_COMMENT(msg)
     end if
   end if

   if (Er%ID==4) then 
     ! === If q-slice of epsilon^-1 has been read then return === 
     !call print_epsilonm1_results(Er)
     RETURN 
   else 
     MSG_ERROR('Wrong Er%ID')
   end if

 case default
   ! ========================
   ! === In-core solution ===
   ! ========================
   MSG_ERROR("you should not be here")
 end select

 DBG_EXIT("COLL")

end subroutine get_pole_epsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/decompose_epsm1
!! NAME
!! decompose_epsm1
!!
!! FUNCTION
!! Decompose the complex symmetrized dielectric 
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine decompose_epsm1(Er,iqibz,eigenvalues)

 use defs_basis
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'decompose_epsm1'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz
 type(Epsilonm1_results),intent(in) :: Er
!arrays
 complex(dpc),intent(out) :: eigenvalues(Er%npwe,Er%nomega)

!Local variables-------------------------------
!scalars
 integer :: info,lwork,iomega,istat,negw,ig1,ig2,idx,sdim
 character(len=500) :: msg                   
!arrays
 real(dp),allocatable :: ww(:),rwork(:)
 complex(dpc),allocatable :: work(:),Adpp(:),eigvec(:,:),Afull(:,:),vs(:,:),wwc(:)
 logical,allocatable :: bwork(:)
 logical :: sortcplx !BUG in abilint

! *********************************************************************

 ABI_CHECK(Er%mqmem/=0,'mqmem==0 not implemented')

 do iomega=1,Er%nomega

   if (ABS(REAL(Er%omega(iomega)))>0.00001) then ! Eigenvalues for a generic complex matrix ===
     !if (.TRUE.) then
     lwork=4*2*Er%npwe
     ABI_ALLOCATE(wwc,(Er%npwe))
     ABI_ALLOCATE(work,(lwork))
     ABI_ALLOCATE(rwork,(Er%npwe))
     ABI_ALLOCATE(bwork,(Er%npwe))
     ABI_ALLOCATE(vs,(Er%npwe,Er%npwe))
     istat = ABI_ALLOC_STAT
     ABI_ALLOCATE(Afull,(Er%npwe,Er%npwe))
     istat = ABI_ALLOC_STAT

     Afull=Er%epsm1(:,:,iomega,iqibz)

     !for the moment no sort, maybe here I should sort using the real part?
     call ZGEES('V','N',sortcplx,Er%npwe,Afull,Er%npwe,sdim,wwc,vs,Er%npwe,work,lwork,rwork,bwork,info)
     if (info/=0) then 
       write(msg,'(2a,i10)')' decompose_epsm1 : Error in ZGEES, diagonalizing complex matrix, info = ',info
       call wrtout(std_out,msg,'COLL') 
     end if

     eigenvalues(:,iomega)=wwc(:)

     ABI_DEALLOCATE(wwc)
     ABI_DEALLOCATE(work)
     ABI_DEALLOCATE(rwork)
     ABI_DEALLOCATE(bwork)
     ABI_DEALLOCATE(vs)
     ABI_DEALLOCATE(Afull)

   else ! Hermitian version.

     lwork=2*Er%npwe-1
     ABI_ALLOCATE(ww,(Er%npwe))
     ABI_ALLOCATE(work,(lwork))
     ABI_ALLOCATE(rwork,(3*Er%npwe-2))
     ABI_ALLOCATE(eigvec,(Er%npwe,Er%npwe))
     ABI_ALLOCATE(Adpp,(Er%npwe*(Er%npwe+1)/2))
     istat = ABI_ALLOC_STAT
     if (istat/=0) STOP ' decompose_epsm1 : out of memory in Adpp'

     idx=0 ! Pack the matrix
     do ig2=1,Er%npwe
       do ig1=1,ig2
         idx=idx+1
         Adpp(idx)=Er%epsm1(ig1,ig2,iomega,iqibz)
       end do 
     end do

     ! For the moment we require also the eigenvectors.
     call ZHPEV('V','U',Er%npwe,Adpp,ww,eigvec,Er%npwe,work,rwork,info)
     if (info/=0) then 
       write(msg,'(a,i5)')' decompose_epsm1 : Error in ZHPEV, diagonalizing matrix, info = ',info
       MSG_ERROR(msg)
     end if

     negw=(COUNT((REAL(ww)<tol6)))
     if (negw/=0) then 
       write(msg,'(a,i5,a,i3,a,f8.4)')&
&        ' Found negative eigenvalues. No. ',negw,' at iqibz= ',iqibz,' minval= ',MINVAL(REAL(ww))
       MSG_WARNING(msg)
     end if

     eigenvalues(:,iomega)=ww(:)

     ABI_DEALLOCATE(ww)
     ABI_DEALLOCATE(work)
     ABI_DEALLOCATE(rwork)
     ABI_DEALLOCATE(eigvec)
     ABI_DEALLOCATE(Adpp)
   end if
 end do !iomega

! contains 
! function sortcplx(carg) result(res)
!  implicit none 
!  complex(dpc),intent(in) :: carg
!  logical :: res
!  res=.TRUE.
! end function sortcplx

end subroutine decompose_epsm1
!!***

!----------------------------------------------------------------------

!!****f* m_screening/outeps
!! NAME
!! outeps
!!
!! FUNCTION
!!  Write the independent matrix elements of epsilon^{-1}_q(G1,G2)(omega) on file
!!
!! INPUTS
!!  Gsphere<Gvectors_type>
!!   %ng=number of G vectors in the sphere
!!   %gvec(3,ng)= G vectors in reduced coordinates
!!   %gmet(3,3)=metric tensor in reciprocal space
!!  Gpairs_q<Gpairs_type>
!!   %ngpi=nuber of independent (G1,G1) pairs   
!!   %ip2fp(2,ngpi)= index of G1 and G2 in the array gvec for each ngpi independent pair
!!  nomega=number of frequencies
!!  title=title describing what is printed
!!  nge=Number of G vectors in epsilon
!!  eps(nge,nge,nomega)=matrix to be printed
!!  omega(nomega)=frequencies
!!  prtvol= if different from 0 first and last 150 independent pairs are written 
!!          if ==0 only first and last 50 pairs are written
!!  unt=unit number of external file
!!
!! OUTPUT
!!  Only writing 
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine outeps(nge,nomega,omega,eps,Gsphere,Gpairs_q,title,unt,prtvol)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outeps'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nge,nomega,prtvol,unt
 character(len=500),intent(in) :: title
 type(Gpairs_type),intent(in) :: Gpairs_q
 type(Gvectors_type),intent(in) :: Gsphere
!arrays
 complex(gwpc),intent(in) :: eps(nge,nge,nomega)
 complex(dpc),intent(in) :: omega(nomega)

!Local variables-------------------------------
!scalars
 integer :: enough,f1,f2,ii,io,ip,ng,ngpi,npmax
 real(dp) :: e1,e2
 character(len=500) :: msg
!arrays
 integer,pointer :: gvec(:,:)
 real(dp) :: gmet(3,3)

!************************************************************************

 ng   =  Gsphere%ng 
 gmet =  Gsphere%gmet
 gvec => Gsphere%gvec(1:3,1:ng)

 ngpi =  Gpairs_q%niggp

 npmax=50 ; if (ngpi<50) npmax=ngpi
 if (prtvol/=0) then 
   npmax=150; if (ngpi<150) npmax=ngpi
 end if

 if (Gpairs_q%niggp<Gpairs_q%ng**2) then 
   !
   ! === For this q-point we used al least 2 symmetries, tables in Gpairs are associated === 
   write(msg,'(3a)')ch10,'# Independent matrix elements of: ',TRIM(title)
   call wrtout(unt,msg,'COLL')

   do io=1,nomega
     ! === Write header ===
     write(msg,'(2a,i3,a,2f9.4,3a)')ch10,&
&     '#',io,'-th omega :',omega(io)*Ha_eV,' [eV] ',ch10,&
&     '#         G1         E1 [H]           G2         E2 [H]         Re        Im '
     call wrtout(unt,msg,'COLL')

     enough=0
     do ip=1,ngpi
       enough=enough+1
       !if (enough>npmax.and.enough<ngpi-npmax.and.prtvol==0) CYCLE
       if (enough>npmax.and.enough<ngpi-npmax) CYCLE
       f1=Gpairs_q%ip2fp(1,ip)
       f2=Gpairs_q%ip2fp(2,ip)
       e1=normv(gvec(:,f1),gmet,'G') ; e1=half*e1**2
       e2=normv(gvec(:,f2),gmet,'G') ; e2=half*e2**2

       write(unt,'(2(3i6,4x,f5.2,2x),4x,2(f8.4,2x))')gvec(:,f1),e1,gvec(:,f2),e2,eps(f1,f2,io)
       if (enough==npmax.and.prtvol==0) then 
         write(unt,*)('-',ii=1,77)
         write(unt,*)' outeps: prtvol=0, stop printing more G1 G2 information' 
         write(unt,*)('-',ii=1,77)
       end if 
     end do !ip
   end do !io

 else 
   !
   ! === We do not have irreducible pairs ===
   write(msg,'(2a,i4,2a)')ch10, '# First and last ',npmax,' matrix elements of ',TRIM(title)
   call wrtout(unt,msg,'COLL')

   do io=1,nomega
     ! === Write header ===
     write(msg,'(2a,i3,a,2f9.4,3a)')ch10,&
&      '#',io,'-th omega :',omega(io)*Ha_eV,' [eV] ',ch10,&
&      '#         G1         E1 [H]           G2         E2 [H]         Re        Im '
     call wrtout(unt,msg,'COLL')

     enough=0
     do f1=1,ng
       e1=normv(gvec(:,f1),gmet,'G') ; e1=half*e1**2
       do f2=1,ng
         enough=enough+1
         !if (enough>npmax.and.enough<ngpi-npmax.and.prtvol==0) CYCLE
         if (enough>npmax.and.enough<ngpi-npmax) CYCLE
         e2=normv(gvec(:,f2),gmet,'G') ; e2=half*e2**2
         write(unt,'(2(3i6,4x,f5.2,2x),4x,2(f8.4,2x))')gvec(:,f1),e1,gvec(:,f2),e2,eps(f1,f2,io)
         if (enough==npmax.and.prtvol==0) then 
           write(unt,*)('-',ii=1,77)
           write(unt,*)' outeps : prtvol=0, stop printing more G-G'' information' 
           write(unt,*)('-',ii=1,77)
         end if 
       end do !f2
     end do !f1
  end do !io

 end if

end subroutine outeps
!!***

!----------------------------------------------------------------------

!!****f* m_screening/symf12
!! NAME
!! symf12
!!
!! FUNCTION 
!!  This subroutine can be used in two different modes:
!!
!!  Mode 1) 
!!   Reconstruct all the reciprocal space matrix elements of a two point function F_{G1,G2}(q,omega) 
!!   for a fixed q-point, using only the set of independent (G1,G2) pairs.
!!   The independent pairs reconstruct the full G \times G space when all the operations of the little group of q are applied. 
!!  Note that F(r1,r2) is supposed to have the same symmetry as the crystal i.e 
!!   F(R^{-1}(r1-t),R^{-1}(r2-t)) = f(r1,r2)  where t is a fractional translation
!!
!!  Mode 2) 
!!   Check if a full matrix satisfies the above mentioned symmetry properties 
!!   (mainly used to debug the code)
!!
!! INPUTS
!!  nsym=number of symmetry operations                                                                             
!!  npwc=number of G vectors in the matrix F_{G1,G2}(q)
!!  symrec(3,3,nsym)= symmetry operations in reciprocal space
!!  Gpairs_q<Gpairs_type>= tables related to the irreducible G1,G2 pairs for this q-point
!!   %ngpi=nuber of independent (G1,G1) pairs   
!!   %ng=Number of G vectors
!!   %ip2fp(2,ngpi)=index of G1 and G2 in the array gvec for each ngpi independent pair 
!!   %fptabo(ng,ng)=index of the symmetry operation S in the array symrec such as (T1,T2)= S (G1,G2)
!!  phgt(npw,nsym)=phase factors associated to non simmorphic operations (e^{-iG \cdot t})
!!
!! OUTPUT
!!  eps(npw,npw,nomega), see side effects
!!
!! SIDE EFFECTS
!!  Mode 1) 
!!   in input: eps contains only the matrix elements corresponding to the independent pair 
!!   in output: full symmetrized eps 
!!  Mode 2) 
!!   Only checking
!!   eps matrix is supposed to contain all the elements, 
!!
!! TODO 
!!  One can take advantage of hermitianess in case of imaginary frequencies
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine symf12(npw,nomega,eps,Gpairs_q,nsym,symrec,phgt,mode)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symf12'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mode,nomega,npw,nsym
 type(Gpairs_type),intent(in) :: Gpairs_q
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 complex(dpc),intent(in) :: phgt(npw,nsym) 
 complex(gwpc),intent(inout) :: eps(npw,npw,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig1,ig2,iid,io,ioimax,iormax,ip,ir1,ir2,isym,itim
 real(dp),parameter :: tol=10d-1
 complex :: diff,phase
 logical :: found_identity
 character(len=500) :: msg
!arrays
 integer :: identity(3,3),iimax(2),irmax(2)
 real(dp) :: den(2),err(2),maxerr(2)

! *************************************************************************

 if (mode/=1.and.mode/=2) then
   write(msg,'(a,i4)')' mode should be 1 or 2 however, mode= ',mode
   MSG_BUG(msg)
 end if
 
 identity(:,:)=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))

 found_identity=.FALSE.
 do isym=1,nsym
   if (ALL(symrec(:,:,isym)==identity)) then 
    iid=isym; found_identity=.TRUE.; EXIT
   end if
 end do 
 
 if (.not.found_identity) then 
   write(msg,'(3a)')&
&    ' Only the inversion was found in the set of symmetries read from the KSS file ',ch10,&
&    ' Likely old version of KSS file '
   MSG_ERROR(msg)
 else 
   write(msg,'(a,i4)')' found identity with index : ',iid
   call wrtout(std_out,msg,'COLL')
 end if 

 ! === The basic property used is ==
 ! If q=Sq then 
 ! $ F_{SG1,SG2}(q)=e^{+iS(G2-G1)\cdot\tau} F_{G1,G2)}(q)$
 ! Besides only at Gamma 
 ! $ F_{-SG1,-SG2)(0)=e^{-iS(G2-G1)\cdot\tau} F_{-G1,-G2}^*(0)$
 !
 SELECT CASE (mode)

 CASE (1) ! symmetrize Matrix.
   do ig1=1,npw
     do ig2=1,npw 

      isym=Gpairs_q%fptabo(ig1,ig2) ; itim=1  
      if (isym<0) then  ! time-reversal only at Gamma, in this case isym is negative     
        isym=-isym; itim=2
        MSG_ERROR("check equation")
      end if

      if (isym==iid.and.itim==1) CYCLE  ! This pair is irreducible
      ir1=Gpairs_q%fp2ip(1,ig1,ig2)     ! Indeces in the array gvec of the irred pair
      ir2=Gpairs_q%fp2ip(2,ig1,ig2)             
      ! TODO can speed up by calculating the phase before entering the loop, but require more memory
      phase=phgt(ig1,isym)*CONJG(phgt(ig2,isym)) 

      do io=1,nomega
        ! might also take into account hermiticity in case of imaginary frequencies
        eps(ig1,ig2,io)=eps(ir1,ir2,io)*phase
      end do 
      if (itim==2) eps(ig1,ig2,:)=CONJG(eps(ig1,ig2,:))
      
     end do 
   end do 

 CASE (2) ! Check symmetries.
   write(msg,'(2a,f5.2,a)')ch10,' symf12 : checking symmetry properties, relative tolerance : ',tol,ch10
   call wrtout(std_out,msg,'COLL')

   !if (ALL(Gpairs_q%fpptabo==iid)) then 
   ! !  A non symmorphyc translation associated to the identity operator means that the 
   ! !  unit cell is not primitive (I think the code should realize and stop)
   ! write(msg,'(2a)')' symff12 : all the pairs are independent, return',ch10
   ! call wrtout(std_out,msg,'COLL')
   ! RETURN
   !end if 

   err=zero ; maxerr(:)=zero

   irmax=0 ;  iimax=0
   iormax=0; ioimax=0

   write(msg,'(a)')' val1 ... val2  percent err'
   call wrtout(std_out,msg,'COLL')

   do ig1=1,npw
     do ig2=1,npw 
       isym=Gpairs_q%fptabo(ig1,ig2) ; itim=1 ! time-reversal only at Gamma, in this case isym is negative
       if (isym<0) then       
         isym=-isym;itim=2
         MSG_ERROR("check equation")
       end if

       if (isym==iid.and.itim==1) CYCLE  ! This pair is irreducible
       ir1=Gpairs_q%fp2ip(1,ig1,ig2)     ! Indeces in the array gvec of the irred pair 
       ir2=Gpairs_q%fp2ip(2,ig1,ig2)        
       phase=phgt(ig1,isym)*CONJG(phgt(ig2,isym))

       do io=1,nomega
         if (itim==1) then 
           diff=eps(ig1,ig2,io)-eps(ir1,ir2,io)*phase
         else 
           diff=eps(ig1,ig2,io)-CONJG(eps(ir1,ir2,io)*phase)
         end if 

         err(1)=ABS(REAL(diff)) ; den(1)=MAX(ABS(REAL(eps(ig1,ig2,io),dp)),tol12)
         err(2)=ABS(AIMAG(diff)); den(2)=MAX(ABS(REAL(AIMAG(eps(ig1,ig2,io)),dp)),tol12)

         if (den(1)*tol<err(1)) then 
           write(std_out,'(3(es18.6,4x),a,es18.6,a,6i4)')&
&            REAL(eps(ig1,ig2,io)),zero,REAL(eps(ir1,ir2,io)*phase),' > ',100*err(1)/den(1),' %',ig1,ig2,ir1,ir2,isym,io
           if (err(1)>maxerr(1)) then 
             maxerr(1)=err(1)
             irmax(1)=ig1
             irmax(2)=ig2
             iormax=io
           end if
         end if 
         
         if (den(2)*tol<err(2)) then 
           write(std_out,'(3(es18.6,4x),a,es18.6,a,6i4)')&
&            zero,AIMAG(eps(ig1,ig2,io)),AIMAG(eps(ir1,ir2,io)*phase),' > ',100*err(2)/den(2),' %',ig1,ig2,ir1,ir2,isym,io
           if (err(2)>maxerr(2)) then 
             maxerr(2)=err(2)
             iimax(1)=ig1
             iimax(2)=ig2
             ioimax=io
           end if
         end if 

       end do !iomega
     end do !ig2 
   end do !ig1

   if (iormax/=0) then 
     write(std_out,'(a,2(2x,f24.6))')' symf12 : max abs error for real part = ',maxerr(1)
     ig1=irmax(1) ; ig2=irmax(2)
     isym=Gpairs_q%fptabo(ig1,ig2)
     ir1=Gpairs_q%ip2fp(1,ip) ; ir2=Gpairs_q%ip2fp(2,ip)
     write(std_out,*)'ig1,ig2,ir1,ir2,iormax=',ig1,ig2,ir1,ir2,iormax
     write(std_out,'(2(2x,es18.6))')&
&      REAL(eps(ig1,ig2,iormax)),REAL(eps(ir1,ir2,iormax)*phase)
   end if 
   if (ioimax/=0) then
     write(std_out,'(a,2(2x,f24.6))')' symf12 : max abs error for imag part = ',maxerr(2)
     ig1=iimax(1) ; ig2=iimax(2)
     isym=Gpairs_q%fptabo(ig1,ig2)
     ir1=Gpairs_q%ip2fp(1,ip) ; ir2=Gpairs_q%ip2fp(2,ip)
     write(std_out,*)'ig1,ig2,ir1,ir2,iormax=',ig1,ig2,ir1,ir2,ioimax
     write(std_out,'(2(2x,es18.6))')AIMAG(eps(ig1,ig2,ioimax)),AIMAG(eps(ir1,ir2,ioimax)*phase)
   end if 

 CASE DEFAULT
  write(msg,'(a,i4)')'Wrong value for mode= ',mode
  MSG_BUG(msg)
 END SELECT

end subroutine symf12
!!***

!----------------------------------------------------------------------

!!****f* m_screening/make_epsm1_driver
!! NAME
!! make_epsm1_driver
!!
!! FUNCTION
!!  Driver routine to calculate the inverse symmetrical dielectric matrix starting
!!  from the irreducible polarizability. The routine considers a single q-point, and 
!!  performs the following tasks:
!!
!!  1) Calculate $\tilde\epsilon^{-1}$ using different approximations:
!!      * RPA
!!      * ALDA within TDDFT
!!
!!  2) Use a special treatment of non-Analytic behavior of heads and wings in reciprocal space
!!     calculating these quantities for different small q-directions specified by the user
!!     (Not yet operative)
!!
!!  3) Output the electron energy loss function and the macroscopic dielectric function with and 
!!     without local field effects (only if non-zero real frequencies are available)
!!
!! INPUTS
!!  iqibz=index of the q-point in the array Vcp%qibz where epsilon^-1 has to be calculated
!!  npwe=Number of G-vectors in chi0.
!!  nI,nJ=Number of rows/columns in chi0_ij (1,1 in collinear case)
!!  nomega=Number of frequencies.
!!  dim_wing=Dimension of the wings (0 or 3 if q-->0)
!!  approx_type=Integer flag defining the type of approximation
!!   == 0 for RPA   ==
!!   == 1 for TDDFT ==
!!  option_test=Only for TDDFT:
!!   == 0 for TESTPARTICLE ==
!!   == 1 for TESTELECTRON ==
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!   %nqibz=Number of q-points.
!!   %qibz(3,nqibz)=q-points in the IBZ.  
!!  nkxc=Integer defining the dimension of the kernel in reciprocal space
!!  kxcg(nfftot,nkxc)=TDDFT kernel in reciprocal space on the FFT mesh. Needed only if approx_type==1
!!  comm=MPI communicator.
!!  ngfft(18)=Info on the FFT mesh.
!!  nfftot=Total number of points in the FFT mesh.
!!  chi0_lwing(npwe*nI,nomega,dim_wing)=Lower wings of chi0 (only for q-->0)
!!  chi0_uwing(npwe*nJ,nomega,dim_wing)=Upper wings of chi0 (only for q-->0)
!!  chi0_head(dim_wing,dim_wing,nomega)=Head of of chi0 (only for q-->0)
!!
!! OUTPUT
!!  Different files are written according to the type of calculation
!!  See also side effects
!!
!! SIDE EFFECTS
!!  chi0(npwe*nI,npwe*nJ,nomega): in input the irreducible polarizability, in output 
!!   the symmetrized inverse dielectric matrix.
!!
!! PARENTS
!!      m_screen,m_screening,screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_epsm1_driver(iqibz,dim_wing,npwe,nI,nJ,nomega,omega,&
& approx_type,option_test,Vcp,nfftot,ngfft,nkxc,kxcg,gvec,chi0_head,&
& chi0_lwing,chi0_uwing,chi0,Spectra,comm,fxc_ADA)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_epsm1_driver'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nI,nJ,npwe,nomega,dim_wing,approx_type,option_test,nkxc,nfftot,comm
 type(vcoul_t),intent(in) :: Vcp
 type(Spectra_type),intent(out) :: Spectra
!arrays
 integer,intent(in) :: ngfft(18),gvec(3,npwe)
 complex(gwpc),intent(in) :: kxcg(nfftot,nkxc)
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(inout) :: chi0_lwing(npwe*nI,nomega,dim_wing)
 complex(dpc),intent(inout) :: chi0_uwing(npwe*nJ,nomega,dim_wing)
 complex(dpc),intent(inout) :: chi0_head(dim_wing,dim_wing,nomega)
 complex(gwpc),intent(inout) :: chi0(npwe*nI,npwe*nJ,nomega)
 complex(gwpc),intent(in), optional :: fxc_ADA(npwe*nI,npwe*nJ)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ig1,ig2,io,istat,ierr,irank,master,iqlwl,my_nqlwl
 integer :: nor,my_rank,nprocs,use_MPI,comm_self,g1mg2_idx
 real(dp) :: ucvol
 logical :: is_qeq0
 character(len=500) :: msg
!arrays
 integer,allocatable :: omega_distrb(:)
 integer,allocatable :: istart(:),istop(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: eelf(:,:) 
 complex(gwpc),allocatable :: chitmp(:,:) 
 complex(dpc),allocatable :: epsm_lf(:,:),epsm_nlf(:,:)
 complex(gwpc),pointer :: vc_sqrt(:)
 complex(gwpc),allocatable :: chi0_save(:,:),kxcg_mat(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 !call lebedev_laikov_int()  

 if (nI/=1.or.nJ/=1) then 
   MSG_ERROR("nI or nJ=/1 not yet implemented")
 end if

 if (.FALSE.) write(std_out,*)chi0_head(1,1,1)

 nprocs  = xcomm_size(comm)
 my_rank = xcomm_rank(comm)
 master=0

 ! MG TODO We use comm_self for the inversion as the single precision version is not yet available
 comm_self = xmpi_self

 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 is_qeq0 = (normv(Vcp%qibz(:,iqibz),gmet,'G')<GW_TOLQ0)

 use_MPI = 0  ! Parallelism is not yet used 
 !if (nprocs>=nomega) use_MPI = 1

 if (use_MPI==1) then ! * Initialize distribution table for frequencies.
   ABI_ALLOCATE(istart,(nprocs))
   ABI_ALLOCATE(istop,(nprocs))
   call xmpi_split_work2_i4b(nomega,nprocs,istart,istop,msg,ierr)
   if (ierr/=0) then
     MSG_WARNING(msg)
   end if 
   omega_distrb(:)=-999
   do irank=0,nprocs-1
     i1 = istart(irank+1)
     i2 = istop (irank+1)
     if (i1<=i2) omega_distrb(i1:i2) = irank
   end do
   ABI_DEALLOCATE(istart)
   ABI_DEALLOCATE(istop)
 end if

 ! * Initialize container for spectral results
 do nor=1,nomega
   if (ABS(AIMAG(omega(nor)))>1.e-3) EXIT
 end do
 nor = nor -1 ; if (nor==0) nor = 1 ! only imag !?

 if (dim_wing==3) then 
   call wrtout(std_out,' Analyzing long wavelength limit for several q','COLL')
   call init_spectra(Spectra,nor,REAL(omega(1:nor)),Vcp%nqlwl,Vcp%qlwl)
   my_nqlwl = 1
   !my_nqlwl = dim_wing ! TODO
   !ABI_CHECK(dim_wing==SIZE(Vcp%vcqlwl_sqrt,DIM=2),"WRONG DIMS")
 else
   call init_spectra(Spectra,nor,REAL(omega(1:nor)),1,Vcp%qibz(:,iqibz))
   my_nqlwl = 1
 end if

!FB: all processors have to perform this operation in order to have the epsm1 matrix when performing
!    a sigma calculation starting with the file _SUS
!if (my_rank==master) then ! presently only master has chi0 in screening

 ! Temporary arrays to store spectra.
 ABI_ALLOCATE(epsm_lf,(nomega,my_nqlwl))
 ABI_ALLOCATE(epsm_nlf,(nomega,my_nqlwl))
 ABI_ALLOCATE(eelf,(nomega,my_nqlwl))
 epsm_lf =czero; epsm_nlf=czero; eelf= zero

 if (my_nqlwl>1)  then
   ABI_ALLOCATE(chi0_save,(npwe*nI,npwe*nJ))
 end if

 SELECT CASE (approx_type)

 CASE (0) ! RPA: \tepsilon=1 - Vc^{1/2} chi0 Vc^{1/2}
   !
   ! * vc_sqrt contains vc^{1/2}(q,G), complex-valued to allow for a possible cutoff.
   do io=1,nomega
     if (my_nqlwl>1)  chi0_save(:,:) = chi0(:,:,io)
     !
     ! * Loop over small q"s (if any) to treat the nonanalytical behavior.
     do iqlwl=my_nqlwl,1,-1

       if (my_nqlwl>1) then 
         chi0(:,:,io) = chi0_save              ! restore pristine polarizability
         chi0(:,1,io) = chi0_lwing(:,io,iqlwl) ! change the wings
         chi0(1,:,io) = chi0_uwing(:,io,iqlwl)
       end if

       if (iqibz==1) then 
         vc_sqrt => Vcp%vcqlwl_sqrt(:,iqlwl)  ! Use Coulomb term for q-->0
       else 
         vc_sqrt => Vcp%vc_sqrt(:,iqibz)  
       end if

       do ig2=1,npwe
!        (XG 090513) The following line has a problem with g95 compiler installed on max
!        chi0(:,ig2,io)=-vc_sqrt(:)*chi0(:,ig2,io)*vc_sqrt(ig2)
         do ig1=1,npwe*nI
           chi0(ig1,ig2,io)=-vc_sqrt(ig1)*chi0(ig1,ig2,io)*vc_sqrt(ig2)
         end do
         chi0(ig2,ig2,io)=one+chi0(ig2,ig2,io) 
       end do

       epsm_nlf(io,iqlwl)=chi0(1,1,io) ! * chi0(io), now contains \tepsilon(io).
       write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
       call wrtout(std_out,msg,'COLL')
       call print_arr(chi0(:,:,io))
       !
       ! === Invert tepsilon and calculate macroscopic dielectric constant ===
       ! * epsm_lf(w)=1/epsm1(G=0,Gp=0,w). 
       ! * Since G=Gp=0 there is no difference btw symmetrical and not symmetrical.
       !
       call xginv(chi0(:,:,io),npwe,comm=comm_self)

       epsm_lf(io,iqlwl) = one/chi0(1,1,io)
       eelf   (io,iqlwl) = -AIMAG(chi0(1,1,io))

       write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
       call wrtout(std_out,msg,'COLL')
       call print_arr(chi0(:,:,io))

       ! Save wings of e^-1 overwriting input values.
       if (dim_wing>0.and..FALSE.) then
         chi0_lwing(:,io,iqlwl) = chi0(:,1,io) 
         chi0_uwing(:,io,iqlwl) = chi0(1,:,io) 
       end if
      
     end do !iqlwl 
   end do ! nomega

 CASE (1) ! Vertex correction from Adiabatic TDDFT. chi_{G1,G2} = [\delta -\chi0 (vc+kxc)]^{-1}_{G1,G3} \chi0_{G3,G2}

   ABI_CHECK(Vcp%nqlwl==1,"nqlwl/=1 not coded")

   ! Make kxcg_mat(G1,G2) = kxcg(G1-G2) from kxcg defined on the FFT mesh.
   ABI_CHECK(nkxc==1,"nkxc/=1 not coded")
   ABI_ALLOCATE(kxcg_mat,(npwe,npwe))
   istat = ABI_ALLOC_STAT
   
   ierr=0
   do ig2=1,npwe
     do ig1=1,npwe
       g1mg2_idx = g2ifft(gvec(:,ig1)-gvec(:,ig2),ngfft)
       if (g1mg2_idx>0) then
         kxcg_mat(ig1,ig2) = kxcg(g1mg2_idx,1)
       else 
         ierr=ierr+1
         kxcg_mat(ig1,ig2) = czero 
       end if
     end do
   end do

   if (ierr/=0) then 
     write(msg,'(a,i4,3a)')&
&     ' Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
&     ' Enlarge the FFT mesh to get rid of this problem. '
     MSG_WARNING(msg)
   end if

   !FIXME "recheck TDDFT code and parallel"

   ABI_ALLOCATE(chitmp,(npwe,npwe))
   istat = ABI_ALLOC_STAT
   if (istat/=0) then 
     write(msg,'(a,f8.2,a)')"out-of-memory, requiring ",npwe**2*gwpc*b2Mb," Mb"
     MSG_ERROR(msg)
   end if

   if (iqibz==1) then 
     !$vc_sqrt => Vcp%vcqlwl_sqrt(:,iqlwl)  ! Use Coulomb term for q-->0
     vc_sqrt => Vcp%vcqlwl_sqrt(:,1)  ! TODO add treatment of non-Analytic behavior
   else 
     vc_sqrt => Vcp%vc_sqrt(:,iqibz)  
   end if

   do io=1,nomega
     ! * Calculate chi0*fxc.
     chitmp(:,:)=MATMUL(chi0(:,:,io),kxcg_mat(:,:)) 
     ! * First calculate the NLF contribution
     do ig1=1,npwe
       do ig2=1,npwe
         chitmp(ig1,ig2)=-chitmp(ig1,ig2)
       end do
       chitmp(ig1,ig1)=chitmp(ig1,ig1)+one
     end do
     call xginv(chitmp,npwe,comm=comm_self)
     chitmp(:,:)=MATMUL(chitmp(:,:),chi0(:,:,io))
     !if (.not. ABS(REAL(omega(io)))> tol3) call hermitianize(chitmp,"All")
     chitmp(1,1)=-vc_sqrt(1)*chitmp(1,1)*vc_sqrt(1)
     chitmp(1,1)=chitmp(1,1)+one

     epsm_nlf(io,1)=chitmp(1,1)     

     chitmp(:,:)=MATMUL(chi0(:,:,io),kxcg_mat(:,:)) 
     ! * Calculate (1-chi0*Vc-chi0*Kxc) and put it in chitmp.
     do ig1=1,npwe
       do ig2=1,npwe
         chitmp(ig1,ig2)=-chitmp(ig1,ig2)-chi0(ig1,ig2,io)*vc_sqrt(ig2)**2
       end do
       chitmp(ig1,ig1)=chitmp(ig1,ig1)+one
     end do

     ! * Invert (1-chi0*Vc-chi0*Kxc) and Multiply by chi0.
     call xginv(chitmp,npwe,comm=comm_self)
     chitmp=MATMUL(chitmp,chi0(:,:,io))

     ! * Save result, now chi0 contains chi.
     chi0(:,:,io)=chitmp

     write(std_out,'(a,i2,a,i1,a)')' chi(q= ',iqibz,',omega= ',io,',G,G")'
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

   SELECT CASE (option_test)

   CASE (0) ! Symmetrized TESTPARTICLE epsilon^-1
     call wrtout(std_out,' Calculating TESTPARTICLE epsilon^-1(G,G") = 1 + Vc*chi','COLL')
     do io=1,nomega
       do ig1=1,npwe
         chi0(ig1,:,io)=(vc_sqrt(ig1)*vc_sqrt(:))*chi0(ig1,:,io)
         chi0(ig1,ig1,io)=one+chi0(ig1,ig1,io)
       end do 
     end do 

   CASE (1) ! Symmetrized TESTELECTRON epsilon^-1
     call wrtout(std_out,' Calculating TESTELECTRON epsilon^-1(G,G") = 1 + (Vc + fxc)*chi',"COLL")
     do io=1,nomega
       chitmp=MATMUL(kxcg_mat(:,:),chi0(:,:,io))
       ! Perform hermitianization, only valid along the imaginary axis.
       if (.not. ABS(REAL(omega(io)))> tol3) call hermitianize(chitmp,"All")
       do ig1=1,npwe
         chi0(ig1,:,io)=(vc_sqrt(ig1)*vc_sqrt(:))*chi0(ig1,:,io)+chitmp(ig1,:)
         chi0(ig1,ig1,io)=one+chi0(ig1,ig1,io)
       end do 
     end do

   CASE DEFAULT 
     write(msg,'(a,i3)')'Wrong value for option_test= ',option_test
     MSG_BUG(msg)
   END SELECT 

   ABI_DEALLOCATE(chitmp)
   ABI_DEALLOCATE(kxcg_mat)
   !
   ! === chi0 now contains symmetrical epsm1 ===
   ! * Calculate macroscopic dielectric constant epsm_lf(w)=1/epsm1(G=0,Gp=0,w) ===
   epsm_lf(:,1) =  one/chi0(1,1,:)
   eelf   (:,1) = -AIMAG(chi0(1,1,:))
   do io=1,nomega
     write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

 CASE(2) ! ADA nonlocal vertex correction contained in fxc_ADA
   MSG_WARNING('Entered fxc_ADA branch: EXPERIMENTAL!')
   ! Test that argument was passed
   if (.NOT.present(fxc_ADA)) then
     MSG_ERROR('make_epsm1_driver was not called with optional argument fxc_ADA')
   end if
   
   ABI_CHECK(Vcp%nqlwl==1,"nqlwl/=1 not coded")

   ABI_ALLOCATE(chitmp,(npwe,npwe))
   istat = ABI_ALLOC_STAT
   if (istat/=0) then 
     write(msg,'(a,f8.2,a)')"out-of-memory, requiring ",npwe**2*gwpc*b2Mb," Mb"
     MSG_ERROR(msg)
   end if

   if (iqibz==1) then 
     !$vc_sqrt => Vcp%vcqlwl_sqrt(:,iqlwl)  ! Use Coulomb term for q-->0
     vc_sqrt => Vcp%vcqlwl_sqrt(:,1)  ! TODO add treatment of non-Analytic behavior
   else 
     vc_sqrt => Vcp%vc_sqrt(:,iqibz)  
   end if

   do io=1,nomega
     ! * Calculate chi0*fxc.
     chitmp(:,:)=MATMUL(chi0(:,:,io),fxc_ADA(:,:))
     ! * First calculate the NLF contribution
     do ig1=1,npwe
       do ig2=1,npwe
         chitmp(ig1,ig2)=-chitmp(ig1,ig2)
       end do
       chitmp(ig1,ig1)=chitmp(ig1,ig1)+one
     end do
     call xginv(chitmp,npwe,comm=comm_self)
     chitmp(:,:)=MATMUL(chitmp(:,:),chi0(:,:,io))
     !if (.not. ABS(REAL(omega(io)))> tol3) call hermitianize(chitmp,"All")
     chitmp(1,1)=-vc_sqrt(1)*chitmp(1,1)*vc_sqrt(1)
     chitmp(1,1)=chitmp(1,1)+one

     epsm_nlf(io,1)=chitmp(1,1)     

     ! Now do normal calculation
     chitmp(:,:)=MATMUL(chi0(:,:,io),fxc_ADA(:,:))
     ! * Calculate (1-chi0*Vc-chi0*fxc_ADA) and put it in chitmp.
     do ig1=1,npwe
       do ig2=1,npwe
         chitmp(ig1,ig2)=-chitmp(ig1,ig2)-chi0(ig1,ig2,io)*vc_sqrt(ig2)**2
       end do
       chitmp(ig1,ig1)=chitmp(ig1,ig1)+one
     end do

     ! * Invert (1-chi0*Vc-chi0*fxc_ADA) and Multiply by chi0.
     call xginv(chitmp,npwe,comm=comm_self)
     chitmp=MATMUL(chitmp,chi0(:,:,io))

     ! * Save result, now chi0 contains chi.
     chi0(:,:,io)=chitmp

     write(std_out,'(a,i2,a,i1,a)')' chi(q= ',iqibz,',omega= ',io,',G,G")'
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

   SELECT CASE (option_test)

   CASE (0) ! Symmetrized TESTPARTICLE epsilon^-1
     call wrtout(std_out,' Calculating TESTPARTICLE epsilon^-1(G,G") = 1 + Vc*chi','COLL')
     do io=1,nomega
       do ig1=1,npwe
         chi0(ig1,:,io)=(vc_sqrt(ig1)*vc_sqrt(:))*chi0(ig1,:,io)
         chi0(ig1,ig1,io)=one+chi0(ig1,ig1,io)
       end do 
     end do 

   CASE (1) ! Symmetrized TESTELECTRON epsilon^-1
     call wrtout(std_out,' Calculating TESTELECTRON epsilon^-1(G,G") = 1 + (Vc + fxc)*chi',"COLL")
     do io=1,nomega
       chitmp=MATMUL(fxc_ADA(:,:),chi0(:,:,io))
       ! Perform hermitianization, only valid along the imaginary axis.
       if (.not. ABS(REAL(omega(io)))> tol3) call hermitianize(chitmp,"All")
       do ig1=1,npwe
         chi0(ig1,:,io)=(vc_sqrt(ig1)*vc_sqrt(:))*chi0(ig1,:,io)+chitmp(ig1,:)
         chi0(ig1,ig1,io)=one+chi0(ig1,ig1,io)
       end do 
     end do

   CASE DEFAULT 
     write(msg,'(a,i3)')'Wrong value for option_test= ',option_test
     MSG_BUG(msg)
   END SELECT 

   ABI_DEALLOCATE(chitmp)
   !
   ! === chi0 now contains symmetrical epsm1 ===
   ! * Calculate macroscopic dielectric constant epsm_lf(w)=1/epsm1(G=0,Gp=0,w) ===
   epsm_lf(:,1)=one/chi0(1,1,:)
   eelf   (:,1) = -AIMAG(chi0(1,1,:))
   do io=1,nomega
     write(msg,'(a,i4,a,2f9.4,a)')' Symmetrical epsilon^-1(G,G'') at the ',io,' th omega',omega(io)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call print_arr(chi0(:,:,io),mode_paral='PERS')
   end do

 CASE DEFAULT
   write(msg,'(a,i3)')'Wrong value for approx_type= ',approx_type 
   MSG_BUG(msg)
 END SELECT

!FB: See comment above
!end if !master

 if (use_MPI==1) then ! * Collect results on each node.
   do io=1,nomega
     if (omega_distrb(io)/=my_rank) then 
       chi0(:,:,io)  = czero_gw
       chi0_lwing(:,io,:) = zero
       chi0_uwing(:,io,:) = zero
       !epsm_lf(io,:) = czero 
       !epsm_nlf(io,:) = czero 
       !eelf(io,:) = zero 
     end if
     call xsum_mpi(chi0(:,:,io), comm,ierr)
     call xsum_mpi(chi0_lwing(:,io,:),comm,ierr)
     call xsum_mpi(chi0_uwing(:,io,:),comm,ierr)
   end do
   call xsum_mpi(epsm_lf, comm,ierr )
   call xsum_mpi(epsm_nlf,comm,ierr)
   call xsum_mpi(eelf,    comm,ierr)
   call xbarrier_mpi(comm)
 end if

 ! * Save results in Spectra%, mind the slicing.
 Spectra%emacro_nlf(:,:) = epsm_nlf(1:nor,:)
 Spectra%emacro_lf (:,:) = epsm_lf (1:nor,:)
 Spectra%eelf      (:,:) = eelf    (1:nor,:)

 ABI_DEALLOCATE(epsm_lf)
 ABI_DEALLOCATE(epsm_nlf)
 ABI_DEALLOCATE(eelf)
 if (allocated(omega_distrb))  then
   ABI_DEALLOCATE(omega_distrb)
 end if
 if (allocated(chi0_save))  then
   ABI_DEALLOCATE(chi0_save)
 end if

 DBG_EXIT("COLL")

end subroutine make_epsm1_driver
!!***

!----------------------------------------------------------------------

!!****f* m_screening/make_W
!! NAME
!! make_W
!!
!! FUNCTION
!!  calculate the symmetrical inverse dielectric matrix starting
!!  from the irreducible polarizability. The routine considers a single q-point, and 
!!  performs the following tasks:
!!
!! INPUTS
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction
!!   %nqibz=Number of q-points.
!!   %qibz(3,nqibz)=q-points in the IBZ.  
!!  Er<Epsilonm1_results>=Data structure containing the inverse dielectric matrix. See also SIDE EFFECTS.
!!
!! SIDE EFFECTS
!!  Er%epsm1(npwe*nI,npwe*nJ,nomega): in input the symmetryzed inverse dieletric matrix, 
!!   in output the screened interaction (including a possible cutoff in real space)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_W(Er,Vcp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_W'
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(vcoul_t),intent(in) :: Vcp
 type(Epsilonm1_results),intent(inout) :: Er

!Local variables-------------------------------
!scalars
 integer :: iq_ibz,ig1,ig2,io
 real(dp) :: ucvol
 logical :: is_qeq0
 character(len=500) :: msg
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 complex(gwpc),pointer :: vc_sqrt(:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (Er%nI/=1.or.Er%nJ/=1) then 
   MSG_ERROR("nI or nJ=/1 not yet implemented")
 end if

 if (Er%ID/=4) then 
   write(msg,'(a,i0,a)')" found Er%ID = ",Er%ID," while it should be 4"
   MSG_ERROR(msg)
 end if

 !TODO mqmem==0, change info ER%, and update Hscr
 call metric(gmet,gprimd,-1,rmet,Vcp%rprimd,ucvol)

 do iq_ibz=1,Er%nqibz
   is_qeq0 = (normv(Er%qibz(:,iq_ibz),gmet,'G')<GW_TOLQ0) ! Check if q==0

   if (is_qeq0) then 
     vc_sqrt => Vcp%vcqlwl_sqrt(:,1)  ! Use Coulomb term for q-->0, only first Q is used, shall we average if nqlwl>1?
   else 
     vc_sqrt => Vcp%vc_sqrt(:,iq_ibz)  
   end if

   do io=1,Er%nomega
     do ig1=1,Er%npwe
       do ig2=1,Er%npwe
         Er%epsm1(ig1,ig2,io,iq_ibz) = Er%epsm1(ig1,ig2,io,iq_ibz) * vc_sqrt(ig1) * vc_sqrt(ig2)
       end do
     end do
   end do
   !
 end do ! nqibz

end subroutine make_W
!!***

!----------------------------------------------------------------------

!!****f* m_screening/mkem1_q0
!! NAME
!! mkem1_q0
!!
!! FUNCTION
!!   This routine construct the microscopic dieletric matrix for q-->0 starting from the heads, wings and the body
!!   of the irreducible polarizability. Afterwards it calculates the symmetrized inverse dieletric matrix 
!!   via a block wise inversion thus obtaining the heads and the wings of e^{-1} that can be 
!!   used to describe the non-analytic behavior for q-->0.
!!
!! INPUTS
!! npwe=Number of Gs used to describe chi0
!! nomega=Number of frequencies in chi0.
!! n1,n2=Factors used to define the same of the chi0 matrix (1,1 if collinear, the typical case)
!! Cryst<crystal_structure>=Structure describing the crystal structure.
!! Vcp<vcoul_t>=datatypes gathering info on the Coulomb term
!! gvec(3,npwe)=G-vector for chi0 in reduced coordinates.
!!
!! OUTPUT
!! eps_head(3,3,nomega)=The macroscopic dieletric tensor in reduced coordinates. 
!!   The dieletric matrix along versor \hat q can be obtained with 
!!     e(\hat q) = \hat q.eps_head \hat q if all quantities are given in Cartesian coordinates.
!!
!! SIDE EFFECTS
!! chi0(npwe*n1,npwe*n2,nomega)= Input: polarizability. output: inverse dieletric matrix (only the body is used)
!! chi0_lwing(npwe*n1,nomega,3) 
!! chi0_uwing(npwe*n2,nomega,3)  Input:  the lower and upper wings of the polarizability
!!                               Output: the "lower" and "upper" wings of the inverse dieletric matrix. See notes below.
!! chi0_head(3,3,nomega)= Input: the polarizability tensor in Cartesian coordinates.
!!                        Ouput: The "head" of the inverse dieletric matrix. See notes below.
!!
!! NOTES
!!  Matrix inversion in block form. 
!!
!!         1  n-1
!!  M =  | c  u^t| 1     ==>   M^{-1} =  |  1/k          -u^t A^{-1}/k                    | 
!!       | v  A  | n-1                   | -A^{-1} v/k    A^{-1} + (A^{-1}v u^t A^{-1})/k |
!!                                       
!!                             where k = c - u^t A^{-1} v
!!
!!  Let q be a versor in reciprocal space, the symmetrized dielectric matrix with bare coulomb interaction
!!  can be written as
!! 
!!  \tilde\epsilon = | q.Eq      q.Wl(G2) |  where   E_ij = \delta_ij -4\pi chi0_head_ij
!!                   | q.Wl(G1)  B(G1,G2  |          Wl(G1) = -4\pi chi0_lwing(G1)
!!                                                   Wu(G2) = -4\pi chi0_uwing(G1)
!!  therefore, in Cartesian coordinates, we have:
!!
!!  1) e^{-1}_{00}(q) = [ q_i q_j (E_{ij} - \sum_{GG'} Wu_i(G)a B_{GG'}^{-1} Wl_j(G')) ]^{-1} = 1/(q.Lq)
!!
!!  2) e^{-1}_{0G'}(q) = -e^{-1}_{00}(q) [ \sum_{iG} q_i Wu_i(G)a B_{GG'}^{-1} ] = (q.Su) /(q.Lq)
!! 
!!  3) e^{-1}_{G0}(q)  = -e^{-1}_{00}(q) [ \sum_{iG'} q_i B_{GG'}^{-1} Wl_i(G') ] = (q.Sl) /(q.Lq)
!!
!!  4) e^{-1}_{GG'}(q) = B_{GG'}^{-1} + 
!!     [ \sum_{ij} q_i q_j ( \sum_T B^{-1}_{GT}^{-1} Wl_i(T)) (\sum_T' Wu_j(T') B^{-1}_{T'G'}^{-1} ] / (q.Lq)
!!
!!  where Su(G,3) and Sl(G,3) are the "upper" and "lower" wings of the inverse dielectric matrix and 
!!  L is the inverse dielectric tensor. Similar equations hold even if vectors and tensors are given in terms
!!  of the reciprocal lattice vectors provided that the metric tensor is taken into account.
!!  The main difference is in the expression for the tensor as only one metric tensor can be 
!!  absorbed in the scalar product, the second metric multiplies one of the wings.
!!
!!  *) The present implementation assumes that no cutoff technique is used in the Coulomb term.
!!
!!  *) Once the tensor in know it is possible to average the quadratic form on the sphere exactly. 
!!  In Cartesian coordinates one obtains.
!!
!!    \dfrac{1}{4\pi} \int v.Tv d\Omega = Trace(T)/3
!!
!!  For the inverse dielectric matrix we have to resort to a numerical integration
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkem1_q0(npwe,n1,n2,nomega,Cryst,Vcp,gvec,chi0_head,chi0_lwing,chi0_uwing,chi0,eps_head)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkem1_q0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwe,nomega,n1,n2
 type(crystal_structure),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
!arrays
 integer,intent(in) :: gvec(3,npwe)
 complex(gwpc),intent(inout) :: chi0(npwe*n1,npwe*n2,nomega)
 complex(dpc),intent(inout) :: chi0_lwing(npwe*n1,nomega,3)
 complex(dpc),intent(inout) :: chi0_uwing(npwe*n2,nomega,3)
 complex(dpc),intent(inout) :: chi0_head(3,3,nomega)
 complex(dpc),intent(out) :: eps_head(3,3,nomega)

!Local variables ------------------------------
!scalars
 integer :: iomega,ig,ig1,ig2,idir,jdir
!arrays 
 real(dp),allocatable :: modg_inv(:)
 complex(dpc),allocatable :: eps_lwing(:,:),eps_uwing(:,:),eps_body(:,:),cvec(:)

!************************************************************************

 ABI_CHECK(npwe/=1,"npwe must be >1")

 ! Precompute 1/|G|.
 ABI_ALLOCATE(modg_inv,(npwe-1))
 do ig=1,npwe-1
   modg_inv(ig) = one/normv(gvec(:,ig+1),Cryst%gmet,'G')
 end do

 ABI_ALLOCATE(eps_uwing,((npwe-1)*n1,3))
 ABI_ALLOCATE(eps_lwing,((npwe-1)*n2,3))
 ABI_ALLOCATE(eps_body,(npwe-1,npwe-1))
 ABI_ALLOCATE(cvec,(npwe-1))

 do iomega=1,nomega
   !
   ! Head and wings of the symmetrized epsilon.
   eps_head(:,:,iomega) = -four_pi*chi0_head(:,:,iomega)
   do idir=1,3
     eps_head(idir,idir,iomega) = one + eps_head(idir,idir,iomega)
     eps_lwing(:,idir) = -four_pi * modg_inv * chi0_lwing(2:,iomega,idir) 
     eps_uwing(:,idir) = -four_pi * modg_inv * chi0_uwing(2:,iomega,idir) 
     !eps_lwing(:,idir) = -chi0_lwing(2:,iomega,idir) * SQRT(four_pi) * Vcp%vcqlwl_sqrt(2:npwe,1)
     !eps_uwing(:,idir) = -chi0_uwing(2:,iomega,idir) * SQRT(four_pi) * Vcp%vcqlwl_sqrt(2:npwe,1)
   end do

   write(std_out,*)" espilon head"
   call print_arr(eps_head(:,:,iomega))
   !
   ! Construct the body of the symmetrized epsilon then invert it.
   do ig2=1,npwe-1
     do ig1=1,npwe-1
       eps_body(ig1,ig2) = -four_pi * modg_inv(ig1)*chi0(ig1+1,ig2+1,iomega )*modg_inv(ig2)
       !eps_body(ig1,ig2) = -Vcp%vcqlwl_sqrt(ig1+1,1)*chi0(ig1+1,ig2+1,iomega)* Vcp%vcqlwl_sqrt(ig2+1,1)
     end do
     eps_body(ig2,ig2) = one + eps_body(ig2,ig2)
   end do

   call xginv(eps_body,npwe-1)
   !
   ! Overwrite chi0_head and chi0_wings with the head and the wings of the inverse dielectric matrix.
   do jdir=1,3
     !
     ! Head.
     cvec=czero
     do idir=1,3
       cvec = cvec + two_pi*Cryst%gmet(jdir,idir)*MATMUL(eps_body,eps_lwing(:,idir)) ! as we work in reciprocal coords.
     end do
     !cvec = MATMUL(eps_body,eps_lwing(:,jdir))
     do idir=1,3
       chi0_head(idir,jdir,iomega) = eps_head(idir,jdir,iomega) - xdotu(npwe-1,eps_uwing(:,idir),1,cvec,1)
     end do
     !
     ! Now the wings.
     chi0_uwing(2:,iomega,jdir) = -MATMUL(eps_uwing(:,jdir),eps_body) 
     chi0_lwing(2:,iomega,jdir) = -MATMUL(eps_body,eps_lwing(:,jdir)) 
     !
   end do !jdir

   write(std_out,*)"espilon^1 head after block inversion"
   call print_arr(chi0_head(:,:,iomega))
   !
   ! Change the body but do not add the corrections due to the head and the wings.
   ! since they can be obtained on the fly from eps_body and the wings of eps^{-1}.
   !$chi0(2:,2:,iomega) = eps_body
 end do !iomega

 ABI_DEALLOCATE(modg_inv)
 ABI_DEALLOCATE(cvec)
 ABI_DEALLOCATE(eps_lwing)
 ABI_DEALLOCATE(eps_uwing)
 ABI_DEALLOCATE(eps_body)

 RETURN
 ABI_UNUSED(Vcp%ng)

end subroutine mkem1_q0
!!***

!----------------------------------------------------------------------

!!****f* m_screening/lebedev_laikov_int
!! NAME
!!  lebedev_laikov_int
!!
!! FUNCTION
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


subroutine lebedev_laikov_int()  

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lebedev_laikov_int'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 
!arrays

!Local variables-------------------------------
!scalars
 integer :: on,npts,ii,ierr,ll,mm,lmax,leb_idx
 real(dp) :: accuracy
 complex(dpc) :: ang_int
!arrays
 real(dp) :: cart_vpt(3),real_pars(0)
 real(dp),allocatable :: vx(:),vy(:),vz(:),ww(:)
 complex(dpc) :: tensor(3,3),cplx_pars(9)
 complex(dpc),allocatable :: ref_func(:),expd_func(:) !tmp_momenta(:)

! *************************************************************************

 !tensor=RESHAPE((/4.0,2.0,4.0,0.5,2.1,0.0,5.4,2.1,5.0/),(/3,3/))
 tensor=RESHAPE((/4.0,0.0,0.0,0.0,4.0,0.0,0.0,0.0,5.0/),(/3,3/))
 !tensor=RESHAPE((/1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/),(/3,3/))

 npts=26
 ABI_ALLOCATE(vx,(npts))
 ABI_ALLOCATE(vy,(npts))
 ABI_ALLOCATE(vz,(npts))
 ABI_ALLOCATE(ww,(npts))

 call LD0026(vx,vy,vz,ww,on)

 ang_int=czero
 do ii=1,npts
   cart_vpt = (/vx(ii),vy(ii),vz(ii)/)
   ang_int = ang_int + ww(ii)*DOT_PRODUCT(cart_vpt,MATMUL(tensor,cart_vpt))
 end do

 !write(std_out,*)"quadratic form associated to tensor=",tensor 
 write(std_out,*)"on ang_int",on,ang_int

 ABI_DEALLOCATE(vx)
 ABI_DEALLOCATE(vy)
 ABI_DEALLOCATE(vz)
 ABI_DEALLOCATE(ww)

 call init_lebedev_gridset()
 cplx_pars = RESHAPE(tensor,(/9/)); accuracy=tol10

 ! This is the function to be expanded evaluated on the lebedev_laikov grid of index leb_idx 
 leb_idx=3; npts=lebedev_npts(leb_idx)
 ABI_ALLOCATE(ref_func,(npts))
 do ii=1,npts
   cart_vpt = Lgridset(leb_idx)%versor(:,ii)
   ref_func(ii) = one/DOT_PRODUCT(cart_vpt,MATMUL(tensor,cart_vpt))
 end do

 ! Calculate the expansion in angular momenta of 1/{q.Tq}. 
 ! Only even l-components contribute thanks to the parity of the integrand.
 ! tol6 seems to be an acceptable error, convergence wrt lmax is very slow even for simple tensors.
 ABI_ALLOCATE(expd_func,(npts))
 expd_func=czero
 lmax=10
 do ll=0,lmax,2
   !allocate(tmp_momenta(-ll:ll))
   do mm=-ll,ll
     call lebedev_quadrature(ylmstar_over_qTq,(/ll,mm/),real_pars,cplx_pars,ang_int,ierr,accuracy)
     write(std_out,*)ll,mm,ang_int
     !tmp_momenta(mm) = ang_int
     do ii=1,npts
       cart_vpt = Lgridset(leb_idx)%versor(:,ii)
       expd_func(ii) = expd_func(ii) + four_pi*ang_int*ylmc(ll,mm,cart_vpt)
     end do
   end do
   !deallocate(tmp_momenta)
   write(std_out,*)"Error in angular expansion at l=",ll," is ",MAXVAL(ABS(expd_func-ref_func))
 end do

!BEGINDEBUG
 do ii=1,npts
   write(777,*)ref_func(ii)
   write(778,*)expd_func(ii)
 end do
!ENDDEBUG

 ABI_DEALLOCATE(expd_func)
 ABI_DEALLOCATE(ref_func)
 call destroy_lebedev_gridset()

 MSG_ERROR("Exiting from lebedev_laikov_int")

end subroutine lebedev_laikov_int
!!***

!----------------------------------------------------------------------

!!****f* m_screening/ylmstar_over_qTq
!! NAME
!!  ylmstar_over_qTq
!!
!! FUNCTION
!!  Return Ylm(q)^*/(q.Tq) where q is a versor in Cartesian coordinates.
!!  and Ylm is a complex spherical Harmonics whose index (l,m) are 
!!  passed via int_pars(1:2). T is a tensore in Cartesian coordinates
!!  passed via cplx_pars(1:9).
!!
!! INPUTS
!!  cart_vers(3)=Cartesian components of the versor
!!  int_pars(1:2)=(l,m) indeces in Ylm. l>=0 and  m \in [-l,l]
!!  cplx_pars(1:9)=Tensor T in Cartesian coordinates.
!!  real_pars=Not used.
!!
!! OUTPUT
!!  Value of Ylm(q)^*/(q.Tq) 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ylmstar_over_qTq(cart_vers,int_pars,real_pars,cplx_pars)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ylmstar_over_qTq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: cart_vers(3)
 integer,intent(in) :: int_pars(:)
 real(dp),intent(in) :: real_pars(:)
 complex(dpc),intent(in) :: cplx_pars(:)
 complex(dpc) :: ylmstar_over_qTq 
!arrays

!Local variables-------------------------------
!scalars
 integer :: ll,mm
!arrays
 complex(dpc) :: tensor(3,3)

! *************************************************************************

 tensor = RESHAPE(cplx_pars(1:9),(/3,3/))
 ll = int_pars(1) ! ll starts from zero. 
 mm = int_pars(2) ! m \in [-l,l]

 ylmstar_over_qTq = CONJG(ylmc(ll,mm,cart_vers))/DOT_PRODUCT(cart_vers,MATMUL(tensor,cart_vers))

 RETURN
 ABI_UNUSED(real_pars(1))
 
end function ylmstar_over_qTq
!!***

!----------------------------------------------------------------------

!!****f* m_screening/ylmstar_wtq_over_qTq
!! NAME
!!  ylmstar_wtq_over_qTq
!!
!! FUNCTION
!!  Return Ylm(q)^* weight(q)/(q.Tq) where q is a versor in Cartesian coordinates.
!!  Ylm is a complex spherical Harmonics whose index (l,m) are 
!!  passed via int_pars(1:2). T is a tensor in Cartesian coordinates
!!  passed via cplx_pars(1:9). weight(q) is the weighting function giving
!!  the length of the vector parallel to versor q that connects the origin
!!  of the lattice to one of the boundaries of the small cell centered at Gamma
!!
!! INPUTS
!!  cart_vers(3)=Cartesian components of the versor
!!  int_pars(1:2)=(l,m) indeces in Ylm. l>=0 and  m \in [-l,l]
!!  cplx_pars(1:9)=Tensor T in Cartesian coordinates.
!!  real_pars(1:9)=The Cartesian vectors defining the small box centered around gamma point
!!    when referred to this vectors the points in the box are given by {(x,y,z) | x,y,z \in [-1,1]}.
!!
!! OUTPUT
!!  Value of Ylm(q)^* weigh(q)/(q.Tq) 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ylmstar_wtq_over_qTq(cart_vers,int_pars,real_pars,cplx_pars)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ylmstar_wtq_over_qTq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: cart_vers(3)
 integer,intent(in) :: int_pars(:)
 real(dp),intent(in) :: real_pars(:)
 complex(dpc),intent(in) :: cplx_pars(:)
 complex(dpc) :: ylmstar_wtq_over_qTq 
!arrays

!Local variables-------------------------------
!scalars
 integer :: ll,mm
 real(dp) :: wtq
!arrays
 real(dp) :: gprimd(3,3),rprimd(3,3),red_vers(3)
 complex(dpc) :: tensor(3,3)

! *************************************************************************

 MSG_ERROR("Work in progress")
 ! box_len has to be tested

 gprimd = RESHAPE(real_pars(1:9),(/3,3/))
 red_vers = MATMUL(rprimd,cart_vers)
 wtq = box_len(red_vers,gprimd)

 tensor = RESHAPE(cplx_pars(1:9),(/3,3/))
 ll = int_pars(1) ! true ll i.e. not shifted
 mm = int_pars(2)

 ylmstar_wtq_over_qTq = CONJG(ylmc(ll,mm,cart_vers))*wtq/DOT_PRODUCT(cart_vers,MATMUL(tensor,cart_vers))
 
end function ylmstar_wtq_over_qTq
!!***

!----------------------------------------------------------------------

!!****f* m_screening/k_fermi
!! NAME
!!  k_fermi
!!
!! FUNCTION
!!  Returns the Fermi wave vector corresponding to the local value of the real space density rhor.
!!
!! INPUTS
!!  rhor=Local density in real space.
!!
!! TODO
!!  Might be defined as elemental function.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function k_fermi(rhor)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'k_fermi'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rhor
 real(dp) :: k_fermi
!arrays

!Local variables-------------------------------
!scalars
 real(dp),parameter :: pisq=pi**2

! *************************************************************************

 k_fermi = (three*pisq*rhor)**third

end function k_fermi
!!***

!----------------------------------------------------------------------

!!****f* m_screening/k_thfermi
!! NAME
!!  k_thfermi
!!
!! FUNCTION
!!  Returns the Thomas-Fermi wave vector corresponding to the local value of the real space density rhor.
!!
!! INPUTS
!!  rhor=Local density in real space.
!!
!! TODO
!!  Might be defined as elemental function.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function k_thfermi(rhor)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'k_thfermi'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: rhor
 real(dp) :: k_thfermi

!Local variables-------------------------------
!scalars
 real(dp),parameter :: pisq=pi**2

! *************************************************************************

 k_thfermi = SQRT(four*k_fermi(rhor)*piinv)

end function k_thfermi
!!***

!----------------------------------------------------------------------

!!****f* m_screening/mdielf_bechstedt      
!! NAME
!!  mdielf_bechstedt
!!
!! FUNCTION
!!  Calculates the model dielectric function for the homogeneous system 
!!  as proposed by F. Bechstedt, in Solid State Commun. 84, 765 1992.
!!
!! INPUTS
!!  eps_inf=Dielectric constant of the material
!!  qnrm=The modulus of the q-point.
!!  rhor=The local value of the density 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function mdielf_bechstedt(eps_inf,qnrm,rhor) result(mdielf)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mdielf_bechstedt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: eps_inf,qnrm,rhor
 real(dp) :: mdielf

! *************************************************************************

 mdielf = one + & 
&  one / ( one/(eps_inf-one) + (qnrm/k_thfermi(rhor))**2 + (three*qnrm**4)/(four*k_fermi(rhor)**2 * k_thfermi(rhor)**2) )

end function mdielf_bechstedt
!!***

!----------------------------------------------------------------------

!!****f* m_screening/screen_mdielf
!! NAME
!!  screen_mdielf
!!
!! FUNCTION
!!  Calculates W_{G,G'}(q,w) for a given q-point in the BZ using a model dielectric function.
!!
!! INPUTS
!!  iq_bz=The index of the q-point in the BZ where W(q) is calculated.
!!  npw=Number of plane waves for W
!!  nomega=Number of frequency points.     
!!  model_type=Flag defining the model.
!!  eps_inf=Dielectric constant of the material. 
!!  Cryst<crystal_structure>=Info on the unit cell
!!  Qmesh<BZ_mesh_type>=Info on the set of q-points.
!!  Vcp<vcoul_t datatype>= containing information on the cutoff technique
!!  Gsph<Gsphere>=The G-sphere for W. 
!!  nspden=Number of spin density components of the density.
!!  nfft=Number of FFT points on the dense FFT mesh
!!  ngfft(18)=contain all needed information about 3D FFT.
!!  rhor(nfft,nspden)=Electron density in real space (The PAW AE term is included)
!!  which= Set it to "EM1" if the symmetrized inverse dielectric matrix is wanted. 
!!   By default the routines returns W.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  w_qbz(npw,npw,nomega)
!!
!! NOTES
!!   W_{G1,G2} =  1/2 { 
!!     v(q+G1) \int em1(|q+G1|,r) e^{-i(G1-G2).r} dr  +
!!     v(q+G2) \int em1(|q+G2|,r) e^{-i(G1-G2).r} dr } / \Omega
!!
!! PARENTS
!!      m_screen,m_shexc
!!
!! CHILDREN
!!
!! SOURCE

subroutine screen_mdielf(iq_bz,npw,nomega,model_type,eps_inf,Cryst,Qmesh,Vcp,Gsph,nspden,nfft,ngfft,rhor,which,w_qbz,comm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screen_mdielf'
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nomega,nfft,nspden,iq_bz,comm,model_type
 real(dp),intent(in) :: eps_inf
 character(len=*),intent(in) :: which
 type(BZ_mesh_type),intent(in) :: Qmesh
 type(crystal_structure),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
 type(gvectors_type),intent(in) :: Gsph
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,nspden)
 complex(gwpc),intent(out) :: w_qbz(npw,npw,nomega)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0,paral_kgb0=0,cplex1=1
 integer :: my_gstart,my_gstop,iq_ibz,ig,itim_q,isym_q,istat
 integer :: ig1,ig2,g1mg2_fft,iw,ii,ierr,nprocs,isg,ifft !,row ,col
 real(dp) :: qpg2_nrm
 complex(dpc) :: ph_mqbzt !,ctmp1,ctmp2
 logical :: is_qeq0,isirred
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq 
!arrays
 integer :: umklp(3)
 integer,allocatable :: igfft(:),g1mg2(:,:)
 real(dp) :: qpg2(3),qpt_bz(3)
 real(dp),allocatable :: em1_qpg2r(:),fofg(:,:)
 complex(gwpc),pointer :: vc_sqrt_ibz(:)
 complex(gwpc),allocatable :: vc_qbz(:),ctmp(:,:)
 logical,allocatable :: mask(:)

! *************************************************************************

 ABI_CHECK(nomega==1,"screen_mdielf does not support nomega>1")

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)

 nprocs = xcomm_size(comm)
 call xmpi_split_work(npw,comm,my_gstart,my_gstop,msg,ierr)

 call get_bz_item(Qmesh,iq_bz,qpt_bz,iq_ibz,isym_q,itim_q,ph_mqbzt,umklp,isirred)

 !if (itim_q/=1.or.isym_q/=1.or.ANY(umklp/=0) ) then
 !  MSG_ERROR("Bug in mdielf_bechstedt")
 !end if
 ! 
 ! Symmetrize Vc in the full BZ.
 is_qeq0 = (normv(qpt_bz,Cryst%gmet,'G')<GW_TOLQ0) ! Check if q==0
 if (is_qeq0) then 
   vc_sqrt_ibz => Vcp%vcqlwl_sqrt(:,1)  ! Use Coulomb term for q-->0, only first Q is used, shall we average if nqlwl>1?
 else 
   vc_sqrt_ibz => Vcp%vc_sqrt(:,iq_ibz)  
 end if

 ABI_ALLOCATE(vc_qbz,(npw))
 do ig=1,npw
   isg = Gsph%rottb(ig,itim_q,isym_q) 
   vc_qbz(isg) = vc_sqrt_ibz(ig)**2
 end do

 ABI_ALLOCATE(igfft,(npw))
 ABI_ALLOCATE(g1mg2,(3,npw))
 ABI_ALLOCATE(fofg,(2,nfft))
 ABI_ALLOCATE(em1_qpg2r,(nfft))
 ABI_ALLOCATE(mask,(npw))

 w_qbz=czero
 do ig2=my_gstart,my_gstop
   !
   ! Compute the index of G-G2 wave in the FFT grid.
   do ii=1,npw
     g1mg2(:,ii) = Gsph%gvec(:,ii) - Gsph%gvec(:,ig2)
   end do
   call kgindex(igfft,g1mg2,mask,MPI_enreg_seq,ngfft,npw)
 
   ! TODO can use zero-padding FFT to speed up the transform.
   !call sphereboundary(gbound,istwfk1,g1mg2,mgfft,npw)

   ! Evaluate em1_qpg2r = \int em1(|q+G2|,r) e^{-i(G1-G2).r} dr }.
   qpg2 = qpt_bz + Gsph%gvec(:,ig2) 
   qpg2_nrm = normv(qpg2,Cryst%gmet,"G")

   do iw=1,nomega
     !
     select case (model_type)
     case (1)
       do ifft=1,nfft
         em1_qpg2r(ifft) = one / mdielf_bechstedt(eps_inf,qpg2_nrm,rhor(ifft,1)) 
       end do
     case default
       write(msg,'(a,i0)')" Unknown value for model_type ",model_type
       MSG_ERROR(msg)
     end select 

     call fourdp(cplex1,fofg,em1_qpg2r,-1,MPI_enreg_seq,nfft,ngfft,paral_kgb0,tim_fourdp0)
     !
     ! Here, unlike the other parts of the code, the unsymmetrized e^{-1} is used.
     do ig1=1,npw
       g1mg2_fft = igfft(ig1)
       w_qbz(ig1,ig2,iw) = DCMPLX(fofg(1,g1mg2_fft), fofg(2,g1mg2_fft)) * vc_qbz(ig2) !/ Cryst%ucvol
     end do
   end do ! iw
   !
 end do ! ig2

 ABI_DEALLOCATE(em1_qpg2r)
 ABI_DEALLOCATE(fofg)
 ABI_DEALLOCATE(igfft)
 ABI_DEALLOCATE(g1mg2)
 ABI_DEALLOCATE(mask)
 !
 ! W = 1/2 * (A + A^H) 
 ! The MPI sum is done inside the loop to avoid problems with the size of the packet.
 ABI_ALLOCATE(ctmp,(npw,npw))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in ctmp")

 do iw=1,nomega
   !% do col=1,npw
   !%   do row=1,col
   !%     ctmp1 = w_qbz(row,col,iw)
   !%     ctmp2 = CONJG(w_qbz(col,row,iw))
   !%     w_qbz(row,col,iw) = half * (ctmp1 + ctmp2) 
   !%     if (col/=row) w_qbz(col,row,iw) = CONJG(w_qbz(row,col,iw))
   !%   end do
   !% end do
   ctmp = TRANSPOSE(CONJG(w_qbz(:,:,iw)))
   w_qbz(:,:,iw) = half * (ctmp + w_qbz(:,:,iw))
   call xsum_mpi(w_qbz(:,:,iw),comm,ierr)
 end do
 !
 ! Calculate the symmetrized Em1. W = vc(G1)^{1/2} \tilde Em1 vc(G2)^{1/2} -------------------------
 if (toupper(which)=="EM1") then 
   do ig=1,npw
     isg = Gsph%rottb(ig,itim_q,isym_q) 
     vc_qbz(isg) = vc_sqrt_ibz(ig)  ! Workspace storing vc*{1/2}(q_BZ,G).
   end do

   do ig2=1,npw
     do ig1=1,npw
       ctmp(ig1,ig2) =  one / (vc_qbz(ig1) * vc_qbz(ig2))
     end do
   end do

   do iw=1,nomega
     w_qbz(:,:,iw) = w_qbz(:,:,iw) * ctmp(:,:)
   end do
 end if

 ABI_DEALLOCATE(vc_qbz)
 if (allocated(ctmp))  then
   ABI_DEALLOCATE(ctmp)
   istat = ABI_ALLOC_STAT
 end if

end subroutine screen_mdielf
!!***

!----------------------------------------------------------------------

!!****f* m_screening/interpolate_w
!! NAME
!! interpolate_w
!!
!! FUNCTION
!!
!! INPUTS
!!  in_fname
!!  out_fname   
!!  new_nkibz
!!  new_kibz(3,new_nkibz)
!!  old_nshiftq
!!  old_qptrlatt(3,3) 
!!  old_shiftq(3,old_nshiftq)
!!
!! OUTPUT
!!   The interpolated matrices are written on file out_fname.
!!
!! NOTES
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine interpolate_w(in_fname,out_fname,old_qptrlatt,old_nshiftq,old_shiftq,new_nkibz,new_kibz,comm)

 use defs_basis

 !use m_crystal
 !use m_crystal_io
 !use m_gsphere
 !use m_bz_mesh
 !use m_screening
 !use m_io_screening

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'interpolate_w'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: new_nkibz,old_nshiftq,comm
 character(len=*),intent(in) :: in_fname,out_fname   
!arrays
 integer,intent(in) :: old_qptrlatt(3,3)
 real(dp),intent(in) :: new_kibz(3,new_nkibz),old_shiftq(3,old_nshiftq)

!Local variables-------------------------------
!scalars
 integer,parameter :: nvert=8,iout2=2,option2=2 !brav1=1,
 integer :: rdwr,master,my_rank,npw,nomega,accesswff,vert,dir !ierr, iq_bz,
 integer :: ii,jj,i1,i2,i3,in_unt,out_unt,ios,old_nqibz,qptopt,fform,timrev,istat,nref_qpt
 integer :: vert_bz
 integer :: max_ond,intq,qref !,vt_bz,vt_ibz,vt_sym,vt_itim !,vt_ph
 integer :: isym_q,itim_q,iq_ibz,iomega,isg1,isg2
 integer :: g1mg0_idx,g2mg0_idx,g0_idx,ig1,ig2
 real(dp) :: tmp_rc,wrap_rc,shift,vfact
 character(len=500) :: msg
 logical :: found,isirred
 type(crystal_structure) :: Cryst
 type(ScrHdr_type) :: Old_Hscr,New_Hscr
 type(BZ_mesh_type) :: Old_Qmesh,New_Qmesh,New_Kmesh
 type(gvectors_type) :: Gsph
!arrays
 integer,pointer :: grottb(:)
 integer :: ondiv(3),bl_idx(3),vert_g0(3),g0(3) !base_g0(3)
 integer,allocatable :: gvec(:,:),iperm(:)
 real(dp) :: newq(3),vert_box(3,nvert),qstep(3),bl_red(3),dummy_qbz(3)
 real(dp),allocatable :: old_qibz(:,:),ref_qbz(:,:) !,spkpt(:,:)
 !integer,allocatable :: bsp2bz(:) !,iperm(:)
 !real(dp),allocatable :: spkpt(:,:),seen(:) 
 real(dp),allocatable :: qvec(:,:) !,yvec(:),zvec(:) !,xknot(:),yknot(:),zknot(:)
 !real(dp),allocatable :: xyzdata(:,:,:),bs_coef(:,:,:),work_bz(:,:,:,:),rout(:,:)
 !complex(dpc),allocatable :: cwork_bz(:,:,:),umat_k(:,:)
 complex(gwpc),pointer :: phmSgt(:)
 complex(gwpc),allocatable :: old_mat(:,:,:,:),cwork(:,:,:),new_mat(:,:,:,:),swap(:,:)

! *************************************************************************

 master    = 0
 my_rank   = xcomm_rank(comm)
 !comm_self = xmpi_self
 accesswff = IO_MODE_FORTRAN

 ! =============================================
 ! ==== Open input file and read the header ====
 ! =============================================
 !
 !if (my_rank==master) then
 call wrtout(std_out," Testing file: "//TRIM(in_fname),"COLL")
 in_unt=get_unit()
 open(unit=in_unt,file=in_fname,status='old',form='unformatted',iostat=ios)
 ABI_CHECK(ios==0,'Opening file '//TRIM(in_fname))

 rdwr=5
 call scr_hdr_io(fform,rdwr,in_unt,comm,master,accesswff,Old_Hscr)
 close(in_unt)
 !
 ! Echo of the header.
 if (my_rank==master) then
   rdwr=4
   call scr_hdr_io(fform,rdwr,std_out,comm,master,accesswff,Old_Hscr)
 end if
 !
 ! Basic dimensions.
 npw       = Old_Hscr%npwe
 nomega    = Old_Hscr%nomega
 old_nqibz = Old_Hscr%nqibz

 ABI_ALLOCATE(old_qibz,(3,old_nqibz))
 old_qibz = Old_Hscr%qibz
 !
 ABI_ALLOCATE(gvec,(3,npw))
 gvec=Old_Hscr%gvec(:,1:npw)
 !
 ! ============================
 ! ==== Init basic objects ====
 ! ============================
 qptopt=3
 timrev=2 ! This information is not reported in the header
          ! 1 => do not use time-reversal symmetry 
          ! 2 => take advantage of time-reversal symmetry
 call init_crystal_from_hdr(Cryst,Old_Hscr%Hdr,timrev)
 !call print_crystal(Cryst)

 call init_gsphere(Gsph,.FALSE.,Cryst,npw,gvec=gvec)
 ABI_DEALLOCATE(gvec)

 ! FIXME: Here I am assuming that kptrlatt is diagonal, no shift is expected!
 ondiv(1) = old_qptrlatt(1,1)
 ondiv(2) = old_qptrlatt(2,2)
 ondiv(3) = old_qptrlatt(3,3)
 max_ond = MAXVAL(ondiv)
 nref_qpt = PRODUCT(ondiv)

 ABI_ALLOCATE(qvec,(max_ond,3))
 ABI_ALLOCATE(iperm,(max_ond))

 ABI_CHECK(old_nshiftq==1,"Multiple shifts not supported")

 do dir=1,3
   qstep(dir) = one/ondiv(dir)
   do ii=1,ondiv(dir)
     tmp_rc = ( (ii-1) + old_shiftq(dir,1) ) * qstep(dir)
     call wrap2_pmhalf(tmp_rc,wrap_rc,shift) ! Wrap the trial values in the interval ]-1/2,1/2] .
     qvec(ii,dir) = wrap_rc
   end do
   !
   ! Sort the reduced coordinates.
   call sort_dp(ondiv(dir),qvec(:,dir),iperm,tol12)
 end do

 ABI_DEALLOCATE(iperm)
 !
 ! Contruct the reference mesh.
 ABI_ALLOCATE(ref_qbz,(3,nref_qpt))
 qref=0
 do i3=1,ondiv(3)
   do i2=1,ondiv(2)
     do i1=1,ondiv(1)
       qref=qref+1
       ref_qbz(:,qref) = (/qvec(i1,1), qvec(i2,2), qvec(i3,3)/)
       write(std_out,*)"ref_qbz ",ref_qbz(:,qref)
     end do
   end do
 end do
 
 ! FIXME here there's a bug when ref_qbz is passed.
 call init_kmesh(Old_Qmesh,Cryst,old_nqibz,old_qibz,qptopt,ref_bz=ref_qbz)
 !call init_kmesh(Old_Qmesh,Cryst,old_nqibz,old_qibz,qptopt)
 call print_bz_mesh(Old_Qmesh,prtvol=100)
 
 ABI_DEALLOCATE(old_qibz)
 ABI_DEALLOCATE(ref_qbz)
 
 ! Find the new mesh for W from the input k-point sampling.
 call init_kmesh(New_Kmesh,Cryst,new_nkibz,new_kibz,qptopt)
 call find_qmesh(New_Qmesh,Cryst,new_Kmesh)
 call destroy_bz_mesh_type(New_Kmesh)
 !  
 !  === Write chi0 on _SUSC file ===
 !  * Master creates and write the header if this is the first q-point calculated.
 if (my_rank==master) then ! * Open file and write header for polarizability  files.
   out_unt = get_unit()
   open(unit=out_unt,file=out_fname,status='unknown',form='unformatted',iostat=ios)
   ABI_CHECK(ios==0,' Opening '//TRIM(out_fname))
   ! Update info on the q-points reported in the header.
   call copy_scrhdr(Old_Hscr,New_Hscr)       
   New_Hscr%nqibz = New_Qmesh%nibz
   ABI_DEALLOCATE(New_Hscr%qibz)
   call deep_copy(New_Qmesh%ibz,New_Hscr%qibz)
   !
   rdwr=2
   call scr_hdr_io(New_Hscr%fform,rdwr,out_unt,xmpi_self,master,accesswff,New_Hscr)
   call free_scrhdr(New_Hscr)
 end if
 !
 ! Read tables from files.
 ABI_ALLOCATE(old_mat,(npw,npw,nomega,old_nqibz))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in old_mat")
                                                                           
 call read_screening(in_fname,npw,old_nqibz,nomega,old_mat,accesswff,comm)
 !
 ! Linear interpolation in 3D.
 !
 ! workspace array used to store F_{G,G'}(w) for the 8 vertices that define
 ! the box that contais the q-point for the interpolation.
 ABI_ALLOCATE(cwork,(npw,npw,nomega))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in cwork")
 !
 ! The interpolated matrix
 ABI_ALLOCATE(new_mat,(npw,npw,nomega,1))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in new_mat")

 do intq=1,New_Qmesh%nibz
   ! FIXME make sure that ibz is in [-half,+half] so that we can bracket it.
   !call wrap2_pmhalf(New_Qmesh%ibz(1,intq),newq(1),base_g0(1)) 
   !call wrap2_pmhalf(New_Qmesh%ibz(2,intq),newq(2),base_g0(2)) 
   !call wrap2_pmhalf(New_Qmesh%ibz(3,intq),newq(3),base_g0(3)) 
   newq = New_Qmesh%ibz(:,intq)
   !
   ! Find the vertex of the box enclosing newq.
   ! Note that we have to take into account the case in which the box 
   ! extends outside the first Brillouin zone BZ since F(q) is not periodic.
   ! For this reason the vertex are calculated explicitly.
   ! This is the difficult part since %bz is not ordered!
   !call interpol3d_indices(newq,ondiv(1),ondiv(2),ondiv(3),bl(1),bl(2),bl(3), tr1,tr2,tr3)
   !write(std_out,'(a,3f5.2,2(a,3i3))')"newq ",newq," is enclosed by bleft",bl(1),bl(2),bl(3)," top right",tr1,tr2,tr3

   bl_idx(1) = bisect(qvec(:,1),newq(1)) 
   bl_idx(2) = bisect(qvec(:,2),newq(2)) 
   bl_idx(3) = bisect(qvec(:,3),newq(3)) 
   !
   do dir=1,3
     if (bl_idx(dir)==0) then ! Bracket on the left side.
       bl_red(dir) = qvec(1,dir) - qstep(dir)
       if (bl_red(dir) > newq(dir)) then
         MSG_ERROR("Cannot bracket newq on the left side")
        end if
     else if (bl_idx(dir)==ondiv(dir)) then 
       bl_red(dir) = qvec(bl_idx(dir), dir)
       ! Make sure we can bracket on the right side.
       if ( bl_red(dir) + qstep(dir) <= newq(dir)) then
         MSG_ERROR("Cannot bracket newq on the right side")
       end if
     else 
       bl_red(dir) = qvec(bl_idx(dir), dir)
     end if
   end do
   !
   ! Now build the vertex of the box.
   vert=0
   do i3=1,2
     do i2=1,2
       do i1=1,2
         vert=vert+1
         vert_box(:,vert) = bl_red + (/(i1-1)*qstep(1), (i2-1)*qstep(2), (i3-1)*qstep(3)/)
         write(std_out,'(a,3f5.2,a,3(2x,f5.2))')"newq ",newq," has vert:",vert_box(:,vert)
       end do
     end do
   end do
   !
   ! For each vertex: symmetrize F and calculate its contribution to the interpolated result.
   ! Here a non-zero umklapp Might be needed!
   ! * Remember the symmetry properties of E_{GG'}(q)
   !  If q_bz=Sq_ibz+G0:
   ! 
   !  $ E_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau} E_{G1,G2)}(q)
   !
   !    The invariance under exchange of the real space position E(1,2) = E(2,1) leads to:
   !      $ E_{-G2,-G1}(-q) = E_{G1,G2)
   !
   new_mat = czero
   do vert=1,nvert
     !
     found = has_bz_item(Old_Qmesh,vert_box(:,vert),vert_bz,vert_g0)
     if (.not.found) then
       write(msg,'(a,3f5.2,a)')" vertex: ",vert_box(:,vert)," not in the BZ"
       MSG_ERROR(msg)
     end if
     write(std_out,*)"vert_g0",vert_g0
     !
     ! * Get iq_ibz, and symmetries from vert_bz.
     !   vert = IS q_ibz + (vert_g0 + g0)
     call get_BZ_item(Old_Qmesh,vert_bz,dummy_qbz,iq_ibz,isym_q,itim_q,umklp=g0,isirred=isirred)
     write(std_out,*)"umklp",g0
     !
     ! Symmetrize the matrix 
     grottb => Gsph%rottb (1:npw,itim_q,isym_q)
     phmSgt => Gsph%phmSGt(1:npw,isym_q) 

     do iomega=1,nomega
       do jj=1,npw
         isg2 = grottb(jj)
         do ii=1,npw
           isg1 = grottb(ii)
           cwork(isg1,isg2,iomega) = old_mat(ii,jj,iomega,iq_ibz) * phmSgt(ii) * CONJG(phmSgt(jj))
         end do
       end do
     end do
     !
     ! Account for time-reversal.
     if (itim_q==2) then
       do iomega=1,nomega
         cwork(:,:,iomega)=TRANSPOSE(cwork(:,:,iomega))
       end do
     end if
     !
     ! Take into account non-zero umklapps
     ! Some G-G0 won't be avaialable since we use gamma-centered spheres.
     ! In this case the G,G' component is set to zero and this might lead
     ! to wrong results in the interpolation. Fortunately this components 
     ! are located at large G.
     g0 = g0 + vert_g0
     if ( ANY(g0/=0) ) then
       !
       MSG_COMMENT("Had to shift the G-sphere")
       write(std_out,*)g0
       g0_idx = gsph_g_idx(Gsph,g0)
       ABI_CHECK(g0_idx>0,"g0 not found!")

       ABI_ALLOCATE(swap,(npw,npw))
       istat = ABI_ALLOC_STAT
       ABI_CHECK(istat==0," out of memory in swap")

       do iomega=1,nomega
         !
         swap = cwork(:,:,iomega)
         cwork(:,:,iomega) = czero
         do ig2=1,npw
           g2mg0_idx = gsph_gmg_idx(Gsph,ig2,g0_idx) 
           !
           if (g2mg0_idx>0) then
             do ig1=1,npw
               g1mg0_idx = gsph_gmg_idx(Gsph,ig1,g0_idx) 
               if (g1mg0_idx>0) cwork(g1mg0_idx,g2mg0_idx,iomega) = swap(ig1,ig2)
             end do
           end if
         end do
         !
       end do
       !
       ABI_DEALLOCATE(swap)
     end if
     !
     ! Accumulate the contribution of this vertex. 
     vfact = vert_fact(vert,newq,ondiv)

     new_mat(:,:,:,1) = new_mat(:,:,:,1) + vfact*cwork
   end do

   call write_screening(out_unt,accesswff,npw,nomega,new_mat)
 end do ! intq

 if (my_rank==master) close(out_unt)

 ABI_DEALLOCATE(old_mat)
 ABI_DEALLOCATE(new_mat)
 ABI_DEALLOCATE(cwork)
 ABI_DEALLOCATE(qvec)
 !
 ! Free memory
 call free_scrhdr(Old_Hscr)
 call destroy_gsphere(Gsph)
 call destroy_crystal(Cryst)
 call destroy_bz_mesh_type(Old_Qmesh)
 call destroy_bz_mesh_type(New_Qmesh)

contains 

 function vert_fact(vert,newq,ondiv)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vert_fact'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: vert
 real(dp) :: vert_fact
!array
 integer,intent(in) :: ondiv(3)
 real(dp) :: newq(3)

!Local variables-------------------------------
!scalars

! *************************************************************************

!Vxyz = V000 (1 - x) (1 - y) (1 - z) +
!       V100 x (1 - y) (1 - z) + 
!       V010 (1 - x) y (1 - z) + 
!       V110 x y (1 - z) + 
!       V001 (1 - x) (1 - y) z +
!       V101 x (1 - y) z + 
!       V011 (1 - x) y z + 
!       V111 x y z

 vert_fact=one

 write(std_out,*)newq
 write(std_out,*)ondiv

 select case (vert)
 case (1)
 case (2)
 case (3)
 case (4)
 case (5)
 case (6)
 case (7)
 case (8)
 case default
   MSG_ERROR("Wrong vert")
 end select

end function vert_fact
 
end subroutine interpolate_w
!!***

!!****f* m_screening/recalculate_epsm1_freq_grid
!! NAME
!! recalculate_epsm1_freq_grid
!!
!! FUNCTION
!!
!!  Recalculate the frequency gridpoints in the Epsilonm1_results structure.
!!  This is useful when a pole-fit screening is used and the grid can be
!!  arbitrarily chaged for the sigma calculation.
!!
!! INPUTS
!!  Er        - The Epsilonm1_results structure
!!  nfreqre   - Number of real frequencies
!!  nfreqim   - Number of imaginary frequencies
!!  freqremax - Maximum real frequency
!!  freqremin - Minimum real frequency
!!
!! OUTPUT
!!   The Epsilonm1_results structure with changed grid parameters
!!
!! NOTES
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine recalculate_epsm1_freq_grid(Er,nfreqre,nfreqim,freqremin,freqremax,ppmfrq)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'recalculate_epsm1_freq_grid'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfreqre,nfreqim
 real(dp),intent(in) :: freqremin,freqremax,ppmfrq
 type(Epsilonm1_results),intent(inout) :: Er
!arrays

!Local variables-------------------------------
!scalars
 integer :: iomega
 real(dp) :: domegareal
!arrays

! *************************************************************************

 Er%nomega   = nfreqre+nfreqim
 Er%nomega_i = nfreqim
 Er%nomega_r = nfreqre
 Er%nomega_c = 0
 if (associated(Er%omega))  then
   ABI_DEALLOCATE(Er%omega)
 end if
 ABI_ALLOCATE(Er%omega,(Er%nomega))

! Real frequencies
 Er%omega(1)=CMPLX(freqremin,zero,kind=dpc)

 if (Er%nomega_r>1) then ! Avoid division by zero.
   domegareal=(freqremax-freqremin)/(Er%nomega_r-1)
   do iomega=2,Er%nomega_r
     Er%omega(iomega)=CMPLX(freqremin+(iomega-1)*domegareal,zero,kind=dpc)
   end do
 end if

! Imaginary frequencies
 do iomega=1,Er%nomega_i
   Er%omega(Er%nomega_r+iomega)=CMPLX(zero,ppmfrq*third*(EXP(two/(Er%nomega_i+1)*LOG(four)*iomega)-one),kind=dpc)
 end do

end subroutine recalculate_epsm1_freq_grid
!!***

END MODULE m_screening
!!***
