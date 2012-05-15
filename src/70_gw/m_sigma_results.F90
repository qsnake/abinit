!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sigma_results
!! NAME
!!  m_sigma_results
!!
!! FUNCTION
!!  This module provides the definition of the sigma_results data type 
!!  used to store results of the GW calculation as well as as 
!!  methods bound to the object.
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
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_sigma_results

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif
 use m_errors 

 use m_numeric_tools,  only : c2r, r2c
 use m_gwdefs,         only : unt_gw, unt_sig, unt_sgr, unt_sgm, sigma_parameters, sigma_needs_w
 use m_crystal,        only : crystal_structure
 use m_crystal_io,     only : abi_crystal_put, init_crystal_from_hdr
 use m_bz_mesh,        only : bz_mesh_type
 use m_ebands,         only : abi_bands_put, abi_bands_read
 use m_screening,      only : epsilonm1_results

 implicit none

 private 
!!***

!----------------------------------------------------------------------

!!****t* m_sigma_results/sigma_results
!! NAME
!! sigma_results
!!
!! FUNCTION
!! For the GW part of ABINIT, the sigma_results structured datatype
!! gather the results of a sigma calculation.
!!
!! SOURCE

 type,public ::  sigma_results

! WARNING : if you modify this datatype, please check there there is no creation/destruction/copy routine,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: b1gw,b2gw      ! min and Max gw band indeces over spin and k-points (used to dimension)
  integer :: gwcalctyp      ! Flag defining the calculation type.
  integer :: nkptgw         ! No. of points calculated
  integer :: nkibz          ! No. of irreducible k-points.
  integer :: nbnds          ! Total number of bands
  integer :: nomega_r       ! No. of real frequencies for the spectral function.
  integer :: nomega_i       ! No. of frequencies along the imaginary axis.
  integer :: nomega4sd      ! No. of real frequencies to evaluate the derivative of $\Sigma(E)$.
  integer :: nsig_ab        ! 1 if nspinor=1,4 for noncollinear case.
  integer :: nsppol         ! No. of spin polarizations.
  integer :: usepawu        ! 1 if we are using LDA+U as starting point (only for PAW)

  real(dp) :: deltae       ! Frequency step for the calculation of d\Sigma/dE
  real(dp) :: maxomega4sd  ! Max frequency around E_ks for d\Sigma/dE.
  real(dp) :: maxomega_r   ! Max frequency for spectral function.
  real(dp) :: scissor_ene  ! Scissor energy value. zero for None.

  integer,pointer :: maxbnd(:,:)   SET2NULL
  ! maxbnd(nkptgw,nsppol)
  ! Max band index considered in GW for this k-point.

  integer,pointer :: minbnd(:,:)   SET2NULL
  ! minbnd(nkptgw,nsppol)
  ! Min band index considered in GW for this k-point.

  !real(dp),pointer :: ame(:,:,:)
  ! ame(nbnds,nkibz,nomega))
  ! Diagonal matrix elements of the spectral function.
  ! Commented out, it can be calculated from the other quantities

  real(dp),pointer :: degwgap(:,:)   SET2NULL
  ! degwgap(nkibz,nsppol)
  ! Difference btw the QP and the KS optical gap.

  real(dp),pointer :: egwgap(:,:)   SET2NULL
  ! egwgap(nkibz,nsppol))
  ! QP optical gap at each k-point and spin.

  real(dp),pointer :: en_qp_diago(:,:,:)   SET2NULL
  ! en_qp_diago(nbnds,nkibz,nsppol))
  ! QP energies obtained from the diagonalization of the Hermitian approximation to Sigma (QPSCGW)

  real(dp),pointer :: e0(:,:,:)    SET2NULL
  ! e0(nbnds,nkibz,nsppol)
  ! KS eigenvalues for each band, k-point and spin. In case of self-consistent?

  real(dp),pointer :: e0gap(:,:)   SET2NULL
  ! e0gap(nkibz,nsppol),
  ! KS gap at each k-point, for each spin.

  real(dp),pointer :: omega_r(:)   SET2NULL
  ! omega_r(nomega_r)
  ! real frequencies used for the self energy.

  real(dp),pointer :: kptgw(:,:)  SET2NULL
  ! kptgw(3,nkptgw)
  ! ! TODO there is a similar array in sigma_parameters
  ! List of calculated k-points.

  real(dp),pointer :: sigxme(:,:,:)  SET2NULL
  ! sigxme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! Diagonal matrix elements of $\Sigma_x$ i.e $\<nks|\Sigma_x|nks\>$

  real(dp),pointer :: vxcme(:,:,:)  SET2NULL
  ! vxcme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! $\<nks|v_{xc}[n_val]|nks\>$ matrix elements of vxc (valence-only contribution).

  real(dp),pointer :: vUme(:,:,:)   SET2NULL
  ! vUme(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! $\<nks|v_{U}|nks\>$ for LDA+U.

  complex(dpc),pointer :: degw(:,:,:)   SET2NULL
  ! degw(b1gw:b2gw,nkibz,nsppol))
  ! Difference between the QP and the KS energies.

  complex(dpc),pointer :: dsigmee0(:,:,:)  SET2NULL
  ! dsigmee0(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! Derivative of $\Sigma_c(E)$ calculated at the KS eigenvalue.

  complex(dpc),pointer :: egw(:,:,:)  SET2NULL
  ! degw(nbnds,nkibz,nsppol))
  ! QP energies, $\epsilon_{nks}^{QP}$.

  complex(dpc),pointer :: eigvec_qp(:,:,:,:)   SET2NULL
  ! eigvec_qp(nbnds,nbnds,nkibz,nsppol))
  ! Expansion of the QP amplitude in the KS basis set.

  complex(dpc),pointer :: hhartree(:,:,:,:)   SET2NULL
  ! hhartree(b1gw:b2gw,b1gw:b2gw,nkibz,nsppol*nsig_ab)
  ! $\<nks|T+v_H+v_{loc}+v_{nl}|mks\>$

  complex(dpc),pointer :: sigcme(:,:,:,:)   SET2NULL
  ! sigcme(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))
  ! $\<nks|\Sigma_{c}(E)|nks\>$ at each nomega_r frequency

  complex(dpc),pointer :: sigmee(:,:,:)  SET2NULL
  ! sigmee(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! $\Sigma_{xc}E_{KS} + (E_{QP}- E_{KS})*dSigma/dE_KS

  complex(dpc),pointer :: sigcmee0(:,:,:)   SET2NULL
  ! sigcmee0(b1gw:b2gw,nkibz,nsppol*nsig_ab))
  ! Diagonal mat. elements of $\Sigma_c(E)$ calculated at the KS energy $E_{KS}$

  complex(dpc),pointer :: sigcmesi(:,:,:,:)   SET2NULL
  ! sigcmesi(b1gw:b2gw,nkibz,nomega_i,nsppol*nsig_ab))
  ! Matrix elements of $\Sigma_c$ along the imaginary axis.
  ! Only used in case of analytical continuation.

  complex(dpc),pointer :: sigcme4sd(:,:,:,:)   SET2NULL
  ! sigcme4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol*nsig_ab))
  ! Diagonal matrix elements of \Sigma_c around the zeroth order eigenvalue (usually KS).

  complex(dpc),pointer :: sigxcme(:,:,:,:)   SET2NULL
  ! sigxme(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))
  ! $\<nks|\Sigma_{xc}(E)|nks\>$ at each real frequency frequency.

  complex(dpc),pointer :: sigxcmesi(:,:,:,:)   SET2NULL
  ! sigxcmesi(b1gw:b2gw,nkibz,nomega_i,nsppol*nsig_ab))
  ! Matrix elements of $\Sigma_{xc}$ along the imaginary axis.
  ! Only used in case of analytical continuation.

  complex(dpc),pointer :: sigxcme4sd(:,:,:,:)   SET2NULL
  ! sigxcme4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol*nsig_ab))
  ! Diagonal matrix elements of \Sigma_xc for frequencies around the zeroth order eigenvalues.

  complex(dpc),pointer :: ze0(:,:,:)   SET2NULL
  ! ze0(b1gw:b2gw,nkibz,nsppol))
  ! renormalization factor. $(1-\dfrac{\partial\Sigma_c} {\partial E_{KS}})^{-1}$

  complex(dpc),pointer :: omega_i(:)  SET2NULL
  ! omegasi(nomega_i)
  ! Frequencies along the imaginary axis used for the analytical continuation.

  complex(dpc),pointer :: omega4sd(:,:,:,:)  SET2NULL
  ! omega4sd(b1gw:b2gw,nkibz,nomega4sd,nsppol).
  ! Frequencies used to evaluate the Derivative of Sigma.

 end type sigma_results
!!***

!----------------------------------------------------------------------

 public :: write_sigma_results_header
 public :: write_sigma_results      
 public :: print_Sigma_perturbative 
 public :: print_Sigma_QPSC         
 public :: nullify_sigma_results    
 public :: init_sigma_results       
 public :: destroy_sigma_results    
 public :: allocate_sigma_results   
 public :: etsf_dump_QP             
 public :: abi_etsf_get_QP

CONTAINS  !========================================================================================
!!***

!!****f* m_sigma_results/write_sigma_results_header
!! NAME
!! write_sigma_results_header
!!
!! FUNCTION
!!  Write basic info and dimensions used during the calculation 
!!  of the QP correctoions (optdriver==4).
!!
!! INPUTS
!!  Sigp=sigma_parameters
!!  Cryst<Crystal_structure>= Info on the Crystal structure
!!  Kmesh<Bz_mesh_type>= Description of the BZ sampling.
!!
!! OUTPUT
!!  (for writing routines, no output) otherwise, should be described
!!
!! NOTES
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE

subroutine write_sigma_results_header(Sigp,Er,Cryst,Kmesh,Qmesh)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_sigma_results_header'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Bz_mesh_type),intent(in) :: Kmesh,Qmesh
 type(Crystal_structure),intent(in) :: Cryst
 type(Epsilonm1_results),intent(in) :: Er
 type(Sigma_parameters),intent(in) :: Sigp

!Local variables-------------------------------
!scalars
 integer :: mod10
 character(len=500) :: msg

! *************************************************************************

 DBG_ENTER("COLL")

 write(msg,'(a)')' SIGMA fundamental parameters:'
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')

 mod10=MOD(Sigp%gwcalctyp,10)
 SELECT CASE (mod10)
 CASE (0)
  write(msg,'(a,i2)')' PLASMON POLE MODEL ',Sigp%ppmodel
 CASE (1)
  write(msg,'(a)')' ANALYTIC CONTINUATION'
 CASE (2)
  write(msg,'(a)')' CONTOUR DEFORMATION'
 CASE (5)
  write(msg,'(a)')' Hartree-Fock'
 CASE (6)
  write(msg,'(a)')' Screened Exchange'
 CASE (7)
  write(msg,'(a)')' COHSEX'
 CASE (8)
  write(msg,'(a,i2)')' MODEL GW with PLASMON POLE MODEL ',Sigp%ppmodel
 CASE (9)
  write(msg,'(a)')' MODEL GW without PLASMON POLE MODEL'
 CASE DEFAULT
  write(msg,'(a,i3)')' Wrong value for Sigp%gwcalctyp = ',Sigp%gwcalctyp 
  MSG_BUG(msg)
 END SELECT
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')

 write(msg,'(a,i12)')' number of plane-waves for SigmaX         ',Sigp%npwx
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of plane-waves for SigmaC and W   ',Sigp%npwc
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of plane-waves for wavefunctions  ',Sigp%npwwfn
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of bands                          ',Sigp%nbnds
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of independent spin polarizations ',Sigp%nsppol
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of spinorial components           ',Sigp%nspinor
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of k-points in IBZ                ',Kmesh%nibz
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of q-points in IBZ                ',Qmesh%nibz
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of symmetry operations            ',Cryst%nsym
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of k-points in BZ                 ',Kmesh%nbz
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of q-points in BZ                 ',Qmesh%nbz
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of frequencies for dSigma/dE      ',Sigp%nomegasrd
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,f12.2)')' frequency step for dSigma/dE [eV]        ',Sigp%deltae*Ha_eV
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,i12)')' number of omega for Sigma on real axis   ',Sigp%nomegasr
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,f12.2)')' max omega for Sigma on real axis  [eV]   ',Sigp%maxomega_r*Ha_eV
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 write(msg,'(a,f12.2)')' zcut for avoiding poles [eV]             ',Sigp%zcut*Ha_eV
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')

 if (Sigp%soenergy>0.1d-4) then 
   write(msg,'(a,f12.2)')' scissor energy [eV]                      ',Sigp%soenergy*Ha_eV
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 end if

 if (MOD(Sigp%gwcalctyp,10)==1) then
   write(msg,'(a,i12)')' number of imaginary frequencies for Sigma',Sigp%nomegasi
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,f12.2)')' max omega for Sigma on imag axis  [eV]   ',Sigp%omegasimax*Ha_eV
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 end if 

 if (sigma_needs_w(Sigp)) then
   write(msg,'(2a)')ch10,' EPSILON^-1 parameters (SCR file):'
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   !write(std_out,*) titem1(2)(1:79)
   write(msg,'(a,i12)')' dimension of the eps^-1 matrix on file   ',Er%Hscr%npwe
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,i12)')' dimension of the eps^-1 matrix used      ',Er%npwe
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,i12)')' number of plane-waves for wavefunctions  ',Er%Hscr%npwwfn_used
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,i12)')' number of bands                          ',Er%Hscr%nbnds_used
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,i12)')' number of q-points in IBZ                ',Qmesh%nibz
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,i12)')' number of frequencies                    ',Er%nomega
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,i12)')' number of real frequencies               ',Er%nomega_r
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,i12)')' number of imag frequencies               ',Er%nomega_i
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 end if

 write(msg,'(3a)')ch10,' matrix elements of self-energy operator (all in [eV])',ch10
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

 if (Sigp%gwcalctyp<10) then
   write(msg,'(a)')' Perturbative Calculation'
 else if (Sigp%gwcalctyp<20) then
   write(msg,'(a)')' Self-Consistent on Energies only'
 else
   write(msg,'(a)')' Self-Consistent on Energies and Wavefunctions'
 end if
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')

 DBG_EXIT("COLL")

end subroutine write_sigma_results_header
!!***

!----------------------------------------------------------------------

!!****f* m_sigma_results/write_sigma_results
!! NAME
!! write_sigma_results
!!
!! FUNCTION
!!  Write the final results of the GW calculation.
!!
!! INPUTS
!!  KS_BSt<Bandstructure_type>=Info on the KS band structure energies.
!!     %eig(mband,nkibz,nsppol)= KS energies
!!  ikibz= index of the k-point in the array kibz, where GW corrections are calculated 
!!  ikcalc= index of the k-point in the array Sigp%kptgw2bz
!!  Sigp=sigma_parameters datatype
!!  sr=sigma results datatype
!!
!! OUTPUT
!!  (for writing routines, no output) otherwise, should be described
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE
!!

subroutine write_sigma_results(ikcalc,ikibz,Sigp,Sr,KS_BSt)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_sigma_results'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikcalc,ikibz
 type(Bandstructure_type),intent(in) :: KS_BSt
 type(Sigma_parameters),intent(in) :: Sigp
 type(sigma_results),intent(in) :: Sr
!arrays

!Local variables-------------------------------
!scalars
 integer :: ib,io,is
 character(len=500) :: msg
!arrays
 character(len=12) :: tag_spin(2)

! *************************************************************************

 !unt_gw  File with GW corrections.
 !unt_sig Self-energy as a function of frequency.
 !unt_sgr Derivative wrt omega of the Self-energy.
 !unt_sgm Sigma on the Matsubara axis.

 tag_spin=(/'            ','            '/); if (Sr%nsppol==2) tag_spin=(/',  SPIN UP  ',',  SPIN DOWN'/)

 do is=1,Sr%nsppol
   write(msg,'(2a,3f8.3,a)')ch10,' k = ',Sigp%kptgw(:,ikcalc),tag_spin(is)
   call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

   msg = ' Band     E0 <VxcLDA>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E'
   if (Sr%usepawu/=0) then
     msg = ' Band     E0 <VxcLDA>   <H_U>  SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E'
   end if

   if (Sigp%gwcalctyp>=10) then
     write(msg,'(2a)')&
&     ' Band     E_lda   <Vxclda>   E(N-1)  <Hhartree>   SigX  SigC[E(N-1)]',&
&     '    Z     dSigC/dE  Sig[E(N)]  DeltaE  E(N)_pert E(N)_diago'
   end if
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')

   write(unt_gw,'(3f10.6)')Sigp%kptgw(:,ikcalc)
   write(unt_gw,'(i4)')Sigp%maxbnd(ikcalc,is)-Sigp%minbnd(ikcalc,is)+1

   write(687,'(3f10.6)')Sigp%kptgw(:,ikcalc)
   write(687,'(i4)')Sigp%maxbnd(ikcalc,is)-Sigp%minbnd(ikcalc,is)+1

   write(unt_sig,'("# k = ",3f10.6)')Sigp%kptgw(:,ikcalc)
   write(unt_sig,'("# b = ",2i10)')Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)

   write(unt_sgr,'("# k = ",3f10.6)')Sigp%kptgw(:,ikcalc)
   write(unt_sgr,'("# b = ",2i10)')Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)

   do ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)
     if (Sigp%gwcalctyp>=10) then
       call print_Sigma_QPSC(Sr,ikibz,ib,is,KS_BSt,unit=ab_out)
       call print_Sigma_QPSC(Sr,ikibz,ib,is,KS_BSt,unit=std_out,prtvol=1)
      
       write(687,'(i6,3f9.4)')          &
       &      ib,                               &
       &      Sr%en_qp_diago(ib,ikibz,is)*Ha_eV,&
       &      (Sr%en_qp_diago(ib,ikibz,is) - KS_BSt%eig(ib,ikibz,is))*Ha_eV,&
       &      zero

     else
       ! If not ppmodel, write out also the imaginary part in ab_out
       SELECT CASE(Sigp%gwcalctyp)
       CASE(1,2)
         call print_Sigma_perturbative(Sr,ikibz,ib,is,unit=ab_out,prtvol=1)
       CASE DEFAULT
         call print_Sigma_perturbative(Sr,ikibz,ib,is,unit=ab_out)
       END SELECT
       call print_Sigma_perturbative(Sr,ikibz,ib,is,unit=std_out,prtvol=1)
     end if

     write(unt_gw,'(i6,3f9.4)')          &
&      ib,                               &
&      REAL (Sr%egw (ib,ikibz,is))*Ha_eV,&
&      REAL (Sr%degw(ib,ikibz,is))*Ha_eV,&
&      AIMAG(Sr%egw (ib,ikibz,is))*Ha_eV
   end do !ib

   if (Sr%e0gap(ikibz,is)**2+Sr%egwgap(ikibz,is)**2+Sr%degwgap(ikibz,is)**2 > tol10) then
     ! Output the direct gap for each spin
     ! If all the gaps are zero, this means that it could not be computed in the calling routine
     write(msg,'(2a,f8.3)')ch10,' E^0_gap       ',Sr%e0gap(ikibz,is)*Ha_eV
     call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
     write(msg,'(a,f8.3)')      ' E^GW_gap      ',Sr%egwgap(ikibz,is)*Ha_eV
     call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
     write(msg,'(a,f8.3,a)')    ' DeltaE^GW_gap ',Sr%degwgap(ikibz,is)*Ha_eV,ch10
     call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
   end if
   !
   ! === Output of the spectral function ===
   do io=1,Sr%nomega_r
     write(unt_sig,'(100(e12.5,2x))')&
&     REAL(Sr%omega_r(io))*Ha_eV,&
&     (REAL(Sr%sigxcme(ib,ikibz,io,is))*Ha_eV,&
&     AIMAG(Sr%sigxcme(ib,ikibz,io,is))*Ha_eV,&
&     one/pi*ABS(AIMAG(Sr%sigcme(ib,ikibz,io,is)))&
&     /( (REAL(Sr%omega_r(io)-Sr%hhartree(ib,ib,ikibz,is)-Sr%sigxcme(ib,ikibz,io,is)))**2&
&       +(AIMAG(Sr%sigcme(ib,ikibz,io,is)))**2) /Ha_eV,&
&     ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is))
   end do
   !
   do ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)
     write(unt_sgr,'("# ik, ib",2i5)')ikibz,ib
     do io=1,Sr%nomega4sd
       write(unt_sgr,'(100(e12.5,2x))')             &
&        REAL (Sr%omega4sd  (ib,ikibz,io,is)) *Ha_eV,&
&        REAL (Sr%sigxcme4sd(ib,ikibz,io,is)) *Ha_eV,&
&        AIMAG(Sr%sigxcme4sd(ib,ikibz,io,is)) *Ha_eV
     end do
   end do
   !
   if (MOD(Sigp%gwcalctyp,10)==1) then ! For AC, write sigma matrix elements along the imaginary axis
     do ib=Sigp%minbnd(ikcalc,is),Sigp%maxbnd(ikcalc,is)
       write(unt_sgm,'("# ik, ib",2i5)')ikibz,ib
       do io=1,Sr%nomega_i
         write(unt_sgm,'(3(e12.5,2x))')              &
&          AIMAG(Sr%omega_i(io))              *Ha_eV,&
&          REAL (Sr%sigxcmesi(ib,ikibz,io,is))*Ha_eV,&
&          AIMAG(Sr%sigxcmesi(ib,ikibz,io,is))*Ha_eV
       end do
     end do
   end if 

 end do !is

end subroutine write_sigma_results
!!***

!----------------------------------------------------------------------

!!****f* m_sigma_results/print_Sigma_perturbative 
!! NAME
!! print_Sigma_perturbative
!!
!! FUNCTION
!!  Write the results of the GW calculation done with the perturbative approach
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwannier,m_sigma_results
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE

subroutine print_Sigma_perturbative(Sr,ik_ibz,iband,isp,unit,prtvol,mode_paral,witheader)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_Sigma_perturbative'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,ik_ibz,isp
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 logical,optional,intent(in) :: witheader
 type(sigma_results),intent(in) :: Sr

!Local variables-------------------------------
!scalars
 integer :: my_unt,verbose
 character(len=4) :: my_mode
 character(len=500) :: msg

! *********************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 verbose=0      ; if (PRESENT(prtvol    )) verbose=prtvol
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

 if (PRESENT(witheader)) then 
   if (witheader) call wrtout(my_unt,' Band     E0 <VxcLDA>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E ',my_mode) 
 end if

 if (Sr%usepawu==0) then 

   if (Sr%nsig_ab/=1) then
     write(msg,'(i5,9f8.3)')                        & 
&           iband,                                  &
&           Sr%e0          (iband,ik_ibz,1)*Ha_eV,  &
&           SUM(Sr%vxcme   (iband,ik_ibz,:))*Ha_eV, &
&           SUM(Sr%sigxme  (iband,ik_ibz,:))*Ha_eV, &
&      REAL(SUM(Sr%sigcmee0(iband,ik_ibz,:)))*Ha_eV,&
&      REAL(Sr%ze0         (iband,ik_ibz,1)),       &
&      REAL(SUM(Sr%dsigmee0(iband,ik_ibz,:))),      &
&      REAL(SUM(Sr%sigmee  (iband,ik_ibz,:)))*Ha_eV,&
&      REAL(Sr%degw        (iband,ik_ibz,1))*Ha_eV, &
&      REAL(Sr%egw         (iband,ik_ibz,1))*Ha_eV
       call wrtout(my_unt,msg,my_mode) 
     if (verbose/=0) then
       write(msg,'(i5,9f8.3)')                         & 
&              iband,                                  &
&              zero,                                   &
&              zero,                                   &
&              zero,                                   &
&        AIMAG(SUM(Sr%sigcmee0(iband,ik_ibz,:)))*Ha_eV,&
&        AIMAG(Sr%ze0         (iband,ik_ibz,1)),       &
&        AIMAG(SUM(Sr%dsigmee0(iband,ik_ibz,:))),      &
&        AIMAG(SUM(Sr%sigmee  (iband,ik_ibz,:)))*Ha_eV,&
&        AIMAG(Sr%degw        (iband,ik_ibz,1))*Ha_eV, &
&        AIMAG(Sr%egw         (iband,ik_ibz,1))*Ha_eV
       call wrtout(my_unt,msg,my_mode) 
     end if
  else
    write(msg,'(i5,9f8.3)')                     & 
&          iband,                               &
&          Sr%e0      (iband,ik_ibz,isp)*Ha_eV, &
&          Sr%vxcme   (iband,ik_ibz,isp)*Ha_eV, &
&          Sr%sigxme  (iband,ik_ibz,isp)*Ha_eV, &
&     REAL(Sr%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
&     REAL(Sr%ze0     (iband,ik_ibz,isp)),      &
&     REAL(Sr%dsigmee0(iband,ik_ibz,isp)),      &
&     REAL(Sr%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
&     REAL(Sr%degw    (iband,ik_ibz,isp))*Ha_eV,&
&     REAL(Sr%egw     (iband,ik_ibz,isp))*Ha_eV
    call wrtout(my_unt,msg,my_mode) 

    if (verbose/=0) then
      write(msg,'(i5,9f8.3)')                       & 
&              iband,                               &
&              zero,                                &
&              zero,                                &
&              zero,                                &
&        AIMAG(Sr%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
&        AIMAG(Sr%ze0     (iband,ik_ibz,isp)),      &
&        AIMAG(Sr%dsigmee0(iband,ik_ibz,isp)),      &
&        AIMAG(Sr%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
&        AIMAG(Sr%degw    (iband,ik_ibz,isp))*Ha_eV,&
&        AIMAG(Sr%egw     (iband,ik_ibz,isp))*Ha_eV
       call wrtout(my_unt,msg,my_mode) 
    end if
  end if

 else  ! PAW+U+GW calculation.
   ABI_CHECK(Sr%nsig_ab==1,'LDA+U with spinor not implemented')
   write(msg,'(i5,10f8.3)')                    & 
&         iband,                               &
&         Sr%e0      (iband,ik_ibz,isp)*Ha_eV, &
&         Sr%vxcme   (iband,ik_ibz,isp)*Ha_eV, &
&         Sr%vUme    (iband,ik_ibz,isp)*Ha_eV, &
&         Sr%sigxme  (iband,ik_ibz,isp)*Ha_eV, &
&    REAL(Sr%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
&    REAL(Sr%ze0     (iband,ik_ibz,isp)),      &
&    REAL(Sr%dsigmee0(iband,ik_ibz,isp)),      &
&    REAL(Sr%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
&    REAL(Sr%degw    (iband,ik_ibz,isp))*Ha_eV,&
&    REAL(Sr%egw     (iband,ik_ibz,isp))*Ha_eV
   call wrtout(my_unt,msg,my_mode) 

   if (verbose/=0) then
     write(msg,'(i5,10f8.3)')                    & 
&           iband,                               &
&           zero,                                &
&           zero,                                &
&           zero,                                &
&           zero,                                &
&     AIMAG(Sr%sigcmee0(iband,ik_ibz,isp))*Ha_eV,&
&     AIMAG(Sr%ze0     (iband,ik_ibz,isp)),      &
&     AIMAG(Sr%dsigmee0(iband,ik_ibz,isp)),      &
&     AIMAG(Sr%sigmee  (iband,ik_ibz,isp))*Ha_eV,&
&     AIMAG(Sr%degw    (iband,ik_ibz,isp))*Ha_eV,&
&     AIMAG(Sr%egw     (iband,ik_ibz,isp))*Ha_eV
      call wrtout(my_unt,msg,my_mode) 
   end if
 end if

end subroutine print_Sigma_perturbative 
!!***

!----------------------------------------------------------------------

!!****f* m_sigma_results/print_Sigma_QPSC
!! NAME
!!  print_Sigma_QPSC
!!
!! FUNCTION
!!  Write the results of the GW calculation in case of self-consistency
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_sigma_results
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE

subroutine print_Sigma_QPSC(Sr,ik_ibz,iband,isp,KS_BSt,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_Sigma_QPSC'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,ik_ibz,isp
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(sigma_results),intent(in) :: Sr
 type(Bandstructure_type),intent(in) :: KS_BSt
!arrays

!Local variables-------------------------------
!scalars
 integer :: my_unt,verbose
 character(len=4) :: my_mode
 character(len=500) :: msg

! *********************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 verbose=0      ; if (PRESENT(prtvol    )) verbose=prtvol
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

! write(msg,'(a)')&
!&   ' Band     E_lda   <Vxclda>   E(N-1)  <Hhartree>   SigX  SigC[E(N-1)]',&
!&   '    Z     dSigC/dE  Sig[E(N)]  DeltaE  E(N)_pert E(N)_diago'

 if (Sr%usepawu==0 .or. .TRUE.) then
   if (Sr%nsig_ab/=1) then
     write(msg,'(i5,12(2x,f8.3))')                        & 
&           iband,                                        &
&           KS_BSt%eig     (iband,ik_ibz,1)*Ha_eV,        &
&           SUM(Sr%vxcme   (iband,ik_ibz,:))*Ha_eV,       &
&           Sr%e0          (iband,ik_ibz,1)*Ha_eV,        &
&      REAL(SUM(Sr%hhartree(iband,iband,ik_ibz,:)))*Ha_eV,&
&           SUM(Sr%sigxme  (iband,ik_ibz,:))*Ha_eV,       &
&      REAL(SUM(Sr%sigcmee0(iband,ik_ibz,:)))*Ha_eV,      &
&      REAL(Sr%ze0         (iband,ik_ibz,1)),             &
&      REAL(SUM(Sr%dsigmee0(iband,ik_ibz,:))),            &
&      REAL(SUM(Sr%sigmee  (iband,ik_ibz,:)))*Ha_eV,      &
&      REAL(Sr%degw        (iband,ik_ibz,1))*Ha_eV,       &
&      REAL(Sr%egw         (iband,ik_ibz,1))*Ha_eV,       &
&           Sr%en_qp_diago (iband,ik_ibz,1)*Ha_eV
     call wrtout(my_unt,msg,my_mode) 

     write(msg,'(i5,12(2x,f8.3))')                         & 
&            iband,                                        &
&            zero,                                         &
&            zero,                                         &
&            zero,                                         &
&      AIMAG(SUM(Sr%hhartree(iband,iband,ik_ibz,:)))*Ha_eV,&
&            zero,                                         &
&      AIMAG(SUM(Sr%sigcmee0(iband,ik_ibz,:)))*Ha_eV,      &
&      AIMAG(Sr%ze0         (iband,ik_ibz,1)),             &
&      AIMAG(SUM(Sr%dsigmee0(iband,ik_ibz,:))),            &
&      AIMAG(SUM(Sr%sigmee  (iband,ik_ibz,:)))*Ha_eV,      &
&      AIMAG(Sr%degw        (iband,ik_ibz,1))*Ha_eV,       &
&      AIMAG(Sr%egw         (iband,ik_ibz,1))*Ha_eV,       &
&            zero
     if (verbose/=0) call wrtout(my_unt,msg,my_mode) 
   else
     write(msg,'(i5,12(2x,f8.3))')                        & 
&           iband,                                        &
&           KS_BSt%eig    (iband,ik_ibz,isp)*Ha_eV,       &
&           Sr%vxcme      (iband,ik_ibz,isp)*Ha_eV,       &
&           Sr%e0         (iband,ik_ibz,isp)*Ha_eV,       &
&      REAL(Sr%hhartree   (iband,iband,ik_ibz,isp))*Ha_eV,&
&           Sr%sigxme     (iband,ik_ibz,isp)*Ha_eV,       &
&      REAL(Sr%sigcmee0   (iband,ik_ibz,isp))*Ha_eV,      &
&      REAL(Sr%ze0        (iband,ik_ibz,isp)),            &
&      REAL(Sr%dsigmee0   (iband,ik_ibz,isp)),            &
&      REAL(Sr%sigmee     (iband,ik_ibz,isp))*Ha_eV,      &
&      REAL(Sr%degw       (iband,ik_ibz,isp))*Ha_eV,      &
&      REAL(Sr%egw        (iband,ik_ibz,isp))*Ha_eV,      &
&           Sr%en_qp_diago(iband,ik_ibz,isp)*Ha_eV
     call wrtout(my_unt,msg,my_mode) 

     write(msg,'(i5,12(2x,f8.3))')                        & 
&            iband,                                       &
&            zero,                                        &
&            zero,                                        &
&            zero,                                        &
&      AIMAG(Sr%hhartree  (iband,iband,ik_ibz,isp))*Ha_eV,&
&            zero,                                        &
&      AIMAG(Sr%sigcmee0   (iband,ik_ibz,isp))*Ha_eV,     &
&      AIMAG(Sr%ze0        (iband,ik_ibz,isp)),           &
&      AIMAG(Sr%dsigmee0   (iband,ik_ibz,isp)),           &
&      AIMAG(Sr%sigmee     (iband,ik_ibz,isp))*Ha_eV,     &
&      AIMAG(Sr%degw       (iband,ik_ibz,isp))*Ha_eV,     &
&      AIMAG(Sr%egw        (iband,ik_ibz,isp))*Ha_eV,     &
&            zero
     if (verbose/=0) call wrtout(my_unt,msg,my_mode) 
   end if

 else ! PAW+U+GW calculation.
   MSG_ERROR("PAW+U+GW not yet implemented")
 end if

end subroutine print_Sigma_QPSC
!!***

!----------------------------------------------------------------------

!!****f* m_sigma_results/nullify_sigma_results
!! NAME
!! nullify_sigma_results
!!
!! FUNCTION
!!  Initialize all pointers defined in the sigma_results data type to NULL.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_sigma_results
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE

subroutine nullify_sigma_results(Sr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_sigma_results'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(sigma_results),intent(inout) :: Sr

! *************************************************************************

 !@sigma_results

!integer 
 nullify(Sr%maxbnd)
 nullify(Sr%minbnd)

!real
 !nullify(Sr%ame       )
 nullify(Sr%degwgap    )
 nullify(Sr%egwgap     )
 nullify(Sr%en_qp_diago)
 nullify(Sr%e0         )
 nullify(Sr%e0gap      )
 nullify(Sr%omega_r    )
 nullify(Sr%kptgw      )
 nullify(Sr%sigxme     )
 nullify(Sr%vxcme      )
 nullify(Sr%vUme       )

!complex
 nullify(Sr%degw      )
 nullify(Sr%dsigmee0  )
 nullify(Sr%egw       )
 nullify(Sr%eigvec_qp )
 nullify(Sr%hhartree  )
 nullify(Sr%sigcme    )
 nullify(Sr%sigmee    )
 nullify(Sr%sigcmee0  )
 nullify(Sr%sigcmesi  )
 nullify(Sr%sigcme4sd )
 nullify(Sr%sigxcme   )
 nullify(Sr%sigxcmesi )
 nullify(Sr%sigxcme4sd)
 nullify(Sr%ze0       )
 nullify(Sr%omega_i   )
 nullify(Sr%omega4sd  )

end subroutine nullify_sigma_results
!!***

!----------------------------------------------------------------------

!!****f* m_sigma_results/init_sigma_results
!! NAME
!! init_sigma_results
!!
!! FUNCTION
!! Main creation method for the sigma_results data type.
!!
!! INPUTS
!! usepawu=1 if we used LDA+U as starting point (only for PAW)
!!
!! OUTPUT
!!
!! TODO
!!  Write documentation.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE

subroutine init_sigma_results(Sigp,nkibz,usepawu,Sr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_sigma_results'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nkibz,usepawu
!scalars
 type(Sigma_parameters),intent(in) :: Sigp
 type(sigma_results),intent(inout) :: Sr

!Local variables-------------------------------
!scalars
 integer :: b1gw,b2gw,mod10

! *************************************************************************

 !@sigma_results
 call nullify_sigma_results(Sr)

 ! === Copy important dimensions ===
 mod10=MOD(Sigp%gwcalctyp,10)

!BEGIN NEW
 Sr%nkptgw     =Sigp%nkptgw
 Sr%gwcalctyp  =Sigp%gwcalctyp
 Sr%deltae     =Sigp%deltae
 Sr%maxomega4sd=Sigp%maxomega4sd
 Sr%maxomega_r =Sigp%maxomega_r
 Sr%scissor_ene=Sigp%soenergy
 !FIXME this should be done in allocate_sigma_results
 ABI_ALLOCATE(Sr%minbnd,(Sr%nkptgw,Sigp%nsppol))
 ABI_ALLOCATE(Sr%maxbnd,(Sr%nkptgw,Sigp%nsppol))
 Sr%minbnd=Sigp%minbnd 
 Sr%maxbnd=Sigp%maxbnd
 ABI_ALLOCATE(Sr%kptgw,(3,Sr%nkptgw))
 Sr%kptgw=Sigp%kptgw
!END NEW

 Sr%b1gw     =Sigp%minbdgw ! * min and Max GW band index over k and spin. 
 Sr%b2gw     =Sigp%maxbdgw !   Used to dimension arrays.
 Sr%nbnds    =Sigp%nbnds
 Sr%nkibz    =nkibz
 Sr%nsppol   =Sigp%nsppol
 Sr%nsig_ab  =Sigp%nsig_ab
 Sr%nomega_r =Sigp%nomegasr  !FIXME change name
 Sr%nomega_i =Sigp%nomegasi
 Sr%nomega4sd=Sigp%nomegasrd
 Sr%usepawu  =usepawu

 !======================================================
 ! === Allocate arrays in the sigma_results datatype ===
 !======================================================
 b1gw=Sr%b1gw  
 b2gw=Sr%b2gw   

 !TODO use this routine
! call allocate_sigma_results(Sr,b1gw,b2gw,Sr%nbnds,Sr%nkibz,Sr%nsppol,&
!& Sr%nsig_ab,Sr%nomega_r,Sr%nomega_i,Sr%nomega4sd,omega_r=Sigp%omega_r,omega_i=Sigp%omegasi)

 ! hhartree(b1,b2,k,s)= <b1,k,s|T+v_{loc}+v_{nl}+v_{H}|b2,k,s>
 ABI_ALLOCATE(Sr%hhartree,(b1gw:b2gw,b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 Sr%hhartree=czero

 ! === QP amplitudes and energies ===
 ABI_ALLOCATE(Sr%en_qp_diago,(Sr%nbnds,Sr%nkibz,Sr%nsppol))
 Sr%en_qp_diago=zero
 ABI_ALLOCATE(Sr%eigvec_qp,(Sr%nbnds,Sr%nbnds,Sr%nkibz,Sr%nsppol))
 Sr%eigvec_qp  =czero

 ! Dont know if it is better to do this here or in the sigma
 ! * Initialize with KS wavefunctions and energies
 !do ib=1,Sr%nbnds
 ! Sr%en_qp_diago(ib,:,:)=en(:,ib,:)
 ! Sr%eigvec_qp(ib,ib,:,:)=cone
 !end do 

 ABI_ALLOCATE(Sr%vxcme   ,(b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_ALLOCATE(Sr%vUme    ,(b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_ALLOCATE(Sr%sigxme  ,(b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))

 ABI_ALLOCATE(Sr%sigcme  ,(b1gw:b2gw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))
 ABI_ALLOCATE(Sr%sigxcme ,(b1gw:b2gw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))

 ABI_ALLOCATE(Sr%sigcmee0,(b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_ALLOCATE(Sr%ze0     ,(b1gw:b2gw,Sr%nkibz,Sr%nsppol))
 ABI_ALLOCATE(Sr%dsigmee0,(b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_ALLOCATE(Sr%sigmee  ,(b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 ABI_ALLOCATE(Sr%degw    ,(b1gw:b2gw,Sr%nkibz,Sr%nsppol))

 ABI_ALLOCATE(Sr%e0 ,(Sr%nbnds,Sr%nkibz,Sr%nsppol))
 ABI_ALLOCATE(Sr%egw,(Sr%nbnds,Sr%nkibz,Sr%nsppol))

 ABI_ALLOCATE(Sr%e0gap  ,(Sr%nkibz,Sr%nsppol))
 ABI_ALLOCATE(Sr%degwgap,(Sr%nkibz,Sr%nsppol))
 ABI_ALLOCATE(Sr%egwgap ,(Sr%nkibz,Sr%nsppol))
 !allocate(Sr%ame(Sr%nbnds,Sr%nkibz,Sr%nomega_r))
 !
 ! === These quantities are used to evaluate $\Sigma(E)$ around the KS\QP eigenvalue ===
 ABI_ALLOCATE(Sr%omega4sd  ,(b1gw:b2gw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol))
 ABI_ALLOCATE(Sr%sigcme4sd ,(b1gw:b2gw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))
 ABI_ALLOCATE(Sr%sigxcme4sd,(b1gw:b2gw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))

 !TODO Find  better treatment
 ! Mesh along the real axis.
 if (Sr%nomega_r>0) then
   ABI_ALLOCATE(Sr%omega_r,(Sr%nomega_r))
   Sr%omega_r(:)=Sigp%omega_r(:)
 end if

 Sr%e0        =zero
 Sr%egw       =czero
 Sr%e0gap     =zero
 Sr%sigcme    =czero
 Sr%sigxme    =czero
 Sr%sigxcme   =czero
 Sr%sigcmee0  =czero
 Sr%ze0       =czero
 Sr%dsigmee0  =czero
 Sr%sigmee    =czero
 Sr%omega4sd  =czero
 Sr%sigcme4sd =czero
 Sr%sigxcme4sd=czero
 Sr%degw      =czero

 ! === Analytical Continuation ===
 if (mod10==1) then 
   ! FIXME omegasi should not be in Sigp% here we should construct the mesh
   ABI_ALLOCATE(Sr%omega_i,(Sr%nomega_i))
   Sr%omega_i=Sigp%omegasi
   ABI_ALLOCATE(Sr%sigcmesi ,(b1gw:b2gw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
   ABI_ALLOCATE(Sr%sigxcmesi,(b1gw:b2gw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
   Sr%sigcmesi =czero
   Sr%sigxcmesi=czero
 end if

end subroutine init_sigma_results
!!***

!----------------------------------------------------------------------

!!****f* m_sigma_results/destroy_sigma_results
!! NAME
!! destroy_sigma_results
!!
!! FUNCTION
!!  Deallocate all associated pointers defined in the sigma_results data type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwannier,sigma
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE

subroutine destroy_sigma_results(Sr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_sigma_results'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(sigma_results),intent(inout) :: Sr

! *************************************************************************

 DBG_ENTER("COLL")

 !@sigma_results

!integer
 if (associated(Sr%maxbnd     ))     then
   ABI_DEALLOCATE(Sr%maxbnd)
 end if
 if (associated(Sr%minbnd     ))     then
   ABI_DEALLOCATE(Sr%minbnd)
 end if

!real
 !if (associated(Sr%ame       ))    deallocate(Sr%ame) 
 if (associated(Sr%degwgap    ))     then
   ABI_DEALLOCATE(Sr%degwgap)
 end if
 if (associated(Sr%egwgap     ))     then
   ABI_DEALLOCATE(Sr%egwgap)
 end if
 if (associated(Sr%en_qp_diago))     then
   ABI_DEALLOCATE(Sr%en_qp_diago)
 end if
 if (associated(Sr%e0         ))     then
   ABI_DEALLOCATE(Sr%e0)
 end if
 if (associated(Sr%e0gap      ))     then
   ABI_DEALLOCATE(Sr%e0gap)
 end if
 if (associated(Sr%omega_r    ))     then
   ABI_DEALLOCATE(Sr%omega_r)
 end if
 if (associated(Sr%kptgw      ))     then
   ABI_DEALLOCATE(Sr%kptgw)
 end if
 if (associated(Sr%sigxme     ))     then
   ABI_DEALLOCATE(Sr%sigxme)
 end if
 if (associated(Sr%vxcme      ))     then
   ABI_DEALLOCATE(Sr%vxcme)
 end if
 if (associated(Sr%vUme       ))     then
   ABI_DEALLOCATE(Sr%vUme)
 end if
 
!complex
 if (associated(Sr%degw       ))     then
   ABI_DEALLOCATE(Sr%degw)
 end if
 if (associated(Sr%dsigmee0   ))     then
   ABI_DEALLOCATE(Sr%dsigmee0)
 end if
 if (associated(Sr%egw        ))     then
   ABI_DEALLOCATE(Sr%egw)
 end if
 if (associated(Sr%eigvec_qp  ))     then
   ABI_DEALLOCATE(Sr%eigvec_qp)
 end if
 if (associated(Sr%hhartree   ))     then
   ABI_DEALLOCATE(Sr%hhartree)
 end if
 if (associated(Sr%sigcme     ))     then
   ABI_DEALLOCATE(Sr%sigcme)
 end if
 if (associated(Sr%sigmee     ))     then
   ABI_DEALLOCATE(Sr%sigmee)
 end if
 if (associated(Sr%sigcmee0   ))     then
   ABI_DEALLOCATE(Sr%sigcmee0)
 end if
 if (associated(Sr%sigcmesi   ))     then
   ABI_DEALLOCATE(Sr%sigcmesi)
 end if
 if (associated(Sr%sigcme4sd  ))     then
   ABI_DEALLOCATE(Sr%sigcme4sd)
 end if
 if (associated(Sr%sigxcme    ))     then
   ABI_DEALLOCATE(Sr%sigxcme)
 end if
 if (associated(Sr%sigxcmesi  ))     then
   ABI_DEALLOCATE(Sr%sigxcmesi)
 end if
 if (associated(Sr%sigxcme4sd ))     then
   ABI_DEALLOCATE(Sr%sigxcme4sd)
 end if
 if (associated(Sr%ze0        ))     then
   ABI_DEALLOCATE(Sr%ze0)
 end if
 if (associated(Sr%omega_i    ))     then
   ABI_DEALLOCATE(Sr%omega_i)
 end if
 if (associated(Sr%omega4sd   ))     then
   ABI_DEALLOCATE(Sr%omega4sd)
 end if

 DBG_EXIT("COLL")

end subroutine destroy_sigma_results
!!***

!----------------------------------------------------------------------

!!****f* m_sigma_results/allocate_sigma_results
!! NAME
!! allocate_sigma_results
!!
!! FUNCTION
!!  Allocates the dynamic arrays in the sigma_results data type starting 
!!  from the knowledge of the basic dimensions used to calculate the QP corrections.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_sigma_results
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE

subroutine allocate_sigma_results(Sr,b1gw,b2gw,nbnds,nkibz,nkptgw,nsppol,nsig_ab,nomega_r,nomega_i,nomega4sd,&
& omega_r,omega_i) ! Optional

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'allocate_sigma_results'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: b1gw,b2gw,nkibz,nsppol,nsig_ab,nbnds
 integer,intent(in) :: nomega_r,nomega_i,nomega4sd,nkptgw
 type(sigma_results),intent(inout) :: Sr
!arrays
 complex(dpc),optional,intent(in) :: omega_r(:),omega_i(:)

!Local variables-------------------------------
! *************************************************************************

 !@sigma_results

 call nullify_sigma_results(Sr)

 Sr%nkptgw=nkptgw
 ABI_ALLOCATE(Sr%minbnd,(Sr%nkptgw,nsppol))
 ABI_ALLOCATE(Sr%maxbnd,(Sr%nkptgw,nsppol))

 ABI_ALLOCATE(Sr%kptgw,(3,Sr%nkptgw))

 ! hhartree(b1,b2,k,s)= <b1,k,s|T+v_{loc}+v_{nl}+v_{H}|b2,k,s>
 ABI_ALLOCATE(Sr%hhartree,(b1gw:b2gw,b1gw:b2gw,nkibz,nsppol*nsig_ab))
 Sr%hhartree=czero

 ! === QP amplitudes and energies ===
 ABI_ALLOCATE(Sr%en_qp_diago,(nbnds,nkibz,nsppol))
 ABI_ALLOCATE(Sr%eigvec_qp,(nbnds,nbnds,nkibz,nsppol))
 Sr%en_qp_diago=zero
 Sr%eigvec_qp  =czero

 ABI_ALLOCATE(Sr%vxcme   ,(b1gw:b2gw,nkibz,nsppol*nsig_ab))
 ABI_ALLOCATE(Sr%vUme    ,(b1gw:b2gw,nkibz,nsppol*nsig_ab))
 ABI_ALLOCATE(Sr%sigxme  ,(b1gw:b2gw,nkibz,nsppol*nsig_ab))

 ABI_ALLOCATE(Sr%sigcme  ,(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))
 ABI_ALLOCATE(Sr%sigxcme ,(b1gw:b2gw,nkibz,nomega_r,nsppol*nsig_ab))

 ABI_ALLOCATE(Sr%sigcmee0,(b1gw:b2gw,nkibz,nsppol*nsig_ab))
 ABI_ALLOCATE(Sr%ze0     ,(b1gw:b2gw,nkibz,nsppol))
 ABI_ALLOCATE(Sr%dsigmee0,(b1gw:b2gw,nkibz,nsppol*nsig_ab))
 ABI_ALLOCATE(Sr%sigmee  ,(b1gw:b2gw,nkibz,nsppol*nsig_ab))
 ABI_ALLOCATE(Sr%degw    ,(b1gw:b2gw,nkibz,nsppol))

 ABI_ALLOCATE(Sr%e0 ,(nbnds,nkibz,nsppol))
 ABI_ALLOCATE(Sr%egw,(nbnds,nkibz,nsppol))

 ABI_ALLOCATE(Sr%e0gap  ,(nkibz,nsppol))
 ABI_ALLOCATE(Sr%degwgap,(nkibz,nsppol))
 ABI_ALLOCATE(Sr%egwgap ,(nkibz,nsppol))
 !allocate(Sr%ame(nbnds,nkibz,nomega_r))

 ! === These quantities are used to evaluate $\Sigma(E)$ around the KS\QP eigenvalue ===
 ABI_ALLOCATE(Sr%omega4sd  ,(b1gw:b2gw,nkibz,nomega4sd,nsppol))
 ABI_ALLOCATE(Sr%sigcme4sd ,(b1gw:b2gw,nkibz,nomega4sd,nsppol*nsig_ab))
 ABI_ALLOCATE(Sr%sigxcme4sd,(b1gw:b2gw,nkibz,nomega4sd,nsppol*nsig_ab))
 
 if (nomega_r>0) then ! Mesh along the real axis.
   ABI_ALLOCATE(Sr%omega_r,(nomega_r))
   if (PRESENT(omega_r)) then 
     ABI_CHECK(SIZE(omega_r)==SIZE(Sr%omega_r),'DIM omega_r=/Sr%omega_r')
     Sr%omega_r(:)=omega_r(:)
   end if
 end if

 ! === Analytical Continuation ===
 !if (mod10==1) then 
 if (nomega_i>0) then
   ! FIXME omegasi should not be in Sigp% here we should construct the mesh
   ABI_ALLOCATE(Sr%omega_i,(nomega_i))
   !; Sr%omega_i=Sigp%omegasi FIXME this has to be done outside
   ABI_ALLOCATE(Sr%sigcmesi ,(b1gw:b2gw,nkibz,nomega_i,nsppol*nsig_ab))
   ABI_ALLOCATE(Sr%sigxcmesi,(b1gw:b2gw,nkibz,nomega_i,nsppol*nsig_ab))
   Sr%omega_i  =czero
   Sr%sigcmesi =czero
   Sr%sigxcmesi=czero
   if (PRESENT(omega_i)) then 
     ABI_CHECK(SIZE(Sr%omega_i)==SIZE(omega_i),'DIM Sr%omega_i /= omega_i')
     Sr%omega_i=omega_i
   end if
 end if

 Sr%e0        =zero
 Sr%egw       =czero
 Sr%e0gap     =zero
 Sr%sigcme    =czero
 Sr%sigxme    =czero
 Sr%sigxcme   =czero
 Sr%sigcmee0  =czero
 Sr%ze0       =czero
 Sr%dsigmee0  =czero
 Sr%sigmee    =czero
 Sr%omega4sd  =czero
 Sr%sigcme4sd =czero
 Sr%sigxcme4sd=czero
 Sr%degw      =czero

end subroutine allocate_sigma_results
!!***

!----------------------------------------------------------------------

!!****f* m_sigma_results/find_wpoles_for_cd
!! NAME
!!  find_wpoles_for_cd
!!
!! FUNCTION
!!  Find the max frequency needed to account for all the poles
!!  of GW used in the contour deformation technique.
!!
!! INPUTS
!!  Sigp=sigma_parameters
!!
!! OUTPUT
!!  omega_max
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE

subroutine find_wpoles_for_cd(Sigp,Sr,Kmesh,BSt,omega_max)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_wpoles_for_cd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(sigma_parameters),intent(in) :: Sigp
 type(sigma_results),intent(in) :: Sr
 type(bandstructure_type),intent(in) :: Bst
 type(BZ_mesh_type),intent(in) :: Kmesh
 real(dp),intent(out) :: omega_max

!Local variables-------------------------------
!scalars
 integer :: spin,ik_ibz,band_gr,bgw_start,bgw_stop,io,ioe0j
 integer :: ikgw,ikgw_ibz,ikgw_bz,band_gw,nomega_tot
 real(dp) :: e_green,e_screen,theta_mu_minus_e0i,e_qp
 real(dp) :: fact_sp
 !character(len=500) :: msg
!arrays
 real(dp),allocatable :: omegame0i(:)

! *************************************************************************

 omega_max = smallest_real
 !
 ! === Normalization of theta_mu_minus_e0i ===
 ! * If nsppol==2, qp_occ $\in [0,1]$
 fact_sp=one
 if (Bst%nsppol==1) then
   fact_sp=half; if (Bst%nspinor==2) fact_sp=one
 end if
 !
 ! Total number of frequencies for sigma (Spectral function + mesh for the derivative).
 nomega_tot=Sr%nomega_r+Sr%nomega4sd
 ABI_ALLOCATE(omegame0i,(nomega_tot))

 ioe0j=Sr%nomega4sd/2+1
 !
 ! Loop over bands used to construct the Green function.
 do spin=1,Bst%nsppol
   do ik_ibz=1,Bst%nkpt
     do band_gr=1,Bst%nband(ik_ibz+(spin-1)*Bst%nkpt)
       e_green           = Bst%eig(band_gr,ik_ibz,spin)
       theta_mu_minus_e0i= Bst%occ(band_gr,ik_ibz,spin)*fact_sp
       !
       ! Loop over GW states.
       do ikgw=1,Sigp%nkptgw
         bgw_start=Sigp%minbnd(ikgw,spin)
         bgw_stop =Sigp%minbnd(ikgw,spin)
         ikgw_bz  =Sigp%kptgw2bz(ikgw_bz)
         ikgw_ibz =Kmesh%tab(ikgw_bz)

         do band_gw=bgw_start,bgw_stop
           e_qp      = Bst%eig(band_gw,ikgw_ibz,spin)
           !
           ! Get frequencies $\omega$-\epsilon_in$ to evaluate  $d\Sigma/dE$, note the spin
           ! subtract e_KS since we have stored e_KS+ Delta \omega in Sr%omega4sd, not required for AC
           if (Sr%nomega_r>0) omegame0i(1:Sr%nomega_r)=DBLE(Sigp%omega_r(1:Sr%nomega_r))-e_green
           do io=Sr%nomega_r+1,nomega_tot
             !omegame0i(io)=DBLE(Sr%omega4sd(band_gw,ikgw_ibz,io-Sr%nomega_r,spin)) - e_green
             !Sr%omega4sd(jb,ik_ibz,io,spin)=Sr%egw(jb,ik_ibz,spin)+Sigp%deltae*(io-ioe0j)
             omegame0i(io) = e_qp + Sigp%deltae*(io-ioe0j) - e_green
           end do

           do io=1,nomega_tot
             e_screen =  ABS(omegame0i(io))
             if (omegame0i(io)>tol12) then
               !ket(spadc+ig,ios)=ket(spadc+ig,ios)+ct*(one-theta_mu_minus_e0i)
               if ( (one-theta_mu_minus_e0i) > tol12 ) omega_max = MAX(omega_max, e_screen)
             end if
             if (omegame0i(io)<-tol12) then
               !ket(spadc+ig,ios)=ket(spadc+ig,ios)-ct*theta_mu_minus_e0i
               if ( theta_mu_minus_e0i > tol12) omega_max = MAX(omega_max, e_screen)
             end if
           end do

         end do
       end do
       !
     end do
   end do
 end do

 ABI_DEALLOCATE(omegame0i)

end subroutine find_wpoles_for_cd
!!***

!----------------------------------------------------------------------

!!****f* m_sigma_results/etsf_dump_QP
!! NAME
!! etsf_dump_QP
!!
!! FUNCTION
!!  Save the data stored in the sigma_results data type on a NETCDF file.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE

subroutine etsf_dump_QP(Sr,QP_BSt,KS_BSt,Hdr,Cryst,filapp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'etsf_dump_QP'
 use interfaces_14_hidewrite
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Bandstructure_type),intent(in) :: QP_BSt,KS_BSt
 type(Crystal_structure),intent(in) :: Cryst
 type(sigma_results),intent(in) :: Sr
 type(Hdr_type),intent(inout) :: Hdr
 character(len=fnlen),intent(in) :: filapp
!arrays

!Local variables ---------------------------------------
#if defined HAVE_TRIO_ETSF_IO
!scalars
 integer :: ncid,nbgw,ndim_sig,b1gw,b2gw,fform,cplex
 logical :: lstat
 character(len=500) :: msg
 character(len=fnlen) :: filetsf
 type(ETSF_io_low_error) :: Error_data
 type(ETSF_gwdata) :: GWdata
!arrays
 real(dp),target,allocatable :: gw_corrections(:,:,:,:) 
 real(dp),allocatable :: rdata2(:,:),rdata4(:,:,:,:),rdata5(:,:,:,:,:)
#endif

! *************************************************************************

 !@sigma_results

#if defined HAVE_TRIO_ETSF_IO

 filetsf=TRIM(filapp)//'-etsf.nc'
 write(msg,'(3a)')ch10,' etsf_dump_QP : about to open file ',TRIM(filetsf)
 call wrtout(std_out,msg,'COLL')

 call abi_crystal_put(Cryst,filetsf)

 call abi_bands_put  (KS_Bst,filetsf)
  !check here that I have problems

 call etsf_io_low_open_modify(ncid,filetsf,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 ! === Write the abinit header ===
 ! Have to update the occupations?
 fform=502
 call hdr_io_etsf(fform,Hdr,2,ncid)

 ! === Write the GW correction ===
 call wrtout(std_out,' etsf_dump_QP: about to write GW corrections','COLL')

 !FIXME this has to be done in a cleaner way use Sr%egw
 ABI_ALLOCATE(gw_corrections,(2,KS_BSt%mband,KS_BSt%nkpt,KS_BSt%nsppol))
 gw_corrections=zero
 gw_corrections(1,:,:,:) = QP_BSt%eig - KS_BSt%eig
 !$gw_corrections = c2r(Sr%degw)
 GWdata%gw_corrections%data4D => gw_corrections

 call etsf_io_gwdata_put(ncid,GWdata,lstat,Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 nullify(GWdata%gw_corrections%data4D)
 ABI_DEALLOCATE(gw_corrections)

 ! === Up to now we have an ETSF file ===
 ! * Now add variables for internal use in abinit.

 call etsf_io_low_set_define_mode(ncid,lstat,Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 cplex=2; b1gw=Sr%b1gw; b2gw=Sr%b2gw; nbgw=b2gw-b1gw+1
 ndim_sig=Sr%nsppol*Sr%nsig_ab

 call etsf_io_low_write_dim(ncid,'cplex',cplex,lstat,Error_data=Error_data)  
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'b1gw',Sr%b1gw,lstat,Error_data=Error_data) 
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'b2gw',Sr%b2gw,lstat,Error_data=Error_data) 
 ETSF_CHECK_MYERROR(lstat,Error_data)

 ! Number of GW bands.
 call etsf_io_low_write_dim(ncid,'nbgw',nbgw,lstat,Error_data=Error_data) 
 ETSF_CHECK_MYERROR(lstat,Error_data)

 ! Number of k-points calculated.
 call etsf_io_low_write_dim(ncid,'nkptgw',Sr%nkptgw,lstat,Error_data=Error_data) 
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_dim(ncid,'ndim_sig',ndim_sig,lstat,Error_data=Error_data)  
 ETSF_CHECK_MYERROR(lstat,Error_data)
 
 if (Sr%nomega_r>0) then ! No. of real frequencies, might be zero.
   call etsf_io_low_write_dim(ncid,'nomega_r',Sr%nomega_r,lstat,Error_data=Error_data) 
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

 if (Sr%nomega_i>0) then ! No. of imaginary frequencies, might be zero.
   call etsf_io_low_write_dim(ncid,'nomega_i',Sr%nomega_i,lstat,Error_data=Error_data) 
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

 ! No. of points for sigma derivative.
 call etsf_io_low_write_dim(ncid,'nomega4sd',Sr%nomega4sd,lstat,Error_data=Error_data) 
 ETSF_CHECK_MYERROR(lstat,Error_data)

 ! No. of components of sigma (1 if collinear, 4 if noncollinear)
 call etsf_io_low_write_dim(ncid,'nsig_ab',Sr%nsig_ab,lstat,Error_data=Error_data) 
 ETSF_CHECK_MYERROR(lstat,Error_data)

 !call etsf_io_low_write_dim(ncid,'usepawu',Sr%usepawu,lstat,Error_data=Error_data) 
 ! 1 if LDA+U TODO changes name to avoid problems with Hdr
 !ETSF_CHECK_MYERROR(lstat,Error_data)

 ! =======================
 ! == Define variables ===
 ! =======================
 ! TODO use more verbose names!

 call etsf_io_low_def_var(ncid,'gwcalctyp',etsf_io_low_integer,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'omegasrdmax',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'deltae',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'omegasrmax',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'scissor_ene',etsf_io_low_double,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'kptgw',etsf_io_low_double,&
& (/pad('number_of_reduced_dimensions'),pad('nkptgw')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'minbnd',etsf_io_low_integer,(/pad('nkptgw'),pad('number_of_spins')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'maxbnd',etsf_io_low_integer,(/pad('nkptgw'),pad('number_of_spins')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 if (Sr%nomega_r>0) then
   call etsf_io_low_def_var(ncid,'omega_r',etsf_io_low_double,(/'nomega_r'/),lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

!here Sr% starts
 call etsf_io_low_def_var(ncid,'degwgap',etsf_io_low_double,&
& (/pad('number_of_kpoints'),pad('number_of_spins')/),lstat,Error_data=Error_data)
  ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'egwgap',etsf_io_low_double,&
& (/pad('number_of_kpoints'),pad('number_of_spins')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'en_qp_diago',etsf_io_low_double,&
& (/pad('max_number_of_states'),pad('number_of_kpoints'),pad('number_of_spins')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'e0',etsf_io_low_double,&
& (/pad('max_number_of_states'),pad('number_of_kpoints'),pad('number_of_spins')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'e0gap',etsf_io_low_double,&
& (/pad('number_of_kpoints'),pad('number_of_spins')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'sigxme',etsf_io_low_double,&
& (/pad('nbgw'),pad('number_of_kpoints'),pad('ndim_sig')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'vxcme',etsf_io_low_double,&
& (/pad('nbgw'),pad('number_of_kpoints'),pad('ndim_sig')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'vUme',etsf_io_low_double,&
& (/pad('nbgw'),pad('number_of_kpoints'),pad('ndim_sig')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'degw',etsf_io_low_double,&
& (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('number_of_spins')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'dsigmee0',etsf_io_low_double,&
& (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('ndim_sig')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'egw',etsf_io_low_double,&
& (/pad('cplex'),pad('max_number_of_states'),pad('number_of_kpoints'),pad('number_of_spins')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'eigvec_qp',etsf_io_low_double,&
& (/pad('cplex'),pad('max_number_of_states'),pad('max_number_of_states'),pad('number_of_kpoints'),pad('number_of_spins')/),&
& lstat,Error_data=Error_data)
  ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'hhartree',etsf_io_low_double,&
& (/pad('cplex'),pad('nbgw'),pad('nbgw'),pad('number_of_kpoints'),pad('ndim_sig')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 if (Sr%nomega_r>0) then
   call etsf_io_low_def_var(ncid,'sigcme',etsf_io_low_double,&
&   (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('nomega_r'),pad('ndim_sig')/),lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

 call etsf_io_low_def_var(ncid,'sigmee',etsf_io_low_double,&
& (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('ndim_sig')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'sigcmee0',etsf_io_low_double,&
& (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('ndim_sig')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'sigcmesi',etsf_io_low_double,&
& (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('ndim_sig')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'sigcme4sd',etsf_io_low_double,&
& (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('nomega4sd'),pad('ndim_sig')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 if (Sr%nomega_r>0) then
   call etsf_io_low_def_var(ncid,'sigxcme',etsf_io_low_double,&
&    (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('nomega_r'),pad('ndim_sig')/),lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

 if (Sr%nomega_i>0) then
   call etsf_io_low_def_var(ncid,'sigxcmesi',etsf_io_low_double,&
&    (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('nomega_i'),pad('ndim_sig')/),lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

 call etsf_io_low_def_var(ncid,'sigxcme4sd',etsf_io_low_double,&
& (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('nomega4sd'),pad('ndim_sig')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_def_var(ncid,'ze0',etsf_io_low_double,&
&  (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('number_of_spins')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 if (Sr%nomega_i>0) then
   call etsf_io_low_def_var(ncid,'omega_i',etsf_io_low_double,& 
&   (/pad('cplex'),pad('nomega_i')/),lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

 call etsf_io_low_def_var(ncid,'omega4sd',etsf_io_low_double,& 
& (/pad('cplex'),pad('nbgw'),pad('number_of_kpoints'),pad('nomega4sd'),pad('number_of_spins')/),lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 ! =====================
 ! === Start writing ===
 ! =====================
 call etsf_io_low_set_write_mode(ncid,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'gwcalctyp',Sr%gwcalctyp,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'omegasrdmax',Sr%maxomega4sd,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'deltae',Sr%deltae,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'omegasrmax',Sr%maxomega_r,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'scissor_ene',Sr%scissor_ene,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'kptgw',Sr%kptgw,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'minbnd',Sr%minbnd,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'maxbnd',Sr%maxbnd,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'degwgap',Sr%degwgap,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'egwgap',Sr%egwgap,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'en_qp_diago',Sr%en_qp_diago,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'e0',Sr%e0,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'e0gap',Sr%e0gap,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 if (Sr%nomega_r>0) then
   call etsf_io_low_write_var(ncid,'omega_r',Sr%omega_r,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

 call etsf_io_low_write_var(ncid,'sigxme',Sr%sigxme,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'vxcme',Sr%vxcme,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_write_var(ncid,'vUme',Sr%vUme,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 ! * Have to transfer complex arrays
 ABI_ALLOCATE(rdata4,(cplex,b1gw:b2gw,Sr%nkibz,Sr%nsppol))
 rdata4=c2r(Sr%degw)
 call etsf_io_low_write_var(ncid,'degw',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata4)

 ABI_ALLOCATE(rdata4,(cplex,b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 rdata4=c2r(Sr%dsigmee0)
 call etsf_io_low_write_var(ncid,'dsigmee0',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata4)

 ABI_ALLOCATE(rdata4,(cplex,Sr%nbnds,Sr%nkibz,Sr%nsppol))
 rdata4=c2r(Sr%egw)
 call etsf_io_low_write_var(ncid,'egw',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata4)

 ABI_ALLOCATE(rdata5,(cplex,Sr%nbnds,Sr%nbnds,Sr%nkibz,Sr%nsppol))
 rdata5=c2r(Sr%eigvec_qp)
 call etsf_io_low_write_var(ncid,'eigvec_qp',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata5)

 ABI_ALLOCATE(rdata5,(cplex,nbgw,nbgw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 rdata5=c2r(Sr%hhartree)
 call etsf_io_low_write_var(ncid,'hhartree',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata5)

 if (Sr%nomega_r>0) then
   ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))
   rdata5=c2r(Sr%sigcme)
   call etsf_io_low_write_var(ncid,'sigcme',rdata5,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
   ABI_DEALLOCATE(rdata5)
 end if

 ABI_ALLOCATE(rdata4,(cplex,nbgw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 rdata4=c2r(Sr%sigmee)
 call etsf_io_low_write_var(ncid,'sigmee',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata4)

 ABI_ALLOCATE(rdata4,(cplex,nbgw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 rdata4=c2r(Sr%sigcmee0)
 call etsf_io_low_write_var(ncid,'sigcmee0',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata4)

 if (Sr%nomega_i>0) then
  ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
  rdata5=c2r(Sr%sigcmesi)
  call etsf_io_low_write_var(ncid,'sigcmesi',rdata5,lstat,Error_data=Error_data)
  ETSF_CHECK_MYERROR(lstat,Error_data)
  ABI_DEALLOCATE(rdata5)
 end if

 ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))
 rdata5=c2r(Sr%sigcme4sd)
 call etsf_io_low_write_var(ncid,'sigcme4sd',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata5)

 if (Sr%nomega_r>0) then
   ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))
   rdata5=c2r(Sr%sigxcme)
   call etsf_io_low_write_var(ncid,'sigxcme',rdata5,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
   ABI_DEALLOCATE(rdata5)
 end if

 if (Sr%nomega_i>0) then
   ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
   rdata5=c2r(Sr%sigxcmesi)
   call etsf_io_low_write_var(ncid,'sigxcmesi',rdata5,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

 ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))
 rdata5=c2r(Sr%sigxcme4sd)
 call etsf_io_low_write_var(ncid,'sigxcme4sd',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata5)

 ABI_ALLOCATE(rdata4,(cplex,nbgw,Sr%nkibz,Sr%nsppol))
 rdata4=c2r(Sr%ze0)
 call etsf_io_low_write_var(ncid,'ze0',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata4)

 if (Sr%nomega_i>0) then
   ABI_ALLOCATE(rdata2,(cplex,Sr%nomega_i))
   rdata2=c2r(Sr%omega_i)
   call etsf_io_low_write_var(ncid,'omega_i',rdata2,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
   ABI_DEALLOCATE(rdata2)
 end if

 ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol))
 rdata5=c2r(Sr%omega4sd)
 call etsf_io_low_write_var(ncid,'omega4sd',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata5)

 ! === Close the file ===
 call etsf_io_low_close(ncid,lstat,Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 DBG_EXIT("COLL")

#else 
  MSG_ERROR(' ETSF-IO support is not activated. ')
#endif

end subroutine etsf_dump_QP
!!***

!----------------------------------------------------------------------

!!****f* m_sigma_results/abi_etsf_get_QP
!! NAME
!! abi_etsf_get_QP
!!
!! FUNCTION
!!  Initializes several structures used for GW calculations from an external NETCDF 
!!  file written following the ETSF-IO specifications.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwannier
!!
!! CHILDREN
!!      abi_bands_read,allocate_sigma_results,etsf_io_dims_get
!!      etsf_io_low_close,etsf_io_low_open_read,etsf_io_low_read_dim
!!      etsf_io_low_read_var,hdr_io_etsf,init_crystal_from_hdr,wrtout
!!
!! SOURCE

subroutine abi_etsf_get_QP(Sr,KS_BSt,Hdr,Cryst,filapp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_etsf_get_QP'
 use interfaces_14_hidewrite
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Bandstructure_type),intent(out) :: KS_BSt
 type(Crystal_structure),intent(out) :: Cryst
 type(sigma_results),intent(out) :: Sr
 type(Hdr_type),intent(out) :: Hdr
 character(len=fnlen),intent(in) :: filapp
!arrays

!Local variables ---------------------------------------

#if defined HAVE_TRIO_ETSF_IO
!scalars
 integer :: ncid,nbgw,b1gw,b2gw,fform,cplex,timrev
 logical :: lstat
 character(len=500) :: msg
 character(len=fnlen) :: filetsf
 type(ETSF_dims) :: Dims
 type(ETSF_io_low_error) :: Error_data
!arrays
 real(dp),allocatable :: rdata2(:,:),rdata4(:,:,:,:),rdata5(:,:,:,:,:)
#endif

! *************************************************************************

 !@sigma_results

#if defined HAVE_TRIO_ETSF_IO
 filetsf=TRIM(filapp)//'-etsf.nc'
 write(msg,'(3a)')ch10,' abi_etsf_get_QP : about to read file ',TRIM(filetsf)
 call wrtout(std_out,msg,'COLL')

 call etsf_io_low_open_read(ncid,filetsf,lstat,Error_data=Error_data,with_etsf_header=.TRUE.)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 ! === Read KS band structure ===
 call abi_Bands_read(KS_Bst,filapp)

 ! === Read the abinit header ===
 call hdr_io_etsf(fform,Hdr,1,ncid)

 timrev=2
 call init_crystal_from_hdr(Cryst,Hdr,timrev,remove_inv=.FALSE.)

 ! === Read dimensions handled by ETSF ===
 call etsf_io_dims_get(ncid,Dims,lstat,Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 ! FIXME: don't handle k_dependent = 1
 !hdr%bantot   = dims%max_number_of_states * dims%number_of_kpoints * dims%number_of_spins
 !hdr%natom    = dims%number_of_atoms
 Sr%nbnds     = Dims%max_number_of_states
 Sr%nkibz     = Dims%number_of_kpoints
 !hdr%nspden   = dims%number_of_components
 !hdr%nspinor  = dims%number_of_spinor_components
 Sr%nsppol     = Dims%number_of_spins
 !hdr%nsym     = dims%number_of_symmetry_operations
 !hdr%ntypat   = dims%number_of_atom_species

 call etsf_io_low_read_dim(ncid,'b1gw',Sr%b1gw,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_dim(ncid,'b2gw',Sr%b2gw,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 !call etsf_io_low_read_dim(ncid,'nbgw',??,lstat,Error_data=Error_data)
 !ETSF_CHECK_MYERROR(lstat,Error_data)

 !FIXME
 call etsf_io_low_read_dim(ncid,'nkptgw',Sr%nkptgw,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 !call etsf_io_low_read_dim(ncid,'ndim_sig',Sr%ndim_sig,lstat,Error_data=Error_data)
 !ETSF_CHECK_MYERROR(lstat,Error_data)

 ! The following dimensions might be not specified
 call etsf_io_low_read_dim(ncid,'nomega_r',Sr%nomega_r,lstat,Error_data=Error_data)
 if (Sr%nomega_r==etsf_no_dimension) lstat=.TRUE.
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_dim(ncid,'nomega_i',Sr%nomega_i,lstat,Error_data=Error_data)
 if (Sr%nomega_i==etsf_no_dimension) lstat=.TRUE.
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_dim(ncid,'nomega4sd',Sr%nomega4sd,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_dim(ncid,'nsig_ab',Sr%nsig_ab,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 !call etsf_io_low_read_dim(ncid,'usepawu',Sr%usepawu,lstat,Error_data=Error_data)
 !ETSF_CHECK_MYERROR(lstat,Error_data)

 ! == Initialize the structure ===
 call allocate_sigma_results(Sr,&
& Sr%b1gw,Sr%b2gw,Sr%nbnds,Sr%nkibz,Sr%nkptgw,Sr%nsppol,Sr%nsig_ab,Sr%nomega_r,Sr%nomega_i,Sr%nomega4sd)

 b1gw=Sr%b1gw 
 b2gw=Sr%b2gw
 nbgw=b2gw-b1gw+1

 ! ======================
 ! === Read variables ===
 ! ======================

 call etsf_io_low_read_var(ncid,'gwcalctyp',Sr%gwcalctyp,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'omegasrdmax',Sr%maxomega4sd,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'deltae',Sr%deltae,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'omegasrmax',Sr%maxomega_r,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'scissor_ene',Sr%scissor_ene,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'kptgw',Sr%kptgw,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'minbnd',Sr%minbnd,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'maxbnd',Sr%maxbnd,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'degwgap',Sr%degwgap,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'egwgap',Sr%egwgap,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'en_qp_diago',Sr%en_qp_diago,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'e0',Sr%e0,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'e0gap',Sr%e0gap,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 if (Sr%nomega_r>0) then
   call etsf_io_low_read_var(ncid,'omega_r',Sr%omega_r,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

 call etsf_io_low_read_var(ncid,'sigxme',Sr%sigxme,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'vxcme',Sr%vxcme,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 call etsf_io_low_read_var(ncid,'vUme',Sr%vUme,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

 cplex=2
 ABI_ALLOCATE(rdata4,(cplex,b1gw:b2gw,Sr%nkibz,Sr%nsppol))
 call etsf_io_low_read_var(ncid,'degw',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%degw=r2c(rdata4)
 ABI_DEALLOCATE(rdata4)

 ABI_ALLOCATE(rdata4,(cplex,b1gw:b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 call etsf_io_low_read_var(ncid,'dsigmee0',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%dsigmee0=r2c(rdata4) 
 ABI_DEALLOCATE(rdata4)

 ABI_ALLOCATE(rdata4,(cplex,Sr%nbnds,Sr%nkibz,Sr%nsppol))
 call etsf_io_low_read_var(ncid,'egw',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%egw=r2c(rdata4)
 ABI_DEALLOCATE(rdata4)

 ABI_ALLOCATE(rdata5,(cplex,Sr%nbnds,Sr%nbnds,Sr%nkibz,Sr%nsppol))
 call etsf_io_low_read_var(ncid,'eigvec_qp',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%eigvec_qp=r2c(rdata5)
 ABI_DEALLOCATE(rdata5)

 ABI_ALLOCATE(rdata5,(cplex,Sr%b1gw:Sr%b2gw,Sr%b1gw:Sr%b2gw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 call etsf_io_low_read_var(ncid,'hhartree',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%hhartree=r2c(rdata5)
 ABI_DEALLOCATE(rdata5)

 if (Sr%nomega_r>0) then
   ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))
   call etsf_io_low_read_var(ncid,'sigcme',rdata5,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
   Sr%sigcme=r2c(rdata5)
   ABI_DEALLOCATE(rdata5)
 end if

 ABI_ALLOCATE(rdata4,(cplex,nbgw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 call etsf_io_low_read_var(ncid,'sigmee',rdata4,lstat,Error_data=Error_data)
 Sr%sigmee=r2c(rdata4)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 ABI_DEALLOCATE(rdata4)

 ABI_ALLOCATE(rdata4,(cplex,nbgw,Sr%nkibz,Sr%nsppol*Sr%nsig_ab))
 call etsf_io_low_read_var(ncid,'sigcmee0',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%sigcmee0=r2c(rdata4) 
 ABI_DEALLOCATE(rdata4)

 if (Sr%nomega_i>0) then
   ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
   call etsf_io_low_read_var(ncid,'sigcmesi',rdata5,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
   Sr%sigcmesi=r2c(rdata5) 
   ABI_DEALLOCATE(rdata5)
 end if

 ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))
 call etsf_io_low_read_var(ncid,'sigcme4sd',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%sigcme4sd=r2c(rdata5)
 ABI_DEALLOCATE(rdata5)

 if (Sr%nomega_r>0) then
   ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))
   call etsf_io_low_read_var(ncid,'sigxcme',rdata5,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
   Sr%sigxcme=r2c(rdata5)
   ABI_DEALLOCATE(rdata5)
 end if

 if (Sr%nomega_i>0) then
   ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
   call etsf_io_low_read_var(ncid,'sigxcmesi',rdata5,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
   Sr%sigxcmesi=r2c(rdata5)
   ABI_DEALLOCATE(rdata5)
 end if

 ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))
 call etsf_io_low_read_var(ncid,'sigcme4sd',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%sigcme4sd=r2c(rdata5)
 ABI_DEALLOCATE(rdata5)

 if (Sr%nomega_r>0) then
   ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_r,Sr%nsppol*Sr%nsig_ab))
   call etsf_io_low_read_var(ncid,'sigxcme',rdata5,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
   Sr%sigxcme=r2c(rdata5)
   ABI_DEALLOCATE(rdata5)
 end if

 if (Sr%nomega_i>0) then
   ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega_i,Sr%nsppol*Sr%nsig_ab))
   call etsf_io_low_read_var(ncid,'sigxcmesi',rdata5,lstat,Error_data=Error_data)
   Sr%sigxcmesi=r2c(rdata5)
   ETSF_CHECK_MYERROR(lstat,Error_data)
 end if

 ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol*Sr%nsig_ab))
 call etsf_io_low_read_var(ncid,'sigxcme4sd',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%sigxcme4sd=r2c(rdata5)
 ABI_DEALLOCATE(rdata5)

 ABI_ALLOCATE(rdata4,(cplex,nbgw,Sr%nkibz,Sr%nsppol))
 call etsf_io_low_read_var(ncid,'ze0',rdata4,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%ze0=r2c(rdata4)
 !Sr%ze0=CMPLX(rdata4(1,:,:,:),rdata4(2,:,:,:))
 !write(std_out,*)rdata4
 !write(std_out,*)Sr%ze0
 ABI_DEALLOCATE(rdata4)

 if (Sr%nomega_i>0) then
   ABI_ALLOCATE(rdata2,(cplex,Sr%nomega_i))
   call etsf_io_low_read_var(ncid,'omega_i',rdata2,lstat,Error_data=Error_data)
   ETSF_CHECK_MYERROR(lstat,Error_data)
   Sr%omega_i=r2c(rdata2)
   ABI_DEALLOCATE(rdata2)
 end if

 ABI_ALLOCATE(rdata5,(cplex,nbgw,Sr%nkibz,Sr%nomega4sd,Sr%nsppol))
 call etsf_io_low_read_var(ncid,'omega4sd',rdata5,lstat,Error_data=Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)
 Sr%omega4sd=r2c(rdata5)
 ABI_DEALLOCATE(rdata5)

 call etsf_io_low_close(ncid,lstat,Error_data)
 ETSF_CHECK_MYERROR(lstat,Error_data)

#else
 MSG_ERROR('ETSF-IO support is not activated.')
#endif

end subroutine abi_etsf_get_QP

!----------------------------------------------------------------------

END MODULE m_sigma_results
!!***
