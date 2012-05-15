!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_dmft
!! NAME
!!  m_paw_dmft
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_dmft

 use m_profiling

 use defs_basis
 use defs_datatypes

 implicit none

 private 

 public :: init_dmft
 public :: init_sc_dmft
 public :: destroy_dmft
 public :: destroy_sc_dmft
 public :: nullify_sc_dmft
 public :: print_dmft
 public :: print_sc_dmft
!!***

!!****t* m_paw_dmft/paw_dmft_type
!! NAME
!!  paw_dmft_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data for the link
!!  between dmft and paw.
!!  occnd(non-diagonal band occupations for self-consistency), band_in
!!  (say which band are taken into account in the calculation), and the
!   dimensions of these arrays.
!!
!! SOURCE

 type, public :: paw_dmft_type

  integer :: dmft_dc
  ! Type of double counting used in DMFT

  integer :: dmft_iter
  ! Nb of iterations for dmft

  integer :: dmft_solv
  ! choice of solver for DMFT

  integer :: dmftcheck
  ! Check various part of the implementation 

  integer :: dmft_nwlo
  ! dmft frequencies

  integer :: dmft_nwr
  ! dmft frequencies

  integer :: dmft_nwli
  ! dmft frequencies

  integer :: dmftqmc_l
  ! qmc related input

  integer :: dmftqmc_n
  ! qmc related input

  integer :: dmftbandi
  ! Number of bands

  integer :: dmftbandf
  ! Number of bands

  integer :: dmft_rslf
  ! Number of bands
      
  integer :: dmft_prgn
  ! Precise the way of printing the green function.
  !  =1   print green  
  !  =2   print self

  integer :: idmftloop
  ! current iteration in the dmft loop

  integer :: maxlpawu         ! Number of correlated atoms 

      
  integer :: mband
  ! Number of bands
      
  integer :: mbandc
  ! Total number of bands in the Kohn-Sham Basis for PAW+DMFT

  integer :: natom
  ! Number of atom

  integer :: natpawu         ! Number of correlated atoms 

  integer :: nkpt
  ! Number of k-point in the IBZ.

  integer :: nspden

  integer :: nspinor

  integer :: nsppol
      
  integer :: prtdos

  integer :: prtvol

  integer  :: lpsichiortho

  integer  :: use_fixed_self

  real(dp) :: edmft

  real(dp) :: dmft_chpr
  ! Precision on charge required for determination of fermi level (fermi_green) with newton method

  real(dp) :: dmft_fepr
  ! Required precision on Fermi level (fermi_green) during the DMFT SCF cycle, (=> ifermie_cv)
  ! used also for self (new_self)  (=> iself_cv).
      
  real(dp) :: dmft_mxsf
  ! Mixing coefficient for Self-Energy during the SCF DMFT cycle.

  real(dp) :: dmft_lcpr
  ! Required precision on local correlated charge  in order to stop SCF
  ! DMFT cycle (integrate_green) => ichargeloc_cv

  real(dp) :: fermie

  real(dp) :: fermie_lda
      
  real(dp) :: nelectval

  character(len=fnlen) :: filapp
      
  real(dp) :: temp

  integer, pointer :: lpawu(:)

  integer, pointer :: include_bands(:)
  ! for each bands included in the calculation (1..mbandc), include_bands
  ! gives the index in the full band index  (1...mband)
      
  integer, pointer :: exclude_bands(:)
  ! gives the bands than are not in the DMFT calculations.

  real(dp), pointer :: occnd(:,:,:,:)
  ! non diagonal band-occupation for each k-point, polarisation.
      
!  real(dp), pointer :: phi0phiiint(:)
!  ! non diagonal band-occupation for each k-point, polarisation.

  logical, pointer :: band_in(:)
  ! true for each band included in the calculation.

  integer :: use_dmft
  ! 1 if non diagonal occupations are used, else 0
      
  integer :: use_sc_dmft
  ! 1 if calculations have to be carried out self-consistently in the
  ! electronic density.

  complex(dpc), pointer :: psichi(:,:,:,:,:,:)


  real(dp), pointer :: eigen_lda(:,:,:)

  real(dp), pointer :: wtk(:)
  real(dp), pointer :: fixed_self(:,:,:,:)
  real(dp), pointer :: omega_lo(:)
  real(dp), pointer :: omega_li(:)
  real(dp), pointer :: omega_r(:)
  real(dp), pointer :: wgt_wlo(:)

 end type paw_dmft_type
!!***

!----------------------------------------------------------------------

CONTAINS  !========================================================================================
!!***

!!****f* m_paw_dmft/init_sc_dmft
!! NAME
!! init_sc_dmft
!!
!! FUNCTION
!!  Allocate variables used in type paw_dmft_type.
!!
!! INPUTS
!! dmftbandi = lower bound for band states included in DMFT calculation
!! dmftbandf = upper bound for band states included in DMFT calculation
!! mband     = max number of bands 
!! nband     = number of bands for each k-point
!! nkpt      = number of k-points
!! nsppol    = number of spin polarisation 
!! occ       =  occupations
!! use_dmft  = ==1 if dmft is activated
!! use_sc_dmft = for self-consistency in dmft
!!
!! OUTPUTS
!! paw_dmft  = structure of data for dmft
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine init_sc_dmft(dmftbandi,dmftbandf,mband,nband,nkpt,nspden,nspinor,nsppol,occ,use_dmft,paw_dmft,use_sc_dmft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_sc_dmft'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: dmftbandi,dmftbandf,mband,nkpt,nspden,nspinor,nsppol,use_dmft,use_sc_dmft
!type
 type(paw_dmft_type),intent(out) :: paw_dmft
! arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)
!Local variables ------------------------------------
 integer :: iband,icb,ikpt,isppol,nband_k,bdtot_index
 character(len=500) :: message

!************************************************************************

 call nullify_sc_dmft(paw_dmft)
 paw_dmft%mband       = mband
 paw_dmft%dmftbandf   = dmftbandf
 paw_dmft%dmftbandi   = dmftbandi
 paw_dmft%nkpt        = nkpt

!  Spin related variables and check
 paw_dmft%nsppol      = nsppol
 paw_dmft%nspinor     = nspinor
 paw_dmft%nspden      = nspden
 if(nspinor==2.and.nspden==1.and.use_dmft/=0) then
   write(message, '(a,a)' )ch10,&
&   ' nspinor==2 and nspden =1 and usedmft=1 is not implemented'
   call wrtout(std_out,message,'PERS')
   call leave_new('COLL')
 endif
 if(nspinor==1.and.nspden==1.and.use_dmft/=0) then
   write(message, '(a,a)' )ch10,&
&   ' nspinor==1 and nspden =1 and usedmft=1 is not implemented'
   call wrtout(std_out,message,'PERS')
   call leave_new('COLL')
 endif

 paw_dmft%use_dmft    = use_dmft
 paw_dmft%use_sc_dmft = use_sc_dmft
 paw_dmft%mbandc  = 0
 ABI_ALLOCATE(paw_dmft%occnd,(mband,mband,nkpt,nsppol*use_dmft))
 ABI_ALLOCATE(paw_dmft%band_in,(mband*use_dmft))
 ABI_ALLOCATE(paw_dmft%include_bands,((dmftbandf-dmftbandi+1)*use_dmft))
 ABI_ALLOCATE(paw_dmft%exclude_bands,(mband*use_dmft))
! allocate(paw_dmft%ph0phiiint()
 paw_dmft%band_in(:)=.false.
 paw_dmft%occnd=zero
 icb=0
 if(use_dmft==1) then
  do iband=1, mband
   if(iband>=paw_dmft%dmftbandi.and.iband<=paw_dmft%dmftbandf) then
    paw_dmft%band_in(iband)=.true.
    paw_dmft%mbandc = paw_dmft%mbandc+1
    paw_dmft%include_bands(paw_dmft%mbandc) = iband
   else
    icb=icb+1
    paw_dmft%exclude_bands(icb)=iband
   endif
  enddo
  bdtot_index=1
  do isppol=1,nsppol
   do ikpt=1,nkpt
    nband_k=nband(ikpt+(isppol-1)*nkpt)
    do iband=1,nband_k
     paw_dmft%occnd(iband,iband,ikpt,isppol)=occ(bdtot_index)
     bdtot_index=bdtot_index+1
    end do
   end do
  end do
 else
  paw_dmft%mbandc = 0
 endif
 if(paw_dmft%use_dmft > 0 .and. paw_dmft%mbandc /= dmftbandf-dmftbandi+1) then
  write(message, '(a,a,a,a,a,a)' )ch10,&
&  ' WARNING init_dmft',ch10,&
&  '  number of bands in dmft is not correctly computed ',ch10, &
&  '  Action : check the code'
  call wrtout(std_out,message,'PERS')
 endif

end subroutine init_sc_dmft
!!***

!!****f* m_paw_dmft/init_dmft
!! NAME
!! init_dmft 
!!
!! FUNCTION
!!  Allocate variables and setup lda hamiltonian and related data
!!  (init_sc_dmft has to been called before)
!!
!! INPUTS
!!  eigen     = LDA eigenvalues
!!  fermie_lda = LDA Fermi level
!!  psichi    = <chi|Psi> projection of KS states over atomic !wavefunction
!!  nkpt      = number of k-points
!!  nsppol    = number of spin polarisation 
!!  nspinor   = number of spinorial component
!!
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE
!!
!! NOTE 
!! The part of the code which deals
!! with the use of logarithmic frequencies
!! is a modification of the GNU GPL
!! code available on http://dmft.rutgers.edu/ and
!! described in the  RMP paper written by
!! G.Kotliar,  S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti.
!!

subroutine init_dmft(dtset, fermie_lda, fnametmp_app, nspinor, paw_dmft, pawtab, psps, typat)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_dmft'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer  :: nspinor
 real(dp) :: fermie_lda
!type
 type(pseudopotential_type), intent(in) :: psps
 type(dataset_type), intent(in) :: dtset
 type(pawtab_type),intent(in)  :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_dmft_type),intent(out) :: paw_dmft
 character(len=fnlen), intent(in) :: fnametmp_app
!arrays
 integer,intent(in) :: typat(dtset%natom)
!Local variables ------------------------------------
 integer :: iatom,ifreq,ifreq_li,isym,itypat
 real(dp) :: deltaomega,expfac,omegamaxmin,prefacexp,sumwtk
 complex(dpc):: ybcbeg,ybcend
 complex(dpc), allocatable :: tospline_lo(:), splined_li(:),ysplin2_lo(:)
 character(len=500) :: message
 real(dp) :: unit_e,step
! *********************************************************************

 write(message,'(6a)') ch10,' ====================================', &
&                      ch10,' =====  Start of DMFT calculation', &
&                      ch10,' ===================================='
 call wrtout(std_out,message,'COLL')

 unit_e=2_dp
!=======================
!==  Check sym
!=======================
 do isym=1,dtset%nsym
   if(dtset%symafm(isym)<0) then
     write(message, '(a,a)' )ch10,&
&     ' WARNING, symafm negative is not implemented in DMFT '
     call wrtout(std_out,message,'COLL') 
     call leave_new('COLL')
   endif
 enddo

!=======================
!==  Define integers
!=======================
 paw_dmft%maxlpawu=maxval(pawtab(:)%lpawu)
 paw_dmft%fermie_lda=fermie_lda ! in Ha
 paw_dmft%fermie= fermie_lda
 paw_dmft%nelectval= dtset%nelect-float(paw_dmft%dmftbandi-1)*two
 paw_dmft%filapp= fnametmp_app
 paw_dmft%natpawu=dtset%natpawu
 paw_dmft%natom=dtset%natom
 paw_dmft%temp=dtset%tsmear!*unit_e
 paw_dmft%dmft_iter=dtset%dmft_iter
 paw_dmft%dmft_dc=dtset%dmft_dc
 paw_dmft%idmftloop=0
 paw_dmft%prtvol = dtset%prtvol
 paw_dmft%prtdos = dtset%prtdos

!=======================
!==  Fixed self for input
!=======================
 paw_dmft%use_fixed_self=dtset%usedmatpu
 paw_dmft%fixed_self=>dtset%dmatpawu

!=======================
!==  Choose solver
!=======================
 paw_dmft%dmft_solv=dtset%dmft_solv
!  0: LDA, no solver
!  1: LDA+U
! -1: LDA+U but LDA values are not renormalized !
 if(paw_dmft%dmft_solv==0.and.paw_dmft%prtvol>4) then
   write(message, '(a,a)') ch10,' DMFT check: no solver and U=J=0'
 else if(paw_dmft%dmft_solv==1) then
   write(message, '(a,a)') ch10,' DMFT check: static solver'
 else if(paw_dmft%dmft_solv==-1) then
   write(message, '(a,a)') ch10,' DMFT check: static solver without renormalization of projectors: should recover LDA+U'
 else if(paw_dmft%dmft_solv==2) then
   write(message, '(a,a)') ch10,' DMFT uses the Hubbard one solver'
 endif
 if((paw_dmft%dmft_solv==0.and.paw_dmft%prtvol>4).or.&
&   (paw_dmft%dmft_solv>=-1.and.paw_dmft%dmft_solv<=2)) then
   call wrtout(std_out,message,'COLL') 
   call wrtout(ab_out,message,'COLL') 
 endif

 if(paw_dmft%dmft_solv==0) then
   do itypat=1,psps%ntypat
     if(pawtab(itypat)%lpawu/=-1) then
       if((pawtab(itypat)%upawu)>tol5.or.(pawtab(itypat)%jpawu)>tol5) then
          write(message, '(2a,i5,2a,2e15.6)' )ch10,&
&          ' WARNING, option dmft_solv=0 requires upaw=jpaw=0 for species',itypat,ch10,&
&          ' Value of upawu and jpawu are here',pawtab(itypat)%upawu,pawtab(itypat)%jpawu
          call wrtout(std_out,message,'COLL')
          call leave_new('COLL')
        endif
     endif
   enddo
 endif
! todo_ab: why upaw and jpawu are not zero (on bigmac) if lpawu==-1 ?
! if(paw_dmft%dmft_solv==0.and.&
!& (maxval(abs(pawtab(:)%upawu))>tol5.or.maxval(abs(pawtab(:)%jpawu))>tol5)) then
!   write(message, '(a,a,2f12.3)' )ch10,&
!&   ' WARNING, option dmft_solv=0 requires upaw=jpaw=0',maxval(abs(pawtab(:)%upawu)),maxval(abs(pawtab(:)%jpawu))
!   call wrtout(std_out,message,'COLL') 
!   call leave_new('COLL')
! endif

 paw_dmft%dmftcheck=dtset%dmftcheck

 if(paw_dmft%dmftcheck==-1) then
   write(message, '(a,a)' )ch10,&
&   ' BUG: init_dmft,   dmftcheck=-1 should not happend here'
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 endif
 paw_dmft%dmft_nwlo=dtset%dmft_nwlo
 paw_dmft%dmft_nwli=dtset%dmft_nwli
 paw_dmft%dmft_nwr=64
 paw_dmft%dmft_rslf=dtset%dmft_rslf
 paw_dmft%dmft_mxsf=dtset%dmft_mxsf
 paw_dmft%dmftqmc_l=128
 paw_dmft%dmftqmc_n=100000
! todo_ba: put input variables in the input of abinit

!=======================
!==  Variables for DMFT itself
!=======================
 ABI_ALLOCATE(paw_dmft%eigen_lda,(paw_dmft%nsppol,paw_dmft%nkpt,paw_dmft%mbandc))
 paw_dmft%eigen_lda=zero

! allocate(paw_dmft%wtk(paw_dmft%nkpt))
 paw_dmft%wtk=>dtset%wtk
 sumwtk=sum(paw_dmft%wtk)
 if(abs(sumwtk-1_dp)>tol13) then
   write(message, '(a,a,f15.10)' )ch10,&
&   ' WARNING sum of k-point is incorrect',sumwtk
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 endif

 ABI_ALLOCATE(paw_dmft%psichi,(dtset%nsppol,dtset%nkpt,paw_dmft%mbandc,nspinor,dtset%natom,(2*paw_dmft%maxlpawu+1)))
 paw_dmft%psichi=cmplx(zero,zero)
 paw_dmft%lpsichiortho=0
! todo_ba: put psichi in m_local (decide..)


 ABI_ALLOCATE(paw_dmft%lpawu,(paw_dmft%natom))
 do iatom=1,paw_dmft%natom
   paw_dmft%lpawu(iatom)=pawtab(typat(iatom))%lpawu
 enddo

!==  Variables for DMFT related to frequencies
! the part of the code which deals
! with the use of logarithmic frequencies
! is a modification of the GNU GPL
! code available on http://dmft.rutgers.edu/ and
! described in the  RMP paper written by
! G.Kotliar, S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti

!=======================
! Imaginary frequencies
!=======================
 ABI_ALLOCATE(paw_dmft%omega_lo,(paw_dmft%dmft_nwlo))
 ABI_ALLOCATE(paw_dmft%omega_li,(paw_dmft%dmft_nwli))
 ABI_ALLOCATE(paw_dmft%omega_r,(2*paw_dmft%dmft_nwr))
 ABI_ALLOCATE(paw_dmft%wgt_wlo,(paw_dmft%dmft_nwlo))
 step=0.0015_dp
 paw_dmft%omega_r(2*paw_dmft%dmft_nwr)=pi*step*(two*float(paw_dmft%dmft_nwr-1)+one)
 do ifreq=1,2*paw_dmft%dmft_nwr-1
  paw_dmft%omega_r(ifreq)=pi*step*(two*float(ifreq-1)+one)-paw_dmft%omega_r(2*paw_dmft%dmft_nwr)
!  write(std_out,*) ifreq,paw_dmft%omega_r(ifreq)
 enddo

 do ifreq=1,paw_dmft%dmft_nwli
  paw_dmft%omega_li(ifreq)=pi*paw_dmft%temp*(two*float(ifreq-1)+one)
 enddo
 omegamaxmin=paw_dmft%omega_li(paw_dmft%dmft_nwli)-paw_dmft%omega_li(1)
 deltaomega=0.5_dp
 expfac=log(omegamaxmin/deltaomega)/(float(paw_dmft%dmft_nwlo-1)/two)
 prefacexp=omegamaxmin/(exp(expfac*float(paw_dmft%dmft_nwlo-1))-one)
! write(std_out,*) "temp",paw_dmft%temp,paw_dmft%omega_li(1)
! write(std_out,*) "omegamaxmin,deltaomega",omegamaxmin,deltaomega
! write(std_out,*) "expfac,prefacexp",expfac,prefacexp
! write(69,*) paw_dmft%dmft_nwlo
 do ifreq=1,paw_dmft%dmft_nwlo
  paw_dmft%omega_lo(ifreq)=prefacexp*(exp(expfac*float(ifreq-1))-one)+paw_dmft%omega_li(1)
!  write(69,*) paw_dmft%omega_lo(ifreq),0.5
 enddo
! call flush(69)
! call flush(68)
 paw_dmft%omega_lo(1)=paw_dmft%omega_li(1)
 paw_dmft%omega_lo(paw_dmft%dmft_nwlo)=paw_dmft%omega_li(paw_dmft%dmft_nwli)

!=======================
!== construct weight for log. freq.
!=======================
 ABI_ALLOCATE(tospline_lo,(paw_dmft%dmft_nwlo))
 ABI_ALLOCATE(splined_li,(paw_dmft%dmft_nwli))
 ABI_ALLOCATE(ysplin2_lo,(paw_dmft%dmft_nwlo))
 do ifreq=1,paw_dmft%dmft_nwlo 
  tospline_lo=cmplx(0_dp,0_dp)
!  do ifreq1=1,paw_dmft%dmft_nwlo
  tospline_lo(ifreq)=cmplx(1_dp,0_dp)
!  tospline_lo(ifreq1)=ifreq1**2-ifreq1
!  enddo
  splined_li=cmplx(0_dp,0_dp)
!  ybcbeg=cmplx(one/tol16**2,zero)
!  ybcend=cmplx(one/tol16**2,zero)
  ybcbeg=czero
  ybcend=czero


!==         spline delta function
  call spline_complex( paw_dmft%omega_lo, tospline_lo, paw_dmft%dmft_nwlo, &
 & ybcbeg, ybcend, ysplin2_lo)
! do ifreq1=1,paw_dmft%dmft_nwlo 
!  write(6588,*) paw_dmft%omega_lo(ifreq1),ysplin2_lo(ifreq1)
! enddo

  call splint_complex( paw_dmft%dmft_nwlo, paw_dmft%omega_lo, tospline_lo,&
 & ysplin2_lo, paw_dmft%dmft_nwli, paw_dmft%omega_li, splined_li)

!==         accumulate weights
   paw_dmft%wgt_wlo(ifreq)=zero
   do ifreq_li=1,paw_dmft%dmft_nwli
   paw_dmft%wgt_wlo(ifreq)=paw_dmft%wgt_wlo(ifreq)+splined_li(ifreq_li)
   enddo
! do ifreq1=1,paw_dmft%dmft_nwlo 
!  write(6688,*) paw_dmft%omega_lo(ifreq1),tospline_lo(ifreq1)
! enddo
! do ifreq1=1,paw_dmft%dmft_nwli 
!  write(6788,*) paw_dmft%omega_li(ifreq1),splined_li(ifreq1)
 enddo
 ABI_DEALLOCATE(tospline_lo)
 ABI_DEALLOCATE(splined_li)
 ABI_DEALLOCATE(ysplin2_lo)
 if(abs(dtset%pawprtvol)>=3) then
   write(message, '(16x,2(2x,a))') "  Log. Freq     weight   "
   call wrtout(std_out,message,'COLL')
   do ifreq=1,paw_dmft%dmft_nwlo
     write(message, '(3x,a,i4,2(2x,e13.5))') "--ifreq--",ifreq,paw_dmft%omega_lo(ifreq),paw_dmft%wgt_wlo(ifreq)
     call wrtout(std_out,message,'COLL')
   enddo
 endif


! todo_ab finish this
 

!************************************************************************


end subroutine init_dmft
!!***

!!****f* m_paw_dmft/destroy_dmft
!! NAME
!! destroy_dmft
!!
!! FUNCTION
!!  deallocate some variables related to paw_dmft
!!
!! INPUTS
!!  paw_dmft
!!
!! OUTPUT
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_dmft(paw_dmft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_dmft'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type),intent(inout) :: paw_dmft

!Local variables-------------------------------

! *********************************************************************

   if (associated(paw_dmft%psichi))  then
     ABI_DEALLOCATE(paw_dmft%psichi)
   end if
!   paw_dmft%wtk is only an explicit pointer =>dtset%wtk 
!   if (associated(paw_dmft%wtk)) deallocate(paw_dmft%wtk)
   nullify(paw_dmft%wtk)
   nullify(paw_dmft%fixed_self)
   if (associated(paw_dmft%eigen_lda))  then
     ABI_DEALLOCATE(paw_dmft%eigen_lda)
   end if
   if (associated(paw_dmft%omega_lo))  then
     ABI_DEALLOCATE(paw_dmft%omega_lo)
   end if
   if (associated(paw_dmft%omega_li))  then
     ABI_DEALLOCATE(paw_dmft%omega_li)
   end if
   if (associated(paw_dmft%omega_r))  then
     ABI_DEALLOCATE(paw_dmft%omega_r)
   end if
   if (associated(paw_dmft%wgt_wlo))  then
     ABI_DEALLOCATE(paw_dmft%wgt_wlo)
   end if
   if (associated(paw_dmft%lpawu))  then
     ABI_DEALLOCATE(paw_dmft%lpawu)
   end if

end subroutine destroy_dmft
!!***

!!****f* m_paw_dmft/nullify_sc_dmft
!! NAME
!! nullify_sc_dmft
!!
!! FUNCTION
!!  nullify paw_dmft
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_dmft
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine nullify_sc_dmft(paw_dmft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_sc_dmft'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type),intent(inout) :: paw_dmft

!Local variables-------------------------------

!*********************************************************************

 nullify(paw_dmft%occnd)
 nullify(paw_dmft%band_in)
 nullify(paw_dmft%include_bands)
 nullify(paw_dmft%exclude_bands)


end subroutine nullify_sc_dmft
!!***

!!****f* m_paw_dmft/destroy_sc_dmft
!! NAME
!! destroy_sc_dmft
!!
!! FUNCTION
!!  deallocate paw_dmft
!!
!! INPUTS
!!  paw_dmft
!!
!! OUTPUT
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_sc_dmft(paw_dmft)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_sc_dmft'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type),intent(inout) :: paw_dmft

!Local variables-------------------------------
 character(len=500) :: message

! *********************************************************************

 if (( .not. associated(paw_dmft%occnd) .or. .not. associated(paw_dmft%band_in) &
&  .or. .not. associated(paw_dmft%include_bands) .or. .not. associated(paw_dmft%exclude_bands)) &
&  .and. paw_dmft%use_dmft == 1 )  then
  write(message, '(a,a,a,a,a,a)' )ch10,&
&  ' WARNING destroy_sc_dmft',ch10,&
&  '  an array is not associated and is not deallocated with use_dmft==1 ',ch10, &
&  '  Action : check the code'
  call wrtout(std_out,message,'PERS')
 endif
 if ( associated(paw_dmft%occnd) )          then
   ABI_DEALLOCATE(paw_dmft%occnd)
 end if
 if ( associated(paw_dmft%band_in) )        then
   ABI_DEALLOCATE(paw_dmft%band_in)
 end if
 if ( associated(paw_dmft%include_bands) )  then
   ABI_DEALLOCATE(paw_dmft%include_bands)
 end if
 if ( associated(paw_dmft%exclude_bands) )  then
   ABI_DEALLOCATE(paw_dmft%exclude_bands)
 end if


end subroutine destroy_sc_dmft
!!***

!!****f* m_paw_dmft/print_dmft
!! NAME
!! print_dmft
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_dmft(paw_dmft,pawprtvol)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_dmft'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(paw_dmft_type),intent(in) :: paw_dmft
 integer :: pawprtvol

!Local variables-------------------------------
 integer :: ikpt,iband,ifreq,isppol
 character(len=500) :: message
! *********************************************************************

 if( abs(pawprtvol) >= 3 )  then
  write(message,'(4a,3(a,2x,e12.5,a))') &
&   "  -----------------------------------------------",ch10,&
&   "  --- Data for DMFT ",ch10,&
&   "  --- paw_dmft%fermie     = ",paw_dmft%fermie    ,ch10,&
&   "  --- paw_dmft%fermie_lda = ",paw_dmft%fermie_lda,ch10,&
&   "  --- paw_dmft%temp       = ",paw_dmft%temp      ,ch10
  call wrtout(std_out,message,'COLL')
  write(message,'(8(a,2x,i8,a),a)') &
&   "  --- paw_dmft%natpawu    = ",paw_dmft%natpawu   ,ch10,&
&   "  --- paw_dmft%dmft_iter  = ",paw_dmft%dmft_iter ,ch10,&
&   "  --- paw_dmft%dmft_solv  = ",paw_dmft%dmft_solv ,ch10,&
&   "  --- paw_dmft%dmft_nwlo  = ",paw_dmft%dmft_nwlo ,ch10,&
&   "  --- paw_dmft%dmft_nwli  = ",paw_dmft%dmft_nwli ,ch10,&
&   "  --- paw_dmft%dmft_dc    = ",paw_dmft%dmft_dc   ,ch10,&
&   "  --- paw_dmft%dmftqmc_l  = ",paw_dmft%dmftqmc_l ,ch10,&
&   "  --- paw_dmft%dmftqmc_n  = ",paw_dmft%dmftqmc_n ,ch10,&
&   "  -----------------------------------------------"
  call wrtout(std_out,message,'COLL')

!  write(message,'(4a,3(a,2x,f8.3,a),8(a,2x,i8,a),a)') "-----------------------------------------------",ch10,&
!&   "--- Data for DMFT ",ch10,&
!&   "--- paw_dmft%fermie     = ",paw_dmft%fermie    ,ch10,&
!&   "--- paw_dmft%fermie_lda = ",paw_dmft%fermie_lda,ch10,&
!&   "--- paw_dmft%temp       = ",paw_dmft%temp      ,ch10,&
!&   "--- paw_dmft%natpawu    = ",paw_dmft%natpawu   ,ch10,&
!&   "--- paw_dmft%dmft_iter  = ",paw_dmft%dmft_iter ,ch10,&
!&   "--- paw_dmft%dmft_solv  = ",paw_dmft%dmft_solv ,ch10,&
!&   "--- paw_dmft%dmft_nwlo  = ",paw_dmft%dmft_nwlo ,ch10,&
!&   "--- paw_dmft%dmft_nwli  = ",paw_dmft%dmft_nwli ,ch10,&
!&   "--- paw_dmft%dmft_dc    = ",paw_dmft%dmft_dc   ,ch10,&
!&   "--- paw_dmft%dmftqmc_l  = ",paw_dmft%dmftqmc_l ,ch10,&
!&   "--- paw_dmft%dmftqmc_n  = ",paw_dmft%dmftqmc_n ,ch10,&
!&   "-----------------------------------------------"
  if(abs(pawprtvol)>10) then
   call wrtout(std_out,message,'COLL')
   write(message, '(a)') " LDA Eigenvalues "
   do isppol=1,paw_dmft%nsppol
    write(message, '(a,i4)') "--isppol--",isppol
    call wrtout(std_out,message,'COLL')
    do ikpt=1,paw_dmft%nkpt
     write(message, '(a,i4,2x,f14.5,a)') "  -k-pt--",ikpt,paw_dmft%wtk(ikpt),"(<-weight(k-pt))"
    
     call wrtout(std_out,message,'COLL')
     do iband=1,paw_dmft%mbandc
      write(message, '(a,i4,f10.5)') "   -iband--",iband,paw_dmft%eigen_lda(isppol,ikpt,iband)
      call wrtout(std_out,message,'COLL')
     enddo
    enddo
   enddo
   write(message, '(3x,a)') "Log. Freq"
   call wrtout(std_out,message,'COLL')
   do ifreq=1,paw_dmft%dmft_nwlo
    write(message, '(3x,a,i4,2(2x,e13.5))') "--ifreq--",ifreq,paw_dmft%omega_lo(ifreq),paw_dmft%wgt_wlo(ifreq)
    call wrtout(std_out,message,'COLL')
   enddo
  endif
 endif

end subroutine print_dmft
!!***

!!****f* m_paw_dmft/print_sc_dmft
!! NAME
!! print_sc_dmft
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_sc_dmft(paw_dmft,pawprtvol)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_sc_dmft'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!type
 type(paw_dmft_type),intent(in) :: paw_dmft
 integer :: pawprtvol

!Local variables-------------------------------
 integer :: iband
 character(len=500) :: message
! *********************************************************************

 if( abs(pawprtvol) >= 3 )  then
   write(message,'(5a,8(a,2x,i5,a),a)')ch10,"-----------------------------------------------",ch10,&
&    "--- Data for SC DMFT ",ch10,&
&    "--- paw_dmft%mband       = ",paw_dmft%mband,ch10,&
&    "--- paw_dmft%dmftbandf   = ",paw_dmft%dmftbandf,ch10,&
&    "--- paw_dmft%dmftbandi   = ",paw_dmft%dmftbandi,ch10,&
&    "--- paw_dmft%nkpt        = ",paw_dmft%nkpt,ch10,&
&    "--- paw_dmft%nsppol      = ",paw_dmft%nsppol,ch10,&
&    "--- paw_dmft%use_dmft    = ",paw_dmft%use_dmft,ch10,&
&    "--- paw_dmft%use_sc_dmft = ",paw_dmft%use_sc_dmft,ch10,&
&    "--- paw_dmft%mbandc      = ",paw_dmft%mbandc,ch10,&
&    "-----------------------------------------------"
   call wrtout(std_out,message,'COLL')
   write(message, '(a)') " paw_dmft%band_in"
   call wrtout(std_out,message,'COLL')
   write(message, '(100i5)') (iband,iband=1,min(paw_dmft%mband,100))
   call wrtout(std_out,message,'COLL')
   write(message, '(100L3)') (paw_dmft%band_in(iband),iband=1,min(paw_dmft%mband,100))
   call wrtout(std_out,message,'COLL')
   do iband=1,paw_dmft%mbandc
     write(message,*) "include_bands",iband,paw_dmft%include_bands(iband)
     call wrtout(std_out,message,'COLL')
   enddo
   write(message, '(a,a,i4,a)' )ch10,&
&    'The',paw_dmft%mband-paw_dmft%dmftbandf+paw_dmft%dmftbandi-1,&
&    '  Following bands are excluded from the DMFT calculation  '
   call wrtout(std_out,message,'COLL')
   write(message,'(100i5)') (paw_dmft%exclude_bands(iband),iband=1,min(paw_dmft%mband-paw_dmft%dmftbandf+paw_dmft%dmftbandi-1,100))
   call wrtout(std_out,message,'COLL')
 endif

end subroutine print_sc_dmft

END MODULE m_paw_dmft
!!***
