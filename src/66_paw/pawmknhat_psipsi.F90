!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmknhat_psipsi
!! NAME
!! pawmknhat_psipsi
!!
!! FUNCTION
!! PAW only:
!! Compute on the fine FFT grid the compensation charge density (and derivatives) associated
!! to the product of two wavefunctions n_{12}(r) = \Psi_1* \Psi_2. Based on pawmknhat.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MG, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  Cprj1(natom,nspinor), Cprj2(natom,nspinor) <type(Cprj_type)>=
!!   projected input wave functions <Proj_i|Cnk> with all NL projectors corresponding to
!!   the \Psi_1 and \Psi_2, respectively.
!!  ider= 0: nhat(r) is computed
!!        1: cartesian derivatives of nhat(r) are computed
!!        2: nhat(r) and derivatives are computed
!!        Note: ider>0 not compatible with ipert>0
!!  izero=if 1, unbalanced components of nhat(g) have to be set to zero
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms on current process, size of PAW arrays
!!  typat(natom)=Type of each atom.
!!  nfft=number of point on the rectangular fft grid
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat12_grdim= 0 if grnhat12 array is not used ; 1 otherwise
!!  ntypat=number of types of atoms in unit cell.
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  === if ider=0 or 2
!!    nhat12(2,nfft,nspinor**2)=nhat on fine rectangular grid
!!  === if ider=1 or 2
!!    grnhat12(nfft,nspinor**2,3)=derivatives of nhat on fine rectangular grid (and derivatives)
!!
!! PARENTS
!!      calc_sigx_me
!!
!! CHILDREN
!!      fourdp,klmn2ijlmn,pawgylm,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawmknhat_psipsi(cprj1,cprj2,ider,izero,mpi_enreg,natom,nfft,ngfft,nhat12_grdim,&
&  nspinor,ntypat,typat,paral_kgb,pawang,pawfgrtab,grnhat12,nhat12,pawtab)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_lmn_indices,  only : klmn2ijlmn

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawmknhat_psipsi'
 use interfaces_53_ffts
 use interfaces_66_paw, except_this_one => pawmknhat_psipsi
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ider,izero,natom,nfft,nhat12_grdim,ntypat,paral_kgb,nspinor
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(out) :: grnhat12(2,nfft,nspinor**2,3*nhat12_grdim)
 real(dp),intent(out) :: nhat12(2,nfft,nspinor**2)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)
 type(Cprj_type),intent(in) :: Cprj1(natom,nspinor),Cprj2(natom,nspinor)

!Local variables ---------------------------------------
!scalars
 integer :: iatom,ic,ils,ilslm,isp1,isp2,isploop,itypat,jc,klm,klmn
 integer :: lmax,lmin,lm_size,mm,optgr0,optgr1,cplex
 integer :: ilmn,jlmn,lmn_size,lmn2_size
 real(dp) :: re_p,im_p
 logical :: compute_grad,compute_nhat
!arrays
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 real(dp) :: rdum(1),cpf(2),cpf_ql(2)
 real(dp),allocatable :: work(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 if (ider>0.and.nhat12_grdim==0) then
   MSG_BUG('Gradients of nhat required but not allocated !')
 end if
 if (nspinor==2) then
   MSG_BUG('nspinor==2 not coded!')
 end if
 if (paral_kgb/=0) then
   MSG_BUG('paral_kgb/=0 not coded!')
 end if

!Initialisations
 compute_nhat=(ider==0.or.ider==2)
 compute_grad=(ider==1.or.ider==2)
 if ((.not.compute_nhat).and.(.not.compute_grad)) return

 if (compute_nhat) nhat12=zero
 if (compute_grad) grnhat12=zero

 if (compute_grad) then
   MSG_BUG('compute_grad not tested!')
 end if

!------------------------------------------------------------------------
!----- Loop over atoms
!------------------------------------------------------------------------
 do iatom=1,natom
   itypat = typat(iatom)
   lm_size=pawfgrtab(iatom)%l_size**2
   lmn_size  = pawtab(itypat)%lmn_size
   lmn2_size = pawtab(itypat)%lmn2_size

!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   if (((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)).or.&
&   ((compute_grad).and.(pawfgrtab(iatom)%gylmgr_allocated==0))) then

     optgr0=0; optgr1=0
     if ((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)) then
       if (associated(pawfgrtab(iatom)%gylm))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(pawfgrtab(iatom)%nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
     end if

     if ((compute_grad).and.(pawfgrtab(iatom)%gylmgr_allocated==0)) then
       if (associated(pawfgrtab(iatom)%gylmgr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,pawfgrtab(iatom)%nfgd,pawfgrtab(iatom)%l_size**2))
       pawfgrtab(iatom)%gylmgr_allocated=2;optgr1=1
     end if

     call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,rdum,&
&     lm_size,pawfgrtab(iatom)%nfgd,optgr0,optgr1,0,pawtab(itypat),&
&     pawfgrtab(iatom)%rfgd,pawfgrtab(iatom)%rfgd_allocated)
   end if

   do isploop=1,nspinor**2    ! Loop over density components of the compensation charge.
!    TODO Here we might take advantage of symmetry relations between the four components if nspinor==2
     isp1=spinor_idxs(1,isploop)
     isp2=spinor_idxs(2,isploop)

     do klmn=1,lmn2_size  ! Loop over ij channels of this atom type.
       klm =pawtab(itypat)%indklmn(1,klmn)
       lmin=pawtab(itypat)%indklmn(3,klmn)  ! abs(il-jl)
       lmax=pawtab(itypat)%indklmn(4,klmn)  ! il+jl

       call klmn2ijlmn(klmn,lmn_size,ilmn,jlmn)  ! This mapping should be stored in pawtab_type

!      Retrieve the factor due to the PAW projections.
       re_p =  Cprj1(iatom,isp1)%cp(1,ilmn) * Cprj2(iatom,isp2)%cp(1,jlmn) &
&       +Cprj1(iatom,isp1)%cp(2,ilmn) * Cprj2(iatom,isp2)%cp(2,jlmn) &
&       +Cprj1(iatom,isp1)%cp(1,jlmn) * Cprj2(iatom,isp2)%cp(1,ilmn) &
&       +Cprj1(iatom,isp1)%cp(2,jlmn) * Cprj2(iatom,isp2)%cp(2,ilmn)

       im_p =  Cprj1(iatom,isp1)%cp(1,ilmn) * Cprj2(iatom,isp2)%cp(2,jlmn) &
&       -Cprj1(iatom,isp1)%cp(2,ilmn) * Cprj2(iatom,isp2)%cp(1,jlmn) &
&       +Cprj1(iatom,isp1)%cp(1,jlmn) * Cprj2(iatom,isp2)%cp(2,ilmn) &
&       -Cprj1(iatom,isp1)%cp(2,jlmn) * Cprj2(iatom,isp2)%cp(1,ilmn)

       cpf(1)=re_p*pawtab(itypat)%dltij(klmn)*half
       cpf(2)=im_p*pawtab(itypat)%dltij(klmn)*half

       if (compute_nhat) then
         do ils=lmin,lmax,2   ! Sum over (L,M)
           do mm=-ils,ils
             ilslm=ils*ils+ils+mm+1
             if (pawang%gntselect(ilslm,klm)>0) then
               cpf_ql(1)=cpf(1)*pawtab(itypat)%qijl(ilslm,klmn)
               cpf_ql(2)=cpf(2)*pawtab(itypat)%qijl(ilslm,klmn)
               do ic=1,pawfgrtab(iatom)%nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 nhat12(1,jc,isploop)=nhat12(1,jc,isploop)+cpf_ql(1)*pawfgrtab(iatom)%gylm(ic,ilslm)
                 nhat12(2,jc,isploop)=nhat12(2,jc,isploop)+cpf_ql(2)*pawfgrtab(iatom)%gylm(ic,ilslm)
               end do
             end if
           end do
         end do
       end if ! compute_nhat

       if (compute_grad) then
         do ils=lmin,lmax,2  ! Sum over (L,M)
           do mm=-ils,ils
             ilslm=ils*ils+ils+mm+1
             if (pawang%gntselect(ilslm,klm)>0) then
               cpf_ql(1)=cpf(1)*pawtab(itypat)%qijl(ilslm,klmn)
               do ic=1,pawfgrtab(iatom)%nfgd
                 jc=pawfgrtab(iatom)%ifftsph(ic)
                 grnhat12(1,jc,isploop,1)=grnhat12(1,jc,isploop,1)+cpf_ql(1)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                 grnhat12(1,jc,isploop,2)=grnhat12(1,jc,isploop,2)+cpf_ql(1)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                 grnhat12(1,jc,isploop,3)=grnhat12(1,jc,isploop,3)+cpf_ql(1)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)

                 grnhat12(2,jc,isploop,1)=grnhat12(2,jc,isploop,1)+cpf_ql(2)*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
                 grnhat12(2,jc,isploop,2)=grnhat12(2,jc,isploop,2)+cpf_ql(2)*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
                 grnhat12(2,jc,isploop,3)=grnhat12(2,jc,isploop,3)+cpf_ql(2)*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)
               end do
             end if
           end do
         end do
       end if ! compute_grad

     end do  ! klmn (ij channels)
   end do ! isploop (density components of the compensation charge)

   if (pawfgrtab(iatom)%gylm_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
     ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(0,0))
     pawfgrtab(iatom)%gylm_allocated=0
   end if
   if (pawfgrtab(iatom)%gylmgr_allocated==2) then
     ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
     ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(0,0,0))
     pawfgrtab(iatom)%gylmgr_allocated=0
   end if

 end do ! iatom

!----- Avoid unbalanced g-components numerical errors -----!
 if (izero==1.and.compute_nhat) then
   ABI_ALLOCATE(work,(2,nfft))
   cplex=2
   do isp1=1,MIN(2,nspinor**2)
     call fourdp(cplex,work,nhat12(:,:,isp1),-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     call zerosym(work,cplex,mpi_enreg,ngfft(1),ngfft(2),ngfft(3))
     call fourdp(cplex,work,nhat12(:,:,isp1),+1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   end do
   ABI_DEALLOCATE(work)
 end if

 DBG_EXIT("COLL")

end subroutine pawmknhat_psipsi
!!***
