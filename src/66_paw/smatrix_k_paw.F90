!{\src2tex{textfont=tt}}
!!****f* ABINIT/smatrix_k_paw
!! NAME
!! smatrix_k_paw
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cprj_k (cprj_type) :: cprj for occupied bands at point k
!!  cprj_kb :: cprj for occupied bands at point k+b
!!  dtefield :: structure referring to all efield and berry's phase variables
!!  kdir :: integer giving direction along which overlap is computed for ket
!!  kfor :: integer indicating whether to compute forward (1) or backward (2)
!!    along kpt string
!!  natom :: number of atoms in cell
!!  typat :: typat(natom) type of each atom
!!  (optional)bdir :: integer giving direction along which overlap is computed for bra
!!  (optional)bfor :: integer indicating whether to compute forward (1) or backward (2)
!!    along kpt string
!!
!! OUTPUT
!! smat_k_paw :: array of the on-site PAW parts of the overlaps between Bloch states at points
!!   k and k+b, for the various pairs of bands, that is, the on-site part of 
!!   <u_nk|u_mk+b>
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This routine assumes that the cprj are not explicitly ordered by 
!! atom type.
!!
!! PARENTS
!!      berry_linemin,berryphase_new,cgwf,make_grad_berry,update_mmat
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine smatrix_k_paw(cprj_k,cprj_kb,dtefield,kdir,kfor,natom,smat_k_paw,typat,bdir,bfor)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'smatrix_k_paw'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: kdir,kfor,natom
 type(efield_type),intent(in) :: dtefield
 type(cprj_type),intent(in) :: cprj_k(natom,dtefield%nband_occ)
 type(cprj_type),intent(in) :: cprj_kb(natom,dtefield%nband_occ)

!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(out) :: smat_k_paw(2,dtefield%nband_occ,dtefield%nband_occ)

!optional arguments
 integer,optional,intent(in) :: bdir,bfor

!Local variables---------------------------
!scalars
 integer :: bdx,bsig,iatom,iband,ilmn,itypat
 integer :: jband,jlmn,kdx,klmn,ksig,twdx
 logical :: need_conjg
 complex(dpc) :: cpk,cpkb,cterm,paw_onsite

! *************************************************************************

!initialize smat_k_paw
 smat_k_paw(:,:,:) = zero

 if(present(bdir)) then
   bsig = -2*bfor+3; ksig=-2*kfor+3
   bdx = 2*bdir-(bsig+1)/2; kdx = 2*kdir-(ksig+1)/2
   twdx = dtefield%twind(bdx,kdx)
   if (twdx < 0) then
     twdx = -twdx; need_conjg=.true.
   else
     need_conjg=.false.
   end if
 end if

 do iatom = 1, natom
   itypat = typat(iatom)

   do ilmn=1,dtefield%lmn_size(itypat)
     do jlmn=1,dtefield%lmn_size(itypat)
       if (ilmn >= jlmn) then
         klmn = ilmn*(ilmn-1)/2+jlmn
       else
         klmn = jlmn*(jlmn-1)/2+ilmn
       end if 
       if (present(bdir)) then
         paw_onsite = cmplx(dtefield%qijb_bk(1,klmn,iatom,twdx),&
&         dtefield%qijb_bk(2,klmn,iatom,twdx))
         if(need_conjg) paw_onsite = conjg(paw_onsite)
       else
         paw_onsite = cmplx(dtefield%qijb_kk(1,klmn,iatom,kdir),&
&         dtefield%qijb_kk(2,klmn,iatom,kdir))
         if (kfor > 1) paw_onsite = conjg(paw_onsite)
       end if
       do iband = 1, dtefield%nband_occ
         cpk=cmplx(cprj_k(iatom,iband)%cp(1,ilmn),cprj_k(iatom,iband)%cp(2,ilmn))
         do jband = 1, dtefield%nband_occ
           cpkb=cmplx(cprj_kb(iatom,jband)%cp(1,jlmn),cprj_kb(iatom,jband)%cp(2,jlmn))
           cterm = conjg(cpk)*paw_onsite*cpkb
           smat_k_paw(1,iband,jband) = smat_k_paw(1,iband,jband)+real(cterm)
           smat_k_paw(2,iband,jband) = smat_k_paw(2,iband,jband)+aimag(cterm)
         end do ! end loop over jband
       end do ! end loop over iband
     end do ! end loop over ilmn
   end do ! end loop over jlmn

 end do ! end loop over atoms

 end subroutine    smatrix_k_paw
!!***

