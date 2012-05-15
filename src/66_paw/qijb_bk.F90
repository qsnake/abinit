!{\src2tex{textfont=tt}}
!!****f* ABINIT/qijb_bk
!! NAME
!! qijb_bk
!!
!! FUNCTION
!! Routine which computes PAW onsite part of wavefunction overlap for Bloch
!! functions at two k-points k+b1 and k+b2,  and also exp(-i (b2-b1).R) at each site. These
!! quantities are stored in an efield structure and used in PAW Berrys Phase 
!! calculations of magnetization.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gprimd(3,3) :: dimensioned primitive translations of reciprocal lattice
!!  natom :: number of atoms in unit cell
!!  nattyp :: number of atoms of each type
!!  ntypat :: number of types of atoms in unit cell
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtefield :: efield structure containing electric field and Berrys Phase data
!!
!! NOTES
!!
!! PARENTS
!!      initberry
!!
!! CHILDREN
!!      initylmr,sbf8,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine qijb_bk(dtefield,gprimd,natom,ntypat,pawang,pawrad,pawtab,typat)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_efield

 use m_radmesh,  only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'qijb_bk'
 use interfaces_32_util
 use interfaces_65_psp
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 type(efield_type),intent(inout) :: dtefield
 type(pawang_type),intent(in) :: pawang

!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: gprimd(3,3)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: bdir,bk_ind,bsig,iatom,ir,isel,itypat,kdir
 integer :: klm,kln,klmn,lbess,lbesslm,lmin,lmax,mbess,mesh_size,mu
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 real(dp) :: arg,bessg,bnorm,intg,rterm
 complex(dpc) :: cterm,etb,ifac
!arrays
 real(dp) :: bb(3),bbn(3),bcart(3),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: ff(:),j_bessel(:,:),ylmb(:)
! the following is (-i)^L mod 4.
 complex(dpc),dimension(0:3) :: il(0:3)=(/(1.0,0.0),(0.0,-1.0),(-1.0,0.0),(0.0,1.0)/)
! *************************************************************************
 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr

 dtefield%qijb_bk(:,:,:,:) = zero

 do iatom = 1, natom

   itypat = typat(iatom)

   mesh_size = pawrad(itypat)%mesh_size

   ABI_ALLOCATE(j_bessel,(mesh_size,pawang%l_size_max))
   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(ylmb,(pawang%l_size_max*pawang%l_size_max))

   do bdir = 1, 3
     do bsig = -1, 1, 2
       kdir = mod(bdir,3)+1
       
       bb(:) = bsig*dtefield%dkvecs(:,bdir) - dtefield%dkvecs(:,kdir) ! will need exp(i(b-k).r), here is b-k
!      get cartesian positions of atom in cell
       do mu=1,3
         bcart(mu)=dot_product(bb(:),gprimd(mu,:))
       end do
!      here is exp(i b.R) for current atom
       bk_ind = 2*bdir-1+(bsig+1)/2
       etb = cmplx(dtefield%expibi(1,iatom,bk_ind),&
&       dtefield%expibi(2,iatom,bk_ind))
       
!      bbn is b-hat (the unit vector in the b direction) 
       bnorm=dsqrt(dot_product(bcart,bcart))
       bbn(:) = bcart(:)/bnorm
       
!      as an argument to the bessel function, need 2pi*b*r = 1 so b is re-normed to two_pi
       bnorm = two_pi*bnorm
       do ir=1,mesh_size
         arg=bnorm*pawrad(itypat)%rad(ir)
         call sbf8(pawang%l_size_max,arg,j_bessel(ir,:)) ! spherical bessel functions at each mesh point
       end do ! end loop over mesh
!      compute Y_LM(b) here
       call initylmr(pawang%l_size_max,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,bbn,ylmb(:),ylmgr)
       
       do klmn = 1, pawtab(itypat)%lmn2_size
         klm =pawtab(itypat)%indklmn(1,klmn)
         kln =pawtab(itypat)%indklmn(2,klmn)
         lmin=pawtab(itypat)%indklmn(3,klmn)
         lmax=pawtab(itypat)%indklmn(4,klmn)
         do lbess = lmin, lmax 
           ifac = il(mod(lbess,4))
           do mbess = -lbess, lbess
             lbesslm = lbess*lbess+lbess+mbess+1
             isel=pawang%gntselect(lbesslm,klm)
             if (isel > 0) then
               bessg = pawang%realgnt(isel)
               ff(1:mesh_size)=(pawtab(itypat)%phiphj(1:mesh_size,kln)&
&               -pawtab(itypat)%tphitphj(1:mesh_size,kln))&
&               *j_bessel(1:mesh_size,lbess+1)
               call simp_gen(intg,ff,pawrad(itypat))
               rterm = four_pi*bessg*intg*ylmb(lbesslm)
               cterm = etb*ifac*rterm
               dtefield%qijb_bk(1,klmn,iatom,bk_ind) = &
               dtefield%qijb_bk(1,klmn,iatom,bk_ind) + real(cterm)
               dtefield%qijb_bk(2,klmn,iatom,bk_ind) = &
               dtefield%qijb_bk(2,klmn,iatom,bk_ind) + aimag(cterm)
               
             end if ! end selection on non-zero Gaunt factors
           end do ! end loop on mbess = -lbess, lbess
         end do ! end loop on lmin-lmax bessel l values
       end do ! end loop on lmn2_size klmn basis pairs
     end do ! end loop over bsig
   end do ! end loop over bdir

   ABI_DEALLOCATE(j_bessel)
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(ylmb)

 end do ! end loop on natom


 end subroutine qijb_bk
!!***
