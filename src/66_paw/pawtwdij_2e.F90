!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawtwdij_2e
!! NAME
!! pawtwdij_2e
!!
!! FUNCTION
!! compute phase-twisted contribution to ${D}^1_{ij}-\tilde{D}^1_{ij}$ due to 
!! Hatree potential of core charge and compensation charge moments.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gprimd(3,3) = primitive translations in recip space
!!  natom = number of atoms in unit cell
!!  ntypat = number of types of atoms in unit cell
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  typat = typat(natom) list of atom types
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtefield <type(efield_type)> = dtefield%twdij0 is updated
!!
!! NOTES
!! This term corresponds to term (2e) of Eq. 49 in Torrent et al.,
!! CMS 42, 337 (2008), including phase shifts as needed for orbital
!! magnetization:
!! $-\sum_{LM}\int_{\Omega_R}e^{i\mathbf{b.r}}v_H[\tilde{n}_{Zc}]\hat{Q}^{LM}_{ij}(\mathbf{r})e^{-i\mathbf{k.r}}$
!! where $\mathbf{b}$ is the bra shift vector and $\mathbf{k}$ is the ket shift vector, both are determined in
!! initberry.F90.
!!
!! PARENTS
!!      initberry
!!
!! CHILDREN
!!      initylmr,jbessel,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! macro to go from row-column indexing to combined indexing
#define RCC(glmn,hlmn) max(glmn,hlmn)*(max(glmn,hlmn)-1)/2+min(glmn,hlmn)

!macro to go from l,m angular momentum indexing to combined indexing
#define LMC(lval,mval) lval*lval+lval+mval+1

#include "abi_common.h"

 subroutine pawtwdij_2e(dtefield,gprimd,natom,ntypat,pawrad,pawtab,psps,typat)

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
#define ABI_FUNC 'pawtwdij_2e'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 type(efield_type),intent(inout) :: dtefield
 type(pseudopotential_type),intent(in) :: psps

!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: gprimd(3,3)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: bdir,bdx,bfor,bl,bm,blm,bln,blmn,bsig,clmn
 integer :: dij_ind,iatom,itypat,ilm2
 integer :: ir,kdir,kfor,kl,km,klm,kln,klmn,ksig,kdx
 integer :: lcmax,ll,mesh_size,mm,mu
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 real(dp) :: bessarg,intgrl,jlx,jldx,jldx2,knorm,xx
 complex(dpc) :: vijsum
!arrays
 real(dp) :: bb(3),bcart(3),kb(3),kcart(3),kbn(3),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: ff(:),ylmk(:)
! the following is (-i)^L mod 4.
 complex(dpc),dimension(0:3) :: iml(0:3)=(/(1.0,0.0),(0.0,-1.0),(-1.0,0.0),(0.0,1.0)/)
! *************************************************************************

 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr


 do itypat = 1, ntypat

   lcmax = pawtab(itypat)%l_size  ! lcmax - 1 is highest angular momentum state used in expansion
   ABI_ALLOCATE(ylmk,(lcmax*lcmax))

   mesh_size = pawrad(itypat)%mesh_size
   ABI_ALLOCATE(ff,(mesh_size))

   do bdir = 1, 3
     do bsig = -1, 1, 2
       bfor = (3-bsig)/2
       bb(:) = bsig*dtefield%dkvecs(:,bdir) ! bra vector

       do mu=1,3
         bcart(mu)=dot_product(bb(:),gprimd(mu,:))
       end do

       do kdir = 1, 3
         do ksig = -1, 1, 2
           if (kdir == bdir) cycle ! never need the kdir // bdir terms
           
           kfor = (3-ksig)/2

!          see efield_type for documentation of kdx
           bdx = 2*bdir-bfor+1; kdx=2*kdir-kfor+1
           dij_ind = dtefield%indhk(bdx,kdx)

           kb(:) = ksig*dtefield%dkvecs(:,kdir)  ! ket vector

           do mu=1,3
             kcart(mu)=dot_product(kb(:),gprimd(mu,:))
           end do

!          form b_k - b_b vector
           kcart(1:3) = kcart(1:3) - bcart(1:3)

           knorm=dsqrt(dot_product(kcart,kcart))
           if (knorm < tol12) then
             kbn(:) = zero
             ylmk(:) = zero; ylmk(1) = one/sqrt(four_pi)
           else
             kbn(:) = kcart(:)/knorm ! unit vector in kb direction
             call initylmr(lcmax,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,kbn,ylmk(:),ylmgr)
           end if

           knorm = two_pi*knorm ! re-normed kb for calls to bessel fnc

           do blmn = 1, pawtab(itypat)%lmn_size
             bl = psps%indlmn(1,blmn,itypat)
             bm = psps%indlmn(2,blmn,itypat)
             blm = psps%indlmn(4,blmn,itypat)
             bln = psps%indlmn(5,blmn,itypat) 

             do klmn = 1, pawtab(itypat)%lmn_size
               kl = psps%indlmn(1,klmn,itypat)
               km = psps%indlmn(2,klmn,itypat)
               klm = psps%indlmn(4,klmn,itypat)
               kln = psps%indlmn(5,klmn,itypat) 
               
               clmn = RCC(blmn,klmn)

               vijsum = cmplx(zero,zero)

               do ll = abs(kl-bl),kl+bl,2
                 do ir = 1, mesh_size
                   xx = pawrad(itypat)%rad(ir)
                   bessarg = xx*knorm
                   call jbessel(jlx,jldx,jldx2,ll,0,bessarg)
                   ff(ir) = xx*xx*jlx*pawtab(itypat)%VHntZC(ir)*pawtab(itypat)%shapefunc(ir,ll+1)
                 end do
                 call simp_gen(intgrl,ff,pawrad(itypat))

                 do mm = -ll, ll
                   ilm2 = LMC(ll,mm)
                   vijsum = vijsum + iml(mod(ll,4))*pawtab(itypat)%qijl(ilm2,clmn)*ylmk(ilm2)*intgrl 
                 end do ! end loop over mm
               end do ! end loop over ll
               
               vijsum = -4.d0*pi*vijsum
               do iatom = 1, natom ! store result
                 if(typat(iatom) == itypat) then
                   dtefield%twdij0(1,blmn,klmn,iatom,dij_ind) = dtefield%twdij0(1,blmn,klmn,iatom,dij_ind) + &
&                   real(vijsum)
                   dtefield%twdij0(2,blmn,klmn,iatom,dij_ind) = dtefield%twdij0(2,blmn,klmn,iatom,dij_ind) + &
&                   aimag(vijsum)
                 end if
               end do

             end do ! end loop over ket states
           end do ! end loop over bra states

         end do ! end loop over ksig
       end do ! end loop over kdir

     end do ! end loop over bsig
   end do ! end loop over bdir

   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(ylmk)

 end do ! end loop on ntypat
 

 end subroutine pawtwdij_2e
!!***
