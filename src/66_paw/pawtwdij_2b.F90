!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawtwdij_2b
!! NAME
!! pawtwdij_2b
!!
!! FUNCTION
!! compute phase-twisted contribution to ${D}^1_{ij}-\tilde{D}^1_{ij}$ due to 
!! Hatree potential of core charge.
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
!! This term corresponds to term (2b) of Eq. 49 in Torrent et al.,
!! CMS 42, 337 (2008), including phase shifts as needed for orbital
!! magnetization:
!! $\langle\Phi_i|e^{i\mathbf{b.r}}v_H[n_{Zc}]e^{-i\mathbf{k.r}}|\Phi_j\rangle -
!!  \langle\tilde{\Phi}_i|e^{i\mathbf{b.r}}v_H[\tilde{n}_{Zc}]e^{-i\mathbf{k.r}}|\tilde{\Phi}_j\rangle $
!! where $\mathbf{b}$ is the bra shift vector and $\mathbf{k}$ is the ket shift vector, both are determined in
!! initberry.F90.
!!
!! PARENTS
!!      initberry
!!
!! CHILDREN
!!      initylmr,jbessel,realgaunt,simp_gen
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

 subroutine pawtwdij_2b(dtefield,gprimd,natom,ntypat,pawrad,pawtab,psps,typat)

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
#define ABI_FUNC 'pawtwdij_2b'
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
 integer :: bdir,bdx,bfor,bgnt,bl,bm,blm,bln,blmn,bsig,dij_ind,expibi_ind
 integer :: iatom,itypat,ilm2
 integer :: ir,kdir,kfor,kl,klm1,km,klm,kln,klmn,ksig,kdx
 integer :: lcmax,ll,mesh_size,mm,mu,ngnt
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 real(dp) :: bessarg,intgrl,jlx,jldx,jldx2,knorm
 complex(dpc) :: afac,vijsum
 logical :: need_conjg
!arrays
 integer,allocatable :: gntselect(:,:)
 real(dp) :: bb(3),bcart(3),kb(3),kcart(3),kbn(3),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: ff(:),realgnt(:),ylmk(:)
! the following is (-i)^L mod 4.
 complex(dpc),dimension(0:3) :: iml(0:3)=(/(1.0,0.0),(0.0,-1.0),(-1.0,0.0),(0.0,1.0)/)
! *************************************************************************

 lcmax = 2*psps%mpsang - 1 ! lcmax - 1 is highest angular momentum state used in expansion

 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr

 ABI_ALLOCATE(ylmk,(lcmax*lcmax))
 ABI_ALLOCATE(gntselect,((2*lcmax-1)**2,lcmax**2*(lcmax**2+1)/2))
 ABI_ALLOCATE(realgnt,((2*lcmax-1)**2*(lcmax)**4))
 call realgaunt(lcmax,ngnt,gntselect,realgnt)

 do itypat = 1, ntypat

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
           expibi_ind = dtefield%twind(bdx,kdx)
           if(expibi_ind < 0) then
             expibi_ind = -expibi_ind
             need_conjg = .true.
           else
             need_conjg = .false.
           end if
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
             ylmk(:)=zero; ylmk(1)=one/sqrt(four_pi)
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

               klm1 = RCC(blm,klm)

               vijsum = cmplx(zero,zero)

               do ll = abs(kl-bl),kl+bl,2
                 do ir = 1, mesh_size
                   bessarg = pawrad(itypat)%rad(ir)*knorm
                   call jbessel(jlx,jldx,jldx2,ll,0,bessarg)
                   ff(ir)=jlx*pawtab(itypat)%phi(ir,bln)*pawtab(itypat)%phi(ir,kln)*&
&                   pawtab(itypat)%VHnZC(ir) - &
&                   jlx*pawtab(itypat)%tphi(ir,bln)*pawtab(itypat)%tphi(ir,kln)*&
&                   pawtab(itypat)%VHntZC(ir)
                 end do
                 call simp_gen(intgrl,ff,pawrad(itypat))

                 do mm = -ll, ll
                   ilm2 = LMC(ll,mm)
                   bgnt = gntselect(ilm2,klm1)
                   if(bgnt /= 0) then
                     vijsum = vijsum + iml(mod(ll,4))*ylmk(ilm2)*realgnt(bgnt)*intgrl 
                   end if
                 end do ! end loop over mm
               end do ! end loop over ll
               
               vijsum = 4.d0*pi*vijsum
               do iatom = 1, natom ! store result multiplied by exp(-i*(k_k-k_b).I)
                 if(typat(iatom) == itypat) then
                   afac = cmplx(dtefield%expibi(1,iatom,expibi_ind),dtefield%expibi(2,iatom,expibi_ind))
                   if(need_conjg) afac=conjg(afac)
                   dtefield%twdij0(1,blmn,klmn,iatom,dij_ind) = dtefield%twdij0(1,blmn,klmn,iatom,dij_ind) + &
&                   real(vijsum*afac)
                   dtefield%twdij0(2,blmn,klmn,iatom,dij_ind) = dtefield%twdij0(2,blmn,klmn,iatom,dij_ind) + &
&                   aimag(vijsum*afac)
                 end if
               end do

             end do ! end loop over ket states
           end do ! end loop over bra states

         end do ! end loop over ksig
       end do ! end loop over kdir

     end do ! end loop over bsig
   end do ! end loop over bdir

   ABI_DEALLOCATE(ff)

 end do ! end loop on ntypat
 
 ABI_DEALLOCATE(gntselect)
 ABI_DEALLOCATE(realgnt)
 ABI_DEALLOCATE(ylmk)

 end subroutine pawtwdij_2b
!!***
