!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawtwdij_2f
!! NAME
!! pawtwdij_2f
!!
!! FUNCTION
!! compute phase-twisted contribution to ${D}^1_{ij}-\tilde{D}^1_{ij}$ due to 
!! Hatree potential of nhat interacting with compensation charge.
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
!!  typat = typat(natom) list of atom types
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtefield <type(efield_type)> = dtefield%tweijkl is updated
!!
!! NOTES
!! This term corresponds to term (2f) of Eq. 49 in Torrent et al.,
!! CMS 42, 337 (2008), including phase shifts as needed for orbital
!! magnetization:
!! $\sum_{LM}\int_{\Omega_R}e^{i\mathbf{b.r}}v_H[\tilde{n}_{Zc}]\hat{Q}^{LM}_{ij}(\mathbf{r})e^{-i\mathbf{k.r}}$
!! where $\mathbf{b}$ is the bra shift vector and $\mathbf{k}$ is the ket shift vector, both are determined in
!! initberry.F90.
!!
!! PARENTS
!!      initberry
!!
!! CHILDREN
!!      initylmr,poisson,realgaunt,sbf8,simp_gen
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

 subroutine pawtwdij_2f(dtefield,gprimd,natom,ntypat,pawrad,pawtab,typat)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_efield
 use m_radmesh,  only : poisson, simp_gen, deducer0

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtwdij_2f'
 use interfaces_32_util
 use interfaces_65_psp
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 type(efield_type),intent(inout) :: dtefield

!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: gprimd(3,3)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: bdir,bigl,bigm,biglm,biglmin,biglmax,bsig,clm,gij,gkl,glm
 integer :: eijkl_ind,iatom,ijlmn,ijlm,ijln,ir,itypat
 integer :: kdir,kllmn,kllm,klln,lcmax
 integer :: litl,litm,litlm,lt,ltmax,ltmin,ltlm,mt,mesh_size,mu,ngnt
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 real(dp) :: bessarg,gjv,knorm,qq
 complex(dpc) :: tijkl
!arrays
 integer,allocatable :: gntselect(:,:)
 real(dp) :: bb(3),bcart(3),kb(3),kcart(3),kbn(3),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: ff(:),j_bessel(:,:),realgnt(:),rvl(:),ylmk(:)
! the following is (-i)^L mod 4.
 complex(dpc),dimension(0:3) :: iml(0:3)=(/(1.0,0.0),(0.0,-1.0),(-1.0,0.0),(0.0,1.0)/)
! *************************************************************************

 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr


 do itypat = 1, ntypat

   lcmax = 2*pawtab(itypat)%l_size + 1 ! lcmax - 1 is highest angular momentum state used in expansion
   ABI_ALLOCATE(ylmk,(lcmax*lcmax))
   ABI_ALLOCATE(gntselect,((2*lcmax-1)**2,lcmax**2*(lcmax**2+1)/2))
   ABI_ALLOCATE(realgnt,((2*lcmax-1)**2*(lcmax)**4))
   call realgaunt(lcmax,ngnt,gntselect,realgnt)

   mesh_size = pawrad(itypat)%mesh_size
   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(j_bessel,(mesh_size,lcmax))
   ABI_ALLOCATE(rvl,(mesh_size))

   do bdir = 1, 3
     do bsig = -1, 1, 2

       eijkl_ind = 2*bdir-1+(bsig+1)/2

       bb(:) = bsig*dtefield%dkvecs(:,bdir) ! bra vector
       do mu=1,3
         bcart(mu)=dot_product(bb(:),gprimd(mu,:))
       end do
       kdir = mod(bdir,3)+1
       kb(:) = dtefield%dkvecs(:,kdir)  ! ket vector
       do mu=1,3
         kcart(mu)=dot_product(kb(:),gprimd(mu,:))
       end do
!      form b_k - b_b vector
       kcart(1:3) = kcart(1:3) - bcart(1:3)
       knorm=dsqrt(dot_product(kcart,kcart))
       if (knorm < tol12) then
         kbn(:) = zero
         ylmk(:) = zero; ylmk(1) = 1.d0/sqrt(four_pi)
       else
         kbn(:) = kcart(:)/knorm ! unit vector in kb direction
         call initylmr(lcmax,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,kbn,ylmk(:),ylmgr)
       end if 
       knorm = two_pi*knorm ! re-normed kb for calls to bessel fnc

!      compute bessel functions here
       do ir = 1, mesh_size
         bessarg = knorm*pawrad(itypat)%rad(ir)
         call sbf8(lcmax,bessarg,j_bessel(ir,:))
       end do

       do ijlmn = 1, pawtab(itypat)%lmn2_size
         ijlm = pawtab(itypat)%indklmn(1,ijlmn)
         ijln = pawtab(itypat)%indklmn(2,ijlmn)
         biglmin = pawtab(itypat)%indklmn(3,ijlmn)
         biglmax = pawtab(itypat)%indklmn(4,ijlmn)

         do kllmn = 1, pawtab(itypat)%lmn2_size
           kllm = pawtab(itypat)%indklmn(1,kllmn)
           klln = pawtab(itypat)%indklmn(2,kllmn)
           ltmin = pawtab(itypat)%indklmn(3,kllmn)
           ltmax = pawtab(itypat)%indklmn(4,kllmn)

           tijkl = cmplx(zero,zero)

           do lt = ltmin, ltmax, 2
             
             ff(1:mesh_size) = four_pi*(pawrad(itypat)%rad(1:mesh_size)**2)*&
&             pawtab(itypat)%shapefunc(1:mesh_size,lt+1)
             call poisson(ff,lt,qq,pawrad(itypat),rvl)

             do bigl = biglmin, biglmax, 2

               do litl = 0, lcmax-1

                 if( lt+bigl-litl >= 0 .and. &
&                 lt-bigl+litl >= 0 .and. &
&                 -lt+bigl+litl >= 0 ) then
                   
!                  get the radial integral 
!                  remember that rvl = r*v_L, need to divide out extra factor of r 
!                  so only one factor of r in Jacobian instead of 2
                   ff(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)*rvl(1:mesh_size)*&
&                   pawtab(itypat)%shapefunc(1:mesh_size,bigl+1)*j_bessel(1:mesh_size,litl+1) 
                   call simp_gen(gjv,ff,pawrad(itypat))
                   
                   do bigm = -bigl, bigl
                     biglm = LMC(bigl,bigm)
                     gij = gntselect(biglm,ijlm)
                     if (gij == 0) cycle

                     do mt = -lt, lt
                       ltlm = LMC(lt,mt)
                       gkl = gntselect(ltlm,kllm)
                       if (gkl == 0) cycle

                       do litm = -litl, litl
                         litlm = LMC(litl,litm)
                         clm = RCC(biglm,ltlm)
                         glm = gntselect(litlm,clm)
                         if (glm == 0) cycle

                         tijkl = tijkl - four_pi*iml(mod(litl,4))*ylmk(litlm)*&
&                         pawtab(itypat)%qijl(biglm,ijlmn)*pawtab(itypat)%qijl(ltlm,kllmn)*&
&                         realgnt(glm)*gjv

                       end do ! end loop over litm
                     end do ! end loop over mt
                   end do ! end loop over bigm
                 end if ! end check on triangle condition
               end do ! end loop over litl
             end do ! end loop over bigl
           end do ! end loop over lt


           do iatom = 1, natom ! store result
             if(typat(iatom) == itypat) then
               dtefield%tweijkl(1,ijlmn,kllmn,iatom,eijkl_ind) = & 
&               dtefield%tweijkl(1,ijlmn,kllmn,iatom,eijkl_ind) + &
&               real(tijkl)
               dtefield%tweijkl(2,ijlmn,kllmn,iatom,eijkl_ind) = & 
&               dtefield%tweijkl(2,ijlmn,kllmn,iatom,eijkl_ind) + &
&               aimag(tijkl)
             end if
           end do ! end loop over atoms in cell

         end do ! end loop over kl states
       end do ! end loop over ij states

     end do ! end loop over bsig
   end do ! end loop over bdir

   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(j_bessel)
   ABI_DEALLOCATE(rvl)
   ABI_DEALLOCATE(ylmk)
   ABI_DEALLOCATE(gntselect)
   ABI_DEALLOCATE(realgnt)

 end do ! end loop on ntypat


 end subroutine pawtwdij_2f
!!***
