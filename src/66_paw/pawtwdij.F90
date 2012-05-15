!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawtwdij
!! NAME
!! pawtwdij
!!
!! FUNCTION
!! Compute the pseudopotential strengths Dij of the PAW non local operator,
!! including phase twists due to different k points as needed in orbital 
!! magnetization.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      initylmr,leave_new,realgaunt,sbf8,simp_gen,wrtout
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

subroutine pawtwdij(dtefield,gprimd,natom,nfftf,nspden,ntypat,&
&                   paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,typat,&
&                   vtrial,vxc)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_efield
 use m_radmesh, only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtwdij'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_65_psp
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: natom,nfftf,nspden,ntypat
 type(efield_type),intent(inout) :: dtefield
 type(pawang_type), intent(in) :: pawang
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: gprimd(3,3),vxc(nfftf,nspden),vtrial(nfftf,nspden)
 type(paw_an_type),intent(in) :: paw_an(natom)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
 type(pawrad_type), intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: bdir,bfor,bif,bigl,biglm,biglmt,biglt,bigltmin,bigltmax,bigm,bigmt,bsig
 integer :: clm,dij_ind,gij,gll
 integer :: iatom,ic,ijlm,ijlmax,ijlmin,ijlmn,ijln,ilmn,ils,ilslm,ir,isel,itypat
 integer :: jc,jlmn,kdir,kfor,kif,klmn,ksig
 integer :: lcmax,litl,litlm,litm,lm0,lm_size,lmn_size,lmn2_size
 integer :: mesh_size,mm,mu,nfgd,ngnt,twind
 integer ::  ylmr_normchoice,ylmr_npts,ylmr_option
 real(dp) :: bessarg,gylm,knorm,vr,vxcij1
 complex(dpc) :: tvxc_3a
 character(len=500) :: message
 logical :: need_conj
!arrays
 integer,allocatable :: gntselect(:,:)
 real(dp) :: bb(3),bcart(3),dijhat(2),eijkl(2),kb(3),kbn(3),kcart(3),phtw(2),ro(2)
 real(dp) :: ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: ff(:),j_bessel(:,:),prod(:,:),realgnt(:),ylmk(:)
! the following is (-i)^L mod 4.
 complex(dpc),dimension(0:3) :: iml(0:3)=(/(1.0,0.0),(0.0,-1.0),(-1.0,0.0),(0.0,1.0)/)

! *************************************************************************

!DBG_ENTER("COLL")

 dtefield%twdij(:,:,:,:,:) = zero
 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr

!load twdij0
 dtefield%twdij(:,:,:,:,:) = dtefield%twdij0(:,:,:,:,:)

 do bdir = 1, 3
   do bsig = -1, 1, 2

     bb(:) = bsig*dtefield%dkvecs(:,bdir) ! bra vector
     do mu=1,3
       bcart(mu)=dot_product(bb(:),gprimd(mu,:))
     end do

     do kdir = 1, 3
       if (kdir /= bdir) then
         do ksig = -1, 1, 2

           bfor = (3-bsig)/2
           kfor = (3-ksig)/2
           
           bif = 2*bdir-bfor+1
           kif = 2*kdir-kfor+1
           dij_ind = dtefield%indhk(bif,kif)
           
           twind = dtefield%twind(bif,kif)
           need_conj = (twind < 0)
           twind = abs(twind)

           kb(:) = dtefield%dkvecs(:,kdir)  ! ket vector
           do mu=1,3
             kcart(mu)=dot_product(kb(:),gprimd(mu,:))
           end do
!          form b_k - b_b vector
           kcart(1:3) = kcart(1:3) - bcart(1:3)
           knorm=dsqrt(dot_product(kcart,kcart))
           
!          loop over atoms
           do iatom = 1, natom
             itypat = typat(iatom)
             lm_size = paw_an(iatom)%lm_size
             ABI_ALLOCATE(prod,(2,lm_size))
             lcmax = pawtab(itypat)%l_size
             lmn_size = pawtab(itypat)%lmn_size
             ABI_ALLOCATE(ylmk,(lcmax*lcmax))
             ABI_ALLOCATE(gntselect,((2*lcmax-1)**2,lcmax**2*(lcmax**2+1)/2))
             ABI_ALLOCATE(realgnt,((2*lcmax-1)**2*(lcmax)**4))
             call realgaunt(lcmax,ngnt,gntselect,realgnt)

             mesh_size = pawrad(itypat)%mesh_size
             ABI_ALLOCATE(ff,(mesh_size))
             ABI_ALLOCATE(j_bessel,(mesh_size,lcmax))

             lmn2_size = pawtab(itypat)%lmn2_size
             nfgd = pawfgrtab(iatom)%nfgd

             if (knorm < tol12) then
               kbn(:) = zero
               ylmk(:) = zero
               ylmk(1) = 1.d0/sqrt(four_pi)
             else
               kbn(:) = kcart(:)/knorm ! unit vector in kb direction
               call initylmr(lcmax,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,kbn,ylmk(:),ylmgr)
             end if
             knorm = two_pi*knorm ! re-normed kb for calls to bessel fnc
             
!            compute bessel functions here
             do ir = 1, mesh_size
               bessarg = knorm*pawrad(itypat)%rad(ir)
               call sbf8(lcmax,bessarg,j_bessel(ir,:))
             end do
             
             do ilmn = 1, lmn_size
               do jlmn = 1, ilmn
                 
                 ijlmn = RCC(ilmn,jlmn)
                 
                 
!                add eijkl terms
                 
                 do isel = 1, pawrhoij(iatom)%nrhoijsel
                   klmn=pawrhoij(iatom)%rhoijselect(isel)
                   if (pawrhoij(iatom)%cplex == 1) then
                     ro(1) = pawrhoij(iatom)%rhoijp(klmn,1)*pawtab(itypat)%dltij(klmn) ! nspden == 1 implicitly
                     eijkl(1) = dtefield%tweijkl(1,ijlmn,klmn,iatom,twind)
                     eijkl(2) = dtefield%tweijkl(2,ijlmn,klmn,iatom,twind)
                     if(need_conj) eijkl(2) = -eijkl(2)
                     dtefield%twdij(1,ilmn,jlmn,iatom,dij_ind) = dtefield%twdij(1,ilmn,jlmn,iatom,dij_ind) &
&                     + ro(1)*eijkl(1)
                     dtefield%twdij(2,ilmn,jlmn,iatom,dij_ind) = dtefield%twdij(2,ilmn,jlmn,iatom,dij_ind) &
&                     + ro(1)*eijkl(2)
                     if (jlmn /= ilmn) then
                       dtefield%twdij(1,jlmn,ilmn,iatom,dij_ind) = dtefield%twdij(1,jlmn,ilmn,iatom,dij_ind)&
&                       + ro(1)*eijkl(1)
                       dtefield%twdij(2,jlmn,ilmn,iatom,dij_ind) = dtefield%twdij(2,jlmn,ilmn,iatom,dij_ind)&
&                       + ro(1)*eijkl(2)
                     end if
                   else
                     write(message,'(a,a,i4,a,a)')ch10,&
&                     ' Bug in pawtwdij : pawrhoij(iatom)%cplex = ',pawrhoij(iatom)%cplex,ch10,&
&                     ' but this case not yet coded '
                     call wrtout(std_out,message,'COLL')
                     call leave_new('COLL')
                   end if
                 end do ! end loop over nrhoijsel
                 
!                add Dij_hat
                 prod(:,:) = zero
                 dijhat(:) = zero
                 do ilslm = 1, lm_size
                   do ic = 1, nfgd
                     jc = pawfgrtab(iatom)%ifftsph(ic)
                     vr = vtrial(jc,1) - vxc(jc,1) ! implicitly nspden = 1 and usexcnhat = 0
                     gylm = pawfgrtab(iatom)%gylm(ic,ilslm)
                     phtw(1) = dtefield%expibr(1,iatom,ic,twind)
                     phtw(2) = dtefield%expibr(2,iatom,ic,twind)
                     prod(1,ilslm) = prod(1,ilslm) + vr*gylm*phtw(1)
                     prod(2,ilslm) = prod(2,ilslm) + vr*gylm*phtw(2)
                   end do ! end loop over nfgd
                   if (need_conj) prod(2,ilslm) = -prod(2,ilslm)
                 end do ! end loop over lm_size
                 ijlm = pawtab(itypat)%indklmn(1,ijlmn)
                 ijln = pawtab(itypat)%indklmn(2,ijlmn)
                 ijlmin = pawtab(itypat)%indklmn(3,ijlmn)
                 ijlmax = pawtab(itypat)%indklmn(4,ijlmn)
                 do ils=ijlmin,ijlmax,2
                   lm0 = ils*ils+ils+1
                   do mm=-ils,ils
                     ilslm=lm0+mm
                     isel=pawang%gntselect(ilslm,ijlm)
                     if(isel > 0) then
                       dijhat(1) = dijhat(1) + prod(1,ilslm)*pawtab(itypat)%qijl(ilslm,ijlmn)
                       dijhat(2) = dijhat(2) + prod(2,ilslm)*pawtab(itypat)%qijl(ilslm,ijlmn)
                     end if
                   end do ! end loop over mm
                 end do ! end loop over ils
                 
                 dtefield%twdij(1,ilmn,jlmn,iatom,dij_ind) = dtefield%twdij(1,ilmn,jlmn,iatom,dij_ind) &
&                 + dijhat(1)
                 dtefield%twdij(2,ilmn,jlmn,iatom,dij_ind) = dtefield%twdij(2,ilmn,jlmn,iatom,dij_ind) &
&                 + dijhat(2)
                 if (jlmn /= ilmn) then
                   dtefield%twdij(1,jlmn,ilmn,iatom,dij_ind) = dtefield%twdij(1,jlmn,ilmn,iatom,dij_ind) &
&                   + dijhat(1)
                   dtefield%twdij(2,jlmn,ilmn,iatom,dij_ind) = dtefield%twdij(2,jlmn,ilmn,iatom,dij_ind) &
&                   + dijhat(2)
                 end if

!                construct twisted vxc terms
!                term 3a is on-site vxc term from Torrent et al, Comp. Mater. Sci. 42, 337 (2008).
!                term 3b is the vxc interaction with nhat but due to implicit usexcnhat = 0 that term vanishes here
!                (see for example pawdij.F90 line 527 and following)

                 tvxc_3a = cmplx(zero,zero)
                 
                 bigltmin = pawtab(itypat)%indklmn(3,ijlmn)
                 bigltmax = pawtab(itypat)%indklmn(4,ijlmn)

                 do biglt = bigltmin, bigltmax, 2
                   do litl = 0, bigltmax
                     do bigl = 0, bigltmax

                       if( biglt + litl - bigl >= 0 .and. &
&                       biglt - litl + bigl >= 0 .and. &
&                       -biglt + litl + bigl >= 0 ) then
                         
                         do bigm = -bigl, bigl
                           biglm = LMC(bigl,bigm)

                           if (paw_an(iatom)%lmselect(biglm)) then
                             ff(1:mesh_size) = j_bessel(1:mesh_size,litl+1)*&
&                             (paw_an(iatom)%vxc1(1:mesh_size,biglm,1)*&
&                             pawtab(itypat)%phiphj(1:mesh_size,ijln) - &
&                             paw_an(iatom)%vxct1(1:mesh_size,biglm,1)*&
&                             pawtab(itypat)%tphitphj(1:mesh_size,ijln))
                             call simp_gen(vxcij1,ff,pawrad(itypat))
                             
                             do bigmt = -biglt, biglt
                               biglmt = LMC(biglt,bigmt)
                               gij = gntselect(biglmt,ijlm)
                               if(gij == 0) cycle
                               
                               do litm = -litl, litl
                                 litlm = LMC(litl,litm)
                                 clm = RCC(litlm,biglm)
                                 gll = gntselect(biglmt,clm)
                                 if(gll == 0) cycle
                                 
                                 tvxc_3a = tvxc_3a+four_pi*iml(mod(litl,4))*ylmk(litlm)*&
&                                 realgnt(gij)*realgnt(gll)*vxcij1
                               end do ! end loop over litm
                             end do ! end loop over bigmt

                           end if ! end selection on nonzero vxc component

                         end do ! end loop on bigm

                       end if ! end selection on litl bigl biglt Gaunt

                     end do ! end loop over bigl
                   end do ! end loop over litl
                 end do ! end loop over biglt

                 dtefield%twdij(1,ilmn,jlmn,iatom,dij_ind) = dtefield%twdij(1,ilmn,jlmn,iatom,dij_ind) &
&                 + real(tvxc_3a)
                 dtefield%twdij(2,ilmn,jlmn,iatom,dij_ind) = dtefield%twdij(2,ilmn,jlmn,iatom,dij_ind) &
&                 + aimag(tvxc_3a)
                 if (jlmn /= ilmn) then
                   dtefield%twdij(1,jlmn,ilmn,iatom,dij_ind) = dtefield%twdij(1,jlmn,ilmn,iatom,dij_ind) &
&                   + real(tvxc_3a)
                   dtefield%twdij(2,jlmn,ilmn,iatom,dij_ind) = dtefield%twdij(2,jlmn,ilmn,iatom,dij_ind) &
&                   + aimag(tvxc_3a)
                 end if

               end do ! end loop over jlmn
             end do ! end loop over ilmn
             
             ABI_DEALLOCATE(ff)
             ABI_DEALLOCATE(gntselect)
             ABI_DEALLOCATE(j_bessel)
             ABI_DEALLOCATE(prod)
             ABI_DEALLOCATE(realgnt)
             ABI_DEALLOCATE(ylmk)
             
           end do ! end loop over iatom
         end do ! end loop over ksig
       end if ! end check that kdir /= bdir
     end do ! end loop over kdir
   end do ! end loop over bsig
 end do ! end loop over bdir

!DBG_EXIT("COLL")

end subroutine pawtwdij
!!***
