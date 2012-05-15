!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawalloc
!! NAME
!! pawalloc
!!
!! FUNCTION
!! Allocate or deallocate datastructures used for PAW
!! at the level of the driving routine (driver.F90)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  idtset=index of the current dataset
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid_vl=dimension of q (or G) grid for Vloc (array vlspl)
!!  npsp=number of pseudopotentials
!!  option=1: allocate PAW datastructures for a new dataset
!!         2: deallocate PAW datastructures depending on paw_size (pawrad,pawtab) for the current dataset
!!         3: deallocate all PAW datastructures (pawang,pawrad,pawtab) for the current dataset
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!                 pseudopotential file header, as well as the psp file name
!!
!! SIDE EFFECTS
!! Allocated/deallocated:
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(paw_size) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(paw_size) <type(pawtab_type)>=paw tabulated starting data
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      destroy_pawang
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawalloc(dtset,idtset,mpsang,mqgrid_vl,npsp,option,paw_size,paw_size_old,&
&                   pawang,pawrad,pawtab,pspheads)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_paw_toolbox

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawalloc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idtset,mpsang,mqgrid_vl,npsp,option,paw_size,paw_size_old
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(inout) :: pawang
!arrays
 type(pawrad_type),intent(inout) :: pawrad(paw_size)
 type(pawtab_type),intent(inout) :: pawtab(paw_size)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 integer,parameter :: mqgrid_shp=300
 integer :: angl_size_new,basis_size_new,ierr,ij_size_new,itypat,gnt_option_new,l_max_new
 integer :: l_size_new,l_size_max_new,lexexch_new,lmn_size_new,lmn2_size_new,lpawu_new
 integer :: mesh_size_new,mqgrid_shp_new,nsym_new,pawspnorb_new,ylm_size_new
 integer :: ngnt_new,ntheta_new,nphi_new
 logical :: need_kij_new,need_nabla_new,need_vhntzc_new,need_vhnzc_new
 logical :: test_alloc,test_angl_size,test_basis_size,test_ij_size,test_gnt_option,test_l_max
 logical :: test_l_size,test_l_size_max,test_lexexch,test_lmn_size,test_lmn2_size,test_kij,test_lpawu
 logical :: test_mesh_size,test_mqgrid,test_mqgrid_shp,test_nabla,test_new,test_nsym,test_spnorb
 logical :: test_vhntzc,test_vhnzc,test_ylm_size
!arrays

! *********************************************************************

 DBG_ENTER("COLL")

!Nothing to do if not PAW
 if (paw_size==0) return

!=====================================================================
!========================= ALLOCATIONS ===============================
!=====================================================================
 if (option==1) then

   test_new=(paw_size/=paw_size_old)
   test_alloc=((idtset/=1).and.(paw_size==paw_size_old))

   l_max_new=mpsang
   l_size_max_new=2*l_max_new-1
   ylm_size_new=0
   if(dtset%pawxcdev==0) then
     if (dtset%xclevel==2) ylm_size_new=(l_size_max_new+1)**2
     if (dtset%xclevel/=2) ylm_size_new=l_size_max_new**2
   end if
   gnt_option_new=1;if (dtset%pawxcdev==2.or.(dtset%pawxcdev==1.and.dtset%positron/=0)) gnt_option_new=2
   angl_size_new=0;if(dtset%pawxcdev==0) angl_size_new=dtset%pawntheta*dtset%pawnphi
   nsym_new=dtset%nsym
   pawspnorb_new=0;if (dtset%pawspnorb>0) pawspnorb_new=1
   mqgrid_shp_new=0
   ngnt_new=0
   ntheta_new=0
   nphi_new=0
   need_kij_new=(dtset%positron/=0)
   need_nabla_new=.false.
   need_vhntzc_new=(abs(dtset%berryopt)==5)
   need_vhnzc_new=(abs(dtset%berryopt)==5)
   if (dtset%optdriver==0.and. &
&   ((dtset%iprcel>=20.and.dtset%iprcel<70).or.dtset%iprcel>=80)) mqgrid_shp_new=mqgrid_shp

   if (idtset==1) then
     pawang%angl_size=-1
     pawang%l_max=-1
     pawang%l_size_max=-1
     pawang%ylm_size=-1
     pawang%nsym=-1
     pawang%gnt_option=-1
     pawang%use_ls_ylm=-1
     nullify(pawang%anginit)
     nullify(pawang%angwgth)
     nullify(pawang%zarot)
     nullify(pawang%gntselect)
     nullify(pawang%realgnt)
     nullify(pawang%ylmr)
     nullify(pawang%ylmrgr)
     nullify(pawang%ls_ylm)
   end if

   test_l_max=(l_max_new/=pawang%l_max)
   test_angl_size=(angl_size_new/=pawang%angl_size)
   test_l_size_max=(l_size_max_new/=pawang%l_size_max)
   test_ylm_size=(ylm_size_new/=pawang%ylm_size)
   test_gnt_option=(gnt_option_new/=pawang%gnt_option)
   test_spnorb=(pawspnorb_new/=pawang%use_ls_ylm)
   test_nsym=(nsym_new/=pawang%nsym)

   do itypat=1,dtset%ntypat

     basis_size_new=pspheads(itypat)%pawheader%basis_size
     lmn_size_new  =pspheads(itypat)%pawheader%lmn_size
     l_size_new    =pspheads(itypat)%pawheader%l_size
     mesh_size_new =pspheads(itypat)%pawheader%mesh_size
     lmn2_size_new =lmn_size_new*(lmn_size_new+1)/2
     ij_size_new   =basis_size_new*(basis_size_new+1)/2
     lpawu_new     =dtset%lpawu(itypat)
     lexexch_new   =dtset%lexexch(itypat)

     if (idtset==1) then
       test_basis_size=.true.;test_ij_size  =.true.
       test_mesh_size =.true.;test_l_size   =.true.
       test_mqgrid_shp=.true.;test_lmn_size =.true.
       test_lmn2_size =.true.;test_lpawu    =.true.
       test_lexexch   =.true.;test_mqgrid   =.true.
       test_kij       =.true.;test_nabla    =.true.
       test_vhntzc    =.true.; test_vhnzc   =.true.
     else
       test_basis_size=(basis_size_new/=pawtab(itypat)%basis_size)
       test_ij_size=(ij_size_new/=pawtab(itypat)%ij_size)
       test_mesh_size=(mesh_size_new/=pawtab(itypat)%mesh_size)
       test_l_size=(l_size_new/=pawtab(itypat)%l_size)
       test_mqgrid_shp=(mqgrid_shp_new/=pawtab(itypat)%mqgrid_shp)
       test_lmn_size=(lmn_size_new/=pawtab(itypat)%lmn_size)
       test_lmn2_size=(lmn2_size_new/=pawtab(itypat)%lmn2_size)
       test_lpawu=(lpawu_new/=pawtab(itypat)%lpawu)
       test_lexexch=(lexexch_new/=pawtab(itypat)%lexexch)
       test_mqgrid=(mqgrid_vl/=pawtab(itypat)%mqgrid)
       test_kij=(need_kij_new.neqv.(pawtab(itypat)%has_kij>0))
       test_nabla=(need_nabla_new.neqv.(pawtab(itypat)%has_nabla>0))
       test_vhntzc=(need_vhntzc_new.neqv.(pawtab(itypat)%has_vhntzc>0))
       test_vhnzc=(need_vhnzc_new.neqv.(pawtab(itypat)%has_vhnzc>0))
     end if

!    Reallocate arrays depending on mesh_size and basis_size
     if (test_new.or.test_ij_size.or.test_basis_size) then
       if(test_alloc)  then
         ABI_DEALLOCATE(pawtab(itypat)%tphi)
         ABI_DEALLOCATE(pawtab(itypat)%phi)
       end if
       ABI_ALLOCATE(pawtab(itypat)%tphi,(mesh_size_new,basis_size_new))
       ABI_ALLOCATE(pawtab(itypat)%phi ,(mesh_size_new,basis_size_new))
     end if

!    Reallocate arrays depending on mesh_size and ij_size
     if (test_new.or.test_mesh_size.or.test_ij_size) then
       if(test_alloc)  then
         ABI_DEALLOCATE(pawtab(itypat)%tphitphj)
         ABI_DEALLOCATE(pawtab(itypat)%phiphj)
       end if
       ABI_ALLOCATE(pawtab(itypat)%tphitphj,(mesh_size_new,ij_size_new))
       ABI_ALLOCATE(pawtab(itypat)%phiphj  ,(mesh_size_new,ij_size_new))
     end if

!    Reallocate arrays depending on mesh_size and l_size
     if (test_new.or.test_mesh_size.or.test_l_size) then
       if(test_alloc)  then
         ABI_DEALLOCATE(pawtab(itypat)%shapefunc)
         ABI_DEALLOCATE(pawtab(itypat)%dshpfunc)
       end if
       ABI_ALLOCATE(pawtab(itypat)%shapefunc,(mesh_size_new,l_size_new))
       ABI_ALLOCATE(pawtab(itypat)%dshpfunc,(mesh_size_new,l_size_new,0))
     end if

!    Reallocate arrays depending on mqgrid_shp and l_size
     if (test_new.or.test_mqgrid_shp.or.test_l_size) then
       if(test_alloc.and.pawtab(itypat)%mqgrid_shp>0)  then
         ABI_DEALLOCATE(pawtab(itypat)%shapefncg)
       end if
       if (mqgrid_shp_new>0)  then
         ABI_ALLOCATE(pawtab(itypat)%shapefncg,(mqgrid_shp_new,2,l_size_new))
       end if
     end if

!    Reallocate arrays depending on l_size and lmn2_size
     if (test_new.or.test_l_size.or.test_lmn2_size) then
       if(test_alloc)  then
         ABI_DEALLOCATE(pawtab(itypat)%qijl)
       end if
       ABI_ALLOCATE(pawtab(itypat)%qijl,(l_size_new*l_size_new,lmn2_size_new))
     end if

!    Reallocate arrays depending on lpawu, lexexch and lmn2_size
     if (test_new.or.test_lmn2_size.or.test_lpawu.or.test_lexexch) then
       if(test_alloc.and.(pawtab(itypat)%lpawu>=0.or.pawtab(itypat)%lexexch>=0))  then
         ABI_DEALLOCATE(pawtab(itypat)%klmntomn)
       end if
       if (lpawu_new>=0.or.lexexch_new>=0)  then
         ABI_ALLOCATE(pawtab(itypat)%klmntomn,(4,lmn2_size_new))
       end if
     end if

!    Reallocate arrays depending on lpawu and lmn2_size
     if (test_new.or.test_lpawu) then
       if(test_alloc.and.(pawtab(itypat)%lpawu>=0))  then
         ABI_DEALLOCATE(pawtab(itypat)%vee)
       end if
       if (lpawu_new>=0)  then
         sz1=2*dtset%lpawu(itypat)+1
         ABI_ALLOCATE(pawtab(itypat)%vee,(sz1,sz1,sz1,sz1))
       end if
     end if

!    Reallocate arrays depending on lexexch and lmn2_size
     if (test_new.or.test_lexexch) then
       if(test_alloc.and.(pawtab(itypat)%lexexch>=0))  then
         ABI_DEALLOCATE(pawtab(itypat)%vex)
         ABI_DEALLOCATE(pawtab(itypat)%fk)
       end if
       if (lexexch_new>=0) then
         sz1=2*dtset%lexexch(itypat)+1
         ABI_ALLOCATE(pawtab(itypat)%vex,(sz1,sz1,sz1,sz1,4))
         ABI_ALLOCATE(pawtab(itypat)%fk,(6,4))
       end if
     end if

!    Reallocate arrays depending on l_size
     if (test_new.or.test_l_size) then
       if(test_alloc)  then
         ABI_DEALLOCATE(pawtab(itypat)%gnorm)
         ABI_DEALLOCATE(pawtab(itypat)%shape_alpha)
         ABI_DEALLOCATE(pawtab(itypat)%shape_q)
       end if
       ABI_ALLOCATE(pawtab(itypat)%gnorm,(l_size_new))
       ABI_ALLOCATE(pawtab(itypat)%shape_alpha,(2,l_size_new))
       ABI_ALLOCATE(pawtab(itypat)%shape_q,(2,l_size_new))
     end if

!    Reallocate arrays depending on lmn_size and has_nabla
     if (test_new.or.test_lmn_size.or.test_nabla) then
       if(test_alloc.and.pawtab(itypat)%has_nabla>0)  then
         ABI_DEALLOCATE(pawtab(itypat)%nabla_ij)
       end if
       if (need_nabla_new)  then
         ABI_ALLOCATE(pawtab(itypat)%nabla_ij,(3,lmn_size_new,lmn_size_new))
       end if
     end if

!    Reallocate arrays depending on lmn2_size and has_kij
     if (test_new.or.test_lmn_size.or.test_kij) then
       if(test_alloc.and.pawtab(itypat)%has_kij>0)  then
         ABI_DEALLOCATE(pawtab(itypat)%kij)
       end if
       if (need_kij_new)  then
         ABI_ALLOCATE(pawtab(itypat)%kij,(lmn2_size_new))
       end if
     end if

!    Reallocate arrays depending on mesh_size and has_vhntzc
     if (test_new.or.test_mesh_size.or.test_vhntzc) then
       if(test_alloc.and.pawtab(itypat)%has_vhntzc>0)  then
         ABI_DEALLOCATE(pawtab(itypat)%VHntZC)
       end if
       if (need_vhntzc_new)  then
         ABI_ALLOCATE(pawtab(itypat)%VHntZC,(mesh_size_new))
       end if
     end if

!    Reallocate arrays depending on mesh_size and has_vhnzc
     if (test_new.or.test_mesh_size.or.test_vhnzc) then
       if(test_alloc.and.pawtab(itypat)%has_vhnzc>0)  then
         ABI_DEALLOCATE(pawtab(itypat)%VHnZC)
       end if
       if (need_vhnzc_new)  then
         ABI_ALLOCATE(pawtab(itypat)%VHnZC,(mesh_size_new))
       end if
     end if

!    Reallocate arrays depending on lmn2_size
     if (test_new.or.test_lmn2_size) then
       if(test_alloc) then
         ABI_DEALLOCATE(pawtab(itypat)%eijkl)
         ABI_DEALLOCATE(pawtab(itypat)%dij0)
         ABI_DEALLOCATE(pawtab(itypat)%dltij)
         ABI_DEALLOCATE(pawtab(itypat)%rhoij0)
         ABI_DEALLOCATE(pawtab(itypat)%sij)
         ABI_DEALLOCATE(pawtab(itypat)%indklmn)
       end if
       if (idtset==1) nullify(pawtab(itypat)%kmix)
       if (associated(pawtab(itypat)%kmix))  then
         ABI_DEALLOCATE(pawtab(itypat)%kmix)
         ierr = ABI_ALLOC_STAT
       end if
       ABI_ALLOCATE(pawtab(itypat)%eijkl,(lmn2_size_new,lmn2_size_new))
       ABI_ALLOCATE(pawtab(itypat)%dij0,(lmn2_size_new))
       ABI_ALLOCATE(pawtab(itypat)%dltij,(lmn2_size_new))
       ABI_ALLOCATE(pawtab(itypat)%rhoij0,(lmn2_size_new))
       ABI_ALLOCATE(pawtab(itypat)%sij,(lmn2_size_new))
       ABI_ALLOCATE(pawtab(itypat)%indklmn,(6,lmn2_size_new))
     end if

!    Reallocate arrays depending on mesh_size
     if (test_new.or.test_mesh_size) then
       if(test_alloc)  then
         ABI_DEALLOCATE(pawtab(itypat)%coredens)
         ABI_DEALLOCATE(pawtab(itypat)%tcoredens)
         ABI_DEALLOCATE(pawrad(itypat)%rad)
         ABI_DEALLOCATE(pawrad(itypat)%radfact)
         ABI_DEALLOCATE(pawrad(itypat)%simfact)
         ABI_DEALLOCATE(pawtab(itypat)%rad_for_spline)
       end if
       ABI_ALLOCATE(pawtab(itypat)%coredens ,(mesh_size_new))
       ABI_ALLOCATE(pawtab(itypat)%tcoredens,(mesh_size_new))
       ABI_ALLOCATE(pawrad(itypat)%rad      ,(mesh_size_new))
       ABI_ALLOCATE(pawrad(itypat)%radfact  ,(mesh_size_new))
       ABI_ALLOCATE(pawrad(itypat)%simfact  ,(mesh_size_new))
       ABI_ALLOCATE(pawtab(itypat)%rad_for_spline,(0))
     end if

!    Reallocate arrays depending on mqgrid_vl
     if (test_new.or.test_mqgrid) then
       if(test_alloc) then
         ABI_DEALLOCATE(pawtab(itypat)%tcorespl)
         if (pawtab(itypat)%usetvale>0)  then
           ABI_DEALLOCATE(pawtab(itypat)%tvalespl)
         end if
       end if
       ABI_ALLOCATE(pawtab(itypat)%tcorespl,(mqgrid_vl,2))
       if (pspheads(itypat)%pawheader%pawver>=4)  then
         ABI_ALLOCATE(pawtab(itypat)%tvalespl,(mqgrid_vl,2))
       end if
     end if

!    Reallocate arrays depending on mqgrid_shp and l_size
     if (test_new.or.test_mqgrid_shp)  then
       if(test_alloc.and.pawtab(itypat)%mqgrid_shp>0)  then
         ABI_DEALLOCATE(pawtab(itypat)%qgrid_shp)
       end if
       if (mqgrid_shp_new>0)  then
         ABI_ALLOCATE(pawtab(itypat)%qgrid_shp,(mqgrid_shp_new))
       end if
     end if

     pawtab(itypat)%basis_size=basis_size_new
     pawtab(itypat)%ij_size=ij_size_new
     pawtab(itypat)%lpawu=lpawu_new
     pawtab(itypat)%lexexch=lexexch_new
     pawtab(itypat)%l_size=l_size_new
     pawtab(itypat)%lmn_size=lmn_size_new
     pawtab(itypat)%lmn2_size=lmn2_size_new
     pawtab(itypat)%mesh_size=mesh_size_new
     pawtab(itypat)%mqgrid=mqgrid_vl
     pawtab(itypat)%mqgrid_shp=mqgrid_shp_new
     pawtab(itypat)%usetvale=0;if (pspheads(itypat)%pawheader%pawver>=4) pawtab(itypat)%usetvale=1
     pawtab(itypat)%has_kij  =0;if (need_kij_new)   pawtab(itypat)%has_kij=1
     pawtab(itypat)%has_nabla=0;if (need_nabla_new) pawtab(itypat)%has_nabla=1
     pawtab(itypat)%has_vhntzc=0; if (need_vhntzc_new) pawtab(itypat)%has_vhntzc=1
     pawtab(itypat)%has_vhnzc=0; if (need_vhntzc_new) pawtab(itypat)%has_vhnzc=1

   end do ! itypat

!  Reallocate arrays depending on angl_size and ylm_size
   if (test_angl_size.or.test_ylm_size) then
     if(idtset/=1) then
       ABI_DEALLOCATE(pawang%ylmr)
       if (associated(pawang%ylmrgr)) ABI_DEALLOCATE(pawang%ylmrgr)
     end if
     ABI_ALLOCATE(pawang%ylmr,(ylm_size_new,angl_size_new))
     if (dtset%xclevel==2) then
       ABI_ALLOCATE(pawang%ylmrgr,(3,ylm_size_new,angl_size_new))
     else
       nullify(pawang%ylmrgr)
     end if
   end if

!  Reallocate arrays depending on nsym, l_size_max and l_max
   if (test_nsym.or.test_l_size_max.or.test_l_max) then
     if(idtset/=1)  then
       ABI_DEALLOCATE(pawang%zarot)
     end if
     ABI_ALLOCATE(pawang%zarot,(l_size_max_new,l_size_max_new,l_max_new,nsym_new))
   end if

!  Reallocate arrays depending on l_max, l_size_max and gnt_option
   if (test_l_max.or.test_l_size_max.or.test_gnt_option) then
     if(idtset/=1)  then
       ABI_DEALLOCATE(pawang%gntselect)
     end if
     if(gnt_option_new/=2) then
       ABI_ALLOCATE(pawang%gntselect,((l_size_max_new)**2,(l_max_new**2)*(l_max_new**2+1)/2))
     else
       ABI_ALLOCATE(pawang%gntselect,((2*l_size_max_new-1)**2,((2*l_max_new-1)**2)*((2*l_max_new-1)**2+1)/2))
     end if
     if (idtset==1) nullify(pawang%realgnt)
     if (associated(pawang%realgnt))  then
       ABI_DEALLOCATE(pawang%realgnt)
       ierr = ABI_ALLOC_STAT
     end if
   end if

!  Reallocate arrays depending on l_max and pawspnorb
   if (test_l_max.or.test_spnorb) then
     if(idtset/=1.and.pawang%use_ls_ylm>0)  then
       ABI_DEALLOCATE(pawang%ls_ylm)
     end if
     if (pawspnorb_new>0)  then
       ABI_ALLOCATE(pawang%ls_ylm,(2,l_max_new**2*(l_max_new**2+1)/2,2))
     end if
   end if

!  Reallocate arrays depending on angl_size
   if (test_angl_size) then
     if(idtset/=1) then
       ABI_DEALLOCATE(pawang%angwgth)
       if (associated(pawang%anginit))  then
         ABI_DEALLOCATE(pawang%anginit)
       end if
     end if
     ABI_ALLOCATE(pawang%angwgth,(angl_size_new))
     if (dtset%xclevel==2) then
       ABI_ALLOCATE(pawang%anginit,(3,angl_size_new))
     else
       nullify(pawang%anginit)
     end if
   end if

   pawang%angl_size=angl_size_new
   pawang%l_max=l_max_new
   pawang%l_size_max=l_size_max_new
   pawang%ylm_size=ylm_size_new
   pawang%nsym=nsym_new
   pawang%gnt_option=gnt_option_new
   pawang%use_ls_ylm=pawspnorb_new
   pawang%ngnt=ngnt_new
   pawang%ntheta=ntheta_new
   pawang%nphi=nphi_new

!  =====================================================================
!  ======================== DEALLOCATIONS ==============================
!  =====================================================================
 else if (option>=2) then

   do itypat=1,dtset%ntypat

     ABI_DEALLOCATE(pawrad(itypat)%rad)
     ABI_DEALLOCATE(pawrad(itypat)%radfact)
     ABI_DEALLOCATE(pawrad(itypat)%simfact)
     ABI_DEALLOCATE(pawtab(itypat)%gnorm)
     ABI_DEALLOCATE(pawtab(itypat)%indklmn)
     if(pawtab(itypat)%lpawu>=0.or.pawtab(itypat)%lexexch>=0)  then
       ABI_DEALLOCATE(pawtab(itypat)%klmntomn)
     end if
     if(pawtab(itypat)%lpawu>=0)  then
       ABI_DEALLOCATE(pawtab(itypat)%vee)
     end if
     if(pawtab(itypat)%lexexch>=0)  then
       ABI_DEALLOCATE(pawtab(itypat)%vex)
     end if
     if(pawtab(itypat)%lexexch>=0)  then
       ABI_DEALLOCATE(pawtab(itypat)%fk)
     end if
     ABI_DEALLOCATE(pawtab(itypat)%shapefunc)
     ABI_DEALLOCATE(pawtab(itypat)%shape_alpha)
     ABI_DEALLOCATE(pawtab(itypat)%shape_q)
     ABI_DEALLOCATE(pawtab(itypat)%tphi)
     ABI_DEALLOCATE(pawtab(itypat)%phi)
     ABI_DEALLOCATE(pawtab(itypat)%tphitphj)
     ABI_DEALLOCATE(pawtab(itypat)%phiphj)
     ABI_DEALLOCATE(pawtab(itypat)%coredens)
     ABI_DEALLOCATE(pawtab(itypat)%tcoredens)
     ABI_DEALLOCATE(pawtab(itypat)%tcorespl)
     if (pawtab(itypat)%usetvale>0)  then
       ABI_DEALLOCATE(pawtab(itypat)%tvalespl)
     end if
     pawtab(itypat)%usetvale=0
     if (pawtab(itypat)%has_kij  >0)  then
       ABI_DEALLOCATE(pawtab(itypat)%kij)
     end if
     if (pawtab(itypat)%has_nabla>0)  then
       ABI_DEALLOCATE(pawtab(itypat)%nabla_ij)
     end if
     if (pawtab(itypat)%has_vhntzc > 0)  then
       ABI_DEALLOCATE(pawtab(itypat)%VHntZC)
     end if
     if (pawtab(itypat)%has_vhnzc > 0)  then
       ABI_DEALLOCATE(pawtab(itypat)%VHnZC)
     end if
     pawtab(itypat)%has_kij=0;pawtab(itypat)%has_nabla=0
     pawtab(itypat)%has_vhntzc = 0
     pawtab(itypat)%has_vhnzc = 0
     ABI_DEALLOCATE(pawtab(itypat)%qijl)
     ABI_DEALLOCATE(pawtab(itypat)%eijkl)
     ABI_DEALLOCATE(pawtab(itypat)%dij0)
     ABI_DEALLOCATE(pawtab(itypat)%dltij)
     ABI_DEALLOCATE(pawtab(itypat)%rhoij0)
     ABI_DEALLOCATE(pawtab(itypat)%sij)
     if (associated(pawtab(itypat)%kmix))  then
       ABI_DEALLOCATE(pawtab(itypat)%kmix)
       ierr = ABI_ALLOC_STAT
     end if
     if (associated(pawtab(itypat)%dshpfunc))  then
       ABI_DEALLOCATE(pawtab(itypat)%dshpfunc)
       ierr = ABI_ALLOC_STAT
     end if
     if (associated(pawtab(itypat)%rad_for_spline))  then
       ABI_DEALLOCATE(pawtab(itypat)%rad_for_spline)
       ierr = ABI_ALLOC_STAT
     end if
     if (pawtab(itypat)%mqgrid_shp>0) then
       pawtab(itypat)%mqgrid_shp=0
       ABI_DEALLOCATE(pawtab(itypat)%qgrid_shp)
       ABI_DEALLOCATE(pawtab(itypat)%shapefncg)
     end if
   end do

   if (option>=3) then
     call destroy_pawang(pawang)
   end if

 end if

 DBG_EXIT("COLL")

end subroutine pawalloc

!!***
