!{\src2tex{textfont=tt}}
!!****f* ABINIT/psps_init_from_dtset
!! NAME
!! psps_init_from_dtset
!!
!! FUNCTION
!! Allocate and initialise all part of psps structure that are dependent
!! of a given dataset.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dtset=<type dataset_type>a given dataset
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! psps=<type pseudopotential_type>the pseudopotentials description
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      getdim_nloc,matr3inv,setmqgrid
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psps_init_from_dtset(dtset, idtset, psps, pspheads)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psps_init_from_dtset'
 use interfaces_32_util
 use interfaces_56_recipspace
 use interfaces_57_iovars
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idtset
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 type(pspheader_type),intent(in) :: pspheads(psps%npsp)

!Local variables-------------------------------
!scalars
 integer,save :: dimekb_old=-1,lmnmax_old=-1,lnmax_old=-1,mqgridff_old=0
 integer,save :: mqgridvl_old=0,ntypat_old=-1,usepaw_old=-1
 integer :: ipsp,lmnmax,lmnmaxso,lnmax,lnmaxso,newmqgrid,newmqgriddg,nptsgvec
 real(dp) :: gprimd_orig(3,3)

! *************************************************************************

 psps%optnlxccc   = dtset%optnlxccc
!Determine the number of points needed in reciprocal space to represent the
!pseudopotentials (either set by hand from input variable or set automatically
!by abinit)
 nptsgvec         = 200 !This has to be chosen one and for all or else ??
 newmqgrid        = dtset%mqgrid
 newmqgriddg      = dtset%mqgriddg
 call matr3inv(dtset%rprimd_orig,gprimd_orig)
 call setmqgrid(newmqgrid,newmqgriddg,dtset%ecut*dtset%dilatmx**2,&
& dtset%pawecutdg*dtset%dilatmx**2,gprimd_orig,nptsgvec,psps%usepaw)
 psps%mqgrid_ff   = newmqgrid
 if (psps%usepaw == 1) then
   psps%mqgrid_vl = newmqgriddg
 else
   psps%mqgrid_vl = newmqgrid
 end if

!Determine the maximum number of projectors, for the set of pseudo atom
 call getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,dtset%mixalch,psps%npsp,dtset%npspalch,&
& dtset%ntypat,dtset%ntypalch,pspheads)

 psps%npspalch = dtset%npspalch
 psps%ntypat   = dtset%ntypat
 psps%ntypalch = dtset%ntypalch
 psps%ntyppure = dtset%ntyppure

!Set the flag for reciprocal space or real space calculations
 psps%vlspl_recipSpace = (dtset%icoulomb /= 1)
!changed by RShaltaf
 psps%positron = dtset%positron
 psps%useylm   = dtset%useylm

 if (idtset > 1) then
   ABI_DEALLOCATE(psps%algalch)
   ABI_DEALLOCATE(psps%mixalch)
 end if
 ABI_ALLOCATE(psps%algalch,(psps%ntypalch))
 ABI_ALLOCATE(psps%mixalch,(psps%npspalch,psps%ntypalch))
 psps%algalch(1:psps%ntypalch)=dtset%algalch(1:psps%ntypalch)
 psps%mixalch(1:psps%npspalch,1:psps%ntypalch)=dtset%mixalch(1:psps%npspalch,1:psps%ntypalch)

!Set mpspso and psps%pspso
!Warning : mpspso might be different for each dataset.
 psps%mpspso=1
 do ipsp=1,dtset%npsp
   if(dtset%nspinor==1)then
     psps%pspso(ipsp)=0
   else
     if(dtset%so_psp(ipsp)/=1)then
       psps%pspso(ipsp)=dtset%so_psp(ipsp)
     else
       psps%pspso(ipsp)=pspheads(ipsp)%pspso
     end if
     if(psps%pspso(ipsp)/=0)psps%mpspso=2
   end if
!  Ideally the following line should not exist, but at present, the space has to be booked
   if(pspheads(ipsp)%pspso/=0)psps%mpspso=2
 end do

!Set mpssoang, lmnmax, lnmax
 if(psps%mpspso==1)then
   psps%mpssoang=psps%mpsang
   psps%lmnmax  =lmnmax
   psps%lnmax   =lnmax
 else
   psps%mpssoang=2*psps%mpsang-1
   psps%lmnmax=lmnmaxso
   psps%lnmax=lnmaxso
 end if
 if (psps%useylm==0) then
   psps%lmnmax=psps%lnmax
 end if

!Set dimekb
 if (psps%usepaw==0) then
   psps%dimekb=psps%lnmax
 else
   psps%dimekb=psps%lmnmax*(psps%lmnmax+1)/2
 end if

!The following arrays are often not deallocated before the end of the dtset loop
!and might keep their content from one dataset to the other,
!if the conditions are fulfilled
 if(dimekb_old/=psps%dimekb .or. ntypat_old/=dtset%ntypat .or. usepaw_old/=psps%usepaw) then
   if(idtset/=1)ABI_DEALLOCATE(psps%ekb)
   ABI_ALLOCATE(psps%ekb,(psps%dimekb,dtset%ntypat*(1-psps%usepaw)))
   dimekb_old=psps%dimekb
 end if
 if(lmnmax_old/=psps%lmnmax .or. ntypat_old/=dtset%ntypat)then
   if(idtset/=1)ABI_DEALLOCATE(psps%indlmn)
   ABI_ALLOCATE(psps%indlmn,(6,psps%lmnmax,dtset%ntypat))
   lmnmax_old=psps%lmnmax
 end if
 if(mqgridff_old/=psps%mqgrid_ff .or. lnmax_old/=psps%lnmax .or. ntypat_old/=dtset%ntypat)then
   if(idtset/=1) then
     ABI_DEALLOCATE(psps%ffspl)
     ABI_DEALLOCATE(psps%qgrid_ff)
   end if
   ABI_ALLOCATE(psps%ffspl,(psps%mqgrid_ff,2,psps%lnmax,dtset%ntypat))
   ABI_ALLOCATE(psps%qgrid_ff,(psps%mqgrid_ff))
   mqgridff_old=psps%mqgrid_ff
   lnmax_old=psps%lnmax
 end if
 if(mqgridvl_old/=psps%mqgrid_vl .or. ntypat_old/=dtset%ntypat)then
   if(idtset/=1) then
     ABI_DEALLOCATE(psps%vlspl)
     ABI_DEALLOCATE(psps%qgrid_vl)
   end if
   if (idtset/=1 .and. .not.psps%vlspl_recipSpace) then
     ABI_DEALLOCATE(psps%dvlspl)
   end if
   ABI_ALLOCATE(psps%vlspl,(psps%mqgrid_vl,2,dtset%ntypat))
   ABI_ALLOCATE(psps%qgrid_vl,(psps%mqgrid_vl))
   if (.not.psps%vlspl_recipSpace) then
     ABI_ALLOCATE(psps%dvlspl,(psps%mqgrid_vl,2,dtset%ntypat))
   end if
   mqgridvl_old=psps%mqgrid_vl
 end if
 if(ntypat_old/=dtset%ntypat.or. usepaw_old/=psps%usepaw)then
   if(idtset/=1)ABI_DEALLOCATE(psps%xccc1d)
   ABI_ALLOCATE(psps%xccc1d,(psps%n1xccc*(1-psps%usepaw),6,dtset%ntypat))
   usepaw_old=psps%usepaw
 end if
 if(ntypat_old/=dtset%ntypat)then
   if(idtset/=1) then
     ABI_DEALLOCATE(psps%xcccrc)
     ABI_DEALLOCATE(psps%ziontypat)
     ABI_DEALLOCATE(psps%znucltypat)
   end if
   ABI_ALLOCATE(psps%xcccrc,(dtset%ntypat))
   ABI_ALLOCATE(psps%znucltypat,(dtset%ntypat))
   ABI_ALLOCATE(psps%ziontypat,(dtset%ntypat))
   ntypat_old=dtset%ntypat
 end if
 psps%ziontypat(:)=dtset%ziontypat(:)
 end subroutine psps_init_from_dtset
!!***
