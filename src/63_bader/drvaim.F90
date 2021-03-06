!{\src2tex{textfont=tt}}
!!****f* ABINIT/drvaim
!! NAME
!! drvaim
!!
!! FUNCTION
!! Main driver for the Bader analysis
!! it looks the values of the input variables
!! and calls corresponding procedures
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! aim_dtset = the structured entity containing all input variables
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  this routine acts primarily on the data contained in the aimprom module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!      addout,aim_follow,cpdrv,critics,graph,initaim,integrho,integvol,plint
!!      rsurf,surf,timein
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine drvaim(aim_dtset)

 use m_profiling

 use defs_basis
 use defs_parameters
 use defs_aimfields
 use defs_aimprom
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'drvaim'
 use interfaces_18_timing
 use interfaces_63_bader, except_this_one => drvaim
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: iat,iatinit,inxat,ipos,iposinit
 integer :: npmax,nstep
 real(dp) :: dstlim,rr,ss,t1,t2,tf,wall
 logical :: debold,sfour,srch,sthree,stwo
!arrays
 real(dp) :: grho(3),xstart(3)

! *********************************************************************

!These input variables might be modified during what follows,
!so, they are copied outside of aim_dtset.
 inxat=aim_dtset%batom
 r0=aim_dtset%atrad
 h0=aim_dtset%folstp
 maxatdst=aim_dtset%maxatd
 maxcpdst=aim_dtset%maxcpd

 dstlim=maxcpdst

!Flags from the old version
!to be remove later
 deb=.false.
 stwo=.true.
 sthree=.true.
 sfour=.false.
 srch=.false.

 npmax=aim_npmaxin

!Main initialisation procedure -
!- it reads ABINIT density file and files
!with core densities and initialises the fields for
!spline interpolation

 call initaim(aim_dtset)


!CP SEARCHING

 if (aim_dtset%crit /= 0) then

   if (aim_dtset%crit==3) then
!    old version of the driver for searching CPs (original code)
     call critics(aim_dtset,inxat,stwo,sthree,sfour,dstlim)
   else
!    driver for searching CPs with Popellier algorithm
     call cpdrv(aim_dtset)
   end if

 end if

!
!BADER SURFACE CALCUL
!

 if (aim_dtset%isurf==1) then
!  driver for determination of the Bader surface
   call surf(aim_dtset)
 end if

!
!CHARGE INTEGRATIOM
!

 if (aim_dtset%irho==1) then
   call integrho(aim_dtset)
 end if

!
!VOLUME INTEGRATION OF THE BADER ATOM
!

 if (aim_dtset%ivol==1) then
   call integvol()
 end if

!
!ONE RADIUS OF THE BADER SURFACE
!

 if (aim_dtset%irsur==1) then
   if (aim_dtset%isurf/=0) srch=.true.
   iat=aim_dtset%batom
   ss=r0
   call timein(t1,wall)
   call rsurf(aim_dtset,rr,grho,aim_dtset%th0,aim_dtset%phi0,ss,iat,npmax,srch)
   call timein(t2,wall)
   t2=t2-t1
   write(unts,'(2F12.8,F15.10)') aim_dtset%th0,aim_dtset%phi0,rr
   write(std_out,'(":RSUR ",2F12.8,2F15.10)') aim_dtset%th0,aim_dtset%phi0,rr,t2
 end if

!
!FOLLOW THE GRADIENT PATH FROM ONE POINT
!

 if (aim_dtset%foll==1) then
   iatinit=aim_dtset%batom
   iposinit=batcell
   if (aim_dtset%isurf/=0) srch=.true.
   debold=deb
   xstart(:)=aim_dtset%foldep(:)
   call timein(t1,wall)
   call aim_follow(aim_dtset,xstart,npmax,srch,iatinit,iposinit,iat,ipos,nstep)
   call timein(t2,wall)
   tf=t2-t1
   write(std_out,'(":TIME in aim_follow:", F12.4)') tf
 end if

 if (aim_dtset%plden == 1) then
!  profile of the density integrated in plane xy
!  belong the z-axes - not finished - cut3d better !
   call plint()
 end if

 if ((aim_dtset%denout > 0).or.(aim_dtset%lapout > 0)) then
!  additional outputs of density and laplacian fields
!  in the plane or line
   call addout(aim_dtset)
 end if

 if (aim_dtset%gpsurf == 1) then
!  script for gnuplot - simple demonstration of the
!  computed surface
   call graph(unts,untg)
 end if

!Deallocation of global variables allocated in initaim
!and declared in defs_aimfields.
 ABI_DEALLOCATE(dig1)
 ABI_DEALLOCATE(dig2)
 ABI_DEALLOCATE(dig3)
 ABI_DEALLOCATE(llg1)
 ABI_DEALLOCATE(llg2)
 ABI_DEALLOCATE(llg3)
 ABI_DEALLOCATE(cdig1)
 ABI_DEALLOCATE(cdig2)
 ABI_DEALLOCATE(cdig3)
 ABI_DEALLOCATE(ddx)
 ABI_DEALLOCATE(ddy)
 ABI_DEALLOCATE(ddz)
 ABI_DEALLOCATE(rrad)
 ABI_DEALLOCATE(crho)
 ABI_DEALLOCATE(sp2)
 ABI_DEALLOCATE(sp3)
 ABI_DEALLOCATE(sp4)
 ABI_DEALLOCATE(corlim)
 ABI_DEALLOCATE(dvl)
 ABI_DEALLOCATE(ndat)
 ABI_DEALLOCATE(rminl)
!Deallocation of global variables allocated in initaim
!and declared in defs_aimprom.
 ABI_DEALLOCATE(typat)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(xatm)

end subroutine drvaim
!!***
