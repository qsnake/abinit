!{\src2tex{textfont=tt}}
!!****f* ABINIT/surf
!! NAME
!! surf
!!
!! FUNCTION
!! Determination of the Bader surface.
!! Use rsurf to determine radius for one direction
!! simple bisection method is used
!! the bassin is tested following the gradient (follow) =
!! = the most time consuming
!! follow stops if the gradient line is near the atom
!! or if it is under already known part of surface - this is why
!! the surface is not computed row by row.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  This routine works primarily on the data contained in the defs_aimprom module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      drvaim
!!
!! CHILDREN
!!      coeffs_gausslegint,rsurf,timein
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine surf(aim_dtset)

 use m_profiling

 use defs_basis
 use defs_aimprom
 use defs_parameters
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'surf'
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
 use interfaces_63_bader, except_this_one => surf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: ii,ijj,iph,iph2,ith,ith2,jj,kk,mm,nn,nph,npmax,nth
 real(dp) :: ct1,ct2,phi,rr,rsmax,rsmin,rthe,rthe0,t1,t2,theta,tt0,vcth,vph,vth
 real(dp) :: wall,xy,xyz
 logical :: srch,stemp
!arrays
 real(dp) :: grho(3),vr(3),vv(3)

!************************************************************************
 ttsrf=zero

 rewind(unts)

 nth=aim_dtset%nth
 nph=aim_dtset%nph

!Coefficients for spherical Gauss quadrature

 ct1=cos(aim_dtset%themin)
 ct2=cos(aim_dtset%themax)
 call coeffs_gausslegint(ct1,ct2,cth,wcth,nth)
 call coeffs_gausslegint(aim_dtset%phimin,aim_dtset%phimax,ph,wph,nph)
 do ijj=1,nth
   th(ijj)=acos(cth(ijj))
   if (aim_dtset%isurf/=-1) then
     do jj=1,nph
       rs(ijj,jj)=zero
     end do
   end if
 end do

 npmax=aim_npmaxin
 rsmax=0.0
 rsmin=100.0
 rthe0=r0
 srch=.false.

 do ijj=1,3
   vv(ijj)=xatm(ijj,aim_dtset%batom)
 end do


 write(unto,*)
 write(unto,*) "BADER SURFACE DETERMINATION"
 write(unto,*) "==========================="
 write(unto,*)

 write(untout,*)
 write(untout,*) "BADER SURFACE DETERMINATION"
 write(untout,*) "==========================="
 write(untout,*)

 write(unto,'(" Atom:  ",i3,3F15.10)') aim_dtset%batom,vv
 write(unto,'(" Theta: ",i3,2F15.10)') nth,aim_dtset%themin,aim_dtset%themax
 write(unto,'(" Phi:   ",i3,2F15.10)') nph,aim_dtset%phimin,aim_dtset%phimax

 write(untout,'(" Atom:  ",i3,3F15.10)') aim_dtset%batom,vv
 write(untout,'(" Theta: ",i3,2F15.10)') nth,aim_dtset%themin,aim_dtset%themax
 write(untout,'(" Phi:   ",i3,2F15.10)') nph,aim_dtset%phimin,aim_dtset%phimax

 write(unts,'(i3,3F15.10)') aim_dtset%batom,vv
 write(unts,'(i3,2F15.10)') nth,aim_dtset%themin,aim_dtset%themax
 write(unts,'(i3,2F15.10)') nph,aim_dtset%phimin,aim_dtset%phimax

!write(unto,*) 'npmax in surf= ',npmax

 ith=0
 iph=0
 tt0=0._dp
 call timein(tt0,wall)

 write(untout,*)
 write(untout,*) "DEVELOPMENT OF THE RADII DETERMINATIONS"
 write(untout,*) "========================================"
 write(untout,*)
 write(untout,*) "Determination near the CPs:"

!Determination of the CP neighbouring radii

 if (aim_dtset%isurf/=-1) then
   srch=.true.

   do ijj=1,nbcp
!    if ((icpc(ijj) == -1)) then

     rthe0=vnorm(pc(:,ijj),0)
     do jj=1,3
       vr(jj)=pc(jj,ijj)-vv(jj)+xatm(jj,aim_dtset%batom)
     end do
     xy=vr(1)*vr(1)+vr(2)*vr(2)
     xyz=xy+vr(3)*vr(3)
     xyz=sqrt(xyz)

     if (xy < aim_xymin) then
       vcth=1._dp
       if (vr(3) < 0._dp) vcth=-vcth
       vph=0._dp
     else
       vcth=vr(3)/xyz
       vph=atan2(vr(2),vr(1))
     end if

     vth=acos(vcth)
     write(untout,'(/," BCP: (index,theta,phi)",I4,2E16.8)') ijj,vth,vph

     if (vth < th(1)) then
       ith=0
     else
       if (vth > th(nth)) then
         ith=nth
       else
         do ii=2,nth
           if (vth < th(ii)) then
             ith=ii-1
             exit
           end if
         end do
       end if
     end if

     if (vph < ph(1)) then
       iph=0
     else
       if (vph > ph(nph)) then
         iph=nph
       else
         do ii=2,nph
           if (vph < ph(ii)) then
             iph=ii-1
             exit
           end if
         end do
       end if
     end if

     write(untout,*) "ATOMIC RADII (ith,iph,theta,phi,radius)"
     do jj=-1,2
       do kk=-1,2
         ith2=ith+jj
         iph2=iph+kk
         stemp=(iph2 > 0).and.(iph2 < nph+1)
         stemp=(stemp.and.((ith2 > 0).and.(ith2 < nth+1)))
         if (stemp) then
           if (abs(rs(ith2,iph2))<1.0d-12) then
             rthe=rthe0
             theta=th(ith2)
             phi=ph(iph2)
             if (deb) write(unto,*) ':CALCULATING NP',theta,phi,rthe,npmax
             if (deb) write(unto,*) ':CALCULATING NP',theta,phi,rthe,npmax
             call timein(t1,wall)
             call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
             call timein(t2,wall)
             t2=t2-t1
             rs(ith2,iph2)=rr
!            write(unts,'(2F12.8,2E16.8)') theta,phi,rr,wcth(ijj)*wph(jj)
             write(unto,'(":RSUR PC ",3i3,4E16.8,F10.4)') ijj,jj,kk,theta,phi,rr,wcth(ith2)*wph(iph2),t2
             write(untout,'(a,2i3,3E16.8)') '-  ',jj,kk,theta,phi,rr
             rthe0=rr
           end if
         end if

       end do ! kk
     end do ! jj

!    end if

   end do ! ijj (loop on BCP)

!  DEBUG
!  write(std_out,*)' surf : near BCP '
!  do ijj=1,nth
!  do jj=1,nph
!  write(std_out,*)ijj,jj,rs(ijj,jj)
!  end do
!  end do
!  ENDDEBUG


   srch=.true.
   do ijj=nbcp+1,nbcp+nrcp     ! Loop on RCP
!    if ((icpc(ijj) == 1)) then
     rthe0=max(rminl(aim_dtset%batom),r0)
     do jj=1,3
       vr(jj)=pc(jj,ijj)-vv(jj)+xatm(jj,aim_dtset%batom)
     end do
     xy=vr(1)*vr(1)+vr(2)*vr(2)
     xyz=xy+vr(3)*vr(3)
     xyz=sqrt(xyz)

     if (xy < aim_xymin) then
       vcth=1._dp
       if (vr(3) < 0._dp) vcth=-vcth
       vph=0._dp
     else
       vcth=vr(3)/xyz
       vph=atan2(vr(2),vr(1))
     end if
     vth=acos(vcth)
     write(untout,'(/,";RCP: (index,theta,phi)",I4,2E16.8)') ijj-nbcp,vth,vph

     if (vth < th(1)) then
       ith=0
     else
       if (vth > th(nth)) then
         ith=nth
       else
         do ii=2,nth
           if (vth < th(ii)) then
             ith=ii-1
             exit
           end if
         end do
       end if
     end if

     if (vph < ph(1)) then
       iph=0
     else
       if (vph > ph(nph)) then
         iph=nph
       else
         do ii=2,nph
           if (vph < ph(ii)) then
             iph=ii-1
             exit
           end if
         end do
       end if
     end if

     write(untout,*) "ATOMIC RADIUS (ith,iph,theta,phi,radius)"
     do jj=-1,2
       do kk=-1,2
         ith2=ith+jj
         iph2=iph+kk
         stemp=(iph2 > 0).and.(iph2 < nph+1)
         stemp=stemp.and.(ith2 > 0).and.(ith2 < nth+1)

         if (stemp) then
           if ((abs(rs(ith2,iph2))<1.0d-12)) then
             rthe=rthe0
             theta=th(ith2)
             phi=ph(iph2)
             if (deb) write(unto,*) ':CALCULATING NP',theta,phi,rthe,npmax
             call timein(t1,wall)
             call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
             call timein(t2,wall)
             t2=t2-t1
             rs(ith2,iph2)=rr
!            write(unts,'(2F12.8,2E16.8)') theta,phi,rr,wcth(ijj)*wph(jj)
             write(unto,'(":RSUR PC ",3i3,4E16.8,F10.4)') ijj,jj,kk,theta,phi,rr,wcth(ith2)*wph(iph2),t2
             write(untout,'(a,2i3,3E16.8)') '-  ',jj,kk,theta,phi,rr
             rthe0=rr
           end if
         end if

       end do ! kk
     end do ! jj
!    end if

   end do ! ijj (Loop on RCP)

!  DEBUG
!  write(std_out,*)' surf : near RCP '
!  do ijj=1,nth
!  do jj=1,nph
!  write(std_out,*)ijj,jj,rs(ijj,jj)
!  end do
!  end do
!  ENDDEBUG

!  Boundary angles
   rthe0=r0
   srch=.true.
   write(untout,*)
   write(untout,*) "The boundary angles:"
   write(untout,*) "===================="
   write(untout,*) "ATOMIC RADIUS (ith,iph,theta,phi,radius)"

!  Must have sufficient angular sampling
   if ((nth > 8).and.(nph > 8)) then
     rthe=r0
     do ijj=1,2
       theta=th(ijj)
       if (ijj==2) rthe=rs(1,1)
       do jj=1,nph
         phi=ph(jj)
         call timein(t1,wall)
         if (abs(rs(ijj,jj))<1.0d-12) then
           if (deb) write(unto,*) ':CALC NP',theta,phi,rthe,npmax
           call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
           rs(ijj,jj)=rr
         end if
         call timein(t2,wall)
         t2=t2-t1
         write(unto,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rs(ijj,jj),wcth(ijj)*wph(jj),t2
         write(untout,'(a,2i3,3E16.8)') '-  ',ijj,jj,theta,phi,rr
         if (rr < rsmin) rsmin=rs(ijj,jj)
         if (rr > rsmax) rsmax=rs(ijj,jj)
         rthe=rs(ijj,jj)
       end do ! jj
     end do ! ijj

     write(untout,*)

     rthe=rs(2,1)
     do jj=1,2
       phi=ph(jj)
       if (jj==2) rthe=rs(2,2)
       do ijj=3,nth
         theta=th(ijj)
         t2=0.0
         call timein(t1,wall)
         if (abs(rs(ijj,jj))<1.0d-12) then
           if (deb) write(unto,*) ':CALC NP',theta,phi,rthe,npmax
           call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
           rs(ijj,jj)=rr
         end if
         call timein(t2,wall)
         t2=t2-t1
         write(unto,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rs(ijj,jj),wcth(ijj)*wph(jj),t2
         write(untout,'(2i3,3E16.8)') ijj,jj,theta,phi,rr
         if (rr < rsmin) rsmin=rs(ijj,jj)
         if (rr > rsmax) rsmax=rs(ijj,jj)
         rthe=rs(ijj,jj)
       end do ! ijj
     end do ! jj

     write(untout,*)

     rthe=rs(nth-1,2)
     do ijj=nth-1,nth
       theta=th(ijj)
       if (ijj==nth) rthe=rs(nth,2)
       do jj=3,nph
         phi=ph(jj)
         call timein(t1,wall)
         if (abs(rs(ijj,jj))<1.0d-12) then
           if (deb) write(unto,*) ':CALC NP',theta,phi,rthe,npmax
           call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
           rs(ijj,jj)=rr
         end if
         call timein(t2,wall)
         t2=t2-t1
         write(unto,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rs(ijj,jj),wcth(ijj)*wph(jj),t2
         write(untout,'(2i3,3E16.8)') ijj,jj,theta,phi,rr
         if (rr < rsmin) rsmin=rs(ijj,jj)
         if (rr > rsmax) rsmax=rs(ijj,jj)
         rthe=rs(ijj,jj)
       end do ! jj
     end do ! ijj

     rthe=rs(2,nph-1)
     do jj=nph-1,nph
       phi=ph(jj)
       if (jj==nph) rthe=rs(2,nph)
       do ijj=3,nth-2
         theta=th(ijj)
         t2=0.0
         call timein(t1,wall)
         if (abs(rs(ijj,jj))<1.0d-12) then
           if (deb) write(unto,*) ':CALC NP',theta,phi,rthe,npmax
           call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
           rs(ijj,jj)=rr
         end if
         call timein(t2,wall)
         t2=t2-t1
         write(unto,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rs(ijj,jj),wcth(ijj)*wph(jj),t2
         write(untout,'(2i3,3E16.8)') ijj,jj,theta,phi,rr
         if (rr < rsmin) rsmin=rs(ijj,jj)
         if (rr > rsmax) rsmax=rs(ijj,jj)
         rthe=rs(ijj,jj)
       end do ! ijj
     end do ! jj
     write(untout,*)

!    Complementary bands for boundary angles
     nn=int(real(nth)/1.4d1)
     if (nn > 1) then
       do ii=1,nn-1
         mm=int(nth/nn)*ii
         do kk=0,1
           mm=mm+kk
           theta=th(mm)
           rthe=rs(mm,2)
           do jj=3,nph-2
             phi=ph(jj)
             call timein(t1,wall)
             if (abs(rs(mm,jj))<1.0d-12) then
               if (deb) write(unto,*) ':CALC NP',theta,phi,rthe,npmax
               call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
               rs(mm,jj)=rr
             end if
             call timein(t2,wall)
             t2=t2-t1
             write(unto,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rs(mm,jj),wcth(mm)*wph(jj),t2
             write(untout,'(2i3,3E16.8)') mm,jj,theta,phi,rr
             if (rr < rsmin) rsmin=rs(mm,jj)
             if (rr > rsmax) rsmax=rs(mm,jj)
             rthe=rs(mm,jj)
           end do ! jj
         end do ! kk
       end do ! ii
     end if ! nn>1

     write(untout,*)

     nn=nint(real(nph)/1.2d1)
     if (nn > 1) then
       do ii=1,nn-1
         mm=int(nph/nn)*ii
         do kk=0,1
           mm=mm+kk
           phi=ph(mm)
           rthe=rs(2,mm)

           do jj=3,nth-2
             theta=th(jj)
             call timein(t1,wall)
             if (abs(rs(jj,mm))<1.0d-12) then
               if (deb) write(unto,*) ':CALC NP',theta,phi,rthe,npmax
               call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
               rs(jj,mm)=rr
             end if
             call timein(t2,wall)
             t2=t2-t1
             write(unto,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rs(jj,mm),wcth(jj)*wph(mm),t2
             write(untout,'(2i3,3E16.8)') jj,mm,theta,phi,rr
             if (rr < rsmin) rsmin=rs(jj,mm)
             if (rr > rsmax) rsmax=rs(jj,mm)
             rthe=rs(jj,mm)
           end do ! jj

         end do ! kk
       end do ! ii
     end if  ! nn>1

   end if ! sufficient sampling to determine boundary angles

   write(untout,*)

!  DEBUG
!  write(std_out,*)' surf : after boundary angles '
!  do ijj=1,nth
!  do jj=1,nph
!  write(std_out,*)ijj,jj,rs(ijj,jj)
!  end do
!  end do
!  ENDDEBUG

!  Output the complete Bader surface

   write(untout,*) "The complete Bader surface:"
   write(untout,*) "==========================="
   write(untout,*) "ATOMIC RADIUS (ith,iph,theta,phi,radius)"
   rthe0=r0
   srch=.true.

!  Systematic scanning of the grid and the determination of the missing radii
   do ijj=1,nth
     theta=th(ijj)
     rthe=rthe0
     do jj=1,nph
       phi=ph(jj)
       t2=0.0
       call timein(t1,wall)
       if (abs(rs(ijj,jj))<1.0d-12) then
         if (deb) write(unto,*) ':CALCULATING NP',ijj,jj,theta,phi,rthe,npmax,rs(ijj,jj)
         call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
         rs(ijj,jj)=rr
       end if
       call timein(t2,wall)
       t2=t2-t1
       write(unts,'(2F12.8,2E16.8)') theta,phi,rs(ijj,jj),wcth(ijj)*wph(jj)
       write(unto,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rs(ijj,jj),wcth(ijj)*wph(jj),t2
       write(untout,'(a,2i3,3E16.8)') '   ',ijj,jj,theta,phi,rs(ijj,jj)
       if (jj == 1) rthe0=rs(ijj,jj)
       if (rs(ijj,jj) < rsmin) rsmin=rs(ijj,jj)
       if (rs(ijj,jj) > rsmax) rsmax=rs(ijj,jj)
       rthe=rs(ijj,jj)
     end do ! jj
   end do ! ijj
   write(unts,'(2F15.10)') rsmin,rsmax
   write(untout,'(/," The minimal and maximal radii:",/,/,"     ",2F15.10)') rsmin,rsmax

!  DEBUG
!  write(std_out,*)' surf : final output '
!  do ijj=1,nth
!  do jj=1,nph
!  write(std_out,*)ijj,jj,rs(ijj,jj)
!  end do
!  end do
!  ENDDEBUG

 end if ! determination of the critical surface

 call timein(ttsrf,wall)
 ttsrf=ttsrf-tt0

end subroutine surf
!!***
