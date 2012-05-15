!{\src2tex{textfont=tt}}
!!****f* ABINIT/cpdrv
!! NAME
!! cpdrv
!!
!! FUNCTION
!! Critical points (CPs) searching driver
!! First Bond CPs are searched for each pair atom-its neighbor
!! (distance cutoff=maxatdst)
!! then Ring CPs for each pair of BCPs
!! and finally Cage CPs for each pair of RCPs.
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
!!
!! SIDE EFFECTS
!!  this routine treat information contained in the aim_prom module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! TODO
!! Should combine parts of code that are similar ...
!!
!! PARENTS
!!      drvaim
!!
!! CHILDREN
!!      aim_follow,critic,sort_dp,timein,vgh_rho
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cpdrv(aim_dtset)

 use m_profiling

 use defs_basis
 use defs_aimprom
 use defs_parameters
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cpdrv'
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
 use interfaces_63_bader, except_this_one => cpdrv
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: iat,iatinit,ii,inxat,inxcell,ipos,iposinit,ires,jj,kk,nb,nb_now
 integer :: nn,nstep,nvs
 real(dp) :: chg,diff1,diff2,diff3,dist,prj,rtdiff,ss,tt0,wall
 logical :: srch=.false.
!arrays
 integer :: ibat(nnpos*natom),inatm(nnpos*natom),incell(nnpos*natom)
 integer :: ipibat(nnpos*natom)
 integer,allocatable :: indexcp(:),nr(:)
 real(dp) :: bmin(natom),dif(3),dists(nnpos*natom),ev(3),evec(3,3),grho(3)
 real(dp) :: hrho(3,3),pom(3),rr(3),uu(3),vv(3),xorig(3)
 real(dp),allocatable :: sortguide(:)
!no_abirules
!Warning : bcp_type should be transformed to cp_type
 type(bcp_type),allocatable :: bcp(:),ccp(:),cp_tmp(:),rcp(:)

!************************************************************************
 inxat=aim_dtset%batom
 slc=-1
 rminl(:)=aim_dtset%rmin
 bmin(:)=0._dp
 ttcp=0._dp

 write(unto,*)
 write(unto,*) "CRITICAL POINTS ANALYSIS"
 write(unto,*) "========================"
 write(unto,*)

 write(untout,*)
 write(untout,*) "CRITICAL POINTS ANALYSIS"
 write(untout,*) "========================"
 write(untout,*)


 xorig(:)=xatm(:,inxat)

 call timein(tt0,wall)

!Searching the neighbouring atoms

 if (aim_dtset%crit > 0) then
   nvs=0
   do ii=1,nnpos
     do jj=1,natom
       dist=0._dp
       dif(:)=xatm(:,inxat)-xatm(:,jj)-atp(:,ii)
       dif(:)=dif(:)/aim_dtset%scal(:)
       dist=vnorm(dif,0)
       if (dist < tol6 ) then
         inxcell=ii
       elseif (dist < maxatdst) then
         nvs=nvs+1
         dists(nvs)=dist
         inatm(nvs)=jj
         incell(nvs)=ii
       end if
     end do
   end do

   write(unto,*) "ATOM:"
   write(unto,*) 'inxat :', inxat, 'inxcell :', inxcell
   write(unto, '(3es16.6)' ) (xorig(ii),ii=1,3)
   write(unto,*)

   write(untout,*) "ATOM:"
   write(untout,*) 'inxat :', inxat, 'inxcell :', inxcell
   write(untout, '(3es16.6)') (xorig(ii),ii=1,3)
   write(untout,*)

   ABI_ALLOCATE(nr,(nvs))
   do ii=1,nvs
     nr(ii)=ii
   end do

!  Ordering of the nearest neighbouring atoms
   call sort_dp(nvs,dists,nr,tol14)

   nb=0
   write(unto,*) "NEIGHBORING ATOMS (atindex,cellindex,distance(in bohr)):"
   write(untout,*) "NEIGHBORING ATOMS (atindex,cellindex,distance(in bohr)):"
   do ii=1,nvs
     nn=nr(ii)
     if (dists(ii) < maxatdst) then
       nb=nb+1
       ibat(nb)=inatm(nn)
       ipibat(nb)=incell(nn)
       write(unto,*) ':NEIG ',inatm(nn),incell(nn),dists(ii)
       write(untout,'("      ",2I6,F16.8)')inatm(nn),incell(nn),dists(ii)
     else
       exit
     end if
   end do

!  SEARCHING BCP
   ABI_ALLOCATE(bcp,(nb))
   nbcp=0
   iatinit=inxat
   iposinit=inxcell
   bcp(:)%iat=0
   bcp(:)%ipos=0

   write(unto,*)
   write(unto,*) "BONDING CRITICAL POINTS (BCP)"
   write(unto,*) "============================="
   write(unto,*)

   write(untout,*)
   write(untout,*) "BONDING CRITICAL POINTS (BCP)"
   write(untout,*) "============================="
   write(untout,*)

   srbcp: do ii=1,nb

!    Start the search for BCP from the midistance between the atom
!    and his neighbor.
     vv(:)=(xatm(:,inxat)+xatm(:,ibat(ii))+atp(:,ipibat(ii)))/2._dp

     call critic(aim_dtset,vv,ev,evec,aim_dmaxcs,ires,-1)

     if (ires==0) then
!      Testing if CP is already known
       if (nbcp > 0) then
         do jj=1,nbcp
           pom(:)=vv(:)-bcp(jj)%rr(:)-xorig(:)
           dist=vnorm(pom,0)
           if (dist < aim_dtset%dpclim) then
             write(unto,*) 'BCP already known  !'
             cycle srbcp
           end if
         end do
       end if
       rr(:)=vv(:)-xorig(:)
       ss=vnorm(rr,0)
       if (ss > maxcpdst) then
         write(unto, '(a,es16.6,a,es16.6)' ) 'BCP distance from atom,',ss,', exceed maxcpdst =',maxcpdst
         cycle srbcp
       end if
       nn=0
       do jj=1,3
         nn=nn+ev(jj)/abs(ev(jj))
       end do
       write(unto, '(a,3es16.6,i4)') ' vv(1:3), nn',(vv(jj), jj=1,3), nn
       write(unto, '(a,3es16.6)') 'ev: ', (ev(jj), jj=1,3)
       if (nn /= -1) then
         write(unto,*) ' The trial critical point is not a BCP !'
         cycle srbcp
       end if
       write(unto, '(a,3es16.6)' ) 'evec(:,1): ',(evec(jj,1), jj=1,3)
       pom(:)=evec(:,1)
       dist=vnorm(pom,0)
       prj=dot_product(evec(:,1),rr)
       write(unto,*) 'prj:', prj, vnorm(evec(:,1),0)
       dist=vnorm(evec(:,1),0)
       uu(:)=vv(:)-sign(aim_epstep,prj)*evec(:,1)/dist

!      Testing whether this BCP "is bonded" to the considered atom
       call aim_follow(aim_dtset,uu,aim_npmaxin,srch,iatinit,iposinit,iat,ipos,nstep)
!      write(unto,*) 'do', iat, ipos
!      if ((iat==0).or.(ipos==0)) cycle
!      write(unto,*) 'APOS: ',(xatm(jj,iat)+atp(jj,ipos), jj=1,3)
       if ((iat/=inxat).or.(inxcell/=ipos)) then
         write(unto,*) ' The trial BCP is not bonded to the Bader atom'
         cycle srbcp
       end if

!      A new BCP has been found !
       nbcp=nbcp+1

!      Searching for the second bonded atom
       ss=vnorm(rr,0)
       diff1=ss
       diff3=dists(ii)
       uu(:)=vv(:)+sign(aim_epstep,prj)*evec(:,1)/dist
       if ((abs(bmin(iat))<1.0d-12).or.( ss<bmin(iat))) then
         bmin(iat)=ss
       end if
       call aim_follow(aim_dtset,uu,aim_npmaxin,srch,iatinit,iposinit,iat,ipos,nstep)
       if ((iat==0).or.(ipos==0)) then
         write(unto,*) ' The trial BCP is not bonded to a bonding atom !'
!        cycle srbcp
       end if
       pom(:)=vv(:)-xatm(:,iat)-atp(:,ipos)
       ss=vnorm(pom,0)
       diff2=ss
       pom(:)=xorig(:)-xatm(:,iat)-atp(:,ipos)
       diff3=vnorm(pom,0)
       rtdiff=diff1/diff3
       if ((abs(bmin(iat))<1.0d-12).or.(ss<bmin(iat))) then
         bmin(iat)=ss
       end if
       pom(:)=xatm(:,iat)+atp(:,ipos)

!      Store more results, for coherent, and portable output
       bcp(nbcp)%iat=iat
       bcp(nbcp)%ipos=ipos
       bcp(nbcp)%chg=chg
       bcp(nbcp)%diff(1)=diff1
       bcp(nbcp)%diff(2)=diff2
       bcp(nbcp)%diff(3)=diff3
       bcp(nbcp)%ev(:)=ev(:)
       bcp(nbcp)%pom(:)=pom(:)
       bcp(nbcp)%rr(:)=rr(:)
       bcp(nbcp)%vec(:,:)=evec(:,:)
       bcp(nbcp)%vv(:)=vv(:)
!      Warning : iat, ipos might be modified by this call
       call vgh_rho(vv,chg,grho,hrho,dist,iat,ipos,0)
       bcp(nbcp)%chg=chg

     end if ! ires==0
   end do srbcp

   if(nbcp>0)then

!    Order the BCP. CPs should appear by increasing values of x,y,z , the latter
!    varying the fastest
     ABI_ALLOCATE(sortguide,(nbcp))
     ABI_ALLOCATE(indexcp,(nbcp))
     ABI_ALLOCATE(cp_tmp,(nbcp))
     do ii=3,1,-1
!      DEBUG
!      write(std_out,*)' cpdrv : sort on index ii=',ii
!      ENDDEBUG

       do jj=1,nbcp
!        DEBUG
!        write(std_out,*)bcp(jj)%vv(:)
!        ENDDEBUG
         sortguide(jj)=bcp(jj)%vv(ii)
         indexcp(jj)=jj
       end do

!      Try to be platform-independent. Might need a larger tolerance.
       call sort_dp(nbcp,sortguide,indexcp,tol3)
       do jj=1,nbcp
         cp_tmp(jj)=bcp(indexcp(jj))
       end do
       do jj=1,nbcp
         bcp(jj)=cp_tmp(jj)
       end do
     end do
!    DEBUG
!    write(std_out,*)' cpdrv : after the sort '
!    do jj=1,nbcp
!    write(std_out,*)bcp(jj)%vv(:)
!    end do
!    ENDDEBUG


!    Output the info about the BCP
     do jj=1,nbcp
       write(untout,'(" Bonded atom (BAT) (indxatm,indxcell,position): ",/,2I6,3F16.8)')&
&       bcp(jj)%iat,bcp(jj)%ipos,bcp(jj)%pom(:)
       write(untout,'("%Bonding CP: ",3F16.8)') bcp(jj)%vv(:)
       write(untout,'("%Eigenval. of Hessian: ",3F16.8)') bcp(jj)%ev(:)
       write(untout,'(a,a,a,3f16.8,a,a,3f16.8,a,a,3f16.8,a)') &
&       ' Eigenvec. of Hessian:',char(10),&
&       '-',bcp(jj)%vec(1,:),char(10),&
&       '-',bcp(jj)%vec(2,:),char(10),&
&       '-',bcp(jj)%vec(3,:),char(10)
       write(untout,'("%Density and laplacian in CP: ",2F16.8)') &
&       bcp(jj)%chg, bcp(jj)%ev(1)+bcp(jj)%ev(2)+bcp(jj)%ev(3)
       write(untout,'("%Relative position of BCP (AT-CP,BAT-CP,AT-BAT,relative(AT): ",/,4F16.8)') &
&       bcp(jj)%diff(:),bcp(jj)%diff(1)/bcp(jj)%diff(3)
       write(untout,*) "********************************************************************"
       write(unto,'(/," BCP: ",3F10.6,3E12.4,E12.4,/)') &
&       bcp(jj)%rr(:),bcp(jj)%ev(:),bcp(jj)%ev(1)+bcp(jj)%ev(2)+bcp(jj)%ev(3)
       write(unto,'(":DISPC ",4F12.6)') bcp(jj)%diff(:),bcp(jj)%diff(1)/bcp(jj)%diff(3)
     end do

     ABI_DEALLOCATE(cp_tmp)
     ABI_DEALLOCATE(indexcp)
     ABI_DEALLOCATE(sortguide)

   end if ! nbcp>0

   if (abs(bmin(inxat))>1.0d-12) then
     rminl(inxat)=aim_dtset%coff1*bmin(inxat)
     r0=bmin(inxat)
   else
     r0=0._dp
   end if

!  !AD-HOC PARAMETER

   do ii=1,natom
     if ((abs(bmin(ii))>1.0d-12).and.(ii /= inxat)) rminl(ii)=aim_dtset%coff2*bmin(ii)
   end do

!  END WARNING

   write(unto,*) ' number of BCP:', nbcp
   write(untout,'(" Number of BCP found: ",I4)') nbcp
   nn=nbcp*(nbcp-1)*(nbcp-2)/6
   if (bit_size(ii) <= nbcp+1) stop ':ERROR b-test !'

!  SEARCHING RCP

   write(unto,*)
   write(unto,*) "RING CRITICAL POINTS (RCP)"
   write(unto,*) "============================="
   write(unto,*)

   write(untout,*)
   write(untout,*) "RING CRITICAL POINTS (RCP)"
   write(untout,*) "============================="
   write(untout,*)

   nrcp=0
   if(aim_dtset%crit==1)nb_now=nbcp
   if(aim_dtset%crit==2)nb_now=nb
!  DEBUG
!  nb_now=nbcp
!  ENDDEBUG
   nn=nb_now*(nb_now-1)/2
   ABI_ALLOCATE(rcp,(nn))

!  Loop on pairs of BCP or atoms
   do ii=1,nb_now-1
     srcp: do jj=ii+1,nb_now
       if (aim_dtset%crit==1) then
         vv(:)=xorig(:)+(bcp(ii)%rr(:)+bcp(jj)%rr(:))/2._dp
       else if (aim_dtset%crit==2) then
         vv(:)=xorig(:)*half+(xatm(:,ibat(ii))+atp(:,ipibat(ii))+xatm(:,ibat(jj))+atp(:,ipibat(jj)))*quarter
       end if

       call critic(aim_dtset,vv,ev,evec,aim_dmaxcs,ires,1)

!      One CP has been found, must make tests to see whether it is a new RCP
       if (ires==0) then
!        Check that it is not the same as a previous one
         if (nrcp > 0) then
           do kk=1,nrcp
             pom(:)=vv(:)-rcp(kk)%rr(:)-xorig(:)
             dist=vnorm(pom,0)
             if (dist < aim_dtset%dpclim) then
               write(unto,*) ':RCP already known'
               cycle srcp
             end if
           end do
         end if
!        Check that it is within the maximum allowed distance for a CP
         rr(:)=vv(:)-xorig(:)
         ss=vnorm(rr,0)
         if (ss > maxcpdst) then
           write(unto,*) 'RCP distance from atom exceed maxcpdst !'
           cycle srcp
         end if
!        Check that it is a RCP
         nn=0
         do kk=1,3
           nn=nn+ev(kk)/abs(ev(kk))
         end do
         if (nn /= 1) then
           write(unto,*) ' the critical point that is found is not a RCP '
           cycle srcp
         end if
!        If crit==2, check that it is on the Bader surface
         if (aim_dtset%crit==2) then
           uu(:)=vv(:)-aim_epstep*rr(:)/ss
           call aim_follow(aim_dtset,uu,aim_npmaxin,srch,iatinit,iposinit,iat,ipos,nstep)
           if ((iat/=inxat).or.(inxcell/=ipos))then
             write(unto,*) ' RCP is not on the Bader surface (outside of it)'
             cycle srcp
           end if
         end if

!        A new RCP has been found !
         nrcp=nrcp+1

         rcp(nrcp)%iat=iat
         rcp(nrcp)%ipos=ipos
         rcp(nrcp)%rr(:)=vv(:)-xorig(:)
         rcp(nrcp)%vec(:,:)=evec(:,:)
         rcp(nrcp)%ev(:)=ev(:)
         rcp(nrcp)%vv(:)=vv(:)
         call vgh_rho(vv,chg,grho,hrho,dist,iat,ipos,0)
         rcp(nrcp)%chg=chg

       end if ! ires==0
     end do srcp ! jj=ii+2,nb_now
   end do ! ii=1,nb_now-1

   if(nrcp>0)then

!    Order the RCP. CPs should appear by increasing values of x,y,z , the latter
!    varying the fastest
     ABI_ALLOCATE(sortguide,(nrcp))
     ABI_ALLOCATE(indexcp,(nrcp))
     ABI_ALLOCATE(cp_tmp,(nrcp))
     do ii=3,1,-1
!      DEBUG
!      write(std_out,*)' cpdrv : sort on index ii=',ii
!      ENDDEBUG
       do jj=1,nrcp

!        DEBUG
!        write(std_out,*)rcp(jj)%vv(:)
!        ENDDEBUG

!        Try to be platform-independent. Might need a larger tolerance.
         sortguide(jj)=rcp(jj)%vv(ii)
         indexcp(jj)=jj
       end do
       call sort_dp(nrcp,sortguide,indexcp,tol3)
       do jj=1,nrcp
         cp_tmp(jj)=rcp(indexcp(jj))
       end do
       do jj=1,nrcp
         rcp(jj)=cp_tmp(jj)
       end do
     end do

!    DEBUG
!    write(std_out,*)' cpdrv : after the sort '
!    do jj=1,nrcp
!    write(std_out,*)rcp(jj)%vv(:)
!    end do
!    ENDDEBUG


!    Write the Ring Critical Point information
     do jj=1,nrcp
       write(untout,'(";Ring CP: ",3F16.8)') rcp(jj)%vv(:)
       write(untout,'("%Eigenval. of Hessian: ",3F16.8)') rcp(jj)%ev(:)
       write(untout,'(a,a,a,3f16.8,a,a,3f16.8,a,a,3f16.8,a)') &
&       ' Eigenvec. of Hessian:',char(10),&
&       '-',rcp(jj)%vec(1,:),char(10),&
&       '-',rcp(jj)%vec(2,:),char(10),&
&       '-',rcp(jj)%vec(3,:),char(10)
       write(untout,'("%Density and laplacian in CP: ",2F16.8)') &
&       rcp(jj)%chg, rcp(jj)%ev(1)+rcp(jj)%ev(2)+rcp(jj)%ev(3)
       write(untout,*) "********************************************************************"
       write(unto,'(/," RCP: ",3F10.6,3E12.4,E12.4,/)') &
&       rcp(jj)%rr(:),rcp(jj)%ev(:),rcp(jj)%ev(1)+rcp(jj)%ev(2)+rcp(jj)%ev(3)
     end do

     ABI_DEALLOCATE(cp_tmp)
     ABI_DEALLOCATE(indexcp)
     ABI_DEALLOCATE(sortguide)

   end if ! nrcp>0

   write(untout,'(" Number of RCP found: ",I4)') nrcp
   write(unto,*) ' Number of RCP:', nrcp

!  SEARCHING CCP

   write(unto,*)
   write(unto,*) "CAGE CRITICAL POINTS (CCP)"
   write(unto,*) "============================="
   write(unto,*)

   write(untout,*)
   write(untout,*) "CAGE CRITICAL POINTS (CCP)"
   write(untout,*) "============================="
   write(untout,*)


   nn=nrcp*(nrcp-1)/2
   ABI_ALLOCATE(ccp,(nn))

   nccp=0
   do ii=1,nrcp-1
     srccp: do jj=ii+1,nrcp
       vv(:)=xorig(:)+(rcp(ii)%rr(:)+rcp(jj)%rr(:))/2._dp
       call critic(aim_dtset,vv,ev,evec,aim_dmaxcs,ires,3)
       if (ires==0) then
         rr(:)=vv(:)-xorig(:)
         ss=vnorm(rr,0)
         if (ss > maxcpdst) then
           write(unto,*) 'CCP distance from atom exceed maxcpdst !'
           cycle srccp
         end if
         nn=0
         do kk=1,3
           nn=nn+ev(kk)/abs(ev(kk))
         end do
         if (nn /= 3) then
           write(unto,*) ' the critical point that is found is not a CCP '
           cycle srccp
         end if

         if (nccp > 0) then
           do kk=1,nccp
             pom(:)=vv(:)-ccp(kk)%rr(:)-xorig(:)
             dist=vnorm(pom,0)
             if (dist < aim_dtset%dpclim) then
               write(unto,*) ':CCP already known'
               cycle srccp
             end if
           end do
         end if
         if (aim_dtset%crit==2) then
           uu(:)=vv(:)-aim_epstep*rr(:)/ss
           call aim_follow(aim_dtset,uu,aim_npmaxin,srch,iatinit,iposinit,iat,ipos,nstep)
           if ((iat/=inxat).or.(inxcell/=ipos)) then
             write(unto,*) ' This CCP is not on the Bader surface (outside of it)'
             cycle srccp
           end if
         end if

         nccp=nccp+1

         ccp(nccp)%iat=iat
         ccp(nccp)%ipos=ipos
         ccp(nccp)%rr(:)=vv(:)-xorig(:)
         ccp(nccp)%vec(:,:)=evec(:,:)
         ccp(nccp)%ev(:)=ev(:)
         ccp(nccp)%vv(:)=vv(:)
         call vgh_rho(vv,chg,grho,hrho,dist,iat,ipos,0)
         ccp(nccp)%chg=chg

       end if
     end do srccp
   end do

   if(nccp>0)then

!    Order the CCP. CPs should appear by increasing values of x,y,z , the latter
!    varying the fastest
     ABI_ALLOCATE(sortguide,(nccp))
     ABI_ALLOCATE(indexcp,(nccp))
     ABI_ALLOCATE(cp_tmp,(nccp))
     do ii=3,1,-1
       do jj=1,nccp
!        Try to be platform-independent. Might need a larger tolerance.
         sortguide(jj)=ccp(jj)%vv(ii)
         indexcp(jj)=jj
       end do
       call sort_dp(nccp,sortguide,indexcp,tol3)
       do jj=1,nccp
         cp_tmp(jj)=ccp(indexcp(jj))
       end do
       do jj=1,nccp
         ccp(jj)=cp_tmp(jj)
       end do
     end do

!    Write the Cage Critical Point information
     do jj=1,nccp
       write(untout,'("%Cage CP: ",3F16.8)') ccp(jj)%vv(:)
       write(untout,'("%Eigenval. of Hessian: ",3F16.8)') ccp(jj)%ev(:)
       write(untout,'(a,a,a,3f16.8,a,a,3f16.8,a,a,3f16.8,a)') &
&       ' Eigenvec. of Hessian:',char(10),&
&       '-',ccp(jj)%vec(1,:),char(10),&
&       '-',ccp(jj)%vec(2,:),char(10),&
&       '-',ccp(jj)%vec(3,:),char(10)
       write(untout,'("%Density and laplacian in CP: ",2F16.8)') &
&       ccp(jj)%chg, ccp(jj)%ev(1)+ccp(jj)%ev(2)+ccp(jj)%ev(3)
       write(untout,*) "********************************************************************"
       write(unto,'(/," CCP: ",3F10.6,3E12.4,E12.4,/)') &
&       ccp(jj)%rr(:),ccp(jj)%ev(:),ccp(jj)%ev(1)+ccp(jj)%ev(2)+ccp(jj)%ev(3)
     end do

     ABI_DEALLOCATE(sortguide)
     ABI_DEALLOCATE(indexcp)
     ABI_DEALLOCATE(cp_tmp)

   end if ! nccp>0

   write(untout,'(" Number of CCP found: ",I4)') nccp
   write(unto,*) 'Number of CCP:', nccp
   write(unto,*)
   write(untout,*)
   write(unto, '(a,3i8)' ) 'BCP-RCP-CCP', nbcp,nrcp,nccp
   write(untout, '(a,3i8)' ) 'BCP-RCP-CCP', nbcp,nrcp,nccp

   write(unto,*)
   write(unto,*) "==============================="
   write(unto,*) "END OF CRITICAL POINTS ANALYSIS"
   write(unto,*)

   write(untout,*)
   write(untout,*) "==============================="
   write(untout,*) "END OF CRITICAL POINTS ANALYSIS"
   write(untout,*)


!  Output of the CPs

   write(untc,'(I4, " :BCP''s, coordinates, laplacian eigs, type of bonding at., sum of lap.eigs., density")') nbcp
   do ii=1,nbcp
     write(untc,'(3F10.6,3E12.4,I4,2E12.4)') &
&     bcp(ii)%rr(:),bcp(ii)%ev(:),bcp(ii)%iat,bcp(ii)%ev(1)+bcp(ii)%ev(2)+bcp(ii)%ev(3),bcp(ii)%chg
   end do

   write(untc,'(I4, " :RCP''s, coordinates, laplacian eigenvalues, sum of these, density")') nrcp
   do ii=1,nrcp
     vv(:)=rcp(ii)%rr(:)+xorig(:)
     call vgh_rho(vv,chg,grho,hrho,dist,iat,ipos,0)
     write(untc,'(3F10.6,3E12.4,2E12.4)') &
&     rcp(ii)%rr(:),rcp(ii)%ev(:),rcp(ii)%ev(1)+rcp(ii)%ev(2)+rcp(ii)%ev(3),rcp(ii)%chg
   end do

   write(untc,'(I4, " :CCP''s coordinates, laplacian eigenvalues, sum of these, density")') nccp
   do ii=1,nccp
     vv(:)=ccp(ii)%rr(:)+xorig(:)
     call vgh_rho(vv,chg,grho,hrho,dist,iat,ipos,0)
     write(untc,'(3F10.6,3E12.4,2E12.4)') &
&     ccp(ii)%rr(:),ccp(ii)%ev(:),ccp(ii)%ev(1)+ccp(ii)%ev(2)+ccp(ii)%ev(3),ccp(ii)%chg
   end do

 end if ! End the condition on aim_dtset%crit > 0

!Reading of the CPs from the file

 if (aim_dtset%crit==-1) then
   read(untc,*) nbcp
   ABI_ALLOCATE(bcp,(nbcp))
   do ii=1,nbcp
     read(untc,*) bcp(ii)%rr(:)
   end do
   read(untc,*) nrcp
   ABI_ALLOCATE(rcp,(nrcp))
   do ii=1,nrcp
     read(untc,*) rcp(ii)%rr(:)
   end do
   read(untc,*) nccp
   ABI_ALLOCATE(ccp,(nccp))
   do ii=1,nccp
     read(untc,*) ccp(ii)%rr(:)
   end do
 end if

 do ii=1,nbcp
   pc(:,ii)=bcp(ii)%rr(:)
   icpc(ii)=-1
 end do
 do ii=1,nrcp
   pc(:,nbcp+ii)=rcp(ii)%rr(:)
   icpc(nbcp+ii)=1
 end do
 do ii=1,nccp
   pc(:,nbcp+nrcp+ii)=ccp(ii)%rr(:)
   icpc(nbcp+nrcp+ii)=3
 end do
 npc=nbcp+nrcp+nccp

!Checking

 if (allocated(bcp)) then
   do ii=1,nbcp
     do jj=1,npc
       iat=bcp(ii)%iat
       ipos=bcp(ii)%ipos
       if ((iat/=0).and.(ipos/=0)) then
         pom(:)=pc(:,jj)+xorig(:)-xatm(:,iat)-atp(:,ipos)
         ss=aim_dtset%coff2*vnorm(pom,0)
         if (rminl(iat) >= ss) rminl(iat)=ss
       end if
     end do
   end do
   ABI_DEALLOCATE(bcp)
 end if
 do ii=1,natom
   write(unto,*) 'atom: ', ii, rminl(ii)
 end do

 if(allocated(rcp))ABI_DEALLOCATE(rcp)
 if(allocated(ccp))ABI_DEALLOCATE(ccp)

!END CP ANALYSIS

 call timein(ttcp,wall)
 ttcp=ttcp-tt0

end subroutine cpdrv
!!***
