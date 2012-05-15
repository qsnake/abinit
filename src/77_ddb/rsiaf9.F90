!{\src2tex{textfont=tt}}
!!****f* ABINIT/rsiaf9
!!
!! NAME
!! rsiaf9
!!
!! FUNCTION
!! Compute the real-space interatomic force constants, including
!! both analytical (short-range) and non-analytical (long-range contribution)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales by which rprim is to be multiplied
!! atifc(natom) =  atifc(ia) equals 1 if the analysis of ifc
!!  has to be done for atom ia; otherwise 0.
!! atmfrc(2,3,natom,3,natom,nrpt)
!!  = Analytical part of the Interatomic Forces in real space.
!!  We used the imaginary part just for debugging
!! dielt(3,3)=dielectric tensor
!! dipdip= if 0, no dipole-dipole interaction was subtracted in atmfrc
!!  if 1, atmfrc has been build without dipole-dipole part
!! dyewq0(3,3,natom)=contraction of the Ewald dynamical matrix at q=0
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! ifcana= 0 => no analysis of ifc ; 1 => full analysis
!! ifcout= Number of interatomic force constants written in the output file
!! iout=unit number for nice output
!! natom=number of atoms in unit cell
!! nrpt= Number of R points in the Big Box
!! nsphere=number of atoms to be included in the cut-off sphere for interatomic
!!  force constant; if = 0 : maximum extent allowed by the grid.
!! rcan(3,natom)=canonical coordinates of atoms
!! rifcsph=radius for cutoff of IFC
!! rprim(3,3)=dimensionless primitive translations in real space
!! rpt(3,nrpt)=canonical coordinates of the points in the BigBox.
!! tcpui,twalli=initial values of cpu and wall clocktime
!! wghatm(natom,natom,nrpt)
!!  = Weights associated to a pair of atoms and to a R vector
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement
!!
!! OUTPUT
!!   written in the output file.
!!
!! NOTES
!! This routine should be executed by one processor only
!!
!! PARENTS
!!      mkifc9
!!
!! CHILDREN
!!      axial9,canct9,dist9,ifclo9,leave_new,matr3inv,sort_dp,timein,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine rsiaf9(acell,atifc,atmfrc,dielt,dipdip,dyewq0,&
& gprim,ifcana,ifcout,iout,natom,nrpt,nsphere,rcan,&
& rifcsph,rprim,rpt,tcpui,twalli,wghatm,zeff)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rsiaf9'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
 use interfaces_77_ddb, except_this_one => rsiaf9
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: dipdip,ifcana,ifcout,iout,natom,nrpt,nsphere
 real(dp),intent(in) :: rifcsph,tcpui,twalli
!arrays
 integer,intent(in) :: atifc(natom)
 real(dp),intent(in) :: acell(3),atmfrc(2,3,natom,3,natom,nrpt),dielt(3,3)
 real(dp),intent(in) :: dyewq0(3,3,natom),gprim(3,3),rcan(3,natom)
 real(dp),intent(in) :: rprim(3,3),rpt(3,nrpt),zeff(3,3,natom)
 real(dp),intent(inout) :: wghatm(natom,natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: flag,ia,ib,ii,index,irpt,jj,kk,mu,nu
 real(dp) :: detdlt,dist1,dist2,rsq,scprod,tcpu,trace1,trace2,trace3
 real(dp) :: twall,yy
 character(len=500) :: message
!arrays
 integer,allocatable :: list(:),sorted(:)
 real(dp) :: ewiaf0(3,3),ewiaf1(3,3),ewloc(3,3),ifcloc(3,3),invdlt(3,3),ra(3)
 real(dp) :: rbpos(3),rcart(3),rdiff(3),rsiaf(3,3),rsloc(3,3),sriaf(3,3)
 real(dp) :: srloc(3,3),vect1(3),vect2(3),vect3(3),work(3),xred(3),xx(3)
 real(dp),allocatable :: dist(:,:,:),wkdist(:)

! *********************************************************************

!call timer(34,1,wall)

!Calculating the inverse (transpose) of the dielectric tensor
 call matr3inv(dielt,invdlt)

!Calculating the determinant of the dielectric tensor
 detdlt=dielt(1,1)*dielt(2,2)*dielt(3,3)+dielt(1,3)*dielt(2,1)*&
& dielt(3,2)+dielt(1,2)*dielt(2,3)*dielt(3,1)-dielt(1,3)*&
& dielt(2,2)*dielt(3,1)-dielt(1,1)*dielt(2,3)*dielt(3,2)-&
& dielt(1,2)*dielt(2,1)*dielt(3,3)

!Compute the distances between atoms
 ABI_ALLOCATE(dist,(natom,natom,nrpt))
 call dist9(acell,dist,gprim,natom,nrpt,rcan,rprim,rpt)
!Now dist(ia,ib,irpt) contains the distance from atom ia
!to atom ib in unit cell irpt.

 write(std_out,'(a)' )' rsiaf9 : analysis of interatomic force constants '
 write(iout, '(/,a,/)' ) &
& ' Analysis of interatomic force constants '
 if(dipdip==1)then
   write(iout, '(a)' )' Are given : column(1-3), the total force constant'
   write(iout, '(a)' )'       then  column(4-6), the Ewald part'
   write(iout, '(a)' )'       then  column(7-9), the short-range part'
   write(iout, '(a)' )' Column 1, 4 and 7 are related to the displacement'
   write(iout, '(a)' )'       of the generic atom along x,               '
   write(iout, '(a)' )' column 2, 5 and 8 are related to the displacement'
   write(iout, '(a)' )'       of the generic atom along y,               '
   write(iout, '(a)' )' column 3, 6 and 9 are related to the displacement'
   write(iout, '(a)')'       of the generic atom along z.               '
 else if(dipdip==0)then
   write(iout, '(a)' )' column 1 is related to the displacement'
   write(iout, '(a)' )'        of the generic atom along x,    '
   write(iout, '(a)' )' column 2 is related to the displacement'
   write(iout, '(a)' )'        of the generic atom along y,    '
   write(iout, '(a)' )' column 3 is related to the displacement'
   write(iout, '(a)' )'        of the generic atom along z,    '
 end if

 ABI_ALLOCATE(list,(natom*nrpt))
 ABI_ALLOCATE(sorted,(natom*nrpt))

!BIG loop on all generic atoms
 do ia=1,natom

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-before big ia loop at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')

!  First transform canonical coordinates to reduced coordinates
   do ii=1,3
     xred(ii)=gprim(1,ii)*rcan(1,ia)+gprim(2,ii)*rcan(2,ia)&
&     +gprim(3,ii)*rcan(3,ia)
   end do
!  Then to cartesian coordinates
   ra(:)=xred(1)*acell(1)*rprim(:,1)+&
&   xred(2)*acell(2)*rprim(:,2)+&
&   xred(3)*acell(3)*rprim(:,3)

!  DEBUG
!  write(std_out,*) ' nsphere, rifcsph, atifc(ia) =', nsphere, rifcsph, atifc(ia)
!  ENDDEBUG

   if( (nsphere/=0.and.nsphere<natom*nrpt) .or. rifcsph > tol10 .or. atifc(ia)==1)then
!    This sorting algorithm is slow ...
     ABI_ALLOCATE(wkdist,(natom*nrpt))
     wkdist(:)=reshape(dist(ia,:,:),(/natom*nrpt/))
     do ii=1,natom*nrpt
       list(ii)=ii
     end do
     call sort_dp(natom*nrpt,wkdist,list,tol14)
   end if

!  In case some cut-off has to be imposed on the IFCs,
!  zero the outside IFCs now : act on wghatm
   if(nsphere/=0.and.nsphere<natom*nrpt)then
     do ii=nsphere+1,natom*nrpt
       index=list(ii)
       irpt=(index-1)/natom+1
       ib=index-natom*(irpt-1)
       wghatm(ia,ib,irpt)=0._dp
     end do
   end if

   if(rifcsph>tol10)then
!    DEBUG
!    write(std_out,*) ' rsiaf9 : wkdist = '
!    write(std_out,'(4E16.6)') wkdist
!    ENDDEBUG
     do ii=nsphere+1,natom*nrpt
       index=list(ii)
!      preserve weights for atoms inside sphere of radius rifcsph
       if (wkdist(ii) < rifcsph) cycle
       irpt=(index-1)/natom+1
       ib=index-natom*(irpt-1)
       wghatm(ia,ib,irpt)=0._dp
     end do
   end if

!  deallocate wkdist if used
   if(allocated(wkdist))  then
     ABI_DEALLOCATE(wkdist)
   end if

   if(atifc(ia)==1)then

     write(iout, '(a)' )
     write(std_out,'(a,i4)' )' generic atom number',ia
     write(iout, '(a,i4)' )' generic atom number',ia
     write(std_out,'(a,3es16.8)' ) ' with cartesian coordinates',ra(1:3)
     write(iout,'(a,3es16.8)' ) ' with cartesian coordinates',ra(1:3)
     write(iout, '(a)' )

     if(ifcana==1)then
!      Generate the local coordinate system for the atom ia
       index=list(2)
       write(std_out,*)index
       call canct9(acell,gprim,ib,index,irpt,natom,nrpt,&
&       rcan,rcart,rprim,rpt)
       dist1=dist(ia,ib,irpt)
       vect2(1)=rcart(1)-ra(1)
       vect2(2)=rcart(2)-ra(2)
       vect2(3)=rcart(3)-ra(3)
       flag=0
       do ii=3,natom*nrpt
         index=list(ii)
         call canct9(acell,gprim,ib,index,irpt,natom,nrpt,&
&         rcan,rcart,rprim,rpt)
         dist2=dist(ia,ib,irpt)
         vect1(1)=(rcart(1)-ra(1))-vect2(1)
         vect1(2)=(rcart(2)-ra(2))-vect2(2)
         vect1(3)=(rcart(3)-ra(3))-vect2(3)
         scprod=0.0_dp
         do jj=1,3
           scprod=scprod+vect1(jj)**2
         end do
         do jj=1,3
           vect1(jj)=vect1(jj)/scprod**0.5
         end do
         scprod=0.0_dp
         do jj=1,3
           scprod=scprod+vect2(jj)*vect1(jj)
         end do
         do jj=1,3
           work(jj)=vect2(jj)-vect1(jj)*scprod
         end do
         scprod=0.0_dp
         do jj=1,3
           scprod=scprod+work(jj)**2
         end do
         if(scprod>1.0d-10)then
           flag=1
         end if
         if(flag==1)exit
       end do
       if(flag==0)then
         write(message, '(5a)' )&
&         ' rsiaf9 : BUG -',ch10,&
&         '  Unable to find a third atom not aligned',ch10,&
&         '  with the two selected ones.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
       vect2(1)=work(1)/scprod**0.5
       vect2(2)=work(2)/scprod**0.5
       vect2(3)=work(3)/scprod**0.5
       vect3(1)=vect1(2)*vect2(3)-vect1(3)*vect2(2)
       vect3(2)=vect1(3)*vect2(1)-vect1(1)*vect2(3)
       vect3(3)=vect1(1)*vect2(2)-vect1(2)*vect2(1)
       write(iout, '(a)' )' Third atom defining local coordinates : '
       write(iout, '(a,i4,a,i4)' )'     ib = ',ib,'   irpt = ',irpt
     end if

!    Analysis and output of
!    force constants, ordered with respect to the distance
!    from atom ia
     do ii=1,ifcout
       index=list(ii)
       call canct9(acell,gprim,ib,index,irpt,natom,nrpt,&
&       rcan,rbpos,rprim,rpt)
       write(iout, '(a)' )
       write(iout, '(i4,a,i6,a,i8)' )&
&       ii,' interaction with atom',ib,' cell',irpt
       write(iout, '(a,3es16.6)' )&
&       ' with coordinates ',rbpos(1:3)*(one+tol8)
       write(iout, '(a,es16.6)' )&
&       ' and distance ',dist(ia,ib,irpt)

       if(ifcana==1.and.ii/=1)then
         dist1=dist(ia,ib,irpt)
         vect1(1)=(rbpos(1)-ra(1))/dist1
         vect1(2)=(rbpos(2)-ra(2))/dist1
         vect1(3)=(rbpos(3)-ra(3))/dist1
       end if

       if(dipdip==0)then

!        Get the "total" force constants (=real space FC)
!        without taking into account the dipole-dipole interaction
         do mu=1,3
           do nu=1,3
             rsiaf(mu,nu)=atmfrc(1,mu,ia,nu,ib,irpt)&
&             *wghatm(ia,ib,irpt)
           end do
         end do

!        Output of the ifcs in cartesian coordinates
         do nu=1,3
           write(iout, '(1x,3f9.5)' )(rsiaf(mu,nu)+tol10,mu=1,3)
         end do

         if(ifcana==1)then
!          Further analysis
           trace1=rsiaf(1,1)+rsiaf(2,2)+rsiaf(3,3)
           write(iout, '(a,f9.5)' ) '  Trace         ',trace1+tol10
           if(ii/=1)then
             call axial9(rsiaf,vect1,vect2,vect3)
           end if
           write(iout, '(a)' )' Transformation to local coordinates '
           write(iout, '(a,3f16.6)' ) ' First  local vector :',vect1
           write(iout, '(a,3f16.6)' ) ' Second local vector :',vect2
           write(iout, '(a,3f16.6)' ) ' Third  local vector :',vect3
           call ifclo9(rsiaf,ifcloc,vect1,vect2,vect3)
           do nu=1,3
             write(iout, '(1x,3f9.5)' )(ifcloc(mu,nu)+tol10,mu=1,3)
           end do
!          Further analysis finished
         end if

       else if(dipdip==1)then

!        Get the Coulomb part
         do jj=1,3
           rdiff(jj)=ra(jj)-rbpos(jj)
         end do
         rsq=0.0_dp
         xx(1:3)=0.0_dp
         do jj=1,3
           do kk=1,3
             ewiaf0(jj,kk)=0.0_dp
             rsq=rsq+rdiff(jj)*invdlt(kk,jj)*rdiff(kk)
             xx(kk)=xx(kk)+invdlt(kk,jj)*rdiff(jj)
           end do
         end do
         yy=sqrt(rsq)
!        Avoid zero denominators in term:
         if (sqrt(rsq)>=tol12) then
           do mu=1,3
             do nu=1,3
               ewiaf0(mu,nu)=(-3*xx(nu)*xx(mu)+invdlt(nu,mu)*yy**2)&
&               /yy**5/sqrt(detdlt)
             end do
           end do
         else
           if (ia/=ib)then
             write(message, '(a,a,a,a,a,a,a,i5,a,i5,a)' )&
&             ' rsiaf9 : ERROR -',ch10,&
&             '  The distance between two atoms vanishes.',ch10,&
&             '  This is not allowed.',ch10,&
&             '  Action : check the input for the atoms number',ia,' and',ib,'.'
             call wrtout(std_out,message,'COLL')
             call leave_new('COLL')
           end if
         end if

!        Take into account the effective charge tensor
         do mu=1,3
           do nu=1,3
             ewiaf1(mu,nu)=0.0_dp
             if(ii==1)then
               ewiaf1(mu,nu)=-dyewq0(mu,nu,ia)
             end if
             do jj=1,3
               do kk=1,3
                 ewiaf1(mu,nu)=ewiaf1(mu,nu)&
&                 +zeff(jj,mu,ia)*zeff(kk,nu,ib)*&
&                 ewiaf0(jj,kk)
               end do
             end do
           end do
         end do

!        Get the short-range force constants and the
!        "total" force constants (=real space FC)
         do mu=1,3
           do nu=1,3
             sriaf(mu,nu)=atmfrc(1,mu,ia,nu,ib,irpt)&
&             *wghatm(ia,ib,irpt)
             rsiaf(mu,nu)=ewiaf1(mu,nu)+sriaf(mu,nu)
           end do
         end do

!        Output of the results
         do nu=1,3
           write(iout, '(1x,3(3f9.5,1x))' )&
&           (rsiaf(mu,nu) +tol10,mu=1,3),&
&           (ewiaf1(mu,nu)+tol10,mu=1,3),&
&           (sriaf(mu,nu) +tol10,mu=1,3)
         end do

         if(ifcana==1)then
!          Further analysis
           write(iout, '(a)' )' Traces (and ratios) :'
           trace1=rsiaf(1,1)+rsiaf(2,2)+rsiaf(3,3)
           trace2=ewiaf1(1,1)+ewiaf1(2,2)+ewiaf1(3,3)
           trace3=sriaf(1,1)+sriaf(2,2)+sriaf(3,3)
           write(iout, '(3(f9.5,17x))' )trace1+tol10,trace2+tol10,trace3+tol10
           write(iout, '(3(f9.5,17x))' )&
&           1.0,trace2/trace1+tol10,trace3/trace1+tol10

           if(ii/=1)then
             call axial9(rsiaf,vect1,vect2,vect3)
           end if
           write(iout, '(a)' )' Transformation to local coordinates '
           write(iout, '(a,3f16.6)' )&
&           ' First  local vector :',vect1
           write(iout, '(a,3f16.6)' )&
&           ' Second local vector :',vect2
           write(iout, '(a,3f16.6)' )&
&           ' Third  local vector :',vect3
           call ifclo9(rsiaf,rsloc,vect1,vect2,vect3)
           call ifclo9(ewiaf1,ewloc,vect1,vect2,vect3)
           call ifclo9(sriaf,srloc,vect1,vect2,vect3)
           do nu=1,3
             write(iout, '(1x,3(3f9.5,1x))' )&
&             (rsloc(mu,nu)+tol10,mu=1,3),&
&             (ewloc(mu,nu)+tol10,mu=1,3),&
&             (srloc(mu,nu)+tol10,mu=1,3)
           end do
           if(ii/=1)then
             write(iout, '(a)' )' Ratio with respect to the longitudinal ifc'
           else
             write(iout, '(a)' )' Ratio with respect to the (1,1) element'
           end if
           do nu=1,3
             write(iout, '(1x,3(3f9.5,1x))' )&
&             (rsloc(mu,nu)/rsloc(1,1)+tol10,mu=1,3),&
&             (ewloc(mu,nu)/rsloc(1,1)+tol10,mu=1,3),&
&             (srloc(mu,nu)/rsloc(1,1)+tol10,mu=1,3)
           end do

!          Further analysis finished
         end if

!        End the condition on dipdip
       end if

!      End loop over all atoms in BigBox:
     end do

!    End Big loop on atoms in the unit cell,
!    and corresponding test .
   end if
 end do

 ABI_DEALLOCATE(dist)
 ABI_DEALLOCATE(list)
 ABI_DEALLOCATE(sorted)

end subroutine rsiaf9
!!***
