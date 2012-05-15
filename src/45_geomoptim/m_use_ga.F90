!{\src2tex{textfont=tt}}
!!****f* ABINIT/predict_ga
!! NAME
!! predict_ga
!!
!! FUNCTION
!! Given a given set of images, which represent a population, it predicts a new sef of images.
!! The implementation is based on a Genetic Algorithm idea, where the best fit candidates are passed
!! to the next generation. Those are chosen from 20% best fit and 80% from Genetic rules
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! itimimage=number of the current time for image propagation (itimimage+1 is to be predicted here)
!! list_dynimage(nimage)=list of dynamical images. 
!! This is quite useful when ground states of the A and B states is known
!! natom=dimension of vel_timimage and xred_timimage
!! ndynimage=number of dynamical images
!! nimage= population size
!! ntimimage=dimension of several arrays
!! results_gs_timimage(ntimimage,nimage)=datastructure that hold all the history of previous computations.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! acell_timimage(3,nimage,ntimimage)
!!   at input, history of the values of acell for all images, up to itimimage
!!   at output, the predicted values of acell for all images, stored in acell_timimage(3,nimage,itimimage+1)
!! rprim_timimage(3,3,nimage,ntimimage)
!!   at input, history of the values of rprim for all images, up to itimimage
!!   at output, the predicted values of rprim for all images, stored in rprim_timimage(3,nimage,itimimage+1)
!! vel_timimage(3,natom,nimage,ntimimage)
!!   at input, history of the values of vel for all images, up to itimimage
!!   at output, the predicted values of vel for all images, stored in vel_timimage(3,natom,nimage,itimimage+1)
!! xred_timimage(3,natom,nimage,ntimimage)
!!   at input, history of the values of xred for all images, up to itimimage
!!   at output, the predicted values of xred for all images, stored in xred_timimage(3,natom,nimage,itimimage+1)
!!
!! PARENTS
!!      predictimg
!!
!! CHILDREN
!!      initialize_perm,metric,sort_dp,swap
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_use_ga

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_results_img, only : results_img_type,gather_array_img
 implicit none


CONTAINS

subroutine predict_ga(itimimage,idum,natom,nimage,ntimimage,results_img)

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'predict_ga'
 use interfaces_28_numeric_noabirule
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)     :: itimimage,natom,nimage,ntimimage
 integer,intent(inout)  :: idum
!arrays
 type(results_img_type) :: results_img(nimage,ntimimage)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk,itmp,iimage,indiv,oper,father1,father2,ndimen,nsurvivor

!arrays
! store the energy and the enthalpy of all elements of the population
 real(dp) :: rprimd(3,3),gprimd(3,3),rmet(3,3),coori(3),gmet(3,3)

 integer,allocatable  :: iperm(:),japerm(:),iaperm(:),ihperm(:)

 real(dp),allocatable :: etotal_img(:),enthalpy_img(:)
 real(dp),allocatable :: acell(:,:),acell_old(:,:),rprim(:,:,:),rprim_old(:,:,:)
 real(dp),allocatable :: fitness(:),zcoori(:),zcoorj(:)
 real(dp),allocatable :: coor(:,:),coor_old(:,:)
 real(dp),allocatable :: distances(:,:)

!real quantities
 real(dp) :: smallacell,sumH,Hmin,Hmax,denom,num_random,ucvol 

!character(len=500)   :: message
!arrays

! *************************************************************************

!DEBUG
!write(std_out,*)' MODULE predict_ga : enter '
!ENDDEBUG

! gen dimension 

 ndimen=3*natom

 ABI_ALLOCATE(coor,(ndimen,nimage))
 ABI_ALLOCATE(acell,(3,nimage))
 ABI_ALLOCATE(rprim,(3,3,nimage))
 ABI_ALLOCATE(coor_old,(ndimen,nimage))
 ABI_ALLOCATE(acell_old,(3,nimage))
 ABI_ALLOCATE(rprim_old,(3,3,nimage))
 ABI_ALLOCATE(iperm,(nimage))
 ABI_ALLOCATE(japerm,(natom))
 ABI_ALLOCATE(iaperm,(natom))
 ABI_ALLOCATE(ihperm,(nimage))
 ABI_ALLOCATE(etotal_img,(nimage))
 ABI_ALLOCATE(enthalpy_img,(nimage))
 ABI_ALLOCATE(fitness,(nimage))
 ABI_ALLOCATE(distances,(natom,nimage))
 ABI_ALLOCATE(zcoori,(natom))
 ABI_ALLOCATE(zcoorj,(natom))

 call initialize_perm(iperm,nimage)
 call initialize_perm(ihperm,nimage)

!from gs_results_image take energies, rprim and acell and build energy and enthalpy vectors reordered from larger to smaller

 write(187,*) '----------------ITERATION = ',itimimage

 do iimage=1,nimage
   etotal_img(iimage)=results_img(iimage,itimimage)%results_gs%etotal
   do ii=1,natom
     coor_old((ii-1)*3+1:ii*3,iimage)=results_img(iimage,itimimage)%xred(:,ii)
   enddo
   acell_old(:,iimage)=results_img(iimage,itimimage)%acell(:)
   rprim_old(:,:,iimage)=results_img(iimage,itimimage)%rprim(:,:)
   do ii=1,3
     do jj=1,3
       rprimd(ii,jj)=rprim_old(ii,jj,iimage)*acell_old(jj,iimage)
     end do
   end do
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
   enthalpy_img(iimage)=etotal_img(iimage)+sum(results_img(iimage,itimimage)%results_gs%strten(1:3))*ucvol
 enddo

!sort energies

 call sort_dp(nimage,etotal_img,iperm,tol9)

!sort enthalpies

 call sort_dp(nimage,enthalpy_img,ihperm,tol9)

! Fitness is calculated

 sumH=0.d0
 Hmin=minval(enthalpy_img)
 Hmax=maxval(enthalpy_img)
 fitness = zero
 
 do iimage=1,nimage
    write(187,*) 'ETOTAL(',iimage,')=',etotal_img(iimage)
    write(187,*) 'ENTHALPY(',iimage,')=',enthalpy_img(iimage)
    write(187,*) 'acell',acell_old(:,iperm(iimage))
    write(187,*) 'rprim'
    write(187,*) rprim_old(1,:,ihperm(iimage))
    write(187,*) rprim_old(2,:,ihperm(iimage))
    write(187,*) rprim_old(3,:,ihperm(iimage))
    if (etotal_img(iimage)< zero) fitness(iimage)=(one-tanh(two*(enthalpy_img(iimage)-Hmax)/(Hmin-Hmax)-one))/two
    sumH=sumH+fitness(iimage)
    if (iimage>1) fitness(iimage)=fitness(iimage-1)+fitness(iimage)
 enddo

 fitness=fitness/sumH

!  do a single boocle of GA

 indiv=0

! Selection over the best 20% of the population

 nsurvivor=0.2*nimage
 if (nsurvivor < one) nsurvivor=1

write(187,*) 'Survivors = ',nsurvivor

! pass coordinates,rprim, acell of survivors to next generation

 do iimage=1,nsurvivor
   indiv=indiv+1
   coor(:,iimage)=coor_old(:,ihperm(iimage))
   acell(:,iimage)=acell_old(:,ihperm(iimage))
   rprim(:,:,iimage)=rprim_old(:,:,ihperm(iimage))
   do jj=1,natom
     coori(:)=coor((jj-1)*3+1:jj*3,iimage)
     distances(jj,iimage)=sqrt(dot_product(coori,coori))
   enddo
   call initialize_perm(japerm,natom)
   call sort_dp(natom,distances(:,iimage),japerm,tol9)
 enddo

! complete the number of individuals of the generation by choosing them through GA

 do while(indiv<nimage)
   oper=9*uniformrandom(idum)+1
   select case(oper)
     case(1) !two point crossover 
       indiv=indiv+1; num_random=uniformrandom(idum)
       father1=choosefather(num_random,fitness,nimage)
       ii=ndimen*uniformrandom(idum)+1; kk=ndimen*uniformrandom(idum)+1
       if(ii>kk) call swap(ii,kk)
       do jj=1,ndimen
             if(ii<=jj.and.jj<=kk)then
                coor(jj,indiv)=coor_old(ii+kk-jj,ihperm(father1))
             else
                coor(jj,indiv)=coor_old(jj,ihperm(father1))
             end if
       enddo
       do jj=1,natom
         coori(:)=coor((jj-1)*3+1:jj*3,indiv)
         distances(jj,indiv)=sqrt(dot_product(coori,coori))
       enddo
       call initialize_perm(japerm,natom)
       call sort_dp(natom,distances(:,indiv),japerm,tol9)
       smallacell=minval(acell_old(:,ihperm(father1)))
       itmp=comp_indiv(distances,indiv,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv),smallacell)
       if (itmp == 2) then
          write(187,*) 'accepted1'
          acell(:,indiv)=acell_old(:,ihperm(father1))
          rprim(:,:,indiv)=rprim_old(:,:,ihperm(father1))
          do jj=1,ndimen,3
            write(187,*) coor(jj:jj+2,indiv)
          enddo
       else
          indiv=indiv-1
          write(187,*) 'Rejected1'
       endif
     case(2) ! random crossing from two parents
       if ((indiv+2)<=nimage) then
         indiv=indiv+2; 
         do
           num_random=uniformrandom(idum); father1=choosefather(num_random,fitness,nimage)
           num_random=uniformrandom(idum); father2=choosefather(num_random,fitness,nimage)
           if (father1 /= father2) EXIT
         enddo
         do jj=1,ndimen
            num_random=uniformrandom(idum)
            if(num_random<=.5)then
              coor(jj,indiv-1)=coor_old(jj,ihperm(father1))
              if(indiv<=nimage) coor(jj,indiv)=coor_old(jj,ihperm(father2))
            else
              coor(jj,indiv-1)=coor_old(jj,ihperm(father2))
              if(indiv<=nimage) coor(jj,indiv)=coor_old(jj,ihperm(father1))
            endif
         end do
         do jj=1,natom
           coori(:)=coor((jj-1)*3+1:jj*3,indiv-1)
           distances(jj,indiv-1)=sqrt(dot_product(coori,coori))
           coori(:)=coor((jj-1)*3+1:jj*3,indiv)
           distances(jj,indiv)=sqrt(dot_product(coori,coori))
         enddo
         call initialize_perm(japerm,natom)
         call sort_dp(natom,distances(:,indiv-1),japerm,tol9)
         smallacell=minval(acell_old(:,iperm(father1)))
         itmp=comp_indiv(distances,indiv-1,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv-1),smallacell)
         if (itmp == 2) then
            acell(:,indiv-1)=acell_old(:,iperm(father1)) 
            rprim(:,:,indiv-1)=rprim_old(:,:,iperm(father1)) 
         else 
            indiv=indiv-1
         endif
         call initialize_perm(japerm,natom)
         call sort_dp(natom,distances(:,indiv),japerm,tol9)
         smallacell=minval(acell_old(:,iperm(father2)))
         itmp=comp_indiv(distances,indiv,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv),smallacell)
         if (itmp == 2) then
            acell(:,indiv)=acell_old(:,iperm(father2))
            rprim(:,:,indiv)=rprim_old(:,:,iperm(father2))
         else
            indiv=indiv-1
         endif
       endif
     case(3)! single-point crossover
       if ((indiv+2)<=nimage) then
         indiv=indiv+2
         do
           num_random=uniformrandom(idum); father1=choosefather(num_random,fitness,nimage)
           num_random=uniformrandom(idum); father2=choosefather(num_random,fitness,nimage)
           if (father1 /= father2) EXIT
         enddo
         ii=ndimen*uniformrandom(idum)+1
         do jj=1,ii
             coor(jj,indiv-1)=coor_old(jj,ihperm(father1))
             if (indiv<=nimage) then
               coor(jj,indiv)=coor_old(jj,ihperm(father2))
             endif
         enddo
         do jj=ii+1,ndimen
             coor(jj,indiv-1)=coor_old(jj,ihperm(father2))
             if (indiv<=nimage) then
               coor(jj,indiv)=coor_old(jj,ihperm(father1))
             endif
         enddo
         do jj=1,natom
           coori(:)=coor((jj-1)*3+1:jj*3,indiv-1)
           distances(jj,indiv-1)=sqrt(dot_product(coori,coori))
           coori(:)=coor((jj-1)*3+1:jj*3,indiv)
           distances(jj,indiv)=sqrt(dot_product(coori,coori))
         enddo
         call initialize_perm(japerm,natom)
         call sort_dp(natom,distances(:,indiv-1),japerm,tol9)
         smallacell=minval(acell_old(:,iperm(father1)))
         itmp=comp_indiv(distances,indiv-1,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv-1),smallacell)
         if (itmp == 2) then
           acell(:,indiv-1)=acell_old(:,iperm(father1))
           rprim(:,:,indiv-1)=rprim_old(:,:,iperm(father1))
         else
           indiv=indiv-1
         endif
         call initialize_perm(japerm,natom)
         call sort_dp(natom,distances(:,indiv),japerm,tol9)
         smallacell=minval(acell_old(:,iperm(father2)))
         itmp=comp_indiv(distances,indiv,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv),smallacell)
         if (itmp == 2) then
            acell(:,indiv)=acell_old(:,iperm(father2))
            rprim(:,:,indiv)=rprim_old(:,:,iperm(father2))
         else
            indiv=indiv-1
         endif
       endif
     case(4) ! coordinates mutation
       indiv=indiv+1; num_random=uniformrandom(idum); father1=choosefather(num_random,fitness,nimage)
       coor(:,indiv)=coor_old(:,iperm(father1))
       itmp=natom/4
       if (itmp<1) itmp=1
       do jj=1,itmp
         ii=natom*uniformrandom(idum)
         do kk=1,3
           coor(3*ii+kk,indiv)=coor_old(3*ii+kk,iperm(father1))+0.15*uniformrandom(idum)
           if (coor(3*ii+kk,indiv)>1.0_dp) coor(3*ii+kk,indiv)=coor(3*ii+kk,indiv)-1.0_dp
           if (coor(3*ii+kk,indiv)<0.0_dp) coor(3*ii+kk,indiv)=coor(3*ii+kk,indiv)+1.0_dp
         enddo
       enddo
       do jj=1,natom
         coori(:)=coor((jj-1)*3+1:jj*3,indiv)
         distances(jj,indiv)=sqrt(dot_product(coori,coori))
       enddo
       smallacell=minval(acell_old(:,iperm(father1)))
       itmp=comp_indiv(distances,indiv,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv),smallacell)
       if (itmp == 2) then
          acell(:,indiv)=acell_old(:,iperm(father1))
          rprim(:,:,indiv)=rprim_old(:,:,iperm(father1))
       else
          indiv=indiv-1
       endif
     case(5) ! random volume scaling
        indiv=indiv+1; num_random=uniformrandom(idum); father1=choosefather(num_random,fitness,nimage)
        denom=1.0_dp+0.2_dp*(2._dp*uniformrandom(idum)-1._dp)
        acell(:,indiv)=acell_old(:,iperm(father1))*denom
        coor(:,indiv)=coor_old(:,iperm(father1))
        do jj=1,natom
          coori(:)=coor((jj-1)*3+1:jj*3,indiv)
          distances(jj,indiv)=sqrt(dot_product(coori,coori))
        enddo
        call initialize_perm(japerm,natom)
        call sort_dp(natom,distances(:,indiv),japerm,tol9)
        smallacell=minval(acell_old(:,iperm(father1)))
        itmp=comp_indiv(distances,indiv,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv),smallacell)
        if (itmp == 2) then
          coor(:,indiv)=coor_old(:,iperm(father1))
          rprim(:,:,indiv)=rprim_old(:,:,iperm(father1))
        else
           indiv=indiv-1
        endif
      case (6) ! randomizing rprimd
!      call flush(187)
       indiv=indiv+1; num_random=uniformrandom(idum); father1=choosefather(num_random,fitness,nimage)
       coor(:,indiv)=coor_old(:,iperm(father1))
       do ii=1,3
         rprim(ii,ii,indiv)=rprim(ii,ii,indiv)+0.5_dp*(two*uniformrandom(idum)-one)
       enddo
       rprim(1,2,indiv)=0.5_dp*(two*uniformrandom(idum)-one)
       rprim(2,1,indiv)=rprim(1,2,indiv)
       rprim(1,3,indiv)=0.5_dp*(two*uniformrandom(idum)-one)
       rprim(3,1,indiv)=rprim(1,3,indiv)
       rprim(2,3,indiv)=0.5_dp*(two*uniformrandom(idum)-one)
       rprim(3,2,indiv)=rprim(2,3,indiv)
!     case(5) ! arithmetic average
!        indiv=indiv+1; 
!        do
!          num_random=uniformrandom(idum); father1=choosefather(num_random,fitness,nimage)
!          num_random=uniformrandom(idum); father2=choosefather(num_random,fitness,nimage)
!          if (father1 /= father2) EXIT
!        enddo
!        do ii=1,ndimen
!            coor(ii,indiv)=(coor_old(ii,iperm(father1))+coor_old(ii,iperm(father2)))/two
!        enddo
!        do jj=1,natom
!          coori(:)=coor((jj-1)*3+1:jj*3,indiv)
!          distances(jj,indiv)=sqrt(dot_product(coori,coori))
!        enddo
!        smallacell=minval(acell(:,iperm(father1)))
!        itmp=comp_indiv(distances,indiv,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv),smallacell)
!        if (itmp>0) then
!           indiv=indiv-1
!          write(187,*) 'Rejected 5 '
!        else
!          write(187,*) 'Accepted 5 '
!          rprim(:,:,indiv)=rprim_old(:,:,iperm(father1))
!          acell(:,indiv)=(acell_old(:,iperm(father1))+acell_old(:,iperm(father2)))/two
!          do jj=1,ndimen,3
!            write(187,*) coor(jj:jj+2,indiv)
!          enddo
!        endif
!     case(6) ! geometric average
!        indiv=indiv+1; 
!        do
!          num_random=uniformrandom(idum); father1=choosefather(num_random,fitness,nimage)
!          num_random=uniformrandom(idum); father2=choosefather(num_random,fitness,nimage)
!          if (father1 /= father2) EXIT
!        enddo
!        do ii=1,ndimen
!            coor(ii,indiv)=sqrt(coor_old(ii,iperm(father1))*coor_old(ii,iperm(father2)))
!        enddo
!        do jj=1,natom
!          coori(:)=coor((jj-1)*3+1:jj*3,indiv)
!          distances(jj,indiv)=sqrt(dot_product(coori,coori))
!        enddo
!        smallacell=minval(acell(:,iperm(father1)))
!        itmp=comp_indiv(distances,indiv,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv),smallacell)
!        if (itmp>0) then
!           indiv=indiv-1
!          write(187,*) 'Rejected 6 '
!        else
!          write(187,*) 'Accepted 6 '
!          rprim(:,:,indiv)=rprim_old(:,:,iperm(father1))
!          acell(:,indiv)=sqrt(acell_old(:,iperm(father1))*acell_old(:,iperm(father2)))
!          do jj=1,ndimen,3
!            write(187,*) coor(jj:jj+2,indiv)
!          enddo
!        endif
!      case (8) ! slicing of crystal structure
!       call flush(187)
!        if ((indiv+2)<=nimage) then
!          indiv=indiv+2
!          do 
!             num_random=uniformrandom(idum); father1=choosefather(num_random,fitness,nimage)
!             num_random=uniformrandom(idum); father2=choosefather(num_random,fitness,nimage)
!             if (father1 /= father2) EXIT
!          enddo
!          idir=3*uniformrandom(idum) ! assuming values between 0 and 2
!          denom=uniformrandom(idum)  ! choose plane "height"
!          ii=0
!          do jj=1,natom
!            zcoori(jj)=coor_old(3*(jj-1)+1+idir,iperm(father1))
!            zcoorj(jj)=coor_old(3*(jj-1)+1+idir,iperm(father2))
!            if (zcoori(jj) < denom) ii=ii+1 
!          enddo
!          write(187,*) 'Estoy en 6'
!          call initialize_perm(japerm,natom)
!          call initialize_perm(iaperm,natom)
!          call sort_dp(natom,zcoori,iaperm,tol9)
!          call sort_dp(natom,zcoorj,japerm,tol9)
!          kk=natom-ii
!          write(187,*) indiv,father1,father2,idir,denom,ii,kk
!          itmp=0
!          do jj=1,ii
!            itmp=itmp+1
!            iatom=iaperm(jj)
!            iatomj=japerm(jj)
!            coor(3*(itmp-1)+1:3*itmp,indiv-1)=coor_old(3*(iatom-1)+1:3*iatom,iperm(father1))
!            coor(3*(itmp-1)+1:3*itmp,indiv)=coor_old(3*(iatomj-1)+1:3*iatomj,iperm(father2))
!          enddo
!          idir=0
!          do jj=1,kk
!            itmp=itmp+1
!            iatom=japerm(natom-jj+1)
!            iatomj=iaperm(natom-jj+1)
!            coor(3*(itmp-1)+1:3*itmp,indiv-1)=coor_old(3*(iatom-1)+1:3*iatom,iperm(father2))
!            coor(3*(itmp-1)+1:3*itmp,indiv)=coor_old(3*(iatomj-1)+1:3*iatomj,iperm(father1))
!          enddo
!          do jj=1,natom
!            coori(:)=coor((jj-1)*3+1:jj*3,indiv)
!            distances(jj,indiv)=sqrt(dot_product(coori,coori))
!            coori(:)=coor((jj-1)*3+1:jj*3,indiv-1)
!            distances(jj,indiv-1)=sqrt(dot_product(coori,coori))
!          enddo
!          call initialize_perm(japerm,natom)
!          call initialize_perm(iaperm,natom)
!          call sort_dp(natom,distances(:,indiv),japerm,tol9)
!          call sort_dp(natom,distances(:,indiv-1),iaperm,tol9)
!          smallacell=minval(acell_old(:,iperm(father1)))
!          itmp=comp_indiv(distances,indiv-1,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv-1),smallacell)
!          if (itmp>0) then
!            coor(:,indiv-1)=coor(:,indiv)
!            acell(:,indiv-1)=acell_old(:,iperm(father2))
!            rprim(:,:,indiv-1)=rprim_old(:,:,iperm(father2))
!            indiv=indiv-1
!          write(187,*) 'Rejected 7'
!           else
!            acell(:,indiv-1)=acell_old(:,iperm(father1))
!            rprim(:,:,indiv-1)=rprim_old(:,:,iperm(father1))
!          write(187,*) checkatomicdist(ndimen,coor(:,indiv),smallacell),smallacell,ndimen,indiv
!          do jj=1,ndimen,3
!            write(187,*) coor(jj:jj+2,indiv)
!          enddo
!         endif
!          smallacell=minval(acell_old(:,iperm(father1)))
!          itmp=comp_indiv(distances,indiv,natom,nimage)+checkatomicdist(ndimen,coor(:,indiv),smallacell)
!          if (itmp>0) then
!             indiv=indiv-1
!          write(187,*) 'Rejected 7 a'
!          else
!             acell(:,indiv)=acell_old(:,iperm(father2))
!            rprim(:,:,indiv)=rprim_old(:,:,iperm(father2))
!          endif
!        endif
     end select
  enddo


 results_img(iimage,itimimage+1)%vel = results_img(iimage,itimimage)%vel 
 do iimage=1,nimage
   results_img(iimage,itimimage+1)%acell = acell(:,iimage)
   results_img(iimage,itimimage+1)%rprim = rprim(:,:,iimage)
   do ii=1,natom
     results_img(iimage,itimimage+1)%xred(1:3,ii) = coor((ii-1)*3+1:ii*3,iimage)
   enddo
 enddo
 ABI_DEALLOCATE(coor)
 ABI_DEALLOCATE(acell)
 ABI_DEALLOCATE(acell_old)
 ABI_DEALLOCATE(rprim)
 ABI_DEALLOCATE(rprim_old)
 ABI_DEALLOCATE(zcoori)
 ABI_DEALLOCATE(zcoorj)
 ABI_DEALLOCATE(coor_old)
 ABI_DEALLOCATE(iperm)
 ABI_DEALLOCATE(japerm)
 ABI_DEALLOCATE(iaperm)
 ABI_DEALLOCATE(etotal_img)
 ABI_DEALLOCATE(fitness)
 ABI_DEALLOCATE(distances)

end subroutine predict_ga

!!*************** local subroutines

INTEGER FUNCTION choosefather(x1,F,n)

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'choosefather'
!End of the abilint section

 implicit none

 integer :: ii
 integer,intent(in) :: n
 real(dp), intent(in) :: x1
 real(dp), dimension(:), intent(in) :: F

 choosefather=1
 do ii=2,n
    if(F(ii-1)<x1.and.x1<=F(ii))then
       choosefather=ii
       exit
    endif
 enddo
end FUNCTION choosefather

SUBROUTINE swap(a,b)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'swap'
!End of the abilint section

 integer, intent(inout) :: a,b
 integer :: dum
   dum=a; a=b;  b=dum
END SUBROUTINE swap

INTEGER FUNCTION comp_indiv(distances,indiv,natom,nimage)

!! comparing individuals from the same generation and check they
!! are not too close.  We compare all individuals with individual: indiv.
!! We assume, all distances for a given individual are ordered and
!! we define a metric from the total difference distances between individuals.
!! if comp_indiv is 0, means that two individuals are two close.

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'comp_indiv'
!End of the abilint section

 implicit none

 integer, intent(in) :: indiv,natom,nimage
 real(dp),intent(in) ::distances(natom,nimage)
 real(dp) :: diff

 integer :: ii,jj
 
 comp_indiv=1

 do ii=1,indiv-1
   diff=0.0_dp
   do jj=1,natom
     diff=(distances(jj,indiv)-distances(jj,ii))**2
   enddo
   if (diff.le.0.001_dp) comp_indiv=0
 enddo

end FUNCTION comp_indiv

SUBROUTINE initialize_perm(iperm,nimage)

!! initialize the vector iperm with corresponding indices.

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initialize_perm'
!End of the abilint section

 implicit none

 integer, intent(in) :: nimage
 integer, intent(inout) :: iperm(nimage)
 integer :: ii


 do ii=1,nimage
   iperm(ii)=ii
 enddo

end SUBROUTINE initialize_perm

! if after a genetic rule, check if two atoms in the same gen
! are too close


INTEGER FUNCTION checkatomicdist(ndim,coord,acell)

!! check if two atoms are two close.
!! if they are, checkatomicdist=0

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'checkatomicdist'
!End of the abilint section

 implicit none

 integer, intent(in) :: ndim
 real(dp), intent(in) :: coord(ndim),acell
 real(dp) :: ci(3),cj(3),cij(3),dist
 integer :: ii,jj

 checkatomicdist=1
 do ii=1,ndim,3
   ci(:)=coord(ii:ii+2)
   do jj=1,ii-3,3
     cj(:)=coord(jj:jj+2)
     cij=ci-cj
     dist=sqrt(dot_product(cij,cij))
!! check if both atoms, ii and jj are too close
     if (dist < 0.001_dp) then
       checkatomicdist=0
       write (187,*) 'DIST=',ii,jj,ci,cj
       EXIT
     endif
!! periodic boundary conditions and closest atoms
     cij=cij-0.5_dp
     cij= cij-ANINT(cij)
     cij=cij+0.5_dp
     dist=acell*sqrt(dot_product(cij,cij))
     if (dist<1.2_dp) then
       checkatomicdist=0
       EXIT
     endif
   enddo
 enddo

END FUNCTION checkatomicdist
 
end MODULE m_use_ga
