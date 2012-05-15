!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawuj_red
!! NAME
!!  pawuj_red
!!
!! FUNCTION
!!  Store atomic occupancies, potential shift, positions in dtpawuj datastructure. 
!!  
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2012 ABINIT group (DJA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! 
!!  fatvshift=factor that multiplies atvshift
!!  natom=number of atoms in cell
!!  ntypat = number of atom types
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawprtvol= printing volume
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type of each atom
!!
!! OUTPUT
!!  dtpawuj(0:ndtpawuj) (initialization of fields vsh, occ, iuj,nnat) 
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      linvmat,prmat,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawuj_red(dtset,dtpawuj,fatvshift,natom,ntypat,paw_ij,pawrad,pawtab,ndtpawuj)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_parameters

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawuj_red'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_66_paw, except_this_one => pawuj_red
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)                :: natom,ntypat,ndtpawuj
 real(dp),intent(in)               :: fatvshift
!arrays
 type(paw_ij_type),intent(in)      :: paw_ij(natom)
 type(pawtab_type),intent(in)      :: pawtab(ntypat)
 type(pawrad_type),intent(in)      :: pawrad(ntypat) 
 type(dataset_type),intent(in)     :: dtset
 type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtpawuj)

!Local variables-------------------------------
!scalars
 integer,parameter           :: natmax=2,ncoeff=3
 integer                     :: iatom,im1,im2,ispden,itypat,ll,nspden,nsppol,iuj 
 integer                     :: natvshift,pawujat,ndtset,typawujat
 logical                     :: usepawu !antiferro,
 character(len=1000)         :: message,hstr
 character(len=500)          :: messg
!arrays
 logical                     :: musk(3,natom)
 integer                     :: nnat,typat(natom),natpawu
 real(dp)                    :: rrtab(ncoeff),wftab(ncoeff),a(ncoeff,ncoeff),b(ncoeff,ncoeff)! ,a(ncoeff,ncoeff)
 real(dp),allocatable        :: nnocctot(:,:) !,magv(:)
 real(dp),allocatable        :: atvshift(:,:,:) ! atvshift(natvshift,2,natom)
 logical,allocatable         :: dmusk(:,:),atvshmusk(:,:,:) !atvshmusk(natvshift,2,natom)

! *********************************************************************

!Initializations

 nspden=paw_ij(1)%nspden ; nsppol=paw_ij(1)%nsppol
 natvshift=dtset%natvshift
 pawujat=dtset%pawujat
 natpawu=dtset%natpawu   ; ndtset=dtset%ndtset
 typat=dtset%typat
 ABI_ALLOCATE(atvshift,(natvshift,nspden,natpawu))
 ABI_ALLOCATE(atvshmusk,(natvshift,nspden,natpawu))
 ABI_ALLOCATE(dmusk,(nspden,natom))
 musk=.false.; dmusk=.false.
 atvshift=fatvshift*dtset%atvshift
 typawujat=typat(pawujat)
 usepawu=(count(pawtab(:)%usepawu>0)>0)

 nnat=0
 if (usepawu) then
   ABI_ALLOCATE(nnocctot,(nspden,natom))
   nnocctot=0
   write(message,'(3a)') ch10, '---------- pawuj_red ------ ',ch10
   call wrtout(std_out,  message,'COLL'); 
   do iatom=1,natom
     itypat=typat(iatom);ll=pawtab(itypat)%lpawu
     if ((ll>=0).and.(pawtab(itypat)%usepawu>0).and.itypat==typawujat) then
       musk(:,iatom)=(/.true., .true., .true. /)
       atvshmusk(:,:,iatom)=reshape((/ (( (im1==1), im1=1,natvshift)  ,im2=1,nspden ) /),(/natvshift,nspden/))
       do ispden=1,nspden
         nnocctot(ispden,iatom)=paw_ij(iatom)%nocctot(ispden)
         dmusk(ispden,iatom)=.true.
       end do
       nnat=nnat+1
     end if
   end do


   iuj=maxval(dtpawuj(:)%iuj)
!  DEBUG
!  write(std_out,*)' pawuj_red: iuj',iuj
!  END DEBUG

!  DEBUG
!  write(std_out,*)' pawuj_red: dtpawuj(:)%iuj ',dtpawuj(:)%iuj
!  END DEBUG

   if (iuj==1.or.iuj==3) then  ! 1 and 3: non-scf steps
     dtpawuj(iuj+1)%iuj=iuj+1
   end if

!  TODO: check that this is correct: this point is passed several times for a given value of iuj - should the stuff be accumulated instead of replaced?
   if(associated(dtpawuj(iuj)%vsh))  then
     ABI_DEALLOCATE(dtpawuj(iuj)%vsh)
   end if
   ABI_ALLOCATE(dtpawuj(iuj)%vsh,(nspden,nnat))
   if(associated(dtpawuj(iuj)%occ))  then
     ABI_DEALLOCATE(dtpawuj(iuj)%occ)
   end if
   ABI_ALLOCATE(dtpawuj(iuj)%occ,(nspden,nnat))
   if(associated(dtpawuj(iuj)%xred))  then
     ABI_DEALLOCATE(dtpawuj(iuj)%xred)
   end if
   ABI_ALLOCATE(dtpawuj(iuj)%xred,(3,nnat))
   
   if (iuj==1) then
     ABI_ALLOCATE(dtpawuj(0)%vsh,(nspden,nnat))
     ABI_ALLOCATE(dtpawuj(0)%occ,(nspden,nnat))
     ABI_ALLOCATE(dtpawuj(0)%xred,(3,nnat))
     dtpawuj(0)%vsh=0
     dtpawuj(0)%occ=0 
     dtpawuj(0)%xred=0 
   end if

   rrtab=(/0.75_dp,0.815_dp,1.0_dp/)*pawtab(typawujat)%rpaw 
   wftab=pawtab(typawujat)%phi(pawrad(typawujat)%mesh_size,pawtab(typawujat)%lnproju(1))

!  DEBUG
!  write(std_out,*)' pawuj_red: rrtab ',rrtab
!  write(std_out,*)' pawuj_red: wftab ',wftab
!  END DEBUG

   do im1=1,ncoeff
!    DEBUG
!    write(std_out,*)' pawuj_red: ncoeff ',ncoeff,' im1 ',im1
!    END DEBUG    
     if (pawrad(typawujat)%mesh_type==1) then
       im2=nint(rrtab(im1)/pawrad(typawujat)%rstep+1)
     else if (pawrad(typawujat)%mesh_type==2) then
       im2=nint(log(rrtab(im1)/pawrad(typawujat)%rstep+1)/pawrad(typawujat)%lstep+1)
     else if (pawrad(typawujat)%mesh_type==3) then
       im2=nint(log(rrtab(im1)/pawrad(typawujat)%rstep)/pawrad(typawujat)%lstep+1)
     else if (pawrad(typawujat)%mesh_type==4) then
       im2=nint(pawrad(typawujat)%mesh_size*(1-exp((-one)*rrtab(im1)/pawrad(typawujat)%rstep))+1)
     end if

!    DEBUG
!    write(std_out,*)' pawuj_red: im2 ',im2
!    END DEBUG

     rrtab(im1)=pawrad(typawujat)%rad(im2)
     wftab(im1)=pawtab(typawujat)%phi(im2,pawtab(typawujat)%lnproju(1))
   end do
   write(message,fmt='(a,i3,a,10f10.5)')' pawuj_red: mesh_type',pawrad(typawujat)%mesh_type,' rrtab:', rrtab 
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,10f10.5)')' pawuj_red: wftab', wftab
   call wrtout(std_out,message,'COLL')
   a=reshape((/ (( rrtab(im2)**(im1-1), im1=1,3)  ,im2=1,3 )/),(/ncoeff,ncoeff/))
   write(messg,fmt='(a)')'A'
   call linvmat(a,b,ncoeff,messg,2,0.0_dp,3) ! linvmat(inmat,oumat,nat,nam,option,gam,prtvol)
   write(std_out,*) 'pawuj_red: a,b ', a,b
   wftab=matmul(wftab,b)
   write(std_out,*) 'pawuj_red: wftab ', wftab
   dtpawuj(iuj)%wfchr(4:6)=wftab 
   
   dtpawuj(iuj)%nat=nnat
   write(std_out,*) 'pawuj_red: m1' 
   dtpawuj(iuj)%vsh=reshape(pack(atvshift,atvshmusk),(/ nspden,nnat /))
!  factor in next line to compensate nocctot contains just occ of 1 spin channel for nspden=1
   write(std_out,*) 'pawuj_red: m2'
   dtpawuj(iuj)%occ=reshape(pack(nnocctot,dmusk),(/nspden,nnat/))*(3-nspden)
   write(std_out,*) 'pawuj_red: m3'
!  dtpawuj(iuj)%occ=dtpawuj(iuj)%occ/pawtab(typawujat)%ph0phiint(1)
   
   write(std_out,*) 'pawuj_red: occ ', dtpawuj(iuj)%occ
   
   dtpawuj(iuj)%xred=reshape(pack(dtset%xred_orig(:,:,1),musk),(/3,nnat/))
   dtpawuj(iuj)%ph0phiint=pawtab(typawujat)%ph0phiint(1)
   dtpawuj(iuj)%wfchr(1:3)=(/ pawtab(typawujat)%zioneff(1)*(dtset%lpawu(typawujat)+2),&
&   one*(dtset%lpawu(typawujat)+1),one*(dtset%lpawu(typawujat))/)
   dtpawuj(iuj)%pawrad=pawtab(typawujat)%rpaw

   write(std_out,*) 'pawuj_red: wfchr ',dtpawuj(iuj)%wfchr


   write (hstr,'(I0)') iuj 
   write(message,'(a,a,I3,I3,a)') ch10, '---------- MARK ------ ',iuj,maxval(dtpawuj(:)%iuj) ,ch10
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a)') 'vsh'//trim(hstr)
   call wrtout(std_out,message,'COLL')
   call prmat(dtpawuj(iuj)%vsh(:,:),1,nnat*nspden,1)
   write(message,fmt='(a)') 'occ'//trim(hstr)
   call wrtout(std_out,message,'COLL')
   call prmat(dtpawuj(iuj)%occ(:,:),1,nnat*nspden,1)
   write(message, '(3a)' )'---------- MARK ---------- ',ch10
   call wrtout(std_out,message,'COLL') 
 end if !usepawu


 ABI_DEALLOCATE(nnocctot)
 ABI_DEALLOCATE(dmusk)
 ABI_DEALLOCATE(atvshift)
 ABI_DEALLOCATE(atvshmusk)

end subroutine pawuj_red

!!***
