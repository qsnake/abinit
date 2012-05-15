!{\src2tex{textfont=tt}}
!!****f* ABINIT/outvar1
!! NAME
!! outvar1
!!
!! FUNCTION
!! Echo variables between acell and natom (by alphabetic order)
!! for the ABINIT code.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  choice= 1 if echo of preprocessed variables, 2 if echo after call driver
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for echoed output
!!  istatr=repetition rate for status file
!!  istatshft=shift of the repetition rate for status file
!!  jdtset_(0:ndtset_alloc)=actual index of the dataset (equal to dtsets(:)%jdtset)
!!  mxvals=maximun size of some arrays along all datasets, including
!!         mxlpawu      =maximal value of input lpawu for all the datasets
!!         mxgw_nqlwl   =maximal value of input gw_nqlwl for all the datasets
!!         mxmband      =maximum number of bands
!!         mxnatom      =maximal value of input natom for all the datasets
!!         mxnatpawu    =maximal value of number of atoms on which +U is applied for all the datasets
!!         mxnatsph     =maximal value of input natsph for all the datasets
!!         mxnatvshift  =maximal value of input natvshift for all the datasets
!!         mxnconeq     =maximal value of input nconeq for all the datasets
!!         mxnimage     =maximal value of input nimage for all the datasets
!!         mxnkptgw     =maximal value of input nkptgw for all the datasets
!!         mxnkpt       =maximal value of input nkpt for all the datasets
!!         mxnnos       =maximal value of input nnos for all the datasets
!!         mxnqptdm     =maximal value of input nqptdm for all the datasets
!!         mxnspinor    =maximal value of input nspinor for all the datasets
!!         mxnsppol     =maximal value of input nsppol for all the datasets
!!         mxnsym       =maximum number of symmetries
!!         mxntypat     =maximum number of type of atoms
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set. Use for most dimensioned arrays.
!!  npsp=number of pseudopotentials
!!  nqptdm=the number of q vectors provided by the user to calculate DM in GW
!!  prtvol_glob= if 0, minimal output volume, if 1, no restriction.
!!  response= 1 if response variables must be output, 0 otherwise.
!!  response_(0:ndtset_alloc)= 1 if response variables must be output, 0 otherwise,
!!   for different datasets
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Note that this routine is called only by the processor me==0 .
!! In consequence, no use of message and wrtout routine.
!! The lines of code needed to output the defaults are preserved
!! (see last section of the routine, but are presently disabled)
!!
!!  Note that acell, occ, rprim, xred and vel might have been modified by the
!!  computation, so that their values if choice=1 or choice=2 will differ.
!!
!! PARENTS
!!      outvars
!!
!! CHILDREN
!!      prttagm,prttagm_images
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outvar1 (choice,dtsets,iout,istatr,istatshft,&
& jdtset_,mxvals, ncid,ndtset,ndtset_alloc,npsp,prtvol_glob,response,&
& response_,results_out,usepaw)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_results_out

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outvar1'
 use interfaces_57_iovars, except_this_one => outvar1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,iout,istatr,istatshft,ndtset
 integer,intent(in) :: ndtset_alloc,npsp,prtvol_glob,response,usepaw,ncid
!arrays
 integer,intent(in) :: jdtset_(0:ndtset_alloc),response_(ndtset_alloc)
 type(ab_maxvals),intent(in) :: mxvals
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
 character(len=*), parameter :: format01110 ="(1x,a16,1x,(t22,8i8) )"
 character(len=*), parameter :: format01150 ="(1x,a16,1x,(t22,3es16.8))"
 character(len=*), parameter :: format01150a="(1x,a16,a,1x,(t22,3es16.8))"
 character(len=*), parameter :: format01155 ="(1x,a16,1x,(t22,10i5))"
 character(len=*), parameter :: format01155a="(1x,a16,a,1x,(t22,10i5))"
 character(len=*), parameter :: format01160 ="(1x,a16,1x,(t22,3es18.10)) "
 character(len=*), parameter :: format01160a="(1x,a16,a,1x,(t22,3es18.10)) "
 character(len=*), parameter :: format01170 ="(1x,a16,1x,(t22,5f11.6)) "
 character(len=*), parameter :: format01170a="(1x,a16,a,1x,(t22,5f11.6)) "
!scalars
 integer,parameter :: nkpt_max=50
 integer :: allowed,first,iatom,idtset,ii,iimage,ikpt,istatr_defo
 integer :: istatshft_defo,kptopt,marr,narr!,jdtset
 integer :: multi_natom,multi_natfix,multi_natfixx
 integer :: multi_natfixy,multi_natfixz,multi_natsph,multi_natvshift,multi_nberry,multi_nimage,multi_nkpt
 integer :: multi_nkptgw,multi_nqptdm,multi_nshiftk,multi_ntypalch
 integer :: multi_nsppol,multi_gw_nqlwl,multi_atsph
 integer :: multi_ntypat,natfix,natfixx,natfixy,natfixz,natom,natsph
 integer :: ndtset_kptopt,nimage,nqpt,nkpt_eff
 integer :: nshiftk,ntypalch,ntypat,tnkpt
 logical :: compute_static_images
 real(dp) :: cpus,kpoint
! character(len=4) :: appen
 character(len=4) :: stringimage
!arrays
 integer,allocatable :: iatfixio_(:,:),iatfixx_(:,:),iatfixy_(:,:)
 integer,allocatable :: iatfixz_(:,:),intarr(:,:),istwfk_2(:,:)
 integer,allocatable :: jdtset_kptopt(:),natfix_(:),natfixx_(:),natfixy_(:)
 integer,allocatable :: natfixz_(:)
 integer,allocatable :: narrm(:)
 integer,allocatable :: nimagem(:),prtimg(:,:)
 real(dp),allocatable :: dprarr(:,:),dprarr_images(:,:,:)
 character(len=8),allocatable :: strimg(:)

! *************************************************************************

!###########################################################
!### 01. Initial allocations

!DEBUG
!write(std_out,*)' outvar1 : enter '
!ENDDEBUG
!
 ABI_ALLOCATE(narrm,(0:ndtset_alloc))
 ABI_ALLOCATE(nimagem,(0:ndtset_alloc))

!Must treat separately the translation of iatfix from the internal
!representation to the input/output representation
 ABI_ALLOCATE(natfix_,(0:ndtset_alloc))
 ABI_ALLOCATE(iatfixio_,(mxvals%mxnatom,0:ndtset_alloc))
 ABI_ALLOCATE(natfixx_,(0:ndtset_alloc))
 ABI_ALLOCATE(iatfixx_,(mxvals%mxnatom,0:ndtset_alloc))
 ABI_ALLOCATE(natfixy_,(0:ndtset_alloc))
 ABI_ALLOCATE(iatfixy_,(mxvals%mxnatom,0:ndtset_alloc))
 ABI_ALLOCATE(natfixz_,(0:ndtset_alloc))
 ABI_ALLOCATE(iatfixz_,(mxvals%mxnatom,0:ndtset_alloc))
 natfix_(0:ndtset_alloc)=0 ; iatfixio_(:,0:ndtset_alloc)=0
 natfixx_(0:ndtset_alloc)=0 ; iatfixx_(:,0:ndtset_alloc)=0
 natfixy_(0:ndtset_alloc)=0 ; iatfixy_(:,0:ndtset_alloc)=0
 natfixz_(0:ndtset_alloc)=0 ; iatfixz_(:,0:ndtset_alloc)=0
 do idtset=1,ndtset_alloc
!  DEBUG
!  write(std_out,*)' outvar1 : iatfix_ for idtset= ',idtset
!  ENDDEBUG
   do iatom=1,dtsets(idtset)%natom
!    First look whether the atom is fixed along the three directions
     if( dtsets(idtset)%iatfix(1,iatom)+ &
&     dtsets(idtset)%iatfix(2,iatom)+ &
&     dtsets(idtset)%iatfix(3,iatom)   ==3 )then
       natfix_(idtset)=natfix_(idtset)+1
!      DEBUG
!      write(std_out,*)' outvar1: iatom,natfix_(idtset)=',iatom,natfix_(idtset)
!      ENDDEBUG
       iatfixio_(natfix_(idtset),idtset)=iatom
     else
!      Now examine each direction, one at a time
       if( dtsets(idtset)%iatfix(1,iatom) ==1)then
         natfixx_(idtset)=natfixx_(idtset)+1
         iatfixx_(natfixx_(idtset),idtset)=iatom
       end if
       if( dtsets(idtset)%iatfix(2,iatom) ==1)then
         natfixy_(idtset)=natfixy_(idtset)+1
         iatfixy_(natfixy_(idtset),idtset)=iatom
       end if
       if( dtsets(idtset)%iatfix(3,iatom) ==1)then
         natfixz_(idtset)=natfixz_(idtset)+1
         iatfixz_(natfixz_(idtset),idtset)=iatom
       end if
     end if
   end do
!  DEBUG
!  write(std_out,*)' natfix ...'
!  write(std_out,*)natfix_(idtset),natfixx_(idtset),natfixy_(idtset),natfixz_(idtset)
!  ENDDEBUG
 end do

 ABI_ALLOCATE(strimg,(mxvals%mxnimage))
 ABI_ALLOCATE(prtimg,(mxvals%mxnimage,0:ndtset_alloc))
 do iimage=1,mxvals%mxnimage
   if(iimage<10)then
     write(stringimage,'(i1)')iimage
   else if(iimage<100)then
     write(stringimage,'(i2)')iimage
   else if(iimage<1000)then
     write(stringimage,'(i3)')iimage
   else if(iimage<10000)then
     write(stringimage,'(i4)')iimage
   end if
   strimg(iimage)='_'//trim(stringimage)//'img'
 end do
 strimg(1)=''

!Maximal size of dprarr and intarr arrays
 marr=max(3*mxvals%mxnatom,3*mxvals%mxnkptgw,mxvals%mxnkpt*mxvals%mxnsppol*mxvals%mxmband,&
& 3*mxvals%mxnkpt,npsp,mxvals%mxntypat,3*mxvals%mxnqptdm,3*mxvals%mxgw_nqlwl,&
& 9*mxvals%mxnsym,mxvals%mxnatsph,mxvals%mxnimage,mxvals%mxnatvshift*mxvals%mxnsppol*mxvals%mxnatom)
 ABI_ALLOCATE(dprarr,(marr,0:ndtset_alloc))
 ABI_ALLOCATE(intarr,(marr,0:ndtset_alloc))
 ABI_ALLOCATE(dprarr_images,(marr,mxvals%mxnimage,0:ndtset_alloc))

!###########################################################
!### 02. Set up dimensions : determine whether these are
!##      different for different datasets. (MULTI)

 multi_gw_nqlwl=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%gw_nqlwl/=dtsets(idtset)%gw_nqlwl)multi_gw_nqlwl=1
   end do
 end if

 multi_nqptdm=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nqptdm/=dtsets(idtset)%nqptdm)multi_nqptdm=1
   end do
 end if

 multi_natfix=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfix_(1)/=natfix_(idtset))multi_natfix=1
   end do
 end if
 if(multi_natfix==0)natfix=natfix_(1)

 multi_natfixx=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfixx_(1)/=natfixx_(idtset))multi_natfixx=1
   end do
 end if
 if(multi_natfixx==0)natfixx=natfixx_(1)

 multi_natfixy=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfixy_(1)/=natfixy_(idtset))multi_natfixy=1
   end do
 end if
 if(multi_natfixy==0)natfixy=natfixy_(1)

 multi_natfixz=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(natfixz_(1)/=natfixz_(idtset))multi_natfixz=1
   end do
 end if
 if(multi_natfixz==0)natfixz=natfixz_(1)

 multi_natom=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%natom/=dtsets(idtset)%natom)multi_natom=1
   end do
 end if
 if(multi_natom==0)natom=dtsets(1)%natom

 multi_natsph=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%natsph/=dtsets(idtset)%natsph)multi_natsph=1
   end do
 end if
 if(multi_natsph==0)natsph=dtsets(1)%natsph

 multi_natvshift=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%natvshift/=dtsets(idtset)%natvshift)multi_natvshift=1
   end do
 end if


 multi_nberry=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nberry/=dtsets(idtset)%nberry)multi_nberry=1
   end do
 end if

 multi_nimage=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nimage/=dtsets(idtset)%nimage)multi_nimage=1
   end do
 end if
 if(multi_nimage==0)nimage=dtsets(1)%nimage

 multi_nkptgw=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nkptgw/=dtsets(idtset)%nkptgw)multi_nkptgw=1
   end do
 end if

 multi_nkpt=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nkpt/=dtsets(idtset)%nkpt)multi_nkpt=1
   end do
 end if

 multi_nshiftk=0
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
   first=0
   do idtset=1,ndtset_alloc
     kptopt=dtsets(idtset)%kptopt
     if(kptopt>=1)then
       if(first==0)then
         first=1
         nshiftk=dtsets(idtset)%nshiftk
       else
         if(nshiftk/=dtsets(idtset)%nshiftk)multi_nshiftk=1
       end if
     end if
   end do
 end if

 multi_nsppol=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nsppol/=dtsets(idtset)%nsppol)multi_nsppol=1
   end do
 end if

 multi_ntypalch=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%ntypalch/=dtsets(idtset)%ntypalch)multi_ntypalch=1
   end do
 end if
 if(multi_ntypalch==0)ntypalch=dtsets(1)%ntypalch

 multi_ntypat=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%ntypat/=dtsets(idtset)%ntypat)multi_ntypat=1
   end do
 end if
 if(multi_ntypat==0)ntypat=dtsets(1)%ntypat

!###########################################################
!### 03. Print all the input variables (From A to N)
!##

!*****************************************
!Templates using generalized prttagm
!
!1. Case without reshape
!
!!VARIABLE
!intarr(1:marr,0)=0             ! default value
!narr=SIZE                      ! default size for all datasets
!do idtset=0,ndtset_alloc       ! especific size for each dataset
!narrm(idtset)=dtsets(idtset)%SIZE
!if (narrm(idtset)>0) then
!intarr(1:narrm(idtset),idtset)=dtsets(idtset)%VARIABLE(1:narrm(idtset))
!end if
!end do
!call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
!& narrm,ncid,ndtset_alloc,'VARIABLE','INT',MULTI)
!
!2. Case with reshape
!
!!VARIBLE
!narr=dtsets(1)%SIZE1*dtsets(1)%SIZE2*... ! default size for all datasets
!do idtset=0,ndtset_alloc       ! especific size for each dataset
!narrm(idtset)=dtsets(idtset)%SIZE1*dtsets(idtset)%SIZE*...
!if (narrm(idtset)>0) then
!dprarr(1:narrm(idtset),idtset)=&
!& reshape(dtsets(idtset)%VARIABLE(1:dtsets(idtset)%SIZE1,&
!& 1:dtsets(idtset)%SIZE2,  ....     ),&
!& (/ narrm(idtset) /) )
!end if
!end do
!call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
!& narrm,ncid,ndtset_alloc,'VARIABLE','DPR',&
!& MULTI)


!###########################################################
!### 03. Print all the input variables (A)
!##

 intarr(1,:)=dtsets(:)%accesswff
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'accesswff','INT',0)

!acell
 if(mxvals%mxnimage==1)then
   do idtset=0,ndtset_alloc
     dprarr(1:3,idtset)=results_out(idtset)%acell(:,1)
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'acell','LEN',0)
 else
   prtimg(:,:)=1
   do idtset=1,ndtset_alloc       ! especific size for each dataset
     nimagem(idtset)=dtsets(idtset)%nimage
     narrm(idtset)=3
     do iimage=1,dtsets(idtset)%nimage
       if (narrm(idtset)>0) then
         dprarr_images(1:narrm(idtset),iimage,idtset)=&
&         results_out(idtset)%acell(1:3,iimage)
       end if
     end do
   end do
   call prttagm_images(dprarr_images,iout,jdtset_,&
&   marr,narrm,ncid,ndtset_alloc,'acell',&
&   mxvals%mxnimage,nimagem,ndtset,prtimg,strimg)

 end if

!algalch
 narr=ntypalch                      ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypalch
   intarr(1:narrm(idtset),idtset)=dtsets(idtset)%algalch(1:narrm(idtset))
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'algalch','INT',multi_ntypalch)

!atvshift
 if(usepaw>0)then
   narr=dtsets(1)%natvshift*dtsets(1)%nsppol*mxvals%mxnatpawu ! default size for all datasets
   do idtset=0,ndtset_alloc       ! especific size for each dataset
     narrm(idtset)=dtsets(idtset)%natvshift*dtsets(idtset)%nsppol*mxvals%mxnatpawu
     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%atvshift(1:dtsets(idtset)%natvshift,&
&     1:dtsets(idtset)%nsppol,1:mxvals%mxnatpawu),&
&     (/ narrm(idtset) /) )
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,5,marr,narr,&
&   narrm,ncid,ndtset_alloc,'atvshift','DPR',&
&   multi_natvshift+multi_nsppol+multi_natom)

 end if

!amu
 narr=ntypat                      ! default size for all datasets
 dprarr(1:narr,0:ndtset_alloc)=zero
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypat
   dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%amu(1:narrm(idtset))
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'amu','DPR',multi_ntypat)

 intarr(1,:)=dtsets(:)%awtr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'awtr','INT',0)

!###########################################################
!### 03. Print all the input variables (B)
!##

 intarr(1,:)=dtsets(:)%bandpp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bandpp','INT',0)

 intarr(1,:)=dtsets(:)%berryopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'berryopt','INT',0)

 intarr(1,:)=dtsets(:)%berrystep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'berrystep','INT',0)

 intarr(1,:)=dtsets(:)%bdberry(1)
 intarr(2,:)=dtsets(:)%bdberry(2)
 intarr(3,:)=dtsets(:)%bdberry(3)
 intarr(4,:)=dtsets(:)%bdberry(4)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,4,narrm,ncid,ndtset_alloc,'bdberry','INT',0)

 if(response==1) then
   intarr(1,:)=dtsets(:)%bdeigrf
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bdeigrf','INT',0)
 end if

!bdgw
 narr=2*dtsets(1)%nkptgw*dtsets(1)%nsppol ! default size for all datasets
 do idtset=1,ndtset_alloc        ! especific size for each dataset
   narrm(idtset)=2*dtsets(idtset)%nkptgw*dtsets(idtset)%nsppol
   if (narrm(idtset)>0)&
&   intarr(1:narrm(idtset),idtset)=&
&   reshape(dtsets(idtset)%bdgw(1:2,1:dtsets(idtset)%nkptgw,1:dtsets(idtset)%nsppol),&
&   (/ narrm(idtset) /) )
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,&
& narrm,ncid,ndtset_alloc,'bdgw','INT',multi_nkptgw)

 dprarr(1,:)=dtsets(:)%bfield(1)
 dprarr(2,:)=dtsets(:)%bfield(2)
 dprarr(3,:)=dtsets(:)%bfield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'bfield','DPR',0)

 dprarr(1,:)=dtsets(:)%bmass
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'bmass','DPR',0)

 dprarr(1,:)=dtsets(:)%boxcenter(1)
 dprarr(2,:)=dtsets(:)%boxcenter(2)
 dprarr(3,:)=dtsets(:)%boxcenter(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'boxcenter','DPR',0)

 dprarr(1,:)=dtsets(:)%boxcutmin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'boxcutmin','DPR',0)

 if (usepaw==1) then

   dprarr(1,:)=dtsets(:)%bxctmindg
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'bxctmindg','DPR',0)

 end if

!###########################################################
!### 03. Print all the input variables (C)
!##

 if (ANY(dtsets(:)%cd_custom_imfrqs/=0)) then
   intarr(1,:)=dtsets(:)%cd_custom_imfrqs
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cd_custom_imfrqs','INT',0)

!  do idtset=1,ndtset_alloc
!  dprarr(1:dtsets(idtset)%cd_custom_imfrqs,idtset)=dtsets(idtset)%cd_imfrqs
!  end do
!  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,MAXVAL(dtsets(:)%cd_custom_imfrqs),narrm,ncid,&
!  &   ndtset_alloc,'cd_imfrqs','ENE',0)
 end if

 dprarr(1,:)=dtsets(:)%cd_halfway_freq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cd_halfway_freq','ENE',0)

 dprarr(1,:)=dtsets(:)%cd_max_freq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cd_max_freq','ENE',0)

 intarr(1,:)=dtsets(:)%cd_use_tangrid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cd_use_tangrid','INT',0)

 intarr(1,:)=dtsets(:)%gw_reconst_scr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_reconst_scr','INT',0)
 
 intarr(1,:)=dtsets(:)%gw_use_pole_scr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_use_pole_scr','INT',0)

 intarr(1,:)=dtsets(:)%gw_npoles
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_npoles','INT',0)

 intarr(1,:)=dtsets(:)%cd_full_grid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'cd_full_grid','INT',0)

 if (ANY(dtsets(:)%cd_subset_freq(1)/=0)) then
   intarr(1,:)=dtsets(:)%cd_subset_freq(1)
   intarr(2,:)=dtsets(:)%cd_subset_freq(2)
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'cd_subset_freq','INT',0)
 end if

 dprarr(1,:)=dtsets(:)%charge
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'charge','DPR',0)

 intarr(1,:)=dtsets(:)%chkexit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'chkexit','INT',0)

 intarr(1,:)=dtsets(:)%chkgwcomp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'chkgwcomp','INT',0)

 intarr(1,:)=dtsets(:)%chksymbreak
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'chksymbreak','INT',0)

 if(dtsets(1)%cpus>one)then
   cpus=dtsets(1)%cpus
   write(iout,'(1x,a16,1x,1p,t22,g10.2,t25,a)') 'cpus',cpus,'(seconds)'
   write(iout,'(1x,a16,1x,1p,t22,g10.2,t25,a)') 'cpum',cpus/60.0_dp,'(minutes)'
   write(iout,'(1x,a16,1x,1p,t22,g10.2,t25,a)') 'cpuh',cpus/3600.0_dp,'(hours)'
 end if

 do idtset=0, ndtset_alloc
   do ii = 1, ntypat
     dprarr(ii,idtset) = dtsets(idtset)%corecs(ii)
   end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,narrm,ncid,ndtset_alloc,'corecs','DPR',0)

!###########################################################
!### 03. Print all the input variables (D)
!##

 intarr(1,:)=dtsets(:)%delayperm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'delayperm','INT',0)

!densty
 narr=ntypat                     ! default size for all datasets
 do idtset=1,ndtset_alloc        ! especific size for each dataset
   narrm(idtset)=dtsets(idtset)%ntypat
!  Only one component of densty is used until now
   dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%densty(1:narrm(idtset),1)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'densty','DPR',multi_ntypat)

 dprarr(1,:)=dtsets(:)%diecut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'diecut','ENE',0)

 dprarr(1,:)=dtsets(:)%diegap
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'diegap','ENE',0)

 dprarr(1,:)=dtsets(:)%dielam
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dielam','DPR',0)

 dprarr(1,:)=dtsets(:)%dielng
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dielng','LEN',0)

 dprarr(1,:)=dtsets(:)%diemac
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'diemac','DPR',0)

 dprarr(1,:)=dtsets(:)%diemix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'diemix','DPR',0)

 if (any(dtsets(1:ndtset_alloc)%diemixmag/=dtsets(1:ndtset_alloc)%diemix)) then
   dprarr(1,:)=dtsets(:)%diemixmag
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'diemixmag','DPR',0)
 end if

 intarr(1,:)=dtsets(:)%diismemory
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'diismemory','INT',0)

 dprarr(1,:)=dtsets(:)%dilatmx
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dilatmx','DPR',0)

 dprarr(1,:)=dtsets(:)%dosdeltae
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dosdeltae','ENE',0)

 dprarr(1,:)=dtsets(:)%dtion
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'dtion','DPR',0)

!dynimage
 intarr(1:marr,0)=1                 ! default value
 narr=nimage                        ! default size for all datasets
 do idtset=1,ndtset_alloc           ! especific size and array for each dataset
   narrm(idtset)=dtsets(idtset)%nimage
   intarr(1:narrm(idtset),idtset)=dtsets(idtset)%dynimage(1:narrm(idtset))
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'dynimage','INT',multi_nimage)


!###########################################################
!### 03. Print all the input variables (E)
!##

 dprarr(1,:)=dtsets(:)%ecut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ecut','ENE',0)

 dprarr(1,:)=dtsets(:)%ecuteps
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ecuteps','ENE',0)

 dprarr(1,:)=dtsets(:)%ecutsigx
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ecutsigx','ENE',0)

 dprarr(1,:)=dtsets(:)%ecutwfn
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ecutwfn','ENE',0)

 dprarr(1,:)=dtsets(:)%ecutsm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ecutsm','ENE',0)

 dprarr(1,:)=dtsets(:)%effmass
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'effmass','DPR',0)

 dprarr(1,:)=dtsets(:)%efield(1)
 dprarr(2,:)=dtsets(:)%efield(2)
 dprarr(3,:)=dtsets(:)%efield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'efield','DPR',0)

 dprarr(1,:)=dtsets(:)%elph2_imagden
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'elph2_imagden','ENE',0)

 intarr(1,:)=dtsets(:)%enunit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'enunit','INT',0)

 dprarr(1,:)=dtsets(:)%eshift
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'eshift','ENE',0)

 dprarr(1,:)=dtsets(:)%esmear
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'esmear','ENE',0)

 dprarr(1,:)=dtsets(:)%exchmix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'exchmix','DPR',0)

!etotal
 if(choice==2)then

   prtimg(:,:)=1
   do idtset=1,ndtset_alloc       ! especific size for each dataset
     compute_static_images=(dtsets(idtset)%istatimg>0)
     nimagem(idtset)=dtsets(idtset)%nimage
     narrm(idtset)=1

     if(dtsets(idtset)%iscf>0 .or. dtsets(idtset)%iscf==-3)then
       do iimage=1,dtsets(idtset)%nimage
         if (narrm(idtset)>0) then
           dprarr_images(1:narrm(idtset),iimage,idtset)=&
&           results_out(idtset)%etotal(iimage)
         end if
         if(.not.(dtsets(idtset)%dynimage(iimage)==1.or.compute_static_images))then
           prtimg(iimage,idtset)=0
         end if
       end do
     else
       narrm(idtset)=0
     end if

   end do
   call prttagm_images(dprarr_images,iout,jdtset_,&
&   marr,narrm,ncid,ndtset_alloc,'etotal',&
&   mxvals%mxnimage,nimagem,ndtset,prtimg,strimg)

 end if

 intarr(1,:)=dtsets(:)%exchn2n3d
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'exchn2n3d','INT',0)

!###########################################################
!### 03. Print all the input variables (F)
!##

 intarr(1,:)=dtsets(:)%pawfatbnd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawfatbnd','INT',0)

 dprarr(1,:)=dtsets(:)%fermie_nest
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'fermie_nest','DPR',0)

 intarr(1,:)=dtsets(:)%ngfft(7)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'fftalg','INT',0)

 intarr(1,:)=dtsets(:)%ngfft(8)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'fftcache','INT',0)

 intarr(1,:)=dtsets(:)%fftgw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'fftgw','INT',0)

 intarr(1,:)=dtsets(:)%fft_opt_lob
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'fft_opt_lob','INT',0)

!force
 if(choice==2)then
   prtimg(:,:)=1
   do idtset=1,ndtset_alloc       ! especific size for each dataset
     compute_static_images=(dtsets(idtset)%istatimg>0)
     nimagem(idtset)=dtsets(idtset)%nimage
     narrm(idtset)=3*dtsets(idtset)%natom

     if(dtsets(idtset)%iscf>0)then
       do iimage=1,dtsets(idtset)%nimage
         if (narrm(idtset)>0) then
           dprarr_images(1:narrm(idtset),iimage,idtset)=&
&           reshape(results_out(idtset)%fcart(1:3,1:dtsets(idtset)%natom,iimage),&
&           (/ narrm(idtset) /) )
         end if
         if(.not.(dtsets(idtset)%dynimage(iimage)==1.or.compute_static_images))then
           prtimg(iimage,idtset)=0
         end if
       end do
     else
       narrm(idtset)=0
     end if

   end do
   call prttagm_images(dprarr_images,iout,jdtset_,&
&   marr,narrm,ncid,ndtset_alloc,'fcart',&
&   mxvals%mxnimage,nimagem,ndtset,prtimg,strimg)

 end if

 dprarr(1,:)=dtsets(:)%fixmom
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'fixmom','DPR',0)

 dprarr(1,:)=dtsets(:)%fxcartfactor
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'fxcartfactor','DPR',0)

 dprarr(1,:)=dtsets(:)%rhoqpmix
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'rhoqpmix','DPR',0)

 dprarr(1,:)=dtsets(:)%freqremin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqremin','ENE',0)

 dprarr(1,:)=dtsets(:)%freqremax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqremax','ENE',0)

 dprarr(1,:)=dtsets(:)%freqspmin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqspmin','ENE',0)

 dprarr(1,:)=dtsets(:)%freqspmax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqspmax','ENE',0)

 dprarr(1,:)=dtsets(:)%freqsusin
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqsusin','DPR',0)

 dprarr(1,:)=dtsets(:)%freqsuslo
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'freqsuslo','DPR',0)

 dprarr(1,:)=dtsets(:)%friction
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'friction','DPR',0)

 intarr(1,:)=dtsets(:)%frzfermi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'frzfermi','INT',0)

!###########################################################
!### 03. Print all the input variables (G)
!##

 intarr(1,:)=dtsets(:)%getbseig
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getbseig','INT',0)

 intarr(1,:)=dtsets(:)%getbsreso
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getbsreso','INT',0)

 intarr(1,:)=dtsets(:)%getbscoup
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getbscoup','INT',0)

 intarr(1,:)=dtsets(:)%gethaydock
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gethaydock','INT',0)

 intarr(1,:)=dtsets(:)%getcell
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getcell','INT',0)

 intarr(1,:)=dtsets(:)%getddk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getddk','INT',0)

 intarr(1,:)=dtsets(:)%getden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getden','INT',0)

 intarr(1,:)=dtsets(:)%getpawden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getpawden','INT',0)

 intarr(1,:)=dtsets(:)%getgam_eig2nkq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getgam_eig2nkq','INT',0)

 intarr(1,:)=dtsets(:)%getqps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getqps','INT',0)

 intarr(1,:)=dtsets(:)%getscr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getscr','INT',0)

 intarr(1,:)=dtsets(:)%getsuscep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getsuscep','INT',0)

 intarr(1,:)=dtsets(:)%getkss
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getkss','INT',0)

 intarr(1,:)=dtsets(:)%getocc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getocc','INT',0)

 intarr(1,:)=dtsets(:)%getvel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getvel','INT',0)

 intarr(1,:)=dtsets(:)%getwfk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getwfk','INT',0)

 intarr(1,:)=dtsets(:)%getwfq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getwfq','INT',0)

 intarr(1,:)=dtsets(:)%getxcart
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getxcart','INT',0)

 intarr(1,:)=dtsets(:)%getxred
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'getxred','INT',0)

 intarr(1,:)=dtsets(:)%get1den
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'get1den','INT',0)

 intarr(1,:)=dtsets(:)%get1wf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'get1wf','INT',0)

 intarr(1,:)=dtsets(:)%goprecon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'goprecon','INT',0)

 dprarr(1,:)=dtsets(:)%goprecprm(1)
 dprarr(2,:)=dtsets(:)%goprecprm(2)
 dprarr(3,:)=dtsets(:)%goprecprm(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'goprecprm','DPR',0)

 intarr(1,:)=dtsets(:)%gwcalctyp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwcalctyp','INT',0)

 intarr(1,:)=dtsets(:)%gwcomp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwcomp','INT',0)

 dprarr(1,:)=dtsets(:)%gwencomp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwencomp','ENE',0)

 intarr(1,:)=dtsets(:)%gwgamma
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwgamma','INT',0)

 if (ANY(dtsets(:)%gw_custom_freqsp/=0)) then
   intarr(1,:)=dtsets(:)%gw_custom_freqsp
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_custom_freqsp','INT',0)
 end if

 intarr(1,:)=dtsets(:)%gw_nqlwl
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_nqlwl','INT',0)

 intarr(1,:)=dtsets(:)%gw_sctype
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_sctype','INT',0)

 intarr(1,:)=dtsets(:)%gw_eet
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_eet','INT',0)

 intarr(1,:)=dtsets(:)%gw_eet_nband
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_eet_nband','INT',0)

 intarr(1,:)=dtsets(:)%gw_eet_inclvkb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_eet_inclvkb','INT',0)

 dprarr(1,:)=dtsets(:)%gw_eet_scale
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_eet_scale','DPR',0)

 intarr(1,:)=dtsets(:)%gw_sigxcore
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_sigxcore','INT',0)

 dprarr(1,:)=dtsets(:)%gw_toldfeig
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'gw_toldfeig','ENE',0)

 intarr(1,:)=dtsets(:)%gwmem
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwmem','INT',0)

 intarr(1,:)=dtsets(:)%gw_nstep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gw_nstep','INT',0)

 intarr(1,:)=dtsets(:)%gwpara
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwpara','INT',0)

 intarr(1,:)=dtsets(:)%gwrpacorr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'gwrpacorr','INT',0)

!gw_qlwl
 narr=3*dtsets(1)%gw_nqlwl ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%gw_nqlwl
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%gw_qlwl(1:3,1:dtsets(idtset)%gw_nqlwl),&
&     (/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'gw_qlwl','DPR',&
& multi_gw_nqlwl)

!gw_custom_freqsp
!It is not output, but it actually overrides the content of nfreqsp (which is forbidden !) in dtset.
!This is to be cleaned ...

!###########################################################
!### 03. Print all the input variables (I)
!##

!iatfix
 narr=natfix                    ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=natfix_(idtset)
   intarr(1:narrm(idtset),idtset)=iatfixio_(1:narrm(idtset),idtset)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'iatfix','INT',multi_natfix)

!iatfixx
 narr=natfixx                   ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=natfixx_(idtset)
   intarr(1:narrm(idtset),idtset)=iatfixx_(1:narrm(idtset),idtset)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'iatfixx','INT',multi_natfixx)

!iatfixy
 narr=natfixy                   ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=natfixy_(idtset)
   intarr(1:narrm(idtset),idtset)=iatfixy_(1:narrm(idtset),idtset)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'iatfixy','INT',multi_natfixy)

!iatfixz
 narr=natfixz                   ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=natfixz_(idtset)
   intarr(1:narrm(idtset),idtset)=iatfixz_(1:narrm(idtset),idtset)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'iatfixz','INT',multi_natfixz)

!iatsph
 multi_atsph=1
 narr=dtsets(1)%natsph          ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=dtsets(idtset)%natsph
!  Need to be printed only if there is some occurence of prtdos==3 or pawfatbnd
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%iatsph(1:narrm(idtset))
   end if
   if(dtsets(idtset)%prtdos==3.or.dtsets(idtset)%pawfatbnd>0)then
     narrm(idtset)=dtsets(idtset)%natsph
   else
     narrm(idtset)=0
   end if
 end do
 if (ndtset_alloc==1.and.sum(narrm(1:ndtset_alloc))==1) multi_atsph=0
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'iatsph','INT',multi_atsph) ! Emulating the case of multiple narr


 intarr(1,:)=dtsets(:)%iboxcut
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iboxcut','INT',0)

 intarr(1,:)=dtsets(:)%icutcoul
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'icutcoul','INT',0)

 intarr(1,:)=dtsets(:)%icoulomb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'icoulomb','INT',0)

 intarr(1,:)=dtsets(:)%idyson
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'idyson','INT',0)

 intarr(1,:)=dtsets(:)%ieig2rf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ieig2rf','INT',0)

 intarr(1,:)=dtsets(:)%ikhxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ikhxc','INT',0)

 intarr(1,:)=dtsets(:)%imgmov
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'imgmov','INT',0)

 intarr(1,:)=dtsets(:)%inclvkb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'inclvkb','INT',0)

 intarr(1,:)=dtsets(:)%intxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'intxc','INT',0)

 intarr(1,:)=dtsets(:)%intexact
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'intexact','INT',0)

 intarr(1,:)=dtsets(:)%ionmov
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ionmov','INT',0)

 intarr(1,:)=dtsets(:)%iextrapwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iextrapwf','INT',0)

 intarr(1,:)=dtsets(:)%iprcch
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iprcch','INT',0)

 intarr(1,:)=dtsets(:)%iprcel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iprcel','INT',0)

 intarr(1,:)=dtsets(:)%iprctfvw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iprctfvw','INT',0)

 intarr(1,:)=dtsets(:)%iprcfc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iprcfc','INT',0)

 intarr(1,:)=dtsets(:)%irandom
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irandom','INT',0)

 intarr(1,:)=dtsets(:)%irdbseig
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdbseig','INT',0)

 intarr(1,:)=dtsets(:)%irdbsreso
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdbsreso','INT',0)

 intarr(1,:)=dtsets(:)%irdbscoup
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdbscoup','INT',0)

 intarr(1,:)=dtsets(:)%irdhaydock
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdhaydock','INT',0)

 intarr(1,:)=dtsets(:)%irdddk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdddk','INT',0)

 intarr(1,:)=dtsets(:)%irdkss
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdkss','INT',0)

 intarr(1,:)=dtsets(:)%irdqps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdqps','INT',0)

 intarr(1,:)=dtsets(:)%irdscr
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdscr','INT',0)

 intarr(1,:)=dtsets(:)%irdsuscep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdsuscep','INT',0)

 intarr(1,:)=dtsets(:)%irdden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdden','INT',0)

 intarr(1,:)=dtsets(:)%irdpawden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdpawden','INT',0)

 intarr(1,:)=dtsets(:)%irdwfk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdwfk','INT',0)

 intarr(1,:)=dtsets(:)%irdwfq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'irdwfq','INT',0)

 intarr(1,:)=dtsets(:)%ird1den
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ird1den','INT',0)

 intarr(1,:)=dtsets(:)%ird1wf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ird1wf','INT',0)

 intarr(1,:)=dtsets(:)%iscf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'iscf','INT',0)

 intarr(1,:)=dtsets(:)%isecur
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'isecur','INT',0)

 intarr(1,:)=dtsets(:)%istatimg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'istatimg','INT',0)

 istatr_defo=49
 if(istatr/=istatr_defo)then
   intarr(1,:)=istatr
   intarr(1,0)=istatr_defo
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'istatr','INT',0)
 end if

 istatshft_defo=1
 if(istatshft/=istatshft_defo)then
   intarr(1,:)=istatshft
   intarr(1,0)=istatshft_defo
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'istatshft','INT',0)
 end if

!istwfk (must first restore the default istwf=0 for non-allowed k points)
 ABI_ALLOCATE(istwfk_2,(mxvals%mxnkpt,0:ndtset_alloc))
 do idtset=1,ndtset_alloc
   nqpt=dtsets(idtset)%nqpt
   do ikpt=1,dtsets(idtset)%nkpt
     allowed=1
     do ii=1,3
!      kpoint=dtsets(idtset)%kptns(ii)
       kpoint=dtsets(idtset)%kpt(ii,ikpt)/dtsets(idtset)%kptnrm
       if(nqpt/=0 .and. response_(idtset)==0)&
&       kpoint=kpoint+dtsets(idtset)%qptn(ii)
       if(abs(kpoint)>1.d-10 .and. abs(kpoint-0.5_dp)>1.e-10_dp )&
&       allowed=0
     end do
     if(allowed==0)then
       istwfk_2(ikpt,idtset)=0
     else
       istwfk_2(ikpt,idtset)=dtsets(idtset)%istwfk(ikpt)
     end if
   end do
 end do

!istwfk
 tnkpt=0
 intarr(1:marr,0)=0 ! default value
 do idtset=1,ndtset_alloc
   nkpt_eff=dtsets(idtset)%nkpt
   if(prtvol_glob==0 .and. nkpt_eff>nkpt_max)then
     nkpt_eff=nkpt_max
     tnkpt=1
   end if
   if((multi_nkpt/=0).and.(sum(istwfk_2(1:nkpt_eff,idtset))==0))&
&   nkpt_eff=0
   narrm(idtset)=nkpt_eff
   intarr(1:narrm(idtset),idtset)=istwfk_2(1:narrm(idtset),idtset)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,nkpt_eff,&
& narrm,ncid,ndtset_alloc,'istwfk','INT',multi_nkpt)
 if(tnkpt==1 .and. sum(istwfk_2(1:nkpt_eff,1:ndtset_alloc))/=0 ) &
& write(iout,'(23x,a,i3,a)' ) &
& 'outvar1 : Printing only first ',nkpt_max,' k-points.'
 ABI_DEALLOCATE(istwfk_2)

!ixc
 intarr(1,:)=dtsets(:)%ixc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ixc','INT',0)

 intarr(1,:)=dtsets(:)%ixcpositron
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ixcpositron','INT',0)

 intarr(1,:)=dtsets(:)%vdw_xc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'vdw_xc','INT',0)

!###########################################################
!### 03. Print all the input variables (J)
!##

 if (ndtset > 0) write(iout,format01155) 'jdtset',jdtset_(1:ndtset)

 intarr(1,:)=dtsets(:)%jellslab
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'jellslab','INT',0)

!###########################################################
!### 03. Print all the input variables (K)
!##

!kberry
 narr=3*dtsets(1)%nberry ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%nberry
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%kberry(1:3,1:dtsets(idtset)%nberry),&
&     (/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'kberry','INT',&
& multi_nberry)

!kpt
 tnkpt=0
 dprarr(:,0)=0
 narr=3*dtsets(1)%nkpt            ! default size for all datasets
 if(prtvol_glob==0 .and. narr>3*nkpt_max)then
   narr=3*nkpt_max
   tnkpt=1
 end if

 do idtset=1,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%nkpt
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=reshape(&
&     dtsets(idtset)%kpt(1:3,1:dtsets(idtset)%nkpt),&
&     (/ narrm(idtset) /) )
   end if

   if(prtvol_glob==0 .and. narrm(idtset)>3*nkpt_max)then
     narrm(idtset)=3*nkpt_max
     tnkpt=1
   end if

 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'kpt','DPR',multi_nkpt)

 if(tnkpt==1) write(iout,'(23x,a,i3,a)' ) &
& 'outvar1 : Printing only first ',nkpt_max,' k-points.'

!kptgw
 narr=3*dtsets(1)%nkptgw ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%nkptgw
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%kptgw(1:3,1:dtsets(idtset)%nkptgw),&
&     (/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'kptgw','DPR',&
& multi_nkptgw)

 dprarr(1,:)=dtsets(:)%kptnrm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'kptnrm','DPR',0)

 dprarr(1,:)=dtsets(:)%kptrlen
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'kptrlen','DPR',0)

 intarr(1,:)=dtsets(:)%kptopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'kptopt','INT',0)

!kptrlatt
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
   ndtset_kptopt=0
   intarr(1:9,0)=reshape( dtsets(0)%kptrlatt(:,:) , (/9/) )
   ABI_ALLOCATE(jdtset_kptopt,(0:ndtset_alloc))
!  Define the set of datasets for which kptopt>0
   do idtset=1,ndtset_alloc
     kptopt=dtsets(idtset)%kptopt
     if(kptopt>0)then
       ndtset_kptopt=ndtset_kptopt+1
       jdtset_kptopt(ndtset_kptopt)=jdtset_(idtset)
       intarr(1:9,ndtset_kptopt)=reshape( dtsets(idtset)%kptrlatt(:,:) , (/9/) )
     end if
   end do
   if(ndtset_kptopt>0)then
     call prttagm(dprarr,intarr,iout,jdtset_kptopt,3,marr,9,&
&     narrm,ncid,ndtset_kptopt,'kptrlatt','INT',0)
   end if
   ABI_DEALLOCATE(jdtset_kptopt)
 end if

 intarr(1,:)=dtsets(:)%kssform
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'kssform','INT',0)

!###########################################################
!### 03. Print all the input variables (L)
!##

 intarr(1,:)=dtsets(:)%ldgapp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ldgapp','INT',0)

 intarr(1,:)=dtsets(:)%localrdwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'localrdwf','INT',0)

!###########################################################
!### 03. Print all the input variables (M)
!##

 dprarr(1,:)=dtsets(:)%mdtemp(1)
 dprarr(2,:)=dtsets(:)%mdtemp(2)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'mdtemp','DPR',0)

 dprarr(1,:)=dtsets(:)%mdwall
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'mdwall','LEN',0)

 intarr(1,:)=dtsets(:)%mffmem
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'mffmem','INT',0)

!mixalch
 narr=dtsets(1)%npspalch*dtsets(1)%ntypalch ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=dtsets(idtset)%npspalch*dtsets(idtset)%ntypalch
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%mixalch(1:dtsets(idtset)%npspalch,&
&     1:dtsets(idtset)%ntypalch),&
&     (/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'mixalch','DPR',&
& multi_ntypalch)

 intarr(1,:)=dtsets(:)%maxnsym
 call prttagm(dprarr,intarr,iout,jdtset_,4,marr,1,narrm,ncid,ndtset_alloc,'maxnsym','INT',0)

 intarr(1,:)=dtsets(:)%mkmem
 call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,narrm,ncid,ndtset_alloc,'mkmem','INT',0)

 if(response==1)then
   intarr(1,:)=dtsets(:)%mkqmem
   call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,narrm,ncid,ndtset_alloc,'mkqmem','INT',0)
 end if

 if(response==1)then
   intarr(1,:)=dtsets(:)%mk1mem
   call prttagm(dprarr,intarr,iout,jdtset_,5,marr,1,narrm,ncid,ndtset_alloc,'mk1mem','INT',0)
 end if

 intarr(1,:)=dtsets(:)%mqgrid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'mqgrid','INT',0)
 if (usepaw==1) then
   intarr(1,:)=dtsets(:)%mqgriddg
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'mqgriddg','INT',0)
 end if

!###########################################################
!### 03. Print all the input variables (N)
!##

 intarr(1,:)=natfix_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natfix','INT',0)

 intarr(1,:)=natfixx_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natfixx','INT',0)

 intarr(1,:)=natfixy_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natfixy','INT',0)

 intarr(1,:)=natfixz_(:)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natfixz','INT',0)

 intarr(1,:)=dtsets(0:ndtset_alloc)%natom
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'natom','INT',0)

!natsph
!Need to be printed only if there is some occurence of prtdos==3 or
!pawfatbnd>0
 narr=1                      ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=1
   intarr(1,idtset)=dtsets(idtset)%natsph

   if(dtsets(idtset)%prtdos==3.or.dtsets(idtset)%pawfatbnd>0)then
     narrm(idtset)=1
   else
     narrm(idtset)=0
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'natsph','INT',multi_atsph) ! Emulating multiple size for narrm


!###########################################################
!### 03. Print all the input variables (BS_)
!##

!BEGIN VARIABLES FOR @Bethe-Salpeter
 intarr(1,:)=dtsets(:)%bs_nstates
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_nstates','INT',0)

 intarr(1,:)=dtsets(:)%bs_algorithm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_algorithm','INT',0)

 intarr(1,:)=dtsets(:)%bs_haydock_niter
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_haydock_niter','INT',0)

 intarr(1,:)=dtsets(:)%bs_hayd_term
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_hayd_term','INT',0)

 intarr(1,:)=dtsets(:)%bs_exchange_term
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_exchange_term','INT',0)

 intarr(1,:)=dtsets(:)%bs_coulomb_term
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_coulomb_term','INT',0)

 intarr(1,:)=dtsets(:)%bs_calctype
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_calctype','INT',0)

 intarr(1,:)=dtsets(:)%bs_coupling
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'bs_coupling','INT',0)

 do idtset=0,ndtset_alloc
   dprarr(1:2,idtset)=dtsets(idtset)%bs_haydock_tol(:)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'bs_haydock_tol','DPR',0)

 intarr(1,:)=dtsets(:)%bs_loband
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'bs_loband','INT',0)

 do idtset=0,ndtset_alloc
   dprarr(1:3,idtset)=dtsets(idtset)%bs_eh_cutoff(1:3)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'bs_eh_cutoff','ENE',0)

 do idtset=0,ndtset_alloc
   dprarr(1:3,idtset)=dtsets(idtset)%bs_freq_mesh(1:3)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'bs_freq_mesh','ENE',0)
!END VARIABLES FOR @Bethe-Salpeter.

!###########################################################
!### 04. Deallocation
!##

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(narrm)
 ABI_DEALLOCATE(nimagem)
 ABI_DEALLOCATE(dprarr_images)
 ABI_DEALLOCATE(natfix_)
 ABI_DEALLOCATE(iatfixio_)
 ABI_DEALLOCATE(natfixx_)
 ABI_DEALLOCATE(iatfixx_)
 ABI_DEALLOCATE(natfixy_)
 ABI_DEALLOCATE(iatfixy_)
 ABI_DEALLOCATE(natfixz_)
 ABI_DEALLOCATE(iatfixz_)
 ABI_DEALLOCATE(strimg)

!DEBUG
!write(std_out,*)' outvar1 : end of subroutine '
!if(.true.)stop
!ENDDEBUG
!
end subroutine outvar1
!!***
