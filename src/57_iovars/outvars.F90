!{\src2tex{textfont=tt}}
!!****f* ABINIT/outvars
!! NAME
!! outvars
!!
!! FUNCTION
!! Echo variables for the ABINIT code.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  choice= 1 if echo of preprocessed variables, 2 if echo after call driver
!!  dmatpuflag=flag controlling the use of an initial density matrix in PAW+U (max. value over datasets)
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for echoed output
!!  istatr=repetition rate for status file
!!  istatshft=shift of the repetition rate for status file
!!  mpi_enreg=informations about MPI parallelization
!!  mxvals=maximun size of some arrays along all datasets, including:
!!         mxgw_nqlwl   =maximal value of input gw_nqlwl for all the datasets
!!         mxlpawu      =maximal value of input lpawu for all the datasets
!!         mxmband      =maximum number of bands
!!         mxnatom      =maximal value of input natom for all the datasets
!!         mxnatpawu    =maximal value of number of atoms on which +U is applied for all the datasets
!!         mxnatsph     =maximal value of input natsph for all the datasets
!!         mxnatvshift  =maximal value of input natvshift for all the datasets
!!         mxnconeq     =maximal value of input nconeq for all the datasets
!!         mxnimage     =maximal value of input nimage for all the datasets
!!         mxnimfrqs    =maximal value of input cd_custom_imfrqs for all the datasets
!!         mxnkpt       =maximal value of input nkpt for all the datasets
!!         mxnkptgw     =maximal value of input nkptgw for all the datasets
!!         mxnnos       =maximal value of input nnos for all the datasets
!!         mxnqptdm     =maximal value of input nqptdm for all the datasets
!!         mxnspinor    =maximal value of input nspinor for all the datasets
!!         mxnsppol     =maximal value of input nsppol for all the datasets
!!         mxnsym       =maximum number of symmetries
!!         mxntypat     =maximum number of type of atoms
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!   one data set. Use for most dimensioned arrays.
!!  npsp=number of atom types
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!  timopt=input variable to modulate the timing
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
!! TODO
!!  MG: What about using modules to store the maximum dimensions as global variables?
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      create_nc_file,etsf_io_low_close,outvar1,prtocc,prttagm,prttagm_images
!!      wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outvars(choice,dmatpuflag,dtsets,filnam4,iout,istatr,&
     & istatshft,mpi_enreg,mxvals,&
     & ndtset,ndtset_alloc,npsp,results_out,timopt)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_results_out
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif
!#if defined HAVE_TRIO_NETCDF
! use netcdf
!#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outvars'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_57_iovars, except_this_one => outvars
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,dmatpuflag,iout,istatr,istatshft
 integer,intent(in) :: ndtset,ndtset_alloc,npsp,timopt
 type(ab_maxvals),intent(in) :: mxvals
 type(MPI_type),intent(in) :: mpi_enreg
 character(len=*),intent(in) :: filnam4
!arrays
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
!scalars
 integer,parameter :: nkpt_max=50
 integer :: first,iat,icount,idtset,ii,iimage,kptopt,lpawu!,jdtset
 integer :: lpawu1,marr,mu,multi_lpawu,multi_mxnsp,multi_natom,multi_natpawu
 integer :: multi_nconeq,multi_nkpt,multi_nnos,multi_nqptdm,multi_nshiftk,multi_kptopt
 integer :: multi_nspinor,multi_nsppol,multi_nsym,multi_ntypat,multi_occopt,mxnsp,natom,natpawu
 integer :: nconeq,ndtset_kptopt,nnos,nshiftk,narr!,nkpt_eff
 integer :: nsym,ntypat,prtvol_glob,response
 integer :: rfddk,rfelfd,rfphon,rfstrs,rfuser
 integer :: timopt_default,tnkpt,usepaw
 integer :: ncid=0 ! Variables for NetCDF output
 logical :: compute_static_images
! character(len=4) :: appen
 character(len=500) :: message
 character(len=4) :: stringimage
 logical :: lstat
!arrays
 integer,allocatable :: intarr(:,:),jdtset_(:),jdtset_kptopt(:),response_(:)
 integer,allocatable :: narrm(:)
 integer,allocatable :: nimagem(:),prtimg(:,:)
 real(dp) :: rprimd(3,3)
 real(dp),allocatable :: dprarr(:,:),xangst(:,:),xangst_(:,:,:,:),xcart(:,:)
 real(dp),allocatable :: xcart_(:,:,:,:),xred(:,:),dprarr_images(:,:,:)
 character(len=8),allocatable :: strimg(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' outvars : enter '
!ENDDEBUG

!XG 080819 Warning : do not remove the following (silly) line :
!this is needed to avoid a compiler bug with g95
 ABI_ALLOCATE(dprarr,(1,1))
 ABI_DEALLOCATE(dprarr)

!Set up a 'global' prtvol value
 prtvol_glob=1
 if(sum((dtsets(:)%prtvol)**2)==0)prtvol_glob=0

!write(std_out,*) 'outvar 00'
!###########################################################
!### 00. Echo of selected default values

 if(choice==1)then
   write(iout, '(10a)' )&
&   '--------------------------------------------------------------------------------',ch10,&
&   '------------- Echo of variables that govern the present computation ------------',ch10,&
&   '--------------------------------------------------------------------------------',ch10,&
&   '-',ch10,&
&   '- outvars: echo of selected default values                                      '
   write(iout, '(3(a,i3),2a)' )&
&   '-   accesswff0 =',dtsets(0)%accesswff,&
&   ' , fftalg0 =',dtsets(0)%ngfft(7),&
&   ' , wfoptalg0 =',dtsets(0)%wfoptalg,&
&   ch10,'-'
   write(iout, '(3a,i5,3a)' )&
&   '- outvars: echo of global parameters not present in the input file              ',ch10,&
&   '-   nproc =',mpi_enreg%nproc,ch10,&
&   '-'
 end if

!write(std_out,*) 'outvar 01'
!###########################################################
!### 01. First line indicating outvars

 if(choice==1)then
   write(iout, '(a)' )&
&   ' -outvars: echo values of preprocessed input variables --------'
 else
   write(iout, '(a)' )&
&   ' -outvars: echo values of variables after computation  --------'
 end if


#if defined HAVE_TRIO_ETSF_IO
 if (iout==std_out)then
   write(iout,*) ch10,' These variables are accessible in NetCDF format (',filnam4//'_OUT.nc',')',ch10
 end if
#endif

!write(std_out,*) 'outvar 02'
!###########################################################
!### 02. Set up dimensions

 marr=max(3*mxvals%mxnatom,&
& 3*mxvals%mxnkptgw,&
& mxvals%mxnkpt*mxvals%mxnsppol*mxvals%mxmband,&
& 3*mxvals%mxnkpt,npsp,&
& mxvals%mxntypat,&
& 9*mxvals%mxnsym,3*8,&
& 3*mxvals%mxnatom*mxvals%mxnconeq,&
& mxvals%mxnnos,&
& 3*mxvals%mxnqptdm,&
& 3*mxvals%mxgw_nqlwl,&
& (2*mxvals%mxlpawu+1)**2*max(mxvals%mxnsppol,mxvals%mxnspinor)*mxvals%mxnatpawu*dmatpuflag)
 ABI_ALLOCATE(dprarr,(marr,0:ndtset_alloc))
 ABI_ALLOCATE(dprarr_images,(marr,mxvals%mxnimage,0:ndtset_alloc))
 ABI_ALLOCATE(intarr,(marr,0:ndtset_alloc))
 ABI_ALLOCATE(narrm,(0:ndtset_alloc))
 ABI_ALLOCATE(nimagem,(0:ndtset_alloc))
 ABI_ALLOCATE(prtimg,(mxvals%mxnimage,0:ndtset_alloc))

!Set up dimensions : determine whether these are different for different
!datasets.
 multi_natom=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%natom/=dtsets(idtset)%natom)multi_natom=1
   end do
 end if
 if(multi_natom==0)natom=dtsets(1)%natom

 multi_nconeq=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nconeq/=dtsets(idtset)%nconeq)multi_nconeq=1
   end do
 end if
 if(multi_nconeq==0)nconeq=dtsets(1)%nconeq

 multi_nnos=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nnos/=dtsets(idtset)%nnos)multi_nnos=1
   end do
 end if
 if(multi_nnos==0)nnos=dtsets(1)%nnos

 multi_nqptdm=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nqptdm/=dtsets(idtset)%nqptdm)multi_nqptdm=1
   end do
 end if

 multi_nshiftk=0
 nshiftk=1
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

 multi_nkpt=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nkpt/=dtsets(idtset)%nkpt)multi_nkpt=1
   end do
 end if

 multi_nsppol=0;multi_nspinor=0;multi_mxnsp=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nsppol/=dtsets(idtset)%nsppol)multi_nsppol=1
     if(dtsets(1)%nspinor/=dtsets(idtset)%nspinor)multi_nspinor=1
     if(dtsets(1)%nsppol*dtsets(1)%nspinor/=dtsets(idtset)%nsppol*dtsets(idtset)%nspinor)multi_mxnsp=1
   end do
 end if

 multi_nsym=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%nsym/=dtsets(idtset)%nsym)multi_nsym=1
   end do
 end if
 if(multi_nsym==0)nsym=dtsets(1)%nsym

 multi_ntypat=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%ntypat/=dtsets(idtset)%ntypat)multi_ntypat=1
   end do
 end if
 if(multi_ntypat==0)ntypat=dtsets(1)%ntypat

 multi_occopt=0
 if(ndtset_alloc>1)then
   do idtset=1,ndtset_alloc
     if(dtsets(1)%occopt/=dtsets(idtset)%occopt)multi_occopt=1
   end do
 end if

 response=0
 ABI_ALLOCATE(response_,(ndtset_alloc))
 response_(:)=0
 do idtset=1,ndtset_alloc
   rfddk=dtsets(idtset)%rfddk
   rfelfd=dtsets(idtset)%rfelfd
   rfphon=dtsets(idtset)%rfphon
   rfstrs=dtsets(idtset)%rfstrs
   rfuser=dtsets(idtset)%rfuser
   if(rfddk/=0 .or. rfelfd/=0 .or. &
&   rfphon/=0 .or. rfstrs/=0 .or. rfuser/=0)then
     response_(idtset)=1 ; response=1
   end if
 end do

!Must compute xangst and xcart
 ABI_ALLOCATE(xangst_,(3,mxvals%mxnatom,mxvals%mxnimage,0:ndtset_alloc))
 ABI_ALLOCATE(xcart_,(3,mxvals%mxnatom,mxvals%mxnimage,0:ndtset_alloc))
 xangst_(:,:,:,:)=0.0_dp ; xcart_(:,:,:,:)=0.0_dp
 do idtset=1,ndtset_alloc
   natom=dtsets(idtset)%natom
   ABI_ALLOCATE(xred,(3,natom))
   ABI_ALLOCATE(xangst,(3,natom))
   ABI_ALLOCATE(xcart,(3,natom))
   do iimage=1,dtsets(idtset)%nimage
     xred(:,1:natom)=results_out(idtset)%xred(:,1:natom,iimage)
     rprimd(:,:)     =results_out(idtset)%rprimd(:,:,iimage)
!    Compute xcart from xred and rprimd
     call xredxcart(natom,1,rprimd,xcart,xred)
!    Compute xangst from xcart
     xangst(:,:)=xcart(:,:)*Bohr_Ang
!    Save the data
     xangst_(1:3,1:natom,iimage,idtset)=xangst(:,:)
     xcart_(1:3,1:natom,iimage,idtset)=xcart(:,:)
   end do
   if(dtsets(idtset)%nimage/=mxvals%mxnimage)then
     xangst_(1:3,1:natom,dtsets(idtset)%nimage+1:mxvals%mxnimage,idtset)=zero
     xcart_(1:3,1:natom,dtsets(idtset)%nimage+1:mxvals%mxnimage,idtset)=zero
   end if
   ABI_DEALLOCATE(xred)
   ABI_DEALLOCATE(xangst)
   ABI_DEALLOCATE(xcart)
 end do

 ABI_ALLOCATE(jdtset_,(0:ndtset_alloc))
 jdtset_(0:ndtset_alloc)=dtsets(0:ndtset_alloc)%jdtset

 ABI_ALLOCATE(strimg,(mxvals%mxnimage))
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

!Determine whether we are in a PAW run
 usepaw=0;if (maxval(dtsets(0:ndtset_alloc)%usepaw)==1) usepaw=1

!DEBUG
!write(std_out,*)' outvars : before outvar1 '
!stop
!ENDDEBUG

!write(std_out,*) 'outvar 03'
!###########################################################
!### 03. Open NetCDF file for export variables

#if defined HAVE_TRIO_ETSF_IO
 call create_nc_file(filnam4//"_OUT.nc",ncid)
#endif

 if (dtsets(1)%prtvol==-1) then
   if (ncid>0)then
     ncid=-ncid
   else
     ncid=-1
   end if
 end if

!write(std_out,*) 'outvar 04'
!###########################################################
!### 04. Print variables from acell to natom

!Print variables between acell and natom (by alphabetic order)
 call outvar1(choice,dtsets,iout,istatr,istatshft,jdtset_,mxvals,&
& ncid,ndtset,ndtset_alloc,npsp,prtvol_glob,response,response_,&
& results_out,usepaw)

!DEBUG
!write(std_out,*)' outvars : after outvar1 '
!stop
!ENDDEBUG

!write(std_out,*) 'outvar 05'
!###########################################################
!### 05. Print variables from nband to znucl

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
!!VARIABLE
!dprarr(1:marr,0)=0             ! default value
!narr=SIZE                      ! default size for all datasets
!do idtset=0,ndtset_alloc       ! especific size for each dataset
!narrm(idtset)=dtsets(idtset)%SIZE
!if (narrm(idtset)>0) then
!dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%VARIABLE(1:narrm(idtset))
!end if
!end do
!call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
!& narrm,ncid,ndtset_alloc,'VARIABLE','DPR',MULTI)
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
!### 03. Print all the input variables (N)
!##

!nband
 if(dtsets(1)%occopt==2)then
   narr=dtsets(1)%nkpt*dtsets(1)%nsppol                      ! default size for all datasets
 else
   narr=1
 end if
 do idtset=0,ndtset_alloc       ! especific size for each dataset

   if(dtsets(idtset)%occopt==2)then
     narrm(idtset)=dtsets(idtset)%nkpt*dtsets(idtset)%nsppol
   else
     narrm(idtset)=1
   end if

   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%nband(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,&
& narrm,ncid,ndtset_alloc,'nband','INT',multi_nkpt+multi_nsppol+multi_occopt)


 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nbandsus
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nbandsus','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nbdblock
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nbdblock','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nbdbuf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nbdbuf','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%nberry
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nberry','INT',0)

 intarr(1,:)=dtsets(:)%nconeq
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nconeq','INT',0)

 intarr(1,:)=dtsets(:)%nctime
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nctime','INT',0)

!ndtset
 if(ndtset>0)then
   intarr(1,:)=ndtset
   intarr(1,0)=0
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ndtset','INT',0)
 end if

 intarr(1,:)=dtsets(:)%ndynimage
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ndynimage','INT',0)

 intarr(1,:)=dtsets(:)%ndyson
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ndyson','INT',0)

 intarr(1,:)=dtsets(:)%nfreqim
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nfreqim','INT',0)

 intarr(1,:)=dtsets(:)%nfreqre
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nfreqre','INT',0)

 intarr(1,:)=dtsets(:)%nfreqsp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nfreqsp','INT',0)

 intarr(1,:)=dtsets(:)%nfreqsus
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nfreqsus','INT',0)

 intarr(1,:)=dtsets(:)%ngfft(1)
 intarr(2,:)=dtsets(:)%ngfft(2)
 intarr(3,:)=dtsets(:)%ngfft(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'ngfft','INT',0)

 if (usepaw==1) then
   intarr(1,:)=dtsets(:)%ngfftdg(1)
   intarr(2,:)=dtsets(:)%ngfftdg(2)
   intarr(3,:)=dtsets(:)%ngfftdg(3)
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'ngfftdg','INT',0)
 end if

 intarr(1,:)=dtsets(:)%ngroup_rf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ngroup_rf','INT',0)

 intarr(1,:)=dtsets(:)%nimage
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nimage','INT',0)

 intarr(1,:)=dtsets(:)%nkptgw
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nkptgw','INT',0)

 intarr(1,:)=dtsets(:)%nkpt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nkpt','INT',0)

 intarr(1,:)=dtsets(:)%nline
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nline','INT',0)

 intarr(1,:)=dtsets(:)%nloalg(1)+10*dtsets(:)%nloalg(5)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nloalg','INT',0)

 intarr(1,:)=dtsets(:)%nnos
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nnos','INT',0)

 intarr(1,:)=dtsets(:)%nnsclo
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nnsclo','INT',0)

 intarr(1,:)=dtsets(:)%nomegasf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nomegasf','INT',0)

 intarr(1,:)=dtsets(:)%nomegasi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nomegasi','INT',0)

 intarr(1,:)=dtsets(:)%nomegasrd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nomegasrd','INT',0)

 dprarr(1,:)=dtsets(:)%noseinert
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'noseinert','DPR',0)

 intarr(1,:)=dtsets(:)%npband
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npband','INT',0)

 intarr(1,:)=dtsets(:)%npfft
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npfft','INT',0)

 intarr(1,:)=dtsets(:)%npimage
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npimage','INT',0)

 intarr(1,:)=dtsets(:)%npkpt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npkpt','INT',0)

 if(multi_ntypat/=0)then
   intarr(1,:)=dtsets(:)%npsp
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npsp','INT',0)
 else if(multi_ntypat==0 .and. ntypat/=npsp)then
   intarr(1,:)=dtsets(:)%npsp
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npsp','INT',0)
 end if

 intarr(1,:)=dtsets(0:ndtset_alloc)%npulayit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npulayit','INT',0)

!intarr(1,:)=dtsets(:)%npspalch
!call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npspalch','INT',0)

 intarr(1,:)=dtsets(:)%npweps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npweps','INT',0)

 intarr(1,:)=dtsets(:)%npwkss
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npwkss','INT',0)

 intarr(1,:)=dtsets(:)%npwsigx
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npwsigx','INT',0)

 intarr(1,:)=dtsets(:)%npwwfn
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'npwwfn','INT',0)

 intarr(1,:)=dtsets(:)%nqpt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nqpt','INT',0)

 intarr(1,:)=dtsets(:)%nqptdm
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'nqptdm','INT',0)

 intarr(1,:)=dtsets(:)%nscforder
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nscforder','INT',0)

 intarr(1,:)=dtsets(:)%nsheps
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nsheps','INT',0)

!nshiftk
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then
   ndtset_kptopt=0
   intarr(1:1,0)=dtsets(0)%nshiftk
   ABI_ALLOCATE(jdtset_kptopt,(0:ndtset_alloc))
!  Define the set of datasets for which kptopt>0
   do idtset=1,ndtset_alloc
     kptopt=dtsets(idtset)%kptopt
     if(kptopt>0)then
       ndtset_kptopt=ndtset_kptopt+1
       jdtset_kptopt(ndtset_kptopt)=jdtset_(idtset)
       intarr(1:1,ndtset_kptopt)=dtsets(idtset)%nshiftk
     end if
   end do
   if(ndtset_kptopt>0)then
     call prttagm(dprarr,intarr,iout,jdtset_kptopt,2,marr,1,&
&     narrm,ncid,ndtset_kptopt,'nshiftk','INT',0)
   end if
   ABI_DEALLOCATE(jdtset_kptopt)
 end if

 intarr(1,:)=dtsets(:)%nshsigx
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nshsigx','INT',0)

 intarr(1,:)=dtsets(:)%nshwfn
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nshwfn','INT',0)

 intarr(1,:)=dtsets(:)%nspden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nspden','INT',0)

 intarr(1,:)=dtsets(:)%nspinor
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nspinor','INT',0)

 intarr(1,:)=dtsets(:)%nsppol
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nsppol','INT',0)

 intarr(1,:)=dtsets(:)%nstep
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nstep','INT',0)

 intarr(1,:)=dtsets(:)%nsym
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nsym','INT',0)

 intarr(1,:)=dtsets(:)%ntime
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ntime','INT',0)

 intarr(1,:)=dtsets(:)%ntimimage
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ntimimage','INT',0)

 intarr(1,:)=dtsets(:)%ntypalch
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ntypalch','INT',0)

!intarr(1,:)=dtsets(:)%ntyppure
!call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ntyppure','INT',0)

 intarr(1,:)=dtsets(:)%ntypat
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ntypat','INT',0)

 intarr(1,:)=dtsets(:)%nwfshist
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'nwfshist','INT',0)

!###########################################################
!### 03. Print all the input variables (O)
!##

!occ
!The use of prttagm for occ if occopt>=2 is not possible because
!the different k-point and spins must be separated on different lines ...
 if(mxvals%mxnimage==1)then
   call prtocc(dtsets,iout,jdtset_,ndtset_alloc,prtvol_glob,results_out)
 end if

 intarr(1,:)=dtsets(:)%occopt
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'occopt','INT',0)

 dprarr(1,:)=dtsets(:)%omegasimax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'omegasimax','ENE',0)

 dprarr(1,:)=dtsets(:)%omegasrdmax
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'omegasrdmax','ENE',0)

 intarr(1,:)=dtsets(:)%optcell
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'optcell','INT',0)

 intarr(1,:)=dtsets(:)%optdriver
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'optdriver','INT',0)

 intarr(1,:)=dtsets(:)%optforces
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'optforces','INT',0)

 intarr(1,:)=dtsets(:)%optstress
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'optstress','INT',0)

 intarr(1,:)=dtsets(:)%optnlxccc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'optnlxccc','INT',0)

 intarr(1,:)=dtsets(:)%ortalg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ortalg','INT',0)

!###########################################################
!### 03. Print all the input variables (P)
!##

 intarr(1,:)=dtsets(:)%paral_rf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'paral_rf','INT',0)

 intarr(1,:)=dtsets(:)%paral_kgb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'paral_kgb','INT',0)


 if (usepaw==1) then

   intarr(1,:)=dtsets(:)%pawcpxocc
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawcpxocc','INT',0)

   intarr(1,:)=dtsets(:)%pawcross
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawcross','INT',0)

   dprarr(1,:)=dtsets(:)%pawecutdg
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'pawecutdg','ENE',0)

   intarr(1,:)=dtsets(:)%pawlcutd
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawlcutd','INT',0)

   intarr(1,:)=dtsets(:)%pawlmix
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawlmix','INT',0)

   intarr(1,:)=dtsets(:)%pawmixdg
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawmixdg','INT',0)

   intarr(1,:)=dtsets(:)%pawnhatxc
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawnhatxc','INT',0)

   intarr(1,:)=dtsets(:)%pawnphi
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawnphi','INT',0)

   intarr(1,:)=dtsets(:)%pawntheta
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawntheta','INT',0)

   intarr(1,:)=dtsets(:)%pawnzlm
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawnzlm','INT',0)

   intarr(1,:)=dtsets(:)%pawoptmix
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawoptmix','INT',0)

   dprarr(1,:)=dtsets(:)%pawovlp
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawovlp','DPR',0)

   intarr(1,:)=dtsets(:)%pawprtden
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprtden','INT',0)

   intarr(1,:)=dtsets(:)%pawprtdos
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprtdos','INT',0)

   intarr(1,:)=dtsets(:)%pawprtvol
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprtvol','INT',0)

   intarr(1,:)=dtsets(:)%pawprtwf
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprtwf','INT',0)

   intarr(1,:)=dtsets(:)%pawprt_k
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprt_k','INT',0)

   intarr(1,:)=dtsets(:)%pawprt_b
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawprt_b','INT',0)

   intarr(1,:)=dtsets(:)%pawspnorb
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawspnorb','INT',0)

   intarr(1,:)=dtsets(:)%pawstgylm
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawstgylm','INT',0)

   intarr(1,:)=dtsets(:)%pawsushat
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawsushat','INT',0)

   intarr(1,:)=dtsets(:)%pawusecp
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawusecp','INT',0)

   intarr(1,:)=dtsets(:)%pawxcdev
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pawxcdev','INT',0)

 end if

!###########################################################
!### 03. Print all the input variables (Q)
!##

!qptdm
 narr=3*dtsets(1)%nqptdm ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%nqptdm
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%qptdm(1:3,&
&     1:dtsets(idtset)%nqptdm),&
&     (/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'qptdm','DPR',&
& multi_nqptdm)

 if (usepaw==1) then

   icount=0
   do idtset=1,ndtset_alloc
     if (dtsets(idtset)%usedmft>0) icount=icount+1
   end do
   if(icount>0) then
!    TO DO : all these variables are not ordered by alphabetical order !! They should
!    be output at the right place in this outvars routine, or in the outvar1 routine
     intarr(1,:)=dtsets(:)%dmft_dc
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_dc','INT',0)
     intarr(1,:)=dtsets(:)%dmft_iter
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_iter','INT',0)
     intarr(1,:)=dtsets(:)%dmft_mxsf
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_mxsf','DPR',0)
     intarr(1,:)=dtsets(:)%dmft_nwli
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_nwli','INT',0)
     intarr(1,:)=dtsets(:)%dmft_nwlo
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_nwlo','INT',0)
     intarr(1,:)=dtsets(:)%dmft_rslf
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_rslf','INT',0)
     intarr(1,:)=dtsets(:)%dmft_solv
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmft_solv','INT',0)
     intarr(1,:)=dtsets(:)%dmftbandi
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmftbandi','INT',0)
     intarr(1,:)=dtsets(:)%dmftbandf
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmftbandf','INT',0)
     intarr(1,:)=dtsets(:)%dmftcheck
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmftcheck','INT',0)
     intarr(1,:)=dtsets(:)%usedmft
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usedmft','INT',0)
   end if

   icount=0
   do idtset=1,ndtset_alloc
     if (dtsets(idtset)%usepawu>0.or.dtsets(idtset)%usedmft>0) icount=icount+1
   end do
   if(icount>0) then
!    TO DO : all these variables are not ordered by alphabetical order !! They should
!    be output at the right place in this outvars routine, or in the outvar1 routine
     intarr(1,:)=dtsets(:)%dmatudiag
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmatudiag','INT',0)
     intarr(1,:)=dtsets(:)%dmatpuopt
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'dmatpuopt','INT',0)
     intarr(1,:)=dtsets(:)%usedmatpu
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usedmatpu','INT',0)
     intarr(1,:)=dtsets(:)%usepawu
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usepawu','INT',0)


!    jpawu
     dprarr(:,0)=0
     narr=ntypat                    ! default size for all datasets
     do idtset=1,ndtset_alloc       ! especific size for each dataset
       narrm(idtset)=dtsets(idtset)%ntypat
       if (narrm(idtset)>0) then
         dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%jpawu(1:narrm(idtset))
       end if
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
&     narrm,ncid,ndtset_alloc,'jpawu','DPR',multi_ntypat)

!    lpawu
     intarr(:,0)=0
     narr=ntypat                    ! default size for all datasets
     do idtset=1,ndtset_alloc       ! especific size for each dataset
       narrm(idtset)=dtsets(idtset)%ntypat
       if (narrm(idtset)>0) then
         intarr(1:narrm(idtset),idtset)=dtsets(idtset)%lpawu(1:narrm(idtset))
       end if
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,&
&     narrm,ncid,ndtset_alloc,'lpawu','INT',multi_ntypat)

!    upawu
     dprarr(:,0)=0
     narr=ntypat                    ! default size for all datasets
     do idtset=1,ndtset_alloc       ! especific size for each dataset
       narrm(idtset)=dtsets(idtset)%ntypat
       if (narrm(idtset)>0) then
         dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%upawu(1:narrm(idtset))
       end if
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
&     narrm,ncid,ndtset_alloc,'upawu','DPR',multi_ntypat)

!    dmatpawu
     if (dmatpuflag==1.and.mxvals%mxnatpawu>0) then
       multi_lpawu=0;multi_natpawu=0
       lpawu=maxval(dtsets(1)%lpawu(:))
       natpawu=dtsets(1)%natpawu
       do idtset=1,ndtset_alloc
         lpawu1=maxval(dtsets(idtset)%lpawu(:))
         if(lpawu/=lpawu1) multi_lpawu=1
         if(dtsets(idtset)%natpawu/=natpawu) multi_natpawu=1
       end do

       mxnsp=max(dtsets(1)%nsppol,dtsets(1)%nspinor)
       narr=(2*lpawu+1)*(2*lpawu+1)*mxnsp*natpawu ! default size for all datasets
       if (lpawu==-1)then
         narr=0
       end if

       dprarr(:,:)=0
       do idtset=0,ndtset_alloc       ! especific size for each dataset
         mxnsp=max(dtsets(idtset)%nsppol,dtsets(idtset)%nspinor)
         lpawu1=maxval(dtsets(idtset)%lpawu(:))
         narrm(idtset)=(2*lpawu1+1)*(2*lpawu1+1)*mxnsp*dtsets(idtset)%natpawu
         if (narrm(idtset)>0) then
           dprarr(1:narrm(idtset),idtset)=&
&           reshape(dtsets(idtset)%dmatpawu(&
&           1:2*lpawu1+1,&
&           1:2*lpawu1+1,&
&           1:mxnsp,&
&           1:dtsets(idtset)%natpawu),&
&           (/ narrm(idtset) /) )
         end if
       end do
       call prttagm(dprarr,intarr,iout,jdtset_,5,marr,narr,&
&       narrm,ncid,ndtset_alloc,'dmatpawu','DPR',&
&       multi_lpawu+multi_mxnsp+multi_natpawu)

     end if
   end if
   icount=0
   do idtset=1,ndtset_alloc
     if (dtsets(idtset)%useexexch>0) icount=icount+1
   end do
   if(icount>0) then
     intarr(1,:)=dtsets(:)%useexexch
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'useexexch','INT',0)


!    lexexch
     narr=ntypat                    ! default size for all datasets
     do idtset=0,ndtset_alloc       ! especific size for each dataset
       narrm(idtset)=dtsets(idtset)%ntypat
       if (narrm(idtset)>0) then
         intarr(1:narrm(idtset),idtset)=dtsets(idtset)%lexexch(1:narrm(idtset))
       end if
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,&
&     narrm,ncid,ndtset_alloc,'lexexch','INT',multi_ntypat)

   end if

   dprarr(1,:)=dtsets(:)%spnorbscl
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'spnorbscl','DPR',0)

!  intarr(1,:)=dtsets(:)%ngfftdg(1)
!  intarr(2,:)=dtsets(:)%ngfftdg(2)
!  intarr(3,:)=dtsets(:)%ngfftdg(3)
!  call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'ngfftdg','INT',0)

 end if

!pimass
 icount=0
 do idtset=0, ndtset_alloc
   do ii = 1, ntypat
     dprarr(ii,idtset) = dtsets(idtset)%pimass(ii)
     if (dtsets(idtset)%pimass(ii)/=dtsets(idtset)%amu(ii)) icount=1
   end do ! end loop over ntypat
 end do ! end loop over datasets
 if (icount/=0) then
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,narrm,ncid,ndtset_alloc,'pimass','DPR',0)
 end if

 intarr(1,:)=dtsets(:)%pitransform
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'pitransform','INT',0)
 intarr(1,:)=dtsets(:)%positron
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'positron','INT',0)
 intarr(1,:)=dtsets(:)%posnstep
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'posnstep','INT',0)
 dprarr(1,:)=dtsets(:)%posocc
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'posocc','DPR',0)
 dprarr(1,:)=dtsets(:)%postoldfe
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'postoldfe','ENE',0)
 dprarr(1,:)=dtsets(:)%postoldff
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'postoldff','DPR',0)

 dprarr(1,:)=dtsets(:)%ppmfrq
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'ppmfrq','ENE',0)

 intarr(1,:)=dtsets(:)%ppmodel
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ppmodel','INT',0)

 do idtset=0, ndtset_alloc
   do ii = 1, ntypat
     dprarr(ii,idtset) = dtsets(idtset)%ptcharge(ii)
   end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,narrm,ncid,ndtset_alloc,'ptcharge','DPR',0)


 intarr(1,:)=dtsets(:)%prepanl
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prepanl','INT',0)

 intarr(1,:)=dtsets(:)%prepgkk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prepgkk','INT',0)

 intarr(1,:)=dtsets(:)%prepscphon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prepscphon','INT',0)

 intarr(1,:)=dtsets(:)%prtbbb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtbbb','INT',0)

 intarr(1,:)=dtsets(:)%prtbltztrp
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtbltztrp','INT',0)

 intarr(1,:)=dtsets(:)%prtcif
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtcif','INT',0)

 intarr(1,:)=dtsets(:)%prtcml
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtcml','INT',0)

 intarr(1,:)=dtsets(:)%prtcs
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtcs','INT',0)

 intarr(1,:)=dtsets(:)%prtden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtden','INT',0)

 intarr(1,:)=dtsets(:)%prtdensph
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtdensph','INT',0)

 intarr(1,:)=dtsets(:)%prtdos
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtdos','INT',0)

 intarr(1,:)=dtsets(:)%prtdosm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtdosm','INT',0)

 intarr(1,:)=dtsets(:)%prtefg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtefg','INT',0)

 intarr(1,:)=dtsets(:)%prteig
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prteig','INT',0)

 intarr(1,:)=dtsets(:)%prtelf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtelf','INT',0)

 intarr(1,:)=dtsets(:)%prtfc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtfc','INT',0)

 intarr(1,:)=dtsets(:)%prtfsurf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtfsurf','INT',0)

 intarr(1,:)=dtsets(:)%prtgden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtgden','INT',0)

 intarr(1,:)=dtsets(:)%prtgeo
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtgeo','INT',0)

 intarr(1,:)=dtsets(:)%prtgkk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtgkk','INT',0)

 intarr(1,:)=dtsets(:)%prtkden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtkden','INT',0)

 intarr(1,:)=dtsets(:)%prtlden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtlden','INT',0)

 intarr(1,:)=dtsets(:)%prtnabla
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtnabla','INT',0)

 intarr(1,:)=dtsets(:)%prtposcar
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtposcar','INT',0)

 intarr(1,:)=dtsets(:)%prtpot
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtpot','INT',0)

 intarr(1,:)=dtsets(:)%prtspcur
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtspcur','INT',0)

 intarr(1,:)=dtsets(:)%prtstm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtstm','INT',0)

 intarr(1,:)=dtsets(:)%prtvha
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvha','INT',0)

 intarr(1,:)=dtsets(:)%prtvhxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvhxc','INT',0)

 intarr(1,:)=dtsets(:)%prtvxc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvxc','INT',0)

 intarr(1,:)=dtsets(:)%prtvol
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvol','INT',0)

 intarr(1,:)=dtsets(:)%prtvolimg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtvolimg','INT',0)

 intarr(1,:)=dtsets(:)%prtwant
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtwant','INT',0)

 intarr(1,:)=dtsets(:)%prtwf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtwf','INT',0)

 intarr(1,:)=dtsets(:)%prtxml
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prtxml','INT',0)

!prtatlist
 if(multi_natom==0)then
   do idtset=0,ndtset_alloc
     intarr(1:natom,idtset)=dtsets(idtset)%prtatlist(1:natom)
   end do
   intarr(1:mxvals%mxnatom,0)=(/ (ii,ii=1,mxvals%mxnatom) /)
   call prttagm(dprarr,intarr,iout,jdtset_,4,marr,natom,narrm,ncid,ndtset_alloc,'prtatlist','INT',0)
 else
!  This thing will disapear with new generalized prttagm
 end if

!prt1dm
 intarr(1,:)=dtsets(:)%prt1dm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'prt1dm','INT',0)

 intarr(1,:)=dtsets(:)%ptgroupma
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'ptgroupma','INT',0)

!qmass
 narr=nnos ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=dtsets(idtset)%nnos
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%qmass(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'qmass','DPR',&
& multi_nnos)

 intarr(1,:)=dtsets(:)%qprtrb(1)
 intarr(2,:)=dtsets(:)%qprtrb(2)
 intarr(3,:)=dtsets(:)%qprtrb(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'qprtrb','INT',0)

 dprarr(1,:)=dtsets(:)%qptn(1)
 dprarr(2,:)=dtsets(:)%qptn(2)
 dprarr(3,:)=dtsets(:)%qptn(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'qpt','DPR',0)

 do idtset=0, ndtset_alloc
   do ii = 1, ntypat
     dprarr(ii,idtset) = dtsets(idtset)%quadmom(ii)
   end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,narrm,ncid,ndtset_alloc,'quadmom','DPR',0)

 do idtset=0, ndtset_alloc
   do ii = 1, ntypat
     dprarr(ii,idtset) = dtsets(idtset)%ratsph(ii)
   end do ! end loop over ntypat
 end do ! end loop over datasets
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,ntypat,narrm,ncid,ndtset_alloc,'ratsph','LEN',0)

 dprarr(1,:)=dtsets(:)%rcut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'rcut','LEN',0)

 intarr(1,:)=dtsets(:)%rdmnb
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rdmnb','INT',0)

!variables used for the random positions in unit cell
 intarr(1,:)=dtsets(:)%random_atpos
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'random_atpos','INT',0)

!Variables used for recursion method
 dprarr(1,:)=dtsets(:)%recefermi
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'recefermi','ENE',0)

 intarr(1,:)=dtsets(:)%recnpath
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'recnpath','INT',0)

 intarr(1,:)=dtsets(:)%recnrec
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'recnrec','INT',0)

 intarr(1,:)=dtsets(:)%recptrott
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'recptrott','INT',0)

 dprarr(1,:)=dtsets(:)%recrcut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'recrcut','LEN',0)

 intarr(1,:)=dtsets(:)%rectesteg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rectesteg','INT',0)

 dprarr(1,:)=dtsets(:)%rectolden
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'rectolden','DPR',0)

 intarr(1,:)=dtsets(:)%restartxf
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'restartxf','INT',0)

 if(response==1)then

   intarr(1,:)=dtsets(:)%rfasr
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfasr','INT',0)

   intarr(1,:)=dtsets(:)%rfatpol(1)
   intarr(2,:)=dtsets(:)%rfatpol(2)
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'rfatpol','INT',0)

   intarr(1,:)=dtsets(:)%rfdir(1)
   intarr(2,:)=dtsets(:)%rfdir(2)
   intarr(3,:)=dtsets(:)%rfdir(3)
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'rfdir','INT',0)

   intarr(1,:)=dtsets(:)%rfddk
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfddk','INT',0)

   intarr(1,:)=dtsets(:)%rfelfd
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfelfd','INT',0)

   intarr(1,:)=dtsets(:)%rfmeth
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfmeth','INT',0)

   intarr(1,:)=dtsets(:)%rfphon
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfphon','INT',0)

   intarr(1,:)=dtsets(:)%rfstrs
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfstrs','INT',0)

   intarr(1,:)=dtsets(:)%rfuser
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rfuser','INT',0)

 end if

 intarr(1,:)=dtsets(:)%rf1atpol(1)
 intarr(2,:)=dtsets(:)%rf1atpol(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'rf1atpol','INT',0)
 intarr(1,:)=dtsets(:)%rf1dir(1)
 intarr(2,:)=dtsets(:)%rf1dir(2)
 intarr(3,:)=dtsets(:)%rf1dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'rf1dir','INT',0)
 intarr(1,:)=dtsets(:)%rf1elfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rf1elfd','INT',0)
 intarr(1,:)=dtsets(:)%rf1phon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rf1phon','INT',0)

 intarr(1,:)=dtsets(:)%rf2atpol(1)
 intarr(2,:)=dtsets(:)%rf2atpol(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'rf2atpol','INT',0)
 intarr(1,:)=dtsets(:)%rf2dir(1)
 intarr(2,:)=dtsets(:)%rf2dir(2)
 intarr(3,:)=dtsets(:)%rf2dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'rf2dir','INT',0)
 intarr(1,:)=dtsets(:)%rf2elfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rf2elfd','INT',0)
 intarr(1,:)=dtsets(:)%rf2phon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rf2phon','INT',0)

 intarr(1,:)=dtsets(:)%rf3atpol(1)
 intarr(2,:)=dtsets(:)%rf3atpol(2)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,2,narrm,ncid,ndtset_alloc,'rf3atpol','INT',0)
 intarr(1,:)=dtsets(:)%rf3dir(1)
 intarr(2,:)=dtsets(:)%rf3dir(2)
 intarr(3,:)=dtsets(:)%rf3dir(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'rf3dir','INT',0)
 intarr(1,:)=dtsets(:)%rf3elfd
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rf3elfd','INT',0)
 intarr(1,:)=dtsets(:)%rf3phon
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'rf3phon','INT',0)

 prtimg(:,:)=1
 if(mxvals%mxnimage==1)then
   do idtset=0,ndtset_alloc
     dprarr(1:9,idtset)= reshape(results_out(idtset)%rprim(:,:,1),(/9/))
     do ii=1,9
       if(abs(dprarr(ii,idtset))<tol12)dprarr(ii,idtset)=zero  ! This is to improve the portability
     end do
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,9,narrm,ncid,ndtset_alloc,'rprim','DPR',0)
 else

   do idtset=1,ndtset_alloc       ! especific size for each dataset
     nimagem(idtset)=dtsets(idtset)%nimage
     narrm(idtset)=9
     do iimage=1,dtsets(idtset)%nimage
       if (narrm(idtset)>0) then
         dprarr_images(1:narrm(idtset),iimage,idtset)=&
&         reshape(results_out(idtset)%rprim(1:3,1:3,iimage),&
&         (/ narrm(idtset) /) )
       end if
     end do
   end do
   call prttagm_images(dprarr_images,iout,jdtset_,&
&   marr,narrm,ncid,ndtset_alloc,'rprim',&
&   mxvals%mxnimage,nimagem,ndtset,prtimg,strimg)

 end if

 dprarr(1,:)=dtsets(:)%scphon_temp
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'scphon_temp','ENE',0)

 if(response==1)then
   dprarr(1,:)=dtsets(:)%sciss
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'sciss','ENE',0)
 end if

!shiftk (printed only when kptopt>0)
 if(sum((dtsets(1:ndtset_alloc)%kptopt)**2)/=0)then

   multi_kptopt=0
   dprarr(:,0)=0.0_dp
   narr=3*dtsets(1)%nshiftk ! default size for all datasets
   do idtset=1,ndtset_alloc       ! especific size for each dataset
     narrm(idtset)=3*dtsets(idtset)%nshiftk
     if (narrm(idtset)>0) then
       dprarr(1:narrm(idtset),idtset)=&
&       reshape(dtsets(idtset)%shiftk(1:3,1:dtsets(idtset)%nshiftk),&
&       (/ narrm(idtset) /) )
     end if
     if(dtsets(idtset)%kptopt<=0)then
       narrm(idtset)=0
       multi_kptopt=1
     end if
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
&   narrm,ncid,ndtset_alloc,'shiftk','DPR',&
&   multi_nshiftk)


!  End of test to see whether kptopt/=0 for some dataset
 end if

 intarr(1,:)=dtsets(:)%signperm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'signperm','INT',0)

 dprarr(1,:)=dtsets(:)%slabwsrad
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'slabwsrad','DPR',0)

 dprarr(1,:)=dtsets(:)%slabzbeg
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'slabzbeg','DPR',0)

 dprarr(1,:)=dtsets(:)%slabzend
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'slabzend','DPR',0)

 dprarr(1,:)=dtsets(:)%smdelta
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'smdelta','INT',0)

 dprarr(1,:)=dtsets(:)%soenergy
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'soenergy','ENE',0)

 dprarr(1,:)=dtsets(:)%spbroad
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'spbroad','ENE',0)

 do idtset=0,ndtset_alloc
   intarr(1:npsp,idtset)=dtsets(idtset)%so_psp(1:npsp)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,npsp,narrm,ncid,ndtset_alloc,'so_psp','INT',0)

 intarr(1,:)=dtsets(:)%spgroup
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'spgroup','INT',0)

!spinat
 dprarr(:,0)=0.0_dp
 narr=3*natom ! default size for all datasets
 do idtset=1,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%natom
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%spinat(1:3,1:dtsets(idtset)%natom),&
&     (/ narrm(idtset) /) )
   end if
   if(sum(abs( dtsets(idtset)%spinat(1:3,1:dtsets(idtset)%natom))) < tol12 )then
     narrm(idtset)=0
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,narr,&
& narrm,ncid,ndtset_alloc,'spinat','DPR',&
& multi_natom)

 intarr(1,:)=dtsets(:)%spmeth
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'spmeth','INT',0)

 dprarr(1,:)=dtsets(:)%stmbias
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'stmbias','DPR',0)

 dprarr(1,:)=dtsets(:)%strfact
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'strfact','DPR',0)

 do ii=1,6
   dprarr(ii,:)=dtsets(:)%strtarget(ii)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,6,narrm,ncid,ndtset_alloc,'strtarget','DPR',0)

!strten
 if(choice==2)then

   prtimg(:,:)=1
   do idtset=1,ndtset_alloc       ! specific size for each dataset
     compute_static_images=(dtsets(idtset)%istatimg>0)
     nimagem(idtset)=dtsets(idtset)%nimage
     narrm(idtset)=6

     if(dtsets(idtset)%iscf>0)then

       do iimage=1,dtsets(idtset)%nimage

         if (narrm(idtset)>0) then
           dprarr_images(1:narrm(idtset),iimage,idtset)=&
&           results_out(idtset)%strten(:,iimage)
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
&   marr,narrm,ncid,ndtset_alloc,'strten',&
&   mxvals%mxnimage,nimagem,ndtset,prtimg,strimg)

 end if

 intarr(1,:)=dtsets(:)%scphon_supercell(1)
 intarr(2,:)=dtsets(:)%scphon_supercell(2)
 intarr(3,:)=dtsets(:)%scphon_supercell(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'scphon_supercell','INT',0)

 intarr(1,:)=dtsets(:)%supercell(1)
 intarr(2,:)=dtsets(:)%supercell(2)
 intarr(3,:)=dtsets(:)%supercell(3)
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'supercell','INT',0)

!symafm
 intarr(:,0)=1
 narr=nsym ! default size for all datasets
 do idtset=1,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=dtsets(idtset)%nsym
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=&
&     dtsets(idtset)%symafm(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'symafm','INT',&
& multi_nsym)

 intarr(1,:)=dtsets(:)%symmorphi
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'symmorphi','INT',0)

!symrel
 intarr(1:9,0)=(/ 1,0,0, 0,1,0, 0,0,1 /)
 narr=9*nsym ! default size for all datasets
 do idtset=1,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=9*dtsets(idtset)%nsym
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%symrel(1:3,1:3,1:dtsets(idtset)%nsym),&
&     (/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,3,marr,narr,&
& narrm,ncid,ndtset_alloc,'symrel','INT',&
& multi_nsym)


 intarr(1,:)=dtsets(:)%symsigma
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'symsigma','INT',0)

 intarr(1,:)=dtsets(:)%symchi
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'symchi','INT',0)

 dprarr(1,:)=dtsets(:)%td_maxene
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'td_maxene','DPR',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%td_mexcit
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'td_mexcit','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%tfkinfunc
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'tfkinfunc','INT',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%recgratio
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'recgratio','INT',0)

!timopt
 timopt_default=1
!MPI parallel case
 if(mpi_enreg%paral_compil==1)then
   timopt_default=0
 end if
 if(timopt/=timopt_default)then
   intarr(1,:)=timopt
   intarr(1,0)=timopt_default
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'timopt','INT',0)
 end if

!WVL - tails related variables
 intarr(1,:)=dtsets(:)%tl_nprccg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'tl_nprccg','INT',0)
 dprarr(1,:)=dtsets(:)%tl_radius
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tl_radius','DPR',0)

!tnons
 dprarr(:,0)=0.0_dp
 narr=3*nsym ! default size for all datasets
 do idtset=1,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%nsym
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%tnons(1:3,1:dtsets(idtset)%nsym),&
&     (/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,-3,marr,narr,&
& narrm,ncid,ndtset_alloc,'tnons','DPR',&
& multi_nsym)

 dprarr(1,:)=dtsets(:)%tolimg
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolimg','ENE',0)

 dprarr(1,:)=dtsets(:)%toldfe
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'toldfe','ENE',0)

 dprarr(1,:)=dtsets(:)%toldff
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'toldff','DPR',0)

 dprarr(1,:)=dtsets(:)%tolrff
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolrff','DPR',0)

 dprarr(1,:)=dtsets(:)%tolmxf
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolmxf','DPR',0)

 dprarr(1,:)=dtsets(:)%tolsym
!DEBUG
!write(std_out,*)' tolsym=',dtsets(:)%tolsym
!ENDDEBUG
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolsym','DPR',0)

 dprarr(1,:)=dtsets(:)%tolvrs
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolvrs','DPR',0)

 dprarr(1,:)=dtsets(:)%tolwfr
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tolwfr','DPR',0)

 dprarr(1,:)=dtsets(:)%tphysel
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tphysel','ENE',0)

 dprarr(1,:)=dtsets(:)%tsmear
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'tsmear','ENE',0)

!typat
 narr=natom                      ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=dtsets(idtset)%natom
   if (narrm(idtset)>0) then
     intarr(1:narrm(idtset),idtset)=dtsets(idtset)%typat(1:narrm(idtset))
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,4,marr,narr,&
& narrm,ncid,ndtset_alloc,'typat','INT',multi_natom)

 intarr(1,:)=dtsets(:)%use_gpu_cuda
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'use_gpu_cuda','INT',0)

 intarr(1,:)=dtsets(:)%use_slk
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'use_slk','INT',0)

 intarr(1,:)=dtsets(:)%usekden
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usekden','INT',0)

 intarr(1,:)=dtsets(:)%useria
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'useria','INT',0)

 intarr(1,:)=dtsets(:)%userib
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'userib','INT',0)

 intarr(1,:)=dtsets(:)%useric
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'useric','INT',0)

 intarr(1,:)=dtsets(:)%userid
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'userid','INT',0)

 intarr(1,:)=dtsets(:)%userie
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'userie','INT',0)

 dprarr(1,:)=dtsets(:)%userra
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'userra','DPR',0)

 dprarr(1,:)=dtsets(:)%userrb
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'userrb','DPR',0)

 dprarr(1,:)=dtsets(:)%userrc
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'userrc','DPR',0)

 dprarr(1,:)=dtsets(:)%userrd
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'userrd','DPR',0)

 dprarr(1,:)=dtsets(:)%userre
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'userre','DPR',0)

 intarr(1,:)=dtsets(:)%usewvl
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usewvl','INT',0)

 if (usepaw==1) then
   intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%usexcnhat
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'usexcnhat','INT',0)
 end if

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%useylm
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'useylm','INT',0)

 dprarr(1,:)=dtsets(:)%vcutgeo(1)
 dprarr(2,:)=dtsets(:)%vcutgeo(2)
 dprarr(3,:)=dtsets(:)%vcutgeo(3)
 call prttagm(dprarr,intarr,iout,jdtset_,3,marr,3,narrm,ncid,ndtset_alloc,'vcutgeo','DPR',0)

!vel
 prtimg(:,:)=1
 if(multi_natom==0 .and. mxvals%mxnimage==1)then
   do idtset=0,ndtset_alloc
     dprarr(1:3*natom,idtset)=reshape(results_out(idtset)%vel(:,1:natom,1),(/3*natom/) )
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3*natom,&
&   narrm,ncid,ndtset_alloc,'vel','DPR',0)
 else
   do idtset=1,ndtset_alloc       ! especific size for each dataset
     nimagem(idtset)=dtsets(idtset)%nimage
     narrm(idtset)=3*dtsets(idtset)%natom
     do iimage=1,dtsets(idtset)%nimage
       if (narrm(idtset)>0) then
         dprarr_images(1:narrm(idtset),iimage,idtset)=&
&         reshape(results_out(idtset)%vel(1:3,1:dtsets(idtset)%natom,iimage),&
&         (/ narrm(idtset) /) )
       end if
       if(sum(abs( results_out(idtset)%vel(:,1:dtsets(idtset)%natom,iimage)- &
&       results_out(0)%vel(:,1:dtsets(idtset)%natom,iimage)      )) < tol12 )then
         prtimg(iimage,idtset)=0
       end if
     end do
   end do
   call prttagm_images(dprarr_images,iout,jdtset_,&
&   marr,narrm,ncid,ndtset_alloc,'vel',&
&   mxvals%mxnimage,nimagem,ndtset,prtimg,strimg)

 end if

 dprarr(1,:)=dtsets(:)%vis
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'vis','DPR',0)

 dprarr(1,:)=dtsets(:)%vprtrb(1)
 dprarr(2,:)=dtsets(:)%vprtrb(2)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,2,narrm,ncid,ndtset_alloc,'vprtrb','ENE',0)

 intarr(1,0:ndtset_alloc)=dtsets(0:ndtset_alloc)%wfoptalg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'wfoptalg','INT',0)

!wtatcon
 narr=3*natom*dtsets(1)%nconeq ! default size for all datasets
 do idtset=0,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=3*dtsets(idtset)%natom*dtsets(idtset)%nconeq
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=&
&     reshape(dtsets(idtset)%wtatcon(1:3,1:dtsets(idtset)%natom,&
&     1:dtsets(idtset)%nconeq),&
&     (/ narrm(idtset) /) )
   end if
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
& narrm,ncid,ndtset_alloc,'wtatcon','DPR',&
& multi_natom+multi_nconeq)

!wtk
 tnkpt=0
 dprarr(:,0)=1
 narr=dtsets(1)%nkpt            ! default size for all datasets
 if(prtvol_glob==0 .and. narr>nkpt_max)then
   narr=nkpt_max
   tnkpt=1
 end if

 do idtset=1,ndtset_alloc       ! especific size for each dataset
   narrm(idtset)=dtsets(idtset)%nkpt
   if (narrm(idtset)>0) then
     dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%wtk(1:narrm(idtset))+tol12
   end if

   if(prtvol_glob==0 .and. narrm(idtset)>nkpt_max)then
     narrm(idtset)=nkpt_max
     tnkpt=1
   end if

 end do
 call prttagm(dprarr,intarr,iout,jdtset_,4,marr,narr,&
& narrm,ncid,ndtset_alloc,'wtk','DPR',multi_nkpt)

 if(tnkpt==1) write(iout,'(23x,a,i3,a)' ) &
& 'outvars : Printing only first ',nkpt_max,' k-points.'

!WVL - wavelets variables
 dprarr(1,:)=dtsets(:)%wvl_crmult
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'wvl_crmult','DPR',0)
 dprarr(1,:)=dtsets(:)%wvl_frmult
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'wvl_frmult','DPR',0)
 dprarr(1,:)=dtsets(:)%wvl_hgrid
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'wvl_hgrid','DPR',0)
 intarr(1,:)=dtsets(:)%wvl_nprccg
 call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'wvl_nprccg','INT',0)

!Wannier90 interface related variables
 if(sum(dtsets(1:ndtset_alloc)%prtwant) >1)then
   intarr(1,:)=dtsets(:)%w90iniprj
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'w90iniprj','INT',0)
   intarr(1,:)=dtsets(:)%w90prtunk
   call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'w90prtunk','INT',0)
!  van der Waals correction with MLWFs related variables
   if(any(dtsets(1:ndtset_alloc)%vdw_xc==10))then
     intarr(1,:)=dtsets(:)%vdw_nfrag
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,1,narrm,ncid,ndtset_alloc,'vdw_nfrag','INT',0)

     intarr(1,:)=dtsets(:)%vdw_supercell(1)
     intarr(2,:)=dtsets(:)%vdw_supercell(2)
     intarr(3,:)=dtsets(:)%vdw_supercell(3)
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,3,narrm,ncid,ndtset_alloc,'vdw_supercell','INT',0)
     do iat=1,mxvals%mxnatom
       intarr(iat,:)=dtsets(:)%vdw_typfrag(iat)
     end do
     call prttagm(dprarr,intarr,iout,jdtset_,2,marr,mxvals%mxnatom,narrm,ncid,ndtset_alloc,'vdw_supercell','INT',0)
   end if !vdw_xc==10
 end if !prtwant>1

!xangst
 prtimg(:,:)=1
 if(multi_natom==0.and.mxvals%mxnimage==1)then
   dprarr(1:3*natom,0:ndtset_alloc)=&
&   reshape(xangst_(1:3,1:natom,1,0:ndtset_alloc),(/3*natom,ndtset_alloc+1/) )
   call prttagm(dprarr,intarr,iout,jdtset_,-2,marr,3*natom,narrm,ncid,ndtset_alloc,'xangst','DPR',0)
 else

   do idtset=1,ndtset_alloc       ! especific size for each dataset
     nimagem(idtset)=dtsets(idtset)%nimage
     narrm(idtset)=3*dtsets(idtset)%natom
     do iimage=1,dtsets(idtset)%nimage
       if (narrm(idtset)>0) then
         dprarr_images(1:narrm(idtset),iimage,idtset)=&
&         reshape(xangst_(1:3,1:dtsets(idtset)%natom,iimage,idtset),&
&         (/ narrm(idtset) /) )
       end if
     end do
   end do
   call prttagm_images(dprarr_images,iout,jdtset_,&
&   marr,narrm,ncid,ndtset_alloc,'xangst',&
&   mxvals%mxnimage,nimagem,ndtset,prtimg,strimg)

 end if

!xcart
 if(multi_natom==0.and.mxvals%mxnimage==1)then
   dprarr(1:3*natom,0:ndtset_alloc)=&
&   reshape(xcart_(1:3,1:natom,1,0:ndtset_alloc),(/3*natom,ndtset_alloc+1/) )
   call prttagm(dprarr,intarr,iout,jdtset_,-2,marr,3*natom,&
&   narrm,ncid,ndtset_alloc,'xcart','DPR',0)
 else

   do idtset=1,ndtset_alloc       ! especific size for each dataset
     nimagem(idtset)=dtsets(idtset)%nimage
     narrm(idtset)=3*dtsets(idtset)%natom
     do iimage=1,dtsets(idtset)%nimage
       if (narrm(idtset)>0) then
         dprarr_images(1:narrm(idtset),iimage,idtset)=&
&         reshape(xcart_(1:3,1:dtsets(idtset)%natom,iimage,idtset),&
&         (/ narrm(idtset) /) )
       end if
     end do
   end do
   call prttagm_images(dprarr_images,iout,jdtset_,&
&   marr,narrm,ncid,ndtset_alloc,'xcart',&
&   mxvals%mxnimage,nimagem,ndtset,prtimg,strimg)

 end if

 dprarr(1,:)=dtsets(:)%xc_denpos
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'xc_denpos','DPR',0)

!xred
 if(multi_natom==0.and.mxvals%mxnimage==1)then
   do idtset=0,ndtset_alloc
     dprarr(1:3*natom,idtset)=&
&     reshape(results_out(idtset)%xred(:,1:natom,1),(/3*natom/) )
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,-2,marr,3*natom,&
&   narrm,ncid,ndtset_alloc,'xred','DPR',0)
 else

   do idtset=1,ndtset_alloc       ! especific size for each dataset
     nimagem(idtset)=dtsets(idtset)%nimage
     narrm(idtset)=3*dtsets(idtset)%natom
     do iimage=1,dtsets(idtset)%nimage
       if (narrm(idtset)>0) then
         dprarr_images(1:narrm(idtset),iimage,idtset)=&
&         reshape(results_out(idtset)%xred(1:3,1:dtsets(idtset)%natom,iimage),&
&         (/ narrm(idtset) /) )
       end if
     end do
   end do
   call prttagm_images(dprarr_images,iout,jdtset_,&
&   marr,narrm,ncid,ndtset_alloc,'xred',&
&   mxvals%mxnimage,nimagem,ndtset,prtimg,strimg)

 end if

 dprarr(1,:)=dtsets(:)%zcut
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,1,narrm,ncid,ndtset_alloc,'zcut','ENE',0)

!zeemanfield
 dprarr(1,:)=dtsets(:)%zeemanfield(1)
 dprarr(2,:)=dtsets(:)%zeemanfield(2)
 dprarr(3,:)=dtsets(:)%zeemanfield(3)
 call prttagm(dprarr,intarr,iout,jdtset_,1,marr,3,narrm,ncid,ndtset_alloc,'zeemanfield','BFI',0)

!ziontypat   ! After all, should always echo this value
 if(sum(dtsets(:)%ntypalch)>0)then   ! After all, should always echo this value ...


   narr=ntypat                      ! default size for all datasets
   do idtset=0,ndtset_alloc       ! especific size for each dataset
     narrm(idtset)=dtsets(idtset)%ntypat
     if (narrm(idtset)>0) then
       dprarr(1:narrm(idtset),idtset)=dtsets(idtset)%ziontypat(1:narrm(idtset))
     end if
   end do
   call prttagm(dprarr,intarr,iout,jdtset_,1,marr,narr,&
&   narrm,ncid,ndtset_alloc,'ziontypat','DPR',multi_ntypat)

 end if

 do idtset=0,ndtset_alloc
   dprarr(1:npsp,idtset)=dtsets(idtset)%znucl(1:npsp)
 end do
 call prttagm(dprarr,intarr,iout,jdtset_,4,marr,npsp,narrm,ncid,ndtset_alloc,'znucl','DPR',0)

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(narrm)
 ABI_DEALLOCATE(nimagem)
 ABI_DEALLOCATE(dprarr_images)
 ABI_DEALLOCATE(jdtset_)
 ABI_DEALLOCATE(response_)
 ABI_DEALLOCATE(xangst_)
 ABI_DEALLOCATE(xcart_)
 ABI_DEALLOCATE(strimg)

 write(message,'(a,80a)')ch10,('=',mu=1,80)
 call wrtout(iout,message,'COLL')

!write(std_out,*) 'outvar 06'
!###########################################################
!### 06. Close NetCDF file

#if defined HAVE_TRIO_ETSF_IO
 call etsf_io_low_close(abs(ncid), lstat)
#endif

!! #if defined HAVE_TRIO_NETCDF
!!  ncerr = nf90_close(ncid)
!! #endif

!**************************************************************************

!DEBUG
!write(std_out,*)' outvars : end of subroutine '
!if(.true.)stop
!ENDDEBUG

end subroutine outvars
!!***
