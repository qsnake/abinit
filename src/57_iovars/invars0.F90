!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars0
!! NAME
!! invars0
!!
!! FUNCTION
!! Initialisation phase : prepare the main input subroutine call by
!! reading most of the NO MULTI variables, as well as natom, nimage, and ntypat,
!! needed for allocating some input arrays in abinit, and also useri
!! and userr. The variable usewvl is also read here for later reading
!! of input path for the atomic orbital file (if required).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lenstr=actual length of string
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!               one data set.
!!  string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here :
!!   cpus,jdtset,natom,nimage,npsp,ntypat,useri*,userr*
!!  istatr=repetition rate for status file
!!  istatshft=shift of the repetition rate for status file
!!  msym=maximal value of input msym for all the datasets
!!  mxnatom=maximal value of input natom for all the datasets
!!  mxnimage=maximal value of input nimage for all the datasets
!!  mxntypat=maximal value of input ntypat for all the datasets
!!  npsp=number of pseudopotentials
!!
!! PARENTS
!!      m_ab6_invars_f90
!!
!! CHILDREN
!!      get_ndevice,intagm,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine invars0(dtsets,istatr,istatshft,lenstr,&
& msym,mxnatom,mxnimage,mxntypat,ndtset,ndtset_alloc,npsp,papiopt,timopt,string)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
#if defined HAVE_GPU_CUDA
 use m_initcuda, only : Get_ndevice
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invars0'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_parser
#if defined HAVE_GPU_CUDA
#endif
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lenstr,ndtset,ndtset_alloc
 integer,intent(out) :: istatr,istatshft,msym,mxnatom,mxnimage,mxntypat,npsp,papiopt
 integer,intent(inout) :: timopt
 character(len=*),intent(in) :: string
!arrays
 type(dataset_type),intent(out) :: dtsets(0:ndtset_alloc)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,idtset,ii,jdtset,marr,tjdtset,tread,treadh,treadm
 integer :: treads,use_gpu_cuda
 real(dp) :: cpus
 character(len=30) :: token
 character(len=500) :: message
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

!******************************************************************

!Set ii to avoid warning of uninitialised variable
 ii = 0

 marr=max(ndtset_alloc,2)
 ABI_ALLOCATE(dprarr,(marr))
 ABI_ALLOCATE(intarr,(marr))

!Set up jdtset
 if(ndtset/=0)then

!  Default values
   dtsets(0)%jdtset = -1 ! unused value
   dtsets(1:ndtset_alloc)%jdtset=(/ (ii,ii=1,ndtset_alloc) /)

!  Read explicitely the jdtset array
   token = 'jdtset'
   call intagm(dprarr,intarr,0,marr,ndtset,string(1:lenstr),&
&   token,tjdtset,'INT')
   if(tjdtset==1) dtsets(1:ndtset)%jdtset=intarr(1:ndtset)

!  Read the udtset array
   token = 'udtset'
   call intagm(dprarr,intarr,0,marr,2,string(1:lenstr),&
&   token,tread,'INT')

!  jdtset and udtset cannot be defined together
   if(tjdtset==1 .and. tread==1)then
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     ' invars0: ERROR -',ch10,&
&     '  jdtset and udtset cannot be defined both in the input file.',ch10,&
&     '  Action : remove one of them from your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Check values of udtset
   if(tread==1)then
     if(intarr(1)<1 .or. intarr(1)>999)then
       write(message, '(a,a,a,a,i5,a,a,a)' )ch10,&
&       ' invars0: ERROR -',ch10,&
&       '  udtset(1) must be between 1 and 999, but it is',intarr(1),'.',ch10,&
&       '  Action : change the value of udtset(1) in your input file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if(intarr(2)<1 .or. intarr(2)>9)then
       write(message, '(a,a,a,a,i5,a,a,a)' )ch10,&
&       ' invars0: ERROR -',ch10,&
&       '  udtset(2) must be between 1 and 9, but it is',intarr(2),'.',ch10,&
&       '  Action : change the value of udtset(2) in your input file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     if(intarr(1)*intarr(2) /= ndtset)then
       write(message, '(a,a,a,a,a,a,i2,a,a,a,i2,a,i3,a,a,a,i3,a,a,a)' )ch10,&
&       ' invars0: ERROR -',ch10,&
&       '  udtset(1)*udtset(2) must be equal to ndtset,',ch10,&
&       '  but it is observed that udtset(1)=',intarr(1),',',ch10,&
&       '  and udtset(2)=',intarr(2),' so that their product is',&
&       intarr(1)*intarr(2),',',ch10,&
&       '  while ndtset is',ndtset,'.',ch10,&
&       '  Action : change udtset or ndtset in your input file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     idtset=0
     do i1=1,intarr(1)
       do i2=1,intarr(2)
         idtset=idtset+1
         dtsets(idtset)%jdtset=i1*10+i2
       end do
     end do
   end if

!  Final check on the jdtset values
   do idtset=1,ndtset
     if(dtsets(idtset)%jdtset<1 .or. dtsets(idtset)%jdtset>9999)then
       write(message, '(a,a,a,a,a,a,i3,a,i3,a,a)' )ch10, &
&       ' invars0: ERROR -',ch10,&
&       '  The components of jdtset must be between 1 and 9999.',ch10,&
&       '  However, the input value of the component',idtset,&
&       ' of jdtset is',dtsets(idtset)%jdtset,ch10,&
&       '  Action : correct jdtset in your input file.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end do

 else
   dtsets(1)%jdtset=0
 end if

 papiopt = 0
 token = 'papiopt'
 call intagm(dprarr,intarr,0,1,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) papiopt=intarr(1)

!Read timopt and pass it to timab
 token = 'timopt'
 call intagm(dprarr,intarr,0,1,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) timopt=intarr(1)

 istatr=49
 token = 'istatr'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) istatr=intarr(1)

 istatshft=1
 token = 'istatshft'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) istatshft=intarr(1)

 cpus=0.0_dp
 token = 'cpus '
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,treads,'DPR')
 if(treads==1) cpus=dprarr(1)
 token = 'cpum '
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,treadm,'DPR')
 if(treadm==1) cpus=dprarr(1)*60.0_dp
 token = 'cpuh '
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,treadh,'DPR')
 if(treadh==1) cpus=dprarr(1)*3600.0_dp
 if(treads+treadm+treadh>1)then
   write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&   ' invars0: ERROR -',ch10,&
&   '  More than one input variable is used to defined the CPU time limit.',&
&   ch10,'  This is not allowed.',ch10,&
&   '  Action : in the input file, suppress either cpus, cpum or cpuh.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 dtsets(:)%cpus=cpus

!Default for natom, nimage, ntypat, useri and userr
 dtsets(:)%natom=1
 dtsets(:)%nimage=1
 dtsets(:)%ntypat=1 ; dtsets(0)%ntypat=0    ! Will always echo ntypat
 dtsets(:)%macro_uj=0
 dtsets(:)%maxnsym=384
 dtsets(:)%useria=0
 dtsets(:)%userib=0
 dtsets(:)%useric=0
 dtsets(:)%userid=0
 dtsets(:)%userie=0
 dtsets(:)%userra=zero
 dtsets(:)%userrb=zero
 dtsets(:)%userrc=zero
 dtsets(:)%userrd=zero
 dtsets(:)%userre=zero
 dtsets(:)%usewvl = 0

!Loop on datasets, to find natom and mxnatom, as well as useri and userr
 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0

!  Read natom from string
   token = 'natom'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
!  Might initialize natom from CML file
   if(tread==0)then
     token = '_natom'
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   end if

   if(tread==1)then
     dtsets(idtset)%natom=intarr(1)
   else
     write(message, '(a,a,a,a,i6,a,a)' ) ch10,&
&     ' invars0: ERROR -',ch10,&
&     '  Input natom must be defined, but was absent for dataset',jdtset,ch10,&
&     '  Action : check the input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
!  Check that natom is greater than 0
   if (dtsets(idtset)%natom<=0) then
     write(message, '(a,a,a,a,i12,a,a,i6,a,a,a)' ) ch10,&
&     ' invars0: ERROR -',ch10,&
&     '  Input natom must be > 0, but was ',dtsets(idtset)%natom,ch10,&
&     '  for dataset',jdtset,'. This is not allowed.',ch10,&
&     '  Action : check the input file.'
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if

   token = 'nimage'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtsets(idtset)%nimage=intarr(1)

!  Check that nimage is greater than 0
   if (dtsets(idtset)%nimage<=0) then
     write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&     ' invars0: ERROR -',ch10,&
&     '  nimage must be > 0, but was ',dtsets(idtset)%nimage,ch10,&
&     '  This is not allowed.',ch10,&
&     '  Action : check the input file.'
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if

   token = 'ntypat'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1)dtsets(idtset)%ntypat=intarr(1)
!  Check that ntypat is greater than 0
   if (dtsets(idtset)%ntypat<=0) then
     write(message, '(a,a,a,a,i12,a,a,i6,a,a,a)' ) ch10,&
&     ' invars0: ERROR -',ch10,&
&     '  Input ntypat must be > 0, but was ',dtsets(idtset)%ntypat,ch10,&
&     '  for dataset',jdtset,'. This is not allowed.',ch10,&
&     '  Action : check the input file.'
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if

   token = 'ntype'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) then
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     ' invars0: ERROR -',ch10,&
&     '  The use of the "ntype" input variable is forbidden since version 4.1 .',ch10,&
&     '  Action : replace "ntype" by "ntypat".'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Read msym from string
   token = 'maxnsym'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1)dtsets(idtset)%maxnsym=intarr(1)
!  Check that maxnsym is greater than 1
   if (dtsets(idtset)%maxnsym<1) then
     write(message, '(a,a,a,a,i12,a,a,i6,a,a,a)' ) ch10,&
&     ' invars0: ERROR -',ch10,&
&     '  Input maxnsym must be > 1, but was ',dtsets(idtset)%maxnsym,ch10,&
&     '  for dataset',jdtset,'. This is not allowed.',ch10,&
&     '  Action : check the input file.'
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if

   token = 'useria'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtsets(idtset)%useria=intarr(1)
   token = 'userib'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtsets(idtset)%userib=intarr(1)
   token = 'useric'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtsets(idtset)%useric=intarr(1)
   token = 'userid'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtsets(idtset)%userid=intarr(1)
   token = 'userie'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtsets(idtset)%userie=intarr(1)

   token = 'userra'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtsets(idtset)%userra=dprarr(1)
   token = 'userrb'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtsets(idtset)%userrb=dprarr(1)
   token = 'userrc'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtsets(idtset)%userrc=dprarr(1)
   token = 'userrd'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtsets(idtset)%userrd=dprarr(1)
   token = 'userre'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtsets(idtset)%userre=dprarr(1)

   token = 'usewvl'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtsets(idtset)%usewvl=intarr(1)

 end do

!mxnatom =maxval(dtsets(1:ndtset_alloc)%natom)
!mxntypat =maxval(dtsets(1:ndtset_alloc)%ntypat)
!msym =maxval(dtsets(1:ndtset_alloc)%maxnsym)
!There is a bug in the HP compiler, the following should execute properly
 mxnatom=dtsets(1)%natom ; mxnimage=dtsets(1)%nimage
 mxntypat=dtsets(1)%ntypat ; msym=dtsets(1)%maxnsym
 if(ndtset_alloc>1)then
   do idtset=2,ndtset_alloc
     mxnatom =max(dtsets(idtset)%natom,mxnatom)
     mxnimage=max(dtsets(idtset)%nimage,mxnimage)
     mxntypat=max(dtsets(idtset)%ntypat,mxntypat)
     msym    =max(dtsets(idtset)%maxnsym,msym)
   end do
 end if

 if(mxnimage>1)then
   do idtset=2,ndtset_alloc
     if(mxnatom/=dtsets(idtset)%natom)then
       write(message,'(8a,i5,a,i5,3a,i5,a)') ch10,&
&       ' invars0: ERROR -',ch10,&
&       '  When there exist one dataset with more than one image,',ch10,&
&       '  the number of atoms in each dataset must be the same.',ch10,&
&       '  However, it has been found that for dataset=',idtset,ch10,&
&       '  natom=',dtsets(idtset)%natom,' differs from the maximum number',ch10,&
&       '  of atoms, mxnatom=',mxnatom,&
&       '  Action : check the input variables natom for different datasets.'
       call wrtout(std_out,  message,'COLL')
       call leave_new('COLL')
     end if
   end do
 end if

!Set up npsp
 npsp=mxntypat   ! Default value
 token = 'npsp'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1)then
   npsp=intarr(1)
 else
   if(ndtset_alloc>1)then
     do idtset=1,ndtset_alloc
       if(dtsets(idtset)%ntypat/=mxntypat)then
         write(message, '(8a,i6,a,i6,a,a,i6,2a)' ) ch10,&
&         ' invars0: ERROR -',ch10,&
&         '  When npsp is not defined, the input variable ntypat must be',ch10,&
&         '  the same for all datasets. However, it has been found that for',ch10,&
&         '  jdtset:',dtsets(idtset)%jdtset,', ntypat=',dtsets(idtset)%ntypat,ch10,&
&         '  differs from the maximum value of ntypat=',mxntypat,ch10,&
&         '  Action : check the input variables npsp and ntypat.'
         call wrtout(std_out,  message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end if
 end if
 dtsets(0)%npsp = mxntypat   ! Default value
 dtsets(1:ndtset_alloc)%npsp = npsp

!BEGIN TF_CHANGES
 dtsets(0)%paral_rf = 0    ! Default value
 dtsets(0)%ngroup_rf = 1   ! Default value
 do idtset=1,ndtset_alloc
   token = 'paral_rf'
   call intagm(dprarr,intarr,dtsets(idtset)%jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) then
     dtsets(idtset)%paral_rf=intarr(1)
   else
     dtsets(idtset)%paral_rf=0
   end if

   token = 'ngroup_rf'
   call intagm(dprarr,intarr,dtsets(idtset)%jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) then
     dtsets(idtset)%ngroup_rf=intarr(1)
   else
     dtsets(idtset)%ngroup_rf=1
   end if

   if(dtsets(idtset)%ngroup_rf<1 .and. dtsets(idtset)%paral_rf==1) then
     write(message, * ) ch10,&
&     ' invars0: ERROR -',ch10,&
&     '  When paral_rf is 1, ngroup_rf must be greater then 0',ch10,&
&     '  However, it has been found that for',ch10,&
&     '  paral_rf=',dtsets(idtset)%paral_rf,ch10,&
&     '  ngroup_rf=',dtsets(idtset)%ngroup_rf,ch10,&
&     '  Action : check the input variables paral_rf and ngroup_rf.'
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if
 end do

!END TF_CHANGES

!Parallel information
 dtsets(:)%paral_kgb=0
 do idtset=1,ndtset_alloc
   token = 'paral_kgb'
   call intagm(dprarr,intarr,idtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1)dtsets(idtset)%paral_kgb=intarr(1)
 end do

 dtsets(:)%use_slk=0
 do idtset=1,ndtset_alloc
   token = 'use_slk'
   call intagm(dprarr,intarr,idtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1)dtsets(idtset)%use_slk=intarr(1)
 end do

!GPU information
 use_gpu_cuda=0
 dtsets(:)%use_gpu_cuda=0
#if defined HAVE_GPU_CUDA && defined HAVE_GPU_CUDA_DP
 call Get_ndevice(ii)
 if (ii>0) then
   do i1=1,ndtset_alloc
     dtsets(i1)%use_gpu_cuda=-1
   end do
 end if
#endif
 do idtset=1,ndtset_alloc
   token = 'use_gpu_cuda'
   call intagm(dprarr,intarr,idtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1)dtsets(idtset)%use_gpu_cuda=intarr(1)
   if (dtsets(idtset)%use_gpu_cuda==1) use_gpu_cuda=1
 end do
 if (use_gpu_cuda==1) then
#if defined HAVE_GPU_CUDA && defined HAVE_GPU_CUDA_DP
   if (ii<=0) then
     write(message,'(6a)') ch10,&
&     ' invars0: ERROR -',ch10,&
&     '   Input variables use_gpu_cuda is on',ch10,&
&     '   but no available GPU device has been detected !',ch10,&
&     '   Action : change the input variable use_gpu_cuda.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
#else
   write(message,'(10a)') ch10,&
&   ' invars0: ERROR -',ch10,&
&   '   Input variables use_gpu_cuda is on but abinit hasn''t been built',ch10,&
&   '   with (double precision) gpu mode enabled !',ch10,&
&   '   Action : change the input variable use_gpu_cuda',ch10,&
&   '            or re-compile ABINIT with double-precision Cuda enabled.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
#endif
 end if

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)

!We allocate the internal array, depending on the computed values.
!WARNING : do not forget to deallocate these arrays in the routine dtsetfree
!WARNING : do not forget to allocate and deallocate these arrays in the routine newsp
!(should make a separate subroutine for allocating/deallocating these records)
 do idtset=0,ndtset_alloc
   ABI_ALLOCATE(dtsets(idtset)%acell_orig,(3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%algalch,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%amu,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%corecs,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%densty,(mxntypat,4))
   ABI_ALLOCATE(dtsets(idtset)%dynimage,(mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%iatfix,(3,mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%jpawu,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%kberry,(3,20))
   ABI_ALLOCATE(dtsets(idtset)%lexexch,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%lpawu,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%mixalch,(npsp,mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%pimass,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%ptcharge,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%prtatlist,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%quadmom,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%ratsph,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%rprim_orig,(3,3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%rprimd_orig,(3,3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%so_psp,(npsp))
   ABI_ALLOCATE(dtsets(idtset)%spinat,(3,mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%shiftk,(3,8))
   ABI_ALLOCATE(dtsets(idtset)%typat,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%upawu,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%vel_orig,(3,mxnatom,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%xred_orig,(3,mxnatom,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%ziontypat,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%znucl,(npsp))
 end do

!DEBUG
!write(std_out,*)' invars0 : nimage, mxnimage = ',dtsets(:)%nimage, mxnimage
!write(std_out,*)' invars0 : natom = ',dtsets(:)%natom
!write(std_out,*)' invars0 : mxnatom = ',mxnatom
!ENDDEBUG

end subroutine invars0
!!***
