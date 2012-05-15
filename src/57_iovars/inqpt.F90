!{\src2tex{textfont=tt}}
!!****f* ABINIT/inqpt
!!
!! NAME
!! inqpt
!!
!! FUNCTION
!! Initialize the q point
!! for one particular dataset, characterized by jdtset.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!! qptn(3)=reduced coordinates of eventual q point (normalisation is already included)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      invars1
!!
!! CHILDREN
!!      getkgrid,intagm,leave_new,metric,symfind,symlatt,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine inqpt(chksymbreak,iout,jdtset,lenstr,natom,qptn,rprimd,spinat,string,typat,vacuum,xred)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inqpt'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_geometry
 use interfaces_42_parser
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chksymbreak,iout,jdtset,lenstr,natom
 character(len=*),intent(in) :: string
!arrays
 integer,intent(in) :: typat(natom),vacuum(3)
 real(dp),intent(out) :: qptn(3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: spinat(3,natom)
 real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ii,iqpt,iscf_fake,marr,msym,nptsym,nqpt,nqpt_computed,nshiftq,nsym_new,qptopt
 integer :: tread,tread_qptrlatt,tread_ngqpt,use_inversion
 real(dp) :: qptnrm,qptrlen,tolsym,ucvol
 character(len=30) :: token
 character(len=500) :: message

!arrays
 integer :: bravais(11)
 integer :: ngqpt(3)
 integer :: qptrlatt(3,3)
 integer, allocatable :: symafm_new(:) 
 integer, allocatable :: ptsymrel(:,:,:),symrel_new(:,:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),qpt(3),rmet(3,3),shiftq(3,8)
 real(dp),allocatable :: qpts(:,:),tnons_new(:,:),wtq(:)


!no_abirules
!Dummy arguments for subroutine 'intagm' to parse input file
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

! *************************************************************************

!DEBUG
 write(std_out,*)' inqpt : enter'
!write(std_out,*)' inqpt : iout=',iout
!ENDDEBUG

!Compute the maximum size of arrays intarr and dprarr (nshiftq is 8 at maximum)
 marr=24
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!Find the method to generate the k points
 qptopt=0
 token = 'qptopt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1)qptopt=intarr(1)

 if(qptopt==0)then

!  Read qpt and qptnrm
   qpt=zero
   token = 'qpt'
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
   if(tread==1) qpt(1:3)=dprarr(1:3)

   qptnrm=1.0_dp
   token = 'qptnrm'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) qptnrm=dprarr(1)
   if(qptnrm<tol10)then
     write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&     ' inqpt: ERROR -',ch10,&
&     '  The input variable qptnrm is lower than 1.0d-10,',ch10,&
&     '  while it must be a positive, non-zero number.   ',ch10,&
&     '  Action : correct the qptnrm in the input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if

   qptn(:)=qpt(:)/qptnrm

 else if(qptopt>=1 .and. qptopt<=4)then

   ngqpt(:)=0
   token = 'ngqpt'
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),&
&   token,tread_ngqpt,'INT')
   if(tread_ngqpt==1)then
     ngqpt(1:3)=intarr(1:3)
     do ii=1,3
       if(ngqpt(ii)<1)then
         write(message,  '(a,a,a,i1,a,a,a,i4,a,a,a)' ) &
&         ' inqpt : ERROR -',ch10,&
&         '  The input variable ngqpt(',ii,') must be strictly positive,',ch10,&
&         '  while it is found to be',ngqpt(ii),'.',ch10,&
&         '  Action : change it in your input file, or change qptopt.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end if

   token = 'qptrlatt'
   call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),&
&   token,tread_qptrlatt,'INT')
   if(tread_qptrlatt==1)&
&   qptrlatt(:,:)=reshape(intarr(1:9), (/3,3/) )

   if(tread_ngqpt==1 .and. tread_qptrlatt==1)then
     write(message,  '(a,a,a,a,a,a,a)' ) &
&     ' inqpt : ERROR -',ch10,&
&     '  The input variables ngqpt and qptrlatt cannot both ',ch10,&
&     '  be defined in the input file.',ch10,&
&     '  Action : change one of ngqpt or qptrlatt in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   else if(tread_ngqpt==1)then
     qptrlatt(:,:)=0
     qptrlatt(1,1)=ngqpt(1)
     qptrlatt(2,2)=ngqpt(2)
     qptrlatt(3,3)=ngqpt(3)
   end if

   nshiftq=1
   token = 'nshiftq'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),&
&   token,tread,'INT')
   if(tread==1)nshiftq=intarr(1)

   if(nshiftq<1 .or. nshiftq>8 )then
     write(message,  '(6a,i4,a,a,a)' )ch10, &
&     ' inqpt : ERROR -',ch10,&
&     '  The only allowed values of nshiftq are between 1 and 8,',ch10,&
&     '  while it is found to be',nshiftq,'.',ch10,&
&     '  Action : change the value of nshiftq in your input file, or change qptopt.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   shiftq=half
   token = 'shiftq'
   call intagm(dprarr,intarr,jdtset,marr,&
&   3*nshiftq,string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     shiftq(:,1:nshiftq)=reshape( dprarr(1:3*nshiftq), (/3,nshiftq/) )
   else
     if(nshiftq/=1)then
       write(message,  '(6a,i4,a,a)' )ch10, &
&       ' inqpt : ERROR -',ch10,&
&       '  When nshiftq is not equal to 1, shiftq must be defined in the input file.',ch10,&
&       '  However, shiftq is not defined, while nshiftq=',nshiftq,ch10,&
&       '  Action : change the value of nshiftq in your input file, or define shiftq.'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end if

!  Re-generate symmetry operations from the lattice and atomic coordinates
!  This is a fundamental difference with respect ot the k point generation.
   tolsym=tol8
   msym=192
   ABI_ALLOCATE(ptsymrel,(3,3,msym))
   ABI_ALLOCATE(symafm_new,(msym))
   ABI_ALLOCATE(symrel_new,(3,3,msym))
   ABI_ALLOCATE(tnons_new,(3,msym))
   call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
   use_inversion=1
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
   call symfind(0,(/zero,zero,zero/),gprimd,0,msym,natom,0,nptsym,nsym_new,&
&   ptsymrel,spinat,symafm_new,symrel_new,tnons_new,tolsym,typat,use_inversion,xred)

!  Prepare to compute the q-point grid in the ZB or IZB
   iscf_fake=0 ! Do not need the weights

!  Compute the maximum number of q points
   nqpt=0
   ABI_ALLOCATE(qpts,(3,nqpt))
   ABI_ALLOCATE(wtq,(nqpt))
   call getkgrid(chksymbreak,0,iscf_fake,qpts,qptopt,qptrlatt,qptrlen,&
&   msym,nqpt,nqpt_computed,nshiftq,nsym_new,rprimd,&
&   shiftq,symafm_new,symrel_new,vacuum,wtq)
   nqpt=nqpt_computed
   ABI_DEALLOCATE(qpts)
   ABI_DEALLOCATE(wtq)

!  Find the index of the q point within the set of q points that will be generated
   iqpt=0
   token = 'iqpt'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1)iqpt=intarr(1)

!  Checks that iqpt is among the computed q points
   if(iqpt<0)then
     write(message, '(4a,i5,3a)' ) ch10,&
&     ' inqpt: ERROR -',ch10,&
&     '  The input variable iqpt,',iqpt,' is negative, while it should be 0 or positive.',ch10,&
&     '  Action : correct iqpt in the input file.'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if

   if(iqpt>nqpt_computed)then
     write(message, '(4a,i5,3a,i5,7a)' ) ch10,&
&     ' inqpt: ERROR -',ch10,&
&     '  The input variable iqpt,',iqpt,' is bigger than the computed number of q-points in the grid,',ch10,&
&     '  which is ',nqpt,'.',ch10,&
&     '  The latter has been computed from the input variables qptrlatt, ngqpt, nshiftq,',ch10,&
&     '  shiftq, as well as qptopt, the symmetries of the lattice, and spinat.',ch10,&
&     '  Action : correct iqpt in the input file, or correct the computed q-point grid.'
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if

!  Compute the q-point grid in the BZ or the IBZ
   ABI_ALLOCATE(qpts,(3,nqpt))
   ABI_ALLOCATE(wtq,(nqpt))
   call getkgrid(chksymbreak,iout,iscf_fake,qpts,qptopt,qptrlatt,qptrlen,&
&   msym,nqpt,nqpt_computed,nshiftq,nsym_new,rprimd,&
&   shiftq,symafm_new,symrel_new,vacuum,wtq)

!  Transfer to qptn, and deallocate
   qptn(:)=zero
   if(iqpt/=0)then
     qptn(:)=qpts(:,iqpt)
   end if

   ABI_DEALLOCATE(ptsymrel)
   ABI_DEALLOCATE(symafm_new)
   ABI_DEALLOCATE(symrel_new)
   ABI_DEALLOCATE(tnons_new)
   ABI_DEALLOCATE(qpts)
   ABI_DEALLOCATE(wtq)

 else

   write(message,  '(a,a,a,a,a,i4,a,a,a)' ) &
&   ' inqpt : ERROR -',ch10,&
&   '  The only values of qptopt allowed are smaller than 4.',ch10,&
&   '  The input value of qptopt is',qptopt,'.',ch10,&
&   '  Action : change qptopt in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')

 end if

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

!DEBUG
 write(std_out,*)' inqpt : exit '
!call flush(6)
!ENDDEBUG

end subroutine inqpt
!!***
