!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars7w
!!
!! NAME
!! invars7w
!!
!! FUNCTION
!! Open input file for the lwf code, then
!! reads or echoes the input information.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! lenstr=actual length of string
!! mqpt=maximum number of q points.
!! natom=number of atoms, needed for xred
!! string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!! grdsize(3)= size of the grid of q points = limit of the shells of the LWF
!! All the other arguments are outputs
!! and are read from the input file.
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      lwf
!!
!! CHILDREN
!!      intagm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine invars7w(allerr,alpha,decflg,enwdmax,enwdmin,frozflg,grdsize,ingss,irwfl,lenstr,localqmode,&
& mqpt,natom,nstom,nwnn,prtvol,rcenter,string,subwdmax,subwdmin,tolomi,trialq,znucl)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invars7w'
 use interfaces_42_parser
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lenstr,mqpt,natom,nwnn
 integer,intent(out) :: decflg,frozflg,irwfl,nstom,prtvol,trialq
 real(dp),intent(out) :: allerr,alpha,enwdmax,enwdmin,subwdmax,subwdmin,tolomi
 character(len=*) :: string
!arrays
 integer,intent(out) :: grdsize(3)
 real(dp),intent(out) :: ingss(nwnn,natom,3,2),localqmode(3),rcenter(3)
 real(dp),intent(out) :: znucl(natom)

!Local variables -------------------------
!Dummy arguments for subroutine 'intagm' to parse input file
!scalars
 integer :: jdtset,marr,mm,tao,tread
 character(len=30) :: token
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

!*********************************************************************

!DEBUG
!write(std_out,*)
!write(std_out,*) ' invars : enter'
!write(std_out,*) ' natom=',natom
!write(std_out,*) ' nwnn=',nwnn
!ENDDEBUG

 marr=4*mqpt
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 jdtset=1

 allerr=zero
 token = 'allerr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) allerr=dprarr(1)

 write(std_out,*) 'allerr=',allerr

 alpha=zero
 token = 'alpha'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) alpha=dprarr(1)

 write(std_out,*) 'alpha=',alpha

 decflg=one
 token = 'decflg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) decflg=intarr(1)

 write(std_out,*) 'decflg=',decflg

 frozflg=one
 token = 'frozflg'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) frozflg=intarr(1)

 write(std_out,*) 'frozflg=',frozflg

 enwdmax=zero
 token = 'enwdmax'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) enwdmax=dprarr(1)

 write(std_out,*) 'enwdmax=',enwdmax

 enwdmin=zero
 token = 'enwdmin'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) enwdmin=dprarr(1)

 write(std_out,*) 'enwdmin=',enwdmin

 grdsize(:)=one
 token = 'grdsize'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
 if(tread==1) grdsize=intarr(1:3)

 write(std_out,*) 'grdsize=',grdsize(:)

 irwfl=one
 token = 'irwfl'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) irwfl=intarr(1)

 write(std_out,*) 'irwfl=',irwfl

 localqmode(:)=zero
 token = 'localqmode'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
 if(tread==1) localqmode=dprarr(1:3)

 write(std_out,*) 'localqmode=',localqmode

 nstom=zero
 token = 'nstom'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) nstom=intarr(1)

 write(std_out,*) 'nstom=',nstom

 prtvol=zero
 token = 'prtvol'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) prtvol=intarr(1)

 write(std_out,*) 'prtvol=',prtvol

 rcenter(:)=zero
 token = 'rcenter'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'DPR')
 if(tread==1) rcenter=dprarr(1:3)

 write(std_out,*) 'rcenter=',rcenter(:)

 subwdmax=zero
 token = 'subwdmax'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) subwdmax=dprarr(1)

 write(std_out,*) 'subwdmax=',subwdmax

 subwdmin=zero
 token = 'subwdmin'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) subwdmin=dprarr(1)

 write(std_out,*) 'subwdmin=',subwdmin

 tolomi=zero
 token = 'tolomi'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
 if(tread==1) tolomi=dprarr(1)

 write(std_out,*) 'tolomi=',tolomi

 trialq=zero
 token = 'trialq'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) trialq=intarr(1)
 if (trialq==1) then
   ingss(:,:,:,:)=0
   token = 'ingss'
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) ingss=reshape(dprarr(1:nwnn*natom*3*2),(/nwnn,natom,3,2/))
 end if

 znucl(:)=zero
 token = 'znucl'
 call intagm(dprarr,intarr,jdtset,marr,natom,string(1:lenstr),token,tread,'DPR')
 if(tread==1) znucl(:)=dprarr(1:natom)

 ingss=zero
 token = 'ingss'
 call intagm(dprarr,intarr,jdtset,marr,nwnn*natom*3*2,string(1:lenstr),token,tread,'DPR')
 if(tread==1) then
   do mm=1,nwnn
     do tao=1,natom
       ingss(mm,tao,:,1)=dprarr((mm-1)*natom*3*2+(tao-1)*3*2+1:(mm-1)*natom*3*2+(tao-1)*3*2+3)
       ingss(mm,tao,:,2)=dprarr((mm-1)*natom*3*2+(tao-1)*3*2+4:(mm-1)*natom*3*2+(tao-1)*3*2+6)
     end do
   end do
 end if


!do mm=1,nwnn
!do tao=1,natom
!write(std_out,'(a,4i4)') 'mm,tao',mm,tao
!write(std_out,'(a,3f12.7)') 'real ing',ingss(mm,tao,:,1)
!write(std_out,'(a,3f12.7)') 'imag ing',ingss(mm,tao,:,2)
!end do
!end do

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

!DEBUG
!write(std_out,*) ' invars : exit'
!ENDDEBUG

end subroutine invars7w
!!***
