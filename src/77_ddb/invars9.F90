!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars9
!!
!! NAME
!! invars9
!!
!! FUNCTION
!! Open input file for the anaddb code, then
!! reads or echoes the input information.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,JCC,CL,MVeithen,XW,MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! lenstr=actual length of string
!! natom=number of atoms, needed for atifc
!! qtol=tolerance for wavevector comparison
!! string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!! anaddb_dtset= (derived datatype) contains all the input variables
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! 27/01/2009: MJV: I have cleaned this routine extensively, putting all
!!  variables in alphabetical order, and in a second segment the dependent
!!  variables which need to be allocated depending on the dimensions read in.
!!  Could be divided into two routines as in abinit.
!!    FIXME: move checks to chkin9?
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      intagm,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine invars9 (anaddb_dtset,lenstr,natom,qtol,string)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invars9'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_parser
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: lenstr,natom
 real(dp),intent(in) :: qtol
 character(len=*),intent(in) :: string
 type(anaddb_dataset_type),intent(out) :: anaddb_dtset

!Local variables -------------------------
!Dummy arguments for subroutine 'intagm' to parse input file
!Set routine version number here:
!scalars
 integer,parameter :: vrsddb=100401
 integer :: ii,iph1,iph2,jdtset,marr,tread
 character(len=30) :: token
 character(len=500) :: message
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

!*********************************************************************

!DEBUG
!write(std_out,*)' invars9 : enter '
!ENDDEBUG

 marr=3
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 jdtset=1

!copy natom to anaddb_dtset
 anaddb_dtset%natom=natom

!=====================================================================
!start reading in dimensions and non-dependent variables
!=====================================================================

!A

!typical value for gaussian smearing of a2F function
 anaddb_dtset%a2fsmear = 0.00002_dp
 token = 'a2fsmear'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'ENERGY')
 if(tread==1) anaddb_dtset%a2fsmear=dprarr(1)
 if (anaddb_dtset%a2fsmear < tol6) then
   write(message,'(a,a,a,f10.3,a,a,a,a,a)' )&
&   ' invars9 : Error -',ch10,&
&   '  a2fsmear is',anaddb_dtset%a2fsmear,', but only values > 1.e-6 ',&
&   ch10,'  are allowed',ch10,&
&   '  Action : correct a2fsmear in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%alphon=0
 token = 'alphon'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%alphon=intarr(1)
!FIXME: need a test on input value

 anaddb_dtset%asr=1
 token = 'asr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%asr=intarr(1)
 if(anaddb_dtset%asr<-2.or.anaddb_dtset%asr>5)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  asr is',anaddb_dtset%asr,', but the only allowed values',ch10,&
&   '  are 0, 1, 2, 3, 4, 5, -1 or -2 .',ch10,&
&   '  Action : correct asr in your input file.'
!  Note : negative values are allowed when the acoustic sum rule
!  is to be applied after the analysis of IFCs
!  3,4 are for rotational invariance (under development)
!  5 is for hermitian imposition of the ASR
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!B

 anaddb_dtset%brav=1
 token = 'brav'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%brav=intarr(1)
 if(anaddb_dtset%brav<=0.or.anaddb_dtset%brav>=5)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  brav is',anaddb_dtset%brav,', but the only allowed values',ch10,&
&   '  are 1,2,3 or 4 .',ch10,&
&   '  Action : correct brav in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!C

 anaddb_dtset%chneut=0
 token = 'chneut'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%chneut=intarr(1)
 if(anaddb_dtset%chneut<0.or.anaddb_dtset%chneut>2)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  chneut is',anaddb_dtset%chneut,', but the only allowed values',&
&   ch10,'  are 0, 1 or 2 .',ch10,&
&   '  Action : correct chneut in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!D

 anaddb_dtset%dieflag=0
 token = 'dieflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%dieflag=intarr(1)
 if(anaddb_dtset%dieflag<0.or.anaddb_dtset%dieflag>4)then
   write(message, '(3a,i8,5a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  dieflag is',anaddb_dtset%dieflag,', but the only allowed values',&
&   ch10,'  are 0, 1, 2, 3 or 4.',ch10,&
&   '  Action : correct dieflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%dipdip=1
 token = 'dipdip'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%dipdip=intarr(1)
 if(anaddb_dtset%dipdip<0.or.anaddb_dtset%dipdip>1)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  dipdip is',anaddb_dtset%dipdip,', but the only allowed values',&
&   ch10,'  are 0 or 1 .',ch10,&
&   '  Action : correct dipdip in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%ep_scalprod = 0
 token = 'ep_scalprod'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ep_scalprod = intarr(1)
 if(anaddb_dtset%ep_scalprod < 0 .or. anaddb_dtset%ep_scalprod > 1) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ep_scalprod is',anaddb_dtset%ep_scalprod,&
&   ', but the only allowed values',ch10,&
&   '  are 0 or 1 .',ch10,&
&   '  Action : correct ep_scalprod in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%dosdeltae=1.0/Ha_cmm1
 token = 'dosdeltae'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%dosdeltae=dprarr(1)
 if(anaddb_dtset%dosdeltae<=zero)then
   write(message, '(a,a,a,es14.4,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  dosdeltae is',anaddb_dtset%dosdeltae,', which is lower than 0 .',&
&   ch10,'  Action : correct dosdeltae in your input file.'
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if

!FIXME : should probably be smaller
 anaddb_dtset%dossmear=5.0/Ha_cmm1
 token = 'dossmear'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%dossmear=dprarr(1)
 if(anaddb_dtset%dossmear<=zero)then
   write(message, '(a,a,a,es14.4,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  dossmear is',anaddb_dtset%dossmear,', which is lower than 0 .',&
&   ch10,'  Action : correct dossmear in your input file.'
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if

 anaddb_dtset%dostol=0.25_dp
 token = 'dostol'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%dostol=dprarr(1)
 if(anaddb_dtset%dostol<zero)then
   write(message, '(a,a,a,es14.4,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  dostol is',anaddb_dtset%dostol,', which is lower than 0 .',ch10,&
&   '  Action : correct dostol in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!AHR

 anaddb_dtset%dossum=0
 token = 'dossum'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%dossum=intarr(1)
 if(anaddb_dtset%dossum < 0 .or. anaddb_dtset%dossum > one)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  dossum is',anaddb_dtset%dossum,', but the only allowed values',&
&   ch10,'  are 0, 1',ch10,&
&   '  Action : correct dossum in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!E

 anaddb_dtset%eivec=0
 token = 'eivec'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%eivec=intarr(1)
 if(anaddb_dtset%eivec<0.or.anaddb_dtset%eivec>4)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  eivec is',anaddb_dtset%eivec,', but the only allowed values',&
&   ch10,'  are 0, 1, 2, 3 or 4.',ch10,&
&   '  Action : correct eivec in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%elaflag=0
 token = 'elaflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%elaflag=intarr(1)
 if(anaddb_dtset%elaflag<0.or.anaddb_dtset%elaflag>5)then
   write(message,'(3a,i8,5a)' )&
&   ' invars9 : Error -',ch10,&
&   '  elaflag is',anaddb_dtset%elaflag,', but the only allowed values',&
&   ch10,'  are 0,1,2,3,4 or 5 .',ch10,&
&   '  Action : correct elaflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!By default use the real fermie (tests for elph_fermie == 0 in the code)
 anaddb_dtset%elph_fermie = zero
 token = 'elph_fermie'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'ENERGY')
 if(tread==1) anaddb_dtset%elph_fermie=dprarr(1)

!extra charge in unit cell (number of electrons) wrt neutral cell
!holes are negative values (reduce number of electrons)
 anaddb_dtset%ep_extrael = zero
 token = 'ep_extrael'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%ep_extrael=dprarr(1)

 anaddb_dtset%elphflag=0
 token = 'elphflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%elphflag=intarr(1)
 if(anaddb_dtset%elphflag<0.or.anaddb_dtset%elphflag>1)then
   write(message,'(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : Error -',ch10,&
&   '  elphflag =',anaddb_dtset%elphflag,', but the allowed values',&
&   ch10,'  are 0, or 1.',ch10,&
&   '  Action : correct elphflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!typical value for gaussian smearing, but can vary sensibly with the metal
 anaddb_dtset%elphsmear = 0.01_dp
 token = 'elphsmear'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'ENERGY')
 if(tread==1) anaddb_dtset%elphsmear=dprarr(1)
 if (anaddb_dtset%elphsmear < tol6) then
   write(message,'(a,a,a,f10.3,a,a,a,a,a)' )&
&   ' invars9 : Error -',ch10,&
&   '  elphsmear is',anaddb_dtset%elphsmear,'. Only values > 1.e-6 ',&
&   ch10,'  are allowed',ch10,&
&   '  Action : correct elphsmear in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%enunit=0
 token = 'enunit'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%enunit=intarr(1)
 if(anaddb_dtset%enunit<0.or.anaddb_dtset%enunit>2)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  enunit is',anaddb_dtset%enunit,', but the only allowed values',&
&   ch10,'  are 0, 1 or 2 . ',ch10,&
&   '  Action : correct enunit in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if


 anaddb_dtset%ep_alter_int_gam = 0
 token = 'ep_alter_int_gam'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ep_alter_int_gam = intarr(1)
 if(anaddb_dtset%ep_alter_int_gam /= 1 .and. anaddb_dtset%ep_alter_int_gam /= 0) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ep_alter_int_gam is',anaddb_dtset%ep_alter_int_gam,&
&   ', but the only allowed values',ch10,&
&   '  are 1 or 0.',ch10,&
&   '  Action : correct ep_alter_int_gam in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Default is 0 - not used unless telphint==2
 anaddb_dtset%ep_b_max = 0
 token = 'ep_b_max'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) then
   anaddb_dtset%ep_b_max = intarr(1)
   if(anaddb_dtset%ep_b_max < 1) then
     write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&     ' invars9 : ERROR -',ch10,&
&     '  ep_b_max is',anaddb_dtset%ep_b_max,&
&     ', but the only allowed values',ch10,&
&     '  are between 1 and nband.',ch10,&
&     '  Action : correct ep_b_max in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if

!Default is 0 - not used unless telphint==2
 anaddb_dtset%ep_b_min = 0
 token = 'ep_b_min'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) then
   anaddb_dtset%ep_b_min = intarr(1)
   if(anaddb_dtset%ep_b_min < 1) then
     write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&     ' invars9 : ERROR -',ch10,&
&     '  ep_b_min is',anaddb_dtset%ep_b_min,&
&     ', but the only allowed values',ch10,&
&     '  are between 1 and nband.',ch10,&
&     '  Action : correct ep_b_min in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if

 anaddb_dtset%ep_keepbands = 0
 token = 'ep_keepbands'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ep_keepbands = intarr(1)
 if(anaddb_dtset%ep_keepbands < 0 .or. anaddb_dtset%ep_keepbands > 1) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ep_keepbands is',anaddb_dtset%ep_keepbands,&
&   ', but the only allowed values',ch10,&
&   '  are 0 or 1 .',ch10,&
&   '  Action : correct ep_keepbands in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%ep_nqpt=0
 token = 'ep_nqpt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ep_nqpt = intarr(1)
 if(anaddb_dtset%ep_nqpt < 0) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ep_nqpt is',anaddb_dtset%ep_nqpt,&
&   ', but the only allowed values',ch10,&
&   '  are > 0.',ch10,&
&   '  Action : correct ep_nqpt in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if (anaddb_dtset%ep_nqpt > 0) then
   ABI_ALLOCATE(anaddb_dtset%ep_qptlist,(3,anaddb_dtset%ep_nqpt))
   if(3*anaddb_dtset%ep_nqpt>marr)then
     marr=3*anaddb_dtset%ep_nqpt
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%ep_qptlist(:,:)=zero
   token = 'ep_qptlist'
   call intagm(dprarr,intarr,jdtset,marr,3*anaddb_dtset%ep_nqpt,&
&   string(1:lenstr),token,tread,'DPR')
   if(tread==1) then
     anaddb_dtset%ep_qptlist(1:3,1:anaddb_dtset%ep_nqpt)=&
&     reshape(dprarr(1:3*anaddb_dtset%ep_nqpt),(/3,anaddb_dtset%ep_nqpt/))
   else
     write(message,'(a,a,a,a,a)')&
&     ' invars9 : Error -',ch10,&
&     '  ep_nqpt is non zero but ep_qptlist is absent ',ch10,&
&     '  Action : specify ep_qptlist in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if


!F
 anaddb_dtset%freeze_displ = zero
 token = 'freeze_displ'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%freeze_displ=dprarr(1)


 anaddb_dtset%frmax=ten
 token = 'frmax'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%frmax=dprarr(1)
 if (anaddb_dtset%frmax < 0) then
   write(message,'(a,a,a,f10.3,a,a,a,a,a)' )&
&   ' invars9 : Error -',ch10,&
&   '  frmax is',anaddb_dtset%frmax,'. Only values > 0 ',&
&   ch10,'  are allowed',ch10,&
&   '  Action : correct frmax in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%frmin=zero
 token = 'frmin'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%frmin=dprarr(1)
 if (anaddb_dtset%frmin < 0) then
   write(message,'(a,a,a,f10.3,a,a,a,a,a)' )&
&   ' invars9 : Error -',ch10,&
&   '  frmin is',anaddb_dtset%frmin,'. Only values > 0 ',&
&   ch10,'  are allowed',ch10,&
&   '  Action : correct frmin in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!G

 anaddb_dtset%gkk2write = 0
!token= 'gkk2write'
!call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
!& tread,'INT')
!if(tread==1) anaddb_dtset%gkk2write = intarr(1)
!if(anaddb_dtset%gkk2write < 0 .or. anaddb_dtset%gkk2write > 1) then
!write(message, '(a,a,a,i8,a,a,a,a,a)' )&
!&  ' invars9 : ERROR -',ch10,&
!&  '  gkk2write is',anaddb_dtset%gkk2write,&
!&  ', but the only allowed values',ch10,&
!&  '  are 0 or 1 .',ch10,&
!&  '  Action : correct gkk2write in your input file.'
!call wrtout(std_out,message,'COLL')
!call leave_new('COLL')
!end if

 anaddb_dtset%gkk_rptwrite = 0
!token= 'gkk_rptwrite'
!call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
!& tread,'INT')
!if(tread==1) anaddb_dtset%gkk_rptwrite = intarr(1)
!if(anaddb_dtset%gkk_rptwrite < 0 .or. anaddb_dtset%gkk_rptwrite > 1) then
!write(message, '(a,a,a,i8,a,a,a,a,a)' )&
!&  ' invars9 : ERROR -',ch10,&
!&  '  gkk_rptwrite is',anaddb_dtset%gkk_rptwrite,&
!&  ', but the only allowed values',ch10,&
!&  '  are 0 or 1 .',ch10,&
!&  '  Action : correct gkk_rptwrite in your input file.'
!call wrtout(std_out,message,'COLL')
!call leave_new('COLL')
!end if

 anaddb_dtset%gkqwrite = 0
 token = 'gkqwrite'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%gkqwrite = intarr(1)
 if(anaddb_dtset%gkqwrite < 0 .or. anaddb_dtset%gkqwrite > 1) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  gkqwrite is',anaddb_dtset%gkqwrite,&
&   ', but the only allowed values',ch10,&
&   '  are 0 or 1 .',ch10,&
&   '  Action : correct gkqwrite in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 
!H

!I

 anaddb_dtset%iavfrq=0
 token = 'iavfrq'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%iavfrq=intarr(1)
 if(anaddb_dtset%iavfrq<0.or.anaddb_dtset%iavfrq>1)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  iavfrq is',anaddb_dtset%iavfrq,', but the only allowed values',&
&   ch10,'  are 0 or 1 .',ch10,&
&   '  Action : correct iavfrq in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%ifcana=0
 token = 'ifcana'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ifcana=intarr(1)
 if(anaddb_dtset%ifcana<0.or.anaddb_dtset%ifcana>1)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ifcana is',anaddb_dtset%ifcana,', but the only allowed values',&
&   ch10,'  are 0 or 1 .',ch10,&
&   '  Action : correct ifcana in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%ifcflag=0
 token = 'ifcflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ifcflag=intarr(1)
 if(anaddb_dtset%ifcflag<0.or.anaddb_dtset%ifcflag>1)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ifcflag is',anaddb_dtset%ifcflag,', but the only allowed values',&
&   ch10,'  are 0 or 1 .',ch10,&
&   '  Action : correct ifcflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%prtsrlr=0
 token = 'prtsrlr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%prtsrlr=intarr(1)
 if(anaddb_dtset%prtsrlr<0.or.anaddb_dtset%prtsrlr>1)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  prtsrlr is',anaddb_dtset%prtsrlr,', but the only allowed values',&
&   ch10,'  are 0 or 1 .',ch10,&
&   '  Action : correct prtsrlr in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%ifcout=0
 token = 'ifcout'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ifcout=intarr(1)
 if(anaddb_dtset%ifcout<0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ifcout is',anaddb_dtset%ifcout,', which is lower than 0 .',ch10,&
&   '  Action : correct ifcout in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%ifltransport = 0
 token = 'ifltransport'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ifltransport = intarr(1)
 if(anaddb_dtset%ifltransport < 0 .or. anaddb_dtset%ifltransport > 1) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ifltransport is',anaddb_dtset%ifltransport,&
&   ', but the only allowed values',ch10,&
&   '  are 0 or 1 .',ch10,&
&   '  Action : correct ifltransport in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 
 anaddb_dtset%instrflag=0
 token = 'instrflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%instrflag=intarr(1)
 if(anaddb_dtset%instrflag<0.or.anaddb_dtset%instrflag>1)then
   write(message,'(3a,i8,5a)' )&
&   'invars9 : Error -',ch10,&
&   ' instrflag is',anaddb_dtset%instrflag,&
&   ', but the only allowed values',&
&   ch10,'  are 0, 1  .',ch10,&
&   '  Action : correct instrflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!J

!K

 anaddb_dtset%kptrlatt(:,:)=0
!why this test on reading in kptrlatt?
 marr = 9
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))
 token = 'kptrlatt'
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1)anaddb_dtset%kptrlatt(1:3,1:3)=reshape(intarr(1:9),(/3,3/))
!
!NOTE: no a priori way to test the validity of the integers in kptrlatt
!

 anaddb_dtset%kptrlatt_fine(:,:)=0
 marr = 9
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))
 token = 'kptrlatt_fine'
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1)anaddb_dtset%kptrlatt_fine(1:3,1:3)=reshape(intarr(1:9),(/3,3/))

 
!L

!M

!typical value for mustar, but can vary sensibly with the metal
 anaddb_dtset%mustar = 0.1_dp
 token = 'mustar'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%mustar=dprarr(1)
 if (anaddb_dtset%mustar < zero) then
   write(message,'(a,a,a,f10.3,a,a,a,a,a)' )&
&   ' invars9 : Error -',ch10,&
&   '  mustar is',anaddb_dtset%mustar,', but only positive values',ch10,&
&   '  are allowed',ch10,&
&   '  Action : correct mustar in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!N

 anaddb_dtset%natfix=0
 token = 'natfix'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%natfix=intarr(1)
 if(anaddb_dtset%natfix > natom)then
   write(message, '(a,a,a,i8,a,a,i4,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  natfix is',anaddb_dtset%natfix,', which is larger than natom',&
&   ' (=',natom,')',ch10,&
&   '  Action : correct natfix in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(anaddb_dtset%natfix < 0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  natfix is',anaddb_dtset%natfix,', which is < 0',ch10,&
&   '  Action : correct natfix in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%natifc=0
 token = 'natifc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%natifc=intarr(1)
 if(anaddb_dtset%natifc<0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  natifc is',anaddb_dtset%natifc,', which is lower than 0 .',ch10,&
&   '  Action : correct natifc in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%natprj_bs=0
 token = 'natprj_bs'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%natprj_bs=intarr(1)
 if(anaddb_dtset%natprj_bs<0 .or. anaddb_dtset%natprj_bs > natom)then
   write(message, '(a,a,a,i8,a,I8,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  natprj_bs is',anaddb_dtset%natprj_bs,', but must be between 0 and natom = ',natom,ch10,&
&   '  Action : correct natprj_bs in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%nchan=800
 token = 'nchan'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%nchan=intarr(1)
!FIXME: check this - it should probably be .ge. 1, not 0
 if(anaddb_dtset%nchan <0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  nchan is',anaddb_dtset%nchan,', which is lower than 0 .',ch10,&
&   '  Action : correct nchan in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%nfreq=1
 token = 'nfreq'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%nfreq=intarr(1)
 if(anaddb_dtset%nfreq<0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  nfreq is',anaddb_dtset%nfreq,', which is lower than 0 .',ch10,&
&   '  Action : correct nfreq in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%ng2qpt(:)=0
 token = 'ng2qpt'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ng2qpt(:)=intarr(1:3)
 do ii=1,3
   if(anaddb_dtset%ng2qpt(ii)<0)then
     write(message, '(a,a,a,i1,a,i8,a,a,a,i1,a)' )&
&     ' invars9 : ERROR -',ch10,&
&     '  ng2qpt(',ii,') is',anaddb_dtset%ng2qpt(ii),&
&     ', which is lower than 0 .',ch10,&
&     '  Action : correct ng2qpt(',ii,') in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end do

 anaddb_dtset%ngqpt(:)=0
 token = 'ngqpt'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ngqpt(1:3)=intarr(1:3)
 do ii=1,3
   if(anaddb_dtset%ngqpt(ii)<0)then
     write(message, '(a,a,a,i1,a,i8,a,a,a,i1,a)' )&
&     ' invars9 : ERROR -',ch10,&
&     '  ngqpt(',ii,') is',anaddb_dtset%ngqpt(ii),&
&     ', which is lower than 0 .',ch10,&
&     '  Action : correct ngqpt(',ii,') in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end do
!DEBUG
!write(std_out,*)' invars9 : ngqpt(:)=',anaddb_dtset%ngqpt(1:3)
!ENDDEBUG

 anaddb_dtset%ngrids=4
 token = 'ngrids'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ngrids=intarr(1)
 if(anaddb_dtset%ngrids<0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ngrids is',anaddb_dtset%ngrids,', which is lower than 0 .',ch10,&
&   '  Action : correct ngrids in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%nlflag=0
 token = 'nlflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%nlflag=intarr(1)
 if(anaddb_dtset%nlflag<0.or.anaddb_dtset%nlflag>2)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  nlflag is',anaddb_dtset%nlflag,', but the only allowed values',&
&   ch10,'  are 0, 1 or 2.',ch10,&
&   '  Action : correct nlflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%nph1l=0
 token = 'nph1l'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%nph1l=intarr(1)
 if(anaddb_dtset%nph1l<0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  nph1l is',anaddb_dtset%nph1l,', which is lower than 0 .',ch10,&
&   '  Action : correct nph1l in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 
 anaddb_dtset%nph2l=0
 token = 'nph2l'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%nph2l=intarr(1)
 if(anaddb_dtset%nph2l<0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  nph2l is',anaddb_dtset%nph2l,', which is lower than 0 .',ch10,&
&   '  Action : correct nph2l in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 
 anaddb_dtset%nqpath=0
 token = 'nqpath'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%nqpath=intarr(1)
 if(anaddb_dtset%nqpath<0)then
   write(message,'(a,a,a,i8,a,a,a)' )&
&   ' invars9 : Error -',ch10,&
&   '  nqpath is',anaddb_dtset%nqpath,', but must be positive',ch10,&
&   '  Action : correct elphflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%nqshft=1
 token = 'nqshft'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%nqshft=intarr(1)
 if(anaddb_dtset%nqshft<0 .or. anaddb_dtset%nqshft==3 .or.&
& anaddb_dtset%nqshft>=5 )then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  nqshft is',anaddb_dtset%nqshft,', but the only allowed values',&
&   ch10,'  are 1, 2 or 4 .',ch10,&
&   '  Action : correct nqshft in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%nsphere=0
 token = 'nsphere'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%nsphere=intarr(1)
 if(anaddb_dtset%nsphere<0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  nsphere is',anaddb_dtset%nsphere,', which is lower than 0',ch10,&
&   '  Action : correct nsphere in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%nstrfix=0
 token = 'nstrfix'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%nstrfix=intarr(1)
 if(anaddb_dtset%nstrfix > 6)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  nstrfix is',anaddb_dtset%nstrfix,', which is larger than 6',ch10,&
&   '  Action : correct nstrfix in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 if(anaddb_dtset%nstrfix < 0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  nstrfix is',anaddb_dtset%nstrfix,', which is < 0',ch10,&
&   '  Action : correct nstrfix in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%ntemper=500
 token = 'ntemper'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ntemper=intarr(1)
 if(anaddb_dtset%ntemper <0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ntemper is',anaddb_dtset%ntemper,', which is lower than 0',ch10,&
&   '  Action : correct ntemper in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%nwchan=10
 token = 'nwchan'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%nwchan=intarr(1)
!FIXME: check this - it should probably be .ge. 1, not 0
 if(anaddb_dtset%nwchan<0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  nwchan is',anaddb_dtset%nwchan,', which is lower than 0 .',ch10,&
&   '  Action : correct nwchan in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!O
 anaddb_dtset%outscphon = 0
 token = 'outscphon'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%outscphon=intarr(1)
 if(anaddb_dtset%outscphon<0.or.anaddb_dtset%outscphon>1)then
   write(message,'(3a,i8,5a)' )&
&   'invars9 : Error -',ch10,&
&   ' outscphon is',anaddb_dtset%outscphon,&
&   ', but the only allowed values',ch10,&
&   '  are 0 or 1  .',ch10,&
&   '  Action : correct outscphon in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if 


!P

 anaddb_dtset%piezoflag=0
 token = 'piezoflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%piezoflag=intarr(1)
 if(anaddb_dtset%piezoflag<0.or.anaddb_dtset%piezoflag>7)then
   write(message,'(3a,i8,5a)' )&
&   'invars9 : Error -',ch10,&
&   ' piezoflag is',anaddb_dtset%piezoflag,&
&   ', but the only allowed values',ch10,&
&   '  are 0, 1,2,3,4,5,6,7  .',ch10,&
&   '  Action : correct piezoflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%polflag=0
 token = 'polflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%polflag=intarr(1)
 if(anaddb_dtset%polflag<0.or.anaddb_dtset%polflag>1)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  polflag is',anaddb_dtset%polflag,&
&   ', but the only allowed values',ch10,&
&   '  are 0 or 1 .',ch10,&
&   '  Action : correct polflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Default is no output for PHDOS
 anaddb_dtset%prtdos=0
 token='prtdos'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%prtdos = intarr(1)
 if(anaddb_dtset%prtdos < 0 .or. anaddb_dtset%prtdos > 2) then
   write(message, '(3a,i8,5a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  prtdos is',anaddb_dtset%prtdos,', but the only allowed values',&
&   ch10,'  are 0 (no output) or 1 (gaussians) or 2 (tetrahedra) ',&
&   ch10,'  Action : correct prtdos in your input file.'
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if
 
!Default is no output for the Fermi Surface
 anaddb_dtset%prtfsurf = 0
 token = 'prtfsurf'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%prtfsurf = intarr(1)
 if(anaddb_dtset%prtfsurf < 0 .or. anaddb_dtset%prtfsurf > 2) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  prtfsurf is',anaddb_dtset%prtfsurf,'. The only allowed values',&
&   ch10,'  are 0 (no output) or 1 (Xcrysden bxsf format)',ch10,  &
&   '  Action : correct prtfsurf in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%prtmbm=0
 token = 'prtmbm'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%prtmbm=intarr(1)
!FIXME: should check whether value of prtmbm is valid

!Default is no output of the nesting factor
 anaddb_dtset%prtnest = 0
 token = 'prtnest'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%prtnest = intarr(1)
 if(anaddb_dtset%prtnest < 0 .or. anaddb_dtset%prtnest > 2) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  prtnest is',anaddb_dtset%prtnest,' The only allowed values',ch10,&
&   '  are 0 (no nesting), 1 (XY format) or 2 (XY + Xcrysden format)',&
&   ch10,'  Action : correct prtnest in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Q

 anaddb_dtset%q2shft(:)=zero
 token = 'q2shft'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%q2shft(:)=dprarr(1:3)
!FIXME: need a test on valid entries for q2shft

 anaddb_dtset%qgrid_type=1 ! default is uniform nqpt(:) grid
 token = 'qgrid_type'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%qgrid_type = intarr(1)
 if(anaddb_dtset%qgrid_type < 1 .or. anaddb_dtset%qgrid_type > 2) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  qgrid_type is',anaddb_dtset%qgrid_type,' The only allowed values',ch10,&
&   '  are 1 (uniform grid from nqpt) or 2 (listed in ep_nqpt, ep_qptlist)',&
&   ch10,'  Action : correct qgrid_type in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%qrefine=1 ! default is no refinement
 token = 'qrefine'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%qrefine = intarr(1)
 if(anaddb_dtset%qrefine < 1) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  qrefine is',anaddb_dtset%qrefine,' The only allowed values',ch10,&
&   '  are integers >= 1 giving the refinement of the ngqpt grid',&
&   ch10,'  Action : correct qrefine in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if



!R

 anaddb_dtset%ramansr=0
 token = 'ramansr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ramansr=intarr(1)
 if(anaddb_dtset%ramansr<0.or.anaddb_dtset%ramansr>2)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ramansr is',anaddb_dtset%ramansr,', but the only allowed values',&
&   ch10,'  are 0, 1 or 2 .',ch10,&
&   '  Action : correct ramansr in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%relaxat=0
 token = 'relaxat'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%relaxat=intarr(1)
 if(anaddb_dtset%relaxat < 0.or.anaddb_dtset%relaxat > 1)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  relaxat is',anaddb_dtset%relaxat,', but the only allowed values',&
&   ch10,'  are 0 or 1 .',ch10,&
&   '  Action : correct relaxat in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%relaxstr=0
 token = 'relaxstr'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%relaxstr=intarr(1)
 if(anaddb_dtset%relaxstr<0.or.anaddb_dtset%relaxstr>1)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  relaxstr is',anaddb_dtset%relaxstr,&
&   '  but the only allowed values',ch10,&
&   '  are 0 or 1 .',ch10,&
&   '  Action : correct relaxstr in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%rfmeth=1
 token = 'rfmeth'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%rfmeth=intarr(1)
 if(anaddb_dtset%rfmeth<1.or.anaddb_dtset%rfmeth>2)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  rfmeth is',anaddb_dtset%rfmeth,', but the only allowed values',&
&   ch10,'  are 1 or 2 . ',ch10,&
&   '  Action : correct rfmeth in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%rifcsph=zero
 token = 'rifcsph'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%rifcsph=dprarr(1)
 if(anaddb_dtset%rifcsph<-tol12)then
   write(message, '(a,a,a,f10.3,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  rifcsph is',anaddb_dtset%rifcsph,', which is lower than zero.',&
&   ch10,'  Action : correct rifcsph in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!S

 anaddb_dtset%selectz=0
 token = 'selectz'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%selectz=intarr(1)
 if(anaddb_dtset%selectz<0.or.anaddb_dtset%selectz>2)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  selectz is',anaddb_dtset%selectz,', but the only allowed values',&
&   ch10,'  are 0, 1 or 2 .',ch10,&
&   '  Action : correct selectz in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%symdynmat=1
 token = 'symdynmat'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%symdynmat=intarr(1)
 if(anaddb_dtset%symdynmat/=0.and.anaddb_dtset%symdynmat/=1)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  symdynmat is',anaddb_dtset%symdynmat,'. The only allowed values',&
&   ch10,'  are 0, or 1.',ch10,&
&   '  Action : correct symdynmat in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!T

 anaddb_dtset%targetpol(:) = 0._dp
 token = 'targetpol'
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%targetpol(1:3) = dprarr(1:3)

!Default is use gaussian integration
 anaddb_dtset%telphint = 1
 token = 'telphint'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%telphint = intarr(1)
 if(anaddb_dtset%telphint < 0 .or. anaddb_dtset%telphint > 3) then
   write(message, '(a,a,a,i8,a,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  telphint is',anaddb_dtset%telphint,'. The only allowed values',&
&   ch10,'  are 0 (tetrahedron) or 1 (gaussian) or ',&
&   '2 (set of bands occupied ep_b_min,ep_b_max) or 3 (Fermi Dirac).',ch10,&
&   '  Action : correct telphint in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%temperinc=two
 token = 'temperinc'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%temperinc=dprarr(1)
 if(anaddb_dtset%temperinc < zero)then
   write(message, '(a,a,a,f10.3,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  temperinc is',anaddb_dtset%temperinc,', which is lower than 0 .',&
&   ch10,'  Action : correct temperinc in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%tempermin=one
 token = 'tempermin'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%tempermin=dprarr(1)
 if(anaddb_dtset%tempermin<-tol12)then
   write(message, '(a,a,a,f10.3,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  tempermin is',anaddb_dtset%tempermin,', which is lower than 0 .',&
&   ch10,'  Action : correct tempermin in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%thmflag=0
 token = 'thmflag'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%thmflag=intarr(1)
 if(anaddb_dtset%thmflag<0.or.anaddb_dtset%thmflag>8)then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  thmflag is',anaddb_dtset%thmflag,', but the only allowed values',&
&   ch10,'  are between 0 to 8 (included).',ch10,&
&   '  Action : correct thmflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%thmtol=0.25_dp
 token = 'thmtol'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'DPR')
 if(tread==1) anaddb_dtset%thmtol=dprarr(1)
 if(anaddb_dtset%thmtol<zero)then
   write(message, '(a,a,a,es14.4,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  thmtol is',anaddb_dtset%thmtol,', which is lower than 0 .',ch10,&
&   '  Action : correct thmtol in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 anaddb_dtset%ep_prt_yambo = 0
 token = 'ep_prt_yambo'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%ep_prt_yambo = intarr(1)
 if(anaddb_dtset%ep_prt_yambo< 0 .or. anaddb_dtset%ep_prt_yambo> 1) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ep_prt_yambo is',anaddb_dtset%ep_prt_yambo,&
&   ', but the only allowed values',ch10,&
&   '  are 0 or 1 .',ch10,&
&   '  Action : correct ep_prt_yambo in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!default means _do_ symmetrize the ep coupling matrices over qpoints
 anaddb_dtset%symgkq = 1
 token = 'symgkq'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,&
& tread,'INT')
 if(tread==1) anaddb_dtset%symgkq = intarr(1)
 if(anaddb_dtset%symgkq< 0 .or. anaddb_dtset%symgkq> 1) then
   write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  symgkq is',anaddb_dtset%symgkq,&
&   ', but the only allowed values',ch10,&
&   '  are 0 or 1.',ch10,&
&   '  Action : correct symgkq in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 else if (anaddb_dtset%symgkq == 0) then
   write (message,'(a)') &
&   ' WARNING: you have turned off el-ph matrix symmetrization over q. Use at own risk'
 end if

!U

!V

!W

!X

!Y

!Z

!=====================================================================
!end non-dependent variables
!=====================================================================

!=======================================================================
!Read in dependent variables (dependent on dimensions above)
!=======================================================================

!A

 ABI_ALLOCATE(anaddb_dtset%atifc,(natom))
 anaddb_dtset%atifc(:)=0
 if(anaddb_dtset%natifc>=1)then
   if(anaddb_dtset%natifc>marr)then
     marr=anaddb_dtset%natifc
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   token = 'atifc'
   call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%natifc,&
&   string(1:lenstr),token,tread,'INT')
   if(tread==1) anaddb_dtset%atifc(1:anaddb_dtset%natifc)=&
&   intarr(1:anaddb_dtset%natifc)
!  check of whether values of atifc are valid is done in chkin9
 end if

!B

!C

!D

!E

!F

!G

!H

!I

 ABI_ALLOCATE(anaddb_dtset%iatfix,(natom))
 anaddb_dtset%iatfix(:) = 0
 if ((anaddb_dtset%relaxat == 1).and.(anaddb_dtset%natfix > 0)) then
   if(natom > marr)then
     marr = natom
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   token = 'iatfix'
   call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%natfix,&
&   string(1:lenstr),token,tread,'INT')
   if(tread==1) anaddb_dtset%iatfix(1:anaddb_dtset%natfix) = &
&   intarr(1:anaddb_dtset%natfix)
 end if
!FIXME: need a test on values of iatfix: are they just 1 or 0?

 if ((anaddb_dtset%relaxstr == 1).and.(anaddb_dtset%nstrfix > 0)) then
   anaddb_dtset%istrfix(:) = 0
   token = 'istrfix'
   call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%nstrfix,&
&   string(1:lenstr),token,tread,'INT')
   if(tread==1) anaddb_dtset%istrfix(1:anaddb_dtset%nstrfix) = &
&   intarr(1:anaddb_dtset%nstrfix)
 end if
!FIXME: need a test on values of istrfix

 nullify (anaddb_dtset%iatprj_bs)
 if (anaddb_dtset%natprj_bs > 0) then
   ABI_ALLOCATE(anaddb_dtset%iatprj_bs,(anaddb_dtset%natprj_bs))
   if(anaddb_dtset%natprj_bs>marr)then
     marr=anaddb_dtset%natprj_bs
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%iatprj_bs(:)=zero
   token = 'iatprj_bs'
   call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%natprj_bs,&
&   string(1:lenstr),token,tread,'INT')
   if(tread==1) then
     anaddb_dtset%iatprj_bs(1:anaddb_dtset%natprj_bs)=&
&     intarr(1:anaddb_dtset%natprj_bs)
   else
     write(message,'(a,a,a,a,a)')&
&     ' invars9 : Error -',ch10,&
&     '  natprj_bs is non zero but iatprj_bs is absent ',ch10,&
&     '  Action : specify iatprj_bs in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if

!J

!K

!L

!M

!N

!O

!P

!Q

 if (anaddb_dtset%nqshft/=0)then
   if(3*anaddb_dtset%nqshft>marr)then
     marr=3*anaddb_dtset%nqshft
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%q1shft(:,:)=zero
   token = 'q1shft'
   call intagm(dprarr,intarr,jdtset,marr,3*anaddb_dtset%nqshft,&
&   string(1:lenstr),token,tread,'DPR')
   if(tread==1) anaddb_dtset%q1shft(1:3,1:anaddb_dtset%nqshft)=&
&   reshape(dprarr(1:3*anaddb_dtset%nqshft),(/3,anaddb_dtset%nqshft/))
 end if

!DEBUG
!write(std_out,*)' invars9 : before allocate qph1l, qnrml1, nph1l=',anaddb_dtset%nph1l
!stop
!ENDDEBUG
 ABI_ALLOCATE(anaddb_dtset%qph1l,(3,anaddb_dtset%nph1l))
 ABI_ALLOCATE(anaddb_dtset%qnrml1,(anaddb_dtset%nph1l))
 if (anaddb_dtset%nph1l/=0)then
   if(4*anaddb_dtset%nph1l>marr)then
     marr=4*anaddb_dtset%nph1l
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%qph1l(:,:)=zero
   anaddb_dtset%qnrml1(:)=zero
   token = 'qph1l'
   call intagm(dprarr,intarr,jdtset,marr,4*anaddb_dtset%nph1l,&
&   string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     do iph1=1,anaddb_dtset%nph1l
       do ii=1,3
         anaddb_dtset%qph1l(ii,iph1)=dprarr(ii+(iph1-1)*4)
       end do
       anaddb_dtset%qnrml1(iph1)=dprarr(4+(iph1-1)*4)
       if(abs(anaddb_dtset%qnrml1(iph1))<qtol)then
         write(message, '(a,a,a,a,a,a,a)' )&
&         ' invars9 : ERROR -',ch10,&
&         '  The first list of wavevectors ',&
&         'should not have non-analytical data.',ch10,&
&         '  Action : correct the first list',&
&         ' of wavevectors in the input file.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end if
 end if

!DEBUG
!write(std_out,*)' invars9 : before allocate qph2l, qnrml2, nph2l=',anaddb_dtset%nph2l
!ENDDEBUG
 ABI_ALLOCATE(anaddb_dtset%qph2l,(3,anaddb_dtset%nph2l))
 ABI_ALLOCATE(anaddb_dtset%qnrml2,(anaddb_dtset%nph2l))
 if (anaddb_dtset%nph2l/=0)then
   if(4*anaddb_dtset%nph2l>marr)then
     marr=4*anaddb_dtset%nph2l
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%qph2l(:,:)=zero
   anaddb_dtset%qnrml2(:)=zero
   token = 'qph2l'
   call intagm(dprarr,intarr,jdtset,marr,4*anaddb_dtset%nph2l,&
&   string(1:lenstr),token,tread,'DPR')
   if(tread==1)then
     do iph2=1,anaddb_dtset%nph2l
       do ii=1,3
         anaddb_dtset%qph2l(ii,iph2)=dprarr(ii+(iph2-1)*4)
       end do
       anaddb_dtset%qnrml2(iph2)=dprarr(4+(iph2-1)*4)
       if(abs(anaddb_dtset%qnrml2(iph2))>qtol)then
         write(message, '(a,a,a,a,a,a,a)' )&
&         ' invars9 : ERROR -',ch10,&
&         '  The second list of wavevectors',&
&         ' should have only non-analytical data.',ch10,&
&         '  Action : correct the second list',&
&         ' of wavevectors in the input file.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end do
   end if
 end if

 if (anaddb_dtset%nqpath > 0) then
   ABI_ALLOCATE(anaddb_dtset%qpath,(3,anaddb_dtset%nqpath))
   if(3*anaddb_dtset%nqpath>marr)then
     marr=3*anaddb_dtset%nqpath
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%qpath(:,:)=zero
   token = 'qpath'
   call intagm(dprarr,intarr,jdtset,marr,3*anaddb_dtset%nqpath,&
&   string(1:lenstr),token,tread,'DPR')
   if(tread==1) then
     anaddb_dtset%qpath(1:3,1:anaddb_dtset%nqpath)=&
&     reshape(dprarr(1:3*anaddb_dtset%nqpath),(/3,anaddb_dtset%nqpath/))
   else
     write(message,'(a,a,a,a,a)')&
&     ' invars9 : Error -',ch10,&
&     '  nqpath is non zero but qpath is absent ',ch10,&
&     '  Action : specify qpath in your input file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if

!R

!S

!T

!U

!V

!W

!X

!Y

!Z

!=======================================================================
!Finished reading in variables - deallocate
!=======================================================================

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)

!=======================================================================
!Check consistency of input variables:
!=======================================================================

 if (anaddb_dtset%frmin > anaddb_dtset%frmax) then
   write(message,'(a,a,a,a,a)' )&
&   ' invars9 : Error -',ch10,&
&   '  frmax should be higher than frmin',ch10,&
&   '  Action : change frmax and/or frmin  in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 
 if (anaddb_dtset%nqpath==0 .and. anaddb_dtset%elphflag==1) then
   write(message,'(a,a,a,a,a,a)' )&
&   ' invars9 : Error -',ch10,&
&   '  elphflag is 1 but no nqpath has been specified',&
&   ' for phonon linewidths',ch10,&
&   '  Action : specify nqpath and qpath(3,nqpath) in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(anaddb_dtset%telphint /= 2 .and. (anaddb_dtset%ep_b_min /= 0 .or.&
& anaddb_dtset%ep_b_max /= 0)) then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' invars9 : WARNING -',ch10,&
&   '  telphint is',anaddb_dtset%telphint,', but ep_b_min or ep_b_max',&
&   ch10,'  are set /= 1. They will not be used'
   call wrtout(std_out,message,'COLL')
 else if(anaddb_dtset%telphint == 2 .and. &
&   (anaddb_dtset%ep_b_min == 0 .or. anaddb_dtset%ep_b_max == 0)) then
   write(message, '(a,a,a,i8,a,a,a,a)' )&
&   ' invars9 : WARNING -',ch10,&
&   '  telphint is',anaddb_dtset%telphint,', but ep_b_min or ep_b_max',&
&   ch10,'  are not both set. ',ch10,&
&   '  Action : set ep_b_min and ep_b_max in your input file.',ch10
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(anaddb_dtset%thmflag < 3) then
   if ((anaddb_dtset%telphint == 0 .or. anaddb_dtset%prtnest == 1 .or. &
&   anaddb_dtset%prtnest == 2 .or. anaddb_dtset%prtfsurf== 1) .and.&
&   sum(anaddb_dtset%kptrlatt) == 0 ) then
     write (message, '(a,a,a)') &
&     ' invars9 :  ERROR : if tetrahedron integration is used, ',&
&     'or the output of the nesting function/Fermi surface is required, ',&
&     'you must specify the kptrlatt'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if

 if(anaddb_dtset%prtdos/=0 .and. anaddb_dtset%ifcflag/=1) then
   write(message, '(5a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ifcflag must be 1 when the calculation of the phonon DOS is required ',ch10,&
&   '  Action : correct ifcflag in your input file.'
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if

 if(anaddb_dtset%prtsrlr/=0 .and. anaddb_dtset%ifcflag/=1) then
   write(message, '(5a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ifcflag must be 1 for the SR/LR decomposition of the phonon frequencies',ch10,&
&   '  Action : correct ifcflag in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(anaddb_dtset%prtdos/=0 .and. sum(abs(anaddb_dtset%ng2qpt(:))) < 3 ) then
   write(message, '(5a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  ng2qpt must be specified when the calculation of the phonon DOS is required ',ch10,&
&   '  Action : correct ng2qpt in your input file.'
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if

 if (anaddb_dtset%ifltransport == 1 .and. anaddb_dtset%ep_keepbands /= 1) then
   write(message, '(5a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  Band dependency of electron phonon matrix elements must be kept for transport ',ch10,&
&   '  Action : set ep_keepbands to 1 in your input file.'
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if

 if (anaddb_dtset%ep_alter_int_gam == 1 .and. anaddb_dtset%ep_keepbands /= 1) then
   write(message, '(5a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  Band dependency of electron phonon matrix elements must be kept for alternative integration mode ',ch10,&
&   '  Action : set ep_keepbands to 1 in your input file.'
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if

 if (anaddb_dtset%ep_alter_int_gam == 1 .and. &
& (anaddb_dtset%kptrlatt(1,1)+anaddb_dtset%kptrlatt(2,2)+anaddb_dtset%kptrlatt(3,3))==0) then
   write(message, '(5a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  kptrlatt must be specified for alternative integration mode ',ch10,&
&   '  Action : set kptrlatt (from ground state output) in your input file.'
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if

 if (anaddb_dtset%ep_alter_int_gam == 1 .and. &
& (anaddb_dtset%kptrlatt_fine(1,1)+anaddb_dtset%kptrlatt_fine(2,2)+anaddb_dtset%kptrlatt_fine(3,3))==0) then
   write(message, '(5a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  kptrlatt_fine must be specified for alternative integration mode ',ch10,&
&   '  Action : set kptrlatt_fine (from dense k-grid run) in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!FIXME: add check that if freeze_displ /= 0 then you need to be doing ifc and phonon interpolation

 if (anaddb_dtset%ifcflag > 0 .and. sum(abs(anaddb_dtset%ngqpt)) == 0) then
   write(message, '(5a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  if you want interatomic force constant output, anaddb needs ngqpt input variable ',ch10,&
&   '  Action : set ngqpt in your input file.'
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
 end if

!check that q-grid refinement is a divisor of ngqpt in each direction
 if(anaddb_dtset%qrefine > 1 .and. sum(abs(dmod(anaddb_dtset%ngqpt/dble(anaddb_dtset%qrefine),one))) > tol10) then
   write(message, '(a,a,a,i8,a,a,a,3i8,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  qrefine is',anaddb_dtset%qrefine,' The only allowed values',ch10,&
&   '  are integers which are divisors of the ngqpt grid', anaddb_dtset%ngqpt,&
&   ch10,'  Action : correct qrefine in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!check that fermie and nelect are not both specified
 if(abs(anaddb_dtset%elph_fermie) > tol10 .and. &
& abs(anaddb_dtset%ep_extrael) > tol10) then
   write(message, '(a,a,a,E10.2,a,E10.2,a,a,a)' )&
&   ' invars9 : ERROR -',ch10,&
&   '  elph_fermie (',anaddb_dtset%elph_fermie,') and ep_extrael (',&
&   anaddb_dtset%ep_extrael, '), may not both be non 0  ', &
&   ch10,'  Action : remove one of the two in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 


!DEBUG
!write(std_out,*)' invars9 : exit '
!ENDDEBUG

end subroutine invars9
!!***
