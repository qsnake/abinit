!{\src2tex{textfont=tt}}
!!****f* ABINIT/importcml
!! NAME
!! importcml
!!
!! FUNCTION
!! Examine the input string, to see whether data from CML
!! file(s) has to be incorporated.
!! For each such CML file, translate the relevant
!! information into intermediate input variables compatible
!! with the usual ABINIT formatting, then append it
!! to the input string.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  string_raw*(strln)=raw string of character from input file (with original case)
!!  strln=maximal number of character of string, as declared in the calling routine
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  lenstr=actual number of character in string
!!  string_upper*(strln)=string of character
!!   the string (with upper case) from the input file, to which the CML data are appended to it
!!
!! PARENTS
!!      m_ab6_invars_f90,newsp,parsefile
!!
!! CHILDREN
!!      append_cml2,incomprs,instrng,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine importcml (lenstr,string_raw,string_upper,strln)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'importcml'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_42_parser
 use interfaces_47_xml, except_this_one => importcml
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: strln
 integer,intent(inout) :: lenstr
 character(len=*),intent(in) :: string_raw
 character(len=*),intent(inout) :: string_upper

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: dtset_len,icml,ii,index_already_done,index_cml_fname
 integer :: index_cml_fname_end,index_cml_token,kk,lenstr_cml,option
 character(len=2) :: dtset_char
 character(len=500) :: message
 character(len=fnlen) :: cml_fname
 character(len=strln) :: string_cml

!************************************************************************

!DEBUG
!write(std_out,*)' importcml : enter '
!stop
!ENDDEBUG

 index_already_done=1
 icml=0

 do    ! Infinite do-loop, to identify the presence of the CMLFILE token

   index_cml_token=index(string_upper(index_already_done:lenstr),"CMLFILE")
   if(index_cml_token==0)exit

   icml=icml+1
   if(icml==1)then
     write(message,'(80a)')('=',ii=1,80)
     call wrtout(ab_out,message,'COLL')
   end if

!  The CMLFILE token has been identified
   index_cml_token=index_already_done+index_cml_token-1

!  write(std_out,*)' index_cml_token =>',string_upper(index_cml_token:index_cml_token+4)

!  Find the related dataset tag, and length
   dtset_char=string_upper(index_cml_token+7:index_cml_token+8)

!  write(std_out,*)' dtset_char =>',string_upper(index_cml_token+7:index_cml_token+8)

   if(dtset_char(1:1)==blank)dtset_char(2:2)=blank
   dtset_len=len_trim(dtset_char)

!  write(std_out,*)' dtset_len =>',dtset_len

!  Find the name of the CML file
   index_cml_fname=index_cml_token+8+dtset_len

!  write(std_out,*)' index_cml_fname =>',string_upper(index_cml_fname:index_cml_fname+10)

   index_cml_fname_end=index(string_upper(index_cml_fname:lenstr),blank)

   if(index_cml_fname_end ==0 )then
     write(message, '(8a,i4,2a)' ) ch10,&
&     ' importcml : ERROR - ',ch10,&
&     '  Could not find the name of the CML file.',ch10,&
&     '  index_cml_fname_end should be non-zero, while it is :',ch10,&
&     '  index_cml_fname_end=',index_cml_fname_end,ch10,&
&     '  Action : check the filename that was provided after the CMLFILE input variable keyword.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   index_cml_fname_end=index_cml_fname_end+index_cml_fname-1

   index_already_done=index_cml_fname_end

   cml_fname=repeat(blank,fnlen)                  ! Initialize cml_fname to a blank line
   cml_fname=string_raw(index_cml_fname:index_cml_fname_end-1)

   write(message, '(3a)') ch10,&
&   ' importcml : Identified token CMLFILE, referring to file ',trim(cml_fname)
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  For the time being, use the common ABINIT parsing routine,
!  but this might lead to bizarre effects. Should modify this ...
   option=0
   call instrng (cml_fname,lenstr_cml,option,strln,string_cml)

   write(message, '(3a)') &
&   ' importcml : Opened file ',trim(cml_fname),'; content stored in string_cml'
   call wrtout(std_out,message,'COLL')

!  write(std_out,*)trim(string_cml)

!  Append the data from the CML file to the string, and update the length of the string
   call append_cml2(dtset_char,lenstr,lenstr_cml,string_upper,string_cml,strln)

 end do ! infinite loop while we find cml file tokens

!append the final ntypat/znucl after all dtset have been accumulated
 if (index_already_done > 1) then
   call append_cml2("-1",lenstr,0,string_upper,string_cml,strln)
 end if

 if(icml/=0)then
   call incomprs(string_upper,lenstr)
!  A blank is needed at the beginning of the string
   do kk=lenstr,1,-1
     string_upper(kk+1:kk+1)=string_upper(kk:kk)
   end do
   string_upper(1:1)=blank
   lenstr=lenstr+1
   write(message,'(a,80a,a)')ch10,('=',ii=1,80),ch10
   call wrtout(ab_out,message,'COLL')
 end if

!DEBUG
!write(std_out,*)' importcml : exit '
!stop
!ENDDEBUG

end subroutine importcml
!!***
