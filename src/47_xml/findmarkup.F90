!{\src2tex{textfont=tt}}
!!****f* ABINIT/findmarkup
!! NAME
!! findmarkup
!!
!! FUNCTION
!! Given a XML mark-up, identifies in a zone of a XML string the places where
!! the markup first occurs  (i.e. find '<markup' ),
!! where the first occurence ends (i.e. the next '>' or '/>' ),
!! and, if the previous was not '/>',  where it is closed   (i.e. find '</markup' ) .
!! Do not treat 'empty' markups yet.
!! Also able to identify a markup that refer to a specific 'title' (often
!! used in CML.
!! Returns a triplet of zero if not found.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  builtin=(string) if non-blank, the routine will select a markup that contains
!!   an attribute title="value_in_title"
!!  index_lower=lower index for search of the markup in the string
!!  index_upper=upper index for search of the markup in the string
!!  markup=the markup
!!  markuplen=length of the markup
!!  strln=maximal number of characters of string
!!  string_xml*(strln)=string of characters to be searched
!!
!! OUTPUT
!!  indices_markup(3)=the three indices : position of the '<' in '<markup',
!!   position of the '>' in '<markup ...>', or '/>' in '<markup .../>', and
!!   position of the '<' in '</markup', if the previous was not '/>'.
!!   Contains a triplet of 0 if the adequate markup was not found
!!   The last index is 0 if '<markup .../>'
!!
!! NOTES
!!  Should translate "builtin" in "title"
!!
!! PARENTS
!!      append_cml2
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine findmarkup(builtin,index_lower,index_upper,indices_markup,&
& markup,markuplen,strln,string_xml)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'findmarkup'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: index_lower,index_upper,markuplen,strln
 character(len=*),intent(in) :: builtin
 character(len=markuplen),intent(in) :: markup
 character(len=strln),intent(in) :: string_xml
!arrays
 integer,intent(out) :: indices_markup(3)

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: found,index_builtin,index_equal,index_lower_trial,index_markup

!************************************************************************

!DEBUG
!write(std_out,*)' findmarkup : enter'
!write(std_out,*)' trim(markup)="',trim(markup),'"'
!write(std_out,*)' trim(builtin)="',trim(builtin),'"'
!write(std_out,*)' string_xml(index_lower:index_upper)='
!write(std_out,'(3a)' )'"',string_xml(index_lower:index_upper),'"'
!stop
!ENDDEBUG

 index_lower_trial=index_lower

!Search for the proper occurence of the markup
 do

!  DEBUG
!  write(std_out,*)' new trial, with index_lower_trial=',index_lower_trial
!  ENDDEBUG

!  Find the critical indices of the mark-up in CML string
!  The character following the markup must be a blank or a '<'
   indices_markup(1)=index(string_xml(index_lower_trial:index_upper),"<"//markup//blank)
   if(indices_markup(1) > 0)then
     index_markup=index(string_xml(index_lower_trial:indices_markup(1)),"<"//markup//">")
     if(index_markup > 0)then
       if(index_markup<indices_markup(1))indices_markup(1)=index_markup
     end if
   else
     indices_markup(1)=index(string_xml(index_lower_trial:index_upper),"<"//markup//">")
   end if

   indices_markup(1)=indices_markup(1)+index_lower_trial-1

   if(indices_markup(1) < index_lower_trial )then
     indices_markup(:)=0
     exit
   end if

   index_markup=index(string_xml(indices_markup(1):index_upper),"/>")&
&   +indices_markup(1)-1

   indices_markup(2)=index(string_xml(indices_markup(1):index_upper),">")&
&   +indices_markup(1)-1

!  Short form or long form
   if(index_markup>indices_markup(1) .and. index_markup<indices_markup(2))then

     indices_markup(2)=index_markup
     indices_markup(3)=0

   else if(indices_markup(2) < indices_markup(1))then

     indices_markup(:)=0
     exit

   else

     indices_markup(3)=index(string_xml(indices_markup(2):index_upper),"</"//markup//">")&
&     +indices_markup(2)-1

     if(indices_markup(3) < indices_markup(2))then
       indices_markup(:)=0
       exit
     end if

   end if ! short or long form

!  DEBUG
!  write(std_out,*)' indices_markup(:)=',indices_markup(:)
!  ENDDEBUG


!  If needed, detect a title
   if(len_trim(builtin)/=0)then

     index_builtin=index(string_xml(indices_markup(1):indices_markup(2)),"title") &
&     +indices_markup(1)-1
     if(index_builtin<indices_markup(1))then
       index_lower_trial=maxval(indices_markup(2:3))
       cycle
     end if

!    DEBUG
!    write(std_out,*)' index_builtin=',index_builtin
!    ENDDEBUG

     index_equal=index(string_xml(index_builtin+5:index_builtin+6),"=") &
&     +index_builtin+5-1
     if(index_equal<index_builtin+5-1)then
       index_lower_trial=maxval(indices_markup(2:3))
       cycle
     end if

!    DEBUG
!    write(std_out,*)' index_equal=',index_equal
!    ENDDEBUG

     found=index(string_xml(index_equal+1:index_equal+3+len_trim(builtin)),&
&     '"'//trim(builtin)//'"')

!    DEBUG
!    write(std_out,*)' found=',found
!    ENDDEBUG

     if(found/=0)exit      ! Succeeded to find the builtin, exit the do-loop search
     index_lower_trial=maxval(indices_markup(2:3))
     cycle

   else ! If no builtin is needed, then the candidate is OK

     exit

   end if

 end do ! End of do-loop on suitable candidates

!DEBUG
!write(std_out,*)' findmarkup : exit ',ch10
!stop
!ENDDEBUG

end subroutine findmarkup
!!***
