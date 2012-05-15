!{\src2tex{textfont=tt}}
!!****f* ABINIT/append_cml2
!! NAME
!! append_cml2
!!
!! FUNCTION
!! Translate the data from a CML2 string (string_cml), of length lenstr_cml,
!! and add it at the end of the usual ABINIT input data string (string),
!! taking into account the dtset (dtset_char)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset_char*2=possible dtset label
!!  lenstr_cml=actual number of characters in string
!!  string_cml*(strln)=string of characters from the CML2 file
!!  strln=maximal number of characters of string, as declared in the calling routine
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  lenstr=actual number of characters in string
!!  string*(strln)=string of characters  (upper case) to which the CML data are appended
!!
!! PARENTS
!!      importcml
!!
!! CHILDREN
!!      findmarkup,getattribute,leave_new,symbol2znucl,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine append_cml2(dtset_char,lenstr,lenstr_cml,string,string_cml,strln)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'append_cml2'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_47_xml, except_this_one => append_cml2
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lenstr_cml,strln
 integer,intent(inout) :: lenstr
 character(len=2),intent(in) :: dtset_char
 character(len=strln),intent(in) :: string_cml
 character(len=strln),intent(inout) :: string

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: found_a,found_acell,iatom,ii,index_lower,index_upper,isym
 integer :: lenstr_new,lenstr_old,markuplen,mu,natom,nsym,nu
 integer :: ntypat
 real(dp) :: amu, rcov, znucl
 logical :: hasattrib,hasattrib_fractional,hasattrib_cartesian
 character(len=2) :: string2
 character(len=20) :: string20
 character(len=3) :: string3
 character(len=5) :: string5
 character(len=500) :: message
 character(len=fnlen) :: builtin
!arrays
 integer :: indices_atomArray(3),indices_crystal(3),indices_markup(3)
 integer :: indices_molecule(3),indices_symmetry(3)
 integer,allocatable :: symrel(:,:,:)
 real(dp),allocatable :: tnons(:,:),xred(:,:),xangst(:,:)
 character(len=2),allocatable :: elementtype(:)
 integer, save :: atomspecies(200) = 0
 character(len=500), save :: znuclstring = ""


!************************************************************************

!DEBUG
!write(std_out,*)' append_cml2 : enter , lenstr=',lenstr
!write(std_out,*)trim(string(1:lenstr))
!string(lenstr+1:lenstr+5)=' TEST'
!write(std_out,*)trim(string(1:lenstr+5))
!stop
!ENDDEBUG

 lenstr_new=lenstr

!Find the critical indices of the first 'molecule' mark-up in CML string
 builtin=blank
 index_lower=1
 index_upper=lenstr_cml
 markuplen=8
 call findmarkup(builtin,index_lower,index_upper,indices_molecule,&
& 'molecule',markuplen,strln,string_cml)

!write(std_out,*)string_cml(indices_molecule(1):indices_molecule(2))

 if(indices_molecule(1)>0)then

   write(message,'(a)') ' Identified CML markup <molecule>'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  ---------------------------------------------------------------------------

!  Find the critical indices of the 'crystal' mark-up,
!  inside the first 'molecule' block.
   builtin=blank
   index_lower=indices_molecule(2)
   index_upper=indices_molecule(3)
   markuplen=7
   call findmarkup(builtin,index_lower,index_upper,indices_crystal,&
&   'crystal',markuplen,strln,string_cml)

   if(indices_crystal(1)>0)then

     write(message,'(a)') ' Identified CML markup <crystal>'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     index_lower=indices_crystal(2)
     index_upper=indices_crystal(3)

!    Find <scalar title="a">,<scalar title="b">,<scalar title="c"> ,
!    <scalar title="alpha">,<scalar title="beta">,<scalar title="gamma"> ,
!    and create adequate append
!    WARNING : if a is given, b and c must be given ! suppose that angstroms are used !
!    WARNING : if alpha is given, beta and gamma must be given ! suppose that degrees are used !
     markuplen=6
     found_a=0 ; found_acell=0

     do ii=1,6

!      If did not found 'a', then should not get 'b' or 'c'
       if(found_a==0)then
         if(ii==2 .or. ii==3)cycle
       end if
!      If did not found 'acell', then should not get 'beta' or 'gamma'
       if(found_acell==0)then
         if(ii==5 .or. ii==6)cycle
       end if

       select case (ii)
         case(1)
           builtin='a'
         case(2)
           builtin='b'
         case(3)
           builtin='c'
         case(4)
           builtin='alpha'
         case(5)
           builtin='beta'
         case(6)
           builtin='gamma'
       end select

       call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&       'scalar',markuplen,strln,string_cml)

       if(indices_markup(1)>0)then

         write(message,'(3a)') ' Identified CML markup <scalar title="',trim(builtin),'">'
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,message,'COLL')

         if(ii==1)found_a=1
         if(ii==4)found_acell=1

         lenstr_old=lenstr_new
         if(ii==1)then
           lenstr_new=lenstr_new+7+len_trim(dtset_char)+1+(indices_markup(3)-indices_markup(2)-1)
           string(lenstr_old+1:lenstr_new)=&
&           " _ACELL"//trim(dtset_char)//blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)
         else if(ii==2 .or. ii==5 .or. ii==6)then
           lenstr_new=lenstr_new+1+(indices_markup(3)-indices_markup(2)-1)
           string(lenstr_old+1:lenstr_new)=&
&           blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)
         else if(ii==3)then
           lenstr_new=lenstr_new+1+(indices_markup(3)-indices_markup(2)-1)+9
           string(lenstr_old+1:lenstr_new)=&
&           blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)//" ANGSTROM"
         else if(ii==4)then
           lenstr_new=lenstr_new+8+len_trim(dtset_char)+1+(indices_markup(3)-indices_markup(2)-1)
           string(lenstr_old+1:lenstr_new)=&
&           " _ANGDEG"//trim(dtset_char)//blank//string_cml(indices_markup(2)+1:indices_markup(3)-1)
         end if

       else if(ii/=1 .and. ii/=4)then

         write(message,'(8a)')ch10,&
&         ' append_cml2 : ERROR -',ch10,&
&         '  Could not identify <scalar title="',builtin,'">',ch10,&
&         '  Action : check your CML file ; it is likely that ABINIT is not yet able to read it.'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')

       end if ! Found <scalar title=

     end do ! ii=1,6

   end if ! Found a crystal markup

!  ---------------------------------------------------------------------------

!  Find the critical indices of the 'symmetry' mark-up,
!  inside the first 'molecule' block.
   builtin=blank
   index_lower=indices_molecule(2)
   index_upper=indices_molecule(3)
   markuplen=8
   call findmarkup(builtin,index_lower,index_upper,indices_symmetry,&
&   'symmetry',markuplen,strln,string_cml)

   if(indices_symmetry(1)>0)then

     write(message,'(a)') ' Identified CML markup <symmetry>'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     index_lower=indices_symmetry(2)
     index_upper=indices_symmetry(3)

!    Count the number of <matrix>
     builtin=blank
     markuplen=6
     nsym=0

     do
       call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&       'matrix',markuplen,strln,string_cml)
       if(indices_markup(1)>0)then
         nsym=nsym+1
         index_lower=indices_markup(3)
       else
         exit
       end if
     end do

     if(nsym/=0)then

       write(message,'(a,i5,a)') ' Found',nsym,' symmetry operations ; translate them.'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')

       write(string3,'(i3)')nsym
       lenstr_old=lenstr_new
       lenstr_new=lenstr_new+6+len_trim(dtset_char)+1+3
       string(lenstr_old+1:lenstr_new)=&
&       " _NSYM"//trim(dtset_char)//blank//string3

       ABI_ALLOCATE(symrel,(3,3,nsym))
       ABI_ALLOCATE(tnons,(3,nsym))

       index_lower=indices_symmetry(2)

!      Read symrel and tnons from CML string
       do isym=1,nsym
         call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&         'matrix',markuplen,strln,string_cml)
         index_lower=indices_markup(3)
         read(string_cml(indices_markup(2)+1:indices_markup(3)-1),*) &
&         symrel(1,1:3,isym),tnons(1,isym),&
&         symrel(2,1:3,isym),tnons(2,isym),&
&         symrel(3,1:3,isym),tnons(3,isym)
       end do

!      Write symrel
       lenstr_old=lenstr_new
       lenstr_new=lenstr_new+8+len_trim(dtset_char)
       string(lenstr_old+1:lenstr_new)=" _SYMREL"//trim(dtset_char)
       do isym=1,nsym
         do mu=1,3
           do nu=1,3
             write(string2,'(i2)')symrel(mu,nu,isym)
             lenstr_old=lenstr_new
             lenstr_new=lenstr_new+3
             string(lenstr_old+1:lenstr_new)=blank//string2
           end do
         end do
       end do

!      Write tnons
       lenstr_old=lenstr_new
       lenstr_new=lenstr_new+7+len_trim(dtset_char)
       string(lenstr_old+1:lenstr_new)=" _TNONS"//trim(dtset_char)
       do isym=1,nsym
         do mu=1,3
           write(string20,'(f20.12)')tnons(mu,isym)
           lenstr_old=lenstr_new
           lenstr_new=lenstr_new+21
           string(lenstr_old+1:lenstr_new)=blank//string20
         end do
       end do

       ABI_DEALLOCATE(symrel)
       ABI_DEALLOCATE(tnons)

     end if ! nsym/=0

   end if ! Found <symmetry>

!  ---------------------------------------------------------------------------

!  Find the critical indices of the 'atomArray' mark-up,
!  inside the first 'molecule' block.
   builtin=blank
   index_lower=indices_molecule(2)
   index_upper=indices_molecule(3)
   markuplen=9
   call findmarkup(builtin,index_lower,index_upper,indices_atomArray,&
&   'atomArray',markuplen,strln,string_cml)

   if(indices_atomArray(1)>0)then

     write(message,'(a)') ' Identified CML markup <atomArray>'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     index_lower=indices_atomArray(2)
     index_upper=indices_atomArray(3)

!    Count the number of <atom>
     builtin=blank
     markuplen=4
     natom=0

!    DEBUG
!    write(std_out,'(a)') 'string_cml(index_lower:index_upper)='
!    write(std_out,'(3a)') '"',string_cml(index_lower:index_upper),'"'
!    ENDDEBUG

     do
       call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&       'atom',markuplen,strln,string_cml)
       if(indices_markup(1)>0)then
         natom=natom+1
         index_lower=indices_markup(2)
       else
         exit
       end if
     end do

     if(natom/=0)then

       write(message,'(a,i5,a)') ' Found',natom,' atoms'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')

       write(string5,'(i5)')natom
       lenstr_old=lenstr_new
       lenstr_new=lenstr_new+7+len_trim(dtset_char)+1+5
       string(lenstr_old+1:lenstr_new)=" _NATOM"//trim(dtset_char)//blank//string5


!      Read once more the <atom ...> markups, and find their attribute
       index_lower=indices_atomArray(2)
!      First see whether coordinates are given in reduced or cartesian
!      coordinates in CML file
       call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&       'atom',markuplen,strln,string_cml)
       call getattribute('xFract',6,hasattrib_fractional,indices_markup,strln,string_cml,string20)
       call getattribute('x3',2,hasattrib_cartesian,indices_markup,strln,string_cml,string20)

       if(hasattrib_fractional)then
         ABI_ALLOCATE(elementtype,(natom))
         ABI_ALLOCATE(xred,(3,natom))
         write(message,'(4a)')ch10,&
&         ' append_cml2 : COMMENT -',ch10,&
&         '  Atomic positions are given in fractional coordinates in CML file'
         call wrtout(std_out,message,'COLL')
         do iatom=1,natom
           call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&           'atom',markuplen,strln,string_cml)

           call getattribute('elementType',11,hasattrib,indices_markup,strln,string_cml,string2)
           if(hasattrib)then
             read(string2,*)elementtype(iatom)
           else
             write(message,'(4a,i4,4a)')ch10,&
&             ' append_cml2 : ERROR -',ch10,&
&             '  attribute "elementType" does not exists for atom number ',iatom,ch10,&
&             '  Action : check your CML file if its correct;',ch10,&
&             '           it is likely that ABINIT is not yet able to read it.'
             call wrtout(std_out,message,'COLL')
             call wrtout(ab_out,message,'COLL')
             call leave_new('COLL')
           end if
           call getattribute('xFract',6,hasattrib,indices_markup,strln,string_cml,string20)
           if(hasattrib)then
             read(string20,*)xred(1,iatom)
           else
             write(message,'(4a,i4,4a)')ch10,&
&             ' append_cml2 : ERROR -',ch10,&
&             '  attribute "xFract" does not exists for atom number ',iatom,ch10,&
&             '  Action : check your CML file if its correct;',ch10,&
&             '           it is likely that ABINIT is not yet able to read it.'
             call wrtout(std_out,message,'COLL')
             call wrtout(ab_out,message,'COLL')
             call leave_new('COLL')
           end if
           call getattribute('yFract',6,hasattrib,indices_markup,strln,string_cml,string20)
           if(hasattrib)then
             read(string20,*)xred(2,iatom)
           else
             write(message,'(4a,i4,4a)')ch10,&
&             ' append_cml2 : ERROR -',ch10,&
&             '  attribute "yFract" does not exists for atom number ',iatom,ch10,&
&             '  Action : check your CML file if its correct;',ch10,&
&             '           it is likely that ABINIT is not yet able to read it.'
             call wrtout(std_out,message,'COLL')
             call wrtout(ab_out,message,'COLL')
             call leave_new('COLL')
           end if
           call getattribute('zFract',6,hasattrib,indices_markup,strln,string_cml,string20)
           if(hasattrib)then
             read(string20,*)xred(3,iatom)
           else
             write(message,'(4a,i4,4a)')ch10,&
&             ' append_cml2 : ERROR -',ch10,&
&             '  attribute "zFract" does not exists for atom number ',iatom,ch10,&
&             '  Action : check your CML file if its correct;',ch10,&
&             '           it is likely that ABINIT is not yet able to read it.'
             call wrtout(std_out,message,'COLL')
             call wrtout(ab_out,message,'COLL')
             call leave_new('COLL')
           end if

           index_lower=indices_markup(2)
         end do
       else if(hasattrib_cartesian)then
         ABI_ALLOCATE(elementtype,(natom))
         ABI_ALLOCATE(xangst,(3,natom))
         write(message,'(4a)')ch10,&
&         ' append_cml2 : COMMENT -',ch10,&
&         '  Atomic positions are given in cartesian coordinates in CML file'
         call wrtout(std_out,message,'COLL')
         do iatom=1,natom
           call findmarkup(builtin,index_lower,index_upper,indices_markup,&
&           'atom',markuplen,strln,string_cml)

           call getattribute('elementType',11,hasattrib,indices_markup,strln,string_cml,string2)
           if(hasattrib)then
             read(string2,*)elementtype(iatom)
           else
             write(message,'(4a,i4,4a)')ch10,&
&             ' append_cml2 : ERROR -',ch10,&
&             '  attribute "elementType" does not exists for atom number ',iatom,ch10,&
&             '  Action : check your CML file if its correct;',ch10,&
&             '           it is likely that ABINIT is not yet able to read it.'
             call wrtout(std_out,message,'COLL')
             call wrtout(ab_out,message,'COLL')
             call leave_new('COLL')
           end if
           call getattribute('x3',2,hasattrib,indices_markup,strln,string_cml,string20)
           if(hasattrib)then
             read(string20,*)xangst(1,iatom)
           else
             write(message,'(4a,i4,4a)')ch10,&
&             ' append_cml2 : ERROR -',ch10,&
&             '  attribute "x3" does not exists for atom number ',iatom,ch10,&
&             '  Action : check your CML file if its correct;',ch10,&
&             '           it is likely that ABINIT is not yet able to read it.'
             call wrtout(std_out,message,'COLL')
             call wrtout(ab_out,message,'COLL')
             call leave_new('COLL')
           end if
           call getattribute('y3',2,hasattrib,indices_markup,strln,string_cml,string20)
           if(hasattrib)then
             read(string20,*)xangst(2,iatom)
           else
             write(message,'(4a,i4,4a)')ch10,&
&             ' append_cml2 : ERROR -',ch10,&
&             '  attribute "y3" does not exists for atom number ',iatom,ch10,&
&             '  Action : check your CML file if its correct;',ch10,&
&             '           it is likely that ABINIT is not yet able to read it.'
             call wrtout(std_out,message,'COLL')
             call wrtout(ab_out,message,'COLL')
             call leave_new('COLL')
           end if
           call getattribute('z3',2,hasattrib,indices_markup,strln,string_cml,string20)
           if(hasattrib)then
             read(string20,*)xangst(3,iatom)
           else
             write(message,'(4a,i4,4a)')ch10,&
&             ' append_cml2 : ERROR -',ch10,&
&             '  attribute "z3" does not exists for atom number ',iatom,ch10,&
&             '  Action : check your CML file if its correct;',ch10,&
&             '           it is likely that ABINIT is not yet able to read it.'
             call wrtout(std_out,message,'COLL')
             call wrtout(ab_out,message,'COLL')
             call leave_new('COLL')
           end if

           index_lower=indices_markup(2)
         end do
       else
         write(message,'(4a,i4,4a)')ch10,&
&         ' append_cml2 : ERROR -',ch10,&
&         '  Could not find if atomic positions are present in CML file,',ch10,&
&         '  neither in fractional nor in cartesian format ...',ch10,&
&         '  Action : check your CML file if its correct;',ch10,&
&         '           it is likely that ABINIT is not yet able to read it.'
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,message,'COLL')
         call leave_new('COLL')
       end if

!      extract znucl for each atom type
       do iatom = 1, natom
         call symbol2znucl(amu,rcov,elementtype(iatom),znucl)
         if (znucl > 200) then
           write (message,'(6a)')ch10,&
&           ' Error: found element beyond Z=200 ', ch10,&
&           ' Solution: increase size of atomspecies in append_cml2', ch10
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if
!        found a new atom type 
         if (atomspecies(int(znucl)) == 0) then
           write(string20,'(f10.2)') znucl
           znuclstring = trim(znuclstring) // " " // trim(string20) // " "
         end if
         atomspecies(int(znucl)) = 1
       end do

!      Write the element type
       lenstr_old=lenstr_new
       lenstr_new=lenstr_new+7+len_trim(dtset_char)+1
       string(lenstr_old+1:lenstr_new)=" _TYPAX"//trim(dtset_char)//blank

       do iatom=1,natom
         lenstr_old=lenstr_new
         lenstr_new=lenstr_new+3
         string(lenstr_old+1:lenstr_new)=elementtype(iatom)//blank
       end do

       lenstr_old=lenstr_new
       lenstr_new=lenstr_new+3
       string(lenstr_old+1:lenstr_new)="XX "

!      Write the fractional or cartesian coordinates (in angstroms)
       if(hasattrib_fractional)then
         lenstr_old=lenstr_new
         lenstr_new=lenstr_new+6+len_trim(dtset_char)+1
         string(lenstr_old+1:lenstr_new)=" _XRED"//trim(dtset_char)//blank

         do iatom=1,natom
           do mu=1,3
             write(string20,'(f20.12)')xred(mu,iatom)
             lenstr_old=lenstr_new
             lenstr_new=lenstr_new+20
             string(lenstr_old+1:lenstr_new)=string20
           end do
         end do

         ABI_DEALLOCATE(elementtype)
         ABI_DEALLOCATE(xred)

       else if(hasattrib_cartesian)then
         lenstr_old=lenstr_new
         lenstr_new=lenstr_new+8+len_trim(dtset_char)+1
         string(lenstr_old+1:lenstr_new)=" _XANGST"//trim(dtset_char)//blank

         do iatom=1,natom
           do mu=1,3
             write(string20,'(f20.12)')xangst(mu,iatom)
             lenstr_old=lenstr_new
             lenstr_new=lenstr_new+20
             string(lenstr_old+1:lenstr_new)=string20
           end do
         end do

         ABI_DEALLOCATE(elementtype)
         ABI_DEALLOCATE(xangst)

       end if ! cartesian or reduced coordinates

     end if ! natom/=0

   end if ! Found <atomArray>

!  ---------------------------------------------------------------------------

 end if ! Found <molecule>

!if we are on the last dataset we can write znucl and ntypat
 if (dtset_char == "-1") then
!  write znucl
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+7+len_trim(znuclstring)+1
   string(lenstr_old+1:lenstr_new)=" ZNUCL"//blank//trim(znuclstring)//blank

!  write ntypat
   ntypat = sum(atomspecies)
   write(string20,'(i10)') ntypat
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+8+len_trim(string20)+1
   string(lenstr_old+1:lenstr_new)=" NTYPAT"//blank//trim(string20)//blank
 end if

!Check the length of the string
 if(lenstr_new>strln)then
   write(message,'(6a)')ch10,&
&   ' append_cml2 : BUG -',ch10,&
&   '  The maximal size of the input variable string has been exceeded.',ch10,&
&   '  The use of a CML file is more character-consuming than the usual input file. Sorry.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Update the length of the string
 lenstr=lenstr_new

!DEBUG
 write(std_out,*)' append_cml2 : exit , lenstr=',lenstr
 write(std_out,*)trim(string(1:lenstr))
!stop
!ENDDEBUG

end subroutine append_cml2
!!***
