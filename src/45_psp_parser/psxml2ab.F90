!{\src2tex{textfont=tt}}
!!****f* ABINIT/psxml2ab
!! NAME
!! psxml2ab
!!
!! FUNCTION
!!  From a SIESTA XML format pseudopotential file which has already been read in,
!!  convert to abinit internal datastructures.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! psxml  = pseudopotential data structure
!!
!! OUTPUT
!! znucl  = atomic charge
!! zion   = valence charge
!! pspcod = index for the pseudopotential flavor type (9 for XML)
!! pspxc  = exchange and correlation functional
!! lmax   = maximum angular momentum shell
!! iwrite = flag that controls whether to dump pseudo info in the output
!!          iwrite = 1, yes ; iwrite = 0, no
!!
!! PARENTS
!!      inpspheads,pspatm
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine psxml2ab( psxml, znucl, zion, pspcod, pspxc, lmax, iwrite )

 use m_profiling

 use defs_basis
 use m_xml_pseudo_types

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psxml2ab'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iwrite
 integer,intent(out) :: lmax,pspcod,pspxc
 real(dp),intent(out) :: zion,znucl
 type(pseudo_t),intent(in) :: psxml

!Local variables-------------------------------
!scalars
 integer :: il,lmm,lshellint,position
 logical :: polarized
 character(len=1) :: ispp,lshell
 character(len=500) :: message
 character(len=70) :: text
!arrays
 real(dp) :: zeld(0:4),zelu(0:4)

! *********************************************************************

 znucl  = psxml%header%atomicnumber
 zion   = psxml%header%zval

 pspcod = 9

 select case(psxml%header%xcfunctionalparametrization)
   case('Ceperley-Alder')
     pspxc = 2
   case('Wigner')
     pspxc = 4
   case('Hedin-Lundqvist')
     pspxc = 5
   case('Gunnarson-Lundqvist')
     write(message, '(a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&     ' psxml2ab : ERROR -',ch10,&
&     '  The exchange and correlation potential by Gunnarson-Lundqvist', ch10,&
&     '  is not implemented in Abinit.',ch10,&
&     '  Action : choose another exchange and correlation potential ',ch10,&
&     '  in the pseudopotential generation. ;'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   case('von Barth-Hedin')
     write(message, '(a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&     ' psxml2ab  : ERROR -',ch10,&
&     '  The exchange and correlation potential by von Barth-Hedin', ch10, &
&     '  is not implemented in Abinit.',ch10,&
&     '  Action : choose another exchange and correlation potential ',ch10,&
&     '  in the pseudopotential generation. ;'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   case('Perdew-Burke-Ernzerhof')
     pspxc = 11
   case('Becke-Lee-Yang-Parr')
     write(message, '(a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&     ' psxml2ab  : ERROR -',ch10,&
&     '  The exchange and correlation potential by Becke-Lee-Yang-Parr', ch10,&
&     '  is not implemented in Abinit.',ch10,&
&     '  Action : choose another exchange and correlation potential ',ch10, &
&     '  in the pseudopotential generation. ;'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
 end select

 lmax = 0

 do il = 1, psxml%npots_down
   lshell = psxml%pot(il)%l
   select case(lshell)
     case('s')
       lshellint = 0
       lmax = max( lmax, lshellint )
     case('p')
       lshellint = 1
       lmax = max( lmax, lshellint )
     case('d')
       lshellint = 2
       lmax = max( lmax, lshellint )
     case('f')
       lshellint = 3
       lmax = max( lmax, lshellint )
   end select
   zeld(lshellint) = psxml%pot(il)%occupation
 end do

 do il = 1, psxml%npots_up
   lshell = psxml%pot(psxml%npots_down+il)%l
   select case(lshell)
     case('s')
       lshellint = 0
       lmax = max( lmax, lshellint )
     case('p')
       lshellint = 1
       lmax = max( lmax, lshellint )
     case('d')
       lshellint = 2
       lmax = max( lmax, lshellint )
     case('f')
       lshellint = 3
       lmax = max( lmax, lshellint )
   end select
   zelu(lshellint) = psxml%pot(psxml%npots_down+il)%occupation
 end do

 select case(psxml%header%relativistic)
   case(.true.)
     ispp      = 'r'
     polarized = .false.
   case(.false.)
     select case(psxml%header%polarized)
       case(.true.)
         ispp      = 's'
         polarized = .true.
       case(.false.)
         ispp      = ' '
         polarized = .false.
     end select
 end select

 text     = ' '
 position = 1
 lmm     = max( psxml%npots_down, psxml%npots_up )
 do il = 1, lmm
   if ( .not. polarized) then
     write(text(position:),9070) &
&     psxml%pot(il)%n,                        &
&     psxml%pot(il)%l,                        &
&     zeld(il-1)+zelu(il-1),                  &
&     ispp,                                   &
&     psxml%pot(il)%cutoff
     9070 format(i1,a1,f5.2,a1,' r=',f5.2,'/')
     position = position + 17
   else
     write(text(position:),9090) &
&     psxml%pot(il)%n,                        &
&     psxml%pot(il)%l,                        &
&     zeld(il-1),                             &
&     zelu(il-1),                             &
&     ispp,                                   &
&     psxml%pot(il)%cutoff
     9090 format(i1,a1,f5.2,',',f5.2,a1,f5.2,'/')
     position = position + 17
   end if
 end do

 write(std_out,'(a,f5.1,a,i4,a,i4)' ) '  read the values zionpsp=',&
& zion,' , pspcod=',pspcod,' , lmax=',lmax

 if ( iwrite .eq. 1 ) then

   write(message,'(a,a)') &
&   '- psxml2ab: Atomic Label:                      ', &
&   psxml%header%symbol
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,f12.5)') &
&   '- psxml2ab: Atomic Number:                     ', &
&   psxml%header%atomicnumber
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,f12.5)') &
&   '- psxml2ab: Valence charge:                    ', &
&   psxml%header%zval
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a)') &
&   '- psxml2ab: Pseudopotential generator code :    ', &
&   psxml%header%creator
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a)') &
&   '- psxml2ab: Date of pseudopotential generation: ', &
&   psxml%header%date
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a)') &
&   '- psxml2ab: Pseudopotential flavor:             ', &
&   psxml%header%flavor
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a)') &
&   '- psxml2ab: Exchange-correlation functional:    ', &
&   psxml%header%xcfunctionaltype
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a)') &
&   '- psxml2ab: Exchange and correlation parametrization: ', &
&   psxml%header%xcfunctionalparametrization
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a,a)') &
&   '- psxml2ab: Reference configuration:            ',ch10,&
&   text
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,l3)') &
&   '- psxml2ab: Relativistic pseudopotential:       ', &
&   psxml%header%relativistic
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,l3)') &
&   '- psxml2ab: Spin-polarized pseudopotential:     ', &
&   psxml%header%polarized
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   select case(psxml%header%core_corrections)
     case("yes")
       write(message, '(a)' ) &
&       '- psxml2ab: XC core correction read in from XML file.'
!      stop
     case("no")
       write(message, '(a)' ) &
&       '- psxml2ab: No core corrections.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
   end select
 end if

end subroutine psxml2ab
!!***
