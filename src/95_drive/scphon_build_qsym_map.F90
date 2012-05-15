!!****f* ABINIT/scphon_build_qsym_map
!! NAME
!! scphon_build_qsym_map
!!
!! FUNCTION
!! Build up map of transformation of qpoints into one another under symmetry
!! operations, and eventually time reversal.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nphononq=number of phonon q-vectors input from anaddb run at equilibrium
!!   geometry
!! nsym_primitive_cell= number of symmetries in the primitive unit cell
!! phononq= phonon q vectors used in anaddb run (reduced coordinates)
!! symrec_primitive_cell= reciprocal space symmetry operations
!!
!! OUTPUT
!! qsym_map= map of qpoints onto one another, by sym ops:
!!   q_{qmap(iq,isym)} = S_{isym} q_{iq}
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!      wrap2_pmhalf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scphon_build_qsym_map(nphononq,nsym_primitive_cell,phononq,&
&    qsym_map,symrec_primitive_cell)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scphon_build_qsym_map'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nphononq,nsym_primitive_cell
!arrays
 integer,intent(in) :: symrec_primitive_cell(3,3,nsym_primitive_cell)
 integer,intent(out) :: qsym_map(nphononq,nsym_primitive_cell,2)
 real(dp),intent(in) :: phononq(3,nphononq)

!Local variables-------------------------------
!scalars
 integer :: dummy_prtvol,found,iq,iq_image,isym,itimrev,maxtimrev
 real(dp) :: shift,timrev_sign
 character(len=500) :: message
!arrays
 real(dp) :: image_qpoint(3),tmpvec(3)

! *************************************************************************

 dummy_prtvol=0
 qsym_map=0

!write(std_out,'(9I6)') symrec_primitive_cell

 do iq=1,nphononq

!  find small group of qpoint: returns
!  symq(4,2,nsym)= (integer) three first numbers define the G vector ;
!  fourth number is zero if the q-vector is not preserved, is 1 otherwise
!  second index is one without time-reversal symmetry, two with
!  call symq3(nsym_primitive_cell,phononq(:,iq),symq,symrec_primitive_cell,qtimrev,dummy_prtvol)

!  write(std_out,*) 'making qmap: ', iq, sum(symq(4,1,:)), qtimrev
!  if (qtimrev==0) maxtimrev=1
!  if (qtimrev==1) maxtimrev=2

   maxtimrev=2
   timrev_sign=one

!  find image of present qpoint through symmetry isym
   do itimrev=1,maxtimrev
     do isym=1,nsym_primitive_cell
!      if (symq(4,itimrev,isym)==0) cycle
       image_qpoint(:) = timrev_sign*&
&       (dble(symrec_primitive_cell(:,1,isym))*phononq(1,iq)&
&       +dble(symrec_primitive_cell(:,2,isym))*phononq(2,iq)&
&       +dble(symrec_primitive_cell(:,3,isym))*phononq(3,iq))
       call wrap2_pmhalf (image_qpoint(1),tmpvec(1),shift)
       call wrap2_pmhalf (image_qpoint(2),tmpvec(2),shift)
       call wrap2_pmhalf (image_qpoint(3),tmpvec(3),shift)
       found=0
       do iq_image=1,nphononq
!        write(std_out,'(3E20.8)') image_qpoint, phononq(:,iq_image)
!        write(std_out,'(E20.8)') sum((phononq(:,iq_image)-image_qpoint(:))**2)
         if (sum((phononq(:,iq_image)-tmpvec(:))**2) < tol10) then
           found=iq_image
           exit
         end if
       end do
       if (found==0) then
         write(message,'(a,I6,3E20.10)') 'Warning: sym qpoint not found ', isym, image_qpoint
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,    message,'COLL')
         write(message,'(a,a)') ' is your qgrid compatible with all symmetries',&
&         ' of the primitive unit cell?'
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,    message,'COLL')
!        call leave_new('COLL')
       end if
       qsym_map(iq,isym,itimrev) = found
     end do
     timrev_sign=-one
   end do
 end do

end subroutine scphon_build_qsym_map
!!***

