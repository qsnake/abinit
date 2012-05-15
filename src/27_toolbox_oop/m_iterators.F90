!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_iterators
!! NAME
!!  m_iterators
!!
!! FUNCTION
!!  This module defines objects (iterators) that are used to facilitate the 
!!  iteration over the elements of an ensemble e.g. set of transitions.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  Iterators are used to obtain an indirect indexing of the indeces used to access array elements 
!!  or array sections. An iterator is an object that contains the list of indeces that should
!!  be treated. Each processor has its own local version initialized according the the distribution 
!!  of the data across the nodes.
!!  Using iterators for looping facilitates the implementation of MPI algorithms since all the 
!!  information on the distribution of the tasks is stored in the iterator itself
!!  For example an MPI loop over spins, k-points and bands can be writted in terms of iterators using:
!!
!!   do isp=1,iter_len(Iter_bks)
!!     spin = iter_yield(Iter_bks,idx3=isp)
!!   
!!     do ikp=1,iter_len(Iter_bks,idx3=isp)
!!       ik_ibz = iter_yield(Iter_bks,idx2=ikp,idx3=isp)
!!   
!!       do ibn=1,iter_len(Iter_bks,idx2=ikp,idx3=isp)    
!!         band = iter_yield(Iter_bks,entry1=ibn,idx2=ikp,idx3=isp)
!!   
!!  iter_len gives the number of non-zero entries 
!!  iter_yield returns the global indeces used to access the data.
!!
!!  The main advantage is seen in sections of code that are MPI parallelized because each node
!!  will have its onw iterator (previously initialized by the programmer). 
!!  Therefore it is very easy to distribute the load among the nodes without having to 
!!  to use "cycle" of "if then" instruction in the inners loops.
!!  Another important advantage is that the MPI implementation will continue to work even 
!!  if the data distribution is changed, only the iterator has to be modified. 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_iterators

 use m_profiling

 use defs_basis
 use m_errors

 implicit none

 private

!!***

!----------------------------------------------------------------------

!Structures
!!***


!!****t* m_iterators/indeces_t
!! NAME
!!  indeces_t
!!
!! FUNCTION
!!  A base datatype for constructing ragged arrays.
!!
!! SOURCE

 type,private :: indeces_t
  integer :: leng
  integer, pointer :: indx(:)   SET2NULL
  ! indx(leng) The set of indeces.
 end type indeces_t
!!***

!----------------------------------------------------------------------

!!****t* m_iterators/iter2_t
!! NAME
!!  iter2_t
!!
!! FUNCTION
!!
!! SOURCE

 type,public :: iter2_t
   private
   integer :: sizes(2) 
   integer :: starts(2)
   type(indeces_t),pointer :: slice(:,:)  SET2NULL
   ! Temporary structure to store the indeces in full storage mode.

#if 0
   integer :: len3
   ! The number of non-zero entries in the last dimension.

   integer,pointer :: map3(:)
   ! map3(len3) 
   ! Indirect indexing packed --> full for the last dimension.

   integer,pointer :: len2(:)  SET2NULL     
   ! len2(len3)
   ! Gives the number of non-zero entries along the second dimension for each non-null index along the 3-dimension.

   type(indeces_t),pointer :: map2(:)  SET2NULL     
   ! map2(len3)%indeces
   ! Indirect indexing packed --> full for the second dimension.

   type(indeces_t),pointer :: len1(:,:)  SET2NULL     
   ! len1(MAX(len2),len3)
   ! Gives the number of non-zero entries along the first dimension.

   type(indeces_t),pointer :: map1(:)  SET2NULL     
   ! map1(MAX(len2),len3)%indeces
   ! Indirect indexing packed --> full for the first dimension.
#endif

 end type iter2_t
!!***

!----------------------------------------------------------------------

!!****t* m_iterators/iter3_t
!! NAME
!!  iter3_t
!!
!! FUNCTION
!!
!! SOURCE
!!
!! type,public :: iter3_t
!!   private
!!   integer :: sizes(3) 
!!   integer :: starts(3)
!!   type(indeces_t),pointer :: slice(:,:,:)  SET2NULL
!! end type iter3_t
!!***

!----------------------------------------------------------------------

 public :: iter_alloc         ! Allocate the iterator.
 public :: iter_push          ! Copy a set of indeces in the iterator 
 public :: iter_destroy       ! Deallocate the iterator.
 public :: OPERATOR(.LBOUND.) ! 
 public :: OPERATOR(.UBOUND.) ! 
 public :: OPERATOR(.SIZE.)   ! 
 public :: iter_len           ! The number of indeces in a slice of the iterator
 public :: yield              ! Return the indeces in of the slice of the iterator.
 public :: iter_print         ! Printout of the iterator, just for debugging purposes.

 interface iter_alloc
   module procedure iter2_alloc
   !module procedure iter3_alloc
 end interface iter_alloc

 interface iter_push
   module procedure iter2_push
   !module procedure iter3_push
 end interface iter_push

 interface iter_destroy
   module procedure iter2_destroy
   !module procedure iter3_destroy
 end interface iter_destroy

 interface operator(.lbound.)
   module procedure iter2_lbound
   !module procedure iter3_lbound
 end interface 

 interface operator(.ubound.)
   module procedure iter2_ubound
   !module procedure iter3_ubound
 end interface

 interface operator(.size.)
   module procedure iter2_size
   !module procedure iter3_size
 end interface 

 interface iter_len
   module procedure iter2_len
   !module procedure iter3_len
 end interface iter_len

 interface yield
   module procedure iter2_yield
   !module procedure iter3_yield
 end interface yield

 interface iter_print
   module procedure iter2_print
   !module procedure iter3_print
 end interface iter_print

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_iterators/indeces_destroy
!! NAME
!!  indeces_destroy
!!
!! FUNCTION
!!
!! INPUTS 
!!
!! OUTPUT
!!
!! PARENTS
!!      m_iterators
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine indeces_destroy(Ids)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'indeces_destroy'
!End of the abilint section

 type(indeces_t),intent(inout) :: Ids

! *************************************************************************
 
 Ids%leng=0
 if (associated(Ids%indx))  then
   ABI_DEALLOCATE(Ids%indx)
 end if

end subroutine indeces_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_iterators/iter2_alloc
!! NAME
!!  iter2_alloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine iter2_alloc(Iter2,sizes,starts)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iter2_alloc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(iter2_t),intent(inout) :: Iter2
!arrays
 integer,intent(in) :: sizes(2)
 integer,optional,intent(in) :: starts(2)

!Local variables ------------------------------
!scalars
 integer :: s1,s2,i1,i2

!************************************************************************
 s1=1; s2=1 
 if (PRESENT(starts)) then 
   s1=starts(1)
   s2=starts(2)
 end if

 Iter2%starts=(/s1,s2/)
 Iter2%sizes =sizes

 ABI_ALLOCATE( Iter2%slice,(s1:s1+sizes(1)-1, s2:s2+sizes(2)-1))

 do i2=LBOUND(Iter2%slice,DIM=2),UBOUND(Iter2%slice,DIM=2)
   do i1=LBOUND(Iter2%slice,DIM=1),UBOUND(Iter2%slice,DIM=1)
     Iter2%slice(i1,i2)%leng=0
     nullify(Iter2%slice(i1,i2)%indx)
   end do
 end do

end subroutine iter2_alloc
!!***

!----------------------------------------------------------------------

!!****f* m_iterators/iter2_push
!! NAME
!!  iter2_push
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine iter2_push(Iter2,i1,i2,list)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iter2_push'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: i1,i2
 type(iter2_t),intent(inout) :: Iter2
!arrays
 integer,intent(in) :: list(:)

!Local variables ------------------------------
!scalars
 integer :: leng

!************************************************************************

 leng = SIZE(list)

 if (associated( Iter2%slice(i1,i2)%indx) ) then
   MSG_ERROR("Iter2%slice already allocated")
 end if

 Iter2%slice(i1,i2)%leng = leng
 ABI_ALLOCATE(Iter2%slice(i1,i2)%indx,(leng))
 Iter2%slice(i1,i2)%indx = list

end subroutine iter2_push
!!***

!----------------------------------------------------------------------

!!****f* m_iterators/iter2_destroy
!! NAME
!!  iter2_destroy
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine iter2_destroy(Iter2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iter2_destroy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(iter2_t),intent(inout) :: Iter2

!Local variables ------------------------------
!scalars
 integer :: i1,i2

!************************************************************************

 do i2=LBOUND(Iter2%slice,DIM=2),UBOUND(Iter2%slice,DIM=2)
   do i1=LBOUND(Iter2%slice,DIM=1),UBOUND(Iter2%slice,DIM=1)
     call indeces_destroy(Iter2%slice(i1,i2))
   end do
 end do

 ABI_DEALLOCATE(Iter2%slice)

end subroutine iter2_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_iterators/iter2_len
!! NAME
!!  iter2_len
!!
!! FUNCTION
!!
!! INPUTS 
!!
!! OUTPUT
!!
!! SOURCE

function iter2_len(Iter2,i1,i2)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iter2_len'
!End of the abilint section

 integer,intent(in) :: i1,i2
 integer :: iter2_len
 type(iter2_t),intent(in) :: Iter2

! *************************************************************************

 iter2_len = Iter2%slice(i1,i2)%leng

end function iter2_len
!!***

!----------------------------------------------------------------------

!!****f* m_iterators/iter2_lbound
!! NAME
!!  iter2_lbound
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function iter2_lbound(Iter2,dim)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iter2_lbound'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dim
 integer :: iter2_lbound
 type(iter2_t),intent(in) :: Iter2

!************************************************************************

 iter2_lbound = Iter2%starts(dim)

end function iter2_lbound
!!***

!----------------------------------------------------------------------

!!****f* m_iterators/iter2_ubound
!! NAME
!!  iter2_ubound
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function iter2_ubound(Iter2,dim)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iter2_ubound'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dim
 integer :: iter2_ubound
 type(iter2_t),intent(in) :: Iter2

!************************************************************************

 iter2_ubound = Iter2%starts(dim) + Iter2%sizes(dim) -1

end function iter2_ubound
!!***

!----------------------------------------------------------------------

!!****f* m_iterators/iter2_size
!! NAME
!!  iter2_size
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function iter2_size(Iter2,dim)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iter2_size'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dim
 integer :: iter2_size
 type(iter2_t),intent(in) :: Iter2

!************************************************************************

 iter2_size = iter2_ubound(Iter2,dim) - iter2_lbound(Iter2,dim) + 1

end function iter2_size
!!***

!----------------------------------------------------------------------

!!****f* m_iterators/iter2_yield
!! NAME
!!  iter2_yield
!!
!! FUNCTION
!!
!! INPUTS 
!!
!! OUTPUT
!!
!! SOURCE

function iter2_yield(Iter2,idx,i1,i2)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iter2_yield'
!End of the abilint section

 integer,intent(in) :: idx,i1,i2
 integer :: iter2_yield
 type(iter2_t),intent(in) :: Iter2

! *************************************************************************

 iter2_yield = Iter2%slice(i1,i2)%indx(idx)

end function iter2_yield
!!***

!----------------------------------------------------------------------

!!****f* m_iterators/iter2_print
!! NAME
!!  iter2_print
!!
!! FUNCTION
!!
!! INPUTS 
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine iter2_print(Iter2,header,unit,mode_paral,prtvol) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iter2_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral 
 character(len=*),optional,intent(in) :: header
 type(iter2_t),intent(in) :: Iter2

!Local variables-------------------------------
 integer :: my_unt,my_prtvol,ntot,i1,i2,idx
 character(len=4) :: my_mode
 character(len=500) :: msg      
! ********************************************************************* 

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol 
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Content of the iter2_t object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 ntot = PRODUCT(Iter2%sizes)
 write(std_out,*)"total number of elements:",ntot
 write(std_out,*)"list of indeces: "

 do i2=iter2_lbound(Iter2,DIM=2),iter2_ubound(Iter2,DIM=2)
   do i1=iter2_lbound(Iter2,DIM=1),iter2_lbound(Iter2,DIM=1)
      write(std_out,*) (yield(Iter2,idx,i1,i2), idx=1,iter_len(Iter2,i1,i2))
   end do
 end do

end subroutine iter2_print

END MODULE m_iterators
!!***
