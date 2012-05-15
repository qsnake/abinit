!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_wfutils
!! NAME
!!  m_wfutils
!!
!! FUNCTION
!!  parameters and function for wave functions copy
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2012 ABINIT group (CS,GZ,FB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_wfutils

 use m_profiling

 use defs_basis
 use m_errors

 implicit none

 private

 public :: setWFParameter ! Definition of wave functions parameters for copy functions.
 public :: wfcopy         ! Copy of wave functions arrays.

!Parameters ------------------------------------
 integer,public  :: x_cplx   ! fortran data type for wave
                             ! functions arrays real or complex.
 integer,public  :: x_me_g0  ! processeur number
 integer,public  :: x_npw_k  ! number of plane waves at this k point
 integer,public  :: x_nspinor! number of spinorial components of the wavefunctions on current proc
 integer,public  :: x_icg    ! shift to be applied on the location of data in the array cg
 integer,public  :: x_igsc   ! shift to be applied on the location of data in the array gsc
 integer,public  :: x_blocksize

contains
!!***

!! NAME
!!
!! setWFParameter
!!
!! FUNCTION
!! Definition of wave functions parameters for copy functions
!!
!! COPYRIGHT
!!   Copyright (C) 1998-2012 ABINIT group (FBottin,CS)
!!   this file is distributed under the terms of the
!!   gnu general public license, see ~abinit/COPYING
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   cplx:         fortran data type for wave functions arrays real or complex
!!   me_g0:        processeur number
!!   npw_k:        number of plane waves at this k point
!!   nspinor:      number of spinorial components of the wavefunctions on current proc
!!   icg:          shift to be applied on the location of data in the array cg
!!   igsc:         shift to be applied on the location of data in the array gsc
!!   blocksize     size of blocks
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine setWFParameter(cplx,me_g0,npw_k,nspinor,icg,igsc,blocksize)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setWFParameter'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: cplx
 integer, intent(in) :: me_g0
 integer, intent(in) :: npw_k
 integer, intent(in) :: nspinor
 integer, intent(in) :: icg
 integer, intent(in) :: igsc
 integer, intent(in) :: blocksize

! *********************************************************************

 x_cplx=cplx
 x_me_g0=me_g0
 x_npw_k=npw_k
 x_nspinor=nspinor
 x_icg=icg
 x_igsc=igsc
 x_blocksize=blocksize
 return

end subroutine setWFParameter
!!***

! correspondence with abinit. here for real wf
! this is the index of a given band in cg array
integer function x_cgindex(iblocksize)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'x_cgindex'
!End of the abilint section

 implicit none
 integer, intent(in) :: iblocksize
 x_cgindex=x_npw_k*x_nspinor*(iblocksize-1)+x_icg+1
 return

end function x_cgindex

! correspondence with abinit. here for real wf
! this is the index of a given band in gsc array
integer function x_gscindex(iblocksize)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'x_gscindex'
!End of the abilint section

 implicit none
 integer, intent(in) :: iblocksize
 x_gscindex=x_npw_k*x_nspinor*(iblocksize-1)+x_igsc+1
 return

end function x_gscindex

integer function x_windex(iblocksize)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'x_windex'
!End of the abilint section

 implicit none
 integer, intent(in) :: iblocksize
 x_windex=x_npw_k*x_nspinor*(iblocksize-1)+1
 return

end function x_windex

integer function wfindex(iblocksize,indtype)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfindex'
!End of the abilint section

 implicit none
 integer, intent(in) :: iblocksize
 character, intent(in) :: indtype
 select case(indtype)
 case ('C')
   wfindex=x_cgindex(iblocksize)
 case ('S')
   wfindex=x_gscindex(iblocksize)
 case ('W')
   wfindex=x_windex(iblocksize)
 end select

end function wfindex

!! NAME
!!
!! wfcopy
!!
!! FUNCTION
!! copy of wave functions arrays.
!! called from lobpcg in REAL or COMPLEX wave fonctions.
!!
!! COPYRIGHT
!!   Copyright (C) 1998-2012 ABINIT group (FBottin,CS)
!!   this file is distributed under the terms of the
!!   gnu general public license, see ~abinit/COPYING
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   direction:    copy direction 'D' direct (global array to local)
!!                'I' indirect (local to global)
!!   size:         number of elements
!!   tsrc:         source array
!!   incsrc:       size incrmentation for tsrc array
!!   tdest:        destination array
!!   incdest:      size incrmentation for tdest array
!!   blockiter:    number of block iteration in case REAL
!!   iblock:       block index
!!   indtype:      indexation type in array
!!   withbbloc:    apply bloc on band for each block
!!
!! PARENTS
!!   lobpcgwf
!!
!! CHILDREN
!!   dcopy, zcopy
!!
!! SOURCE
!!
subroutine wfcopy(direction,size,tsrc,incsrc,tdest,incdest,blockiter,iblock,indtype,&
&                withbbloc,timopt,tim_wfcopy) ! optional arguments

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfcopy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character, intent(in) :: direction
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 real(dp), dimension(:,:) :: tsrc
 real(dp), dimension(:,:) :: tdest
 integer, intent(in) :: blockiter
 integer, intent(in) :: iblock
 character, intent(in) :: indtype
 logical, optional, intent(in) :: withbbloc
 integer, intent(in), optional :: timopt,tim_wfcopy
!Local variables ------------------------------------
 logical :: bblock=.false.
 integer :: lig
 integer :: vectsize
 integer :: rvectsize
 integer :: blocksize
 integer :: bblocksize
 integer :: iblocksize
 integer :: iband
 real(dp) :: factor
 real(dp) :: tsec(2)

! *********************************************************************

 if (present(tim_wfcopy).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_wfcopy,1,tsec)
 end if
 if ( present(withbbloc) ) then
   bblock=withbbloc
 endif
 if ( indtype == 'C' ) then
   lig=x_icg
 else if ( indtype == 'S' ) then
   lig=x_igsc
 else
   lig=0
 endif
 rvectsize=x_npw_k*x_nspinor
 if ( x_me_g0 == 1 ) then
   vectsize=2*rvectsize-1
 else
   vectsize=2*rvectsize
 endif
 if ( x_cplx == 2 ) then
   vectsize=x_npw_k*x_nspinor
 endif
 blocksize  = x_blocksize
 bblocksize=(iblock-1)*blocksize
 if ( direction == 'D' ) then
   if ( x_cplx == 1 ) then
     factor=sqrt(2.0_dp)
     do iblocksize=1,blockiter
       iband=iblocksize
       if ( bblock ) then
         iband=iblocksize+bblocksize
       endif
       if ( x_me_g0 == 1 ) then
         !call dcopy(1,tsrc(1,wfindex(iband,indtype)),1,tdest(1,iblocksize),1)
         tdest(1 ,iblocksize)=tsrc(1,wfindex(iband,indtype))
         call dcopy(rvectsize-1,tsrc(1,wfindex(iband,indtype)+1:wfindex(iband+1,indtype)-1)*factor,&
&                   1,tdest(2:rvectsize         ,iblocksize),1)
         call dcopy(rvectsize-1,tsrc(2,wfindex(iband,indtype)+1:wfindex(iband+1,indtype)-1)*factor,&
&                   1,tdest(rvectsize+1:vectsize,iblocksize),1)
       else
         call dcopy(rvectsize,tsrc(1,wfindex(iband,indtype):wfindex(iband+1,indtype)-1)*factor,&
&                   1,tdest(1:rvectsize         ,iblocksize),1)
         call dcopy(rvectsize,tsrc(2,wfindex(iband,indtype):wfindex(iband+1,indtype)-1)*factor,&
&                   1,tdest(rvectsize+1:vectsize,iblocksize),1)
       endif
     enddo
   else
     if ( indtype == 'C' ) then
       if ( bblock ) then
         call zcopy(size,&
 &        tsrc(:,vectsize*((iblock-1)*blocksize)+lig+1:vectsize*(iblock*blocksize)+lig),&
 &        incsrc,tdest(:,1:blocksize),incdest)
       else
         call zcopy(size,&
 &        tsrc(:,lig+1:vectsize*((iblock-1)*blocksize)+lig),incsrc,tdest(:,1:bblocksize),incdest)
       endif
     else if ( indtype == 'S' ) then
       call zcopy(size,&
 &      tsrc(:,lig+1:vectsize*((iblock-1)*blocksize)+lig),&
 &      incsrc,tdest(:,1:bblocksize),incdest)
     else
       call zcopy(size,tsrc,incsrc,tdest,incdest)
     endif
   endif
 else if ( direction == 'I' ) then
   if ( x_cplx == 1 ) then
     factor=1/sqrt(2.0_dp)
     do iblocksize=1,blockiter
       iband=iblocksize
       if ( bblock ) then
         iband=iblocksize+(iblock-1)*blocksize
       endif
       if (x_me_g0 == 1) then
         !call dcopy(1,tsrc(1,iblocksize),1,tdest(1,wfindex(iband,indtype)),1)
         tdest(1,wfindex(iband,indtype))=tsrc(1,iblocksize)
         tdest(2,wfindex(iband,indtype))=0._DP
         call dcopy(rvectsize-1,tsrc(2:rvectsize         ,iblocksize)*factor,&
&                   1,tdest(1,wfindex(iband,indtype)+1:wfindex(iband+1,indtype)-1),1)
         call dcopy(rvectsize-1,tsrc(rvectsize+1:vectsize,iblocksize)*factor,&
&                   1,tdest(2,wfindex(iband,indtype)+1:wfindex(iband+1,indtype)-1),1)
       else
         call dcopy(rvectsize,tsrc(1:rvectsize         ,iblocksize)*factor,&
&                   1,tdest(1,wfindex(iband,indtype):wfindex(iband+1,indtype)-1),1)
         call dcopy(rvectsize,tsrc(rvectsize+1:vectsize,iblocksize)*factor,&
&                   1,tdest(2,wfindex(iband,indtype):wfindex(iband+1,indtype)-1),1)
       end if
     enddo
   else
     if ( indtype == 'C' ) then
       call zcopy(size,tsrc,1,tdest(:,vectsize*((iblock-1)*blocksize)+lig+1:vectsize*(iblock*blocksize)+lig),1)
     else if ( indtype == 'S' ) then
       call zcopy(size,tsrc,1,&
 &      tdest(:,vectsize*((iblock-1)*blocksize)+lig+1:vectsize*(iblock*blocksize)+lig),1)
     else
       call zcopy(size,tsrc,1,tdest,1)
     endif
   endif
 endif
 if (present(tim_wfcopy).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_wfcopy,2,tsec)
 end if
end subroutine wfcopy
!!***

end module m_wfutils
