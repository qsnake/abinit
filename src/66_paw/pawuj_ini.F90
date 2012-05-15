!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawuj_ini
!! NAME
!! pawuj_ini 
!!
!! FUNCTION
!!  Initialize dtpawuj datastructure for one SCF run
!!  Relevant only for automatic determination of U in PAW+U context
!! 
!! COPYRIGHT 
!! Copyright (C) 1998-2012 ABINIT group (DJA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtpawuj(0:ndtpawuj) (initialization of fields vsh, occ, iuj,nnat)
!!
!! PARENTS
!!      pawuj_drive,ujdet
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
subroutine pawuj_ini(dtpawuj,ndtset)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_parameters
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawuj_ini'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer                           :: ndtset
 type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtset) 

!Local variables -------------------------
!Variables for partial dos calculation
!scalars
 integer, parameter :: nwfchr=6
 integer            :: iuj,im1 
!arrays
! *********************************************************************
 
 DBG_ENTER("COLL")

!DEBUG
 write(std_out,*)'pawuj_ini enter'
!END DEBUG

 do iuj=0,ndtset
!  DEBUG
!  write(std_out,*)'pawuj_ini iuj ',iuj
!  END DEBUG
   dtpawuj(iuj)%diemix=0.45_dp
   dtpawuj(iuj)%iuj=0
   dtpawuj(iuj)%nat=0
   dtpawuj(iuj)%ndtset=1
   dtpawuj(iuj)%nspden=1
   dtpawuj(iuj)%macro_uj=0
   dtpawuj(iuj)%option=1
   dtpawuj(iuj)%pawujat=1 
   dtpawuj(iuj)%pawujga=one 
   dtpawuj(iuj)%pawprtvol=1
   dtpawuj(iuj)%ph0phiint=one
   dtpawuj(iuj)%dmatpuopt=3
   dtpawuj(iuj)%pawujrad=3.0_dp
   dtpawuj(iuj)%pawrad=20.0_dp 
!  Allocate arrays
!  DEBUG
!  write(std_out,*)'pawuj_ini before arrays'
!  END DEBUG
   ABI_ALLOCATE(dtpawuj(iuj)%rprimd,(3,3))
   ABI_ALLOCATE(dtpawuj(iuj)%scdim,(3))
   ABI_ALLOCATE(dtpawuj(iuj)%wfchr,(nwfchr))
   dtpawuj(iuj)%rprimd=reshape((/ 1,0,0,0,1,0,0,0,1/),(/ 3,3 /))
   dtpawuj(iuj)%scdim=reshape((/ 250,0,0 /),(/3 /))
   dtpawuj(iuj)%wfchr=(/ (0,im1=1,nwfchr) /)
   if (iuj>0) then 
     dtpawuj(iuj)%iuj=-1
   end if
 end do

!DEBUG
!write(std_out,*)'pawuj_ini leave'
!END DEBUG

 DBG_EXIT("COLL")
end subroutine pawuj_ini 
!!***
