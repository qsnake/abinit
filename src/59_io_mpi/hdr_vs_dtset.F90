!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_vs_dtset
!! NAME
!! hdr_vs_dtset
!!
!! FUNCTION
!!  Check the compatibility of the Abinit header with respect to the
!!  input variables defined in the input file. Mainly used in case 
!!  of GW calculation
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Dtset<type(dataset_type)>=all input variables for this dataset
!!  Hdr <type(hdr_type)>=the header structured variable
!!
!! OUTPUT
!!  Only check
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      rdm,setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine hdr_vs_dtset(Hdr,Dtset)

 use m_profiling
    
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_crystal, only : print_symmetries

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_vs_dtset'
 use interfaces_14_hidewrite
 use interfaces_59_io_mpi, except_this_one => hdr_vs_dtset
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(Hdr_type),intent(in) :: Hdr
 type(Dataset_type),intent(inout) :: Dtset

!Local variables-------------------------------
 logical,parameter :: NO_FATAL=.FALSE.,IS_FATAL=.TRUE.
 integer :: ik,jj,ierr
 logical :: test
 logical :: tsymrel,ttnons,tsymafm
 character(len=500) :: msg      
! *************************************************************************

!=== Check basic dimensions ===
 ierr=0
 call compare_and_change_int(Hdr%natom,  Dtset%natom,  'natom',  IS_FATAL,ierr)
 call compare_and_change_int(Hdr%nkpt,   Dtset%nkpt,   'nkpt',   IS_FATAL,ierr)
 call compare_and_change_int(Hdr%npsp,   Dtset%npsp,   'npsp',   IS_FATAL,ierr)
 call compare_and_change_int(Hdr%nspden, Dtset%nspden, 'nspden', IS_FATAL,ierr)
 call compare_and_change_int(Hdr%nspinor,Dtset%nspinor,'nspinor',IS_FATAL,ierr)
 call compare_and_change_int(Hdr%nsppol, Dtset%nsppol, 'nsppol', IS_FATAL,ierr)
 call compare_and_change_int(Hdr%nsym,   Dtset%nsym,   'nsym',   IS_FATAL,ierr)
 call compare_and_change_int(Hdr%ntypat, Dtset%ntypat, 'ntypat', IS_FATAL,ierr)
 call compare_and_change_int(Hdr%usepaw, Dtset%usepaw, 'usepaw', IS_FATAL,ierr)
 call compare_and_change_int(Hdr%usewvl, Dtset%usewvl, 'usewvl', IS_FATAL,ierr)

!=== The number of fatal errors must be zero ===
 if (ierr/=0) then 
   write(msg,'(3a)')&
&   ' Cannot continue, basic dimensions reported in the header do not agree with input file. ',ch10,&
&   ' Check consistency between the content of the external file and the input file. '
   MSG_ERROR(msg)
 end if

 test=ALL(ABS(Hdr%xred-Dtset%xred_orig(:,1:Dtset%natom,1))<tol6)
 ABI_CHECK(test,'Mismatch in xred')

 test=ALL(Hdr%typat==Dtset%typat(1:Dtset%natom)) 
 ABI_CHECK(test,'Mismatch in typat')
!
!* Check if the lattice from the input file agrees with that read from the KSS file
 if ( (ANY(ABS(Hdr%rprimd-Dtset%rprimd_orig(1:3,1:3,1))>tol6)) ) then
   write(msg,'(6a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   ' real lattice vectors read from Header ',ch10,&
&   ' differ from the values specified in the input file'
   call wrtout(std_out,msg,'COLL')
   write(msg,'(3a,3(3es16.6),3a,3(3es16.6),3a)')ch10,&
&   ' rprimd from Hdr file   = ',ch10,(Hdr%rprimd(:,jj),jj=1,3),ch10,&
&   ' rprimd from input file = ',ch10,(Dtset%rprimd_orig(:,jj,1),jj=1,3),ch10,ch10,&
&   '  Modify the lattice vectors in the input file '
   call wrtout(std_out,msg,'COLL') 
   MSG_ERROR("")
 end if 

!=== Check symmetry operations ===
 tsymrel=(ALL(Hdr%symrel==Dtset%symrel(:,:,1:Dtset%nsym)))
 if (.not.tsymrel) then
   write(msg,'(6a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   ' real space symmetries read from Header ',ch10,&
&   ' differ from the values inferred from the input file'
   call wrtout(std_out,msg,'COLL')
   tsymrel=.FALSE.
 end if 

 ttnons=ALL(ABS(Hdr%tnons-Dtset%tnons(:,1:Dtset%nsym))<tol6)
 if (.not.ttnons) then
   write(msg,'(6a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   ' fractional translations read from Header ',ch10,&
&   ' differ from the values inferred from the input file'
   call wrtout(std_out,msg,'COLL')
   ttnons=.FALSE.
 end if 

 tsymafm=ALL(Hdr%symafm==Dtset%symafm(1:Dtset%nsym))
 if (.not.tsymafm) then
   write(msg,'(6a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   ' AFM symmetries read from Header ',ch10,&
&   ' differ from the values inferred from the input file'
   call wrtout(std_out,msg,'COLL')
   tsymafm=.FALSE.
 end if

 if (.not.(tsymrel.and.ttnons.and.tsymafm)) then
   write(msg,'(a)')' Header ' 
   call wrtout(std_out,msg,'COLL') 
   call print_symmetries(Hdr%nsym,Hdr%symrel,Hdr%tnons,Hdr%symafm)
   write(msg,'(a)')' Dtset  ' 
   call wrtout(std_out,msg,'COLL') 
   call print_symmetries(Dtset%nsym,Dtset%symrel,Dtset%tnons,Dtset%symafm)
   MSG_ERROR('Check symmetry operations ')
 end if

!* Check if the k-points from the input file agrees with that read from the KSS file
 if ( (ANY(ABS(Hdr%kptns(:,:)-Dtset%kpt(:,1:Dtset%nkpt))>tol6)) ) then
   write(msg,'(9a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   '  k-points read from Header ',ch10,&
&   '  differ from the values specified in the input file',ch10,&
&   '  k-points from Hdr file                        | k-points from input file ',ch10
   call wrtout(std_out,msg,'COLL') 
   do ik=1,Dtset%nkpt
     write(msg,'(3(3es16.6,3x))')Hdr%kptns(:,ik),Dtset%kpt(:,ik)
     call wrtout(std_out,msg,'COLL') 
   end do
   write(msg,'(a)')' Modify the k-mesh in the input file '
   MSG_ERROR(msg)
 end if 

 if (ANY(ABS(Hdr%wtk(:)-Dtset%wtk(1:Dtset%nkpt))>tol6)) then
   write(msg,'(9a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   '  k-point weights read from Header ',ch10,&
&   '  differ from the values specified in the input file',ch10,&
&   '  Hdr file  |  File ',ch10
   call wrtout(std_out,msg,'COLL') 
   do ik=1,Dtset%nkpt
     write(msg,'(2(f11.5,1x))')Hdr%wtk(ik),Dtset%wtk(ik)
     call wrtout(std_out,msg,'COLL') 
   end do
   write(msg,'(a)')' Check the k-mesh and the symmetries of the system. '
   MSG_ERROR(msg)
 end if 

!Check istwfk storage
 if ( (ANY(Hdr%istwfk(:)/=Dtset%istwfk(1:Dtset%nkpt))) ) then
   write(msg,'(9a)')ch10,&
&   ' hdr_vs_dtset : ERROR - ',ch10,&
&   '  istwfk read from Header ',ch10,&
&   '  differ from the values specified in the input file',ch10,&
&   '  Hdr | input ',ch10
   call wrtout(std_out,msg,'COLL') 
   do ik=1,Dtset%nkpt
     write(msg,'(i5,3x,i5)')Hdr%istwfk(ik),Dtset%istwfk(ik)
     call wrtout(std_out,msg,'COLL') 
   end do
   write(msg,'(a)')' Modify istwfk in the input file '
   MSG_ERROR(msg)
 end if 

 CONTAINS  !===========================================================
!!***

!!****f* hdr_vs_dtset/compare_and_change_int
!! NAME
!! compare_and_change_int
!!
!! FUNCTION
!!  Compare to int value and may raise an exception on error.
!!
!! INPUTS
!!  int_exp= expected value.
!!  info_string= some information to put in the error if any.
!!  isfatal= if .true., increase ierr and left int_found unchanged. 
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ierr=if isfatal is .true. and int_exp is different from int_found,
!!       then increase ierr of one.
!!  int_found=the read value, may be changed to int_exp if isfatal is .false..
!!
!! PARENTS
!!      hdr_vs_dtset
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

 subroutine compare_and_change_int(int_exp,int_found,info_string,isfatal,ierr)

 use m_profiling

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compare_and_change_int'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: int_exp
 integer,intent(inout) :: int_found,ierr
 character(len=*),intent(in) :: info_string
 logical,intent(in) :: isfatal

!Local variables-------------------------------
 logical :: leq                       
 character(len=500) :: msg                                              
! *************************************************************************

   leq=(int_exp==int_found)

   if (.not.leq) then 
     write(msg,'(4a,i6,a,i6)')ch10,&
&     ' hdr_vs_dtset : WARNING - Mismatch in '//TRIM(info_string),ch10,&
&     '  Expected = ',int_exp,' Found = ',int_found
     call wrtout(std_out,msg,'COLL') 
     if (isfatal) then
!      Increase ierr to signal we should stop in the caller.
       ierr=ierr+1 
     else 
!      Change the value found to reflect the expected one.
       int_found=int_exp  
       write(msg,'(a)')' The value found has been set to the expected one.'
       call wrtout(std_out,msg,'COLL') 
     end if
   end if

 end subroutine compare_and_change_int
!!***

end subroutine hdr_vs_dtset
