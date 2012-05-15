!{\src2tex{textfont=tt}}
!!****f* ABINIT/testfi
!!
!! NAME
!! testfi
!!
!! FUNCTION
!! Routine "Final test" for generation of the test report in the status file :
!! if it appears that the run was a "Build-in Test", then
!! compare the final values of different quantities to the reference
!! values, here hard-coded.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA,XG,GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  etotal=total energy (sum of 7 contributions) (hartree)
!!  filnam(5)=names of the files
!!  filstat=name of the status file
!!  fred(3,natom)=symmetrized gradient of etotal with respect to tn
!!  natom=number of atoms in cell.
!!  strten(6)=components of the stress tensor (hartree/bohr^3)
!!  xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine testfi(etotal,filnam,filstat,fred,natom,strten,xred)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'testfi'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 real(dp),intent(in) :: etotal
 character(len=fnlen),intent(in) :: filstat
!arrays
 real(dp),intent(in) :: fred(3,natom),strten(6),xred(3,natom)
 character(len=fnlen),intent(in) :: filnam(5)

!Local variables-------------------------------
 !character(len=9), parameter :: chara = '123456789'
 character(len=fnlen) :: testname(7)='         '
 character(len=*), parameter :: format01000="(a,d22.14,a,d12.4)"
!scalars
 integer,parameter :: mtest=7
 integer :: iatom,ii,testin,tok
 real(dp) :: etot_mxdev,etot_ref,fred_mxdev,strten_mxdev,xred_mxdev
 character(len=fnlen) :: nam
!arrays
 integer,parameter :: natom_test(mtest)=(/2,1,2,1,1,1,2/)
 real(dp) :: fred_ref(3,2),strten_ref(6),xred_ref(3,2)

! ***********************************************************************

!DEBUG
 write(std_out,*)' testfi : enter '
 write(std_out,*)'filnam(1)="',trim(filnam(1)),'"'
!ENDDEBUG
 testname(1)='in_1'
 testname(2)='in_2'
 testname(3)='in_3'
 testname(4)='in_bigdft'
 testname(5)='in_etsf_io'
 testname(6)='in_libxc'
 testname(7)='in_wannier90'

!Determine whether the run was a build-in test
 testin=0
 do ii=1,mtest
   nam='test'//trim(testname(ii))//'.in   '
!  DEBUG
   write(std_out,*)'nam="',trim(nam),'"'
!  ENDDEBUG
   if(trim(filnam(1))/=trim(nam)    )cycle
   nam='test'//trim(testname(ii))//'.out  '
   if(trim(filnam(2))/=trim(nam)    )cycle
   nam='test'//trim(testname(ii))//'i     '
   if(trim(filnam(3))/=trim(nam)    )cycle
   nam='test'//trim(testname(ii))//'o     '
   if(trim(filnam(4))/=trim(nam)    )cycle
   nam='test'//trim(testname(ii))//'      '
   if(trim(filnam(5))/=trim(nam)    )cycle
   if( natom /= natom_test(ii) ) cycle
   testin=ii
   exit
 end do

!DEBUG
 write(std_out,*)' testfi : testin= ',testin
!ENDDEBUG

!---------------------------------------------------------

!Now, open the status file, and either delete it, or produce a report
 open (tmp_unit,file=trim(filstat),form='formatted',status='unknown')

 if(testin==0)then

   close (tmp_unit,status='delete')

 else

!  Note : all processors have their own file, so no special
!  attention must be paid to the parallel case.

   write(tmp_unit,*)
   write(tmp_unit,*)'Status file, reporting on test ',trim(testname(ii))
   write(tmp_unit,*)

!  Define reference values, as well as maximum tolerable deviation
   if(testin==1)then

     etot_ref=-1.05814441948188d+00
     etot_mxdev=1.0d-9
     xred_ref(1:3,1)=(/ -0.65048430042634D-01 , 0.0_dp , 0.0_dp /)
     xred_ref(1:3,2)=(/  0.65048430042634D-01 , 0.0_dp , 0.0_dp /)
!    xred(*,*) are reduced coordinates
     xred_mxdev=1.0d-6
     fred_ref(1:3,1:2)= 0.0_dp
!    fred(*,*) are gradients with respect to reduced coordinates
     fred_mxdev=5.0d-4
     strten_ref(1:6)=(/ 0.149D-04  , 0.560D-04 , 0.560D-04 ,&
&     0.0_dp , 0.0_dp , 0.0_dp /)
     strten_mxdev=1.0d-5

   else if(testin==2)then
     
!    This value of etot is accurate to about 1.0d-12
     etot_ref=-.70811958266295D+02 
!    Initialisation conditions might give fluctuations
     etot_mxdev=3.0d-7
     xred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-12
!    This value of strten is accurate to at least 1.0d-8
     strten_ref(1:3)= 5.09324870E-03
     strten_ref(4:6)= 0.0_dp
     strten_mxdev=1.0d-6

   else if(testin==3)then

     etot_ref=-.892746696311772D+01
     etot_mxdev=1.0d-8
     xred_ref(1:3,1)=(/ -0.125_dp , 0.0_dp , 0.0_dp /)
     xred_ref(1:3,2)=(/  0.125_dp , 0.0_dp , 0.0_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ -.140263620278D+00 , 0.0_dp , 0.0_dp /)
     fred_ref(1:3,2)=(/  .140013483725D+00 , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-3
     strten_ref(1:6)=(/  1.3949d-3 ,  1.3643d-3 , 1.3643d-3 ,&
&     .0_dp ,  .0_dp  ,  .0_dp    /)
     strten_mxdev=1.0d-5

!    Bigdft
   else if(testin==4)then

!    This value of etot is accurate to about 1.0d-12
     etot_ref=-0.56231990141D+00
!    Initialisation conditions might give fluctuations
     etot_mxdev=3.0d-7
     xred_ref(1:3,1)=(/ 0.5_dp , 0.5_dp , 0.5_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-12
     strten_ref(1)=-0.22537238382D-01
     strten_ref(2)=-0.22536232141D-01
     strten_ref(3)=-0.22529038043D-01
     strten_ref(4)= 0.39104928264D-06
     strten_ref(5)= 0.62823846408D-06
     strten_ref(6)= 0.38414125548E-09
     strten_mxdev=1.0d-7

!    ETSF_IO
   else if(testin==5)then

!    This value of etot is accurate to about 1.0d-12
     etot_ref=-0.33307825915795D+02
!    Initialisation conditions might give fluctuations
     etot_mxdev=3.0d-7
     xred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-7
!    This value of strten is accurate to at least 1.0d-8
     strten_ref(1)= 0.13569940015175D-01
     strten_ref(2)= 0.33822108610352D-01
     strten_ref(3)= 0.37991262607028D-01
     strten_ref(4:6)= 0.0_dp
     strten_mxdev=1.0d-6

!    LibXC
   else if(testin==6)then

     etot_ref=-0.55196688523897D+01
     etot_mxdev=5.0d-7
     xred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-12
     strten_ref(1:3)= 0.13246699375127D-04
     strten_ref(4:6)= 0.0_dp
     strten_mxdev=1.0d-8

!    Wannier90
   else if(testin==7)then

!    This value of etot is accurate to about 1.0d-12
     etot_ref=-0.10620085133544D+02
!    Initialisation conditions might give fluctuations
     etot_mxdev=3.0d-7
     xred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     xred_ref(1:3,2)=(/ 0.25_dp , 0.25_dp , 0.25_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_ref(1:3,2)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-12
!    This value of strten is accurate to at least 1.0d-8
     strten_ref(1:3)=0.31413922197317D-03
     strten_ref(4:6)= 0.0_dp
     strten_mxdev=1.0d-6

!    End of the reference value set up, for different tests
   end if

!  Compare reference and actual values

!  Take care of total energy
   if(abs(etot_ref-etotal)<etot_mxdev)then
     write(tmp_unit,'(a)')' OK for total energy '
   else
     write(tmp_unit,'(a)')' Error for total energy : '
     write(tmp_unit,format01000)'        expected ',etot_ref,&
&     '  with maximum   deviation',etot_mxdev
     write(tmp_unit,format01000)'        computed ',etotal,&
&     '  with effective deviation',abs(etotal-etot_ref)
   end if

!  Take care of nuclei positions
   tok=1
   do iatom=1,natom
     do ii=1,3
       if(abs(xred_ref(ii,iatom)-xred(ii,iatom))>xred_mxdev)then
         tok=0
         write(tmp_unit, '(a,i1,a,i1,a)' )&
&         ' Error for nuclei position xred(',ii,',',iatom,')'
         write(tmp_unit,format01000)'        expected ',xred_ref(ii,iatom),&
&         '  with maximum   deviation',xred_mxdev
         write(tmp_unit,format01000)'        computed ',xred(ii,iatom),&
&         '  with effective deviation',&
&         abs( xred(ii,iatom)-xred_ref(ii,iatom) )
       end if
     end do
   end do
   if(tok==1)write(tmp_unit,'(a)')' OK for nuclei positions '

!  Take care of forces
   tok=1
   do iatom=1,natom
     do ii=1,3
       if(abs(fred_ref(ii,iatom)-fred(ii,iatom))&
&       >fred_mxdev)then
         tok=0
         write(tmp_unit, '(a,i1,a,i1,a)' )&
&         ' Error for force fred(',ii,',',iatom,')'
         write(tmp_unit,format01000)'        expected ',fred_ref(ii,iatom),&
&         '  with maximum   deviation',fred_mxdev
         write(tmp_unit,format01000)'        computed ',fred(ii,iatom),&
&         '  with effective deviation',&
&         abs( fred(ii,iatom)-fred_ref(ii,iatom) )
       end if
     end do
   end do
   if(tok==1)write(tmp_unit,'(a)')' OK for forces '

!  Take care of stress
   tok=1
   do ii=1,6
     if(abs(strten_ref(ii)-strten(ii))>strten_mxdev)then
       tok=0
       write(tmp_unit,'(a,i1,a)')' Error for stress strten(',ii,')'
       write(tmp_unit,format01000)'        expected ',strten_ref(ii),&
&       '  with maximum   deviation',strten_mxdev
       write(tmp_unit,format01000)'        computed ',strten(ii),&
&       '  with effective deviation',&
&       abs( strten(ii)-strten_ref(ii) )
     end if
   end do
   if(tok==1)write(tmp_unit,'(a)')' OK for stresses '

   write(tmp_unit,*)

   close(tmp_unit)

!  End of the choice between produce a report, and produce no report
 end if

!DEBUG
!write(std_out,*)' testfi : exit '
!ENDDEBUG

end subroutine testfi

!!***
