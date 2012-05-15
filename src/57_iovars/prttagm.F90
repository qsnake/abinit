!{\src2tex{textfont=tt}}
!!****f* ABINIT/prttagm
!!
!! NAME
!! prttagm
!!
!! FUNCTION
!! Eventually print the content of dprarr (if typevarphys='DPR','LEN', 'ENE' and 'BFI'),
!! or intarr (if typevarphys='INT'), arrays of effective dimensions narr and 0:ndtset_alloc
!! For the second dimension, the 0 index relates to a default.
!! Print the array only if the content for at least one value of the second
!! index is different from the default.
!! Print a generic value if the non-default values are all equal.
!! Print the detail of all values otherwise.
!! The input variable 'length' controls the print format, and, in the case
!! of the real(dp) variable, the way two numbers are determined to be
!! different or not.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  intarr(1:marr,0:ndtset_alloc), dprarr(1:marr,0:ndtset_alloc)
!!   integer or real(dp) arrays, respectively,
!!   containing the data to be printed. Use these arrays even for scalars.
!!   For the first index, only the range 1:narr is relevant.
!!  iout=unit number for echoed output
!!  jdtset_(0:ndtset_alloc)=list of dataset indices.
!!  length= if 1, short format for printing, if 2, long format for printing
!!     special formats: if 3, INT : for symrel
!!                      if 4, INT : for type
!!                      if 5, INT : for mkmem, mkqmem, mk1mem
!!                      if 3, DPR : for tnons
!!                      if 4, DPR : for wtk and znucl
!!                      if 5, DPR : for atvshift
!!     If the typevarphys is 'DPR', a negative value of 'length' will request that
!!        the equality of real(dp) numbers is determined by an ABSOLUTE
!!        difference criterion only. The absolute value of length is used
!!        to determine the format, as above.
!!
!!  marr=first dimension of the intarr and dprarr arrays, as declared in the
!!   calling subroutine.
!!  narr=actual first dimension of intarr and dprarr.
!!  narrm=used when the effective first dimension of intarr is variable
!!        in this case narrm(0:ndtset_alloc)
!!  ncid= NETCDF id
!!  ndtset_alloc=govern second dimension of intarr and dprarr
!!  token=character string for 'tag'.  Assumed no longer than 9 characters
!!  typevarphys=physical variable type (might indicate the physical meaning of
!!   for dimensionality purposes)
!!   'INT'=>integer
!!   'DPR'=>real(dp) (no special treatment)
!!   'LEN'=>real(dp) (output in bohr and angstrom)
!!   'ENE'=>real(dp) (output in hartree and eV)
!!   'BFI'=>real(dp) (output in Tesla)
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      outvar1,outvars,pawuj_det
!!
!! CHILDREN
!!      appdig,leave_new,write_var_netcdf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prttagm(dprarr,intarr,iout,jdtset_,length,&
& marr,narr,narrm,ncid,ndtset_alloc,token,typevarphys,&
& use_narrm)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prttagm'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_57_iovars, except_this_one => prttagm
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,length,marr,narr,ndtset_alloc,ncid,use_narrm
 character(len=*),intent(in) :: token
 character(len=3),intent(in) :: typevarphys
!arrays
 integer,intent(in) :: intarr(marr,0:ndtset_alloc)
 integer,intent(in) :: jdtset_(0:ndtset_alloc)
 integer,intent(in) :: narrm(0:ndtset_alloc)
 real(dp),intent(in) :: dprarr(marr,0:ndtset_alloc)


!!! marr=size(intarr,DIM=1)

!Local variables-------------------------------
 character(len=*), parameter :: short_int="10i5))"
 character(len=*), parameter :: long_int ="8i8) )"
 character(len=*), parameter :: short_beg="(a,a16,a,1x,(t20,"
 character(len=*), parameter :: long_beg="(a,a16,a,1x,(t22,"
 character(len=*), parameter :: short_dpr="es16.8))"
 character(len=*), parameter :: long_dpr="es18.10))"
 character(len=*), parameter :: short_dim="es16.8),a)"
 character(len=*), parameter :: long_dim="es18.10),a)"
 character(len=*), parameter :: f_symrel="3(3i3,1x),4x,3(3i3,1x)))"
 character(len=*), parameter :: f_type  ="20i3))"
 character(len=*), parameter :: f_mem   ="8i8))"
 character(len=*), parameter :: f_tnons ="3f11.7,3x,3f11.7))"
 character(len=*), parameter :: f_wtk   ="6f11.5))"
 character(len=*), parameter :: f_atvshift   ="5f11.5))"
!scalars
 integer :: iarr,idtset,jdtset,multi,ndtset_eff,print,narr_eff
 real(dp),parameter :: tol21=1.0d-21
 real(dp) :: diff,scale_factor,sumtol
 character(len=1) :: digit,first_column
 character(len=4) :: appen
 character(len=8) :: out_unit
 character(len=50) :: format_dp,format_int,full_format
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*)'prttagm start'
!END DEBUG

!###########################################################
!### 01. Check consistency of input

 if(len_trim(token)>16)then
   write(message, '(a,a,a,a,i6,a,a,a,a,a)' ) ch10,&
&   ' prttagm : BUG -',ch10,&
&   '  The length of the name of the input variable ',trim(token),' is ',len_trim(token),ch10,&
&   '  This exceeds 16 characters, the present maximum in routine prttagm.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(ndtset_alloc<1)then
   write(message, '(a,a,a,a,i6,a,a,a,a,a)' ) ch10,&
&   ' prttagm : BUG -',ch10,&
&   '  ndtset_alloc=',ndtset_alloc,', while it should be >= 1.',ch10,&
&   '  This happened for token=',token,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(ndtset_alloc>9999)then
   write(message, '(a,a,a,a,i6,a,a,a,a,a)' ) ch10,&
&   ' prttagm : BUG -',ch10,&
&   '  ndtset_alloc=',ndtset_alloc,', while it must be lower than 10000.',ch10,&
&   '  This happened for token=',token,'.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if(narr>99 .and. (typevarphys=='ENE'.or.typevarphys=='LEN'))then
   write(message, '(6a,i6,a)' ) ch10,&
&   ' prttagm : BUG -',ch10,&
&   '  typevarphys=',typevarphys,' with narr=',narr,'  is not allowed.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 if ((narr>0).or.(use_narrm/=0)) then

   print=1
   multi=0

!  ###########################################################
!  ### 02. Treatment of integer 'INT'

   if(typevarphys=='INT')then

!    Determine whether the different non-default occurences are all equal

     if (use_narrm==0) then ! use of scalar 'narr' instead of array 'narrm'
       if(ndtset_alloc>1)then
         do idtset=1,ndtset_alloc
           do iarr=1,narr
             if(intarr(iarr,1)/=intarr(iarr,idtset))multi=1
           end do
         end do
       end if
     else
!      If the sizes of the arrays are different we can not compare them
!      So we have to asume they are different
       multi=1
     end if

!    If they are all equal, then determine whether they are equal to the default
     if(multi==0)then
       print=0
       do iarr=1,narr
         if(intarr(iarr,1)/=intarr(iarr,0))print=1
       end do
     end if

!    Print only if the values differ from the default
     if((print==1).or.(ncid<0))then
       ndtset_eff=ndtset_alloc
       if((multi==0).or.(ncid<0)) ndtset_eff=1
       do idtset=1,ndtset_eff
!        Initialize the character in the first column
         first_column=' '
         if(abs(length)==5)first_column='P'
!        Initialize the format
         if(abs(length)==1)format_int=trim(short_int)
         if(abs(length)==2)format_int=trim(long_int)
         if(abs(length)==3)format_int=trim(f_symrel)
         if(abs(length)==4)format_int=trim(f_type)
         if(abs(length)==5)format_int=trim(f_mem)
!        Initialize the dataset number string, and print
         if((multi==0).or.(ncid<0))then
           appen=' '
         else
           jdtset=jdtset_(idtset)
           call appdig(jdtset,'',appen)
         end if
         full_format=trim(long_beg)//trim(format_int)

!        narr_eff could be narr or narrm(idtset)
!        It depends if the size is variable for different datasets
         if (use_narrm==0)then
           narr_eff=narr
         else
           narr_eff=narrm(idtset)
         end if

         if (narr_eff/=0) then

           if (print==1) write(iout,full_format)first_column,token,trim(appen),&
&           intarr(1:narr_eff,idtset)
#if defined HAVE_TRIO_ETSF_IO
           call write_var_netcdf(intarr(1:narr_eff,idtset),&
&           dprarr(1:narr_eff,idtset),&
&           marr,narr_eff,abs(ncid),typevarphys,token//appen)
#endif
         end if

       end do
     end if !(print==1)

!    ###########################################################
!    ### 03. Treatment of real 'DPR', 'LEN', 'ENE', 'BFI'

   else if (typevarphys=='DPR' .or. typevarphys=='LEN' .or. typevarphys=='ENE' .or. typevarphys=='BFI') then

     if((ndtset_alloc>1).and.(use_narrm==0))then
       do idtset=1,ndtset_alloc
         do iarr=1,narr
!          The determination of effective equality is more difficult than in the
!          integer case :
!          - if length > 0, ask for a relative accuracy, and also include
!          the case of zero values, thanks to tol21.
!          - if length < 0, ask for absolute accuracy.
           diff=abs( dprarr(iarr,1)-dprarr(iarr,idtset) )
           if(length>0)then
             sumtol=abs(dprarr(iarr,1))+abs(dprarr(iarr,idtset))+10*tol21
             if(diff>sumtol*tol11)multi=1
           else
             if(diff>tol14)multi=1
           end if
         end do
       end do
     elseif (use_narrm/=0) then
       multi=1 ! Assume that values could not be compare with a reference
     end if

     if(multi==0)then
       print=0
       do iarr=1,narr
         diff=abs( dprarr(iarr,1)-dprarr(iarr,0) )
         if(length>0)then
           sumtol=abs(dprarr(iarr,1))+abs(dprarr(iarr,0))+10*tol21
           if(diff>sumtol*tol11)print=1
         else
           if(diff>tol14)print=1
         end if
       end do
     end if

     if((print==1).or.(ncid<0))then
!      Select the proper format
       if(abs(length)==1 .or. abs(length)==2)then
         if(typevarphys=='DPR')then
           digit='3'
           if(abs(length)==1)format_dp=digit//short_dpr
           if(abs(length)==2)format_dp=digit//long_dpr
         else if(typevarphys=='ENE' .or. typevarphys=='LEN' .or. typevarphys=='BFI')then
           if (narr<10) write(digit,'(i1)')narr
           if (narr> 9) write(digit,'(i2)')narr
           if(abs(length)==1)format_dp=digit//short_dim
           if(abs(length)==2)format_dp=digit//long_dim
         end if
       else
         if(abs(length)==3)format_dp=f_tnons
         if(abs(length)==4)format_dp=f_wtk
         if(abs(length)==5)format_dp=f_atvshift
       end if
       ndtset_eff=ndtset_alloc
       if((multi==0).or.(ncid<0))ndtset_eff=1
       do idtset=1,ndtset_eff

!        narr_eff could be narr or narrm(idtset)
!        It depends if the size is variable for different datasets
         if (use_narrm==0)then
           narr_eff=narr
         else
           narr_eff=narrm(idtset)
         end if

         if (narr_eff/=0) then

!          Initialize the character in the first column
           first_column=' '
!          Define scale_factor
           scale_factor=one
           if(typevarphys=='BFI')scale_factor=one/BField_Tesla
!          Define out_unit
           if(typevarphys=='ENE')out_unit=' Hartree'
           if(typevarphys=='LEN')out_unit=' Bohr   '
           if(typevarphys=='BFI')out_unit=' Tesla  '
!          Format, according to the length of the dataset string
           if((multi==0).or.(ncid<0))then
             appen=' '
           else
             jdtset=jdtset_(idtset)
             call appdig(jdtset,'',appen)
           end if
           full_format=trim(long_beg)//trim(format_dp)
           if(typevarphys=='DPR')then
             if (print==1) write(iout,full_format)&
&             first_column,token,trim(appen),dprarr(1:narr_eff,idtset)*scale_factor
           else
             if (print==1) write(iout,full_format)&
&             first_column,token,trim(appen),dprarr(1:narr_eff,idtset)*scale_factor,trim(out_unit)
           end if
#if defined HAVE_TRIO_ETSF_IO
           call write_var_netcdf(intarr(1:narr_eff,idtset),dprarr(1:narr_eff,idtset),&
&           marr,narr_eff,abs(ncid),'DPR',token//trim(appen))
#endif

         end if

       end do
     end if

!    ###########################################################
!    ### 04. The type is neither 'INT' nor 'DPR','ENE','LEN','BFI'
   else

     write(message, '(a,a,a,a,a,a)' ) ch10,&
&     ' prttagm : BUG -',ch10,&
&     '  Disallowed typevarphys=',typevarphys,'.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')

   end if

!  End condition of narr>0
 end if

end subroutine prttagm
!!***
