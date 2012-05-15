!{\src2tex{textfont=tt}}
!!****f* ABINIT/read_wfrspa
!! NAME
!! read_wfrspa
!!
!! FUNCTION
!!  Internal subroutine used to fill wfrspa from the scratch file when we
!!  have mkmem==0.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  fnametmp_tdwf=name of the file appended with _TDWF
!!  state=the current state is an id to know if this band has already been read or not
!!  (to be completed)
!!
!! OUTPUT
!!  iband=index of the band after reading
!!  imkmem=
!!  isppol=index of the polarisation after reading
!!  (to be completed)
!!
!! SIDE EFFECTS
!!  eigbnd=the eigenvalue of this band
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  wfrspa=read values
!!
!! PARENTS
!!      tddft
!!
!! CHILDREN
!!      int2char4
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine read_wfrspa(state,fnametmp_tdwf,eigbnd,iband,isppol,imkmem,ndiel4,ndiel5,ndiel6,wfrspa_extract)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_wfrspa'
 use interfaces_27_toolbox_oop
!End of the abilint section

 implicit none

!Arguments -----------------------------
!scalars
 integer,intent(out) :: iband,isppol
 integer,intent(in) :: imkmem,ndiel4,ndiel5,ndiel6,state
 real(dp),intent(inout) :: eigbnd
 character(len=fnlen),intent(in) :: fnametmp_tdwf
!arrays
 real(dp),intent(inout) :: wfrspa_extract(ndiel4,ndiel5,ndiel6)

!Local variables---------------------------------
!scalars
 character(len=4) :: tag
 character(len=fnlen) :: fnametmp_tdwf_tag
!arrays
 integer,save :: states(4)=(/0,0,0,0/)
!real(dp),allocatable :: work2(:,:,:)

! ***********************************************************************

 if (state == states(imkmem)) then
!  We have already read wfrspa(:,:,:,imkmem) last time
   return
 else
!  We save the band number for future reference and continue reading
   states(imkmem)=state
 end if

 call int2char4(state,tag)
 fnametmp_tdwf_tag=trim(fnametmp_tdwf)//tag

!DEBUG
!write(message,'(a,a,a,i3)')'reading phi : ',fnametmp_tdwf_tag, &
!& ' by proc ',me_loc
!call wrtout(std_out,message,'PERS')
!ENDDEBUG

 open(tmp_unit,file=fnametmp_tdwf_tag,form='unformatted',status='unknown')
 read(tmp_unit)iband,isppol,eigbnd
 read(tmp_unit)wfrspa_extract
 close(tmp_unit)

end subroutine read_wfrspa
!!***
