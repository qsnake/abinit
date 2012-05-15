!{\src2tex{textfont=tt}}
!!****f* ABINIT/read_md_hist
!!
!! NAME
!! read_md_hist
!!
!! FUNCTION
!! Read the history file in netcdf format and store it
!! into a hist dataset structure
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filename = Filename of the NetCDF to read
!!
!! OUTPUT
!! hist<type ab_movehistory>=Historical record of positions, forces
!!      |                    acell, stresses, and energies,
!!      |                    contains:
!!      | mxhist:  Maximun number of records
!!      | ihist:   Index of present record of hist
!!      | histA:   Historical record of acell(A) and rprimd(R)
!!      | histE:   Historical record of energy(E)
!!      | histR:   Historical record of rprimd(R)
!!      | histS:   Historical record of strten(S)
!!      | histV:   Historical record of velocity(V)
!!      | histXF:  Historical record of positions(X) and forces(F)
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      handle_ncerr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine read_md_hist(filename,hist)

 use m_profiling

! define dp,sixth,third,etc...
 use defs_basis
! type(ab_movetype), type(ab_movehistory)
 use defs_mover
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_md_hist'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 type(ab_movehistory),intent(out) :: hist
 character(len=fnlen),intent(in) :: filename
!Local variables-------------------------------
!scalars
 integer :: ncerr,ncid
 integer :: natom,time
 integer :: xyz_id,natom_id,time_id,six_id
 integer :: xcart_id,xred_id,fcart_id,fred_id,ekin_id,mdtime_id
 integer :: vel_id,etotal_id,acell_id,rprimd_id,strten_id
 character(len=5) :: char_tmp
!arrays

! *************************************************************************

!fname_hist=trim(ab_mover%filnam_ds(4))//'_HIST'

#if defined HAVE_TRIO_NETCDF

!#####################################################################
!### Reading of NetCDF file

!1. Open netCDF file

 ncerr=nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ncid)

 if(ncerr /= NF90_NOERR) then
   write(std_out,*) 'Could no open ',trim(filename),', starting from scratch'
   hist%ihist=0
   hist%mxhist=0
   return
 else
   write(std_out,*) 'Succesfully open ',trim(filename),' for reading'
   write(std_out,*) 'Extracting information from NetCDF file...'
   hist%ihist=0
   hist%mxhist=0
 end if

!2. Inquire dimensions IDs

 ncerr = nf90_inq_dimid(ncid,"natom",natom_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," inquire dimension ID for natom")
 ncerr = nf90_inq_dimid(ncid,"xyz",xyz_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," inquire dimension ID for xyz")
 ncerr = nf90_inq_dimid(ncid,"time",time_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," inquire dimension ID for time")
 ncerr = nf90_inq_dimid(ncid,"six",six_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," inquire dimension ID for six")

!3. Inquire dimensions lenghts

 ncerr = nf90_inquire_dimension(ncid,natom_id,char_tmp,natom)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," inquire dimension ID for natom")
 ncerr = nf90_inquire_dimension(ncid,time_id,char_tmp,time)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," inquire dimension ID for time")

 write(std_out,*) 'Number of atoms readed:',natom
 write(std_out,*) 'Number of iterations recorded:',time

!4. Allocate hist structure 

 hist%ihist=1
 hist%mxhist=time
 ABI_ALLOCATE(hist%histA,(3,hist%mxhist))
 ABI_ALLOCATE(hist%histE,(hist%mxhist))
 ABI_ALLOCATE(hist%histEk,(hist%mxhist))
 ABI_ALLOCATE(hist%histT,(hist%mxhist))
 ABI_ALLOCATE(hist%histR,(3,3,hist%mxhist))
 ABI_ALLOCATE(hist%histS,(6,hist%mxhist))
 ABI_ALLOCATE(hist%histV,(3,natom,hist%mxhist))
 ABI_ALLOCATE(hist%histXF,(3,natom,4,hist%mxhist))

!5. Get the ID of a variables from their name

 ncerr = nf90_inq_varid(ncid, "xcart", xcart_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for xcart")
 ncerr = nf90_inq_varid(ncid, "xred", xred_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for xred")
 ncerr = nf90_inq_varid(ncid, "fcart", fcart_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for fcart")
 ncerr = nf90_inq_varid(ncid, "fred", fred_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for fred")
 ncerr = nf90_inq_varid(ncid, "vel", vel_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for vel")
 ncerr = nf90_inq_varid(ncid, "acell", acell_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for acell")
 ncerr = nf90_inq_varid(ncid, "rprimd", rprimd_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for rprimd")
 ncerr = nf90_inq_varid(ncid, "etotal", etotal_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for etotal")
 ncerr = nf90_inq_varid(ncid, "ekin", ekin_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for ekin")
 ncerr = nf90_inq_varid(ncid, "mdtime", mdtime_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for mdtime")
 ncerr = nf90_inq_varid(ncid, "strten", strten_id)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," get the id for strten")

!#####################################################################
!### Read variables from the dataset and write them into hist

!XCART,XRED,FCART,FRED,VEL
 ncerr = nf90_get_var(ncid, xcart_id,hist%histXF(:,:,1,:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable xcart")
 ncerr = nf90_get_var(ncid, xred_id,hist%histXF(:,:,2,:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable xred")
 ncerr = nf90_get_var(ncid, fcart_id,hist%histXF(:,:,3,:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable fcart")
 ncerr = nf90_get_var(ncid, fred_id,hist%histXF(:,:,4,:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable fred")
 ncerr = nf90_get_var(ncid, vel_id,hist%histV(:,:,:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable vel")

!RPRIMD
 ncerr = nf90_get_var(ncid, rprimd_id,hist%histR(:,:,:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable rprimd")

!ACELL
 ncerr = nf90_get_var(ncid, acell_id,hist%histA(:,:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable acell")

!STRTEN
 ncerr = nf90_get_var(ncid, strten_id,hist%histS(:,:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable strten")

!ETOTAL
 ncerr = nf90_get_var(ncid, etotal_id,hist%histE(:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable etotal")

!Ekin
 ncerr = nf90_get_var(ncid, ekin_id,hist%histEk(:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable ekin")

!MDTime
 ncerr = nf90_get_var(ncid, mdtime_id,hist%histT(:))
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," read variable mdtime")

!#####################################################################
!### Close NetCDF file
 ncerr = nf90_close(ncid)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," close netcdf history file")

#endif

end subroutine read_md_hist
!!***
