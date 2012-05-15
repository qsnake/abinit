!{\src2tex{textfont=tt}}
!!****f* ABINIT/write_md_hist
!!
!! NAME
!! write_md_hist
!!
!! FUNCTION
!! Write the history into a netcdf dataset
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
!! filname = Filename of the file where the history will be stored
!! natom = Number of atoms
!! hist<type ab_movehistory>=Historical record of positions, forces
!!      |                    acell, stresses, and energies,
!!      |                    contains:
!!      | mxhist:  Maximun number of records
!!      | histA:   Historical record of acell(A) and rprimd(R)
!!      | histE:   Historical record of energy(E)
!!      | histEk:  Historical record of Ionic kinetic energy(Ek)
!!      | histT:   Historical record of time(T) (For MD or iteration for GO)
!!      | histR:   Historical record of rprimd(R)
!!      | histS:   Historical record of strten(S)
!!      | histV:   Historical record of velocity(V)
!!      | histXF:  Historical record of positions(X) and forces(F)
!! itime=index of the present iteration 
!! icycle=index of the present cycle
!!
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      ab_define_var,handle_ncerr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine write_md_hist(filename,hist,icycle,itime,natom)

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
#define ABI_FUNC 'write_md_hist'
 use interfaces_61_ionetcdf, except_this_one => write_md_hist
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filename 
 type(ab_movehistory),intent(in) :: hist
 integer,intent(in) :: itime,icycle,natom
!Local variables-------------------------------
!scalars
 integer :: ncerr,ncid
 integer :: xyz_id,natom_id,time_id,six_id
 integer :: xcart_id,xred_id,fcart_id,fred_id
 integer :: vel_id,etotal_id,acell_id,rprimd_id,strten_id,ekin_id,mdtime_id
!arrays
 integer :: dimXFids(3),dimAids(2),dimEids(1),dimRids(3),dimSids(2)
 integer :: count1(1),start1(1)
 integer :: count2(2),start2(2)
 integer :: count3(3),start3(3)

! *************************************************************************

!fname_hist=trim(ab_mover%filnam_ds(4))//'_HIST'

 write(std_out,*) 'ihist @ write_md_hist',hist%ihist
 write(std_out,*) 'mxhist @ write_md_hist',hist%mxhist

#if defined HAVE_TRIO_NETCDF

!#####################################################################
!### Creation of NetCDF file

 if (itime==1.and.icycle==1) then
!  1. Create netCDF file
   ncerr = nf90_create(path=trim(filename),&
&   cmode=NF90_CLOBBER, ncid=ncid)
   if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," create netcdf history file")
!  2. Define dimensions
   ncerr = nf90_def_dim(ncid,"natom",natom,natom_id)
   if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," define dimension natom")
   ncerr = nf90_def_dim(ncid,"xyz",3,xyz_id)
   if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," define dimension xyz")
   ncerr = nf90_def_dim(ncid,"time",NF90_UNLIMITED,time_id)
   if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," define dimension time")
   ncerr = nf90_def_dim(ncid,"six",6,six_id)
   if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," define dimension six")

!  Dimensions for xcart,xred,fcart,fred and vel
   dimXFids = (/ xyz_id, natom_id, time_id /)
!  Dimensions for acell
   dimAids = (/ xyz_id, time_id /)
!  Dimensions for etotal
   dimEids = (/ time_id /)
!  Dimensions for rprimd
   dimRids = (/ xyz_id, xyz_id, time_id /)
!  Dimensions for strten
   dimSids = (/ six_id, time_id /)

!  3. Define variables and their attributes (units and mnemonics)
   call ab_define_var(ncid, dimXFids, xcart_id, NF90_DOUBLE,&
&   "xcart",&
&   "vectors (X) of atom positions in CARTesian coordinates",&
&   "bohr" )
   call ab_define_var(ncid, dimXFids, xred_id, NF90_DOUBLE,&
&   "xred",&
&   "vectors (X) of atom positions in REDuced coordinates",&
&   "dimensionless" )
   call ab_define_var(ncid, dimXFids, fcart_id, NF90_DOUBLE,&
&   "fcart",&
&   "atom Forces in CARTesian coordinates",&
&   "Ha/bohr" )
   call ab_define_var(ncid, dimXFids, fred_id, NF90_DOUBLE,&
&   "fred",&
&   "atom Forces in REDuced coordinates",&
&   "dimensionless" )
   call ab_define_var(ncid, dimXFids, vel_id, NF90_DOUBLE,&
&   "vel",&
&   "VELocity",&
&   "bohr*Ha/hbar" )
   call ab_define_var(ncid, dimAids, acell_id, NF90_DOUBLE,&
&   "acell",&
&   "CELL lattice vector scaling",&
&   "bohr" )
   call ab_define_var(ncid, dimRids, rprimd_id, NF90_DOUBLE,&
&   "rprimd",&
&   "Real space PRIMitive translations, Dimensional",&
&   "bohr" )
   call ab_define_var(ncid, dimEids, etotal_id, NF90_DOUBLE,&
&   "etotal",&
&   "TOTAL Energy",&
&   "Ha" )
   call ab_define_var(ncid, dimEids, ekin_id, NF90_DOUBLE,&
&   "ekin",&
&   "Energy KINetic ionic",&
&   "Ha" )
   call ab_define_var(ncid, dimEids, mdtime_id, NF90_DOUBLE,&
&   "mdtime",&
&   "Molecular Dynamics TIME",&
&   "hbar/Ha" )
   call ab_define_var(ncid, dimSids, strten_id, NF90_DOUBLE,&
&   "strten",&
&   "STRess tensor",&
&   "Ha/bohr^3" )

!  4. End define mode
   ncerr = nf90_enddef(ncid)
   if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," end define mode")

!  #####################################################################
!  ### Open the NetCDF file for write new iterations
 else
   write(std_out,*) 'OPEN NETCDF FILE'
!  1. Open netCDF file
   ncerr = nf90_open(path=trim(filename),&
&   mode=NF90_WRITE, ncid=ncid)
   if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," open netcdf history file")
!  2. Get the ID of a variables from their name
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

 end if

!#####################################################################
!### Write variables into the dataset

 count3 = (/ 3, natom, 1 /)
 start3 = (/ 1, 1, 1 /)
 start3(3)=hist%ihist

 ncerr = nf90_put_var(ncid, xcart_id,&
& hist%histXF(:,:,1,hist%ihist),&
& start = start3,count = count3)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable xcart")
 ncerr = nf90_put_var(ncid, xred_id,&
& hist%histXF(:,:,2,hist%ihist),&
& start = start3,count = count3)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable xred")
 ncerr = nf90_put_var(ncid, fcart_id,&
& hist%histXF(:,:,3,hist%ihist),&
& start = start3,count = count3)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable fcart")
 ncerr = nf90_put_var(ncid, fred_id,&
& hist%histXF(:,:,4,hist%ihist),&
& start = start3,count = count3)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable fred")
 ncerr = nf90_put_var(ncid, vel_id,&
& hist%histV(:,:,hist%ihist),&
& start = start3,count = count3)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable vel")

!RPRIMD
 count3(2)=3
 ncerr = nf90_put_var(ncid, rprimd_id,&
& hist%histR(:,:,hist%ihist),&
& start = start3,count = count3)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable rprimd")

!ACELL
 count2 = (/ 3, 1 /)
 start2 = (/ 1, 1 /)
 start2(2)=hist%ihist
 ncerr = nf90_put_var(ncid, acell_id,&
& hist%histA(:,hist%ihist),&
& start = start2,count = count2)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable acell")

!STRTEN
 count2(1)=6
 ncerr = nf90_put_var(ncid, strten_id,&
& hist%histS(:,hist%ihist),&
& start = start2,count = count2)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable strten")

!ETOTAL
 count1 = (/ 1 /)
 start1 = (/ 1 /)
 start1(1)=hist%ihist
 ncerr = nf90_put_var(ncid, etotal_id,&
& hist%histE(hist%ihist), start = (/ hist%ihist /) )
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable etotal")

!Ekin
 count1 = (/ 1 /)
 start1 = (/ 1 /)
 start1(1)=hist%ihist
 ncerr = nf90_put_var(ncid, ekin_id,&
& hist%histEk(hist%ihist), start = (/ hist%ihist /) )
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable ekin")

!MDTIME
 count1 = (/ 1 /)
 start1 = (/ 1 /)
 start1(1)=hist%ihist
 ncerr = nf90_put_var(ncid, mdtime_id,&
& hist%histT(hist%ihist), start = (/ hist%ihist /) )
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," write variable mdtime")

!#####################################################################
!### Close NetCDF file
 ncerr = nf90_close(ncid)
 if(ncerr /= NF90_NOERR) call handle_ncerr(ncerr," close netcdf history file")

#endif

end subroutine write_md_hist
!!***
