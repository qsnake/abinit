!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_io_netcdf_wfftype
!! NAME
!! hdr_io_netcdf_wfftype
!!
!! FUNCTION
!! This subroutine deals with the I/O of the hdr_type
!! structured variables (read/write/echo).
!! It has been modified to output netCDF format files for the
!! binary calls, as an alternative to standard hdr_io.f
!! According to the value of rdwr, it reads the header
!! of a file, writes it, or echo the value of the structured
!! variable to a file.
!! Note that, when reading, different records of hdr
!! are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated
!! correctly by a call to hdr_clean when hdr is not used anymore.
!! Two instances of the hdr_io_netcdf routines are defined :
!!  hdr_io_netcdf_int to which only the unit number is given
!!  hdr_io_netcdf_wfftype to which a wffil datatype is given
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (XG,Mver)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  rdwr= if 1, read the hdr structured variable from the header of the netCDF file,
!!        if 2, write the header to unformatted netCDF file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!        if 5, read the hdr without rewinding (unformatted), identical to 1 for netCDF
!!        if 6, read the hdr without rewinding (unformatted), identical to 2 for netCDF
!!  unitfi=unit number of the file (unformatted if rdwr=1, 2, 5 or 6 formatted if rdwr=3,4)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!  hdr <type(hdr_type)>=the header structured variable
!!   if rdwr=1,5 : will be output
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!
!! NOTES
!! In all cases, the file is supposed to be open already
!! When reading (rdwr=1) or writing (rdwr=2), rewind the file
!! When echoing (rdwr=3) does not rewind the file.
!! When reading (rdwr=5) or writing (rdwr=6), DOES NOT rewind the file
!!
!! NetCDF: the netcdf file should be in data mode (not define) coming
!!   in to and leaving this routine.
!!   To read:
!!     1) inq the dimension id
!!     2) Inquire the value
!!     1) inq the variable id
!!     2) get the variable value
!!   To write:
!!     if the object is undefined
!!      1) define the dimension or variable, and get its id.
!!         Dimensions are put here.
!!     otherwise
!!      1) inq the dimension or variable id
!!     end if
!!     2) put the variable values
!!
!! PARENTS
!!
!! CHILDREN
!!      handle_ncerr,hdr_io_int,leave_new,rhoij_alloc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine hdr_io_netcdf_wfftype(fform,hdr,rdwr,wff)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_wffile

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_io_netcdf_wfftype'
 use interfaces_16_hideleave
 use interfaces_59_io_mpi, except_this_one => hdr_io_netcdf_wfftype
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(inout) :: fform
 integer,intent(in) :: rdwr
 type(hdr_type),intent(inout) :: hdr
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: ierr
#endif

! *************************************************************************

!DEBUG
!write(std_out,*)' hdr_io_netcdf : enter hdr_io_netcdf_wfftype '
!call flush(6)
!ENDDEBUG
#if !defined HAVE_TRIO_NETCDF
 write(std_out,*) 'hdr_io_netcdf Error : NETCDF undefined at compile time.'
 call leave_new('COLL')
#endif

 if( wff%accesswff==IO_MODE_NETCDF .and. wff%master==wff%me ) then
   call hdr_io_netcdf_int(fform,hdr,rdwr,wff%unwff)
 end if

!NORMALLY THIS SHOULD NOT BE HERE: netCDF not in parallel yet
#if defined HAVE_MPI
 write(std_out,*) 'hdr_io_netcdf Error : MPI defined at compile time. Netcdf does not work with MPI yet.'
 stop
!In the parallel case, if the files were not local, need to bcast the data
 if(rdwr==1)then
   if (wff%accesswff==IO_MODE_FORTRAN_MASTER .or. wff%accesswff==IO_MODE_MPI) then
     call MPI_BCAST(fform,1,MPI_INTEGER,wff%master,wff%spaceComm,ierr)
     call hdr_comm(hdr,wff%master,wff%me,wff%spaceComm)
!    XG 040121 : I am not sure that the next line is correct
     if(wff%accesswff==IO_MODE_MPI .and. wff%master/=wff%me)then
       call hdr_skip_wfftype(wff,ierr)
     end if
   end if
 end if
#endif

end subroutine hdr_io_netcdf_wfftype

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

 subroutine hdr_io_netcdf_int(fform,hdr,rdwr,unitfi)

 use m_profiling

 use defs_basis
 use defs_abitypes
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_io_netcdf_int'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_44_abitypes_defs
 use interfaces_59_io_mpi, except_this_one => hdr_io_netcdf_int
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(inout) :: fform
 integer,intent(in) :: rdwr,unitfi
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
!no_abirules
#if defined HAVE_TRIO_NETCDF
 integer :: cplex,headform,iatom,ikpt,nband_kptsppol
 integer :: isppol,iocc

!DEBUG
!character(len=500) :: tmpname
!integer :: ivar
!ENDDEBUG

 integer :: ilmn,irhoij,ispden,itypat,lmn2_size,lmn_size_max,nselect,rhoijdim1,rhoijdim2

 character(len=500) :: message
 character(len=6) :: codvsn

!netCDF file id and error code
 integer :: ncid_hdr, ncerr

!id for dimensions
 integer :: bantot_id,codvsnlen_id,dimr3_id,dimg3_id,dim1scalar_id
 integer :: mband_id,natom_id,nkpt_id,nspden_id,nspinor_id
 integer :: nsppol_id,nsym_id,npsp_id,ntypat_id
 integer :: psptitlen_id
!headform >= 44
 integer :: rhoijdim1_id,rhoijdim2_id

!id for scalar variables:
 integer :: codvsn_id,date_id,ecut_id,ecut_eff_id,ecutsm_id
 integer :: headform_id,fform_id,intxc_id,ixc_id,ngfft_id
 integer :: occopt_id,pertcase_id,residm_id,stmbias_id,tphysel_id,tsmear_id
!headform >= 44
 integer :: usepaw_id, ecutdg_id

!id for allocatable variables
 integer :: qptn_id,rprimd_id,istwfk_id,nband_id,npwarr_id
 integer :: so_psp_id,symafm_id,symrel_id,typat_id
 integer :: kptns_id,occ_id,tnons_id,znucltypat_id
 integer :: title_id,znuclpsp_id,zionpsp_id,pspso_id,pspdat_id,pspcod_id,pspxc_id
 integer :: xred_id,etot_id,fermie_id
!headform >= 44
 integer :: lmn_size_id,rhoij_id
!headform >= 53
 integer :: wtk_id


!temp variables
 integer :: mband_dummy
!integer, allocatable :: nband_tmp(:,:)
 real(dp), allocatable :: occ_1kpt(:)
 real(dp), allocatable :: rhoij(:,:,:)

 character(len=500) :: dummy_name
#endif

!*************************************************************************

!DEBUG
!write(std_out,*)' hdr_io_netcdf : enter hdr_io_netcdf_int '
!call flush(6)
!ENDDEBUG

#if defined HAVE_TRIO_NETCDF
!-------------------------------------------------------------------------
!Reading the header of an unformatted file
!-------------------------------------------------------------------------

 if(rdwr==1 .or. rdwr==5)then

!  Reading the first record of the file ------------------------------------

!  PRESUME FILE IS OPEN AND unitfi IS ACTUALLY ncid
!  
!  !  should we pass the name of the file to this subroutine,
!  !    or do the nf90_open in the same place as the normal file open?
!  ncerr = nf90_open(path="tmpfile.nc", mode=NF90_NOWRITE, ncid=ncid_hdr)
!  if (ncerr /= nf90_noerr)  write(std_out,*) trim(nf90_strerror(ncerr))

   ncid_hdr = unitfi

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="codvsn",varid=codvsn_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire codvsn")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=codvsn_id,values=codvsn)
   call handle_ncerr(ncerr," hdr_io_netcdf : get codvsn")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="fform",varid=fform_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire fform")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=fform_id,values=fform)
   call handle_ncerr(ncerr," hdr_io_netcdf : get fform")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="headform",varid=headform_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire headform")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=headform_id,values=headform)
   call handle_ncerr(ncerr," hdr_io_netcdf : get headform")

   if(headform<=42)then
     write(message, '(4a,i3,3a,i8,3a)' ) ch10,&
&     ' hdr_io_netcdf : ERROR -',ch10,&
&     '  The first records of the netCDF (WF, DEN or POT) file read in unit ',unitfi,&
&     ' is erroneous.',ch10,&
&     '  headform is ',headform,', while it should be <= 42.',ch10,&
&     '  Action : check the correctness of your file.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   hdr%codvsn=codvsn
   hdr%headform=headform
!  fform is not a record of hdr_type

!  DEBUG
!  write(std_out,*)' hdr_io_netcdf : debug '
!  write(std_out,*)' hdr_io_netcdf : codvsn,headform,fform',codvsn,headform,fform
!  ENDDEBUG

!  =================================================
!  get main dimension attributes:
!  =================================================

!  get dimid number from nf90_inq_dimid,
!  then length from call nf90_Inquire_Dimension

   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="bantot",dimid=bantot_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire bantot")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=bantot_id,&
&   name=dummy_name,len=hdr%bantot)
   call handle_ncerr(ncerr," hdr_io_netcdf : get bantot")

   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="mband",dimid=mband_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire mband")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=mband_id,&
&   name=dummy_name,len=mband_dummy)
   call handle_ncerr(ncerr," hdr_io_netcdf : get mband")

   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="natom",dimid=natom_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire natom")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=natom_id,&
&   name=dummy_name,len=hdr%natom)
   call handle_ncerr(ncerr," hdr_io_netcdf : get natom")

   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nkpt",dimid=nkpt_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire nkpt")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nkpt_id,&
&   name=dummy_name,len=hdr%nkpt)
   call handle_ncerr(ncerr," hdr_io_netcdf : get nkpt")

   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="npsp",dimid=npsp_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire npsp")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=npsp_id,&
&   name=dummy_name,len=hdr%npsp)
   call handle_ncerr(ncerr," hdr_io_netcdf : get npsp")

   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nspden",dimid=nspden_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire nspden")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nspden_id,&
&   name=dummy_name,len=hdr%nspden)
   call handle_ncerr(ncerr," hdr_io_netcdf : get nspden")

   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nspinor",dimid=nspinor_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire nspinor")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nspinor_id,&
&   name=dummy_name,len=hdr%nspinor)
   call handle_ncerr(ncerr," hdr_io_netcdf : get nspinor")

   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nsppol",dimid=nsppol_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire nsppol")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nsppol_id,&
&   name=dummy_name,len=hdr%nsppol)
   call handle_ncerr(ncerr," hdr_io_netcdf : get nsppol")

   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nsym",dimid=nsym_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire nsym")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nsym_id,&
&   name=dummy_name,len=hdr%nsym)
   call handle_ncerr(ncerr," hdr_io_netcdf : get nsym")

   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="ntypat",dimid=ntypat_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire ntypat")
   ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=ntypat_id,&
&   name=dummy_name,len=hdr%ntypat)
   call handle_ncerr(ncerr," hdr_io_netcdf : get ntypat")

   if (headform >= 44) then
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="rhoijdim1",dimid=rhoijdim1_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire rhoijdim1")
     ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=rhoijdim1_id,&
&     name=dummy_name,len=rhoijdim1)
     call handle_ncerr(ncerr," hdr_io_netcdf : get rhoijdim1")
   end if
   if (headform >= 57) then
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="rhoijdim2",dimid=rhoijdim2_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire rhoijdim2")
     ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=rhoijdim2_id,&
&     name=dummy_name,len=rhoijdim2)
     call handle_ncerr(ncerr," hdr_io_netcdf : get rhoijdim2")
   else
     rhoijdim2=hdr%nspden
   end if

!  =================================================
!  variables which dont need to be allocated
!  alphabetical order
!  =================================================

!  get varid from nf90_inq_varid (and check for existence in passing)
!  get values from nf90_get_var

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="date",varid=date_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire date")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=date_id,values=hdr%date)
   call handle_ncerr(ncerr," hdr_io_netcdf : get date")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="ecut",varid=ecut_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire ecut")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=ecut_id,values=hdr%ecut)
   call handle_ncerr(ncerr," hdr_io_netcdf : get ecut")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="ecut_eff",varid=ecut_eff_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire ecut_eff")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=ecut_eff_id,values=hdr%ecut_eff)
   call handle_ncerr(ncerr," hdr_io_netcdf : get ecut_eff")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="ecutsm",varid=ecutsm_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire ecutsm")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=ecutsm_id,values=hdr%ecutsm)
   call handle_ncerr(ncerr," hdr_io_netcdf : get ecutsm")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="etot",varid=etot_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire etot")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=etot_id,values=hdr%etot)
   call handle_ncerr(ncerr," hdr_io_netcdf : get etot")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="fermie",varid=fermie_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire fermie")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=fermie_id,values=hdr%fermie)
   call handle_ncerr(ncerr," hdr_io_netcdf : get fermie")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="intxc",varid=intxc_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire intxc")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=intxc_id,values=hdr%intxc)
   call handle_ncerr(ncerr," hdr_io_netcdf : get intxc")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="ixc",varid=ixc_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire ixc")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=ixc_id,values=hdr%ixc)
   call handle_ncerr(ncerr," hdr_io_netcdf : get ixc")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="ngfft",varid=ngfft_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire ngfft")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=ngfft_id,values=hdr%ngfft)
   call handle_ncerr(ncerr," hdr_io_netcdf : get ngfft")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="occopt",varid=occopt_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire occopt")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=occopt_id,values=hdr%occopt)
   call handle_ncerr(ncerr," hdr_io_netcdf : get occopt")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="pertcase",varid=pertcase_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire pertcase")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=pertcase_id,values=hdr%pertcase)
   call handle_ncerr(ncerr," hdr_io_netcdf : get pertcase")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="qptn",varid=qptn_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire qptn")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=qptn_id,values=hdr%qptn)
   call handle_ncerr(ncerr," hdr_io_netcdf : get qptn")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="residm",varid=residm_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire residm")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=residm_id,values=hdr%residm)
   call handle_ncerr(ncerr," hdr_io_netcdf : get residm")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="rprimd",varid=rprimd_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire rprimd")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=rprimd_id,values=hdr%rprimd)
   call handle_ncerr(ncerr," hdr_io_netcdf : get rprimd")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="stmbias",varid=stmbias_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire stmbias")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=stmbias_id,values=hdr%stmbias)
   call handle_ncerr(ncerr," hdr_io_netcdf : get stmbias")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="tphysel",varid=tphysel_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire tphysel")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=tphysel_id,values=hdr%tphysel)
   call handle_ncerr(ncerr," hdr_io_netcdf : get tphysel")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="tsmear",varid=tsmear_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire tsmear")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=tsmear_id,values=hdr%tsmear)
   call handle_ncerr(ncerr," hdr_io_netcdf : get tsmear")

!  Compared to v4.2, add usepaw and ecutdg
   if (headform >= 44) then
     ncerr = nf90_inq_varid(ncid=ncid_hdr,name="usepaw",varid=usepaw_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire usepaw")
     ncerr = nf90_get_var(ncid=ncid_hdr,varid=usepaw_id,values=hdr%usepaw)
     call handle_ncerr(ncerr," hdr_io_netcdf : get usepaw")

     ncerr = nf90_inq_varid(ncid=ncid_hdr,name="ecutdg",varid=ecutdg_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire ecutdg")
     ncerr = nf90_get_var(ncid=ncid_hdr,varid=ecutdg_id,values=hdr%ecutdg)
     call handle_ncerr(ncerr," hdr_io_netcdf : get ecutdg")
   end if

!  test for old wavefunction style
   if(hdr%ecutsm>tol6 .and. headform<44 .and. .not.(fform==51.or.fform==52.or.fform==101.or.fform==102))then
     write(message, '(4a,es16.6,9a)' ) ch10,&
&     ' hdr_io_netcdf : ERROR -',ch10,&
&     '  The value of ecutsm is',hdr%ecutsm,', while the file has been produced prior to v4.4 .',ch10,&
&     '  The definition of the smearing function has changed, so that you are not allowed',ch10,&
&     '  to restart from a old wavefunction file. By contrast, you can restart from an old',ch10,&
&     '  potential or density file, and perform a self-consistent cycle with a new ABINIT version.',ch10,&
&     '  Action : produce a density or potential file using the old version of ABINIT, and restart from it.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if


!  DEBUG
!  write(std_out,*)' hdr_io_netcdf : before allocate '
!  write(std_out,*)' hdr_io_netcdf : bantot,natom,nkpt,npsp,nsppol,nsym,ntypat',&
!  &  hdr%bantot,hdr%natom,hdr%nkpt,hdr%npsp,hdr%nsppol,hdr%nsym,hdr%ntypat
!  ENDDEBUG

!  Allocate all parts of hdr that need to be --------------------------------

   ABI_ALLOCATE(hdr%istwfk,(hdr%nkpt))
   ABI_ALLOCATE(hdr%lmn_size,(hdr%npsp))
   ABI_ALLOCATE(hdr%nband,(hdr%nkpt*hdr%nsppol))
   ABI_ALLOCATE(hdr%npwarr,(hdr%nkpt))
   ABI_ALLOCATE(hdr%pspcod,(hdr%npsp))
   ABI_ALLOCATE(hdr%pspdat,(hdr%npsp))
   ABI_ALLOCATE(hdr%pspso,(hdr%npsp))
   ABI_ALLOCATE(hdr%pspxc,(hdr%npsp))
   ABI_ALLOCATE(hdr%so_psp,(hdr%npsp))
   ABI_ALLOCATE(hdr%symafm,(hdr%nsym))
   ABI_ALLOCATE(hdr%symrel,(3,3,hdr%nsym))
   ABI_ALLOCATE(hdr%typat,(hdr%natom))
   ABI_ALLOCATE(hdr%kptns,(3,hdr%nkpt))
   ABI_ALLOCATE(hdr%occ,(hdr%bantot))
   ABI_ALLOCATE(hdr%tnons,(3,hdr%nsym))
   ABI_ALLOCATE(hdr%wtk,(hdr%nkpt))
   ABI_ALLOCATE(hdr%xred,(3,hdr%natom))
   ABI_ALLOCATE(hdr%znuclpsp,(hdr%npsp))
   ABI_ALLOCATE(hdr%znucltypat,(hdr%ntypat))
   ABI_ALLOCATE(hdr%zionpsp,(hdr%npsp))
   ABI_ALLOCATE(hdr%title,(hdr%npsp))
   if(hdr%usepaw==1)  then
     ABI_ALLOCATE(hdr%pawrhoij,(hdr%natom))
   end if


!  DEBUG
!  do ivar=1,46
!  ncerr = nf90_Inquire_Variable(ncid=ncid_hdr,varid=ivar,name=tmpname)
!  write(std_out,*) 'hdr_io_netcdf : ivar,name = ', ivar,trim(tmpname)
!  end do
!  write(std_out,*)' hdr_io_netcdf : after allocate '
!  ENDDEBUG

!  ===========================================================
!  variables which need to be allocated (dimensions above)
!  alphabetical order
!  ===========================================================

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="istwfk",varid=istwfk_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire istwfk")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=istwfk_id,values=hdr%istwfk)
   call handle_ncerr(ncerr," hdr_io_netcdf : get istwfk")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="kptns",varid=kptns_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire kptns")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=kptns_id,values=hdr%kptns)
   call handle_ncerr(ncerr," hdr_io_netcdf : get kptns")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="wtk",varid=wtk_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire wtk")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=wtk_id,values=hdr%wtk)
   call handle_ncerr(ncerr," hdr_io_netcdf : get wtk")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="nband",varid=nband_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire nband")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=nband_id,values=hdr%nband,&
&   start=(/1,1/),count=(/hdr%nkpt,hdr%nsppol/))
   call handle_ncerr(ncerr," hdr_io_netcdf : get nband")

!  !allocate (nband_tmp(hdr%nsppol,hdr%nkpt))
!  allocate (nband_tmp(hdr%nkpt,hdr%nsppol))
!  ncerr = nf90_get_var(ncid=ncid_hdr,varid=nband_id,values=nband_tmp)
!  call handle_ncerr(ncerr," hdr_io_netcdf : get nband")
!  
!  write(std_out,*) 'read in nband_tmp = ', nband_tmp
!  hdr%nband = reshape(nband_tmp,(/hdr%nkpt*hdr%nsppol/))
!  write(std_out,*) 'allocated(hdr%nband) = ', allocated(hdr%nband)
!  write(std_out,*) 'shape(hdr%nband) = ', shape(hdr%nband)
!  write(std_out,*) 'kind(hdr%nband) = ', kind(hdr%nband)
!  write(std_out,*) 'hdr%nband = ', hdr%nband
!  deallocate (nband_tmp)

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="npwarr",varid=npwarr_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire npwarr")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=npwarr_id,values=hdr%npwarr)
   call handle_ncerr(ncerr," hdr_io_netcdf : get npwarr")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="occ",varid=occ_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire occ")
   ABI_ALLOCATE(occ_1kpt,(maxval(hdr%nband)))
   iocc = 0
   do isppol=1,hdr%nsppol
     do ikpt=1,hdr%nkpt
       nband_kptsppol = hdr%nband((isppol-1)*hdr%nkpt + ikpt)
       ncerr = nf90_get_var(ncid=ncid_hdr,varid=occ_id,values=occ_1kpt,&
&       start=(/1,ikpt,isppol/),count=(/nband_kptsppol,1,1/))
       call handle_ncerr(ncerr," hdr_io_netcdf : get occ")
       hdr%occ(iocc+1:iocc+nband_kptsppol) = occ_1kpt(1:nband_kptsppol)
       iocc = iocc + nband_kptsppol
     end do
   end do
   ABI_DEALLOCATE(occ_1kpt)


   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="pspcod",varid=pspcod_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire pspcod")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=pspcod_id,values=hdr%pspcod)
   call handle_ncerr(ncerr," hdr_io_netcdf : get pspcod")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="pspdat",varid=pspdat_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire pspdat")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=pspdat_id,values=hdr%pspdat)
   call handle_ncerr(ncerr," hdr_io_netcdf : get pspdat")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="pspso",varid=pspso_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire pspso")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=pspso_id,values=hdr%pspso)
   call handle_ncerr(ncerr," hdr_io_netcdf : get pspso")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="pspxc",varid=pspxc_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire pspxc")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=pspxc_id,values=hdr%pspxc)
   call handle_ncerr(ncerr," hdr_io_netcdf : get pspxc")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="so_psp",varid=so_psp_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire so_psp")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=so_psp_id,values=hdr%so_psp)
   call handle_ncerr(ncerr," hdr_io_netcdf : get so_psp")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="symafm",varid=symafm_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire symafm")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=symafm_id,values=hdr%symafm)
   call handle_ncerr(ncerr," hdr_io_netcdf : get symafm")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="symrel",varid=symrel_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire symrel")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=symrel_id,values=hdr%symrel)
   call handle_ncerr(ncerr," hdr_io_netcdf : get symrel")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="title",varid=title_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire title")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=title_id,values=hdr%title)
   call handle_ncerr(ncerr," hdr_io_netcdf : get title")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="typat",varid=typat_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire typat")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=typat_id,values=hdr%typat)
   call handle_ncerr(ncerr," hdr_io_netcdf : get typat")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="tnons",varid=tnons_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire tnons")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=tnons_id,values=hdr%tnons)
   call handle_ncerr(ncerr," hdr_io_netcdf : get tnons")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="xred",varid=xred_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire xred")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=xred_id,values=hdr%xred)
   call handle_ncerr(ncerr," hdr_io_netcdf : get xred")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="zionpsp",varid=zionpsp_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire zionpsp")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=zionpsp_id,values=hdr%zionpsp)
   call handle_ncerr(ncerr," hdr_io_netcdf : get zionpsp")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="znuclpsp",varid=znuclpsp_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire znuclpsp")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=znuclpsp_id,values=hdr%znuclpsp)
   call handle_ncerr(ncerr," hdr_io_netcdf : get znuclpsp")

   ncerr = nf90_inq_varid(ncid=ncid_hdr,name="znucltypat",varid=znucltypat_id)
   call handle_ncerr(ncerr," hdr_io_netcdf : inquire znucltypat")
   ncerr = nf90_get_var(ncid=ncid_hdr,varid=znucltypat_id,values=hdr%znucltypat)
   call handle_ncerr(ncerr," hdr_io_netcdf : get znucltypat")

!  Compared to 4.2, add lmn_size and
   if (headform>=44) then
     ncerr = nf90_inq_varid(ncid=ncid_hdr,name="lmn_size",varid=lmn_size_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire lmn_size")
     ncerr = nf90_get_var(ncid=ncid_hdr,varid=lmn_size_id,values=hdr%lmn_size)
     call handle_ncerr(ncerr," hdr_io_netcdf : get lmn_size")

     if (hdr%usepaw==1) then
       ABI_ALLOCATE(rhoij,(rhoijdim1,rhoijdim2,hdr%natom))
       ncerr = nf90_inq_varid(ncid=ncid_hdr,name="rhoij",varid=rhoij_id)
       call handle_ncerr(ncerr," hdr_io_netcdf : inquire rhoij")
       ncerr = nf90_get_var(ncid=ncid_hdr,varid=rhoij_id,values=rhoij)
       call handle_ncerr(ncerr," hdr_io_netcdf : get rhoij")
       lmn_size_max=maxval(hdr%lmn_size)
       lmn2_size=lmn_size_max*(lmn_size_max+1)/2
       cplex=min(rhoijdim1/lmn2_size,2)
       call rhoij_alloc(cplex,hdr%lmn_size,rhoijdim2,hdr%nspinor,hdr%nsppol,hdr%pawrhoij,hdr%typat)
       do iatom=1,hdr%natom
         itypat=hdr%typat(iatom)
         lmn2_size=hdr%lmn_size(itypat)*(hdr%lmn_size(itypat)+1)/2
         nselect=0
         if (cplex==1) then
           do ilmn=1,lmn2_size
             if (any(abs(rhoij(ilmn,:,iatom))>tol10)) then
               nselect=nselect+1
               hdr%pawrhoij(iatom)%rhoijselect(nselect)=ilmn
               do ispden=1,rhoijdim2
                 hdr%pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoij(ilmn,ispden,iatom)
               end do
             end if
           end do
         else
           do ilmn=1,lmn2_size
             if (any(abs(rhoij(2*ilmn-1:2*ilmn,:,iatom))>tol10)) then
               nselect=nselect+1
               hdr%pawrhoij(iatom)%rhoijselect(nselect)=ilmn
               do ispden=1,rhoijdim2
                 hdr%pawrhoij(iatom)%rhoijp(2*nselect-1,ispden)=rhoij(2*ilmn-1,ispden,iatom)
                 hdr%pawrhoij(iatom)%rhoijp(2*nselect  ,ispden)=rhoij(2*ilmn  ,ispden,iatom)
               end do
             end if
           end do
         end if
         if (nselect<lmn2_size) then
           hdr%pawrhoij(iatom)%rhoijselect(nselect+1:lmn2_size)=0
           do ispden=1,rhoijdim2
             hdr%pawrhoij(iatom)%rhoijp(cplex*nselect+1:cplex*lmn2_size,ispden)=zero
           end do
         end if
         hdr%pawrhoij(iatom)%nrhoijsel=nselect
       end do
       ABI_DEALLOCATE(rhoij)
     end if
   end if

!  ! close netCDF file
!  ncerr = nf90_close(ncid_hdr)
!  if (ncerr /= nf90_noerr) then
!  write(std_out,*) ' Error closing netCDF file'
!  stop
!  end if

!  DEBUG
!  write(std_out,*)' hdr_io_netcdf : read mode, hdr%so_psp(:), hdr%symafm(:)=',&
!  &                                 hdr%so_psp(:), hdr%symafm(:)
!  ENDDEBUG


!  -------------------------------------------------------------------------
!  Writing the header of an unformatted file
!  -------------------------------------------------------------------------

 else if(rdwr==2 .or. rdwr==6)then

!  natom,nkpt,npsp,ntypat... are not defined in this section :
!  always address then from hdr

!  should not be a problem for netCDF files: access is not sequential
!  if(rdwr==2) rewind(unitfi)


!  PRESUME FILE IS ALREADY OPEN (FUNCTION LIKE HDR_IO.F)
!  
!  ! Add CLOBBER?
!  ncerr = nf90_create(path="tmpfile.nc", cmode=NF90_NOCLOBBER, ncid=ncid_hdr)

   ncid_hdr = unitfi

!  -------------------------------------------------------------------------
!  It is possible that the file has already been written to before
!  in that case, the dimensions and variables do not have to be
!  defined, and we skip to the filling immediately
!  Should perhaps check if the dimensions have not changed... Should not happen!
!  go to definition mode
!  -------------------------------------------------------------------------
   ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="bantot",dimid=bantot_id)

   if (ncerr /= NF90_NOERR) then
     ncerr = nf90_redef (ncid=ncid_hdr)
     call handle_ncerr(ncerr," going back to define mode ")

!    Writing always use format 44
     headform=57

!    
!    types = NF90_DOUBLE, NF90_INT  NF90_CHAR
!    NF90_UNLIMITED
     write(std_out,*) 'NF90_UNLIMITED = ', NF90_UNLIMITED
     write(std_out,*) 'natom = ', hdr%natom

!    =================================================
!    dimensions
!    =================================================

!    bantot is used in arrays
     write(std_out,*) 'bantot = ', hdr%bantot
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="bantot",len=hdr%bantot,dimid=bantot_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define bantot")

!    add dimension for length of codvsn string
     write(std_out,*) 'codvsnlen = ', len(hdr%codvsn)
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="codvsnlen",len=len(hdr%codvsn),&
&     dimid=codvsnlen_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define codvsnlen")

!    version 44 add first and second dimensions for rhoij
     rhoijdim1 = maxval(hdr%lmn_size)
     rhoijdim1 = rhoijdim1*(rhoijdim1+1)/2
     if (hdr%usepaw==1) rhoijdim1=rhoijdim1*hdr%pawrhoij(1)%cplex
     rhoijdim1 = max (rhoijdim1,1) ! impose rhoijdim1 >= 1 : if 0, it defaults to NF90_UNLIMITED
     write(std_out,*) 'rhoijdim1 = ', rhoijdim1
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="rhoijdim1",len=rhoijdim1,&
&     dimid=rhoijdim1_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define rhoijdim1")
     rhoijdim2 = hdr%nspden
     if (hdr%usepaw==1) rhoijdim2=hdr%pawrhoij(1)%nspden
     rhoijdim2 = max (rhoijdim2,1) ! impose rhoijdim2 >= 1 : if 0, it defaults to NF90_UNLIMITED
     write(std_out,*) 'rhoijdim2 = ', rhoijdim2
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="rhoijdim2",len=rhoijdim2,&
&     dimid=rhoijdim2_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define rhoijdim2")

!    dimensions for 3x3 matrix = 3 in recip and real space
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="dimr3",len=3,dimid=dimr3_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define dimr3")
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="dimg3",len=3,dimid=dimg3_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define dimg3")
!    reference dimension for scalars. Could probably be done away with
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="dim1scalar",len=1,dimid=dim1scalar_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define dim1scalar")

!    Now on to real physical dimensions used in header
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="mband",len=maxval(hdr%nband),dimid=mband_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define mband")
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="natom",len=hdr%natom,dimid=natom_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define natom")
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="nkpt",len=hdr%nkpt,dimid=nkpt_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define nkpt")
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="npsp",len=hdr%npsp,dimid=npsp_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define npsp")
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="nspden",len=hdr%nspden,dimid=nspden_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define nspden")
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="nspinor",len=hdr%nspinor,dimid=nspinor_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define nspinor")
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="nsppol",len=hdr%nsppol,dimid=nsppol_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define nsppol")
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="nsym",len=hdr%nsym,dimid=nsym_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define nsym")
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="ntypat",len=hdr%ntypat,dimid=ntypat_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define ntypat")

!    add dimension for length of psp titles
     ncerr = nf90_def_dim(ncid=ncid_hdr,name="psptitlen",len=fnlen,&
&     dimid=psptitlen_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define psptitlen")


!    ==========================================================
!    DEFINE constant scalars and strings
!    (some could be changed in favor of attributes)
!    ==========================================================

     ncerr = nf90_def_var(ncid=ncid_hdr,name="date",xtype=NF90_INT,&
&     dimids=(/dim1scalar_id/),varid=date_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define date")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="codvsn",xtype=NF90_CHAR,&
&     dimids=(/codvsnlen_id/),varid=codvsn_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define codvsn")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="ecut",xtype=NF90_DOUBLE,&
&     dimids=(/dim1scalar_id/),varid=ecut_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define ecut")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="ecut_eff",xtype=NF90_DOUBLE,&
&     dimids=(/dim1scalar_id/),varid=ecut_eff_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define ecut_eff")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="ecutsm",xtype=NF90_DOUBLE,&
&     dimids=(/dim1scalar_id/),varid=ecutsm_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define ecutsm")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="etot",xtype=NF90_DOUBLE,&
&     dimids=(/dim1scalar_id/),varid=etot_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define etot")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="fermie",xtype=NF90_DOUBLE,&
&     dimids=(/dim1scalar_id/),varid=fermie_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define fermie")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="headform",xtype=NF90_INT,&
&     dimids=(/dim1scalar_id/),varid=headform_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define headform")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="fform",xtype=NF90_INT,&
&     dimids=(/dim1scalar_id/),varid=fform_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define fform")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="intxc",xtype=NF90_INT,&
&     dimids=(/dim1scalar_id/),varid=intxc_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define intxc")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="ixc",xtype=NF90_INT,&
&     dimids=(/dim1scalar_id/),varid=ixc_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define ixc")

!    !!! May need to make this into 3 individual dimensions to be able
!    to use them in the definition of the wfk
     ncerr = nf90_def_var(ncid=ncid_hdr,name="ngfft",xtype=NF90_INT,&
&     dimids=(/dimg3_id/),varid=ngfft_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define ngfft")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="occopt",xtype=NF90_INT,&
&     dimids=(/dim1scalar_id/),varid=occopt_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define occopt")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="pertcase",xtype=NF90_INT,&
&     dimids=(/dim1scalar_id/),varid=pertcase_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define pertcase")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="qptn",xtype=NF90_DOUBLE,&
&     dimids=(/dimg3_id/),varid=qptn_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define qptn")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="residm",xtype=NF90_DOUBLE,&
&     dimids=(/dim1scalar_id/),varid=residm_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define residm")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="rprimd",xtype=NF90_DOUBLE,&
&     dimids=(/dimr3_id,dimr3_id/),varid=rprimd_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define rprimd")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="stmbias",xtype=NF90_DOUBLE,&
&     dimids=(/dim1scalar_id/),varid=stmbias_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define stmbias")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="tphysel",xtype=NF90_DOUBLE,&
&     dimids=(/dim1scalar_id/),varid=tphysel_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define tphysel")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="tsmear",xtype=NF90_DOUBLE,&
&     dimids=(/dim1scalar_id/),varid=tsmear_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define tsmear")

!    Version 44 add ecutdg and usepaw
     ncerr = nf90_def_var(ncid=ncid_hdr,name="ecutdg",xtype=NF90_DOUBLE,&
&     dimids=(/dim1scalar_id/),varid=ecutdg_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define ecutdg")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="usepaw",xtype=NF90_INT,&
&     dimids=(/dim1scalar_id/),varid=usepaw_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define usepaw")


!    =============================================================
!    DEFINE multidimensional variables
!    eventually add map between dimensions for eg occ or wfk
!    probably ok as is
!    =============================================================

     ncerr = nf90_def_var(ncid=ncid_hdr,name="istwfk",xtype=NF90_INT,&
&     dimids=(/nkpt_id/),varid=istwfk_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define istwfk")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="kptns",xtype=NF90_DOUBLE,&
&     dimids=(/dimg3_id,nkpt_id/),varid=kptns_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define kptns")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="wtk",xtype=NF90_DOUBLE,&
&     dimids=(/nkpt_id/),varid=wtk_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define wtk")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="nband",xtype=NF90_INT,&
&     dimids=(/nkpt_id,nsppol_id/),varid=nband_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define nband")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="npwarr",xtype=NF90_INT,&
&     dimids=(/nkpt_id/),varid=npwarr_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define npwarr")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="occ",xtype=NF90_DOUBLE,&
&     dimids=(/mband_id,nkpt_id,nsppol_id/),varid=occ_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define occ")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="pspcod",xtype=NF90_INT,&
&     dimids=(/npsp_id/),varid=pspcod_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define pspcod")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="pspdat",xtype=NF90_INT,&
&     dimids=(/npsp_id/),varid=pspdat_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define pspdat")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="pspso",xtype=NF90_INT,&
&     dimids=(/npsp_id/),varid=pspso_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define pspso")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="pspxc",xtype=NF90_INT,&
&     dimids=(/npsp_id/),varid=pspxc_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define pspxc")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="so_psp",xtype=NF90_INT,&
&     dimids=(/ntypat_id/),varid=so_psp_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define so_psp")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="symafm",xtype=NF90_INT,&
&     dimids=(/nsym_id/),varid=symafm_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define symafm")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="symrel",xtype=NF90_INT,&
&     dimids=(/dimr3_id,dimr3_id,nsym_id/),varid=symrel_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define symrel")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="title",xtype=NF90_CHAR,&
&     dimids=(/psptitlen_id,npsp_id/),varid=title_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define title")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="typat",xtype=NF90_INT,&
&     dimids=(/natom_id/),varid=typat_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define typat")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="tnons",xtype=NF90_DOUBLE,&
&     dimids=(/dimr3_id,nsym_id/),varid=tnons_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define tnons")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="xred",xtype=NF90_DOUBLE,&
&     dimids=(/dimr3_id,natom_id/),varid=xred_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define xred")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="zionpsp",xtype=NF90_DOUBLE,&
&     dimids=(/npsp_id/),varid=zionpsp_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define zionpsp")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="znuclpsp",xtype=NF90_DOUBLE,&
&     dimids=(/npsp_id/),varid=znuclpsp_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define znuclpsp")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="znucltypat",xtype=NF90_DOUBLE,&
&     dimids=(/ntypat_id/),varid=znucltypat_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define znucltypat")

!    Version 44 add lmn_size and rhoij
     ncerr = nf90_def_var(ncid=ncid_hdr,name="lmn_size",xtype=NF90_INT,&
&     dimids=(/ntypat_id/),varid=lmn_size_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define lmn_size")

     ncerr = nf90_def_var(ncid=ncid_hdr,name="rhoij",xtype=NF90_DOUBLE,&
&     dimids=(/rhoijdim1_id,nspden_id,natom_id/),varid=rhoij_id)
     call handle_ncerr(ncerr," hdr_io_netcdf : define rhoij")


!    End define mode and go to data mode
     ncerr = nf90_enddef(ncid=ncid_hdr)
     call handle_ncerr(ncerr," enddef call ")

!    if the dimensions and variables have already been defined,
!    aquire all the variable ids
   else

     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=date_id,name="date")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire date")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=codvsn_id,name="codvsn")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire codvsn")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=ecut_id,name="ecut")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire ecut")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=ecut_eff_id,name="ecut_eff")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire ecut_eff")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=ecutsm_id,name="ecutsm")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire ecutsm")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=etot_id,name="etot")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire etot")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=fermie_id,name="fermie")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire fermie")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=headform_id,name="headform")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire headform")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=fform_id,name="fform")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire fform")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=intxc_id,name="intxc")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire intxc")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=ixc_id,name="ixc")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire ixc")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=ngfft_id,name="ngfft")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire ngfft")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=occopt_id,name="occopt")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire occopt")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=pertcase_id,name="pertcase")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire pertcase")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=qptn_id,name="qptn")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire qptn")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=residm_id,name="residm")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire residm")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=rprimd_id,name="rprimd")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire rprimd")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=stmbias_id,name="stmbias")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire stmbias")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=tphysel_id,name="tphysel")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire tphysel")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=tsmear_id,name="tsmear")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire tsmear")

!    Version 44 add usepaw ecutdg
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=ecutdg_id,name="ecutdg")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire ecutdg")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=usepaw_id,name="usepaw")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire usepaw")

     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=istwfk_id,name="istwfk")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire istwfk")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=kptns_id,name="kptns")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire kptns")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=wtk_id,name="wtk")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire wtk")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=nband_id,name="nband")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire nband")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=npwarr_id,name="npwarr")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire npwarr")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=occ_id,name="occ")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire occ")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=pspcod_id,name="pspcod")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire pspcod")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=pspdat_id,name="pspdat")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire pspdat")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=pspso_id,name="pspso")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire pspso")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=pspxc_id,name="pspxc")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire pspxc")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=so_psp_id,name="so_psp")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire so_psp")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=symafm_id,name="symafm")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire symafm")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=symrel_id,name="symrel")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire symrel")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=title_id,name="title")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire title")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=typat_id,name="typat")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire typat")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=tnons_id,name="tnons")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire tnons")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=xred_id,name="xred")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire xred")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=zionpsp_id,name="zionpsp")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire zionpsp")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=znuclpsp_id,name="znuclpsp")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire znuclpsp")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=znucltypat_id,name="znucltypat")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire znucltypat")

!    Version 44 add lmn_size and rhoij
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=lmn_size_id,name="lmn_size")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire lmn_size")
     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=rhoij_id,name="rhoij")
     call handle_ncerr(ncerr," hdr_io_netcdf : inquire rhoij")

   end if
!  end if on dimensions and variables being defined already


!  ==========================================================
!  FILL constant scalars and strings
!  (some could be changed in favor of attributes)
!  ==========================================================

   ncerr = nf90_put_var(ncid=ncid_hdr,varid=date_id,values=hdr%date)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill date")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=codvsn_id,values=hdr%codvsn)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill codvsn")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=ecut_id,values=hdr%ecut)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill ecut")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=ecut_eff_id,values=hdr%ecut_eff)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill ecut_eff")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=ecutsm_id,values=hdr%ecutsm)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill ecutsm")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=etot_id,values=hdr%etot)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill etot")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=fermie_id,values=hdr%fermie)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill fermie")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=headform_id,values=headform)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill headform")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=fform_id,values=fform)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill fform")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=intxc_id,values=hdr%intxc)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill intxc")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=ixc_id,values=hdr%ixc)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill ixc")
!  !!! May need to make this into 3 individual dimensions to be able
!  to use them in the definition of the wfk
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=ngfft_id,values=hdr%ngfft(1:3))
   call handle_ncerr(ncerr," hdr_io_netcdf : fill ngfft")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=occopt_id,values=hdr%occopt)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill occopt")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=pertcase_id,values=hdr%pertcase)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill pertcase")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=qptn_id,values=hdr%qptn)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill qptn")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=residm_id,values=hdr%residm)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill residm")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=rprimd_id,values=hdr%rprimd)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill rprimd")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=stmbias_id,values=hdr%stmbias)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill stmbias")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=tphysel_id,values=hdr%tphysel)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill tphysel")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=tsmear_id,values=hdr%tsmear)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill tsmear")

!  Version 44 add usepaw ecutdg
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=ecutdg_id,values=hdr%ecutdg)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill ecutdg")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=usepaw_id,values=hdr%usepaw)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill usepaw")

!  =============================================================
!  FILL variables
!  eventually add map between dimensions for eg occ, nband
!  or wfk
!  probably ok as is
!  =============================================================

   ncerr = nf90_put_var(ncid=ncid_hdr,varid=istwfk_id,values=hdr%istwfk)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill istwfk")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=kptns_id,values=hdr%kptns)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill kptns")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=wtk_id,values=hdr%wtk)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill wtk")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=nband_id,values=hdr%nband,&
&   start=(/1,1/),count=(/hdr%nkpt,hdr%nsppol/))
   call handle_ncerr(ncerr," hdr_io_netcdf : fill nband")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=npwarr_id,values=hdr%npwarr)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill npwarr")
   ABI_ALLOCATE(occ_1kpt,(maxval(hdr%nband)))
   iocc = 0
   do isppol=1,hdr%nsppol
     do ikpt=1,hdr%nkpt
       nband_kptsppol = hdr%nband((isppol-1)*hdr%nkpt + ikpt)
       occ_1kpt(1:nband_kptsppol) = hdr%occ(iocc+1:iocc+nband_kptsppol)
       ncerr = nf90_put_var(ncid=ncid_hdr,varid=occ_id,values=occ_1kpt,&
&       start=(/1,ikpt,isppol/),count=(/nband_kptsppol,1,1/))
       call handle_ncerr(ncerr," hdr_io_netcdf : fill occ")
       iocc = iocc + nband_kptsppol
     end do
   end do
   ABI_DEALLOCATE(occ_1kpt)
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=pspcod_id,values=hdr%pspcod)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill pspcod")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=pspdat_id,values=hdr%pspdat)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill pspdat")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=pspso_id,values=hdr%pspso)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill pspso")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=pspxc_id,values=hdr%pspxc)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill pspxc")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=so_psp_id,values=hdr%so_psp)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill so_psp")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=symafm_id,values=hdr%symafm)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill symafm")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=symrel_id,values=hdr%symrel)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill symrel")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=title_id,values=hdr%title)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill title")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=typat_id,values=hdr%typat)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill typat")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=tnons_id,values=hdr%tnons)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill tnons")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=xred_id,values=hdr%xred)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill xred")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=zionpsp_id,values=hdr%zionpsp)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill zionpsp")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=znuclpsp_id,values=hdr%znuclpsp)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill znuclpsp")
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=znucltypat_id,values=hdr%znucltypat)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill znucltypat")

!  Version 44 add lmn_size and rhoij
   ncerr = nf90_put_var(ncid=ncid_hdr,varid=lmn_size_id,values=hdr%lmn_size)
   call handle_ncerr(ncerr," hdr_io_netcdf : fill lmn_size")

   if (hdr%usepaw == 1) then
     ABI_ALLOCATE(rhoij,(rhoijdim1,rhoijdim2,hdr%natom))
     do iatom=1,hdr%natom
       itypat=hdr%typat(iatom)
       lmn2_size = hdr%lmn_size(itypat)*(hdr%lmn_size(itypat)+1)/2
       cplex=hdr%pawrhoij(iatom)%cplex
       do ispden=1,rhoijdim2
         rhoij(1:cplex*lmn2_size,ispden,iatom)=zero
         if (cplex==1) then
           do irhoij=1,hdr%pawrhoij(iatom)%nrhoijsel
             ilmn=hdr%pawrhoij(iatom)%rhoijselect(irhoij)
             rhoij(ilmn,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(irhoij,ispden)
           end do
         else
           do irhoij=1,hdr%pawrhoij(iatom)%nrhoijsel
             ilmn=hdr%pawrhoij(iatom)%rhoijselect(irhoij)
             rhoij(2*ilmn-1,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(2*irhoij-1,ispden)
             rhoij(2*ilmn  ,ispden,iatom)=hdr%pawrhoij(iatom)%rhoijp(2*irhoij  ,ispden)
           end do
         end if
       end do
     end do
     ncerr = nf90_put_var(ncid=ncid_hdr,varid=rhoij_id,values=rhoij)
     call handle_ncerr(ncerr," hdr_io_netcdf : fill rhoij")
     ABI_DEALLOCATE(rhoij)
   end if

!  ! Normal code
!  
!  write(unitfi) hdr%codvsn, headform, fform
!  
!  write(unitfi) hdr%bantot, hdr%date, hdr%intxc, hdr%ixc, &
!  &  hdr%natom, hdr%ngfft(1:3), hdr%nkpt, &
!  &  hdr%nspden, hdr%nspinor, &
!  &  hdr%nsppol, hdr%nsym, hdr%npsp, hdr%ntypat, hdr%occopt, hdr%pertcase,&
!  &  hdr%ecut, hdr%ecutsm, hdr%ecut_eff, &
!  &  hdr%qptn, hdr%rprimd, hdr%stmbias, hdr%tphysel, hdr%tsmear
!  
!  write(unitfi) hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:),&
!  &   hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
!  &   hdr%typat(:), hdr%kptns(:,:), hdr%occ(:), &
!  &   hdr%tnons(:,:), hdr%znucltypat(:)

!  DEBUG
!  write(std_out,*)' hdr_io_netcdf : write psp record, headform= ',hdr%headform
!  ENDDEBUG

!  =================================================
!  psp related variables
!  =================================================
!  DEBUG
!  write(std_out,*)' hdr_io_netcdf : write mode, hdr%so_psp(:), hdr%symafm(:)=',&
!  &                                  hdr%so_psp(:), hdr%symafm(:)
!  ENDDEBUG


!  variables


!  ! close netCDF file
!  ncerr = nf90_close(ncid_hdr)
!  if (ncerr /= nf90_noerr) then
!  write(std_out,*) ' Error closing netCDF file'
!  stop
!  end if
!  -------------------------------------------------------------------------
!  Writing the header of a formatted file
!  unchanged for netCDF
!  -------------------------------------------------------------------------

 else if(rdwr==3 .or. rdwr==4)then

   call hdr_io_int(fform,hdr,rdwr,unitfi)

 end if ! choice read/write/echo

#endif


 return
 fform=0 ; return   ! This is to allow treatment of old epsm1 format

end subroutine hdr_io_netcdf_int
!!***
