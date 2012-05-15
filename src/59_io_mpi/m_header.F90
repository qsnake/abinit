!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_header
!! NAME
!! m_header
!!
!! FUNCTION
!! This module contains the definition of the abinit header (TODO)
!! and methods acting on the data type.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (XG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

MODULE m_header

 use m_profiling

 use defs_basis
#if defined HAVE_MPI2
 use mpi
#endif
 use m_xmpi
 use m_errors

 use m_copy,          only : deep_copy
 use defs_wvltypes,   only : wvl_internal_type
 use defs_datatypes,  only : bandstructure_type, pseudopotential_type, pawtab_type
 use defs_abitypes,   only : hdr_type, dataset_type

 implicit none

 private
!!**

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 public :: m_header_init           ! Initialize the global variable HDR_fforms stored in this module.
 public :: hdr_fform2ftype         ! Return the filetype flag, a string with the name of the file and the patch level from fform.
 public :: hdr_ftype2fform         ! Returns the last fform associated to the filetype ftype.
 public :: hdr_init                ! Initialize the header structured datatype and most of its content from dtset and psps.
 public :: hdr_init_lowlvl         ! Low level initialization method for Hdr (no dtset).
 public :: hdr_nullify             ! Set all pointers to null.
 public :: hdr_clean               ! Deallocates the components of the header.
 public :: hdr_copy                ! Deep copy of the Header.
 public :: hdr_get_nelect_byocc    ! Returns the number of electrons calculated from Hdr%occ
 public :: isknown_headform        ! Returns .TRUE. is headform is in HDR_KNOWN_HEADFORMS.
 public :: hdr_mpio_skip           ! Skip the abinit header using MPI-IO routines and returns the offset of the 
                                   ! first Fortran record after the header.

 integer,public,parameter :: HDR_NHEADFORMS=9
 ! Number of abinit header formats used so far.

 integer,public,parameter :: HDR_KNOWN_HEADFORMS( HDR_NHEADFORMS ) = (/23,34,40,41,42,44,53,56,57/)
 ! The list of headforms used so far.

 integer,public,parameter :: HDR_LATEST_HEADFORM = HDR_KNOWN_HEADFORMS( HDR_NHEADFORMS )
 ! The latest headform to be used for writing.

 ! Each filetype is denoted by the filetype flag HDR_* and a patch level.
 ! To add a new filetype do the following:
 !   1) add a new integer filetype flag
 !   2) increase HDR_NFTYPES 
 !   3) Modify m_header_init to add a new entry in the array HDR_fforms.

 integer,public,parameter ::   &
&  HDR_WF_PW      = 1,         &
&  HDR_DENSITY    = 2,         &
&  HDR_POTENTIAL  = 3,         &
&  HDR_WF_WVL     = 4,         &
&  HDR_SCREENING  = 5,         &
&  HDR_BS_HAM     = 6,         &
&  HDR_BS_EIG     = 7,         &
&  HDR_BS_HAYDOCK = 8,         &
&  HDR_NFTYPES    = HDR_BS_HAYDOCK

 type,private :: fform_t
   !integer :: ftype
   character(len=500) :: fname
   integer,pointer :: fforms(:)
 end type fform_t

 type(fform_t),private,save,allocatable :: HDR_fforms(:)
 ! HDR_fforms(HDR_NFTYPES)
 ! The internal databases with the correspondence ftype ==> (fform, patch_level)

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_header/m_header_init
!! NAME
!! m_header_init
!!
!! FUNCTION
!!  Initialize the global variable HDR_fforms stored in this module.
!!
!! OUTPUT
!!  ierr=A nonzero values signals an inconsistency in the internal tables.
!!
!! SIDE EFFECTS 
!!   HDR_fforms(HDR_NFTYPES): Database completely initialized.
!!
!! PARENTS
!!      abi_init_globals
!!
!! CHILDREN
!!      mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine m_header_init(ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'm_header_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ierr
!Local variables-------------------------------
!scalars
 integer :: ft,nft

! *************************************************************************

 !@fform_t
 nft = HDR_NFTYPES
 ABI_ALLOCATE(HDR_fforms,(nft))

 do ft=1,nft
   select case(ft)

   case (HDR_WF_PW)
     HDR_fforms(ft)%fname = "wf_planewave"
     ABI_ALLOCATE(HDR_fforms(ft)%fforms,(2))
     HDR_fforms(ft)%fforms = (/1,2/)

   case (HDR_DENSITY)
     HDR_fforms(ft)%fname = "density"
     ABI_ALLOCATE(HDR_fforms(ft)%fforms,(2))
     HDR_fforms(ft)%fforms = (/51,52/)

   case (HDR_POTENTIAL)
     HDR_fforms(ft)%fname = "potential"
     ABI_ALLOCATE(HDR_fforms(ft)%fforms,(2))
     HDR_fforms(ft)%fforms = (/101,102/)

   case (HDR_WF_WVL)
     HDR_fforms(ft)%fname = "wf_wavelet"
     ABI_ALLOCATE(HDR_fforms(ft)%fforms,(1))
     HDR_fforms(ft)%fforms = (/202/)

   ! FIXME Screening part
   case (HDR_SCREENING)
     HDR_fforms(ft)%fname = "screening"
     ABI_ALLOCATE(HDR_fforms(ft)%fforms,(2))
     HDR_fforms(ft)%fforms = (/1002,1102/)

   case (HDR_BS_HAM)
     HDR_fforms(ft)%fname = "bs_hamiltonian"
     ABI_ALLOCATE(HDR_fforms(ft)%fforms,(1))
     HDR_fforms(ft)%fforms = (/301/)

   case (HDR_BS_EIG)
     HDR_fforms(ft)%fname = "bs_eigenstates"
     ABI_ALLOCATE(HDR_fforms(ft)%fforms,(1))
     HDR_fforms(ft)%fforms = (/401/)
    
   case (HDR_BS_HAYDOCK)
     HDR_fforms(ft)%fname = "bs_haydock"
     ABI_ALLOCATE(HDR_fforms(ft)%fforms,(1))
     HDR_fforms(ft)%fforms = (/501/)

   case default
     ierr=ierr+1
   end select

 end do

end subroutine m_header_init
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_fform2ftype
!! NAME
!! hdr_fform2ftype
!!
!! FUNCTION
!!  Return the filetype flag, a string with the name of the file and the patch level from the
!!  input fform.
!!
!! INPUTS
!!  fform=The value of fform read from the header.
!!
!! OUTPUT
!!  fname=The name of the file.
!!  ftype=The integer flag giving the filetype.
!!  patch_level=The 
!!  ierr=A nonzero value signals that fform is not allowed 
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine hdr_fform2ftype(fform,fname,ftype,patch_level,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_fform2ftype'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fform
 integer,intent(out) :: ftype,patch_level,ierr
 character(len=500),intent(out):: fname

!Local variables-------------------------------
!scalars
 integer :: ii,patch

! *************************************************************************

 !@fform_t
 ierr=0
 do ii=1,SIZE(HDR_fforms)
   if ( ANY( fform == HDR_fforms(ii)%fforms) ) then
     ierr = ierr+1
     ftype = ii
     fname = HDR_fforms(ii)%fname
     do patch=1,SIZE(HDR_fforms(ii)%fforms)
       if (fform == HDR_fforms(ii)%fforms(patch)) then 
         patch_level=patch-1
         EXIT
       end if
     end do
   end if
 end do

 if (ierr==0) then ! Invalid fform, raise the error.
   ierr=1
 else if (ierr>1) then
   MSG_ERROR("Internal database is inconsistent")
 else
  ierr=0
 end if

end subroutine hdr_fform2ftype
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_ftype2fform
!! NAME
!! hdr_ftype2fform
!!
!! FUNCTION
!!  Returns the last fform associated to the filetype ftype.
!!
!! INPUTS
!!  ftype
!! 
!! OUTPUT
!!  fform
!!  ierr
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine hdr_ftype2fform(ftype,fform,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_ftype2fform'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ftype
 integer,intent(out) :: fform,ierr

!Local variables-------------------------------
!scalars
 integer :: last

! *************************************************************************

 !@fform_t
 ierr=0
 if ( ftype<1 .or. ftype>SIZE(HDR_fforms) ) then ! Wrong ftype.
  ierr=1; RETURN
 end if
 !
 ! Return the last fform associated to this ftype.
 last = SIZE(HDR_fforms(ftype)%fforms)
 fform = HDR_fforms(ftype)%fforms(last) 

end subroutine hdr_ftype2fform
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_init
!! NAME
!! hdr_init
!!
!! FUNCTION
!! This subroutine initializes the header structured datatype
!! and most of its content from dtset and psps, and put default values for
!! evolving variables.
!!
!! INPUTS
!! bstruct <type(bandstructure_type)>=band structure information
!!  including Brillouin zone description
!! codvsn=code version
!! dtset <type(dataset_type)>=all input variables for this dataset
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! pertcase=index of the perturbation, or 0 if GS calculation
!! psps <type(pseudopotential_type)>=all the information about psps
!!
!! OUTPUT
!! hdr <type(hdr_type)>=the header, initialized, and for most part of
!!   it, contain its definite values, except for evolving variables
!!
!! PARENTS
!!      gstate,loper3,newsp,nonlinear,respfn,setup_bse,setup_screening
!!      setup_sigma
!!
!! CHILDREN
!!      mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine hdr_init(bstruct,codvsn,dtset,hdr,pawtab,pertcase,psps,wvl)

 use defs_basis
 use defs_wvltypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: pertcase
 character(len=6),intent(in) :: codvsn
 type(bandstructure_type),intent(in) :: bstruct
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(out) :: hdr
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type),intent(in) :: wvl
!arrays
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************

 !@hdr_type

! More checking would be needed ...
 if (dtset%ntypat/=psps%ntypat) then
   write(msg,'(a,2i0)')' dtset%ntypat and psps%ntypat differs. They are :',dtset%ntypat,psps%ntypat
   MSG_ERROR(msg)
 end if

 if (dtset%npsp/=psps%npsp) then
   write(msg,'(a,2i0)')' dtset%npsp and psps%npsp differs. They are :',dtset%npsp,psps%npsp
   MSG_ERROR(msg)
 end if

 call hdr_init_lowlvl(hdr,bstruct,psps,pawtab,wvl,&
&  codvsn,pertcase,dtset%natom,dtset%nsym,dtset%nspden,dtset%ecut,dtset%pawecutdg,dtset%ecutsm,dtset%dilatmx,&
&  dtset%intxc,dtset%ixc,dtset%stmbias,dtset%usewvl,dtset%pawcpxocc,dtset%ngfft,dtset%ngfftdg,dtset%so_psp,dtset%qptn,&
&  dtset%rprimd_orig(:,:,1),dtset%xred_orig(:,:,1),dtset%symrel,dtset%tnons,dtset%symafm,dtset%typat)

end subroutine hdr_init
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_nullify
!! NAME
!! hdr_nullify
!!
!! FUNCTION
!! This subroutine set all the pointers to null()
!!
!! INPUTS
!! hdr <type(hdr_type)>=the header
!!
!! PARENTS
!!      elphon,m_gamma,m_io_gkk,m_io_screening,read_el_veloc,read_gkk
!!
!! CHILDREN
!!      mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine hdr_nullify(hdr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(hdr_type),intent(inout) :: hdr

! *************************************************************************
 !@hdr_type

! integer
 nullify(hdr%istwfk)
 nullify(hdr%lmn_size)
 nullify(hdr%nband)
 nullify(hdr%npwarr)
 nullify(hdr%pspcod)
 nullify(hdr%pspdat)
 nullify(hdr%pspso)
 nullify(hdr%pspxc)
 nullify(hdr%so_psp)
 nullify(hdr%symafm)
 nullify(hdr%symrel)
 nullify(hdr%typat)

! real
 nullify(hdr%kptns)
 nullify(hdr%occ)
 nullify(hdr%tnons)
 nullify(hdr%wtk)
 nullify(hdr%xred)
 nullify(hdr%zionpsp)
 nullify(hdr%znuclpsp)
 nullify(hdr%znucltypat)

! string arrays
 nullify(hdr%title)

! types
 nullify(hdr%pawrhoij)

end subroutine hdr_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_clean
!! NAME
!! hdr_clean
!!
!! FUNCTION
!! This subroutine deallocates the components of the header structured datatype
!!
!! INPUTS
!! hdr <type(hdr_type)>=the header
!!
!! OUTPUT
!!  (only deallocate)
!!
!! PARENTS
!!      bethe_salpeter,compare_interpol,conducti_nc,conducti_paw
!!      conducti_paw_core,cut3d,elphon,emispec_paw,gstate,gw_tools,initaim
!!      inpgkk,inwffil,ioarr,kss2wfk,linear_optics_paw,loper3,m_bse_io,m_ebands
!!      m_gamma,m_gwannier,m_io_gkk,m_io_kss,m_io_screening,m_wannier2abinit
!!      m_wfs,macroave,mrggkk,newsp,nonlinear,optic,rdm,read_el_veloc,read_gkk
!!      respfn,screening,sigma,suscep
!!
!! CHILDREN
!!      mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine hdr_clean(hdr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_clean'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------

! *************************************************************************

 DBG_ENTER("COLL")

 !@hdr_type

 !integer
 if (associated(hdr%istwfk  ))   then
   ABI_DEALLOCATE(hdr%istwfk)
 end if
 if (associated(hdr%lmn_size))   then
   ABI_DEALLOCATE(hdr%lmn_size)
 end if
 if (associated(hdr%nband   ))   then
   ABI_DEALLOCATE(hdr%nband)
 end if
 if (associated(hdr%npwarr  ))   then
   ABI_DEALLOCATE(hdr%npwarr)
 end if
                            
 if (associated(hdr%pspcod))   then
   ABI_DEALLOCATE(hdr%pspcod)
 end if
 if (associated(hdr%pspdat))   then
   ABI_DEALLOCATE(hdr%pspdat)
 end if
 if (associated(hdr%pspso ))   then
   ABI_DEALLOCATE(hdr%pspso)
 end if
 if (associated(hdr%pspxc ))   then
   ABI_DEALLOCATE(hdr%pspxc)
 end if
 if (associated(hdr%so_psp))   then
   ABI_DEALLOCATE(hdr%so_psp)
 end if
 if (associated(hdr%symafm))   then
   ABI_DEALLOCATE(hdr%symafm)
 end if
 if (associated(hdr%symrel))   then
   ABI_DEALLOCATE(hdr%symrel)
 end if
 if (associated(hdr%typat ))   then
   ABI_DEALLOCATE(hdr%typat)
 end if
                            
 !real
 if (associated(hdr%kptns     ))   then
   ABI_DEALLOCATE(hdr%kptns)
 end if
 if (associated(hdr%occ       ))   then
   ABI_DEALLOCATE(hdr%occ)
 end if
 if (associated(hdr%tnons     ))   then
   ABI_DEALLOCATE(hdr%tnons)
 end if
 if (associated(hdr%wtk       ))   then
   ABI_DEALLOCATE(hdr%wtk)
 end if
 if (associated(hdr%xred      ))   then
   ABI_DEALLOCATE(hdr%xred)
 end if
 if (associated(hdr%zionpsp   ))   then
   ABI_DEALLOCATE(hdr%zionpsp)
 end if
 if (associated(hdr%znuclpsp  ))   then
   ABI_DEALLOCATE(hdr%znuclpsp)
 end if
 if (associated(hdr%znucltypat))   then
   ABI_DEALLOCATE(hdr%znucltypat)
 end if
                            
 !string arrays
 if(associated(hdr%title))   then
   ABI_DEALLOCATE(hdr%title)
 end if

 if (hdr%usepaw==1 .and. associated(hdr%pawrhoij) ) then
  call rhoij_free(hdr%pawrhoij)
  ABI_DEALLOCATE(hdr%pawrhoij)
 end if

 DBG_EXIT("COLL")

end subroutine hdr_clean
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_copy
!! NAME
!! hdr_copy
!!
!! FUNCTION
!! Make a deep copy of the abinit header.
!!
!! INPUTS
!!  Hdr_in=The header to be copied.
!!
!! OUTPUT
!!  Hdr_cp=The deep copy of Hdr_in.
!!
!! NOTES
!!  The present version deals with versions of the header up to 56.
!!
!! PARENTS
!!      m_io_kss,m_io_screening,m_wannier2abinit
!!
!! CHILDREN
!!      mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine hdr_copy(Hdr_in,Hdr_cp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_copy'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(hdr_type),intent(in) :: Hdr_in
 type(hdr_type),intent(inout) :: Hdr_cp

!Local variables-------------------------------
!scalars
 integer :: cplex
 character(len=500) :: msg
!arrays
 integer,pointer :: nlmn(:)
! *************************************************************************

!TODO add method to nullify header, deep_copy might crash if
!a pointer in the Hdr is not initialized to null(), tipically
!this happens for PAW quantities but we always check Hdr%usepaw

 !@hdr_type
 if (Hdr_in%headform>57) then
  write(msg,'(3a,i3,2a)')&
&  ' hdr_copy deals with versions of the header only up to 57. ',ch10,&
&  ' However headform = ',Hdr_in%headform,ch10,&
&  ' Change the source to add the changes done in the new version. '
  MSG_ERROR(msg)
 end if

!=== Integer values ===
 Hdr_cp%bantot   = Hdr_in%bantot
 Hdr_cp%date     = Hdr_in%date
 Hdr_cp%headform = Hdr_in%headform
 Hdr_cp%intxc    = Hdr_in%intxc
 Hdr_cp%ixc      = Hdr_in%ixc
 Hdr_cp%natom    = Hdr_in%natom
 Hdr_cp%nkpt     = Hdr_in%nkpt
 Hdr_cp%npsp     = Hdr_in%npsp
 Hdr_cp%nspden   = Hdr_in%nspden
 Hdr_cp%nspinor  = Hdr_in%nspinor
 Hdr_cp%nsppol   = Hdr_in%nsppol
 Hdr_cp%nsym     = Hdr_in%nsym
 Hdr_cp%ntypat   = Hdr_in%ntypat
 Hdr_cp%occopt   = Hdr_in%occopt
 Hdr_cp%pertcase = Hdr_in%pertcase
 Hdr_cp%usepaw   = Hdr_in%usepaw
 Hdr_cp%usewvl   = Hdr_in%usewvl

 ! === Integer arrays ===
 Hdr_cp%ngfft(:)   = Hdr_in%ngfft(:)
 Hdr_cp%nwvlarr(:) = Hdr_in%nwvlarr(:)

!=== Integer pointers ====
 call deep_copy( Hdr_in%istwfk,  Hdr_cp%istwfk   )
 call deep_copy( Hdr_in%lmn_size,Hdr_cp%lmn_size )
 call deep_copy( Hdr_in%nband,   Hdr_cp%nband    )
 call deep_copy( Hdr_in%npwarr,  Hdr_cp%npwarr   )
 call deep_copy( Hdr_in%pspcod,  Hdr_cp%pspcod )
 call deep_copy( Hdr_in%pspdat,  Hdr_cp%pspdat )
 call deep_copy( Hdr_in%pspso ,  Hdr_cp%pspso  )
 call deep_copy( Hdr_in%pspxc ,  Hdr_cp%pspxc  )
 call deep_copy( Hdr_in%so_psp,  Hdr_cp%so_psp )
 call deep_copy( Hdr_in%symafm,  Hdr_cp%symafm )
 call deep_copy( Hdr_in%symrel,  Hdr_cp%symrel )
 call deep_copy( Hdr_in%typat ,  Hdr_cp%typat  )

!=== Real variables ====
 Hdr_cp%ecut        = Hdr_in%ecut
 Hdr_cp%ecutdg      = Hdr_in%ecutdg
 Hdr_cp%ecutsm      = Hdr_in%ecutsm
 Hdr_cp%ecut_eff    = Hdr_in%ecut_eff
 Hdr_cp%etot        = Hdr_in%etot
 Hdr_cp%fermie      = Hdr_in%fermie
 Hdr_cp%residm      = Hdr_in%residm
 Hdr_cp%stmbias     = Hdr_in%stmbias
 Hdr_cp%tphysel     = Hdr_in%tphysel
 Hdr_cp%tsmear      = Hdr_in%tsmear

 Hdr_cp%qptn(:)     = Hdr_in%qptn(:)
 Hdr_cp%rprimd(:,:) = Hdr_in%rprimd(:,:)

!=== Real pointers ===
 call deep_copy( Hdr_in%kptns     ,Hdr_cp%kptns     )
 call deep_copy( Hdr_in%occ       ,Hdr_cp%occ       )
 !write(std_out,*)"DEBUG: ",Hdr_in%occ
 !write(std_out,*)"DEBUG: ",Hdr_in%bantot
 !if (associated(Hdr_cp%occ)) write(std_out,*)"DEBUG: ",Hdr_cp%occ
 call deep_copy( Hdr_in%tnons     ,Hdr_cp%tnons     )
 call deep_copy( Hdr_in%wtk       ,Hdr_cp%wtk       )
 call deep_copy( Hdr_in%xred      ,Hdr_cp%xred      )
 call deep_copy( Hdr_in%zionpsp   ,Hdr_cp%zionpsp   )
 call deep_copy( Hdr_in%znuclpsp  ,Hdr_cp%znuclpsp  )
 call deep_copy( Hdr_in%znucltypat,Hdr_cp%znucltypat)

!=== Character pointers ===
 Hdr_cp%codvsn = Hdr_in%codvsn
! THIS DOES NOT WORK ON XLF: Hdr_cp%title string length becomes huge and segfaults
! call deep_copy( Hdr_in%title,Hdr_cp%title )
 ABI_ALLOCATE(Hdr_cp%title,(Hdr_cp%npsp))
 Hdr_cp%title = Hdr_in%title

!=== For PAW have to copy Pawrhoij ====
!* TODO alchemy requires a different treatment but for the moment it is not available within PAW.
 if (Hdr_in%usepaw==1) then
  cplex = Hdr_in%Pawrhoij(1)%cplex
  nlmn => Hdr_in%lmn_size(1:Hdr_in%ntypat)
  ABI_ALLOCATE(Hdr_cp%Pawrhoij,(Hdr_in%natom))
  call rhoij_alloc(cplex,nlmn,Hdr_in%nspden,Hdr_in%nspinor,Hdr_in%nsppol,Hdr_cp%Pawrhoij,Hdr_in%typat)
  call rhoij_copy(Hdr_in%Pawrhoij,Hdr_cp%Pawrhoij)
 end if

end subroutine hdr_copy
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_get_nelect_byocc
!! NAME
!! hdr_get_nelect_byocc
!!
!! FUNCTION
!!  Return the number of electrons from the occupation numbers
!!  thus taking into account a possible additional charge or alchemy.
!!
!! INPUTS
!!  Hdr<hdr_type>
!!
!! OUTPUT
!!  nelect=Number of electrons in the unit cell.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function hdr_get_nelect_byocc(Hdr) result(nelect)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_get_nelect_byocc'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 type(hdr_type),intent(in) :: Hdr
 real(dp) :: nelect

!Local variables ---------------------------------------
!scalars
 integer :: idx,isppol,ikibz,nband_k

! *************************************************************************

!* Cannot use znucl because we might have additional charge or alchemy.
 nelect=zero ; idx=0
 do isppol=1,Hdr%nsppol
   do ikibz=1,Hdr%nkpt
     nband_k=Hdr%nband(ikibz+(isppol-1)*Hdr%nkpt)
     nelect = nelect + Hdr%wtk(ikibz)*SUM(Hdr%occ(idx+1:idx+nband_k))
     idx=idx+nband_k
   end do
 end do

!Might also check also Hdr%znuclpsp(:) to avoid round off errors

end function hdr_get_nelect_byocc
!!***

!----------------------------------------------------------------------

!!****f* m_header/isknown_headform
!! NAME
!!  isknown_headform
!!
!! FUNCTION
!!  Returns .TRUE. if headform is one of the allowed values.
!!
!! INPUTS
!!  headform
!!
!! SOURCE

function isknown_headform(headform) result(ans)

 use defs_basis

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isknown_headform'
!End of the abilint section

 integer,intent(in) :: headform
 logical :: ans

! *************************************************************************

 ans = ANY(headform == HDR_KNOWN_HEADFORMS)

end function isknown_headform
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_init_lowlvl
!! NAME
!! hdr_init_lowlvl
!!
!! FUNCTION
!! This subroutine initializes the header structured datatype
!! and most of its content from psps and other input variables that 
!! are passed explicitly. It also use default values for evolving variables. 
!! Note that Dtset is not required thus rendering the initialization of the header
!! much easier.
!!
!! INPUTS
!! bstruct <type(bandstructure_type)>=band structure information
!!  including Brillouin zone description
!! codvsn=code version
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! pertcase=index of the perturbation, or 0 if GS calculation
!! psps <type(pseudopotential_type)>=all the information about psps
!! For the meaning of the other varialble see the definition of dataset_type.
!!
!! OUTPUT
!! hdr <type(hdr_type)>=the header, initialized, and for most part of
!!   it, contain its definite values, except for evolving variables
!!
!! PARENTS
!!      m_header,m_wfs
!!
!! CHILDREN
!!      mpi_file_read_at,xmpio_read_frm
!!
!! SOURCE

subroutine hdr_init_lowlvl(hdr,bstruct,psps,pawtab,wvl,&
&  codvsn,pertcase,natom,nsym,nspden,ecut,pawecutdg,ecutsm,dilatmx,&
&  intxc,ixc,stmbias,usewvl,pawcpxocc,ngfft,ngfftdg,so_psp,qptn,&
&  rprimd,xred,symrel,tnons,symafm,typat)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_init_lowlvl'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym,nspden,intxc,ixc,usewvl,pawcpxocc,pertcase
 real(dp),intent(in) :: ecut,ecutsm,dilatmx,stmbias,pawecutdg
 character(len=6),intent(in) :: codvsn
 type(bandstructure_type),intent(in) :: bstruct
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type),intent(in) :: wvl
 type(hdr_type),intent(out) :: hdr
!arrays
 integer,intent(in) :: typat(natom)
 integer,intent(in) ::  so_psp(psps%npsp)
 integer,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 integer,intent(in) :: ngfft(18),ngfftdg(18)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(in) :: qptn(3) ! the wavevector, in case of a perturbation
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: bantot,date,itypat,nkpt,npsp,ntypat,nsppol,nspinor
 integer :: idx,isppol,ikpt,iband,ipsp
#ifndef HAVE_DFT_BIGDFT
 character(len=500) :: msg
#endif
 character(len=8) :: date_time
!arrays
 integer,allocatable :: nlmn(:)

! *************************************************************************

 !@hdr_type
 call date_and_time(date_time)
 read(date_time,'(i8)')date

 npsp   = psps%npsp
 ntypat = psps%ntypat
 nkpt   = bstruct%nkpt
 nsppol = bstruct%nsppol
 nspinor= bstruct%nspinor
 bantot = bstruct%bantot

!Transfer dimensions and other scalars to hdr.
 hdr%intxc    =intxc
 hdr%ixc      =ixc
 hdr%natom    =natom
 hdr%npsp     =npsp
 hdr%nspden   =nspden
 hdr%nspinor  =nspinor
 hdr%nsym     =nsym
 hdr%ntypat   =ntypat
 hdr%bantot   =bantot
 hdr%nkpt     =nkpt
 hdr%nsppol   =nsppol
 hdr%usepaw   =psps%usepaw
 hdr%usewvl   =usewvl !hdr%nwvlarr will be set later since the number !of wavelets have not yet been computed.
 hdr%occopt   =bstruct%occopt
 hdr%codvsn   =codvsn
 hdr%date     =date
 hdr%headform =HDR_LATEST_HEADFORM ! Initialize with the latest headform
 hdr%pertcase =pertcase
 hdr%ecut     =ecut
 hdr%ecutsm   =ecutsm
 hdr%ecut_eff =ecut * (dilatmx)**2
 hdr%stmbias  =stmbias
 hdr%tphysel  =bstruct%tphysel
 hdr%tsmear   =bstruct%tsmear
 hdr%qptn     =qptn
 hdr%rprimd   =rprimd      ! Evolving data

!Default for other data  (all evolving data)
 hdr%etot     =1.0d20
 hdr%fermie   =1.0d20
 hdr%residm   =1.0d20

!Allocate all components of hdr

!Transfer data from bstruct
 ABI_ALLOCATE(hdr%istwfk,(nkpt))
 hdr%istwfk(1:nkpt)      =bstruct%istwfk(1:nkpt)
 ABI_ALLOCATE(hdr%kptns,(3,nkpt))
 hdr%kptns(:,:)          =bstruct%kptns(:,:)
 ABI_ALLOCATE(hdr%nband,(nkpt*nsppol))
 hdr%nband(1:nkpt*nsppol)=bstruct%nband(1:nkpt*nsppol)
 ABI_ALLOCATE(hdr%npwarr,(nkpt))
 hdr%npwarr(:)           =bstruct%npwarr(:)
 ABI_ALLOCATE(hdr%wtk,(nkpt))
 hdr%wtk(:)=bstruct%wtk(:)

!Transfer data from psps
 ABI_ALLOCATE(hdr%pspcod,(npsp))
 hdr%pspcod    =psps%pspcod
 ABI_ALLOCATE(hdr%pspdat,(npsp))
 hdr%pspdat    =psps%pspdat         
 ABI_ALLOCATE(hdr%pspso,(npsp))
 hdr%pspso     =psps%pspso          
 ABI_ALLOCATE(hdr%pspxc,(npsp))
 hdr%pspxc     =psps%pspxc          
 ABI_ALLOCATE(hdr%znuclpsp,(npsp))
 hdr%znuclpsp  =psps%znuclpsp
 ABI_ALLOCATE(hdr%znucltypat,(ntypat))
 hdr%znucltypat=psps%znucltypat
 ABI_ALLOCATE(hdr%zionpsp,(npsp))
 hdr%zionpsp   =psps%zionpsp     
 ABI_ALLOCATE(hdr%title,(npsp))
 do ipsp=1,psps%npsp
   write(hdr%title(ipsp), "(A)") psps%title(ipsp)(1:132)
 end do

 ABI_ALLOCATE(hdr%so_psp,(npsp))
 hdr%so_psp=so_psp

 ABI_ALLOCATE(hdr%symafm,(nsym))
 hdr%symafm(1:min(size(symafm),size(hdr%symafm)))=symafm(1:min(size(symafm),size(hdr%symafm)))

 ABI_ALLOCATE(hdr%symrel,(3,3,nsym))
 hdr%symrel(:,:,1:min(size(symrel,3),size(hdr%symrel,3))) =symrel(:,:,1:min(size(symrel,3),size(hdr%symrel,3)))

 ABI_ALLOCATE(hdr%tnons,(3,nsym))
 hdr%tnons(:,1:min(size(tnons,2),size(hdr%tnons,2)))=tnons(:,1:min(size(tnons,2),size(hdr%tnons,2)))

 ABI_ALLOCATE(hdr%typat,(natom))
 hdr%typat(1:natom) =typat(1:natom)  ! PMA : in tests/v2/t11 size(dtset%typat) is bigger dtset%natom
 ABI_ALLOCATE(hdr%xred,(3,natom))
 hdr%xred(:,1:natom)=xred(:,1:natom) ! Evolving data

 if(psps%usepaw==1)then
   ABI_ALLOCATE(hdr%pawrhoij,(natom))
   ABI_ALLOCATE(nlmn,(ntypat))
   do itypat=1,ntypat
     nlmn(itypat)=pawtab(itypat)%lmn_size
   end do
   !Values of nspden/nspinor/nsppol are dummy ones; they are overwritten later (by hdr_update)
   call rhoij_alloc(pawcpxocc,nlmn,nspden,nspinor,nsppol,hdr%pawrhoij,typat)
   ABI_DEALLOCATE(nlmn)
 end if

 if (psps%usepaw==1) then
   hdr%ngfft(:) =ngfftdg(1:3)
 else if (usewvl==1) then
#if defined HAVE_DFT_BIGDFT
   hdr%ngfft(:) = (/ wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i /)
#else
 write(msg, '(a,a,a,a)' ) ch10,&
& ' hdr_init_lowlvl : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(std_out,msg,'COLL')
 call leave_new('COLL')
#endif
 else
   hdr%ngfft(:) =ngfft(1:3)
 end if

!Transfer paw data
 ABI_ALLOCATE(hdr%lmn_size,(npsp))
 if(psps%usepaw==1) then
   hdr%ecutdg   =pawecutdg
   hdr%lmn_size(1:npsp)=pawtab(1:npsp)%lmn_size
 else
   hdr%ecutdg=hdr%ecut
   hdr%lmn_size(:)=psps%lmnmax
 end if

 !call get_eneocc_vect(bstruct,'occ',hdr%occ)  ! Evolving data
 ABI_ALLOCATE(hdr%occ,(bantot))
 hdr%occ(:)=zero; idx=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do iband=1,hdr%nband(ikpt+(isppol-1)*nkpt)
       idx=idx+1
       hdr%occ(idx)=bstruct%occ(iband,ikpt,isppol)
     end do
   end do
 end do

end subroutine hdr_init_lowlvl
!!***

!----------------------------------------------------------------------

!!****f* m_header/hdr_mpio_skip
!! NAME
!!  hdr_mio_skip
!!
!! FUNCTION
!!   Skip the abinit header in MPI-IO mode. This routine uses local MPI-IO calls hence 
!!   it can be safely called by master node only. Note however that in this case the
!!   offset has to be communicated to the other nodes.
!!
!! INPUTS
!!  mpio_fh=MPI-IO file handler
!!
!! OUTPUT
!!  fform=kind of the array in the file
!!  offset=The offset of the Fortran record located immediately below the Abinit header.
!!
!! SOURCE

subroutine hdr_mpio_skip(mpio_fh,fform,offset) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_mpio_skip'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: mpio_fh
 integer,intent(out) :: fform
 integer(kind=XMPI_OFFSET_KIND),intent(out) :: offset

!Local variables-------------------------------
!scalars
 integer :: bsize_frm,mpi_type_frm
#ifdef HAVE_MPI_IO
 integer :: headform,ierr,mu,usepaw,npsp
!arrays
 integer(kind=MPI_OFFSET_KIND) :: fmarker,positloc
 integer :: statux(MPI_STATUS_SIZE)
#endif

! *************************************************************************

 offset = 0
 fform  = 0

 bsize_frm    = xmpio_bsize_frm    ! bsize_frm= Byte length of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

#ifdef HAVE_MPI_IO
!Reading the first record of the file -------------------------------------
!read (unitfi)   codvsn,headform,..............
 positloc = bsize_frm + 6*xmpi_bsize_ch
 call MPI_FILE_READ_AT(mpio_fh,positloc,fform,1,MPI_INTEGER,statux,ierr)
 !write(std_out,*)" MPI_FILE_READ_AT test fform: ",fform

 if (ANY(fform == (/1,2,51,52,101,102/) )) then 
   headform=22  ! This is the old format !read (unitfi) codvsn,fform
 else 
   !read (unitfi)codvsn,headform,fform
   call MPI_FILE_READ_AT(mpio_fh,positloc,headform,1,MPI_INTEGER,statux,ierr)
   positloc = positloc + xmpi_bsize_int
   call MPI_FILE_READ_AT(mpio_fh,positloc,fform,1,MPI_INTEGER,statux,ierr)
   !write(std_out,*)" MPI_FILE_READ_AT: headform, fform: ",headform,fform
 end if

 ! Skip first record.
 call xmpio_read_frm(mpio_fh,offset,xmpio_at,fmarker,ierr)

!Reading the second record of the file: read(unitfi) bantot, hdr%date, hdr%intxc.................
!Read npsp and usepaw.
 positloc  = offset + bsize_frm + xmpi_bsize_int*13
 call MPI_FILE_READ_AT(mpio_fh,positloc,npsp,1,MPI_INTEGER,statux,ierr)

 usepaw=0
 if (headform >= 44) then
   positloc = positloc +  xmpi_bsize_int*4
   call MPI_FILE_READ_AT(mpio_fh,positloc,usepaw,1,MPI_INTEGER,statux,ierr)
 end if

 ! Skip second record.
 call xmpio_read_frm(mpio_fh,offset,xmpio_at,fmarker,ierr)

 ! Skip the rest of the file ---------------------------------------------
 do mu=1,2+npsp
   call xmpio_read_frm(mpio_fh,offset,xmpio_at,fmarker,ierr)
 end do

 if (headform>=44.and.usepaw==1) then ! skip rhoij records.
   call xmpio_read_frm(mpio_fh,offset,xmpio_at,fmarker,ierr)
   call xmpio_read_frm(mpio_fh,offset,xmpio_at,fmarker,ierr)
 end if

#else
 MSG_ERROR("hdr_mpio_skip cannot be used when MPI-IO is not enabled")
 ABI_UNUSED(mpio_fh) 
#endif

end subroutine hdr_mpio_skip
!!***

!----------------------------------------------------------------------

END MODULE m_header
!!***
