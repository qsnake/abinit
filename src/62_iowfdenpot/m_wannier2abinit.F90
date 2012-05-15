!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_wannier2abinit
!! NAME
!!  m_wannier2abinit
!!
!! FUNCTION
!!  This module defines an interface to perform the interpolation of k-dependent
!!  matrix elements using a representation in terms of Wannier functions.
!!  The unitary transformation defining the Wannier gauge is supposed to be
!!  stored in the WAN file produced by Abinit after having performed the wannierization 
!!  by calling wannier90 in library-mode.
!!  This module provides a procedure to read the WAN file initializing a WannierData object.
!!  The object, in turn, provides a high-level interface that can be used to perform 
!!  the interpolation of user-provided matrix elements.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_wannier2abinit

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

 use m_fstrings,     only : toupper
 use m_geometry,     only : wigner_seitz, normv
 use m_io_tools,     only : get_unit
 use m_header,       only : hdr_copy, hdr_clean
 use m_bz_mesh,      only : bz_mesh_type

 implicit none

 private

!!***

!!****t* m_wannier2abinit/WannierData
!! NAME
!! WannierData
!!
!! FUNCTION
!! Object used to store and handle Wannier90 results inside Abinit
!!
!! SOURCE

 ! TODO a better integration of Wannier90 and improvements of the library-mode

 type,public ::  WannierData

  integer :: mband            ! Total number of bands to be processed.
  integer :: mwan             ! Max number of Wannier functions over spin, i.e MAXVAL(nwan) (to dimension arrays).
  integer :: nntot            ! Number of k-point neighbour.
  integer :: nkpt             ! Number of k-points.
  integer :: nsppol           ! Number of independent spin polarizations (presently only nsppol=1 is implemented).
  integer :: WDversion        ! Version of the WannierData file.

  real(dp) :: W90version      ! Wannier90 version.

  character(len=800) :: title

  !from W90
  !if (have_disentangled) then
   ! write(chk_unit) omega_invariant     ! Omega invariant
   ! lwindow, ndimwin and U_matrix_opt 
   !OK write(chk_unit) ((lwindow(i,nkp),i=1,num_bands),nkp=1,WanData%nkpt)
   !?  write(chk_unit) (ndimwin(nkp),nkp=1,WanData%nkpt)
   !OK write(chk_unit) (((u_matrix_opt(i,j,nkp),i=1,num_bands),j=1,num_wann),nkp=1,WanData%nkpt)
  !endif

!Arrays
  integer,pointer :: nwan(:)        SET2NULL
   ! nwan(nsppol)
   ! Number of wannier functions (read in wannier90.win).

  real(dp),pointer :: eigen(:,:,:)  SET2NULL
   ! eigen(mband,nkpt,nsppol)

!TODO convert Everything to Bohr to be consistent with Abinit internal conventions.
! inside the creation method
  real(dp),pointer :: spreadw(:,:)  SET2NULL
   ! spreadw(3,nsppol)

  real(dp),pointer :: wann_centres(:,:,:)   SET2NULL
   ! wann_centres(3,mwan,nsppol)

  real(dp),pointer :: wann_spreads(:,:)   SET2NULL
   ! wann_spreads(mwan,nsppol)

  complex(dp),pointer :: U_matrix(:,:,:,:)   SET2NULL
   ! U_matrix(mwan,mwan,nkpt,nsppol)

  complex(dp),pointer :: U_matrix_opt(:,:,:,:)   SET2NULL
   ! U_matrix_opt(mband,mwan,nkpt,nsppol) 

  logical,pointer :: band_in(:,:)   SET2NULL
   ! band_in(mband,nsppol)
   ! .TRUE. if the band is included in the calculation.   

  logical,pointer :: lwindow(:,:,:)  SET2NULL
   ! lwindow(mband,nkpt,nsppol)
   ! Only if disentanglement, .TRUE. if this band at this k-point lies within the outer window

  logical,pointer :: have_disentangled(:)   SET2NULL
   ! have_disentangled(nsppol)
   ! Whether a disentanglement has been performed

  character(len=100) :: cut_mode='None'
   ! Wheter the Hamiltonian in real space in the Wannier gauge has to be truncated.

  type(Hdr_type) :: Hdr 
   ! The abinit header

  ! ==================================================================
  ! Variable and arrays used to perform the Wannier interpolation
  ! * NB: The quantities below are initialized and defined during the 
  !  interpolation, They are _not_ stored in the WAN file but are 
  !  calculated starting from the Wannierization results.
  ! ==================================================================

  integer :: nrpts
  ! Number of points in the Wigner-Seitz cell

  integer,pointer :: irvec(:,:)   SET2NULL
  ! irvec(3,nrpts) 
  ! Lattice vectors in the WS cell in the basis of the lattice vectors 
  ! defining the unit cell

  integer,pointer :: ndegen(:)   SET2NULL
  ! ndegen(nrpts)
  ! Degeneracy of each point. It will be weighted using 1/ndegen(ii)

  integer,pointer :: ndimwin(:,:)  SET2NULL
  ! ndimwin(nkpt,nsppol)
  ! Number of bands inside outer window at nkpt-th k point

  complex(dpc),pointer :: hamWR(:,:,:,:)   SET2NULL
  ! hamWR(mwan,mwan,nrpts,nsppol))
  ! Hamiltonian in k-space (ab-initio grid) in the Wannier gauge.
 end type WannierData

! === List of available public routines and functions ===
 public ::  InitWanData       ! Creation method.
 public ::  DestroyWanData    ! Destruction method.
 public ::  DumpWanData       ! Dump the object on file.
 public ::  ReadWanData       ! Read the object from file.
 public ::  MakeWannierHR     ! Construct matrix elements of H in the Wannier representation.
 public ::  PlotWannierHR     ! Plot the above matrix elements.
 public ::  PrintWanData      ! Printout of the content of the object
 public ::  WanMatInterpol    ! Public interface to interpolate matrix elements.
 public ::  WannierInterpol   ! Perform the Wannier interpolation of the energies stored in the object.
 public ::  WrapInitWanData   ! Wrapper for the creation method, needed to bypass the problem with nsppol=2

 interface WanMatInterpol
  module procedure WanMatInterpol_dpc
  module procedure WanMatInterpol_dp
 end interface 

!----------------------------------------------------------------------

 integer,private,parameter :: LATEST_WD_VERSION=0

 integer,private,parameter :: LMAX(3)=(/2,2,2/)
 real(dp),private,parameter :: gamma_point(3)=(/zero,zero,zero/)

CONTAINS  !===========================================================
!!***

!!****f* m_wannier2abinit/InitWanData
!! NAME
!! InitWanData
!!
!! FUNCTION
!!  Creation method for the WannierData object.
!!  The most important results of wannier90 are passed in input to the routine which
!!  allocates and initializes the structure.
!!
!! INPUTS
!! nkpt=Number of k-points.
!! mband=Maximum number of bands.
!! nntot=Number of k-point neighbour
!! nsppol=Number of independent spin polarizations.
!! mwan=Max number of wannier functions (over spin)
!! W90version=Version of the Wannier90 code.
!! Hdr<Hdr_type>=The abinit header.
!! title=Informative title.
!! nwan(nsppol)=Number of Wannier functions for each spin.
!! spreadw(3,nsppol)=
!! wann_centres(3,mwan,nsppol)=Centres of the wannier functions (output of Wannier90)
!! wann_spreads(mwan,nsppol)=Spread of the wannier functions (output of Wannier90)
!! eigen(mband*nkpt*nsppol)=Kohn-Sham Electronic eigenvalues.
!! U_matrix(mwan,mwan,nkpt,nsppol)=Matrix defining the Wannier gauge.
!! U_matrix_opt(:,:,:,:)=Matrix defining the optimal subspace. Used only if we have used the disentanglement procedure.
!! have_disentangled(nsppol)=Defines whether a disentanglement has been done (for each spin, separately).
!! band_in(mband,nsppol)=Number of bands used to construct the Wannier functions.
!! lwindow(mband,nkpt,nsppol)=TRUE if the state falls within the energy window used for the disentanglement.
!!
!! OUTPUT
!! WanData<WannierData>=The object containing all data needed to perform the interpolation.
!!
!! PARENTS
!!      m_wannier2abinit
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine InitWanData(WanData,mband,nkpt,nsppol,nntot,mwan,nwan,have_disentangled,eigen,&
& band_in,lwindow,wann_centres,wann_spreads,spreadw,U_matrix,U_matrix_opt,Hdr,W90version,title)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'InitWanData'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,mband,nntot,nsppol,mwan
 real(dp),intent(in) :: W90version
 type(Hdr_type),intent(in) :: Hdr
 type(WannierData),intent(out) :: WanData
 character(len=*),intent(in) :: title
!arrays
 integer,intent(in) :: nwan(nsppol)
 real(dp),intent(in) :: spreadw(3,nsppol),wann_centres(3,mwan,nsppol),wann_spreads(mwan,nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 complex(dpc),intent(in) :: U_matrix(mwan,mwan,nkpt,nsppol)
 !complex(dpc),intent(in) :: U_matrix_opt(nwan,nwan,nkpt)
 complex(dpc),intent(in) :: U_matrix_opt(:,:,:,:)
 logical,intent(in) :: have_disentangled(nsppol)
 logical,intent(in) :: band_in(mband,nsppol)
 logical,intent(in) :: lwindow(mband,nkpt,nsppol)

!Local variables---------------------------
!scalars
 integer :: ikpt,isp,nband_k,idx,iband
! *************************************************************************

 call NullifyWanData_(WanData)

 WanData%mband     =mband
 WanData%mwan      =mwan
 WanData%nntot     =nntot
 WanData%nkpt      =nkpt
 WanData%nsppol    =nsppol
 WanData%W90version=W90version
 WanData%WDversion =LATEST_WD_VERSION

 WanData%title=title

 ! === Arrays ===
 ABI_ALLOCATE(WanData%nwan,(nsppol))
 WanData%nwan=nwan

 ABI_ALLOCATE(WanData%eigen,(mband,nkpt,nsppol))
 idx=0
 do isp=1,nsppol
  do ikpt=1,nkpt
   nband_k=Hdr%nband(ikpt+(isp-1)*nkpt)
   do iband=1,nband_k
    idx=idx+1
    WanData%eigen(iband,ikpt,isp)=eigen(idx)
   end do
  end do
 end do

 ABI_ALLOCATE(WanData%spreadw,(3,nsppol))
 WanData%spreadw=spreadw

 ABI_ALLOCATE(WanData%wann_centres,(3,mwan,nsppol))
 WanData%wann_centres=wann_centres

 ABI_ALLOCATE(WanData%wann_spreads,(mwan,nsppol))
 WanData%wann_spreads=wann_spreads

 ABI_ALLOCATE(WanData%U_matrix,(mwan,mwan,nkpt,nsppol))
 WanData%U_matrix=U_matrix

 ABI_ALLOCATE(WanData%have_disentangled,(nsppol))
 WanData%have_disentangled=have_disentangled

 ABI_ALLOCATE(WanData%lwindow,(mband,nkpt,nsppol))
 WanData%lwindow=lwindow(:,:,:)

 if (ANY(have_disentangled)) then 
  ABI_CHECK(mband==SIZE(U_matrix_opt,DIM=1),'Wrong mband size')
  ABI_CHECK(mwan ==SIZE(U_matrix_opt,DIM=2),'Wrong nwan size')
  ABI_CHECK(nkpt ==SIZE(U_matrix_opt,DIM=3),'Wrong nkpt size')
  ABI_CHECK(nsppol==SIZE(U_matrix_opt,DIM=4),'Wrong nsppol size')
  ABI_ALLOCATE(WanData%ndimwin,(nkpt,nsppol))
  do isp=1,nsppol
   do ikpt=1,nkpt
    WanData%ndimwin(ikpt,isp)=COUNT(WanData%lwindow(:,ikpt,isp))
   end do
  end do
  !TODO initialize U_matrix_opt if one of the spin has not been disentangled
  ABI_ALLOCATE(WanData%U_matrix_opt,(mband,mwan,nkpt,nsppol))
  WanData%U_matrix_opt=U_matrix_opt
 end if

 ABI_ALLOCATE(WanData%band_in,(mband,nsppol))
 WanData%band_in=band_in

 call hdr_copy(Hdr,WanData%Hdr)

end subroutine InitWanData
!!***

!!****f* m_wannier2abinit/WrapInitWanData
!! NAME
!! WrapInitWanData
!!
!! FUNCTION
!!  Just a wrapper for the real creation method to avoid problems with the extra dimension
!!  allocated for the spin. Wannier90, for the moment, accepts only a single spin for run!
!!  The WannierData, instead, takes into account this possibility.
!!
!! INPUTS
!! nkpt=Number of k-points.
!! mband=Maximum number of bands.
!! nntot=Number of k-point neighbour
!! nsppol=Number of independent spin polarizations.
!! mwan=Max number of wannier functions (over spin)
!! W90version=Version of the Wannier90 code.
!! Hdr<Hdr_type>=The abinit header.
!! title=Informative title.
!! nwan(nsppol)=Number of Wannier functions for each spin.
!! spreadw(3,nsppol)=
!! wann_centres(3,mwan,nsppol)=Centres of the wannier functions (output of Wannier90)
!! wann_spreads(mwan,nsppol)=Spread of the wannier functions (output of Wannier90)
!! eigen(mband*nkpt*nsppol)=Kohn-Sham Electronic eigenvalues.
!! U_matrix(mwan,mwan,nkpt,nsppol)=Matrix defining the Wannier gauge.
!! U_matrix_opt(:,:,:,:)=Matrix defining the optimal subspace. Used only if we have used the disentanglement procedure.
!! have_disentangled(nsppol)=Defines whether a disentanglement has been done (for each spin, separately).
!! band_in(mband,nsppol)=Number of bands used to construct the Wannier functions.
!! lwindow(mband,nkpt,nsppol)=TRUE if the state falls within the energy window used for the disentanglement.
!!
!! OUTPUT
!!  WanData<WannierData>=The object containing all data needed to perform the interpolation.
!!
!! PARENTS
!!      mlwfovlp
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine WrapInitWanData(WanData,mband,nkpt,nsppol,nntot,mwan,nwan_W,have_disentangled_W,eigen,&
& band_in_W,lwindow_W,wann_centres_W,wann_spreads_W,spreadw_W,U_matrix_W,U_matrix_opt_W,Hdr,W90version,title)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WrapInitWanData'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,mband,nntot,nsppol,mwan
 real(dp),intent(in) :: W90version
 type(Hdr_type),intent(in) :: Hdr
 type(WannierData),intent(out) :: WanData
 character(len=*),intent(in) :: title
!arrays
 integer,intent(in) :: nwan_W
 real(dp),intent(in) :: spreadw_W(3)
 real(dp),intent(in) :: wann_centres_W(3,mwan) 
 real(dp),intent(in) :: wann_spreads_W(mwan)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 complex(dpc),intent(in) :: U_matrix_W(mwan,mwan,nkpt)
 complex(dpc),intent(in) :: U_matrix_opt_W(mband,mwan,nkpt)
 logical,intent(in) :: have_disentangled_W
 logical,intent(in) :: band_in_W(mband)
 logical,intent(in) :: lwindow_W(mband,nkpt)

!Local variables-------------------------------
! Used to call the correct creation method
 integer :: nwan(nsppol)
 real(dp) :: spreadw(3,nsppol),wann_centres(3,mwan,nsppol),wann_spreads(mwan,nsppol)
 complex(dpc) :: U_matrix(mwan,mwan,nkpt,nsppol)
 complex(dpc) :: U_matrix_opt(mband,mwan,nkpt,nsppol)
 logical :: have_disentangled(nsppol)
 logical :: band_in(mband,nsppol)
 logical :: lwindow(mband,nkpt,nsppol)

! *************************************************************************

 nwan(:)              =nwan_W
 spreadw(:,:)         =SPREAD(spreadw_W(:),DIM=2,NCOPIES=nsppol) 
 wann_centres(:,:,:)  =SPREAD(wann_centres_W(:,:),DIM=3,NCOPIES=nsppol)
 wann_spreads(:,:)    =SPREAD(wann_spreads_W(:),DIM=2,NCOPIES=nsppol)
 U_matrix(:,:,:,:)    =SPREAD(U_matrix_W(:,:,:),DIM=4,NCOPIES=nsppol)
 U_matrix_opt(:,:,:,:)=SPREAD(U_matrix_opt_W(:,:,:),DIM=4,NCOPIES=nsppol)
 have_disentangled(:) =have_disentangled_W
 band_in(:,:)         =SPREAD(band_in_W(:)  ,DIM=2,NCOPIES=nsppol)
 lwindow(:,:,:)       =SPREAD(lwindow_W(:,:),DIM=3,NCOPIES=nsppol)

 call InitWanData(WanData,mband,nkpt,nsppol,nntot,mwan,nwan,have_disentangled,eigen,&
& band_in,lwindow,wann_centres,wann_spreads,spreadw,U_matrix,U_matrix_opt,Hdr,W90version,title)

end subroutine WrapInitWanData
!!***

!!****f* m_wannier2abinit/NullifyWanData_
!! NAME
!! NullifyWanData_  
!!
!! FUNCTION
!!  Nullify all pointers defined in the WannierData data type [Private]
!!
!! INPUTS
!!  WanData<WannierData>=The structure. See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  All pointers set to NULL.
!!
!! PARENTS
!!      m_wannier2abinit
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine NullifyWanData_(WanData)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'NullifyWanData_'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(WannierData),intent(inout) :: WanData
!Local variables-------------------------------

! *************************************************************************

 nullify(WanData%nwan        )

 nullify(WanData%eigen       )
 nullify(WanData%spreadw     )
 nullify(WanData%wann_centres)
 nullify(WanData%wann_spreads)

 nullify(WanData%U_matrix    )
 nullify(WanData%U_matrix_opt)

 nullify(WanData%band_in     )
 nullify(WanData%lwindow     )
 nullify(WanData%have_disentangled)

 !TODO
 !call nullify(WanData%Hdr)

 ! === Variables used to interpolate ===
 nullify(WanData%irvec  )
 nullify(WanData%ndegen )
 nullify(WanData%ndimwin)
 nullify(WanData%hamWR  )

end subroutine NullifyWanData_
!!***

!!****f* m_wannier2abinit/DestroyWanData
!! NAME
!! DestroyWanData
!!
!! FUNCTION
!! Free the dynamic memory in the object.
!!
!! INPUTS
!!  WanData<WannierData>=The structure. See also SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  All associated pointers are deallocated.
!!
!! PARENTS
!!      m_gwannier,mlwfovlp
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine DestroyWanData(WanData)
 
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'DestroyWanData'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(WannierData),intent(inout) :: WanData
!Local variables-------------------------------

! *************************************************************************

 if (associated(WanData%nwan             ))   then
   ABI_DEALLOCATE(WanData%nwan)
 end if
 if (associated(WanData%eigen            ))   then
   ABI_DEALLOCATE(WanData%eigen)
 end if
 if (associated(WanData%spreadw          ))   then
   ABI_DEALLOCATE(WanData%spreadw)
 end if
 if (associated(WanData%wann_centres     ))   then
   ABI_DEALLOCATE(WanData%wann_centres)
 end if
 if (associated(WanData%wann_spreads     ))   then
   ABI_DEALLOCATE(WanData%wann_spreads)
 end if

 if (associated(WanData%U_matrix         ))   then
   ABI_DEALLOCATE(WanData%U_matrix)
 end if
 if (associated(WanData%U_matrix_opt     ))   then
   ABI_DEALLOCATE(WanData%U_matrix_opt)
 end if

 if (associated(WanData%band_in          ))   then
   ABI_DEALLOCATE(WanData%band_in)
 end if
 if (associated(WanData%lwindow          ))   then
   ABI_DEALLOCATE(WanData%lwindow)
 end if
 if (associated(WanData%have_disentangled))   then
   ABI_DEALLOCATE(WanData%have_disentangled)
 end if

 call hdr_clean(WanData%Hdr)

 ! === Variables used to interpolate
 if (associated(WanData%irvec  ))  then
   ABI_DEALLOCATE(WanData%irvec)
 end if
 if (associated(WanData%ndegen ))  then
   ABI_DEALLOCATE(WanData%ndegen)
 end if
 if (associated(WanData%ndimwin))  then
   ABI_DEALLOCATE(WanData%ndimwin)
 end if
 if (associated(WanData%hamWR  ))  then
   ABI_DEALLOCATE(WanData%hamWR)
 end if

end subroutine DestroyWanData
!!***

!!****f* m_wannier2abinit/DumpWanData
!! NAME
!! DumpWanData
!!
!! FUNCTION
!! Dump the object on an external binary file.
!!
!! INPUTS
!!  WanData<WannierData>=The structure. See also SIDE EFFECTS
!!  fname=Name of the external file.
!!  accesswff=Flag defining the file format. 
!!   * 0 for plain Fortran file
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      mlwfovlp
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine DumpWanData(WanData,fname,accesswff)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'DumpWanData'
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff 
 !character(len=*)  This might create problems with other abinit routines!
 character(len=fnlen),intent(in) :: fname
 type(WannierData),intent(inout) :: WanData

!Local variables-------------------------------
!scalars
 integer :: unt,fform,rdwr,isp
 character(len=500) :: msg

! *************************************************************************

 fform=502 ; rdwr=2
 ! Add date and time?

 SELECT CASE (accesswff)

 CASE (IO_MODE_FORTRAN)
  unt=get_unit()
  open(unt,file=fname,form='unformatted')
  
  ! === Write the (in)famous Abinit Header
  call hdr_io(fform,WanData%Hdr,rdwr,unt)

  ! === Info ===
  write(unt)WanData%title, WanData%W90version, WanData%WDversion

  ! === Dimensions and parameters of the calculation  ===
  write(unt)WanData%mband, WanData%nkpt, WanData%nntot, WanData%nsppol, WanData%mwan
  do isp=1,WanData%nsppol
   write(unt)WanData%nwan(isp),WanData%have_disentangled(isp)
   write(unt)WanData%band_in(:,isp) 
   write(unt)WanData%lwindow(:,:,isp) !this is only for disentangled!
  end do

  do isp=1,WanData%nsppol
   ! * Energies
   write(unt)WanData%eigen(:,:,isp)
   ! * Results for Wan spread.
   write(unt)WanData%spreadw(:,isp), WanData%wann_centres(:,:,isp), WanData%wann_spreads(:,isp)

   write(unt)WanData%U_matrix(:,:,:,isp)
   if (WanData%have_disentangled(isp)) then 
    write(unt)WanData%U_matrix_opt(:,:,:,isp)
   end if
  end do

  close(unt)

 CASE DEFAULT 
  write(msg,'(a,i4,a)')' accesswff= ',accesswff,'. Wrong value or not implemented error.'
  MSG_BUG(msg)
 END SELECT

end subroutine DumpWanData
!!***

!!****f* m_wannier2abinit/ReadWanData
!! NAME
!! ReadWanData
!!
!! FUNCTION
!! Read the object from file.
!!
!! INPUTS
!!  fname=Name of the external file.
!!  accesswff=Flag defining the file format. 
!!   * 0 for plain Fortran file
!!
!! OUTPUT
!!  WanData<WannierData>=The structure initialized with all values needed for the Wannier interpolation. 
!!
!! PARENTS
!!      m_gwannier
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine ReadWanData(WanData,fname,accesswff)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ReadWanData'
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff 
 !character(len=*)  This might create problems with other abinit routines!
 character(len=fnlen),intent(in) :: fname
 type(WannierData),intent(out) :: WanData

!Local variables-------------------------------
!scalars
 integer :: iband,ikpt,isp,nband_k,unt,fform,rdwr,idx,mwan
 integer :: mband,nntot,nkpt,nsppol,WDversion
 real(dp) :: W90version
 character(len=500) :: msg
 character(len=800) :: title
 type(Hdr_type) ::  Hdr
!Arrays 
 integer,allocatable ::  nwan(:)
 real(dp),allocatable :: spreadw(:,:)
 real(dp),allocatable :: wann_centres(:,:,:),wann_spreads(:,:),eigen(:,:,:),eigpack(:)
 logical,allocatable ::  band_in(:,:),lwindow(:,:,:),have_disentangled(:)
 complex(dp),allocatable :: U_matrix(:,:,:,:),U_matrix_opt(:,:,:,:) 
! *************************************************************************

 rdwr=1 ! Add date and time?

 SELECT CASE (accesswff)
 CASE (IO_MODE_FORTRAN)
  unt=get_unit()
  open(unt,file=fname,form='unformatted')
  
  ! === Read the (in)famous Abinit Header
  call hdr_io(fform,Hdr,rdwr,unt)

  ! === Read Info ===
  read(unt)title, W90version, WDversion

  ! === Read Dimensions and parameters of the calculation  ===
  read(unt)mband, nkpt, nntot, nsppol, mwan

  ABI_ALLOCATE(nwan,(nsppol))
  ABI_ALLOCATE(have_disentangled,(nsppol))
  ABI_ALLOCATE(band_in,(mband,nsppol))
  ABI_ALLOCATE(lwindow,(mband,nkpt,nsppol))
  
  do isp=1,nsppol
   read(unt)nwan(isp), have_disentangled(isp)
   read(unt)band_in(:,isp) 
   read(unt)lwindow(:,:,isp) !this is only for disentangled!
  end do

  ABI_ALLOCATE(eigen,(mband,nkpt,nsppol))
  ABI_ALLOCATE(spreadw,(3,nsppol))
  ABI_ALLOCATE(wann_centres,(3,mwan,nsppol))
  ABI_ALLOCATE(wann_spreads,(mwan,nsppol))
  ABI_ALLOCATE(U_matrix,(mwan,mwan,nkpt,nsppol))
  if (ANY(have_disentangled)) then
   ABI_ALLOCATE(U_matrix_opt,(mband,mwan,nkpt,nsppol))
  end if

  do isp=1,nsppol
   ! * Read Energies
   read(unt)eigen(:,:,isp)
   ! * Read Results for Wan spread.
   read(unt)spreadw(:,isp), wann_centres(:,:,isp), wann_spreads(:,isp)

   read(unt)U_matrix(:,:,:,isp)
   if (have_disentangled(isp)) then 
    read(unt)U_matrix_opt(:,:,:,isp)
   end if
  end do

  close(unt)

  ! I have to temporarily pack the energies to follow the funny Abinit convention.
  ! TODO call pack_eneocc, once the build-sytem will treat correctly dependencies.
  ABI_ALLOCATE(eigpack,(mband*nkpt*nsppol))
  idx=0
  do isp=1,nsppol
   do ikpt=1,nkpt
    nband_k=Hdr%nband(ikpt+(isp-1)*nkpt)
    do iband=1,nband_k
     idx=idx+1
     eigpack(idx)=eigen(iband,ikpt,isp)
    end do
   end do
  end do

  call InitWanData(WanData,mband,nkpt,nsppol,nntot,mwan,nwan,have_disentangled,eigpack,&
&  band_in,lwindow,wann_centres,wann_spreads,spreadw,U_matrix,U_matrix_opt,Hdr,W90version,title)

  ! Brand-new WanData object, deallocate temporary stuff.
  ABI_DEALLOCATE(eigpack)

  call hdr_clean(Hdr)
  ABI_DEALLOCATE(nwan)
  ABI_DEALLOCATE(band_in)
  ABI_DEALLOCATE(lwindow)
  ABI_DEALLOCATE(eigen)
  ABI_DEALLOCATE(spreadw)
  ABI_DEALLOCATE(wann_centres)
  ABI_DEALLOCATE(wann_spreads)
  ABI_DEALLOCATE(U_matrix)
  if (ANY(have_disentangled))  then
    ABI_DEALLOCATE(U_matrix_opt)
  end if
  ABI_DEALLOCATE(have_disentangled)

 CASE DEFAULT 
  write(msg,'(a,i4,a)')' accesswff= ',accesswff,'. Wrong value or not implemented error.'
  MSG_BUG(msg)
 END SELECT

end subroutine ReadWanData
!!***

!!****f* m_wannier2abinit/PrintWanData
!! NAME
!! PrintWanData
!!
!! FUNCTION
!! Printout of the content of the object (useful for debugging)
!!
!! INPUTS
!! [unit}=Unit number for output, defaults to std_out.
!! [prtvol]=Verbosity level, defaults to 0
!! [mode_paral]=Either "COLL" or "PERS", defaults to "COLL"
!! WanData<WannierData>=The object whose basic info have to be printed. 
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_gwannier
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine PrintWanData(WanData,unit,mode_paral,prtvol)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'PrintWanData'
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 type(WannierData),intent(inout) :: WanData

!Local variables-------------------------------
!scalars
 integer :: rdwr,fform,unt,verbose,isp,ikpt
 character(len=4) :: mode

! *************************************************************************

 unt=std_out ; if (PRESENT(unit))       unt=unit
 verbose=0   ; if (PRESENT(prtvol))     verbose=prtvol
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 ! TODO define fform and version, solve well known problems with build system
 rdwr=4; fform=502 ! Define fform for WAN HDR?
 call hdr_io(fform,WanData%Hdr,rdwr,unt) !this call requires INOUT intent.

 ! TODO cleaning and beautification
 write(unt,*)' Wannier90 version ................... ',WanData%W90version 
 write(unt,*)' Wannier file version ................ ',WanData%WDversion 
 write(unt,*)' Maximun number of bands processed ... ',WanData%mband
 write(unt,*)' Number of spin polarizations ........ ',WanData%nsppol
 write(unt,*)' Number of Wannier functions ......... ',WanData%nwan(:)
 write(unt,*)' Number of k-points and neighbours ... ',WanData%nkpt,WanData%nntot
 write(unt,*)' Have disentangled ................... ',WanData%have_disentangled(:)
 write(unt,*)' spreadw ',WanData%spreadw     
 write(unt,*)' Wannier centres ..................... ',WanData%wann_centres
 write(unt,*)' Wannier spreads ..................... ',WanData%wann_spreads

 do isp=1,WanData%nsppol
  write(unt,*)' For spin ',isp,'; ',COUNT(WanData%band_in(:,isp)),' bands are included '
  if (WanData%have_disentangled(isp)) then
   do ikpt=1,WanData%nkpt
    write(unt,*)'At k-point ',ikpt,'; ',COUNT(WanData%lwindow(:,ikpt,isp)),' bands are included'
   end do
  end if
 end do

 if (verbose>0) then
  do isp=1,WanData%nsppol
   write(unt,*)' Energies ',WanData%eigen(:,:,isp)
   !write(unt,*)' U matrix ',WanData%U_matrix(:,:,:,isp)
   if (WanData%have_disentangled(isp)) then
    !write(unt,*)'U_matrix_opt',WanData%U_matrix_opt(:,:,:,isp)
   end if
  end do
 end if

end subroutine PrintWanData
!!***

!!****f* m_wannier2abinit/MakeWannierHR
!! NAME
!! MakeWannierHR
!!
!! FUNCTION
!! Construct the matrix elements of the KS Hamiltonian in real space starting
!! from the data stored in the WanData type. 
!!
!! INPUTS
!!  WanData<WannierData>=Results of the Wannierization. See also SIDE EFFECTS
!!  kptrlatt(3,3)=Array defining the BZ sampling, see also notes below.
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  The following quantities defined in the WannierData data type are allocated and/or calculated.
!!  %nrpts=Number of points in the Wigner-Seitz cell (WS)
!!  %irvec(3,nrpts)=Lattice vectors in the WS cell in the basis of the lattice vectors defining the unit cell
!!  %ndegen(nrpts)=Degeneracy of each point. It will be weighted using 1/ndegen(ii)
!!  %hamWR(mwan,mwan,nrpts,nsppol)=The Hamiltonian in the Wannier gauge.
!!
!! TODO 
!! kptrlatt is passed in input but it should be reported in the abinit header as well as shiftk and kptopt
!! Add additional options to cut the matrix.
!!
!! PARENTS
!!      m_gwannier
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine MakeWannierHR(WanData,kptrlatt)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MakeWannierHR'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(WannierData),target,intent(inout) :: WanData
!scalars
 integer,intent(in) :: kptrlatt(3,3)

!Local variables-------------------------------
!scalars
 integer :: isp,ikpt,idx,ib,nwan,nkpt,mband,ii,jj,kk,iR,mwan,nrpts
 real(dp) :: rdotk
 complex(dpc) :: phase
 logical :: ass_irvec,ass_ndegen

!Arrays 
 real(dp) :: irvec_tmp(3)
 integer,allocatable :: shift_vec(:,:)
 integer,pointer :: irvec(:,:),ndegen(:)
 real(dp) :: rmet(3,3)
 real(dp),pointer :: kpt(:,:),rprimd(:,:)
 real(dp),allocatable :: eigval_opt(:,:)
 real(dp),allocatable :: eigval2(:,:)
 complex(dpc),pointer :: U_matrix(:,:,:),U_matrix_opt(:,:,:)
 complex(dpc),allocatable :: ham_k(:,:,:),ham_r(:,:,:)

! *************************************************************************

 ! === If not yet done calculate Wigner-Seitz cell ===
 ! * Here we should consider the possibility of a different real mesh.
 ass_irvec =associated(WanData%irvec )
 ass_ndegen=associated(WanData%ndegen)

 if (.not.ass_irvec.or..not.ass_ndegen) then
  ABI_CHECK(.not.ass_irvec,'irvec should not be associated')
  ABI_CHECK(.not.ass_ndegen,'ndegen should not be associated')
  rprimd => WanData%Hdr%rprimd
  do ii=1,3 ! Compute real space metrics
   rmet(ii,:)=rprimd(1,ii)*rprimd(1,:)+&
&             rprimd(2,ii)*rprimd(2,:)+&
&             rprimd(3,ii)*rprimd(3,:)
  end do

  call wigner_seitz(gamma_point,LMAX,kptrlatt,rmet,nrpts,irvec,ndegen,prtvol=0)
  WanData%nrpts=nrpts
  ABI_ALLOCATE(WanData%irvec,(3,nrpts))
  WanData%irvec =irvec
  ABI_ALLOCATE(WanData%ndegen,(nrpts))
  WanData%ndegen=ndegen
  ABI_DEALLOCATE(irvec)
  ABI_DEALLOCATE(ndegen)
 end if

 ! === Allocate Hamiltonian in the Wannier gauge on the ab initio grid ===
 ! * The Max number of Wan bands over spin is used to allocate sufficient space.
 mband =  WanData%mband
 nkpt  =  WanData%nkpt
 mwan  =  WanData%mwan 
 kpt   => WanData%Hdr%kptns(1:3,1:nkpt)

 if (associated(WanData%hamWR))  then
   ABI_DEALLOCATE(WanData%hamWR)
 end if
 ABI_ALLOCATE(WanData%hamWR,(mwan,mwan,WanData%nrpts,WanData%nsppol))
 WanData%hamWR=czero

 do isp=1,WanData%nsppol
  nwan = WanData%nwan(isp)
  ABI_ALLOCATE(ham_k,(nwan,nwan,nkpt))
  ham_k = czero

  ABI_ALLOCATE(eigval2,(nwan,nkpt))
  eigval2 = zero

  if (WanData%have_disentangled(isp)) then
   !
   ! === Slim down eigval to contain states within the outer window ===
   ABI_ALLOCATE(eigval_opt,(mband,nkpt))
   eigval_opt=zero
   do ikpt=1,WanData%nkpt
    idx=0
    do ib=1,WanData%mband
     if (WanData%lwindow(ib,ikpt,isp)) then
      idx=idx+1
      eigval_opt(idx,ikpt)=WanData%eigen(ib,ikpt,isp)
     end if
    end do
   end do
   
   ! === Rotate eigval into the optimal subspace ===
   ! * In general eigval would be a matrix at each kpoint
   !   but we choose u_matrix_opt such that the Hamiltonian is
   !   diagonal at each kpoint. (I guess we should check it here)
   u_matrix_opt => WanData%U_matrix_opt(:,:,:,isp)
   
   do ikpt=1,WanData%nkpt
    do jj=1,nwan
     do kk=1,WanData%ndimwin(ikpt,isp)
      eigval2(jj,ikpt)=eigval2(jj,ikpt)+ & 
&     eigval_opt(kk,ikpt)*real(conjg(u_matrix_opt(kk,jj,ikpt))*u_matrix_opt(kk,jj,ikpt),dp)
     end do
    end do
   end do
   ABI_DEALLOCATE(eigval_opt)

  else
   !FIXME here there is a bug if we dont start from 1
   eigval2(1:nwan,:)=WanData%eigen(1:nwan,:,isp)
  end if

  ! At this point eigval2 contains nwan values belonging to the Wannier subspace.
  ! Rotate the Hamiltonian into the basis of smooth Bloch states.
  !          H(k)=U^{dagger}(k).H_0(k).U(k)
  ! Note: we enforce hermiticity here
  u_matrix => WanData%U_matrix(:,:,:,isp)

  do ikpt=1,WanData%nkpt
   do jj=1,nwan
    do ii=1,jj
     do kk=1,nwan
      ham_k(ii,jj,ikpt) = ham_k(ii,jj,ikpt)+ &
&      eigval2(kk,ikpt)*CONJG(u_matrix(kk,ii,ikpt))*u_matrix(kk,jj,ikpt)
     end do
     if (ii<jj) ham_k(jj,ii,ikpt)=CONJG(ham_k(ii,jj,ikpt))
    end do
   end do
  end do

  ABI_DEALLOCATE(eigval2)

  ! === Fourier transform rotated hamiltonian into WF basis ===
  ! * H_ij(k) --> H_ij(R) = (1/N_kpts) sum_k e^{-ikR} H_ij(k)
  ABI_ALLOCATE(ham_r,(nwan,nwan,WanData%nrpts))
  ham_r=czero

  !TODO
  !if (.not.use_translation) then
  if (.TRUE.) then
   do iR=1,WanData%nrpts
    do ikpt=1,WanData%nkpt
     rdotk=two_pi*DOT_PRODUCT(kpt(:,ikpt),REAL(WanData%irvec(:,iR),dp))
     phase=EXP(-j_dpc*rdotk) !/real(num_kpts,dp)
     ham_r(:,:,iR)=ham_r(:,:,iR)+phase*ham_k(:,:,ikpt)
    end do
   end do
   ham_r=ham_r/WanData%nkpt
   !have_translated = .false.

  else
   ABI_ALLOCATE(shift_vec,(3,nwan))
   stop "not implemented" !TODO
   !call internal_translate_centres()
   do iR=1,nrpts
    do ikpt=1,WanData%nkpt
     do ii=1,nwan
      do jj=1,nwan
       ! ham_r(j,i,iR) interaction btw j at 0 and i at irvec(:,iR)
       irvec_tmp(:)=WanData%irvec(:,iR)+shift_vec(:,ii)-shift_vec(:,jj)   
       rdotk=two_pi*DOT_PRODUCT(kpt(:,ikpt),real(irvec_tmp(:),dp))
       phase=EXP(-j_dpc*rdotk) !/real(num_kpts,dp)
       ham_r(jj,ii,iR)=ham_r(jj,ii,iR)+phase*ham_k(jj,ii,ikpt)
      end do
     end do
    end do
   end do
   ham_r=ham_r/WanData%nkpt
   !have_translated = .true.
   ABI_DEALLOCATE(shift_vec)
  end if

  ! === Save Hamiltonian in real space ===
  WanData%hamWR(1:nwan,1:nwan,:,isp)=ham_r(:,:,:)
  ABI_DEALLOCATE(ham_r)
 end do !isp

end subroutine MakeWannierHR
!!***

!!****f* m_wannier2abinit/PlotWannierHR
!! NAME
!! PlotWannierHR
!!
!! FUNCTION
!!  Plot the matrix elements in real space of the Hamitonian in the Wannier Gauge.
!!
!! INPUTS
!!  fname=File name for output.
!!  [cut_mode]=If present, it defines which type of cutoff in real space has to be used.
!!  WanData<WannierData>=The object containing the quantities used for the Wannier interpolation.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine PlotWannierHR(WanData,fname,cut_mode)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'PlotWannierHR'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: fname
 character(len=*),optional,intent(in) :: cut_mode
 type(WannierData),target,intent(in) :: WanData
!scalars

!Local variables-------------------------------
!scalars
 integer :: unt,iR,irad,isp
 integer :: nwan,ii
 real(dp) :: datum(WanData%nsppol)
 character(len=500) :: msg
!Arrays 
 integer,allocatable :: sph2arr(:)
 real(dp) :: rmet(3,3)
 real(dp),pointer :: rprimd(:,:)
 real(dp),allocatable :: length(:)

! *************************************************************************

 if (PRESENT(cut_mode)) then
  if (WanData%cut_mode/=cut_mode) then 
   msg="WanData%cut_mode/=cut_mode not implemented"
   MSG_ERROR(msg)
  end if
 end if

 rprimd => WanData%Hdr%rprimd
 do ii=1,3 ! Compute real space metrics
  rmet(ii,:)=rprimd(1,ii)*rprimd(1,:)+&
&            rprimd(2,ii)*rprimd(2,:)+&
&            rprimd(3,ii)*rprimd(3,:)
 end do

 ! Sort points by length
 ABI_ALLOCATE(length,(WanData%nrpts))
 ABI_ALLOCATE(sph2arr,(WanData%nrpts))
 do iR=1,WanData%nrpts
  length (iR)=normv(WanData%irvec(:,iR),rmet,'R')
  sph2arr(iR)=iR
 end do
 call sort_dp(WanData%nrpts,length,sph2arr,tol6)
 ABI_DEALLOCATE(length)

 unt=get_unit()
 open(file=fname,unit=unt,form='formatted')

 write(unt,'(a)')'# '
 write(unt,'(a)')'# Matrix elements of H in the Wannier representation '
 write(unt,'(a)')'# Length(ii), Lattice point Rws(1:3,ii), (MAXVAL |hamW_ij(R,isp)|, isp=1,nsppol)'
 write(unt,'(a)')'# '
 do irad=1,WanData%nrpts
  ii=sph2arr(irad)
  do isp=1,WanData%nsppol
   nwan=WanData%nwan(isp)
   datum(isp)=MAXVAL(ABS(WanData%hamWR(1:nwan,1:nwan,ii,isp)))
  end do
  write(unt,'(es14.6,3i3,(es14.6))')length(ii),WanData%irvec(:,ii),datum(:)
 end do

 ABI_DEALLOCATE(sph2arr)
 close(unt)

end subroutine PlotWannierHR
!!***


!!****f* m_wannier2abinit/WanMatInterpol_dpc
!! NAME
!! WanMatInterpol_dpc
!!
!! FUNCTION
!! Interpolate matrix elements on a denser k-mesh using Wannier representation
!! Version for double precision complex values.
!!
!! INPUTS
!! nk4intp=Number of interpolating k-points.
!! nkab=Number of ab-initio calculated k-points.
!! nsppol=Number of spin polarizations.
!! mb4intp=Max number of bands in the ab-inition calculated matrix elements.
!! ishermitian=TRUE if the operator to be interpolated is Hermitian
!! WanData<WannierData>=The object storing the quantities used for the Wannier interpolation.
!! k4intp(3,nk4intp)=Reduced coordinated of the interpolating k-points.
!! matrix_in(mb4intp,mb4intp,nkab,nsppol)=Ab-initio calculated matrix elements.
!!
!! OUTPUT
!! matrix_out(mb4intp,mb4intp,nk4intp,nsppol)=Interpolated matrix elements.
!!
!! PARENTS
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine WanMatInterpol_dpc(WanData,nkab,nk4intp,k4intp,ishermitian,nsppol,mb4intp,matrix_in,matrix_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WanMatInterpol_dpc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nk4intp,nkab,nsppol,mb4intp
 logical,intent(in) :: ishermitian
 type(WannierData),intent(in) :: WanData
!Arrays
 real(dp),intent(in) :: k4intp(3,nk4intp)
 complex(dpc),intent(in)  :: matrix_in (mb4intp,mb4intp,nkab,   nsppol)
 complex(dpc),intent(out) :: matrix_out(mb4intp,mb4intp,nk4intp,nsppol)

!Local variables-------------------------------
!scalars
! *********************************************************************

 ! === Consistency check === 
 ABI_CHECK(WanData%nkpt==nkab,'WanData%nkpt/=nkab')
 ABI_CHECK(WanData%nsppol==nsppol,'WanData%nsppol/=nsppol')
 ABI_CHECK(mb4intp==WanData%mwan,'WanData%mwan/=mb4intp')

 ! === Call the private method to do the dirty job ===
 call DoWanMatInterpol_(WanData,nk4intp,k4intp,ishermitian,matrix_in,matrix_out)

end subroutine WanMatInterpol_dpc
!!***

!!****f* ABINIT/WanMatInterpol_dp
!! NAME
!! WanMatInterpol_dp
!!
!! FUNCTION
!! Interpolate matrix elements on a denser k-mesh using Wannier representation
!! Version for double precision real values.
!!
!! INPUTS
!! nk4intp=Number of interpolating k-points.
!! nkab=Number of ab-initio calculated k-points.
!! nsppol=Number of spin polarizations.
!! mb4intp=Max number of bands in the ab-inition calculated matrix elements.
!! ishermitian=TRUE if the operator to be interpolated is Hermitian
!! WanData<WannierData>=The object storing the quantities used for the Wannier interpolation.
!! k4intp(3,nk4intp)=Reduced coordinated of the interpolating k-points.
!! matrix_in(mb4intp,mb4intp,nkab,nsppol)=Ab-initio calculated matrix elements.
!!
!! OUTPUT
!! matrix_out(mb4intp,mb4intp,nk4intp,nsppol)=Interpolated matrix elements.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE
subroutine WanMatInterpol_dp(WanData,nkab,nk4intp,k4intp,ishermitian,nsppol,mb4intp,matrix_in,matrix_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WanMatInterpol_dp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nk4intp,nkab,nsppol,mb4intp
 logical,intent(in) :: ishermitian
 type(WannierData),intent(in) :: WanData
!Arrays
 real(dp),intent(in) :: k4intp(3,nk4intp)
 real(dp),intent(in)  :: matrix_in (mb4intp,mb4intp,nkab,   nsppol)
 real(dp),intent(out) :: matrix_out(mb4intp,mb4intp,nk4intp,nsppol)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
!Arrays
 complex(dpc) :: c_matrix_in (mb4intp,mb4intp,nkab,   nsppol)
 complex(dpc) :: c_matrix_out(mb4intp,mb4intp,nk4intp,nsppol)
! *********************************************************************

 ! === Consistency check === 
 ABI_CHECK(WanData%nkpt==nkab,'%nkpt/=nkab')
 ABI_CHECK(WanData%nsppol==nsppol,'%nsppol/=nsppol')
 ABI_CHECK(mb4intp==WanData%mwan,'%mwan/=mb4intp')
 if (.not.ishermitian) then 
  msg='ishermitian should be true for the real version'
  MSG_ERROR(msg)
 end if

 ! === Transfer input data to complex array ===
 c_matrix_in=CMPLX(matrix_in)

 ! === Call the private method to do the dirty job ===
 call DoWanMatInterpol_(WanData,nk4intp,k4intp,ishermitian,c_matrix_in,c_matrix_out)

 matrix_out=REAL(c_matrix_out)

end subroutine WanMatInterpol_dp
!!***

!!****f* m_wannier2abinit/DoWanMatInterpol_
!! NAME
!! DoWanMatInterpol_ [PRIVATE]
!!
!! FUNCTION
!! Interpolate matrix elements on a denser k-mesh using Wannier representation.
!!
!! INPUTS
!! nk4intp=Number of k-points of the interpolating set.
!! ishermitian=TRUE if the operator is Hermitian.
!! WanData<WannierData>=Object storing values needed for the interpolation.
!! k4intp(3,nk4intp)=Reduced k-points of the interpolating set.
!! matrix_in(WanData%mwan,WanData%mwan,WanData%nkpt,WanData%nsppol)=Ab-initio calculated matrix elements.
!!
!! OUTPUT
!! matrix_out(WanData%mwan,WanData%mwan,nk4intp,WanData%nsppol)=Interpolated matrix elements.
!!
!! PARENTS
!!      m_wannier2abinit
!!
!! NOTES
!!  Do not use this method directly (it is private). Use the public wrapper 
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine DoWanMatInterpol_(WanData,nk4intp,k4intp,ishermitian,matrix_in,matrix_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'DoWanMatInterpol_'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nk4intp
 logical,intent(in) :: ishermitian
 type(WannierData),target,intent(in) :: WanData
!Arrays
 real(dp),intent(in) :: k4intp(3,nk4intp)
 complex(dpc),target,intent(in)  :: matrix_in (WanData%mwan,WanData%mwan,WanData%nkpt,WanData%nsppol)
 complex(dpc),target,intent(out) :: matrix_out(WanData%mwan,WanData%mwan,nk4intp,     WanData%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ikab,ikint,nwan,nrpts,isp
 integer :: ierr,iR,ii,jj,kk,ll,info,nfound
 real(dp) :: rdotk,ABSTOL
 complex(dpc) :: fact,phase
 character(len=500) :: msg
!arrays
 integer,pointer :: irvec(:,:),ndegen(:)
 integer,allocatable :: iwork(:),ifail(:)
 real(dp),allocatable :: eig_interp(:,:),rwork(:)     
 complex(dpc),pointer :: Hwan_R(:,:,:),Ub2w(:,:,:)
 complex(dpc),allocatable  :: Ham_pack(:),Ham_kint(:,:),Owan_kint(:,:)
!complex(dpc),allocatable  :: Ham_r_cut(:,:,:)
 complex(dpc),allocatable  :: Ub2w_intp(:,:),cwork(:),Owan_kab(:,:,:),Owan_R(:,:,:) 
 complex(dpc),pointer :: Oham_kab(:,:,:)
 !integer, allocatable :: irvec_cut(:,:)
 !integer              :: irvec_max(3)
 !integer              :: nrpts_cut
! *********************************************************************

 ! =====================================
 ! === Perform Wannier interpolation ===
 ! ======================================

 if (.not.associated(WanData%irvec).or..not.associated(WanData%ndegen)) then
  msg='irvec and ndegen are not associated'
  MSG_ERROR(msg)
 end if

 nrpts  =  WanData%nrpts
 irvec  => WanData%irvec
 ndegen => WanData%ndegen
 
 matrix_out=czero

 do isp=1,WanData%nsppol

  ! === Rotate the Operator into the basis of smooth Bloch states ===
  ! * Owan(k)=U^*(k).0ham(k).U(k)

  nwan     =  WanData%nwan(isp) 
  Ub2w     => WanData%U_matrix(:,:,:,isp)
  Oham_kab => matrix_in(1:nwan,1:nwan,:,isp)

  ABI_ALLOCATE(Owan_kab,(nwan,nwan,WanData%nkpt))
  Owan_kab=czero

  if (ishermitian) then 
   ! * Impose Hermiticity of Owan_kab
   do ikab=1,WanData%nkpt

    do jj=1,nwan
     do ii=1,jj
      do kk=1,nwan
       do ll=1,nwan
        Owan_kab(ii,jj,ikab) = Owan_kab(ii,jj,ikab)+ &
&        Oham_kab(kk,ll,ikab)*CONJG(Ub2w(kk,ii,ikab))*Ub2w(ll,jj,ikab)
       end do
      end do
      if (ii< jj) Owan_kab(jj,ii,ikab)=CONJG(Owan_kab(ii,jj,ikab))
      if (ii==jj) Owan_kab(jj,ii,ikab)=half*(Owan_kab(ii,jj,ikab)+CONJG(Owan_kab(ii,jj,ikab)))
     end do
    end do

   end do

  else
   do ikab=1,WanData%nkpt
    Owan_kab(:,:,ikab)= &
&    MATMUL( CONJG(TRANSPOSE(Ub2w(:,:,ikab))),MATMUL(Oham_kab(:,:,ikab),Ub2w(:,:,ikab)))
    !eigval2(kk,ikab)*CONJG(Ub2w(kk,ii,ikab))*Ub2w(kk,jj,ikab)
   end do
  end if

  ! === Fourier transform rotated hamiltonian into Wannier basis ===
  ! * H_ij(k) --> H_ij(R) = (1/N_kpts) sum_k e^{-ikR} H_ij(k)
  ABI_ALLOCATE(Owan_R,(nwan,nwan,WanData%nrpts))
  Owan_R=czero

  do iR=1,WanData%nrpts
   do ikab=1,WanData%nkpt
    rdotk=two_pi*DOT_PRODUCT(WanData%Hdr%kptns(:,ikab),real(WanData%irvec(:,iR),dp))
    phase=EXP(-j_dpc*rdotk) 
    Owan_R(:,:,iR) = Owan_R(:,:,iR) + phase*Owan_kab(:,:,ikab)
   end do
  end do
  Owan_R=Owan_R/WanData%nkpt

  ABI_DEALLOCATE(Owan_kab)

  ! Pay attention here. Point only the Wannier functions for this spin i.e up to nwan.
  Hwan_R => WanData%hamWR(1:nwan,1:nwan,1:WanData%nrpts,isp)
  ! TODO add cut case at this level !if (cut) then ...

  ABI_ALLOCATE(eig_interp,(nwan,nk4intp))
  eig_interp=zero 

  ABI_ALLOCATE(Ham_pack,((nwan*(nwan+1))/2))
  ierr = ABI_ALLOC_STAT
  ABI_ALLOCATE(Ham_kint,(nwan,nwan))
  ierr = ABI_ALLOC_STAT
  ABI_ALLOCATE(Owan_kint,(nwan,nwan))
  ierr = ABI_ALLOC_STAT
  ABI_ALLOCATE(Ub2w_intp,(nwan,nwan))
  ierr = ABI_ALLOC_STAT
  ABI_ALLOCATE(cwork,(2*nwan))
  ABI_ALLOCATE(rwork,(7*nwan))
  ABI_ALLOCATE(iwork,(5*nwan))
  ABI_ALLOCATE(ifail,(nwan))
  !
  ! Cut H matrix in real-space
  !if (index(bands_plot_mode,'cut').ne.0)  call plot_cut_hr()

  ! === Interpolate the Hamiltonian at each kpoint 
  ! * Get Owan_kint on the interpolating grid then rotate to get Oham_kint. 

  do ikint=1,nk4intp

   Ham_kint(:,:)=czero
   Owan_kint(:,:)=czero
   !if (index(bands_plot_mode,'s-k').ne.0) then
    do iR=1,nrpts
     rdotk=two_pi*DOT_PRODUCT(k4intp(:,ikint),irvec(:,iR))
     fact=EXP(j_dpc*rdotk)/REAL(ndegen(iR),dp)
     Ham_kint  =Ham_kint  +fact*Hwan_R(:,:,iR)
     Owan_kint =Owan_kint +fact*Owan_R(:,:,iR)
    end do
   !TODO
   !else if (index(bands_plot_mode,'cut').ne.0) then
   ! do iR=1,nrpts_cut
   !  rdotk=two_pi*DOT_PRODUCT(k4intp(:,ikint),irvec_cut(:,iR))
   !  !!$[aam] check divide by ndegen?
   !  fact=EXP(j_dpc*rdotk)
   !  Ham_kint=Ham_kint+fact*Ham_r_cut(:,:,iR)
   !  Owan_kint =Owan_kint +fact*Owan_R(:,:,iR)
   ! end do
   !end if

   ! === Diagonalise H_k (->basis of eigenstates) ===
   do jj=1,nwan
    do ii=1,jj
     Ham_pack(ii+((jj-1)*jj)/2)=Ham_kint(ii,jj)
    end do
   end do

   ! TODO Do not need to allocate all the k-point in eig_interp
   ABSTOL=-one !; ABSTOL= 2*DLAMCH('S')
   call ZHPEVX('V','A','U',nwan,ham_pack,zero,zero,0,0,ABSTOL,&
    nfound,eig_interp(1,ikint),Ub2w_intp,nwan,cwork,rwork,iwork,ifail,info)

   if (info<0) then
    write(msg,'(a,i3,a)') 'THE ',-info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
    MSG_ERROR(msg)
   end if
   if (info>0) then
    write(msg,'(i3,a)')info,' EIGENVECTORS FAILED TO CONVERGE'
    MSG_ERROR(msg)
   end if

   ! * Rotate and save interpolated matrix elements.
   if (ishermitian) then 
    do jj=1,nwan
     do ii=1,jj
      do kk=1,nwan
       do ll=1,nwan
        matrix_out(ii,jj,ikint,isp) = matrix_out(ii,jj,ikint,isp) + &
&        Owan_kint(kk,ll)*CONJG(Ub2w_intp(kk,ii))*Ub2w_intp(ll,jj)
       end do
      end do
      if (ii< jj) matrix_out(jj,ii,ikint,isp) = CONJG(matrix_out(ii,jj,ikint,isp))
      if (ii==jj) matrix_out(jj,ii,ikint,isp) = half*(matrix_out(ii,jj,ikint,isp)+CONJG(matrix_out(ii,jj,ikint,isp)))
     end do
    end do

   else 
    matrix_out(1:nwan,1:nwan,ikint,isp) = MATMUL( CONJG(TRANSPOSE(Ub2w_intp)),MATMUL(Owan_kint,Ub2w_intp))
   end if

   !write(77,'(7f8.4)')k4intp(:,ikint),eig_interp(1:nwan,ikint)

  end do !ikint

  ! Interpolation for this spin done.
  ABI_DEALLOCATE(eig_interp)
  ABI_DEALLOCATE(Ham_pack)
  ABI_DEALLOCATE(Ham_kint)
  ABI_DEALLOCATE(Owan_kint)
  ABI_DEALLOCATE(Ub2w_intp)
  ABI_DEALLOCATE(cwork)
  ABI_DEALLOCATE(rwork)
  ABI_DEALLOCATE(iwork)
  ABI_DEALLOCATE(ifail)
  !if (allocated(Ham_r_cut)) deallocate(Ham_r_cut,stat=ierr)
  !if (allocated(irvec_cut)) deallocate(irvec_cut,stat=ierr)

 end do !isp

end subroutine DoWanMatInterpol_
!!***

!!****f* m_electrons/WannierInterpol
!! NAME
!! WannierInterpol
!!
!! FUNCTION
!!  Perform the Wannier interpolation of the eigenvalues stored in a 
!!  bandstructure_type. See also m_wannier2abinit.
!!  Taken from the Wannier90 code. Modified by MG to have an OO interface. 
!!
!! Copyright (C) 2007 Jonathan Yates, Arash Mostofi,          
!!  Young-Su Lee, Nicola Marzari, Ivo Souza, David Vanderbilt 
!! This file is distributed under the terms of the GNU        
!! General Public License. See the file `LICENSE' in          
!! the root directory of the present distribution, or         
!! http://www.gnu.org/copyleft/gpl.txt .                      
!!                                                            
!! INPUTS
!!  WData<WanData>=Structure containing quantities used for the Wannier interpolation.
!!  Kmesh<bz_mesh_type>=Structure defining the k-point sampling for the interpolating mesh.
!! 
!! OUTPUT
!!  BSt<Bandstructure_type>=The interpolated band structure.
!!
!! PARENTS
!!      m_gwannier
!!
!! CHILDREN
!!      metric,zhpevx
!!
!! SOURCE

subroutine WannierInterpol(WData,Kmesh,BSt)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'WannierInterpol'
 use interfaces_28_numeric_noabirule
 use interfaces_42_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(WannierData),intent(inout) :: WData
 type(Bandstructure_type),intent(out)  :: BSt
 type(BZ_mesh_type),intent(in) :: Kmesh
!arrays

!Local variables-------------------------------
!scalars
 integer :: ikpt,nwan,nrpts,isp,prtvol
 integer :: ierr,iR,ii,jj,info,nfound,spad
 real(dp) :: rdotk,ABSTOL,ucvol
 complex(dpc) :: fact
 character(len=500) :: msg

!arrays
 integer,pointer :: irvec(:,:),ndegen(:)
 integer,allocatable :: iwork(:),ifail(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: eig_interp(:,:),rwork(:)     
!real(dp),allocatable :: bands_proj(:,:)
 complex(dpc),pointer :: ham_r(:,:,:)
!complex(dpc),allocatable  :: ham_r_cut(:,:,:)
 complex(dpc),allocatable  :: ham_pack(:)
 complex(dpc),allocatable  :: ham_kprm(:,:)
 complex(dpc),allocatable  :: U_int(:,:)
 complex(dpc),allocatable  :: cwork(:)     
 !integer, allocatable :: irvec_cut(:,:)
 !integer              :: irvec_max(3)
 !integer              :: nrpts_cut
! *********************************************************************

 prtvol=0
 call metric(gmet,gprimd,-1,rmet,WData%Hdr%rprimd,ucvol)

 ! === Bst will contain only the interpolated bands ===
 !FIXME We hack a bit Bst but this has to be done in a much cleaner way to have a working object
 Bst%nkpt=Kmesh%nibz 
 ABI_ALLOCATE(Bst%kptns,(3,Bst%nkpt))
 ABI_ALLOCATE(Bst%wtk,(Bst%nkpt))
 Bst%kptns=Kmesh%ibz 
 Bst%wtk  =Kmesh%wt

 Bst%nsppol=WData%nsppol 
 Bst%mband=MAXVAL(WData%nwan)

 !FIXME this has to be initialized correctly
 Bst%nspinor=1

 ABI_ALLOCATE(Bst%nband,(Bst%nkpt*Bst%nsppol))
 do isp=1,Bst%nsppol
  spad=(isp-1)*Bst%nkpt
  Bst%nband(1+spad:Bst%nkpt+spad)=WData%nwan(isp)
 end do

 ABI_ALLOCATE(Bst%eig,(Bst%mband,Bst%nkpt,Bst%nsppol))

 ! =======================================
 ! === Start the Wannier interpolation ===
 ! =======================================
 if (.not.associated(WData%irvec).or..not.associated(WData%ndegen)) then
  msg='irvec and ndegen must be associated'
  MSG_ERROR(msg)
 end if
 nrpts  =  WData%nrpts
 irvec  => WData%irvec
 ndegen => WData%ndegen

 do isp=1,WData%nsppol

  ! Pay attention here: point the Wannier functions for this spin
  nwan=WData%nwan(isp) 
  ham_r => WData%hamWR(1:nwan,1:nwan,1:WData%nrpts,isp)

  ! TODO add cut case at this level !if (cut) then ...

  ABI_ALLOCATE(eig_interp,(nwan,BSt%nkpt))
  eig_interp=zero 
  !allocate(bands_proj(nwan,BSt%nkpt),stat=ierr)
  !bands_proj=zero

  ABI_ALLOCATE(ham_pack,((nwan*(nwan+1))/2))
  ierr = ABI_ALLOC_STAT
  ABI_ALLOCATE(ham_kprm,(nwan,nwan))
  ierr = ABI_ALLOC_STAT
  ABI_ALLOCATE(U_int,(nwan,nwan))
  ierr = ABI_ALLOC_STAT
  ABI_ALLOCATE(cwork,(2*nwan))
  ABI_ALLOCATE(rwork,(7*nwan))
  ABI_ALLOCATE(iwork,(5*nwan))
  ABI_ALLOCATE(ifail,(nwan))
  !
  ! Cut H matrix in real-space
  !if (index(bands_plot_mode,'cut').ne.0)  call plot_cut_hr()

  ! === Interpolate the Hamiltonian at each kpoint ===
  do ikpt=1,BSt%nkpt

   ham_kprm(:,:)=czero
   !if (index(bands_plot_mode,'s-k').ne.0) then
    do iR=1,nrpts
     rdotk=two_pi*DOT_PRODUCT(BSt%kptns(:,ikpt),irvec(:,iR))
     fact=EXP(j_dpc*rdotk)/REAL(ndegen(iR),dp)
     ham_kprm=ham_kprm+fact*ham_r(:,:,iR)
    end do
   !TODO
   !else if (index(bands_plot_mode,'cut').ne.0) then
   ! do iR=1,nrpts_cut
   !  rdotk=two_pi*DOT_PRODUCT(BSt%kptns(:,ikpt),irvec_cut(:,iR))
   !  !!$[aam] check divide by ndegen?
   !  fact=EXP(j_dpc*rdotk)
   !  ham_kprm=ham_kprm+fact*ham_r_cut(:,:,iR)
   ! end do
   !end if

   ! === Diagonalise H_k (->basis of eigenstates) ===
   do jj=1,nwan
    do ii=1,jj
     ham_pack(ii+((jj-1)*jj)/2)=ham_kprm(ii,jj)
    end do
   end do

   ABSTOL=-one !; ABSTOL= 2*DLAMCH('S')

   call ZHPEVX('V','A','U',nwan,ham_pack,zero,zero,0,0,ABSTOL,&
    nfound,eig_interp(:,ikpt),U_int,nwan,cwork,rwork,iwork,ifail,info)

   if (info<0) then
    write(msg,'(a,i3,a)') 'THE ',-info, ' ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE'
    MSG_ERROR(msg)
   end if
   if (info>0) then
    write(msg,'(i3,a)')info,' EIGENVECTORS FAILED TO CONVERGE'
    MSG_ERROR(msg)
   end if

   ! Compute projection onto WF if requested
   !if (num_bands_project>0) then
   ! do loop_w=1,nwan
   !  do loop_p=1,nwan 
   !   if (any(bands_plot_project==loop_p)) then
   !    bands_proj(loop_w,ikpt)=bands_proj(loop_w,ikpt)+abs(U_int(loop_p,loop_w))**2
   !   end if 
   !  end do
   ! end do
   !end if

   ! * eig_interp contains the interpolated nwan bands for this k-points
   ! TODO here be careful since first index might refer to an arbitrary band!
   BSt%eig(1:nwan,ikpt,isp)=eig_interp(:,ikpt)
   !write(97,'(7f8.4)')BSt%kptns(:,ikpt),eig_interp(1:nwan,ikpt)
  end do !ikpt

  ! Interpolation Finished.
  ABI_DEALLOCATE(eig_interp)
  ABI_DEALLOCATE(ham_pack)
  ABI_DEALLOCATE(ham_kprm)
  ABI_DEALLOCATE(U_int)
  ABI_DEALLOCATE(cwork)
  ABI_DEALLOCATE(rwork)
  ABI_DEALLOCATE(iwork)
  ABI_DEALLOCATE(ifail)

  !if (allocated(ham_r_cut)) deallocate(ham_r_cut,stat=ierr)
  !if (allocated(irvec_cut)) deallocate(irvec_cut,stat=ierr)
 end do !isp

end subroutine WannierInterpol

END MODULE m_wannier2abinit
!!***
