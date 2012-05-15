!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ebands
!! NAME
!!  m_ebands
!!
!! FUNCTION
!!  This module contains utilities to analyze and retrieve information
!!  from the bandstructure_type.
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

MODULE m_ebands

 use m_profiling

 use defs_basis
 use defs_datatypes, only : bandstructure_type
 use defs_abitypes,  only : hdr_type
 use m_errors
#if defined HAVE_TRIO_ETSF_IO
 use etsf_io
#endif

 use m_copy,           only : deep_copy
 use m_fstrings,       only : tolower
 use m_io_tools,       only : get_unit
 use m_numeric_tools,  only : arth, imin_loc, imax_loc
 use m_geometry,       only : normv
 use m_crystal,        only : crystal_structure
 use m_header,         only : hdr_clean, hdr_get_nelect_byocc
 use m_bz_mesh,        only : bz_mesh_type, init_kmesh, find_qmesh, isamek, get_BZ_item, get_BZ_diff, destroy_BZ_mesh_type,&
&                             little_group, nullify_little_group, destroy_little_group

 implicit none

 private

 public :: bstruct_init            ! Main creation method.
 public :: bst_init_from_hdr       ! Init from abinit header.
 public :: bstruct_clean           ! Destruction method.
 public :: copy_bandstructure      ! Deep copy of the bandstructure_type.
 public :: print_bandstructure     ! Reports basic info on the data type.
 public :: unpack_eneocc           ! Helper function for reshaping (energies|occupancies|derivate of occupancies).
 public :: pack_eneocc             ! Helper function for reshaping (energies|occupancies|derivate of occupancies).
 public :: get_eneocc_vect         ! Reshape (ene|occ|docdde) returning a matrix instead of a vector.
 public :: put_eneocc_vect         ! Put (ene|occ|doccde) in vectorial form into the data type doing a reshape.
 public :: get_bandenergy          ! Returns the band energy of the system.
 public :: get_valence_idx         ! Gives the index of the (valence|bands at E_f).
 public :: apply_scissor           ! Apply a scissor operator (no k-dependency)
 public :: get_occupied            ! Returns band indeces after wich occupations are less than an input value.
 public :: enclose_degbands        ! Adjust band indeces such that all degenerate states are treated.
 public :: get_minmax              ! Returns min and Max value of (eig|occ|doccde).
 public :: get_FS                  ! Returns the Fermi surface (k-points and bands at E_f).
 public :: bst_ismetal             ! .TRUE. if metallic occupation scheme is used.
 public :: bst_print_fs            ! Dumps file to be used for Fermi surface visualization. TODO
 public :: get_dos                 ! Calculates the electronic DOS using (gaussian|tetrahedrons).
 public :: update_occ              ! Update the occupation numbers.
 public :: ReportGap               ! Print info on the fundamental and optical gap.
 public :: bst_plot_bands          ! Plot bands in the XMGRACE format.
 public :: SelectBands             ! Extracts a subset of bands from the structure.
 public :: ExpandBands             ! Returns a new object defined in the BZ starting from the IBZ.
 public :: joint_dos               ! Calculate the joint density of states.
 public :: abi_bands_put           ! Dump the object into NETCDF file.
 public :: abi_Bands_read          ! Initialize the object from a NETCDF file.

CONTAINS  !=========================================================================================================================
!!***

!!****f* m_ebands/bstruct_init
!! NAME
!! bstruct_init
!!
!! FUNCTION
!! This subroutine initializes the bandstructure structured datatype
!!
!! INPUTS
!! bantot=total number of bands (=sum(nband(:))
!! doccde(bantot)=derivative of the occupation numbers with respect to the energy (Ha)
!! eig(bantot)=eigenvalues (hartree)
!! istwfk(nkpt)=parameter that describes the storage of wfs.
!! kptns(3,nkpt)=k points in terms of recip primitive translations
!! nband(nkpt*nsppol)=number of bands
!! nelect=Number of electrons.
!! nkpt=number of k points
!! npwarr(nkpt)=number of planewaves at each k point
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nspinor=Number of spinor components
!! occopt=Occupation options (see input variable)
!! occ(bantot)=occupation numbers
!! tphysel=Physical temperature (input variable)
!! tsmear=Temperature of smearing.
!! wtk(nkpt)=weight assigned to each k point
!!
!! OUTPUT
!! bstruct<bstruct_type>=the bandstructure datatype
!!
!! SIDE EFFECTS
!!  %entropy and %fermie initialized to zero.
!!
!! PARENTS
!!      gstate,loper3,m_ebands,m_gwannier,m_shirley,mlwfovlp_qp,newsp,nonlinear
!!      respfn,setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine bstruct_init(bantot,bstruct,nelect,doccde,eig,istwfk,kptns,&
& nband,nkpt,npwarr,nsppol,nspinor,tphysel,tsmear,occopt,occ,wtk)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bstruct_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot,nkpt,nsppol,nspinor,occopt
 real(dp),intent(in) :: nelect,tphysel,tsmear
 type(bandstructure_type),intent(out) :: bstruct
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt)
 real(dp),intent(in) :: doccde(bantot),eig(bantot),kptns(3,nkpt),occ(bantot)
 real(dp),intent(in) :: wtk(nkpt)
! *************************************************************************

 ! Copy the scalars
 ! MG TODO here there is a inconsistency in the way occ are treated in the header
 ! (only the states used, bantot. are saved, and the way occ. and energies
 ! are passed to routines (mband,nkpt,nsppol). It might happen that bantot<mband*nktp*nsppol
 ! this should not lead to problems since arrays are passed by reference
 ! anyway the treatment of these arrays have to be rationalized
 bstruct%bantot =bantot
 bstruct%mband  =MAXVAL(nband(1:nkpt*nsppol))
 bstruct%nkpt   =nkpt
 bstruct%nspinor=nspinor
 bstruct%nsppol =nsppol
 bstruct%occopt =occopt
 
 bstruct%entropy=zero  ! Initialize results
 bstruct%fermie =zero  ! Initialize results
 bstruct%nelect =nelect
 bstruct%tphysel=tphysel
 bstruct%tsmear =tsmear

 ! Allocate the components
 ABI_ALLOCATE(bstruct%nband,(nkpt*nsppol))
 ABI_ALLOCATE(bstruct%istwfk,(nkpt))
 ABI_ALLOCATE(bstruct%npwarr,(nkpt))
 ABI_ALLOCATE(bstruct%kptns,(3,nkpt))

 ! Copy the arrays
 bstruct%nband(1:nkpt*nsppol)=nband(1:nkpt*nsppol)
 bstruct%istwfk(1:nkpt)      =istwfk(1:nkpt)
 bstruct%npwarr(1:nkpt)      =npwarr(1:nkpt)
 bstruct%kptns(1:3,1:nkpt)   =kptns(1:3,1:nkpt)

 ! In bstruct, energies and occupations are stored in a matrix (mband,nkpt,nsppol).
 ! put_eneocc_vect is used to reshape the values stored in vectorial form.
 ABI_ALLOCATE(bstruct%eig   ,(bstruct%mband,nkpt,nsppol))
 bstruct%eig   = HUGE(one)
 ABI_ALLOCATE(bstruct%occ   ,(bstruct%mband,nkpt,nsppol))
 bstruct%occ   =zero        
 ABI_ALLOCATE(bstruct%doccde,(bstruct%mband,nkpt,nsppol))
 bstruct%doccde=zero

 call put_eneocc_vect(bstruct,'eig',   eig   ) 
 call put_eneocc_vect(bstruct,'occ',   occ   ) 
 call put_eneocc_vect(bstruct,'doccde',doccde) 

 ABI_ALLOCATE(bstruct%wtk,(nkpt))
 bstruct%wtk(1:nkpt)=wtk(1:nkpt)

end subroutine bstruct_init
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/bst_init_from_hdr
!! NAME
!! bst_init_from_hdr
!!
!! FUNCTION
!! This subroutine initializes the bandstructure datatype from the abinit header by
!! calling the main creation method.
!!
!! INPUTS
!!  Hdr<hdr_type>=Abinit header.
!!  mband=Maximun number of bands.
!!  ene3d(mband,Hdr%nkpt,Hdr%nsppol)=Energies.
!!  [nelect]=Number of electrons per unit cell.
!!    Optional argument that can be used for performing a ridid shift of the fermi level. 
!!    in the case of metallic occupancies.
!!    If not specified, nelect will be initialized from Hdr.
!!
!! OUTPUT
!!  Bst<bstruct_type>=The bandstructure datatype completely initialized.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!
!! SOURCE

subroutine bst_init_from_hdr(Bst,Hdr,mband,ene3d,nelect)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bst_init_from_hdr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband
 type(hdr_type),intent(in) :: Hdr
 type(bandstructure_type),intent(out) :: Bst
 real(dp),optional,intent(in) :: nelect
!arrays
 real(dp),intent(in) :: ene3d(mband,Hdr%nkpt,Hdr%nsppol)

!Local variables-------------------------------
!scalars
 real(dp) :: my_nelect
!arrays
 real(dp),allocatable :: ugly_doccde(:),ugly_ene(:)
! *************************************************************************

 if (PRESENT(nelect)) then
   my_nelect = nelect
 else 
   ! TODO Have to add nelect to the header
   my_nelect = hdr_get_nelect_byocc(Hdr) 
 end if
 !
 ! Have to use ugly 1d vectors to call bstruct_init
 ABI_ALLOCATE(ugly_doccde,(Hdr%bantot))
 ugly_doccde=zero

 ABI_ALLOCATE(ugly_ene,(Hdr%bantot))
 call pack_eneocc(Hdr%nkpt,Hdr%nsppol,mband,Hdr%nband,Hdr%bantot,ene3d,ugly_ene)

 call bstruct_init(Hdr%bantot,Bst,my_nelect,ugly_doccde,ugly_ene,Hdr%istwfk,Hdr%kptns,Hdr%nband,Hdr%nkpt,&
&  Hdr%npwarr,Hdr%nsppol,Hdr%nspinor,Hdr%tphysel,Hdr%tsmear,Hdr%occopt,Hdr%occ,Hdr%wtk)

 ABI_DEALLOCATE(ugly_doccde)
 ABI_DEALLOCATE(ugly_ene)

end subroutine bst_init_from_hdr
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/bstruct_clean
!! NAME
!! bstruct_clean
!!
!! FUNCTION
!! Deallocates the components of the bandstructure structured datatype
!!
!! INPUTS
!!  bstruct<bandstructure_type>=The data type to be deallocated.
!!
!! OUTPUT
!!  Deallocate the dynamic arrays in the bandstructure type.
!!  (only deallocate)
!!
!! PARENTS
!!      bethe_salpeter,elphon,exc_interp_ham,gstate,loper3,m_gwannier
!!      mlwfovlp_qp,newsp,nonlinear,respfn,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine bstruct_clean(bstruct)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bstruct_clean'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(bandstructure_type),intent(inout) :: bstruct
! *************************************************************************

 DBG_ENTER("COLL")

!Deallocate all components of bstruct
 ABI_DEALLOCATE(bstruct%istwfk)
 ABI_DEALLOCATE(bstruct%nband)
 ABI_DEALLOCATE(bstruct%npwarr)
 ABI_DEALLOCATE(bstruct%kptns)
 ABI_DEALLOCATE(bstruct%eig)
 ABI_DEALLOCATE(bstruct%occ)
 ABI_DEALLOCATE(bstruct%doccde)
 ABI_DEALLOCATE(bstruct%wtk)

 DBG_EXIT("COLL")

end subroutine bstruct_clean
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/copy_bandstructure
!! NAME
!!  copy_bandstructure
!!
!! FUNCTION
!! This subroutine performs a deep copy of a bandstructure datatype.
!! All the associated pointers in the input object will be copied preserving the shape.
!! If a pointer in BSt_in happens to be not associated, the corresponding
!! pointer in the copied object will be nullified.
!!
!! INPUTS
!!  BSt_in<bandstructure_type>=The data type to be copied.
!!
!! OUTPUT
!!  BSt_cp<bandstructure_type>=The copy.
!!
!! TODO 
!!  To be on the safe side one should nullify all pointers in the bandstructure_type 
!!  in the creation method. We have to follow F90 specifications and the initial status 
!!  of a pointer is not defined. This might lead to problem in deep_copy.
!!
!! PARENTS
!!      m_gwannier,screening,setup_bse,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_bandstructure(BSt_in,BSt_cp)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_bandstructure'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(bandstructure_type),intent(in)  :: BSt_in
 type(bandstructure_type),intent(out) :: BSt_cp

! *********************************************************************

 ! === Copy scalars ===
 BSt_cp%bantot  = BSt_in%bantot  
 BSt_cp%mband   = BSt_in%mband
 BSt_cp%nkpt    = BSt_in%nkpt    
 BSt_cp%nspinor = BSt_in%nspinor 
 BSt_cp%nsppol  = BSt_in%nsppol  
 BSt_cp%occopt  = BSt_in%occopt  

 BSt_cp%entropy = BSt_in%entropy 
 BSt_cp%fermie  = BSt_in%fermie  
 BSt_cp%nelect  = BSt_in%nelect  
 BSt_cp%tphysel = BSt_in%tphysel 
 BSt_cp%tsmear  = BSt_in%tsmear  

 ! === Copy pointers ===
 call deep_copy( BSt_in%istwfk, BSt_cp%istwfk)
 call deep_copy( BSt_in%nband , BSt_cp%nband )     
 call deep_copy( BSt_in%npwarr, BSt_cp%npwarr)    

 call deep_copy( BSt_in%kptns , BSt_cp%kptns ) 
 call deep_copy( BSt_in%eig   , BSt_cp%eig   )  
 call deep_copy( BSt_in%occ   , BSt_cp%occ   )   
 call deep_copy( BSt_in%doccde, BSt_cp%doccde)   
 call deep_copy( BSt_in%wtk   , BSt_cp%wtk   )

end subroutine copy_bandstructure
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/print_bandstructure
!! NAME
!! print_bandstructure
!!
!! FUNCTION
!! Print the content of the object.
!!
!! INPUTS
!!  BSt<bandstructure_type>The type containing the data.
!!  [unit]=Unit number (std_out if None)
!!  [header]=title for info
!!  [prtvol]=Verbosity level (0 if None)
!!  [mode_paral]=Either 'COLL' or 'PERS' ('COLL' if None).
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_gwannier,setup_bse
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_bandstructure(BSt,header,unit,prtvol,mode_paral)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_bandstructure'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=*),optional,intent(in) :: header
 character(len=4),optional,intent(in) :: mode_paral
 type(bandstructure_type),intent(in) :: BSt

!Local variables-------------------------------
 integer :: isppol,ikpt,my_unt,my_prtvol,ii
 character(len=4) :: my_mode
 character(len=500) :: msg
! *************************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the bandstructure_type ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(5(a,i5,a))')&
&  '  Number of spinorial components ...... ',BSt%nspinor,ch10,&
&  '  Number of spin polarizations ........ ',BSt%nsppol,ch10,&
&  '  Number of k-points in the IBZ ....... ',BSt%nkpt,ch10,&
&  '  Maximum number of bands ............. ',BSt%mband,ch10,&
&  '  Occupation option ................... ',BSt%occopt,ch10
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(a,f14.2,a,4(a,f14.6,a))')&
&  '  Number of valence electrons ......... ',BSt%nelect,ch10,&
&  '  Fermi level  ........................ ',BSt%fermie,ch10,&
&  '  Entropy ............................. ',BSt%entropy,ch10,&
&  '  Tsmear value ........................ ',BSt%tsmear,ch10,&
&  '  Tphysel value ....................... ',BSt%tphysel,ch10
 call wrtout(my_unt,msg,my_mode)

 if (my_prtvol>0) then

   if (BSt%nsppol==1)then
     write(msg,'(a,i4,a)')' New occ. numbers for occopt= ',BSt%occopt,' , spin-unpolarized case. '
     call wrtout(my_unt,msg,my_mode)
   end if

   do isppol=1,BSt%nsppol
     if (BSt%nsppol==2) then
       write(msg,'(a,i4,a,i2)')' New occ. numbers for occopt= ',BSt%occopt,' spin ',isppol
       call wrtout(my_unt,msg,my_mode)
     end if

     do ikpt=1,BSt%nkpt
       write(msg,'(2a,i4,a,3f12.6,a,f6.3)')ch10,&
&        ' k-point number ',ikpt,') ',BSt%kptns(:,ikpt),'; weight: ',BSt%wtk(ikpt)
       call wrtout(my_unt,msg,my_mode)
       do ii=1,BSt%nband(ikpt+(isppol-1)*BSt%nkpt)
         write(msg,'(3(f7.3,1x))')BSt%eig(ii,ikpt,isppol)*Ha_eV,BSt%occ(ii,ikpt,isppol),BSt%doccde(ii,ikpt,isppol)
         call wrtout(my_unt,msg,my_mode)
       end do
     end do !ikpt

   end do !isppol

   !TODO add additional info useful for debugging)
   !istwfk(:), nband(:)
 end if !my_prtvol

end subroutine print_bandstructure
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/unpack_eneocc
!! NAME
!! unpack_eneocc
!!
!! FUNCTION
!!  Helper function to do a reshape of (energies|occupancies|derivate of occupancies)
!!  initially stored in a vector. Return a 3D array index by (iband,ikpt,isppol) 
!!
!! INPUTS
!!  nkpt=number of k-points
!!  nsppol=number of spin polarizations
!!  mband=Max number of bands over k-points (just to dimension the output)
!!  nbands(nkpt*nsppol)=Number of bands at eack k and spin
!!  bantot=Total number of bands
!!  vect(bantot)=The input values to reshape
!!
!! OUTPUT
!!  array3d(mband,nkpt,nsppol)=Arrays containing the values of vect. 
!!   Note that the first dimension is usually larger than the 
!!   number of bands really used for a particular k-point and spin.
!!
!! TODO 
!!  Switch to function just to allow inlining
!!
!! PARENTS
!!      cchi0q0_intraband,kss2wfk,m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine unpack_eneocc(nkpt,nsppol,mband,nband,bantot,vect,array3d)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unpack_eneocc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsppol,mband,bantot
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: vect(bantot)
 real(dp),intent(out) :: array3d(mband,nkpt,nsppol)

!Local variables-------------------------------
 integer :: isppol,ikpt,iband,idx
! *************************************************************************

 array3d(:,:,:)=zero
 idx=0
 ! elements in vect are packed in the first positions.
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do iband=1,nband(ikpt+(isppol-1)*nkpt)
      idx=idx+1
      array3d(iband,ikpt,isppol)=vect(idx)
     end do
   end do
 end do

end subroutine unpack_eneocc
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/pack_eneocc
!! NAME
!! pack_eneocc
!!
!! FUNCTION
!!  Helper function to do a reshape of (energies|occupancies|derivate of occupancies)
!!  initially stored in a 3D arrays returning a vector. 
!!
!! INPUTS
!!  nkpt=number of k-points
!!  nsppol=number of spin polarizations
!!  mband=Max number of bands over k-points (just to dimension the output)
!!  nbands(nkpt*nsppol)=Number of bands at eack k and spin
!!  bantot=Total number of bands
!!  array3d(mband,nkpt,nsppol)=Arrays containing the values to reshape.
!!
!! OUTPUT
!!  vect(bantot)=The input values stored in vector mode. Only the values really
!!   considered at each k-point and spin are copied.
!!
!! PARENTS
!!      cchi0q0_intraband,m_ebands,m_shirley
!!
!! CHILDREN
!!
!! SOURCE

subroutine pack_eneocc(nkpt,nsppol,mband,nband,bantot,array3d,vect)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pack_eneocc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsppol,mband,bantot
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: array3d(mband,nkpt,nsppol)
 real(dp),intent(out) :: vect(bantot)

!Local variables-------------------------------
 integer :: isppol,ikpt,iband,idx

! *************************************************************************

 vect(:)=zero
 idx=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do iband=1,nband(ikpt+(isppol-1)*nkpt)
       idx=idx+1
       vect(idx)=array3d(iband,ikpt,isppol)
     end do
   end do
 end do

end subroutine pack_eneocc 
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_eneocc_vect
!! NAME
!! get_eneocc_vect
!!
!! FUNCTION
!!  Retrieve energies or occupations from a Bandstructure structure accessing by name. 
!!  Results are reported in a vector to facilitate the interface with other abinit routines.
!!
!! INPUTS
!!  BSt<bandstructure_type>The type containing the data.
!!  arr_name=The name of the quantity to retrieve. Allowed values are
!!   == "eig"    == For the eigenvalues.
!!   == "occ"    == For the occupation numbers.
!!   == "doccde" == For the derivative of the occupancies wrt the energy.
!!
!! OUTPUT
!!  vect(BSt%bantot)=The values required.
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_eneocc_vect(BSt,arr_name,vect)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_eneocc_vect'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arr_name
 type(bandstructure_type),intent(in) :: BSt
 real(dp),intent(out) :: vect(BSt%bantot)

!Local variables-------------------------------
 integer :: nkpt,nsppol,mband,bantot
 character(len=500) :: msg
! *************************************************************************

 mband =BSt%mband 
 bantot=BSt%bantot
 nkpt  =BSt%nkpt
 nsppol=BSt%nsppol

 SELECT CASE (arr_name)
 CASE ('occ')
  call pack_eneocc(nkpt,nsppol,mband,BSt%nband,bantot,BSt%occ,vect)
 CASE ('eig')
  call pack_eneocc(nkpt,nsppol,mband,BSt%nband,bantot,BSt%eig,vect)
 CASE ('doccde')
  call pack_eneocc(nkpt,nsppol,mband,BSt%nband,bantot,BSt%doccde,vect)
 CASE DEFAULT
  write(msg,'(2a)')' Wrong value of arr_name= ',TRIM(arr_name)
  MSG_BUG(msg)
 END SELECT

end subroutine get_eneocc_vect
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/put_eneocc_vect
!! NAME
!! put_eneocc_vect
!!
!! FUNCTION
!!  Update the energies or the occupations stored a Bandstructure structure. 
!!  The input values are stored in a vector according to the abinit convention
!!  In the data type, on the contrary,  we use 3D arrays (mband,nkpt,nsspol) 
!!  which are much easier to use inside loops.
!!
!! INPUTS
!!  vect(BSt%bantot)=The new values to be stored in the structure.
!!  arr_name=The name of the quantity to be saved (CASE insensitive). 
!!  Allowed values are
!!   == "eig"    == For the eigenvalues.
!!   == "occ"    == For the occupation numbers.
!!   == "doccde" == For the derivative of the occupancies wrt the energy.
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  BSt<bandstructure_type>=The object with updated values depending on the value of arr_name
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine put_eneocc_vect(BSt,arr_name,vect) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'put_eneocc_vect'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arr_name
 type(bandstructure_type),intent(inout) :: BSt
 real(dp),intent(in) :: vect(BSt%bantot)

!Local variables-------------------------------
 integer :: nkpt,nsppol,mband,bantot
 character(len=500) :: msg
! *************************************************************************

 mband =BSt%mband 
 bantot=BSt%bantot
 nkpt  =BSt%nkpt
 nsppol=BSt%nsppol

 SELECT CASE (tolower(arr_name))
 CASE ('occ')
  call unpack_eneocc(nkpt,nsppol,mband,BSt%nband,bantot,vect,BSt%occ)
 CASE ('eig')
  call unpack_eneocc(nkpt,nsppol,mband,BSt%nband,bantot,vect,BSt%eig)
 CASE ('doccde')
  call unpack_eneocc(nkpt,nsppol,mband,BSt%nband,bantot,vect,BSt%doccde)
 CASE DEFAULT 
  write(msg,'(2a)')' Wrong value of arr_name= ',TRIM(arr_name)
  MSG_BUG(msg)
 END SELECT

end subroutine put_eneocc_vect
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_bandenergy
!! NAME
!! get_bandenergy
!!
!! FUNCTION
!!  Return the band energy (weighted sum of occupied eigenvalues)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!! TODO Likely this expression is not accurate since it is not variatonal
!!  One should use 
!!   band_energy = \int e N(e) de   for e<Ef , where N(e) is the e-DOS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_bandenergy(BSt) result(band_energy)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_bandenergy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(bandstructure_type),intent(in) :: BSt
 real(dp) :: band_energy

!Local variables-------------------------------
 integer :: isppol,ikibz,nband_k
 real(dp) :: wtk
! *********************************************************************

 band_energy=zero
 do isppol=1,BSt%nsppol
   do ikibz=1,BSt%nkpt
     wtk=BSt%wtk(ikibz)
     nband_k=BSt%nband(ikibz+(isppol-1)*BSt%nkpt) 
     band_energy = band_energy + wtk*SUM( BSt%eig(1:nband_k,ikibz,isppol)*BSt%occ(1:nband_k,ikibz,isppol) )
   end do
 end do

end function get_bandenergy
!!***

!!****f* m_ebands/get_valence_idx
!! NAME
!!  get_valence_idx
!!
!! FUNCTION
!!  For each k-point and spin polarisation, report: 
!!   The index of the valence in case of Semiconductors.
!!   The index of the band at the Fermi energy+toldfe
!!
!! INPUTS
!!  BSt<bandstructure_type>=The object describing the band structure.
!!  tol_fermi[optional]
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_valence_idx(BSt,tol_fermi) result(val_idx)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_valence_idx'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: tol_fermi
 type(bandstructure_type),intent(in) :: BSt
!arrays
 integer :: val_idx(BSt%nkpt,BSt%nsppol)

!Local variables-------------------------------
 integer :: iband,ikpt,isppol,idx,nband_k 
 real(dp) :: tol_

! *************************************************************************

 tol_=tol6 ; if (PRESENT(tol_fermi)) tol_=tol_fermi

 do isppol=1,BSt%nsppol
   do ikpt=1,BSt%nkpt
     nband_k=BSt%nband(ikpt+(isppol-1)*BSt%nkpt)

     idx=0
     do iband=1,nband_k
       if (BSt%eig(iband,ikpt,isppol) > BSt%fermie+ABS(tol_)) then
         idx=iband ; EXIT
       end if
     end do
     val_idx(ikpt,isppol)=idx-1
     if (idx==1) val_idx(ikpt,isppol)=idx
     if (idx==0) val_idx(ikpt,isppol)=nband_k

   end do
 end do

end function get_valence_idx
!!***

!!****f* m_ebands/apply_scissor
!! NAME
!!  apply_scissor
!!
!! FUNCTION
!!  Apply a scissor operator of amplitude scissor_energy.
!!
!! INPUTS
!!  scissor_energy=The energy shift
!!
!! OUTPUT
!!
!! SIDE EFFECT
!!  BSt<bandstructure_type>=The following quantities are modified:
!!   %eig(mband,nkpt,nsppol)=The band structure after the application of the scissor operator
!!   %fermi_energy
!!
!! PARENTS
!!      screening,setup_bse
!!
!! CHILDREN
!!
!! SOURCE

subroutine apply_scissor(BSt,scissor_energy)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'apply_scissor'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: scissor_energy
 type(bandstructure_type),intent(inout) :: BSt

!Local variables-------------------------------
 integer :: ikpt,isppol,ival,nband_k 
 real(dp) :: fixmom_
 character(len=500) :: msg
!arrays
 integer :: val_idx(BSt%nkpt,BSt%nsppol)
! *************************************************************************

 ! === Get the valence band index for each k and spin ===
 val_idx(:,:) = get_valence_idx(BSt)

 do isppol=1,BSt%nsppol
  if (ANY(val_idx(:,isppol)/=val_idx(1,isppol))) then
   write(msg,'(a,i2,a)')&
&   ' Trying to apply a scissor operator on a metallic band structure for spin = ',isppol,&
&   ' Assuming you know what you are doing, continuing anyway! '
   MSG_COMMENT(msg)
   !Likely newocc will stop, unless the system is semimetallic ?
  end if
 end do

 ! === Apply the scissor ===
 do isppol=1,BSt%nsppol
  do ikpt=1,BSt%nkpt
   nband_k=BSt%nband(ikpt+(isppol-1)*BSt%nkpt)
   ival=val_idx(ikpt,isppol)

   if (nband_k>=ival+1) then
    BSt%eig(ival+1:,ikpt,isppol) = BSt%eig(ival+1:,ikpt,isppol)+scissor_energy
   else 
    write(msg,'(2a,4(a,i4))')&
&    ' Not enough bands to apply the scissor operator. ',ch10,&
&    ' spin = ',isppol,' ikpt = ',ikpt,' nband_k = ',nband_k,' but valence index = ',ival
    MSG_COMMENT(msg)
   end if

  end do
 end do

 ! === Recalculate the fermi level and occ. factors ===
 ! * For Semiconductors only the Fermi level is changed (in the middle of the new gap) 
 fixmom_=-99.99_dp !?; if (PRESENT(fixmom)) fixmom_=fixmom
 call update_occ(BSt,fixmom_)

end subroutine apply_scissor
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_occupied
!! NAME
!!  get_occupied
!!
!! FUNCTION
!!  For each k-point and spin polarisation, report the band index
!!  after which the occupation numbers are less than tol_occ.
!!
!! INPUTS
!!  BSt<bandstructure_type>=The object describing the band structure.
!!  tol_occ[Optional]=Tollerance on the occupation factors.
!!
!! OUTPUT
!!
!! NOTES
!!  We assume that the occupation factors monotonically decrease as a function of energy.
!!  This is not always true for eavery smearing technique implemented in Abinit.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_occupied(BSt,tol_occ) result(occ_idx)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_occupied'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: tol_occ
 type(bandstructure_type),intent(in) :: BSt
!arrays
 integer :: occ_idx(BSt%nkpt,BSt%nsppol)

!Local variables-------------------------------
 integer :: iband,ikpt,isppol,idx,nband_k 
 real(dp) :: tol_

! *************************************************************************

 tol_=tol8 ; if (PRESENT(tol_occ)) tol_=tol_occ

 do isppol=1,BSt%nsppol
  do ikpt=1,BSt%nkpt
   nband_k=BSt%nband(ikpt+(isppol-1)*BSt%nkpt)

   idx=0
   do iband=1,nband_k
    if (BSt%occ(iband,ikpt,isppol)<ABS(tol_)) then
     idx=iband; EXIT
    end if
   end do
   occ_idx(ikpt,isppol)=idx-1
   if (idx==1) occ_idx(ikpt,isppol)=idx
   if (idx==0) occ_idx(ikpt,isppol)=nband_k

  end do
 end do

end function get_occupied
!!***

!!****f* m_ebands/enclose_degbands
!! NAME
!!  enclose_degbands
!!
!! FUNCTION
!!  Adjust ibmin and ibmax such that all the degenerate states are enclosed
!!  between ibmin and ibmax. The routine works for a given k-point a spin.
!!
!! INPUTS
!!  BSt<bandstructure_type>=The object describing the band structure.
!!  ikibz=Index of the k-point.
!!  isppol=Spin index.
!!  tol_enedif=Tolerance on the energy difference.
!!
!! OUTPUT 
!!  changed=.TRUE. if ibmin or ibmax has been changed.
!!
!! SIDE EFFECTS
!!  ibmin,ibmax=
!!    Input: initial guess for the indeces
!!    Output: All the denerate states are between ibmin and ibmax 
!!
!! PARENTS
!!      setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine enclose_degbands(BSt,ikibz,isppol,ibmin,ibmax,changed,tol_enedif)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'enclose_degbands'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikibz,isppol
 integer,intent(inout) :: ibmin,ibmax
 real(dp),intent(in) :: tol_enedif
 logical,intent(out) :: changed
 type(bandstructure_type),intent(in) :: BSt

!Local variables-------------------------------
!scalars
 integer :: ib,ibmin_bkp,ibmax_bkp
 real(dp) :: emin,emax

! *************************************************************************

 ibmin_bkp = ibmin
 ibmax_bkp = ibmax

 emin =  BSt%eig(ibmin,ikibz,isppol)
 do ib=ibmin-1,1,-1
   if ( ABS(BSt%eig(ib,ikibz,isppol) - emin) > tol_enedif) then
     ibmin = ib +1 
     EXIT 
   else 
     ibmin = ib
   end if
 end do

 emax =  BSt%eig(ibmax,ikibz,isppol)
 do ib=ibmax+1,BSt%nband(ikibz+(isppol-1)*BSt%nkpt)
   if ( ABS(BSt%eig(ib,ikibz,isppol) - emax) > tol_enedif) then
     ibmax = ib - 1 
     EXIT 
   else 
     ibmax = ib
   end if
 end do

 changed = (ibmin /= ibmin_bkp) .or. (ibmax /= ibmax_bkp)

end subroutine enclose_degbands
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_minmax
!! NAME
!!  get_minmax
!!
!! FUNCTION
!!  Report the min and max value over k-points and bands of (eig|occ|doccde) for each 
!!  spin. Cannot use F90 array syntax due to the internal storage used in abinit.
!!
!! INPUTS
!!  BSt<bandstructure_type>=The object describing the band structure.
!!  arr_name=The name of the array whose min and Max value has to be calculated.
!!   Possible values: 'occ', 'eig' 'doccde'
!!
!! OUTPUT
!! minmax(2,BSt%nsppol)=For each spin the min and max value of the quantity specified by "arr_name"
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_minmax(BSt,arr_name) result(minmax)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_minmax'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(bandstructure_type),intent(in) :: BSt
 character(len=*),intent(in) :: arr_name
!arrays
 real(dp) :: minmax(2,BSt%nsppol)

!Local variables-------------------------------
!scalars
 integer :: iband,ikpt,isppol,nband_k 
 real(dp) :: datum
 character(len=500) :: msg
!arrays
 real(dp),pointer :: rdata(:,:,:)

! *************************************************************************

 SELECT CASE (tolower(arr_name))
 CASE ('occ')
  rdata => BSt%occ
 CASE ('eig')
  rdata => BSt%eig
 CASE ('doccde')
  rdata => BSt%doccde
 CASE DEFAULT
  write(msg,'(2a)')' Wrong value of arr_name = ',TRIM(arr_name)
  MSG_BUG(msg)
 END SELECT

 minmax(1,:)=greatest_real
 minmax(2,:)=smallest_real
 
 do isppol=1,BSt%nsppol
  do ikpt=1,BSt%nkpt
   nband_k=BSt%nband(ikpt+(isppol-1)*BSt%nkpt)
   do iband=1,nband_k
    datum=rdata(iband,ikpt,isppol)
    minmax(1,isppol)=MIN(minmax(1,isppol),datum)
    minmax(2,isppol)=MAX(minmax(2,isppol),datum)
   end do
  end do
 end do

end function get_minmax
!!***

!!****f* m_ebands/get_FS
!! NAME
!!  get_FS
!!
!! FUNCTION
!!  Returns the indeces of the k-points belonging to the Fermi surface as well
!!  as the minimum and Maximum index of the bands crossing the Fermi level.
!!
!! INPUTS
!!  BSt<bandstructure_type>=The object describing the band structure.
!!  tolweight=Tolerange on the weighs for FS integrations. A k-point belongs
!!    to the Fermi surface if for some band its weight is > tolweight
!!
!! OUTPUT
!!  nFSkpt=Number of points on the Fermi surface.
!!  fs2ibz(1:nFSkpt)=Index of the FS points in the input IBZ array.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_FS(BSt,isppol,tolweight,fs_weight,bmin,bmax,nFSkpt,fs2ibz)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_FS'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol
 integer,intent(out) :: bmin,bmax
 real(dp),intent(in) :: tolweight
 type(bandstructure_type),intent(in) :: BSt
 integer,intent(out) ::  nFSkpt
!arrays
 integer,intent(out) :: fs2ibz(BSt%nkpt)
 real(dp),intent(in) :: fs_weight(BSt%mband,BSt%nkpt)

!Local variables-------------------------------
 integer :: iband,ikpt,nband_k,i1,i2
 logical :: seen
!arrays

! *************************************************************************
 
 nFSkpt=0
 bmin = HUGE(1)
 bmax = 0
 fs2ibz(:)=0

 do ikpt=1,BSt%nkpt
   nband_k=BSt%nband(ikpt+(isppol-1)*BSt%nkpt)

   seen=.FALSE.
   do iband=1,nband_k
    if (fs_weight(iband,ikpt) > tolweight) then 
      if (.not.seen) then
        seen=.TRUE.
        i1=iband
        i2=iband
      end if
    else 
      i2=iband
    end if
   end do

   if (seen) then 
     nFSkpt = nFSkpt+1 
     fs2ibz(nFSkpt) = ikpt
     bmin = MIN(bmin,i1)
     bmax = MAX(bmax,i2)
   end if
 end do !ikpt

end subroutine get_FS
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/bst_ismetal
!! NAME
!! bst_ismetal
!!
!! FUNCTION
!! Returns .TRUE. if metallic occupation scheme is used.
!!
!! INPUTS
!! bstruct<bstruct_type>=The bandstructure datatype
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function bst_ismetal(Bst)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bst_ismetal'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical :: bst_ismetal
 type(bandstructure_type),intent(in) :: Bst

! *************************************************************************

 bst_ismetal = ( ANY(Bst%occopt == (/3,4,5,6,7,8/)) )

end function bst_ismetal
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/bst_print_fs
!! NAME
!!  bst_print_fs
!!
!! FUNCTION
!!
!! INPUTS
!!  BSt<bandstructure_type>=The object describing the band structure.
!!  Cryst<crystal_structure>=Info on unit cell and symmetries.
!!  fname=File name for output.
!!  kptrlatt(3,3)=Matrix partly defining the mesh. See related input variable.
!!  nshiftk=Number of shifts in K-mesh (usually 1)
!!  shiftk(3,nshiftk)=The shifts of the mesh. Xcrysden requires Gamma-centered k-meshes.
!!  [fformat]=Integer flags specifying the format to be used:
!!     1 => BXSF Xcrysden file format [DEFAULT]
!!
!! OUTPUT
!!  ierr=Status error.
!!  BXSF file.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine bst_print_fs(BSt,Cryst,kptrlatt,nshiftk,shiftk,fname,ierr,fformat)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bst_print_fs'
 use interfaces_62_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nshiftk
 integer,intent(out) :: ierr
 integer,optional,intent(in) :: fformat
 character(len=*),intent(in) :: fname
 type(bandstructure_type),intent(in) :: BSt
 type(crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)

!Local variables-------------------------------
 integer :: fform
 logical :: use_timrev
 character(len=500) :: msg
!arrays

! *************************************************************************

 ABI_UNUSED(fname)
 !ABI_UNUSED(BSt%mband)

 ierr=0
 fform=1; if (PRESENT(fformat)) fform=fformat

 use_timrev = (Cryst%timrev==2)

 select case (fform)

 case (1)
   call printbxsf(BSt%eig,zero,BSt%fermie,Cryst%gprimd,kptrlatt,BSt%mband,BSt%nkpt,Bst%kptns,&
&    Cryst%nsym,Cryst%use_antiferro,Cryst%symrec,Cryst%symafm,use_timrev,Bst%nsppol,shiftk,nshiftk,fname,ierr)

 case default
   ierr = ierr+1
   write(msg,'(a,i0,a)')" Unsupported value for fform: ",fform," Fermi surface file won't be produced."
   MSG_WARNING(msg)
 end select

end subroutine bst_print_fs
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_dos
!! NAME
!!  get_dos
!!
!! FUNCTION
!!  Calculate the electronic density of states.
!!  Ideally the routine should report a DOS object 
!!  that can be written on a file using a method dump_dos 
!!
!! INPUTS
!!  BSt<Bandstructure_type>=The object describing the band structure.
!!  Kmesh<BZ_mesh_type>=Strcuture defining the BZ sampling.
!!  method
!!  fildos
!!  broad
!!  dosdeltae
!!
!! OUTPUT
!!  For the moment results are printed on an external file.
!!
!! PARENTS
!!      exc_interp_ham,m_gwannier
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_dos(BSt,Kmesh,method,fildos,broad,dosdeltae)

 use defs_basis
 use m_tetrahedron

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_dos'
 use interfaces_14_hidewrite
 use interfaces_62_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method
 real(dp),intent(in) :: dosdeltae,broad
 character(len=*),intent(in) :: fildos
 type(bandstructure_type),intent(in)  :: BSt
 type(BZ_mesh_type),intent(in)  :: Kmesh

!Local variables-------------------------------
!scalars
 integer,parameter :: option4dos=2
 integer :: unitdos,isppol,iband,ikpt
 integer :: nene,iene
 real(dp),parameter :: buffer=0.01_dp
 real(dp) :: nelect_out,entropy_out,max_occ,rcvol,enex
 real(dp) :: enemax,enemin,deltaene,integral_DOS
 character(len=500) :: msg
 logical :: ltest
!arrays
 integer :: G0(3)
 real(dp) :: minmax(2,BSt%nsppol)
 real(dp) :: gprimd(3,3),k1(3),k2(3)
 real(dp),allocatable :: eigen(:),doccde(:),occ(:)
 real(dp),allocatable :: dtweightde(:,:) 
 !real(dp),allocatable :: integ_dos(:,:),partial_dos(:,:)
 real(dp),allocatable :: tmp_eigen(:),total_dos(:)
 real(dp),allocatable :: total_integ_dos(:)
 real(dp),allocatable :: tweight(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 SELECT CASE (method)

 CASE (1) ! DOS with broadening.
   ! The unit is closed in getnel, not so elegant!
   unitdos=get_unit()
   open(unit=unitdos,file=fildos,status='unknown',form='formatted')

   ! To be consistent with the interface of getnel 
   max_occ=two/(BSt%nspinor*BSt%nsppol)  ! Will not work in the fixed moment case ????? MG should check why
   ABI_ALLOCATE(doccde,(BSt%mband*BSt%nkpt*BSt%nsppol))
   ABI_ALLOCATE(eigen ,(BSt%mband*BSt%nkpt*BSt%nsppol))
   ABI_ALLOCATE(occ   ,(BSt%mband*BSt%nkpt*BSt%nsppol))
   doccde(:)=zero
   call get_eneocc_vect(BSt,'eig',eigen)
   call get_eneocc_vect(BSt,'occ',occ  )

   deltaene=dosdeltae ; if (ABS(deltaene)<tol10) deltaene=0.001_dp

   ! Use Gaussian smearing. TODO should be input
   call getnel(doccde,deltaene,eigen,entropy_out,BSt%fermie,max_occ,BSt%mband,BSt%nband,&
&   nelect_out,BSt%nkpt,BSt%nsppol,occ,3,option4dos,BSt%tphysel,broad,unitdos,BSt%wtk)

   !if (nelect_out/=BSt%nelect) STOP
   ABI_DEALLOCATE(doccde)
   ABI_DEALLOCATE(eigen)
   ABI_DEALLOCATE(occ)

 CASE (2) ! Tetrahedron Method.
  !
  ! * Check input.
  ABI_CHECK(BSt%nkpt>=2,'At least 2 points are needed for tetrahedrons')

  if (Kmesh%nshift>1) then 
    write(msg,'(a,i0)')" For tetrahedrons, nshift must be (0,1) but found: ",Kmesh%nshift
    MSG_ERROR(msg)
  end if

  if ( ANY(BSt%nband/=BSt%nband(1)) ) then
    MSG_ERROR('For tetrahedrons, nband(:) must be constant')
  end if

  gprimd = Kmesh%gprimd
  rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
&  -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
&  +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

  ! Choose the lower and upper energies and Extend the range to a nicer value
  minmax=get_minmax(BSt,'eig')
  enemin=MINVAL(minmax(1,:))-buffer
  enemax=MAXVAL(minmax(2,:))+buffer
  enemax=0.1_dp*CEILING(enemax*ten)
  enemin=0.1_dp*FLOOR(enemin*ten)

  !Choose the energy increment
  deltaene=dosdeltae
  if (ABS(deltaene)<tol10) deltaene=0.0005_dp ! Higher resolution possible (and wanted) for tetrahedron
  nene=NINT((enemax-enemin)/deltaene)+1

  ! === kpoints must be the same and with the same ordering ===
  ltest=(BSt%nkpt==Kmesh%nibz)
  ABI_CHECK(ltest,'Mismatch in number of k-points')
  ltest=.TRUE.
  do ikpt=1,BSt%nkpt
    k1=BSt%kptns(:,ikpt)
    k2=Kmesh%ibz(:,ikpt)
    ltest=(ltest.and.isamek(k1,k2,G0))
  end do
  ABI_CHECK(ltest,'k-points in Kmesh and BSt are not equivalent.')

  ltest=(associated(Kmesh%tetrahedra%tetra_mult).and.associated(Kmesh%tetrahedra%tetra_full))
  !TODO recalculate tetra on-the-fly 
  ABI_CHECK(ltest,'tetra_full or tetra_mult not associated.')

  ! Here ndosfraction has been removed, should provide a low 
  ! level routine to perform integrations with tetrahedrons of matrix elements provided by the user

  max_occ=two/(BSt%nspinor*BSt%nsppol)  ! Will not work in the fixed moment case ????? MG should check why
  ABI_ALLOCATE(tweight,(BSt%nkpt,nene))
  ABI_ALLOCATE(dtweightde,(BSt%nkpt,nene))
  ABI_ALLOCATE(total_dos,(nene))
  ABI_ALLOCATE(total_integ_dos,(nene))

  unitdos=get_unit()
  open(unit=unitdos,file=fildos,status='unknown',form='formatted')

  ! For each spin and band, interpolate over kpoints, 
  ! calculate integration weights and DOS contribution.
  do isppol=1,BSt%nsppol

    total_dos(:)=zero
    total_integ_dos(:)= zero

    if (BSt%nsppol==2) then
      if (isppol==1) write(msg,'(a,16x,a)')  '#','Spin-up DOS'
      if (isppol==2) write(msg,'(2a,16x,a)')  ch10,'#','Spin-dn DOS'
      call wrtout(unitdos,msg,'COLL')
    end if

    ABI_ALLOCATE(tmp_eigen,(BSt%nkpt))

    do iband=1,BSt%nband(1)
      ! For each band get its contribution
      tmp_eigen(:)=Bst%eig(iband,:,isppol)
      !  
      ! === Calculate integration weights at each irred kpoint ===
      ! * Blochl et al PRB 49 16223 ===
      call get_tetra_weight(tmp_eigen,enemin,enemax,&
&             max_occ,nene,Bst%nkpt,Kmesh%tetrahedra,&
&            tweight,dtweightde)

      ! === Calculate DOS and integrated DOS projected with the input dos_fractions ===
!      call get_dos_1band (dos_fractions(:,iband,isppol,:),enemin,enemax,&
!&      integ_dos(:,:,iband),nene,Bst%nkpt,ndosfraction,partial_dos(:,:,iband),tweight,dtweightde)

      do iene=1,nene
        do ikpt=1,Bst%nkpt
          total_dos(iene)       = total_dos(iene) + dtweightde(ikpt,iene) 
          total_integ_dos(iene) = total_integ_dos(iene) + tweight(ikpt,iene) 
        end do
      end do
    end do ! iband

    ABI_DEALLOCATE(tmp_eigen)

    ! === Write DOS values ===
    call wrtout(unitdos,'#  energy(Ha)     DOS  integrated DOS','COLL')

    enex=enemin
    do iene=1,nene
      ! Print the data for this energy. Note the upper limit, to be
      ! consistent with the format. The use of "E" format is not adequate,
      ! for portability of the self-testing procedure.
      ! write(msg, '(i5,f9.4,f14.6)' ) iene-1,enex,total_dos(iene,:)
      write(msg,'(f11.5,2f10.4)')enex,min(total_dos(iene),9999.9999_dp),total_integ_dos(iene)
      call wrtout(unitdos,msg,'COLL')
      enex=enex+deltaene
    end do

    !integral_DOS=deltaene*SUM(total_dos(iene,:))
    integral_DOS=total_integ_dos(nene)
    write(msg,'(a,es16.8)')' tetrahedron : integrate to',integral_DOS
    call wrtout(std_out,msg,'COLL')
  end do !isppol

  close(unitdos)

  ! === Free memory ===
  ABI_DEALLOCATE(total_dos)
  ABI_DEALLOCATE(total_integ_dos)
  ABI_DEALLOCATE(tweight)
  ABI_DEALLOCATE(dtweightde)

 CASE DEFAULT
   write(msg,'(a,i0)')' Wrong value for method= ',method
   MSG_BUG(msg)
 END SELECT

 DBG_EXIT("COLL")

end subroutine get_dos
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/update_occ
!! NAME
!! update_occ
!!
!! FUNCTION
!! Calculate new occupation numbers, the Fermi level and the Max occupied band index 
!! for each spin channel starting from the the knowledge of eigenvalues.
!!
!! INPUTS
!!  fixmom=if differ from -99.99d0, fix the magnetic moment (in Bohr magneton)
!!  [stmbias]=
!!  [prtvol]=Verbosity level (0 for lowest level)
!!  BSt<bandstructure_type>=Info on the band structure, the smearing technique and the physical temperature used.
!!
!! OUTPUT
!!  see also SIDE EFFECTS.
!!
!! SIDE EFFECTS
!!  === For metallic occupation the following quantites are recalculated ===
!!   %fermie=the new Fermi energy
!!   %entropy=the new entropy associated with the smearing.
!!   %occ(mband,nkpt,nsppol)=occupation numbers
!!   %doccde(mband,nkpt,nsppol)=derivative of occupancies wrt the energy for each band and k point
!!  === In case of semiconductors ===
!!   All the quantitities in BSt are left unchanged with the exception of:
!!   %fermie=Redefined so that it is in the middle of the gap
!!   %entropy=Set to zero
!!
!! PARENTS
!!      bethe_salpeter,elphon,m_ebands,screening,setup_bse,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine update_occ(BSt,fixmom,stmbias,prtvol)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'update_occ'
 use interfaces_14_hidewrite
 use interfaces_62_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(bandstructure_type),intent(inout) :: BSt
 integer,optional,intent(in) :: prtvol
 real(dp),intent(in) :: fixmom
 real(dp),optional,intent(in) :: stmbias
!arrays

!Local variables-------------------------------
!scalars
 integer :: iband,mband,ikibz,nkpt,isppol,nsppol,my_prtvol,nband_k
 real(dp) :: entropy,fermie,stmbias_local,ndiff,cbot,vtop,maxocc
 character(len=500) :: msg
!arrays
 real(dp) :: nelect_spin(BSt%nsppol),condbottom(BSt%nsppol),valencetop(BSt%nsppol)
 real(dp),allocatable :: doccdet(:),occt(:),eigent(:)

! *************************************************************************

 my_prtvol    =0   ; if (PRESENT(prtvol )) my_prtvol    =prtvol
 stmbias_local=zero; if (PRESENT(stmbias)) stmbias_local=stmbias

 if (BSt%occopt>=3.and.BSt%occopt<=8) then !  If occupation is metallic have to compute new occupation numbers.
   write(msg,'(a,f9.5)')' metallic scheme, calling newocc with fixmom = ',fixmom
   call wrtout(std_out,msg,'COLL')

   mband  = BSt%mband ! to make the interface of newocc happy.
   nkpt   = BSt%nkpt
   nsppol = BSt%nsppol

   ABI_ALLOCATE(eigent,(mband*nkpt*nsppol))
   call get_eneocc_vect(BSt,'eig',eigent)

   ABI_ALLOCATE(occt,(mband*nkpt*nsppol))
   ABI_ALLOCATE(doccdet,(mband*nkpt*nsppol))

   call newocc(doccdet,eigent,entropy,fermie,fixmom,mband,BSt%nband,&
&    BSt%nelect,BSt%nkpt,BSt%nspinor,BSt%nsppol,occt,BSt%occopt,&
&    my_prtvol,stmbias_local,BSt%tphysel,BSt%tsmear,BSt%wtk)
   !
   ! Save output in BSt%. 
   BSt%entropy = entropy
   BSt%fermie  = fermie
   call put_eneocc_vect(BSt,'occ'   ,occt   ) 
   call put_eneocc_vect(BSt,'doccde',doccdet) 
   ABI_DEALLOCATE(eigent)
   ABI_DEALLOCATE(occt)
   ABI_DEALLOCATE(doccdet)

 else  !  Semiconductor or Insulator.
   ! 
   ! FIXME here there is an inconsistency btw GW and Abinit 
   ! In abinit Fermi is set to HOMO while in GW fermi is in the middle
   ! of Gap. In case of crystal systems, the later convention should be preferable.
   ! Anyway we have to decide and follow a unique convention to avoid problems.
   !
   ! occupation factors MUST be initialized
   if (ALL(ABS(Bst%occ) < tol6)) then
     msg = "occupation factors are not initialized, likely due to the use of iscf=-2"
     MSG_ERROR(msg)
   end if

   maxocc=two/(BSt%nsppol*BSt%nspinor)

   ! * Calculate the valence index for each spin channel.
   do isppol=1,BSt%nsppol
     valencetop(isppol)= smallest_real
     condbottom(isppol)= greatest_real

     do ikibz=1,BSt%nkpt
       nband_k=BSt%nband(ikibz+(isppol-1)*BSt%nkpt) 
       do iband=1,nband_k
         if (BSt%occ(iband,ikibz,isppol)/maxocc>one-tol6 .and. valencetop(isppol)<BSt%eig(iband,ikibz,isppol)) then 
           valencetop(isppol)=BSt%eig(iband,ikibz,isppol)
         end if
         if (BSt%occ(iband,ikibz,isppol)/maxocc<tol6 .and. condbottom(isppol)>BSt%eig(iband,ikibz,isppol)) then 
           condbottom(isppol)=BSt%eig(iband,ikibz,isppol)
         end if
       end do
     end do 

   end do 

   vtop=MAXVAL(valencetop)
   cbot=MINVAL(condbottom)
   write(msg,'(a,f6.2,2a,f6.2)')&
&    ' top of valence       [eV] ',vtop*Ha_eV,ch10,&
&    ' bottom of conduction [eV] ',cbot*Ha_eV
   call wrtout(std_out,msg,'COLL')
   if (BSt%nsppol==2) then 
     if (ABS(vtop-MINVAL(valencetop))>tol6) then 
       write(msg,'(a,i2)')' top of valence is spin ',MAXLOC(valencetop)
       call wrtout(std_out,msg,'COLL')
     end if
     if (ABS(cbot-MAXVAL(condbottom))>tol6) then 
       write(msg,'(a,i2)')' bottom of conduction is spin ',MINLOC(condbottom)
       call wrtout(std_out,msg,'COLL')
     end if
   end if

   ! === Save output === 
   ! Here I dont know if it is better to be consistent with the abinit convention i.e fermi=vtop
   BSt%entropy=zero
   BSt%fermie=(vtop+cbot)/2 
   if (ABS(cbot-vtop)<1.d-4) BSt%fermie=vtop ! To avoid error on the last digit FIXME is it really needed
 end if

 write(msg,'(a,f6.2,a)')' Fermi energy         [eV] ',BSt%fermie*Ha_eV,ch10
 call wrtout(std_out,msg,'COLL')
 !
 ! === Compute number of electrons for each spin channel ===
 nelect_spin(:)=zero 
 do isppol=1,BSt%nsppol
   do ikibz=1,BSt%nkpt
     nband_k=BSt%nband(ikibz+(isppol-1)*BSt%nkpt)
     nelect_spin(isppol)= nelect_spin(isppol) + BSt%wtk(ikibz)*SUM(BSt%occ(1:nband_k,ikibz,isppol))
   end do
 end do

 ndiff=BSt%nelect-SUM(nelect_spin)
 if (my_prtvol>0) then
   write(msg,'(2a,f6.2,2a,f7.4)')ch10,&
&    ' total number of electrons = ',SUM(nelect_spin),ch10,&
&    ' input and calculated no. of electrons differ by ',ndiff 
   call wrtout(std_out,msg,'COLL')
 end if

 if (ABS(ndiff)>5.d-2*BSt%nelect) then
   write(msg,'(2a,2(a,f6.2))')&
&    ' Too large difference in no. of electrons:,',ch10,&
&    ' Expected= ',BSt%nelect,' Calculated= ',SUM(nelect_spin)
   MSG_ERROR(msg)
 end if

end subroutine update_occ
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ReportGap
!! NAME
!! ReportGap
!!
!! FUNCTION
!!  Print info on the fundamental and optical gap.
!!
!! INPUTS
!!  BSt<bandstructure_type>=Info on the band structure, the smearing technique and the physical temperature used.
!!  [header]=Optional title.
!!  [kmask]=Logical mask used to exclude k-points.
!!  [unit]=Optional unit for output (std_out if not specified)
!!  [mode_paral]=Either "COLL" or "PERS", former is default.
!!
!! OUTPUT
!!  writing.
!!  [gaps(2,nsppol)]=Fundamental and optical gaps.
!!
!! PARENTS
!!      exc_diago_driver,gstate,setup_bse,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ReportGap(BSt,header,kmask,unit,mode_paral,gaps)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ReportGap'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: unit
 character(len=4),intent(in),optional :: mode_paral
 character(len=*),intent(in),optional :: header
 type(bandstructure_type),intent(in)  :: BSt
!arrays
 real(dp),optional,intent(out) :: gaps(2,BSt%nsppol)
 logical,optional,intent(in) ::  kmask(BSt%nkpt)

!Local variables-------------------------------
!scalars
 integer :: ikibz,nband_k,isppol,ikopt,ivk,ick,ivb,icb,my_unt,first
 real(dp),parameter :: tol_fermi=tol6
 real(dp) :: gap,optical_gap
 logical :: ismetal
 character(len=4) :: my_mode
 character(len=500) :: msg
!arrays
 integer :: val_idx(BSt%nkpt,BSt%nsppol)
 real(dp) :: top_valence(BSt%nkpt),bot_conduct(BSt%nkpt) 
 logical :: my_kmask(BSt%nkpt)
! *********************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral
 my_kmask=.TRUE.; if (PRESENT(kmask     )) my_kmask=kmask

 if (PRESENT(gaps)) gaps=zero

 val_idx(:,:) = get_valence_idx(BSt,tol_fermi)
 first=0

 do isppol=1,BSt%nsppol

  ismetal=ANY(val_idx(:,isppol)/=val_idx(1,isppol)) 
  if (ismetal) CYCLE ! No output if system i metallic
  first=first+1
  if (first==1) then
    msg=ch10
    if (PRESENT(header)) msg=ch10//' === '//TRIM(ADJUSTL(header))//' === '
    call wrtout(my_unt,msg,my_mode) 
  end if

  ivb=val_idx(1,isppol)
  icb=ivb+1

  do ikibz=1,BSt%nkpt
    if (.not.my_kmask(ikibz)) CYCLE
    nband_k=BSt%nband(ikibz+(isppol-1)*BSt%nkpt)
    top_valence(ikibz)=BSt%eig(ivb,ikibz,isppol)
    if (icb>nband_k) GOTO 10 ! Only occupied states are present, no output!
    bot_conduct(ikibz)=BSt%eig(icb,ikibz,isppol)
  end do

  ! === Get minimum of the optical Gap ===
  ikopt= imin_loc(bot_conduct-top_valence,MASK=my_kmask)
  optical_gap=bot_conduct(ikopt)-top_valence(ikopt)

  ! === Get fundamental Gap ===
  ick = imin_loc(bot_conduct,MASK=my_kmask)
  ivk = imax_loc(top_valence,MASK=my_kmask)
  gap = BSt%eig(icb,ick,isppol)-BSt%eig(ivb,ivk,isppol)

  write(msg,'(a,i2,a,2(a,f8.4,a,3f8.4,a),33x,a,3f8.4)')&
&  '  >>>> For spin ',isppol,ch10,&
&  '   Minimum optical gap = ',optical_gap*Ha_eV,' [eV], located at k-point      : ',BSt%kptns(:,ikopt),ch10,&
&  '   Fundamental gap     = ',gap*Ha_eV,        ' [eV], Top of valence bands at : ',BSt%kptns(:,ivk),ch10,  &
&                                                '       Bottom of conduction at : ',BSt%kptns(:,ick)
  call wrtout(my_unt,msg,my_mode) 

  if (PRESENT(gaps)) then
    gaps(:,isppol) = (/gap,optical_gap/)
  end if

 end do !isppol

 return

10 msg = " ReportGap : WARNING - not enough states to calculate the band gap. "
 call wrtout(my_unt,msg,my_mode) 

end subroutine ReportGap
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/bst_write_bands
!! NAME
!! bst_write_bands
!!
!! FUNCTION
!! Plots the interpolated band structure in Xmgrace format  
!! Based on plot_interpolate_xmgrace, routine of the Wannier90 code.
!!
!! INPUTS
!!  BSt<bandstructure_type>The type containing the data.
!!  gmet(3,3)=Metric in G-space.
!!  fname=The name of the file.
!!  [bands_label(BSt%nkpt)]=Optional string identifying the k-points on the path.
!!  [format]=integer flag specifying the format to be used.
!!    1 for xmgrace format [default]
!!
!! OUTPUT
!!  ierr=Status error.
!!  Results are written on file.
!!
!! PARENTS
!!      bloch_interp,m_gwannier
!!
!! CHILDREN
!!
!! SOURCE

subroutine bst_plot_bands(BSt,gmet,fname,ierr,bands_label,format)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bst_plot_bands'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: format
 integer,intent(out) :: ierr
 character(len=*),intent(in) :: fname
 type(bandstructure_type),intent(in) :: BSt
!arrays
 real(dp),intent(in) :: gmet(3,3)
 character(len=1),optional,intent(in) :: bands_label(BSt%nkpt)

!Local variables-------------------------------
 integer :: fform,xmgr_ierr
 character(len=500) :: msg
!arrays 
 character(len=1) :: my_bands_label(BSt%nkpt)
! *************************************************************************

 ierr=0

 fform=1; if (PRESENT(format)) fform=format
 my_bands_label=""; if (PRESENT(bands_label)) my_bands_label = bands_label

 if (ANY(BSt%nband /= BSt%nband(1))) then
   MSG_WARNING("Variable nband not supported")
   ierr=1; RETURN
 end if

 select case (fform) 

 case (1)
  call xmgrace_bandplot(BSt%mband,BSt%nkpt,BSt%nsppol,BSt%kptns,BSt%eig,gmet,fname,my_bands_label,xmgr_ierr)
  ierr = ierr + xmgr_ierr

 case default
   write(msg,'(a,i0)')" Unsupported value for fform: ",fform
   MSG_ERROR(msg)
 end select

 if (ierr/=0) then
   MSG_WARNING("Band plot file has not been generated")
 end if

end subroutine bst_plot_bands 
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/SelectBands
!! NAME
!! SelectBands
!!
!! FUNCTION
!!  Return a new object of type bandstructure_type corresponding to 
!!  to the logical masks specified in input.
!!
!! INPUTS
!!  BSt<bandstructure_type>=Info on the band structure 
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function SelectBands(BSt,bandselect,kselect,spinselect) result(NewBSt)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'SelectBands'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(bandstructure_type),intent(in)  :: BSt 
 logical,optional,intent(in) :: bandselect(BSt%mband,BSt%nkpt,BSt%nsppol)
 logical,optional,intent(in) :: kselect(BSt%nkpt)
 logical,optional,intent(in) :: spinselect(BSt%nsppol)
 type(bandstructure_type) :: NewBSt

!Local variables-------------------------------
!scalars
 integer :: ikpt,nband_k,iband,isppol,new_nband_k,idx,isp_new,ik_new,ib_new
!arrays
 logical :: band_in(BSt%mband,BSt%nkpt,BSt%nsppol)
 logical :: kpt_in(BSt%nkpt)
 logical :: spin_in(BSt%nsppol)
! *********************************************************************

 band_in=.TRUE.; if (PRESENT(bandselect)) band_in=bandselect
 kpt_in =.TRUE.; if (PRESENT(kselect   )) kpt_in =kselect
 spin_in=.TRUE.; if (PRESENT(spinselect)) spin_in=spinselect

 NewBSt%nkpt   =COUNT(kpt_in)
 NewBSt%nspinor=BSt%nspinor      
 NewBSt%nsppol =COUNT(spin_in)  
 NewBSt%occopt =BSt%occopt        

 NewBSt%entropy = BSt%entropy 
 NewBSt%fermie  = BSt%fermie  
 NewBSt%nelect  = BSt%nelect  
 NewBSt%tphysel = BSt%tphysel 
 NewBSt%tsmear  = BSt%tsmear  

 ABI_ALLOCATE(NewBSt%nband,(NewBSt%nkpt*NewBSt%nsppol))
 idx=0
 do isppol=1,BSt%nsppol
   if (.not.spin_in(isppol)) CYCLE
   do ikpt=1,BSt%nkpt
     if (.not.kpt_in(ikpt)) CYCLE
     nband_k=BSt%nband(ikpt+(isppol-1)*BSt%nkpt)
     idx=idx+1
     new_nband_k=COUNT(band_in(1:nband_k,ikpt,isppol))
     NewBSt%nband(idx)=new_nband_k
   end do
 end do

 NewBSt%mband =MAXVAL(NewBSt%nband)
 NewBSt%bantot=SUM   (NewBSt%nband)

 ABI_ALLOCATE(NewBSt%istwfk,(NewBSt%nkpt))
 ABI_ALLOCATE(NewBSt%npwarr,(NewBSt%nkpt))
 ABI_ALLOCATE(NewBSt%kptns,(3,NewBSt%nkpt))
 ABI_ALLOCATE(NewBSt%eig   ,(NewBSt%mband,NewBSt%nkpt,NewBSt%nsppol))
 ABI_ALLOCATE(NewBSt%occ   ,(NewBSt%mband,NewBSt%nkpt,NewBSt%nsppol))
 ABI_ALLOCATE(NewBSt%doccde,(NewBSt%mband,NewBSt%nkpt,NewBSt%nsppol))
 ABI_ALLOCATE(NewBSt%wtk,(NewBSt%nkpt))

 ! === Copy arrays that depend only on k ===
 idx=0
 do ikpt=1,BSt%nkpt
   if (.not.kpt_in(ikpt)) CYCLE
   idx=idx+1
   NewBSt%istwfk(idx) =BSt%istwfk(ikpt)      
   NewBSt%npwarr(idx) =BSt%npwarr(ikpt)      
   NewBSt%kptns(:,idx)=BSt%kptns(:,ikpt)    
   NewBSt%wtk(idx)    =BSt%wtk(ikpt)        
 end do

 ! * Renormalize weights
 NewBSt%wtk=NewBSt%wtk/SUM(NewBSt%wtk)

 ! === Copy arrays that depend on b-k-s ===
 isp_new=0
 do isppol=1,BSt%nsppol
   if (.not.spin_in(isppol)) CYCLE
   isp_new=isp_new+1

   ik_new=0
   do ikpt=1,BSt%nkpt
     if (.not.kpt_in(ikpt)) CYCLE
     ik_new=ik_new+1

     nband_k=BSt%nband(ikpt+(isppol-1)*BSt%nkpt)
     ib_new=0
     do iband=1,nband_k
       if (.not.band_in(iband,ikpt,isppol)) CYCLE
       ib_new=ib_new+1
       NewBSt%eig   (ib_new,ik_new,isp_new)=BSt%eig   (iband,ikpt,isppol)        
       NewBSt%occ   (ib_new,ik_new,isp_new)=BSt%occ   (iband,ikpt,isppol)        
       NewBSt%doccde(ib_new,ik_new,isp_new)=BSt%doccde(iband,ikpt,isppol)     
     end do !iband
   end do !ikpt
 end do !isppol

end function SelectBands
!!***

!!****f* m_ebands/ExpandBands
!! NAME
!! ExpandBands
!!
!! FUNCTION
!!  Return a new object of type bandstructure_type corresponding to a list of k-points 
!!  specified in input. Symmetry properties of the eigenvectors are used to 
!!  symmetrize energies and occupation numbers.
!!
!! INPUTS
!!  BSt<Bandstructure_type>=Info on the band structure 
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  In info/=0, some of the k-points in klist have no corresponding image in Bst_in%kptns
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ExpandBands(BSt_in,nklist,klist,use_tr,use_afm,nsym,symrec,symafm,info) result(BSt_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ExpandBands'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym,nklist
 integer,intent(out) :: info
 logical,intent(in) :: use_tr,use_afm
 type(bandstructure_type),intent(in) :: BSt_in 
 type(bandstructure_type) :: BSt_out
!arrays
 real(dp),intent(in) :: klist(3,nklist)
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(in) :: symafm(nsym)

!Local variables-------------------------------
!scalars
 integer :: iklist,itim,timrev,isym,ikibz,nkfound,ii,jj
 integer :: nband_k,isppol,ikdx
 character(len=500) :: msg
 logical :: found
!arrays
 integer :: G0(3)
 real(dp) :: krot(3)
 integer,allocatable :: klist2ibz(:)
! *********************************************************************

 do ii=1,nklist-1
   do jj=ii+1,nklist
     if (isamek(klist(:,ii),klist(:,jj),G0)) then
       write(msg,'(a,i4,a,i4,5a)')&
&       ' Points ',ii,' and ',jj,' in klist are equivalent ',ch10,&
&       ' This is not allowed in the present implementation.',ch10,&
&       ' Change the input values to avoid duplicated k-points. '
       MSG_ERROR(msg)
     end if
   end do
 end do

 info=0
 timrev=1 ; if (use_tr) timrev=2
 ABI_ALLOCATE(klist2ibz,(nklist))
 klist2ibz=0

 do iklist=1,nklist
   found=.FALSE.

ibzl: do ikibz=1,BSt_in%nkpt
        do itim=1,timrev
          do isym=1,nsym
            if (use_afm.and.symafm(isym)==-1) CYCLE
            ! * Form IS k ===
            krot(:)=(3-2*itim)*MATMUL(symrec(:,:,isym),BSt_in%kptns(:,ikibz))

            ! * Check whether it is equal to klist(:,ilist) within a RL vector.
            ! FIXME see notes below related to this function
            if (isamek(krot,klist(:,iklist),G0)) then 
              found=.TRUE.
              klist2ibz(iklist)=ikibz
              EXIT ibzl
            end if

          end do !isym 
        end do !itim
      end do ibzl

      if (.not.found) info=info+1
 end do !iklist

 nkfound=COUNT(klist2ibz/=0)

 Bst_out%nkpt    = nkfound
 Bst_out%nspinor = BSt_in%nspinor      
 Bst_out%nsppol  = Bst_in%nsppol
 Bst_out%occopt  = Bst_in%occopt        

 Bst_out%entropy = Bst_in%entropy 
 Bst_out%fermie  = Bst_in%fermie  
 Bst_out%nelect  = Bst_in%nelect  
 Bst_out%tphysel = Bst_in%tphysel 
 Bst_out%tsmear  = Bst_in%tsmear  

 ABI_ALLOCATE(Bst_out%nband,(Bst_out%nkpt*Bst_out%nsppol))
 ikdx=0
 do isppol=1,Bst_out%nsppol
   do iklist=1,nklist
     ikibz=klist2ibz(iklist)
     if (ikibz==0) CYCLE
     ikdx=ikdx+1
     nband_k=Bst_in%nband(ikibz+(isppol-1)*Bst_in%nkpt)
     Bst_out%nband(ikdx)=nband_k
   end do
 end do

 Bst_out%mband =MAXVAL(Bst_out%nband)
 Bst_out%bantot=SUM   (Bst_out%nband)

 ABI_ALLOCATE(Bst_out%istwfk,(Bst_out%nkpt))
 ABI_ALLOCATE(Bst_out%nband,(Bst_out%nkpt*Bst_out%nsppol))
 ABI_ALLOCATE(Bst_out%npwarr,(Bst_out%nkpt))
 ABI_ALLOCATE(Bst_out%kptns,(3,Bst_out%nkpt))
 ABI_ALLOCATE(Bst_out%eig   ,(Bst_out%mband,Bst_out%nkpt,Bst_out%nsppol))
 ABI_ALLOCATE(Bst_out%occ   ,(Bst_out%mband,Bst_out%nkpt,Bst_out%nsppol))
 ABI_ALLOCATE(Bst_out%doccde,(Bst_out%mband,Bst_out%nkpt,Bst_out%nsppol))
 ABI_ALLOCATE(Bst_out%wtk,(Bst_out%nkpt))
 !
 ! === Copy arrays that depend only on k ===
 ikdx=0
 do iklist=1,nklist
   ikibz=klist2ibz(iklist)
   if (ikibz==0) CYCLE
   ikdx=ikdx+1
   Bst_out%istwfk(ikdx) =Bst_in%istwfk(ikibz)      
   Bst_out%npwarr(ikdx) =Bst_in%npwarr(ikibz)      
   Bst_out%kptns(:,ikdx)=klist(:,iklist) !Use klist, not Bst_in%kptns
   Bst_out%wtk(ikdx)    =Bst_in%wtk(ikibz)        
 end do
 !
 ! * Renormalize weights
 Bst_out%wtk=Bst_out%wtk/SUM(Bst_out%wtk)

 ! === Copy arrays that depend on b-k-s ===
 do isppol=1,Bst_in%nsppol
   ikdx=0
   do iklist=1,nklist
     ikibz=klist2ibz(iklist)
     if (ikibz==0) CYCLE
     ikdx=ikdx+1
     Bst_out%eig   (:,ikdx,isppol)=Bst_in%eig   (:,ikibz,isppol)        
     Bst_out%occ   (:,ikdx,isppol)=Bst_in%occ   (:,ikibz,isppol)        
     Bst_out%doccde(:,ikdx,isppol)=Bst_in%doccde(:,ikibz,isppol)     
   end do !ikpt
 end do !isppol

 ABI_DEALLOCATE(klist2ibz)

end function ExpandBands
!!***

!!****f* m_ebands/joint_dos
!! NAME
!! joint_dos
!!
!! FUNCTION
!!  Calculate the joint density of states 
!!   $ J(\omega,q)=\sum_{k,b1,b2} \delta(\omega-\epsilon_{k-q,b2}+\epsilon_{k,b1}) $
!!  for a given set of external q-points. 
!!
!!  Two different quadrature methods are employed according to the input variable method:
!!   method=1 => simple gaussian broadenig 
!!   method=2 => tetrahedron method
!!
!! INPUTS
!!  fname=Name of the output file.
!!  Cryst<Crystal_structure>=Info on unit cell and its symmetries
!!    %nsym=number of symmetry operations
!!    %symrec(3,3,nsym)=symmetry operations in reciprocal space (reduced coordinates)
!!    %tnons(3,nsym)=fractional translations
!!    %timrev=2 if time-reversal holds, 1 otherwise
!!  qibz(3,:)=q-points in reduced coordinates. If size is 0, calculate all the q-points
!!  BSt<Bandstructure_type>=Type gathering info on the electronic band structure
!!   %kpt(3,nkpt)=the irreducible k-points (reduced coordinates)
!!   %nkpt=Number of irreducible k-points
!!   %nsppol=Number of Independendent spin polarizations
!!   %eig(mband,nkpt,nsppol)=energies
!!   %occ(mband,nkpt,nsppol)=occupation numbers
!!   %nband(nkpt*nsppol)=number of bands for each k-point and spin
!!  method=1 for gaussian broadening, 2 for tetrahedron 
!!  step=freqency step for JDOS
!!  broad=only for gaussian method, the broadening in Ha 
!!  
!! OUTPUT
!!  Only write
!!
!! TODO
!!  1) This should be a method of Bandstructure_type, waiting for restructuring
!!  of the build systems to solve dependencies.
!!  Ideally this method should report a class (J)DOS, 
!!
!!  2) Symmetries are not yet used to reduce the number of k-points (it should be easy)
!!
!!  3) tetrahedron method not yet implemented.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine joint_dos(qibz,Cryst,BSt,fname,method,step,broad)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'joint_dos'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method
 real(dp),intent(in) :: broad,step
 character(len=fnlen),intent(in) :: fname
 type(Crystal_structure),intent(in) :: Cryst
 type(bandstructure_type),intent(in) :: BSt
!arrays
 real(dp),target,intent(in) :: qibz(:,:)

!Local variables-------------------------------
!scalars
 integer :: ib1,ib2,iend,ikbz,ikibz,ikmq_bz,ikmq_ibz,iomega,nqibz
 integer :: iqibz,isppol,istart,isym_k,isym_kmq,itim_k,itim_kmq,nband_k,nband_kmq
 integer :: nfound,nomega,prtvol,unt,kptopt
 !integer :: npwe,npwvec,use_umklp
 real(dp),parameter :: TOL_OCC=tol6
 real(dp) :: dossmear,ene1,ene2,gaussfactor,gaussprefactor,gaussval,max_transition,occ1
 real(dp) :: occ2,trans,xx
 logical :: ltest
 character(len=500) :: msg
 type(BZ_mesh_type) :: Kmesh,Qmesh
 type(little_group) :: Ltg_q
!arrays
 integer :: G0(3)
 real(dp) :: kbz(3),kmq(3),qq(3)
 real(dp) :: eminmax(2,BSt%nsppol)
 real(dp),allocatable :: jdos(:,:,:),omega(:) 
 real(dp),pointer :: qcalc(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 MSG_ERROR("Work in progress")
 !
 ! === Check input ===
 ltest=(Cryst%timrev==1.or.Cryst%timrev==2)
 ABI_CHECK(ltest,'timrev should be 1 or 2,')
 ltest=(method==1.or.method==2)
 ABI_CHECK(ltest,'method should be 1 or 2,')
 !
 ! === Initialize the k-mesh structure ===
 prtvol=0; kptopt=1
 call init_kmesh(Kmesh,Cryst,BSt%nkpt,BSt%kptns,kptopt)

 nqibz=SIZE(qibz,DIM=2)
 if (nqibz/=0) then 
   qcalc => qibz(:,:)
 else
   call find_qmesh(Qmesh,Cryst,Kmesh)
   nqibz=Qmesh%nibz
   qcalc => Qmesh%ibz(:,:)
 end if 

 ! === Linear mesh starting at zero ===
 eminmax=get_minmax(BSt,'eig')
 max_transition=MAXVAL(eminmax(2,:)-eminmax(1,:)) 
 if (method==1) max_transition=max_transition+five*broad+tol6

 nomega=NINT(max_transition/step)+1
 ABI_ALLOCATE(omega,(nomega))
 omega(:)=arth(zero,step,nomega)

 ABI_ALLOCATE(jdos,(nomega,BSt%nsppol,nqibz))
 jdos=zero

 SELECT CASE (method) 

 CASE (1) ! Gaussian method.
   gaussprefactor=one/(dossmear*SQRT(two_pi))
   gaussfactor=one/(sqrt2*dossmear)
   write(msg,'(4a,f8.5,2a,f8.5)')ch10,&
&    ' joint_dos : calculating Joint DOS using gaussian method :',ch10,&
&    ' gaussian smearing [eV] = ',broad*Ha_eV,ch10,&
&    ' energy step       [eV] = ',step*Ha_eV

 CASE (2) ! Tetrahedron method.
   write(msg,'(2a)')ch10,&
&   ' mkphdos : calculating joint DOS using tetrahedron method '
   MSG_ERROR("not implemented yet")

 CASE DEFAULT 
   write(msg,'(a,i3)')'Wrong value for method= ',method
   MSG_BUG(msg)
 END SELECT
 call wrtout(std_out,msg,'COLL')

 call nullify_little_group(Ltg_q)

 do iqibz=1,nqibz
   qq(:)=qcalc(:,iqibz) 

   !TODO add little group
   !use_umklp=1 ; npwvec=0 ; npwe=0
   !call setup_little_group(qq,Kmesh,Cryst,use_umklp,prtvol,Ltg_q,0)

   do ikbz=1,Kmesh%nbz
     do isppol=1,BSt%nsppol
       !
       ! === Get k and k-q in BZ and their symmetric in the IBZ ===
       call get_BZ_item(Kmesh,ikbz,kbz,ikibz,isym_k,itim_k)
       call get_BZ_diff(Kmesh,kbz,qq,ikmq_bz,G0,nfound)
       if (nfound/=1) then 
         msg='Multiple k1-k2 points found in the BZ zone'
         MSG_ERROR(msg)
       end if

       call get_BZ_item(Kmesh,ikmq_bz,kmq,ikmq_ibz,isym_kmq,itim_kmq)

       nband_k  =BSt%nband(ikibz   +(isppol-1)*BSt%nkpt)
       nband_kmq=BSt%nband(ikmq_ibz+(isppol-1)*BSt%nkpt)

       do ib2=1,nband_kmq
         occ2 = BSt%occ(ib2,ikmq_ibz,isppol)
         ene2 = BSt%eig(ib2,ikmq_ibz,isppol)

         do ib1=1,nband_k
           occ1 = BSt%occ(ib1,ikibz,isppol); if (ABS(occ1*(one-occ2))<TOL_OCC) CYCLE !occ1 might be not monotonic
           ene1 = BSt%eig(ib1,ikibz,isppol)

           ! three*broad should be enough, five should give smooth curves for coarse k-meshes.
           trans=ene2-ene1
           iend  =NINT((trans+five*broad)/step) ; if (iend>nomega) iend  =nomega
           istart=NINT((trans-five*broad)/step) ; if (istart<=0)   istart=1
           !
           ! === Accumulate ===
           do iomega=istart,iend
             xx=(omega(iomega)-trans)*gaussfactor
             gaussval=gaussprefactor*EXP(-xx*xx)
             jdos(iomega,isppol,iqibz)=jdos(iomega,isppol,iqibz)+gaussval!+Kmesh%wtk(ikibz)*gaussval
           end do

         end do !ib1
       end do !ib2
     end do !isppol
   end do !ikbz

 end do !iqibz
 !
 ! === Write results ===
 !fnam='JDOS' ; call isfile(fname,'new')
 unt=get_unit()
 open(file=fname,unit=unt,form='formatted')
 write(unt,'(a)')         '# Joint density of states. All in eV units '
 write(unt,'(a,es16.8,a)')'# Frequency step ',step*Ha_eV,' [eV]'
 write(unt,'(a,es16.8,a)')'# Max Frequency  ',max_transition*Ha_eV,' [eV]'
 write(unt,'(a,i2)')      '# Number of Independendent polarization',BSt%nsppol
 write(unt,'(a,i5)')      '# Number of k-points in BZ ',Kmesh%nbz
 write(unt,'(a,i5)')      '# Number of analysed q-points',nqibz
 do iqibz=1,nqibz
   write(unt,'(a,i4,a,3es16.8)')'# ',iqibz,') ',qcalc(:,iqibz)
 end do
 write(unt,'(a)')'#'
 write(unt,'(a)')'# do iomega=1,nomega '
 write(unt,'(a)')'#  write(std_out,*) omega(iomega),((jdos(iomega,isppol,iqibz),isppol=1,nsppol),iqibz=1,nqibz)'
 write(unt,'(a)')'# end do '
 write(unt,'(a)')'#'
 do iomega=1,nomega
   write(unt,'(es16.8)')omega(iomega),((jdos(iomega,isppol,iqibz),isppol=1,BSt%nsppol),iqibz=1,nqibz)
 end do
 close(unt)

 nullify(qcalc)

 ! === Free memory ===
 ABI_DEALLOCATE(omega)
 ABI_DEALLOCATE(jdos)

 call destroy_Little_group(Ltg_q)
 call destroy_BZ_mesh_type(Kmesh)
 if (SIZE(qibz,DIM=2)==0) then
   call destroy_BZ_mesh_type(Qmesh)
 end if

 DBG_EXIT("COLL")

end subroutine joint_dos
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/abi_bands_put
!! NAME
!! abi_bands_put
!!
!! FUNCTION
!!  Writes the content of a bandstructure_type object to a NETCDF file 
!!  according to the ETSF-IO specifications.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_sigma_results
!!
!! CHILDREN
!!
!! SOURCE

subroutine abi_bands_put(Bst,fname)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_bands_put'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: fname
 type(bandstructure_type),target,intent(in) :: Bst
!arrays

!Local variables-------------------------------
!scalars
 integer :: untbst,usewvl,formeig,idx,isppol,nkpt,iband,ikpt
 integer,target :: nelect_int
 real(dp) :: fixmom_
 logical :: lstat
 character(len=500) :: msg
#if defined HAVE_TRIO_ETSF_IO
 character(len=etsf_io_low_error_len) :: errmess
 !type(ETSF_dims) :: Dims
 type(ETSF_io_low_error) :: Error_data
 type(ETSF_basisdata) :: Basisdata
 type(ETSF_groups) :: GroupFolder
 type(ETSF_electrons),target :: Electrons
 type(ETSF_kpoints),target :: Kpoints
#endif
!arrays
 real(dp),allocatable,target :: eig_vec(:),occ_vec(:)

! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_TRIO_ETSF_IO
 ! === Open the file ===
 call wrtout(std_out,' abi_bands_put : about to modify file '//TRIM(fname),'COLL')

 call etsf_io_low_open_modify(untbst,fname,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 ! === Write dimensions handled by ETSF ===
 ! * Use low level procedures to write dimensions, since 
 !   etsf_io_dims_def correctly complains that some ETSF dims are missing.

 !FIXME: do not handle k_dependent = 1
 !Dims%max_number_of_states        = Bst%mband    
 !Dims%number_of_kpoints           = Bst%nkpt
 !Dims%number_of_spinor_components = Bst%nspinor  
 !Dims%number_of_spins             = Bst%nsppol  
 !!Dims%number_of_components       = Bst%nspden 
 !call etsf_io_dims_def(untbst,Dims,lstat,Error_data)
 !if (.not.lstat) goto 1000

 call etsf_io_low_write_dim(untbst,'max_number_of_states',Bst%mband,lstat,Error_data=Error_data)  
 if (.not.lstat) goto 1000  
 call etsf_io_low_write_dim(untbst,'number_of_kpoints',Bst%nkpt,lstat,Error_data=Error_data)  
 if (.not.lstat) goto 1000  
 call etsf_io_low_write_dim(untbst,'number_of_spinor_components',Bst%nspinor,lstat,Error_data=Error_data)  
 if (.not.lstat) goto 1000  
 call etsf_io_low_write_dim(untbst,'number_of_spins',Bst%nsppol,lstat,Error_data=Error_data)  
 if (.not.lstat) goto 1000  


 usewvl=0 !FIXME
 ! Get all variables included in ETSF
 if (usewvl==0) then
  BasisData%number_of_coefficients => Bst%npwarr
  call etsf_io_basisdata_put(untbst,Basisdata,lstat,Error_data)
  if (.not.lstat) goto 1000
 end if

 ! === Close the file, due to call to etsf_io_data_write  ===
 call etsf_io_low_close(untbst,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 Electrons%fermi_energy            => Bst%fermie
 nelect_int=NINT(Bst%nelect)
 Electrons%number_of_electrons     => nelect_int  !FIXME in ETSF this is integer
 Electrons%smearing_width          => Bst%tsmear 
 Electrons%number_of_states%data1D => Bst%nband 

 ! TODO DFPT not treated, we allocate correctly but the object should be modified a bit. 
 formeig=0
 ABI_ALLOCATE(eig_vec,((2*Bst%mband)**formeig*Bst%mband*Bst%nkpt*Bst%nsppol))
 ABI_ALLOCATE(occ_vec,(Bst%mband*Bst%nkpt*Bst%nsppol))
 idx=0 !call my helper function once modules will be supported.
 do isppol=1,Bst%nsppol
  do ikpt=1,Bst%nkpt
   do iband=1,Bst%nband(ikpt+(isppol-1)*nkpt)
    idx=idx+1
    eig_vec(idx)=Bst%eig(iband,ikpt,isppol)
    occ_vec(idx)=Bst%occ(iband,ikpt,isppol)
   end do
  end do
 end do

 Electrons%eigenvalues%data1D => eig_vec
 Electrons%occupations%data1D => occ_vec

 Kpoints%reduced_coordinates_of_kpoints => Bst%kptns
 Kpoints%kpoint_weights                 => Bst%wtk

 GroupFolder%Electrons => Electrons
 GroupFolder%Kpoints   => Kpoints

 call etsf_io_data_write(fname,GroupFolder,lstat,Error_data)

 ABI_DEALLOCATE(eig_vec)
 ABI_DEALLOCATE(occ_vec)

 ! === Write additional stuff contained in the the abinit header ===
 ! * Define additional variables.
 call etsf_io_low_open_modify(untbst,fname,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 call etsf_io_low_def_var(untbst,'tphysel',etsf_io_low_double,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 call etsf_io_low_def_var(untbst,'occopt',etsf_io_low_integer,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 call etsf_io_low_def_var(untbst,'istwfk',etsf_io_low_integer,(/'number_of_kpoints'/),lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 ! === Write data ===
 call etsf_io_low_set_write_mode(untbst,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 call etsf_io_low_write_var(untbst,'occopt',Bst%occopt,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 call etsf_io_low_write_var(untbst,'tphysel',Bst%tphysel,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 call etsf_io_low_write_var(untbst,'istwfk',Bst%istwfk,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 ! === Close the file ===
 call etsf_io_low_close(untbst,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 1000 continue
 ! === Handle the error ===
 if (.not.lstat) then
  call etsf_io_low_error_to_str(errmess,Error_data)
  msg=errmess(1:min(500,len(errmess)))
  MSG_ERROR(msg)
 end if

 ! === Finalize the object ===
 !FIXME At the moment this has to be done outside the routine
 ! Moreover I have to solve the problem with the convention for the Fermi level
 ! Maybe it makes sense if fixmom is included in the object!
 fixmom_=-99.99_dp 
 !call update_occ(BSt,fixmom_)

#else 
 MSG_ERROR("ETSF-IO support is not activated. ")
#endif

 DBG_EXIT("COLL")

end subroutine abi_bands_put
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/abi_Bands_read
!! NAME
!! abi_Bands_read
!!
!! FUNCTION
!!  Initialize an instance of the bandstructure_type from a NETCDF file
!!  written according to the ETSF-IO specifications.
!!
!! INPUTS
!!
!! OUTPUT
!!  Data written in file whose name is fname.
!!
!! PARENTS
!!      m_sigma_results
!!
!! CHILDREN
!!
!! SOURCE

subroutine abi_bands_read(Bst,fname)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_bands_read'
 use interfaces_14_hidewrite
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: fname
 type(bandstructure_type),target,intent(out) :: Bst
!arrays

!Local variables-------------------------------
!scalars
 integer :: ncid,fform,usewvl,formeig,idx,isppol,nkpt,iband,ikpt
 integer,target :: nelect_int
 real(dp) :: fixmom_
 logical :: lstat
#if defined HAVE_TRIO_ETSF_IO
 character(len=etsf_io_low_error_len) :: errmess
 type(ETSF_dims) :: Dims
 type(ETSF_io_low_error) :: Error_data
 type(ETSF_basisdata) :: Basisdata
 type(ETSF_groups) :: GroupFolder
 type(ETSF_electrons),target :: Electrons
 type(ETSF_kpoints),target :: Kpoints
#endif
 type(Hdr_type) :: Hdr
 !character(len=500) :: msg
!arrays
 real(dp),allocatable,target :: eig_vec(:),occ_vec(:)

! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_TRIO_ETSF_IO
 ! === Open the file ===
 call wrtout(std_out,' abi_Bands_read : about to read file '//TRIM(fname),'COLL')

 call etsf_io_low_open_read(ncid,fname,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 ! === Read dimensions handled by ETSF ===
 call etsf_io_dims_get(ncid,Dims,lstat,Error_data)
 if (.not.lstat) goto 1000

 ! TODO : do not handle k_dependent = 1
 Bst%bantot   = Dims%max_number_of_states * Dims%number_of_kpoints * Dims%number_of_spins
 Bst%mband    = Dims%max_number_of_states
 Bst%nkpt     = Dims%number_of_kpoints
 Bst%nspinor  = Dims%number_of_spinor_components
 Bst%nsppol   = Dims%number_of_spins
 !? Bst%nspden   = Dims%number_of_components

 ABI_ALLOCATE(Bst%istwfk,(Bst%nkpt))
 ABI_ALLOCATE(Bst%nband,(Bst%nkpt*Bst%nsppol))
 ABI_ALLOCATE(Bst%npwarr,(Bst%nkpt))
 ABI_ALLOCATE(Bst%kptns,(3,Bst%nkpt))
 ABI_ALLOCATE(Bst%eig,(Bst%mband,Bst%nkpt,Bst%nsppol))
 ABI_ALLOCATE(Bst%occ,(Bst%mband,Bst%nkpt,Bst%nsppol))
 ABI_ALLOCATE(Bst%doccde,(Bst%mband,Bst%nkpt,Bst%nsppol))
 ABI_ALLOCATE(Bst%wtk,(Bst%nkpt))

 usewvl=0
 call etsf_io_low_read_var(ncid,'usewvl',usewvl,lstat,error_data=Error_data)
 if (.not.lstat) goto 1000

 ! Get all variables included in ETSF
 if (usewvl==0) then
  BasisData%number_of_coefficients => Bst%npwarr
  call etsf_io_basisdata_get(ncid,Basisdata,lstat,Error_data)
  if (.not.lstat) goto 1000
 end if

 Electrons%fermi_energy            => Bst%fermie
 !FIXME this is integer, therefore alchemy or charged cells won"t work!
 Electrons%number_of_electrons     => nelect_int 
 Electrons%smearing_width          => Bst%tsmear 
 Electrons%number_of_states%data1D => Bst%nband 

 !TODO DFPT not treated, we read correctly but the object should be modified a bit. 
 formeig=0
 ABI_ALLOCATE(eig_vec,((2*Bst%mband)**formeig*Bst%mband*Bst%nkpt*Bst%nsppol))
 ABI_ALLOCATE(occ_vec,(Bst%mband*Bst%nkpt*Bst%nsppol))
 Electrons%eigenvalues%data1D      => eig_vec  !then we have to unpack
 Electrons%occupations%data1D      => occ_vec

 Kpoints%reduced_coordinates_of_kpoints => Bst%kptns
 Kpoints%kpoint_weights                 => Bst%wtk

 GroupFolder%Electrons => Electrons
 GroupFolder%Kpoints   => Kpoints

 call etsf_io_data_read(fname,GroupFolder,lstat,Error_data)

 idx=0 ! TODO call my helper function once modules will be supported.
 do isppol=1,Bst%nsppol
  do ikpt=1,Bst%nkpt
   do iband=1,Bst%nband(ikpt+(isppol-1)*nkpt)
    idx=idx+1
    Bst%eig(iband,ikpt,isppol)=eig_vec(idx)
    Bst%occ(iband,ikpt,isppol)=occ_vec(idx)
   end do
  end do
 end do
 ABI_DEALLOCATE(eig_vec)
 ABI_DEALLOCATE(occ_vec)

 ! === Read the abinit header ===
 ! * Fill missing quantities in Bst using Hdr.
 call hdr_io_etsf(fform,Hdr,1,ncid)

 Bst%tphysel = Hdr%tphysel
 Bst%istwfk  = Hdr%istwfk
 Bst%occopt  = Hdr%occopt
 call hdr_clean(Hdr)

 ! === Close the file ===
 call etsf_io_low_close(ncid,lstat,Error_data=Error_data)
 if (.not.lstat) goto 1000

 ! === Finalize the object ===
 !FIXME At the moment this has to be done outside the routine
 ! Moreover I have to solve the problem with the convention for the Fermi level
 fixmom_=99.99_dp 
 !call update_occ(BSt,fixmom_)

 ! === Handle the error ===
 1000 continue
 if (.not.lstat) then
  call etsf_io_low_error_to_str(errmess,Error_data)
  MSG_ERROR(errmess)
 end if

#else 
 MSG_ERROR("ETSF-IO support is not activated.")
#endif
 
 DBG_EXIT("COLL")

end subroutine abi_Bands_read
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/xmgrace_bandplot
!! NAME
!! xmgrace_bandplot
!!
!! FUNCTION
!! Plots the interpolated band structure in Xmgrace format  
!! Based on plot_interpolate_xmgrace, routine of the Wannier90 code.
!!
!! INPUTS
!!  fname=The name of the file.
!!
!! Copyright (C) 2007 Jonathan Yates, Arash Mostofi,          
!!  Young-Su Lee, Nicola Marzari, Ivo Souza, David Vanderbilt 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine xmgrace_bandplot(mband,nkpt,nsppol,kptns,eig,gmet,fname,bands_label,ierr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xmgrace_bandplot'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,mband,nsppol
 integer,intent(out) :: ierr
 character(len=*),intent(in) :: fname
!arrays
 real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt),eig(mband,nkpt,nsppol)
 character(len=1),optional,intent(in) :: bands_label(nkpt)

!Local variables-------------------------------
 integer :: ikpt,ikpath,ib,xunit 
 integer :: nspec_k
 real(dp) :: emin,emax
 !character(len=9) :: cdate,ctime
 character(len=fnlen) :: xname
!arrays 
 real(dp) :: kdiff(3)
 real(dp) :: xval(nkpt)
 character(len=10),allocatable  :: ctemp(:) !,xlabel(:)
! *************************************************************************

 ABI_CHECK(nsppol==1,"nsppol=2 not coded")

 !call io_date(cdate, ctime)
 ierr=0

 ! === Axis labels ===
 ! TODO doesnt work yet
 ABI_ALLOCATE(ctemp,(nkpt))
 ctemp='0'

 ! Switch any G to Gamma
 if (PRESENT(bands_label)) then
   do ikpt=1,nkpt 
    if (bands_label(ikpt)/='0') then
     ctemp(ikpt)=bands_label(ikpt)
     if (ctemp(ikpt)=='G') ctemp(ikpt)='\xG\0'
    end if
   end do
 end if

 nspec_k=COUNT(ctemp/='0')

 !allocate(xlabel(nspec_k))
 !xlabel(1)=' '//TRIM(ctemp(1))//' '
 !do i=2,num_paths
 ! if(ctemp(2*(i-1))/=ctemp(2*(i-1)+1)) then
 !  xlabel(i)=ctemp(2*(i-1))//'/'//ctemp(2*(i-1)+1)
 ! else
 !  xlabel(i)=ctemp(2*(i-1))
 ! end if
 !end do
 !xlabel(num_spts)=ctemp(bands_num_spec_points)
 !deallocate(ctemp)

 ! Get min and Max energy in Ev
 !minmax_ene = get_minmax(BSt,'eig')*Ha_eV
 emin = MINVAL(eig)*Ha_eV - one
 emax = MAXVAL(eig)*Ha_eV + one

 xval(1)=zero
 do ikpath=2,nkpt
   kdiff=kptns(:,ikpath)-kptns(:,ikpath-1)
   xval(ikpath) = xval(ikpath-1)+ normv(kdiff,gmet,'G')
 end do

 ! Xmgrace format
 xunit=get_unit() 
 xname=TRIM(fname)//'_band.agr'
 open(xunit,file=xname,form='formatted')

 write(xunit,'(a)') '# Grace project file                      '
 write(xunit,'(a)') '# written using Abinit www.abinit.org     '
 write(xunit,'(a)') '@version 50113                            '
 write(xunit,'(a)') '@page size 792, 612                       '
 write(xunit,'(a)') '@page scroll 5%                           '
 write(xunit,'(a)') '@page inout 5%                            '
 write(xunit,'(a)') '@link page off                            '

 !write(xunit,'(a)') '@timestamp def "'//cdate//' at '//ctime//'" ' 
 write(xunit,'(a)') '@with g0'                                  
 write(xunit,'(a)') '@    world xmin 0.00'
 !write(xunit,'(a,f10.5)') '@    world xmax ',xval(nkpt)
 write(xunit,'(a,i5)') '@    world xmax ',nkpt-1
 write(xunit,'(a,f10.5)') '@    world ymin ',emin
 write(xunit,'(a,f10.5)') '@    world ymax ',emax
 write(xunit,'(a)') '@default linewidth 1.5'
 write(xunit,'(a)') '@    xaxis  tick on'
 write(xunit,'(a)') '@    xaxis  tick major 1'
 write(xunit,'(a)') '@    xaxis  tick major color 1'
 write(xunit,'(a)') '@    xaxis  tick major linestyle 3'
 write(xunit,'(a)') '@    xaxis  tick major grid on'
 write(xunit,'(a)') '@    xaxis  tick spec type both'
 !write(xunit,'(a,i0)') '@    xaxis  tick spec ',1+bands_num_spec_points/2
 write(xunit,'(a)') '@    xaxis  tick major 0, 0'
 !do i=1,bands_num_spec_points/2
 ! write(xunit,'(a,i0,a,a)') '@    xaxis  ticklabel ',i-1,',', '"'//trim(adjustl(xlabel(i)))//'"'
 ! write(xunit,'(a,i0,a,f10.5)') '@    xaxis  tick major ',i,' , ',sum(kpath_len(1:i))
 !end do
 !write(xunit,'(a,i0,a)') '@    xaxis  ticklabel ',bands_num_spec_points/2 &
 ! ,',"'//trim(adjustl(xlabel(1+bands_num_spec_points/2)))//'"'
 write(xunit,'(a)') '@    xaxis  ticklabel char size 1.500000'
 write(xunit,'(a)') '@    yaxis  tick major 10'
 write(xunit,'(a)') '@    yaxis  label "Band Energy (eV)"'
 write(xunit,'(a)') '@    yaxis  label char size 1.500000'
 write(xunit,'(a)') '@    yaxis  ticklabel char size 1.500000'

 !here I suppose no dependency on k and spin
 do ib=1,mband !num_wann
   write(xunit,'(a,i4,a)') '@    s',ib-1,' line color 1'
 end do

 do ib=1,mband !num_wann
   write(xunit,'(a,i4)') '@target G0.S',ib-1
   write(xunit,'(a)') '@type xy'
   do ikpt=1,nkpt
     !write(xunit,'(2e16.8)')xval(ikpt),eig(ib,ikpt,1)*Ha_eV !FIXME doesnt work for nsppol=2
     write(xunit,'(i5,e16.8)')ikpt-1,eig(ib,ikpt,1)*Ha_eV !FIXME doesnt work for nsppol=2
   end do
   write(xunit,'(a)')'&'
 end do

 close(xunit)

end subroutine xmgrace_bandplot 
!!***

END MODULE m_ebands
!!***
