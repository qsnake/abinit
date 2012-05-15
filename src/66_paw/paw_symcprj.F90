!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_symcprj
!! NAME
!!  paw_symcprj
!!
!! FUNCTION
!!
!! INPUTS
!!  ik_ibz=The index of the k-point in the full BZ where the matrix elements have to be symmetrized.
!!  nspinor=Number of spinorial components
!!  nband_k=Number of bands stored in cprjnk_kibz for this k-point.
!!  Cryst<Crystal_structure>=data type gathering information on unit cell and symmetries.
!!     %ntypat=number of type of atoms
!!     %natom=number of atoms in the unit cell
!!     %typat(natom)=type of each atom
!!     %indsym(4,nsym,natom)=indirect indexing array:
!!      for each isym,iatom, fourth element is label of atom into which iatom is sent by the INVERSE of the
!!      symmetry operation symrel(isym); first three elements are the primitive translations that must be subtracted
!!      after the transformation to get back to the original unit cell.
!!  Kmesh<bz_mesh_type>: datatype gathering information on the k-point sampling.
!!     %nbz=number of k-points in the full Brillouin zone
!!     %nibz=number of k-points in the irreducible wedge
!!     %tab(nkbz)=table giving for each k-point in the BZ (array kbz), the corresponding irred. point in the IBZ.
!!       i.e k_BZ = (IS) kIBZ where S is one of the symrec operations and I is the inversion or the identity
!!     %tabi(nkbz)=for each k-point in the BZ defines whether inversion has to be considered in the
!!       relation k_BZ=(IS) k_IBZ (1 => only S; -1 => -S)
!!     %tabo(nkbz)= the symmetry operation S that takes k_IBZ to each k_BZ
!!  Pawtab(Cryst%ntypat) <type(pawtab_type)>=paw tabulated starting data.
!!  Pawang <type(pawang_type)>=paw angular mesh and related data
!!     %lmax=Max angular momentum included in the PAW datasets used. mentioned at the second line of the psp file
!!     %zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)=coefficients of the transformation of real spherical
!!      harmonics under the symmetry operations.
!!  Psps<pseudopotential_type>
!!     %indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!
!! NOTES  
!!  Derivatives are not symmetrized.
!!
!! OUTPUT
!! 
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,cchi0q0_intraband
!!      check_completeness,cohsex_me,debug_tools,m_shirley
!!
!! CHILDREN
!!      cprj_copy,get_bz_item,get_ibz_item
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine paw_symcprj(ik_bz,nspinor,nband_k,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj_bz) 

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

 use m_crystal,   only : crystal_structure
 use m_bz_mesh,   only : bz_mesh_type, get_ibz_item, get_bz_item

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_symcprj'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nspinor,nband_k,ik_bz
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)
 type(Cprj_type),intent(inout) :: Cprj_bz(Cryst%natom,nspinor*nband_k)

!Local variables-------------------------------
!scalars
 integer :: iatom,iat_sym,iband,ibsp_bz
 integer :: ibsp_ibz,ik_ibz,indexj,ispinor,isym,itim
 integer :: itypat,jl,jl0,jlmn,jln,jln0,jlpm,jm,jn,lmax,mm,ncpgr
 real(dp) :: arg,wtk
 logical :: isirred
!arrays
 integer :: r0(3),nlmn_atom(Cryst%natom)
 real(dp) :: dum(2,nspinor),kbz(3),kirr(3),phase(2),swp(2),tmp(2,nspinor)
 real(dp),allocatable :: DS_mmpl(:,:,:)
 type(Cprj_type) :: Cprjnk_kibz(Cryst%natom,nspinor*nband_k)

! *********************************************************************

 ncpgr = Cprj_bz(1,1)%ncpgr
 ABI_CHECK(ncpgr==0,"Derivatived of cprj are not coded")

!Get the index of the IBZ image associated to the BZ k-point ik_bz and related simmetry.
 call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym,itim,isirred=isirred)

 if (isirred) RETURN  ! It is a point in the IBZ, Symmetrization is not needed.
!
!The corresponding point kirr in the IBZ.
 call get_IBZ_item(Kmesh,ik_ibz,kirr,wtk)

!Local copy.
 do iatom=1,Cryst%natom
   nlmn_atom(iatom)=Pawtab(Cryst%typat(iatom))%lmn_size
 end do
 
 call cprj_alloc(Cprjnk_kibz,ncpgr,nlmn_atom)
 call cprj_copy(Cprj_bz,Cprjnk_kibz)
!
!=== DS_mmpl is the rotation matrix for real spherical harmonics associated to symrec(:,:,isym) ===
!* Note the convention used by Blanco in Eq. 27 : DS_mmp multiply spherical harmonics as row vectors
 lmax=Pawang%l_max-1 ! l_max is Max l+1
 ABI_ALLOCATE(DS_mmpl,(2*lmax+1,2*lmax+1,lmax+1))
 DS_mmpl=Pawang%zarot(:,:,:,isym)
!
!===========================================
!==== Loop over atoms to be symmetrized ====
!===========================================
 do iatom=1,Cryst%natom
   itypat=Cryst%typat(iatom)
   iat_sym=Cryst%indsym(4,isym,iatom)
   r0=Cryst%indsym(1:3,isym,iatom) ! R^{-1} (xred(:,iatom)-tnons) = xred(:,iat_sym) + r0.
   arg=two_pi*dot_product(kirr,r0)
   phase(1)=COS(arg)
   phase(2)=SIN(arg)
!  
!  Loop over the (jl,jm,jn) components to be symmetrized.
   jl0=-1; jln0=-1; indexj=1
   do jlmn=1,Pawtab(itypat)%lmn_size
     jl  =Psps%indlmn(1,jlmn,itypat)
     jm  =Psps%indlmn(2,jlmn,itypat)
     jn  =Psps%indlmn(3,jlmn,itypat)
     jln =Psps%indlmn(5,jlmn,itypat)
     jlpm=1+jl+jm
     if (jln/=jln0) indexj=indexj+2*jl0+1
!    
!    === For each band, calculate contribution due to rotated real spherical harmonics ===
!    FIXME check this expression; according to Blanco I should have D(S^-1} but it seems D(S) is correct
!    Recheck spinorial case, presently is wrong
     ibsp_ibz=0
     ibsp_bz=0
     do iband=1,nband_k

       tmp(:,:)=zero
       do ispinor=1,nspinor
         ibsp_ibz=ibsp_ibz+1
         do mm=1,2*jl+1
           tmp(1,ispinor)=tmp(1,ispinor)+DS_mmpl(mm,jlpm,jl+1)*Cprjnk_kibz(iat_sym,ibsp_ibz)%cp(1,indexj+mm)
           tmp(2,ispinor)=tmp(2,ispinor)+DS_mmpl(mm,jlpm,jl+1)*Cprjnk_kibz(iat_sym,ibsp_ibz)%cp(2,indexj+mm)
         end do
       end do !ispinor
!      
!      * Apply the phase to account if the symmetric atom belongs to a different unit cell.
       do ispinor=1,nspinor
         dum(1,ispinor)=tmp(1,ispinor)*phase(1)-tmp(2,ispinor)*phase(2)
         dum(2,ispinor)=tmp(1,ispinor)*phase(2)+tmp(2,ispinor)*phase(1)
       end do
!      
!      * If required, apply time-reversal symmetry to retrieve the correct point in the BZ.
       if (itim==2) then
         if (nspinor==1) then
           dum(2,1)=-dum(2,1)
         else if (nspinor==2) then ! TODO rotate wavefunction in spinor space.
           swp(:)=dum(:,1)
           dum(1,1)= dum(1,2)
           dum(2,1)=-dum(2,2)
           dum(1,2)=-swp(1)
           dum(2,2)= swp(2)
         end if
       end if
!      
!      ==== Save values ====
       do ispinor=1,nspinor
         ibsp_bz=ibsp_bz+1
         Cprj_bz(iatom,ibsp_bz)%cp(1,jlmn)=dum(1,ispinor)
         Cprj_bz(iatom,ibsp_bz)%cp(2,jlmn)=dum(2,ispinor)
       end do
     end do !iband

     jl0=jl; jln0=jln
   end do !jlmn
 end do !iatom

 call cprj_free(Cprjnk_kibz)
 ABI_DEALLOCATE(DS_mmpl)

end subroutine paw_symcprj
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/paw_symcprj_op
!! NAME
!!  paw_symcprj_op
!!
!! FUNCTION
!!
!! INPUTS
!!  ik_ibz=The index of the k-point in the full BZ where the matrix elements have to be symmetrized.
!!  nspinor=Number of spinorial components
!!  nband_k=Number of bands stored in cprjnk_kibz for this k-point.
!!  Cryst<Crystal_structure>=data type gathering information on unit cell and symmetries.
!!     %ntypat=number of type of atoms
!!     %natom=number of atoms in the unit cell
!!     %typat(natom)=type of each atom
!!     %indsym(4,nsym,natom)=indirect indexing array:
!!      for each isym,iatom, fourth element is label of atom into which iatom is sent by the INVERSE of the
!!      symmetry operation symrel(isym); first three elements are the primitive translations that must be subtracted
!!      after the transformation to get back to the original unit cell.
!!  Kmesh<bz_mesh_type>: datatype gathering information on the k-point sampling.
!!     %nbz=number of k-points in the full Brillouin zone
!!     %nibz=number of k-points in the irreducible wedge
!!     %tab(nkbz)=table giving for each k-point in the BZ (array kbz), the corresponding irred. point in the IBZ.
!!       i.e k_BZ = (IS) kIBZ where S is one of the symrec operations and I is the inversion or the identity
!!     %tabi(nkbz)=for each k-point in the BZ defines whether inversion has to be considered in the
!!       relation k_BZ=(IS) k_IBZ (1 => only S; -1 => -S)
!!     %tabo(nkbz)= the symmetry operation S that takes k_IBZ to each k_BZ
!!  Pawtab(Cryst%ntypat) <type(pawtab_type)>=paw tabulated starting data.
!!  Pawang <type(pawang_type)>=paw angular mesh and related data
!!     %lmax=Max angular momentum included in the PAW datasets used. mentioned at the second line of the psp file
!!     %zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)=coefficients of the transformation of real spherical
!!      harmonics under the symmetry operations.
!!  Psps<pseudopotential_type>
!!     %indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!   in_Cprj(Cryst%natom,nspinor*nband_k)<Cprj_type>=Input cprj
!!
!! OUTPUT
!!   out_Cprj(Cryst%natom,nspinor*nband_k)<Cprj_type>=Symmetrized cprj matrix elements.
!!
!! NOTES  
!!  Derivatives are not symmetrized.
!! 
!! PARENTS
!!      exc_build_block
!!
!! CHILDREN
!!      cprj_copy,get_bz_item,get_ibz_item
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine paw_symcprj_op(ik_bz,nspinor,nband_k,Cryst,Kmesh,Psps,Pawtab,Pawang,in_Cprj,out_Cprj) 

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

 use m_crystal,   only : crystal_structure
 use m_bz_mesh,   only : bz_mesh_type, get_ibz_item, get_bz_item

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_symcprj_op'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nspinor,nband_k,ik_bz
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)
 type(Cprj_type),intent(in) :: in_Cprj(Cryst%natom,nspinor*nband_k)
 type(Cprj_type),intent(out) :: out_Cprj(Cryst%natom,nspinor*nband_k)

!Local variables-------------------------------
!scalars
 integer :: iatom,iat_sym,iband,ibsp_bz
 integer :: ibsp_ibz,ik_ibz,indexj,ispinor,isym,itim
 integer :: itypat,jl,jl0,jlmn,jln,jln0,jlpm,jm,jn,lmax,mm,ncpgr
 real(dp) :: arg,wtk
 logical :: isirred
!arrays
 integer :: r0(3) !,nlmn_atom(Cryst%natom)
 real(dp) :: dum(2,nspinor),kbz(3),kirr(3),phase(2),swp(2),tmp(2,nspinor)
 real(dp),allocatable :: DS_mmpl(:,:,:)

! *********************************************************************

 ncpgr = in_Cprj(1,1)%ncpgr
 ABI_CHECK(ncpgr==0,"Derivatived of cprj are not coded")

!Get the index of the IBZ image associated to the BZ k-point ik_bz and related simmetry.
 call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym,itim,isirred=isirred)

 if (isirred) then  ! It is a point in the IBZ, Symmetrization is not needed.
   call cprj_copy(in_Cprj,out_Cprj)
   RETURN  
 end if
!
!The corresponding point kirr in the IBZ.
 call get_IBZ_item(Kmesh,ik_ibz,kirr,wtk)
!
!=== DS_mmpl is the rotation matrix for real spherical harmonics associated to symrec(:,:,isym) ===
!* Note the convention used by Blanco in Eq. 27 : DS_mmp multiply spherical harmonics as row vectors
 lmax=Pawang%l_max-1 ! l_max is Max l+1
 ABI_ALLOCATE(DS_mmpl,(2*lmax+1,2*lmax+1,lmax+1))
 DS_mmpl=Pawang%zarot(:,:,:,isym)

!Local copy.
!do iatom=1,Cryst%natom
!nlmn_atom(iatom)=Pawtab(Cryst%typat(iatom))%lmn_size
!end do
!call cprj_alloc(out_Cprj,ncpgr,nlmn_atom)
!
!===========================================
!==== Loop over atoms to be symmetrized ====
!===========================================
 do iatom=1,Cryst%natom
   itypat=Cryst%typat(iatom)
   iat_sym=Cryst%indsym(4,isym,iatom)
   r0=Cryst%indsym(1:3,isym,iatom) ! R^{-1} (xred(:,iatom)-tnons) = xred(:,iat_sym) + r0.
   arg=two_pi*dot_product(kirr,r0)
   phase(1)=COS(arg)
   phase(2)=SIN(arg)
!  
!  Loop over the (jl,jm,jn) components to be symmetrized.
   jl0=-1; jln0=-1; indexj=1
   do jlmn=1,Pawtab(itypat)%lmn_size
     jl  =Psps%indlmn(1,jlmn,itypat)
     jm  =Psps%indlmn(2,jlmn,itypat)
     jn  =Psps%indlmn(3,jlmn,itypat)
     jln =Psps%indlmn(5,jlmn,itypat)
     jlpm=1+jl+jm
     if (jln/=jln0) indexj=indexj+2*jl0+1
!    
!    === For each band, calculate contribution due to rotated real spherical harmonics ===
!    FIXME check this expression; according to Blanco I should have D(S^-1} but it seems D(S) is correct
!    Recheck spinorial case, presently is wrong
     ibsp_ibz=0
     ibsp_bz=0
     do iband=1,nband_k

       tmp(:,:)=zero
       do ispinor=1,nspinor
         ibsp_ibz=ibsp_ibz+1
         do mm=1,2*jl+1
           tmp(1,ispinor)=tmp(1,ispinor)+DS_mmpl(mm,jlpm,jl+1)*in_Cprj(iat_sym,ibsp_ibz)%cp(1,indexj+mm)
           tmp(2,ispinor)=tmp(2,ispinor)+DS_mmpl(mm,jlpm,jl+1)*in_Cprj(iat_sym,ibsp_ibz)%cp(2,indexj+mm)
         end do
       end do !ispinor
!      
!      * Apply the phase to account if the symmetric atom belongs to a different unit cell.
       do ispinor=1,nspinor
         dum(1,ispinor)=tmp(1,ispinor)*phase(1)-tmp(2,ispinor)*phase(2)
         dum(2,ispinor)=tmp(1,ispinor)*phase(2)+tmp(2,ispinor)*phase(1)
       end do
!      
!      * If required, apply time-reversal symmetry to retrieve the correct point in the BZ.
       if (itim==2) then
         if (nspinor==1) then
           dum(2,1)=-dum(2,1)
         else if (nspinor==2) then ! TODO rotate wavefunction in spinor space.
           swp(:)=dum(:,1)
           dum(1,1)= dum(1,2)
           dum(2,1)=-dum(2,2)
           dum(1,2)=-swp(1)
           dum(2,2)= swp(2)
         end if
       end if
!      
!      ==== Save values ====
       do ispinor=1,nspinor
         ibsp_bz=ibsp_bz+1
         out_Cprj(iatom,ibsp_bz)%cp(1,jlmn)=dum(1,ispinor)
         out_Cprj(iatom,ibsp_bz)%cp(2,jlmn)=dum(2,ispinor)
       end do
     end do !iband

     jl0=jl; jln0=jln
   end do !jlmn
 end do !iatom

 ABI_DEALLOCATE(DS_mmpl)

end subroutine paw_symcprj_op
!!***
