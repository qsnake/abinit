!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_cprj_bspline
!! NAME
!! m_cprj_bspline
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_cprj_bspline

 use m_profiling

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_bspline
 use m_errors 

 use m_numeric_tools,   only : arth
 use m_crystal,         only : crystal_structure

 implicit none

 private
!!***

!----------------------------------------------------------------------
!!***

!!****t* m_cprj_bspline/cprj_bspl_t
!! NAME
!!  cprj_bspl_t
!! 
!! FUNCTION
!! 
!! SOURCE

 type,public :: cprj_bspl_t 
   integer :: natom
   integer :: nspinor
   integer :: blc_size
   integer :: nkx,nky,nkz,nk_tot
   integer :: kxord, kyord, kzord
   integer :: nxknot, nyknot,nzknot

   integer,pointer :: nlmn_sort(:)

   real(dp),pointer :: xvec(:),yvec(:),zvec(:),kvec(:,:)

   real(dp),pointer :: xknot(:),yknot(:),zknot(:)

   type(coeff5_type),pointer :: Bscoef(:,:,:)
   ! Bscoef(ntypat,nspinor,blc_size)

   type(cprj_type),pointer :: Pblc_aki(:,:,:,:)
 end type cprj_bspl_t
!!***

!----------------------------------------------------------------------

 public  :: cprj_bspline_init
 public  :: cprj_bspline_free
 public  :: cprj_bspline_mktensor
 public  :: cprj_bspline_eval
 public  :: pbloch

CONTAINS  !====================================================================

!----------------------------------------------------------------------

!!****f* m_cprj_bspline/cprj_bspline_init
!! NAME
!!  cprj_bspline_init
!!
!! FUNCTION
!!  Creation method for the cprj_bspl_t structured datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      initmpi_seq,initylmg,mkffnl,mkkpg,nonlop_ylm,ph1d3d
!!
!! SOURCE

subroutine cprj_bspline_init(Cp_bspl,Cryst,Psps,Pawtab,Pawang,nspinor,ngfft,kdiv,kord,npw,blc_size,kg_k,blc_ug,nlmn_sort) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_bspline_init'
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nspinor,blc_size,npw
 type(Crystal_structure),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawang_type),intent(in) :: Pawang
 type(cprj_bspl_t),intent(inout) :: Cp_bspl
!arrays
 integer,intent(in) :: kdiv(3),nlmn_sort(Cryst%natom),kord(3),ngfft(18),kg_k(3,npw)
 complex(gwpc),intent(in) :: blc_ug(npw*nspinor,blc_size)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)

!Local variables ------------------------------
!scalars
 integer :: nxknot,nyknot,nzknot,nk_tot,nkx,nky,nkz,kxord,kyord,kzord,natom
 integer :: ik_mesh,ik_bmsh,ik1,ik2,ik3,blc,size5,spinor,iatom
 integer :: istwf_k,npw_k,ii,mgfft
 real(dp) :: k3,k2,k1,start,step,shift
!arrays
 real(dp) :: kpoint(3)
 real(dp),allocatable :: ph1d(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_UNUSED(Pawang%l_max)
 ABI_UNUSED(Pawtab%basis_size)

 natom            = Cryst%natom
 mgfft            = MAXVAL(ngfft(1:3))

 Cp_bspl%natom    = natom
 Cp_bspl%nspinor  = nspinor
 Cp_bspl%blc_size = blc_size

 ABI_ALLOCATE(Cp_bspl%nlmn_sort,(natom))
 Cp_bspl%nlmn_sort = nlmn_sort

 nkx = kdiv(1) + 2 
 nky = kdiv(2) + 2
 nkz = kdiv(3) + 2
 nk_tot         = nkx*nky*nkz 

 Cp_bspl%nk_tot = nk_tot
 Cp_bspl%nkx    = nkx 
 Cp_bspl%nky    = nky 
 Cp_bspl%nkz    = nkz 

 ABI_ALLOCATE(Cp_bspl%xvec,(nkx))
 ABI_ALLOCATE(Cp_bspl%yvec,(nky))
 ABI_ALLOCATE(Cp_bspl%zvec,(nkz))
 ABI_ALLOCATE(Cp_bspl%kvec,(3,nk_tot))

 if (nkx<2.or.nky<2.or.nkz<2) then
   MSG_ERROR("At least two point are needed for B-spline")
 end if

 !Cp_bspl%xvec = arth(-half,one/(nkx-1),nkx) 
 !Cp_bspl%yvec = arth(-half,one/(nky-1),nky)
 !Cp_bspl%zvec = arth(-half,one/(nkz-1),nkz)
 
 step=one/(nkx-1); start = -half-step; step=(one+2*step)/(nkx-1)
 Cp_bspl%xvec = arth(start,step,nkx) 

 step=one/(nky-1); start = -half-step; step=(one+2*step)/(nky-1)
 Cp_bspl%yvec = arth(start,step,nky)

 step=one/(nkz-1); start = -half-step; step=(one+2*step)/(nkz-1)
 Cp_bspl%zvec = arth(start,step,nkz)

 ik_bmsh=0
 do ik3=1,nkz
   k3 = Cp_bspl%zvec(ik3) 
   do ik2=1,nky
     k2 = Cp_bspl%yvec(ik2) 
     do ik1=1,nkx
       k1 = Cp_bspl%xvec(ik1) 
       ik_bmsh=ik_bmsh+1
       Cp_bspl%kvec(:,ik_bmsh) = (/k1,k2,k3/)
       !write(std_out,*) (/k1,k2,k3/)
     end do
   end do
 end do

 kxord = kord(1); kxord=MIN(kxord,nkx)
 kyord = kord(2); kyord=MIN(kyord,nky)
 kzord = kord(3); kzord=MIN(kzord,nkz)

 Cp_bspl%kxord = kxord
 Cp_bspl%kyord = kyord
 Cp_bspl%kzord = kzord

 nxknot = nkx + kxord
 nyknot = nky + kyord
 nzknot = nkz + kzord

 Cp_bspl%nxknot = nxknot
 Cp_bspl%nyknot = nyknot
 Cp_bspl%nzknot = nzknot

 ABI_ALLOCATE(Cp_bspl%xknot,(nxknot))
 ABI_ALLOCATE(Cp_bspl%yknot,(nyknot))
 ABI_ALLOCATE(Cp_bspl%zknot,(nzknot))
 !
 ! Generate knots
 call dbsnak (nkx, Cp_bspl%xvec, kxord, Cp_bspl%xknot)
 call dbsnak (nky, Cp_bspl%yvec, kyord, Cp_bspl%yknot)
 call dbsnak (nkz, Cp_bspl%zvec, kzord, Cp_bspl%zknot)

 write(std_out,*)"xnot",Cp_bspl%xknot
 write(std_out,*)"ynot",Cp_bspl%yknot
 write(std_out,*)"znot",Cp_bspl%zknot

 ABI_ALLOCATE(Cp_bspl%Pblc_aki,(natom,nspinor,blc_size,nk_tot))

 do ik_mesh=1,nk_tot
   do blc=1,blc_size
     call cprj_alloc(Cp_bspl%Pblc_aki(:,:,blc,ik_mesh),0,nlmn_sort)
   end do
 end do

 ABI_ALLOCATE(Cp_bspl%Bscoef,(natom,nspinor,blc_size))

 do blc=1,Cp_bspl%blc_size
   do spinor=1,Cp_bspl%nspinor
     do iatom=1,Cp_bspl%natom
       size5 = nlmn_sort(iatom)
       ABI_ALLOCATE(Cp_bspl%Bscoef(iatom,spinor,blc)%value,(nkx,nky,nkz,2,size5))
       Cp_bspl%Bscoef(iatom,spinor,blc)%value=zero 
     end do
   end do
 end do

#if 1
 ABI_ALLOCATE(ph1d,(2,3*(2*mgfft+1)*Cryst%natom))
 call getph(Cryst%atindx,Cryst%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,Cryst%xred)
 !
 ! Evaluate projections on the B-spline k-mesh.
 !
 do ik_mesh=1,Cp_bspl%nk_tot
   istwf_k = 1
   !kg_k    => Wfd%Kdata(1)%kg_k
   npw_k   = npw

   do ii=1,3 ! Wrap in the first BZ thus enforcing traslational invariance.
     call wrap2_pmhalf(Cp_bspl%kvec(ii,ik_mesh),kpoint(ii),shift)  ! TODO overloaded interface.
   end do
   !
   ! FIXME complete the implementation of NC pseudos.
   call pbloch(Cryst,Psps,Cryst%natom,nspinor,npw_k,&
&    blc_size,blc_ug,kpoint,kg_k,ngfft,mgfft,ph1d,Cp_bspl%Pblc_aki(:,:,:,ik_mesh))
 end do 

 ABI_DEALLOCATE(ph1d)

 call cprj_bspline_mktensor(Cp_bspl) 
#endif

 DBG_EXIT("COLL")

end subroutine cprj_bspline_init
!!***

!----------------------------------------------------------------------

!!****f* m_cprj_bspline/cprj_bspline_mktensor
!! NAME
!!  cprj_bspline_mktensor
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_cprj_bspline
!!
!! CHILDREN
!!      initmpi_seq,initylmg,mkffnl,mkkpg,nonlop_ylm,ph1d3d
!!
!! SOURCE

subroutine cprj_bspline_mktensor(Cp_bspl) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_bspline_mktensor'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(cprj_bspl_t),intent(inout) :: Cp_bspl
!arrays

!Local variables ------------------------------
!scalars
 integer :: nxknot,nyknot,nzknot,ikx,iky,ikz,ik_tot
 integer :: nkx,nky,nkz,ldf,mdf,kxord,kyord,kzord
 integer :: spinor,iatom,ilmn,reim,blc
!arrays
 real(dp),pointer :: xvec(:),yvec(:),zvec(:),bscoef(:,:,:) 
 real(dp),allocatable :: xyzdata(:,:,:)
 real(dp),pointer :: xknot(:),yknot(:),zknot(:)

!************************************************************************

 DBG_ENTER("COLL")

 nkx = Cp_bspl%nkx 
 nky = Cp_bspl%nky 
 nkz = Cp_bspl%nkz 
 ldf = nkx
 mdf = nky

 kxord = Cp_bspl%kxord
 kyord = Cp_bspl%kyord 
 kzord = Cp_bspl%kzord 

 nxknot = Cp_bspl%nxknot
 nyknot = Cp_bspl%nyknot
 nzknot = Cp_bspl%nzknot

 xvec => Cp_bspl%xvec 
 yvec => Cp_bspl%yvec 
 zvec => Cp_bspl%zvec 

 xknot => Cp_bspl%xknot
 yknot => Cp_bspl%yknot
 zknot => Cp_bspl%zknot

 !write(std_out,*)nkx,nky,nkz
 !write(std_out,*)nxknot,nyknot,nzknot
 !write(std_out,*)kxord,kyord,kzord
 !write(std_out,*)"xvec",xvec
 !write(std_out,*)"yvec",yvec
 !write(std_out,*)"zvec",yvec
 !write(std_out,*)"xknot", xknot 
 !write(std_out,*)"yknot", yknot 
 !write(std_out,*)"zknot", zknot 

 ABI_ALLOCATE(xyzdata,(nkx,nky,nkz))

 do blc=1,Cp_bspl%blc_size
   do spinor=1,Cp_bspl%nspinor
     do iatom=1,Cp_bspl%natom
       do ilmn=1,Cp_bspl%nlmn_sort(iatom)
         do reim=1,2
           !
           ik_tot=0 ! Load data
           do ikz=1,nkz
             do iky=1,nky
               do ikx=1,nkx
                  ik_tot = ik_tot + 1
                  xyzdata(ikx,iky,ikz) = Cp_bspl%Pblc_aki(iatom,spinor,blc,ik_tot)%cp(reim,ilmn) 
               end do
             end do
           end do
           ! 
           ! Evaluate 3D tensor for B-spline.
           bscoef => Cp_bspl%bscoef(iatom,spinor,blc)%value(:,:,:,reim,ilmn)

           call dbs3in(nkx,xvec,nky,yvec,nkz,zvec,xyzdata,ldf,mdf,kxord,kyord,kzord,&
&            Cp_bspl%xknot,Cp_bspl%yknot,Cp_bspl%zknot,bscoef)

         end do
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(xyzdata)

 DBG_EXIT("COLL")

end subroutine cprj_bspline_mktensor
!!***

!----------------------------------------------------------------------

!!****f* m_cprj_bspline/cprj_bspline_free
!! NAME
!!  cprj_bspline_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      initmpi_seq,initylmg,mkffnl,mkkpg,nonlop_ylm,ph1d3d
!!
!! SOURCE

subroutine cprj_bspline_free(Cp_bspl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_bspline_free'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(cprj_bspl_t),intent(inout) :: Cp_bspl
!arrays

!Local variables ------------------------------
!scalars
 integer :: ik_mesh,blc,i3,i2,i1
!arrays

!************************************************************************

 if (associated(Cp_bspl%nlmn_sort))  then
   ABI_DEALLOCATE(Cp_bspl%nlmn_sort)
 end if

 if (associated(Cp_bspl%xvec))  then
   ABI_DEALLOCATE(Cp_bspl%xvec)
 end if
 if (associated(Cp_bspl%yvec))  then
   ABI_DEALLOCATE(Cp_bspl%yvec)
 end if
 if (associated(Cp_bspl%zvec))  then
   ABI_DEALLOCATE(Cp_bspl%zvec)
 end if
 if (associated(Cp_bspl%kvec))  then
   ABI_DEALLOCATE(Cp_bspl%kvec)
 end if

 if (associated(Cp_bspl%xknot))  then
   ABI_DEALLOCATE(Cp_bspl%xknot)
 end if
 if (associated(Cp_bspl%yknot))  then
   ABI_DEALLOCATE(Cp_bspl%yknot)
 end if
 if (associated(Cp_bspl%zknot))  then
   ABI_DEALLOCATE(Cp_bspl%zknot)
 end if

 if (associated(Cp_bspl%Bscoef)) then      ! TODO Write a module for ragged arrays. With free and nullify methods.
   do i3=1,SIZE(Cp_bspl%Bscoef,DIM=3)
     do i2=1,SIZE(Cp_bspl%Bscoef,DIM=2)
       do i1=1,SIZE(Cp_bspl%Bscoef,DIM=2)
         ABI_DEALLOCATE(Cp_bspl%Bscoef(i1,i2,i3)%value)
       end do
     end do
   end do
   ABI_DEALLOCATE(Cp_bspl%Bscoef)
 end if

 if (associated(Cp_bspl%Pblc_aki)) then
   do ik_mesh=1,SIZE(Cp_bspl%Pblc_aki,DIM=4)
     do blc=1,SIZE(Cp_bspl%Pblc_aki,DIM=3)
       call cprj_free(Cp_bspl%Pblc_aki(:,:,blc,ik_mesh))
     end do
   end do
   ABI_DEALLOCATE(Cp_bspl%Pblc_aki)
 end if

end subroutine cprj_bspline_free
!!***

!----------------------------------------------------------------------

!!****f* m_cprj_bspline/cprj_bspline_eval
!! NAME
!!  cprj_bspline_eval
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      initmpi_seq,initylmg,mkffnl,mkkpg,nonlop_ylm,ph1d3d
!!
!! SOURCE

subroutine cprj_bspline_eval(Cp_bspl,blc,k4intp,Cp_lk) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_bspline_eval'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blc
 real(dp),intent(in) :: k4intp(3)
 type(cprj_bspl_t),intent(in) :: Cp_bspl
 type(cprj_type),intent(inout) :: Cp_lk(Cp_bspl%natom,Cp_bspl%nspinor)

!Local variables ------------------------------
!scalars
 integer :: iatom,spinor,reim,ilmn
 integer :: nkx,nky,nkz,kxord,kyord,kzord
 real(dp) :: bval,x,y,z
!arrays
 real(dp),pointer :: xknot(:),yknot(:),zknot(:)
 real(dp),pointer :: xvec(:),yvec(:),zvec(:),bscoef(:,:,:) 

!************************************************************************

 DBG_ENTER("COLL")
 !
 ! B-spline interpolation.
 !
 nkx = Cp_bspl%nkx
 nky = Cp_bspl%nky
 nkz = Cp_bspl%nkz

 kxord = Cp_bspl%kxord
 kyord = Cp_bspl%kyord 
 kzord = Cp_bspl%kzord 

 xknot => Cp_bspl%xknot
 yknot => Cp_bspl%yknot
 zknot => Cp_bspl%zknot

 xvec => Cp_bspl%xvec 
 yvec => Cp_bspl%yvec 
 zvec => Cp_bspl%zvec 

 x = k4intp(1)
 y = k4intp(2)
 z = k4intp(3)

 do spinor=1,Cp_bspl%nspinor
   do iatom=1,Cp_bspl%natom
     do ilmn=1,Cp_bspl%nlmn_sort(iatom)
       do reim=1,2

           bscoef => Cp_bspl%Bscoef(iatom,spinor,blc)%value(:,:,:,reim,ilmn)

           bval = dbs3vl(x,y,z,kxord,kyord,kzord,xknot,yknot,zknot,nkx,nky,nkz,bscoef)
           Cp_lk(iatom,spinor)%cp(reim,ilmn) = bval

       end do
     end do
   end do
 end do

 DBG_EXIT("COLL")

end subroutine cprj_bspline_eval
!!***

!----------------------------------------------------------------------

!!****f* m_cprj_bspline/pbloch
!! NAME
!! pbloch
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Cryst<Crystal_structure>= data type gathering info on symmetries and unit cell
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang<pawang_type> angular mesh discretization and related data:
!! nspinor
!! npw
!! blc_size2
!! blc_ug(npw*nspinor,blc_size2)
!! kpoint(3)
!! kg_k(3,npw)
!! ngfft(18)=Information about 3D FFT.
!!
!! OUTPUT
!!  Pblc
!!
!! PARENTS
!!      m_cprj_bspline
!!
!! CHILDREN
!!      initmpi_seq,initylmg,mkffnl,mkkpg,nonlop_ylm,ph1d3d
!!
!! SOURCE

subroutine pbloch(Cryst,Psps,natom,nspinor,npw,blc_size2,blc_ug,kpoint,kg_k,ngfft,mgfft,ph1d,Pblc)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pbloch'
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,nspinor,npw,blc_size2,natom
 type(Crystal_structure),intent(in) :: Cryst
 !type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: ngfft(18),kg_k(3,npw)
 real(dp),intent(in) :: kpoint(3)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom) 
 complex(gwpc),intent(in) :: blc_ug(npw*nspinor,blc_size2)
 !type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Psps%usepaw)
 type(Cprj_type),intent(inout) :: Pblc(natom,nspinor,blc_size2)

!Local variables ------------------------------
!scalars
 integer,parameter :: idir0=0,ider0=0
 integer :: istat,optder,mkmem_,nkpg,dimffnl
 integer :: choice,cpopt,matblk,paw_opt,signs
 integer :: cplex_dij,lmnmax_ !,itypat !,ispden !,isp !,ilmn
 integer :: blc,istwf_k,useylm_
 integer :: iat,iatom,dimenl1,dimenl2,nnlout
 real(dp),parameter :: lambda0=zero
 real(dp) :: arg
 !character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: nloalg(5)
 integer,pointer :: indlmn_(:,:,:)
 real(dp) :: kptns_(3,1),ylmgr_dum(1),dum_enlout(0)
 real(dp),allocatable :: ylm_k(:,:),dum_ylm_gr_k(:,:,:),ffnl(:,:,:,:),kpg_k(:,:)
 real(dp),allocatable :: ph3d(:,:,:),vectin(:,:),enl(:,:,:),sij(:,:),phkxred(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 call initmpi_seq(MPI_enreg_seq) ! Fake MPI_type for sequential part.
 !
 ! The reason of this routine is that <p_i|ug> are needed for the B-spline of Vnl.
 ! and I didn't find any easy way to extract these quantities from nonlop_pl.
 ! Therefore I've decided to wrap nonlop_ylm so that B-spline can be used also
 ! for NC pseudos. Mind however that some dimensions in Psps cannot be used 
 ! since they are not correctly when useylm=0.

 ABI_CHECK(Psps%useylm==1,"useylm==0 not coded") 

 useylm_=1       ! FIXME Check Psps%lmnmax, I think it is wrong!
 lmnmax_ = Psps%lmnmax
 indlmn_ => Psps%indlmn

 !allocate(indlmn(6,psps%lmnmax))
 !call make_indlmn(ln_size,lmn_size,orbitals,indlmn)

 ! nloalg is also a special case
 nloalg(1)=4
 nloalg(2)=4
 nloalg(3)=199
 nloalg(4)=10
 nloalg(5)=1
 !nloalg(5)=psps%usepaw
 !
 ! Here I assume that the G-sphere is gamma-centered.
 ! Real Spherical Harmonics are always used to apply the non-local part even for NC pseudos.
 ! I did non find any easy way to extract only <p_nl|psi> from nonlop_pl.
 !
 ! ================================
 ! ==== Allocate NL strengths  ====
 ! ================================
 ! * Not used here, just avoiding problems with boundary checking.
 !
 if (Psps%usepaw==0) then ! Norm-conserving => constant Kleimann-Bylander energies function of (n,l).
   dimenl1=Psps%dimekb
   dimenl2=Cryst%ntypat
   ABI_ALLOCATE(enl,(Psps%dimekb,dimenl2,nspinor**2))
   ABI_ALLOCATE(sij,(0,0))
   enl(:,:,1)=Psps%ekb(:,:)    ! have to fill the value since paw_opt==0
   if (nspinor==2) then
     enl(:,:,2)=Psps%ekb(:,:)
     enl(:,:,3:4)=zero
   end if
                                                                                                   
 else ! PAW: store overlap coefficients and allocate memory for Dij coefficients (spin dependent) 
   !cplex_dij=Paw_ij(1)%cplex_dij
   cplex_dij=1
   dimenl1=Psps%dimekb*cplex_dij
   dimenl2=Cryst%natom
   ABI_ALLOCATE(enl,(dimenl1,dimenl2,nspinor**2))
   ABI_ALLOCATE(sij,(dimenl1,Cryst%ntypat))
 end if
 !
 ! ==============================================================
 ! ==== Tabulate Vnlk_ij on the homogeneous mesh of k-points ==== 
 ! ==============================================================
 !
 ABI_ALLOCATE(ylm_k,(npw,Psps%mpsang**2*useylm_))

 if (useylm_==1) then
   kptns_(:,1)=kpoint; optder=0; mkmem_=1
   ABI_ALLOCATE(dum_ylm_gr_k,(npw,3+6*(optder/2),Psps%mpsang**2))

   !  Here mband is not used if paral_compil_kpt=0
   call initylmg(Cryst%gprimd,kg_k,kptns_,mkmem_,MPI_enreg_seq,Psps%mpsang,npw,(/1/),1,&
&    (/npw/),1,optder,Cryst%rprimd,0,0,ylm_k,dum_ylm_gr_k)

   ABI_DEALLOCATE(dum_ylm_gr_k)
   istat = ABI_ALLOC_STAT
 end if
 !
 ! Compute (k+G) vectors (only if useylm_=1)
 nkpg=3*nloalg(5)  
 ABI_ALLOCATE(kpg_k,(npw,nkpg))
 if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw)
 !
 ! ========================================================
 ! ==== Compute nonlocal form factors ffnl at all (k+G) ====
 ! ========================================================
 !
 dimffnl=1+ider0 ! Derivatives are not needed. 
 ABI_ALLOCATE(ffnl,(npw,dimffnl,lmnmax_,Psps%ntypat))

 call mkffnl(Psps%dimekb,dimffnl,Psps%ekb,ffnl,Psps%ffspl,Cryst%gmet,Cryst%gprimd,ider0,idir0,indlmn_,&
&  kg_k,kpg_k,kpoint,lmnmax_,Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw,& 
&  Psps%ntypat,Psps%pspso,Psps%qgrid_ff,Cryst%rmet,Psps%usepaw,useylm_,ylm_k,ylmgr_dum)

 ABI_DEALLOCATE(ylm_k)
 istat = ABI_ALLOC_STAT

 !Allocate the arrays phkxred and ph3d, compute phkxred and eventually ph3d.
 ABI_ALLOCATE(phkxred,(2,Cryst%natom))
 do iat=1,Cryst%natom
   iatom=Cryst%atindx(iat)
   arg=two_pi*DOT_PRODUCT(kpoint,Cryst%xred(:,iat))
   phkxred(1,iatom)=DCOS(arg)
   phkxred(2,iatom)=DSIN(arg)
 end do
 !
 ! ============================================
 ! === Computation or <p_lmn|e^(ikr)U_i>  =====
 ! ============================================
 !
 matblk=nloalg(4); if (nloalg(1)>0) matblk=Cryst%natom
 istwf_k=1
 choice =0     ! only compute WF projected with NL projectors.
 cpopt  =0     ! <p_lmn|in> are computed and saved
 nnlout =0
 paw_opt=Psps%usepaw
 signs=1   ! NOT used. Applies the non-local operator to a function in reciprocal space
           ! In the case signs=1, the array vectout is not used, nor modified
           ! so that the same array as vectin can be used as a dummy argument;
           ! the same is true for the pairs npwin-npwout, ffnlin-ffnlout,
           ! kgin-kgout, ph3din-ph3dout, phkredin-phkxredout).

 ABI_ALLOCATE(ph3d,(2,npw,matblk))
 istat = ABI_ALLOC_STAT
 if (nloalg(1)>0) then ! Allocation as well as precomputation
   if (MPI_enreg_seq%mode_para/='b') then
     call ph1d3d(1,Cryst%natom,kg_k,matblk,Cryst%natom,npw,ngfft(1),ngfft(2),ngfft(3),phkxred,ph1d,ph3d)
   else 
     MSG_ERROR("Stop not coded")
   end if
 end if

 ABI_ALLOCATE(vectin,(2,npw*nspinor))

 do blc=1,blc_size2

   vectin(1,:) = DBLE (blc_ug(:,blc))
   vectin(2,:) = AIMAG(blc_ug(:,blc))

   call nonlop_ylm(Cryst%atindx1,choice,cpopt,Pblc(:,:,blc),dimenl1,dimenl2,dimffnl,dimffnl,&
&    enl,dum_enlout,ffnl,ffnl,Cryst%gprimd,idir0,indlmn_,istwf_k,&
&    kg_k,kg_k,kpg_k,kpg_k,kpoint,kpoint,lambda0,lmnmax_,matblk,mgfft,&
&    MPI_enreg_seq,Cryst%natom,Cryst%nattyp,ngfft,nkpg,nkpg,nloalg,nnlout,&
&    npw,npw,nspinor,nspinor,Cryst%ntypat,paw_opt,phkxred,phkxred,ph1d,&
&    ph3d,ph3d,signs,sij,vectin,Cryst%ucvol,vectin,vectin)
 end do

 ABI_DEALLOCATE(vectin)
 ABI_DEALLOCATE(ffnl)
 ABI_DEALLOCATE(kpg_k)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(phkxred)
 ABI_DEALLOCATE(ph3d)
 ABI_DEALLOCATE(enl)
 ABI_DEALLOCATE(sij)

 DBG_EXIT("COLL")

end subroutine pbloch       
!!***

!----------------------------------------------------------------------

END MODULE m_cprj_bspline
!!***
