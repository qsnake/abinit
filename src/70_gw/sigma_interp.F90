!{\src2tex{textfont=tt}}
!!****f* ABINIT/interpolate_sigmak
!! NAME
!! interpolate_sigmak
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
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      dbs3in,dbsnak,set2unit,smpbz,sort_dp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine interpolate_sigmak(Cryst,Kmesh,kptrlatt,nshiftk,shiftk,lbd,ubd,isigk_ij,onkpt,okpt,osigk_ij,ierr)

 use m_profiling

 use defs_basis
 use m_bspline
 use m_errors

 use m_numeric_tools,  only : set2unit
 use m_crystal,        only : crystal_structure
 use m_bz_mesh,        only : bz_mesh_type, isequalk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'interpolate_sigmak'
 use interfaces_28_numeric_noabirule
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lbd,ubd,onkpt,nshiftk
 integer,intent(out) :: ierr
 type(Crystal_structure),intent(in) :: Cryst
 type(bz_mesh_type),intent(in) :: Kmesh
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: okpt(3,onkpt) 
 real(dp),intent(in) :: shiftk(3,nshiftk)
 complex(dpc),intent(in) :: isigk_ij(lbd:ubd,lbd:ubd,Kmesh%nibz)
 complex(dpc),intent(out) :: osigk_ij(lbd:ubd,lbd:ubd,onkpt)

!Local variables ------------------------------
!scalars
 integer :: kxord,kyord,kzord
 integer :: nxknot,nyknot,nzknot
 integer :: b1,b2,reim,ik_bz,ik_ibz,nkz,nky,nkx,ii,nkptlatt
 integer :: dir,mkpt,brav,nseen,nvec,nmiss,smp_nk
 integer :: ladd,radd,ikpt,ik1,ik2,ik3,bs_idx
 real(dp),parameter :: KPT_TOL=tol12
 logical :: found
 character(len=500) :: msg

!arrays
 integer,allocatable :: bsp2bz(:),iperm(:)
 real(dp) :: bs_kpt(3)
 real(dp),allocatable :: spkpt(:,:),seen(:),xvec(:),yvec(:),zvec(:),xknot(:),yknot(:),zknot(:)
 real(dp),allocatable :: xyzdata(:,:,:),bs_coef(:,:,:),work_bz(:,:,:,:),rout(:,:)
 complex(dpc),allocatable :: cwork_bz(:,:,:),umat_k(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_UNUSED(Cryst%nsym)  

 ierr=0
 if (Kmesh%nbz==1) then
   ierr=1
   MSG_WARNING("Cannot interpolate with a single k-point") 
   RETURN
 end if

 ABI_ALLOCATE(cwork_bz,(Kmesh%nbz,lbd:ubd,lbd:ubd))
 ABI_ALLOCATE(umat_k,(lbd:ubd,lbd:ubd))

 call set2unit(umat_k)
 !
 ! Reconstruct matrix elements in the full BZ, then rotate to the new basis.
 ! TODO: Here pay attention to time reversal symmetry.
 do ik_bz=1,Kmesh%nbz
   ik_ibz = Kmesh%tab(ik_bz)
   cwork_bz(ik_bz,:,:) =  isigk_ij(:,:,ik_ibz)
   cwork_bz(ik_bz,:,:) =  MATMUL( TRANSPOSE(CONJG(umat_k)), MATMUL(cwork_bz(ik_bz,:,:),umat_k) ) 
 end do

 ABI_DEALLOCATE(umat_k)

 ABI_ALLOCATE(work_bz,(2,Kmesh%nbz,lbd:ubd,lbd:ubd))
 work_bz(1,:,:,:) = DBLE(cwork_bz)
 work_bz(2,:,:,:) = AIMAG(cwork_bz)
 ABI_DEALLOCATE(cwork_bz)
 !
 ! Prepare B-spline interpolation.
 ! FIXME: Here I am assuming that kptrlatt is diagonal, only one shift is expected!
 !
 !kptrlatt = Kmesh%kptrlatt  ! FIXME these quantities have to be inititialized!
 !nshiftk  = Kmesh%nshift
 !shiftk   => Kmesh%shift

 ! Compute the number of k points in the G-space unit cell (will be multiplied by nshiftk later).
 nkptlatt = &
&    kptrlatt(1,1)*kptrlatt(2,2)*kptrlatt(3,3) &
&   +kptrlatt(1,2)*kptrlatt(2,3)*kptrlatt(3,1) &
&   +kptrlatt(1,3)*kptrlatt(2,1)*kptrlatt(3,2) &
&   -kptrlatt(1,2)*kptrlatt(2,1)*kptrlatt(3,3) &
&   -kptrlatt(1,3)*kptrlatt(2,2)*kptrlatt(3,1) &
&   -kptrlatt(1,1)*kptrlatt(2,3)*kptrlatt(3,2)

 mkpt=nkptlatt*nshiftk; brav=1 !brav=1 is able to treat all bravais lattices.

 ABI_ALLOCATE(spkpt,(3,mkpt))
 call smpbz(brav,std_out,kptrlatt,mkpt,smp_nk,nshiftk,0,shiftk,spkpt)

 ABI_ALLOCATE(seen,(smp_nk))
 ABI_ALLOCATE(iperm,(smp_nk))

 do dir=1,3
   nseen=1; seen(nseen) = spkpt(dir,1)
   do ik_bz=2,smp_nk
      if (ALL( ABS(spkpt(dir,ik_bz)-seen(1:nseen)) > KPT_TOL ) ) then
        nseen = nseen + 1
        seen(nseen) = spkpt(dir,ik_bz) 
      end if
   end do

   call sort_dp(nseen,seen,iperm, KPT_TOL)
   nvec = nseen
   radd = 0; ladd = 0
   if (ABS(seen(1)+one - seen(nseen)) > KPT_TOL) then 
     radd=1
     nvec = nvec + 1
   end if
   if (ABS(seen(nseen)-one - seen(1)) > KPT_TOL) then
     ladd=1
     nvec = nvec + 1
   end if

   if (dir==1) then
     ABI_ALLOCATE(xvec,(nvec))
     xvec(1+ladd:nvec-radd) = seen(1:nseen)
     if (ladd==1) xvec(ladd) = seen(nseen)-one
     if (radd==1) xvec(nvec) = seen(1)+one

   else if (dir==2) then
     ABI_ALLOCATE(yvec,(nvec))
     yvec(1+ladd:nvec-radd) = seen(1:nseen)
     if (ladd==1) yvec(ladd) = seen(nseen)-one
     if (radd==1) yvec(nvec) = seen(1)+one

   else if (dir==3) then
     ABI_ALLOCATE(zvec,(nvec))
     zvec(1+ladd:nvec-radd) = seen(1:nseen)
     if (ladd==1) zvec(ladd) = seen(nseen)-one
     if (radd==1) zvec(nvec) = seen(1)+one
   end if
 end do

 ABI_DEALLOCATE(seen)
 ABI_DEALLOCATE(iperm)

 ! TODO Add check on final partion. 
 ABI_DEALLOCATE(spkpt)
 ! 
 !
 ! Map B-spline mesh onto BZ.
 !
 nkx = SIZE(xvec)
 nky = SIZE(yvec)
 nkz = SIZE(zvec)

 write(std_out,*)" B-spline mesh: ",nkx,nky,nkz

 ABI_ALLOCATE(bsp2bz,(nkx*nky*nkz))

 bs_idx=0; nmiss=0
 do ik3=1,nkz
   do ik2=1,nky
     do ik1=1,nkx
       bs_idx = bs_idx+1
       bs_kpt = (/xvec(ik1), yvec(ik2), zvec(ik3)/)
       
       ik_bz=0; found=.FALSE.
       do while (ik_bz<Kmesh%nbz .and. .not.found)
         ik_bz = ik_bz + 1
         found = isequalk(bs_kpt,Kmesh%bz(:,ik_bz))
       end do

       if (found) then
         bsp2bz(bs_idx) = ik_bz
       else 
         nmiss = nmiss+1
       end if

     end do
   end do
 end do

 if (nmiss>0) then
   ierr=ierr+1
   write(msg,'(a,i0,a)')"Found ",nmiss," k-points in B-spline mesh that does not belong to the BZ"
   MSG_ERROR(msg)
 end if
 !
 ! Generate knots (Order should be selected in a more careful way)
 kxord = nkx
 kyord = nky
 kzord = nkz

 nxknot = nkx + kxord
 nyknot = nky + kyord
 nzknot = nkz + kzord

 ABI_ALLOCATE(xknot,(nxknot))
 ABI_ALLOCATE(yknot,(nyknot))
 ABI_ALLOCATE(zknot,(nzknot))

 call dbsnak (nkx, xvec, kxord, xknot)
 call dbsnak (nky, yvec, kyord, yknot)
 call dbsnak (nkz, zvec, kzord, zknot)
 !
 ABI_ALLOCATE(xyzdata,(nkx,nky,nkz))
 ABI_ALLOCATE(bs_coef,(nkx,nky,nkz))
 ABI_ALLOCATE(rout,(2,onkpt))

 do b2=lbd,ubd
   do b1=lbd,ubd
     do reim=1,2
       !
       ! Load data.
       bs_idx = 0  
       do ik3=1,nkz
         do ik2=1,nky
           do ik1=1,nkx
             bs_idx = bs_idx+1
             ik_bz = bsp2bz(bs_idx)
             xyzdata(ik1,ik2,ik3) = work_bz(reim,ik_bz,b1,b2)
           end do
         end do
       end do
       !
       ! Construct 3D tensor for B-spline.
       call dbs3in(nkx,xvec,nky,yvec,nkz,zvec,xyzdata,nkx,nky,kxord,kyord,kzord,xknot,yknot,zknot,bs_coef)
       !
       do ikpt=1,onkpt ! B-spline interpolation.
         rout(reim,ikpt) = dbs3vl(okpt(1,ikpt),okpt(2,ikpt),okpt(3,ikpt),kxord,kyord,kzord,xknot,yknot,zknot,nkx,nky,nkz,bs_coef)
       end do
     end do ! reim

     osigk_ij(b1,b2,:) = DCMPLX(rout(1,:),rout(2,:))

   end do
 end do

 do ikpt=1,onkpt
   write(77,'(a,3es16.8,a)')"# kpt= (",okpt(:,ikpt),") "
   write(77,'(1x,(10f9.5))')(REAL(osigk_ij(ii,ii,ikpt))*Ha_eV,ii=lbd,ubd)
 end do

 ABI_DEALLOCATE(rout)
 ABI_DEALLOCATE(xyzdata)
 ABI_DEALLOCATE(bs_coef)

 ABI_DEALLOCATE(xknot)
 ABI_DEALLOCATE(yknot)
 ABI_DEALLOCATE(zknot)
 ABI_DEALLOCATE(xvec)
 ABI_DEALLOCATE(yvec)
 ABI_DEALLOCATE(zvec)

 ABI_DEALLOCATE(bsp2bz)
 ABI_DEALLOCATE(work_bz)

 DBG_EXIT("COLL")

end subroutine interpolate_sigmak
!!***


