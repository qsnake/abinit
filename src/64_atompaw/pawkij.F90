!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawkij
!! NAME
!! pawkij
!!
!! FUNCTION
!! PAW: deduce kinetic part of psp strength (Dij) from the knowledge of frozen Dij (Dij0)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT,GJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  indlmn(6,lmnmax)= array giving l,m,n,lm,ln,s for i=lmn
!!  lmnmax=max number of (l,m,n) components over all type of psps
!!  ncore(radmesh_core%mesh_size)=atomic core density
!!  opt_init=flag defining the storage of PAW atomic data
!!           0: PAW atomic data have not been initialized (in pawtab)
!!           1: PAW atomic data have been initialized (in pawtab)
!!  opt_vhnzc=flag defining the inclusion of VH(nZc) in computation
!!            0: VH(nZc) is not taken into account
!!            1: VH(nZc) is taken into account
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  radmesh <type(pawrad_type)>=paw radial mesh (and related data)
!!  radmesh_core <type(pawrad_type)>=radial mesh (and related data) for the core densities
!!  radmesh_vloc <type(pawrad_type)>=radial mesh (and related data) for the local potential (VH(tnZc))
!!  vhtnzc(radmesh_vloc%mesh_size)= local potential VH(tnZc)
!!  znucl= valence and total charge of the atomic species
!!
!! OUTPUT
!!  kij(pawtab%lmn2_size)= kinetic part of Dij
!!
!! PARENTS
!!      psp7in
!!
!! CHILDREN
!!      bound_deriv,pawshpfun,pawvhnzc,simp_gen,spline,splint
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawkij(indlmn,kij,lmnmax,ncore,opt_init,opt_vhnzc,pawtab,radmesh,radmesh_core,radmesh_vloc,vhtnzc,znucl)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_splines

 use m_radmesh,   only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawkij'
 use interfaces_32_util
 use interfaces_64_atompaw, except_this_one => pawkij
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lmnmax,opt_init,opt_vhnzc
 real(dp),intent(in) :: znucl
 type(pawrad_type),intent(in) :: radmesh,radmesh_core,radmesh_vloc
 type(pawtab_type),intent(in) :: pawtab
!arrays
 integer,intent(in) :: indlmn(6,lmnmax)
 real(dp),intent(out) :: kij(pawtab%lmn2_size)
 real(dp),intent(in) :: ncore(radmesh_core%mesh_size)
 real(dp),intent(in) :: vhtnzc(radmesh_vloc%mesh_size)

!Local variables ---------------------------------------
 integer :: il,ilm,iln,ilmn,j0lmn,jl,jlm,jln,jlmn,klmn,lmn2_size,meshsz
 real(dp) :: intg,intvh,yp1,ypn
 real(dp),allocatable :: ff(:),shpf(:),vhnzc(:),vhtnzc_sph(:),work1(:),work2(:)

! *********************************************************************

 DBG_ENTER("COLL")

 lmn2_size=pawtab%lmn2_size
 meshsz=radmesh%mesh_size
 ABI_ALLOCATE(ff,(meshsz))

!Retrieve VH(tnZc) on the correct radial mesh
 ABI_ALLOCATE(vhtnzc_sph,(meshsz))
 if ((radmesh%mesh_type/=radmesh_vloc%mesh_type).or.&
& (radmesh%rstep    /=radmesh_vloc%rstep)    .or.&
& (radmesh%lstep    /=radmesh_vloc%lstep)) then
   call bound_deriv(vhtnzc(1:radmesh_vloc%mesh_size),radmesh_vloc,radmesh_vloc%mesh_size,yp1,ypn)
   ABI_ALLOCATE(work1,(radmesh_vloc%mesh_size))
   ABI_ALLOCATE(work2,(radmesh_vloc%mesh_size))
   call spline(radmesh_vloc%rad,vhtnzc,radmesh_vloc%mesh_size,yp1,ypn,work1)
   call splint(radmesh_vloc%mesh_size,radmesh_vloc%rad,vhtnzc,work1,meshsz,radmesh%rad(1:meshsz),vhtnzc_sph)
   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
 else
   vhtnzc_sph(1:meshsz)=vhtnzc(1:meshsz)
 end if

!Initiialize Kij with Dij0
!=========================
 kij(1:lmn2_size)=pawtab%dij0(1:lmn2_size)

!Substraction of <phi_i|vh(nZc)|phi_j> on the PAW sphere
!=======================================================
 if (opt_vhnzc/=0) then
   ABI_ALLOCATE(vhnzc,(radmesh_core%mesh_size))
   call pawvhnzc(ncore,radmesh_core,vhnzc,znucl)
   do jlmn=1,pawtab%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
       if (jlm==ilm) then
         ff(1:meshsz)=pawtab%phi(1:meshsz,iln)*pawtab%phi(1:meshsz,jln)*vhnzc(1:meshsz)
         call simp_gen(intg,ff,radmesh)
         kij(klmn)=kij(klmn)-intg
       end if
     end do
   end do
   ABI_DEALLOCATE(vhnzc)
 end if

!Substraction of -<tphi_i|vh(tnZc)|tphi_j> on the PAW sphere
!===========================================================
 do jlmn=1,pawtab%lmn_size
   j0lmn=jlmn*(jlmn-1)/2
   jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
   do ilmn=1,jlmn
     klmn=j0lmn+ilmn
     ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
     if (jlm==ilm) then
       ff(1:meshsz)=pawtab%tphi(1:meshsz,iln)*pawtab%tphi(1:meshsz,jln)*vhtnzc_sph(1:meshsz)
       call simp_gen(intg,ff,radmesh)
       kij(klmn)=kij(klmn)+intg
     end if
   end do
 end do

!Substraction of -int[vh(tnzc)*Qijhat(r)dr]
!==========================================
 if (opt_init==0) then
   ABI_ALLOCATE(shpf,(meshsz))
   call pawshpfun(0,radmesh,intg,pawtab,shpf)
   ff(1:meshsz)=vhtnzc_sph(1:meshsz)*shpf(1:meshsz)*radmesh%rad(1:meshsz)**2
   ABI_DEALLOCATE(shpf)
   call simp_gen(intvh,ff,radmesh)
   do jlmn=1,pawtab%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     jl=indlmn(1,jlmn);jln=indlmn(5,jlmn);jlm=indlmn(4,jlmn)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       il=indlmn(1,ilmn);iln=indlmn(5,ilmn);ilm=indlmn(4,ilmn)
       if (ilm==jlm) then
         ff(1:meshsz)=(pawtab%phi (1:meshsz,iln)*pawtab%phi (1:meshsz,jln)&
&         -pawtab%tphi(1:meshsz,iln)*pawtab%tphi(1:meshsz,jln))
         call simp_gen(intg,ff,radmesh)
         kij(klmn)=kij(klmn)+intvh*intg
       end if
     end do
   end do
 else
   ff(1:meshsz)=vhtnzc_sph(1:meshsz)*pawtab%shapefunc(1:meshsz,1)*radmesh%rad(1:meshsz)**2
   call simp_gen(intvh,ff,radmesh)
   do jlmn=1,pawtab%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     jl=indlmn(1,jlmn);jln=indlmn(5,jlmn);jlm=indlmn(4,jlmn)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       il=indlmn(1,ilmn);iln=indlmn(5,ilmn);ilm=indlmn(4,ilmn)
       if (ilm==jlm) then
         intg=pawtab%qijl(1,klmn)*sqrt(four_pi)
         kij(klmn)=kij(klmn)+intvh*intg
       end if
     end do
   end do
 end if

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(vhtnzc_sph)

 DBG_EXIT("COLL")

 end subroutine pawkij
!!***
