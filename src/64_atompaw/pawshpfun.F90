!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawshpfun
!! NAME
!! pawshpfun
!!
!! FUNCTION
!! Compute shape function used in the definition
!! of compensation density (PAW)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ll= l quantum number
!!  mesh <type(pawrad_type)>=data containing radial grid information
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  norm= factor for shape function normalization
!!
!! SIDE effects
!!  shapefunc(mesh%mesh_size)=shape function g(r)
!!    In case of numerical shape function (shape_type=-1), shapefunc
!!    array contains the shape function read in psp file at input.
!!
!! NOTES
!!  Types of shape functions:
!!   type -1: numerical shape function, given in psp file
!!   type  1: g(r)=k(r).r^l; k(r)=exp(-(r/sigma)^lambda)
!!   type  2: g(r)=k(r).r^l; k(r)=[sin(Pi.r/rshp)/(Pi.r/rshp)]^2
!!   type  3: g(r)=alpha1.jl(q1.r)+alpha2.jl(q2.r)
!!
!! PARENTS
!!      pawdij0,pawinit,pawkij,psp7in
!!
!! CHILDREN
!!      jbessel,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawshpfun(ll,mesh,norm,pawtab,shapefunc)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_radmesh,          only : simp_gen, ifromr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawshpfun'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll
 real(dp),intent(out) :: norm
 type(pawrad_type),intent(in) :: mesh
 type(pawtab_type),intent(in) :: pawtab
!arrays
 real(dp),intent(inout) :: shapefunc(mesh%mesh_size)

!Local variables ------------------------------
!scalars
 integer :: ir,ishp
 real(dp) :: arg,argj1,argj2,besp,bespp,jbes1,jbes2,shapefunc1,shapefunc2
 real(dp) :: shapefunc3
!arrays
 real(dp) :: alpha(2),qq(2)
 real(dp),allocatable :: r2k(:)
!no_abirules
!Statement functions -----------------------------------
 shapefunc1(arg)= exp(-(arg/pawtab%shape_sigma)**pawtab%shape_lambda)
 shapefunc2(arg)= (sin(pi*arg/pawtab%rshp)/(pi*arg/pawtab%rshp))**2
 shapefunc3(argj1,argj2)= alpha(1)*argj1+alpha(2)*argj2

!***************************************************************************

 DBG_ENTER("COLL")

!Index for shape function cut-off radius
 ishp=ifromr(mesh,pawtab%rshp)-1

!Computation of non-normalized shape function
 if (pawtab%shape_type==-1) then
   shapefunc(1:ishp)=pawtab%shapefunc(1:ishp,1+ll)
 else if (pawtab%shape_type==1) then
   if (ll==0) then
     shapefunc(1)=one
     do ir=2,ishp
       shapefunc(ir)=shapefunc1(mesh%rad(ir))
     end do
   else
     shapefunc(1)=zero
     do ir=2,ishp
       shapefunc(ir)=shapefunc1(mesh%rad(ir))*mesh%rad(ir)**ll
     end do
   end if
 else if (pawtab%shape_type==2) then
   if (ll==0) then
     shapefunc(1)=one
     do ir=2,ishp
       shapefunc(ir)=shapefunc2(mesh%rad(ir))
     end do
   else
     shapefunc(1)=zero
     do ir=2,ishp
       shapefunc(ir)=shapefunc2(mesh%rad(ir))*mesh%rad(ir)**ll
     end do
   end if
 else if (pawtab%shape_type==3) then
   alpha(1:2)=pawtab%shape_alpha(1:2,1+ll)
   qq(1:2)=pawtab%shape_q(1:2,1+ll)
   do ir=1,ishp
     call jbessel(jbes1,besp,bespp,ll,0,qq(1)*mesh%rad(ir))
     call jbessel(jbes2,besp,bespp,ll,0,qq(2)*mesh%rad(ir))
     shapefunc(ir)=shapefunc3(jbes1,jbes2)
   end do
 end if

 if (ishp<mesh%mesh_size) shapefunc(ishp+1:mesh%mesh_size)=zero

!Shape function normalization
 if (pawtab%shape_type==-1.or.pawtab%shape_type==1.or.pawtab%shape_type==2) then
   ABI_ALLOCATE(r2k,(mesh%mesh_size))
   r2k=zero
   r2k(2:ishp)=shapefunc(2:ishp)*mesh%rad(2:ishp)**(2+ll)
   call simp_gen(norm,r2k,mesh);norm=one/norm
   shapefunc(1:ishp)=shapefunc(1:ishp)*norm
   if (pawtab%shape_type==-1) norm=one
   ABI_DEALLOCATE(r2k)
 else if (pawtab%shape_type==3) then
   norm=one
 end if

 DBG_EXIT("COLL")

end subroutine pawshpfun
!!***
