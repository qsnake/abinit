!{\src2tex{textfont=tt}}
!!****f* ABINIT/opernld_ylm
!! NAME
!! opernld_ylm
!!
!! FUNCTION
!! * Operate with the non-local part of the hamiltonian,
!!   in order to get contributions to energy/forces/stress/dyn.matrix/elst tens.
!!   from projected scalars
!! * Operate with the non-local projectors and the overlap matrix Sij
!!   in order to get contributions to <c|S|c>
!!   from projected scalars
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  choice=chooses possible output
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  cplex_fac=1 if gxfac scalars are real, 2 if gxfac scalars are complex
!!  d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)=2nd gradients of projected scalars
!!  dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)=gradients of projected scalars
!!  dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)=gradients of reduced projected scalars
!!  gx(cplex,nlmn,nincat,nspinor)= projected scalars
!!  gxfac(cplex_fac,nlmn,nincat,nspinor)= reduced projected scalars related to Vnl (NL operator)
!!  gxfac_sij(cplex,nlmn,nincat,nspinor)= reduced projected scalars related to Sij (overlap)
!!  ia3=gives the absolute number of the first atom in the subset presently treated
!!  natom=number of atoms in cell
!!  nd2gxdt=second dimension of d2gxdt
!!  ndgxdt=second dimension of dgxdt
!!  ndgxdtfac=second dimension of dgxdtfac
!!  nincat=number of atoms in the subset here treated
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nnlout=dimension of enlout
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  ucvol=unit cell volume (bohr^3)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! --If (paw_opt==0, 1 or 2)
!!    enlout(nnlout)= contribution to the non-local part of the following properties:
!!      if choice=1 : enlout(1)               -> the energy
!!      if choice=2 : enlout(1:3*natom)       -> the forces
!!      if choice=3 : enlout(1:6)             -> the stresses
!!      if choice=23: enlout(1:6+3*natom)     -> the forces and the stresses
!!      if choice=4 : enlout(1:6*natom)       -> the frozen wf part of dyn. mat.
!!      if choice=4 : enlout(3)       -> the frozen wf part of dyn. mat.
!!      if choice=24: enlout(1:9*natom)       -> the forces and the frozen wf part of dyn. mat.
!!      if choice=6 : enlout(1:6*(3*natom+6)) -> the frozen wf part of elastic tensor
!!    if (choice==3.or.choice==6.or.choice==23)
!!      enlk=contribution to the non-local part of the energy
!!    if (choice==5)
!!      nlout(nnlout)==contribution to the derivatives of nl part of the energy wrt to k
!!    if (choice==6)
!!      fnlk(3*natom)=contribution to the non-local part of the forces
!!      strnlk(6)=contribution to the non-local part of the stresses
!! --If (paw_opt==3)
!!    if choice=1 : enlout(nnlout)= contribution to <c|S|c>  (nnlout=1)
!! --If (paw_opt==4)
!!    not available
!!
!! NOTES
!! Operate for one type of atom, and within this given type of atom,
!! for a subset of at most nincat atoms.
!!
!! PARENTS
!!      nonlop_ylm
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine opernld_ylm(choice,cplex,cplex_fac,d2gxdt,dgxdt,dgxdtfac,enlk,enlout,fnlk,&
&                      gx,gxfac,gxfac_sij,ia3,natom,nd2gxdt,ndgxdt,ndgxdtfac,nincat,&
&                      nlmn,nnlout,nspinor,paw_opt,strnlk,ucvol)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'opernld_ylm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,cplex_fac,ia3,natom,nd2gxdt,ndgxdt
 integer,intent(in) :: ndgxdtfac,nincat,nlmn,nnlout,nspinor,paw_opt
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: enlk
!arrays
 real(dp),intent(in) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: gx(cplex,nlmn,nincat,nspinor),gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(inout) :: enlout(nnlout),fnlk(3*natom),strnlk(6)

!Local variables-------------------------------
!scalars
 integer :: ia,iashift,ilmn,iplex,ishift,ispinor,mu,mua,mub,mushift,mut,muu,nu,nushift
 real(dp) :: factvol
 complex(dpc) :: cft, cfu
!arrays
 integer :: twist_dir(6)=(/2,3,3,1,1,2/)
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 real(dp) :: enlj(6),gxfacj(cplex)

! *************************************************************************

 if (cplex_fac<cplex) stop "opernld_ylm: BUG - invalid cplex_fac (<cplex) !"

 if (paw_opt==0.or.paw_opt==1.or.paw_opt==2) then

!  ============== Accumulate the non-local energy ===============
   if (choice==1) then
     do ispinor=1,nspinor
       do ia=1,nincat
         do ilmn=1,nlmn
           do iplex=1,cplex
             enlout(1)=enlout(1)+gxfac(iplex,ilmn,ia,ispinor)*gx(iplex,ilmn,ia,ispinor)
           end do
         end do
       end do
     end do
   end if

!  ============ Accumulate the forces contributions =============
   if (choice==2.or.choice==23.or.choice==24) then
     ishift=0;if (choice==23) ishift=6
     factvol=two*ucvol
     do ispinor=1,nspinor
       do ia=1,nincat
         enlj(1:3)=zero
         iashift=3*(ia+ia3-2)+ishift
         do ilmn=1,nlmn
           do mu=1,3
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*dgxdt(iplex,mu+ishift,ilmn,ia,ispinor)
             end do
           end do
         end do
         enlout(iashift+1:iashift+3)=enlout(iashift+1:iashift+3)+factvol*enlj(1:3)
       end do
     end do
   end if

!  ======== Accumulate the stress tensor contributions ==========
   if (choice==3.or.choice==23) then
     enlj(1:6)=zero
     do ispinor=1,nspinor
       do ia=1,nincat
         do ilmn=1,nlmn
           gxfacj(1:cplex)=gxfac(1:cplex,ilmn,ia,ispinor)
           do iplex=1,cplex
             enlk=enlk+gxfacj(iplex)*gx(iplex,ilmn,ia,ispinor)
           end do
           do mu=1,6
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfacj(iplex)*dgxdt(iplex,mu,ilmn,ia,ispinor)
             end do
           end do
         end do
       end do
     end do
     enlout(1:6)=enlout(1:6)+two*enlj(1:6)
   end if

!  ====== Accumulate the dynamical matrix contributions =========
   if (choice==4.or.choice==24) then
     ishift=0;if (choice==24) ishift=3*natom
     factvol=two*ucvol
     do ispinor=1,nspinor
       do ia=1,nincat
         enlj(1:6)=zero
         iashift=6*(ia+ia3-2)+ishift
         do ilmn=1,nlmn
           do mu=1,6
             mua=alpha(mu);mub=beta(mu)
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*d2gxdt(iplex,mu,ilmn,ia,ispinor)&
&               +dgxdtfac(iplex,mub,ilmn,ia,ispinor)*dgxdt(iplex,mua,ilmn,ia,ispinor)
             end do
           end do
         end do
         enlout(iashift+1:iashift+6)=enlout(iashift+1:iashift+6)+factvol*enlj(1:6)
       end do
     end do
   end if

!  ======== Accumulate the contributions derivatives of E wrt to k ==========
   if (choice==5) then
     enlj(1:3)=zero
     do ispinor=1,nspinor
       do ia=1,nincat
         do ilmn=1,nlmn
           do mu=1,3
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*dgxdt(iplex,mu,ilmn,ia,ispinor)
             end do
           end do
         end do
       end do
     end do
     enlout(1:3)=enlout(1:3)+two*enlj(1:3)
   end if

!  ======== Accumulate the contributions of twist derivatives of E wrt to k ==========
   if (choice==53) then
     enlj(1:3)=zero
     do ispinor=1,nspinor
       do ia=1,nincat
         do ilmn=1,nlmn
           do mu=1,3
             mut = twist_dir(2*mu-1)
             muu = twist_dir(2*mu)

             cft = cmplx(dgxdt(1,mut,ilmn,ia,ispinor),dgxdt(2,mut,ilmn,ia,ispinor))
             cfu = cmplx(dgxdtfac(1,muu,ilmn,ia,ispinor),dgxdtfac(2,muu,ilmn,ia,ispinor))
             enlj(mu) = enlj(mu) + aimag(conjg(cft)*cfu)

             cfu = cmplx(dgxdt(1,muu,ilmn,ia,ispinor),dgxdt(2,muu,ilmn,ia,ispinor))
             cft = cmplx(dgxdtfac(1,mut,ilmn,ia,ispinor),dgxdtfac(2,mut,ilmn,ia,ispinor))
             enlj(mu) = enlj(mu) - aimag(conjg(cfu)*cft)

           end do
         end do
       end do
     end do
     enlout(1:3)=enlout(1:3)+enlj(1:3)
   end if


!  ======= Accumulate the elastic tensor contributions ==========
   if (choice==6) then
     do ispinor=1,nspinor
       do ia=1,nincat
         iashift=3*(ia+ia3-2)
         do ilmn=1,nlmn
           do iplex=1,cplex
             enlk=enlk+gxfac(iplex,ilmn,ia,ispinor)*gx(iplex,ilmn,ia,ispinor)
           end do
           enlj(1:3)=zero
           do mu=1,3
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*dgxdt(iplex,6+mu,ilmn,ia,ispinor)
             end do
           end do
           fnlk(iashift+1:iashift+3)=fnlk(iashift+1:iashift+3)+two*enlj(1:3)
           enlj(1:6)=zero
           do mu=1,6
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac(iplex,ilmn,ia,ispinor)*dgxdt(iplex,mu,ilmn,ia,ispinor)
             end do
           end do
           strnlk(1:6)=strnlk(1:6)+two*enlj(1:6)
           do mub=1,6
             mushift=6*(mub-1);nushift=(3*natom+6)*(mub-1)
             do mua=1,6
               mu=mushift+mua;nu=nushift+mua
               do iplex=1,cplex
                 enlout(nu)=enlout(nu)+two* &
&                 (gxfac(iplex,ilmn,ia,ispinor)*d2gxdt(iplex,mu,ilmn,ia,ispinor)&
&                 +dgxdtfac(iplex,mua,ilmn,ia,ispinor)*dgxdt(iplex,mub,ilmn,ia,ispinor))
               end do
             end do
             mushift=36+3*(mub-1);nushift=6+iashift+(3*natom+6)*(mub-1)
             do mua=1,3
               mu=mushift+mua;nu=nushift+mua
               do iplex=1,cplex
                 enlout(nu)=enlout(nu)+two* &
&                 (gxfac(iplex,ilmn,ia,ispinor)*d2gxdt(iplex,mu,ilmn,ia,ispinor)&
&                 +dgxdtfac(iplex,mub,ilmn,ia,ispinor)*dgxdt(iplex,6+mua,ilmn,ia,ispinor))
               end do
             end do
           end do
         end do
       end do
     end do
   end if
 end if

 if (paw_opt==3) then

!  ============== Accumulate contribution to <c|S|c> ===============
   if (choice==1) then
     do ispinor=1,nspinor
       do ia=1,nincat
         do ilmn=1,nlmn
           do iplex=1,cplex
             enlout(1)=enlout(1)+gxfac_sij(iplex,ilmn,ia,ispinor)*gx(iplex,ilmn,ia,ispinor)
           end do
         end do
       end do
     end do
   end if

!  ============== Accumulate contribution to <c|dS/atm_pos|c> ===============
   if (choice==2.or.choice==23.or.choice==24) then
     ishift=0;if (choice==23) ishift=6
     factvol=two*ucvol
     do ispinor=1,nspinor
       do ia=1,nincat
         enlj(1:3)=zero
         iashift=3*(ia+ia3-2)+ishift
         do ilmn=1,nlmn
           do mu=1,3
             do iplex=1,cplex
               enlj(mu)=enlj(mu)+gxfac_sij(iplex,ilmn,ia,ispinor)*dgxdt(iplex,mu+ishift,ilmn,ia,ispinor)
             end do
           end do
         end do
         enlout(iashift+1:iashift+3)=enlout(iashift+1:iashift+3)+factvol*enlj(1:3)
       end do
     end do
   end if

 end if

end subroutine opernld_ylm
!!***
