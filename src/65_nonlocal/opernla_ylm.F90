!{\src2tex{textfont=tt}}
!!****f* ABINIT/opernla_ylm
!! NAME
!! opernla_ylm
!!
!! FUNCTION
!! For a given wave-function |c>, get all projected scalars
!! <p_lmn|c> where |p_lmn> are non-local projectors
!!   With:
!!   <p_lmn|c>=4pi/sqrt(vol) (i)^l Sum_g[c(g).f_nl(g).Y_lm(g).exp(2pi.i.g.R)]
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  choice=chooses possible output (see below):
!!         if choice>=0: compute projected scalars
!!         if choice<0: use already computed projected scalars
!!         if ABS(choice)>1, then compute additional quantities:
!!         2 : compute projected scalars and derivatives wrt atm pos.
!!         3 : compute projected scalars and derivatives wrt strains
!!         23: compute projected scalars, derivatives wrt atm pos. and derivatives wrt strains
!!         4 : compute projected scalars, derivatives wrt atm pos. and 2nd derivatives wrt atm pos.
!!         24: compute projected scalars, derivatives wrt atm pos. and 2nd derivatives wrt atm pos.
!!         5, 51, 52 : compute projected scalars and derivatives wrt wave vector k
!!         53: computed projected scalars and derivatives wrt wave vector k in direction idir+1
!!                 and idir-1
!!         6 : compute projected scalars, derivatives wrt atm pos., 2nd derivatives wrt 2 atm pos.,
!!           and 2nd derivatives wrt atm. and strains
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  dimffnl=second dimension of ffnl
!!  ffnl(npw,dimffnl,nlmn)= nonlocal quantities containing nonlocal form factors
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,signs=2)
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  kpg(npw,nkpg)=(k+G) components          for ikpg=1...3   (if nkpg=3 or 9)
!!       [(k+G)_a].[(k+G)_b] quantities for ikpg=4...9   (if nkpg=9)
!!  matblk=dimension of the array ph3d
!!  mpi_enreg=informations about MPI parallelization
!!  ndgxdt=second dimension of dgxdt
!!  nd2gxdt=second dimension of d2gxdt
!!  nincat=number of atoms in the subset here treated
!!  nkpg=second dimension of array kpg (0, 3 or 9)
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npw=number of plane waves in reciprocal space
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!  signs=chooses possible output (see below)
!!  ucvol=unit cell volume (bohr^3)
!!  vect(2,npw*my_nspinor)=starting vector in reciprocal space
!!
!! OUTPUT
!!  if (choice>1) dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)=
!!     gradients of projected scalars wrt coords  (choice=2, 4, 6 or 23)
!!                                    wrt strains (choice=3 or 23)
!!                                    wrt k       (choice=5, 51, 52, 53)
!!  if (choice=4, 24 or 6) d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)=
!!     2nd grads of projected scalars wrt 2 coords         (choice=4 or 24)
!!                                    wrt coords & strains (choice=6)
!!
!! SIDE EFFECTS
!!  gx(cplex,nlmn,nincat,nspinor)= projected scalars - input if choice<0, output if choice>=0
!!
!! NOTES
!! Operate for one type of atom, and within this given type of atom,
!! for a subset of at most nincat atoms.
!!
!! PARENTS
!!      getcprj,nonlop_ylm
!!
!! CHILDREN
!!      timab,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine opernla_ylm(choice,cplex,dimffnl,d2gxdt,dgxdt,ffnl,gx,ia3,idir,&
&       indlmn,istwf_k,kpg,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpg,nlmn,&
&       nloalg,npw,nspinor,ph3d,signs,ucvol,vect)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'opernla_ylm'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,dimffnl,ia3,idir,istwf_k,matblk,nd2gxdt
 integer,intent(in) :: ndgxdt,nincat,nkpg,nlmn,npw,nspinor,signs
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: indlmn(6,nlmn),nloalg(5)
 real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn),kpg(npw,nkpg),ph3d(2,npw,matblk)
 real(dp),intent(in) :: vect(2,npw*nspinor)
 real(dp),intent(out) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)
 real(dp),intent(out) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor),gx(cplex,nlmn,nincat,nspinor)

!Local variables-------------------------------
!no_abirules
 integer :: choice_,ffnl_dir,gama,gamb,gamc,gamd,ia,iaph3d,ierr,il,ilmn
 integer :: ipw,ipw0,ipwshft,ishift,ispinor,jpw,mu
 integer :: mua,mub,nua1,nua2,nub1,nub2,old_paral_level,spaceComm
 real(dp), parameter :: two_pi2=two_pi*two_pi
 real(dp) :: aux_i,aux_i1,aux_i2,aux_i3,aux_i4
 real(dp) :: aux_r,aux_r1,aux_r2,aux_r3,aux_r4
 real(dp) :: buffer_i,buffer_i1,buffer_i2,buffer_i3,buffer_i4,buffer_i5,buffer_i6
 real(dp) :: buffer_ia,buffer_ib,buffer_ic,buffer_id,buffer_ie,buffer_if
 real(dp) :: buffer_r,buffer_r1,buffer_r2,buffer_r3,buffer_r4,buffer_r5,buffer_r6
 real(dp) :: buffer_ra,buffer_rb,buffer_rc,buffer_rd,buffer_re,buffer_rf
 real(dp) :: kpga,kpgb,kpgc,kpgd,scale,scale2,scale3,scale4,scale5,wt
 logical :: parity
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 integer :: ffnl_dir_dat(6)=(/3,4,4,2,2,3/)
 real(dp) :: tsec(2)
 real(dp),allocatable :: scali(:),scalr(:),scalari(:,:),scalarr(:,:)

! *************************************************************************

 if (choice==-1) return

!Useful variables
 choice_=abs(choice)
 wt=four_pi/sqrt(ucvol);if (cplex==1) wt=2.d0*wt
 ipw0=1;if (istwf_k==2.and.mpi_enreg%me_g0==1) ipw0=2

!Allocate work space
 ABI_ALLOCATE(scali,(npw))
 ABI_ALLOCATE(scalr,(npw))

!Loop on spinorial components
 do ispinor =1,nspinor

   ipwshft=(ispinor-1)*npw

!  Loops (blocking)
!  $OMP PARALLEL DEFAULT(PRIVATE) &
!  $OMP&SHARED(cplex,nincat,nloalg,npw,nlmn,indlmn,ia3)
!  $OMP&SHARED(idir,ipw0,ipwshft,choice_,signs,ndgxdt,nd2gxdt,alpha,beta)
!  $OMP&SHARED(vect,ph3d,ffnl,scalr,scali)
!  $OMP&SHARED(gx,dgxdt,d2gxdt,kpg)
!  $OMP DO

!  Loop on atoms
   do ia=1,nincat
     iaph3d=ia;if (nloalg(1)>0) iaph3d=ia+ia3-1

!    Compute Sum_g[c(g).exp(2pi.i.g.R)]
     do ipw=ipw0,npw
       jpw=ipw+ipwshft
       scalr(ipw)=(vect(1,jpw)*ph3d(1,ipw,iaph3d)-vect(2,jpw)*ph3d(2,ipw,iaph3d))
       scali(ipw)=(vect(2,jpw)*ph3d(1,ipw,iaph3d)+vect(1,jpw)*ph3d(2,ipw,iaph3d))
     end do
     if (ipw0==2) then
       scalr(1)=half*vect(1,1+ipwshft)*ph3d(1,1,iaph3d)
       scali(1)=half*vect(1,1+ipwshft)*ph3d(2,1,iaph3d)
     end if

!    --------------------------------------------------------------------
!    ALL CHOICES:
!    Accumulate Gx
!    --------------------------------------------------------------------
     if (choice>=0) then ! JWZ to check: I dont think 53 needs this
       do ilmn=1,nlmn
         il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
         scale=wt;if (il>1) scale=-scale
         if (cplex==2) then
           buffer_r = zero ; buffer_i = zero
           do ipw=1,npw
             buffer_r = buffer_r + scalr(ipw)*ffnl(ipw,1,ilmn)
             buffer_i = buffer_i + scali(ipw)*ffnl(ipw,1,ilmn)
           end do
           if (parity) then
             gx(1,ilmn,ia,ispinor) = scale*buffer_r ; gx(2,ilmn,ia,ispinor) = scale*buffer_i
           else
             gx(1,ilmn,ia,ispinor) =-scale*buffer_i ; gx(2,ilmn,ia,ispinor) = scale*buffer_r
           end if
         else
           if (parity) then
             buffer_r =  zero
             do ipw=1,npw
               buffer_r = buffer_r + scalr(ipw)*ffnl(ipw,1,ilmn)
             end do
             gx(1,ilmn,ia,ispinor) = scale*buffer_r
           else
             buffer_i = zero
             do ipw=1,npw
               buffer_i = buffer_i + scali(ipw)*ffnl(ipw,1,ilmn)
             end do
             gx(1,ilmn,ia,ispinor) =-scale*buffer_i
           end if
         end if
       end do
     end if

!    --------------------------------------------------------------------
!    CHOICE 2:
!    Accumulate dGxdt --- derivative wrt atm pos. ---
!    --------------------------------------------------------------------
     if (choice_==2) then
       if (signs==1) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_r2 = zero
             buffer_r3 = zero ; buffer_i1 = zero
             buffer_i2 = zero ; buffer_i3 = zero
             do ipw=1,npw
               aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
               aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
               buffer_r1 = buffer_r1 - aux_i*kpg(ipw,1)
               buffer_r2 = buffer_r2 - aux_i*kpg(ipw,2)
               buffer_r3 = buffer_r3 - aux_i*kpg(ipw,3)
               buffer_i1 = buffer_i1 + aux_r*kpg(ipw,1)
               buffer_i2 = buffer_i2 + aux_r*kpg(ipw,2)
               buffer_i3 = buffer_i3 + aux_r*kpg(ipw,3)
             end do
             if (parity) then
               dgxdt(1,1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(2,1,ilmn,ia,ispinor) = scale2*buffer_i1
               dgxdt(1,2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(2,2,ilmn,ia,ispinor) = scale2*buffer_i2
               dgxdt(1,3,ilmn,ia,ispinor) = scale2*buffer_r3
               dgxdt(2,3,ilmn,ia,ispinor) = scale2*buffer_i3
             else
               dgxdt(1,1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(2,1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(1,2,ilmn,ia,ispinor) =-scale2*buffer_i2
               dgxdt(2,2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(1,3,ilmn,ia,ispinor) =-scale2*buffer_i3
               dgxdt(2,3,ilmn,ia,ispinor) = scale2*buffer_r3
             end if
           else
             if (parity) then
               buffer_r1 = zero ; buffer_r2 = zero ; buffer_r3 = zero
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_r1 = buffer_r1 - aux_i*kpg(ipw,1)
                 buffer_r2 = buffer_r2 - aux_i*kpg(ipw,2)
                 buffer_r3 = buffer_r3 - aux_i*kpg(ipw,3)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(1,2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(1,3,ilmn,ia,ispinor) = scale2*buffer_r3
             else
               buffer_i1 = zero ; buffer_i2 = zero ; buffer_i3 = zero
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_i1 = buffer_i1 + aux_r*kpg(ipw,1)
                 buffer_i2 = buffer_i2 + aux_r*kpg(ipw,2)
                 buffer_i3 = buffer_i3 + aux_r*kpg(ipw,3)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(1,2,ilmn,ia,ispinor) =-scale2*buffer_i2
               dgxdt(1,3,ilmn,ia,ispinor) =-scale2*buffer_i3
             end if
           end if
         end do
       end if

!      ------------------------------------------
       if (signs==2) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_i1 = zero
             do ipw=1,npw
               aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
               aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
               buffer_r1 = buffer_r1 - aux_i*kpg(ipw,idir)
               buffer_i1 = buffer_i1 + aux_r*kpg(ipw,idir)
             end do
             if (parity) then
               dgxdt(1,1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(2,1,ilmn,ia,ispinor) = scale2*buffer_i1
             else
               dgxdt(1,1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(2,1,ilmn,ia,ispinor) = scale2*buffer_r1
             end if
           else
             if (parity) then
               buffer_r1 = zero
               do ipw=1,npw
                 buffer_r1 = buffer_r1 - scali(ipw)*ffnl(ipw,1,ilmn)*kpg(ipw,idir)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) = scale2*buffer_r1
             else
               buffer_i1 = zero
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scalr(ipw)*ffnl(ipw,1,ilmn)*kpg(ipw,idir)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) =-scale2*buffer_i1
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE 3:
!      Accumulate dGxdt --- derivative wrt strain ---
!      --------------------------------------------------------------------
     else if (choice_==3) then
       if (signs==1) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=half*scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_r2 = zero
             buffer_r3 = zero ; buffer_r4 = zero
             buffer_r5 = zero ; buffer_r6 = zero
             buffer_i1 = zero ; buffer_i2 = zero
             buffer_i3 = zero ; buffer_i4 = zero
             buffer_i5 = zero ; buffer_i6 = zero
             do ipw=1,npw
               aux_r2 = scalr(ipw)*ffnl(ipw,2,ilmn)
               aux_r3 = scalr(ipw)*ffnl(ipw,3,ilmn)
               aux_r4 = scalr(ipw)*ffnl(ipw,4,ilmn)
               aux_i2 = scali(ipw)*ffnl(ipw,2,ilmn)
               aux_i3 = scali(ipw)*ffnl(ipw,3,ilmn)
               aux_i4 = scali(ipw)*ffnl(ipw,4,ilmn)
               buffer_r1 = buffer_r1 + aux_r2*kpg(ipw,1)
               buffer_r2 = buffer_r2 + aux_r3*kpg(ipw,2)
               buffer_r3 = buffer_r3 + aux_r4*kpg(ipw,3)
               buffer_r4 = buffer_r4 + aux_r4*kpg(ipw,2) + aux_r3*kpg(ipw,3)
               buffer_r5 = buffer_r5 + aux_r4*kpg(ipw,1) + aux_r2*kpg(ipw,3)
               buffer_r6 = buffer_r6 + aux_r3*kpg(ipw,1) + aux_r2*kpg(ipw,2)
               buffer_i1 = buffer_i1 + aux_i2*kpg(ipw,1)
               buffer_i2 = buffer_i2 + aux_i3*kpg(ipw,2)
               buffer_i3 = buffer_i3 + aux_i4*kpg(ipw,3)
               buffer_i4 = buffer_i4 + aux_i4*kpg(ipw,2) + aux_i3*kpg(ipw,3)
               buffer_i5 = buffer_i5 + aux_i4*kpg(ipw,1) + aux_i2*kpg(ipw,3)
               buffer_i6 = buffer_i6 + aux_i3*kpg(ipw,1) + aux_i2*kpg(ipw,2)
             end do
             if (parity) then
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_r1
               dgxdt(2,1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(1,2,ilmn,ia,ispinor) =-scale*buffer_r2
               dgxdt(2,2,ilmn,ia,ispinor) =-scale*buffer_i2
               dgxdt(1,3,ilmn,ia,ispinor) =-scale*buffer_r3
               dgxdt(2,3,ilmn,ia,ispinor) =-scale*buffer_i3
               dgxdt(1,4,ilmn,ia,ispinor) =-scale2*buffer_r4
               dgxdt(2,4,ilmn,ia,ispinor) =-scale2*buffer_i4
               dgxdt(1,5,ilmn,ia,ispinor) =-scale2*buffer_r5
               dgxdt(2,5,ilmn,ia,ispinor) =-scale2*buffer_i5
               dgxdt(1,6,ilmn,ia,ispinor) =-scale2*buffer_r6
               dgxdt(2,6,ilmn,ia,ispinor) =-scale2*buffer_i6
             else
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(2,1,ilmn,ia,ispinor) =-scale*buffer_r1
               dgxdt(1,2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(2,2,ilmn,ia,ispinor) =-scale*buffer_r2
               dgxdt(1,3,ilmn,ia,ispinor) = scale*buffer_i3
               dgxdt(2,3,ilmn,ia,ispinor) =-scale*buffer_r3
               dgxdt(1,4,ilmn,ia,ispinor) = scale2*buffer_i4
               dgxdt(2,4,ilmn,ia,ispinor) =-scale2*buffer_r4
               dgxdt(1,5,ilmn,ia,ispinor) = scale2*buffer_i5
               dgxdt(2,5,ilmn,ia,ispinor) =-scale2*buffer_r5
               dgxdt(1,6,ilmn,ia,ispinor) = scale2*buffer_i6
               dgxdt(2,6,ilmn,ia,ispinor) =-scale2*buffer_r6
             end if
           else
             if (parity) then
               buffer_r1 = zero ; buffer_r2 = zero
               buffer_r3 = zero ; buffer_r4 = zero
               buffer_r5 = zero ; buffer_r6 = zero
               do ipw=1,npw
                 aux_r2 = scalr(ipw)*ffnl(ipw,2,ilmn)
                 aux_r3 = scalr(ipw)*ffnl(ipw,3,ilmn)
                 aux_r4 = scalr(ipw)*ffnl(ipw,4,ilmn)
                 buffer_r1 = buffer_r1 + aux_r2*kpg(ipw,1)
                 buffer_r2 = buffer_r2 + aux_r3*kpg(ipw,2)
                 buffer_r3 = buffer_r3 + aux_r4*kpg(ipw,3)
                 buffer_r4 = buffer_r4 + aux_r4*kpg(ipw,2) + aux_r3*kpg(ipw,3)
                 buffer_r5 = buffer_r5 + aux_r4*kpg(ipw,1) + aux_r2*kpg(ipw,3)
                 buffer_r6 = buffer_r6 + aux_r3*kpg(ipw,1) + aux_r2*kpg(ipw,2)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_r1
               dgxdt(1,2,ilmn,ia,ispinor) =-scale*buffer_r2
               dgxdt(1,3,ilmn,ia,ispinor) =-scale*buffer_r3
               dgxdt(1,4,ilmn,ia,ispinor) =-scale2*buffer_r4
               dgxdt(1,5,ilmn,ia,ispinor) =-scale2*buffer_r5
               dgxdt(1,6,ilmn,ia,ispinor) =-scale2*buffer_r6
             else
               buffer_i1 = zero ; buffer_i2 = zero
               buffer_i3 = zero ; buffer_i4 = zero
               buffer_i5 = zero ; buffer_i6 = zero
               do ipw=1,npw
                 aux_i2 = scali(ipw)*ffnl(ipw,2,ilmn)
                 aux_i3 = scali(ipw)*ffnl(ipw,3,ilmn)
                 aux_i4 = scali(ipw)*ffnl(ipw,4,ilmn)
                 buffer_i1 = buffer_i1 + aux_i2*kpg(ipw,1)
                 buffer_i2 = buffer_i2 + aux_i3*kpg(ipw,2)
                 buffer_i3 = buffer_i3 + aux_i4*kpg(ipw,3)
                 buffer_i4 = buffer_i4 + aux_i4*kpg(ipw,2) + aux_i3*kpg(ipw,3)
                 buffer_i5 = buffer_i5 + aux_i4*kpg(ipw,1) + aux_i2*kpg(ipw,3)
                 buffer_i6 = buffer_i6 + aux_i3*kpg(ipw,1) + aux_i2*kpg(ipw,2)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(1,3,ilmn,ia,ispinor) = scale*buffer_i3
               dgxdt(1,4,ilmn,ia,ispinor) = scale2*buffer_i4
               dgxdt(1,5,ilmn,ia,ispinor) = scale2*buffer_i5
               dgxdt(1,6,ilmn,ia,ispinor) = scale2*buffer_i6
             end if
           end if
         end do
       end if

!      ------------------------------------------
       if (signs==2) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_i1 = zero
             do ipw=1,npw
               buffer_r1 = buffer_r1 - scalr(ipw)*ffnl(ipw,2,ilmn)
               buffer_i1 = buffer_i1 - scali(ipw)*ffnl(ipw,2,ilmn)
             end do
             if (parity) then
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,1,ilmn,ia,ispinor) = scale*buffer_i1
             else
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,1,ilmn,ia,ispinor) = scale*buffer_r1
             end if
           else
             if (parity) then
               buffer_r1 = zero
               do ipw=1,npw
                 buffer_r1 = buffer_r1 - scalr(ipw)*ffnl(ipw,2,ilmn)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_r1
             else
               buffer_i1 = zero
               do ipw=1,npw
                 buffer_i1 = buffer_i1 - scali(ipw)*ffnl(ipw,2,ilmn)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_i1
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE 23:
!      Accumulate dGxdt --- derivative wrt atm pos. ---
!      Accumulate dGxdt --- derivative wrt strain ---
!      --------------------------------------------------------------------
     else if (choice_==23) then
       if (signs==1) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=half*scale;scale3=two_pi*scale
           if (cplex==2) then
             buffer_r1 =zero
             buffer_r2 = zero ; buffer_r3 = zero
             buffer_ra = zero ; buffer_rb = zero
             buffer_rc = zero ; buffer_rd = zero
             buffer_re = zero ; buffer_rf = zero
             buffer_i1 = zero
             buffer_i2 = zero ; buffer_i3 = zero
             buffer_ia = zero ; buffer_ib = zero
             buffer_ic = zero ; buffer_id = zero
             buffer_ie = zero ; buffer_if = zero
             do ipw=1,npw
               aux_r1 = scalr(ipw)*ffnl(ipw,1,ilmn)
               aux_r2 = scalr(ipw)*ffnl(ipw,2,ilmn)
               aux_r3 = scalr(ipw)*ffnl(ipw,3,ilmn)
               aux_r4 = scalr(ipw)*ffnl(ipw,4,ilmn)
               aux_i1 = scali(ipw)*ffnl(ipw,1,ilmn)
               aux_i2 = scali(ipw)*ffnl(ipw,2,ilmn)
               aux_i3 = scali(ipw)*ffnl(ipw,3,ilmn)
               aux_i4 = scali(ipw)*ffnl(ipw,4,ilmn)
               buffer_r1 = buffer_r1 - aux_i1*kpg(ipw,1)
               buffer_r2 = buffer_r2 - aux_i1*kpg(ipw,2)
               buffer_r3 = buffer_r3 - aux_i1*kpg(ipw,3)
               buffer_ra = buffer_ra + aux_r2*kpg(ipw,1)
               buffer_rb = buffer_rb + aux_r3*kpg(ipw,2)
               buffer_rc = buffer_rc + aux_r4*kpg(ipw,3)
               buffer_rd = buffer_rd + aux_r4*kpg(ipw,2) + aux_r3*kpg(ipw,3)
               buffer_re = buffer_re + aux_r4*kpg(ipw,1) + aux_r2*kpg(ipw,3)
               buffer_rf = buffer_rf + aux_r3*kpg(ipw,1) + aux_r2*kpg(ipw,2)
               buffer_i1 = buffer_i1 + aux_r1*kpg(ipw,1)
               buffer_i2 = buffer_i2 + aux_r1*kpg(ipw,2)
               buffer_i3 = buffer_i3 + aux_r1*kpg(ipw,3)
               buffer_ia = buffer_ia + aux_i2*kpg(ipw,1)
               buffer_ib = buffer_ib + aux_i3*kpg(ipw,2)
               buffer_ic = buffer_ic + aux_i4*kpg(ipw,3)
               buffer_id = buffer_id + aux_i4*kpg(ipw,2) + aux_i3*kpg(ipw,3)
               buffer_ie = buffer_ie + aux_i4*kpg(ipw,1) + aux_i2*kpg(ipw,3)
               buffer_if = buffer_if + aux_i3*kpg(ipw,1) + aux_i2*kpg(ipw,2)
             end do
             if (parity) then
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_ra
               dgxdt(2,1,ilmn,ia,ispinor) =-scale*buffer_ia
               dgxdt(1,2,ilmn,ia,ispinor) =-scale*buffer_rb
               dgxdt(2,2,ilmn,ia,ispinor) =-scale*buffer_ib
               dgxdt(1,3,ilmn,ia,ispinor) =-scale*buffer_rc
               dgxdt(2,3,ilmn,ia,ispinor) =-scale*buffer_ic
               dgxdt(1,4,ilmn,ia,ispinor) =-scale2*buffer_rd
               dgxdt(2,4,ilmn,ia,ispinor) =-scale2*buffer_id
               dgxdt(1,5,ilmn,ia,ispinor) =-scale2*buffer_re
               dgxdt(2,5,ilmn,ia,ispinor) =-scale2*buffer_ie
               dgxdt(1,6,ilmn,ia,ispinor) =-scale2*buffer_rf
               dgxdt(2,6,ilmn,ia,ispinor) =-scale2*buffer_if
               dgxdt(1,7,ilmn,ia,ispinor) = scale3*buffer_r1
               dgxdt(2,7,ilmn,ia,ispinor) = scale3*buffer_i1
               dgxdt(1,8,ilmn,ia,ispinor) = scale3*buffer_r2
               dgxdt(2,8,ilmn,ia,ispinor) = scale3*buffer_i2
               dgxdt(1,9,ilmn,ia,ispinor) = scale3*buffer_r3
               dgxdt(2,9,ilmn,ia,ispinor) = scale3*buffer_i3
             else
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_ia
               dgxdt(2,1,ilmn,ia,ispinor) =-scale*buffer_ra
               dgxdt(1,2,ilmn,ia,ispinor) = scale*buffer_ib
               dgxdt(2,2,ilmn,ia,ispinor) =-scale*buffer_rb
               dgxdt(1,3,ilmn,ia,ispinor) = scale*buffer_ic
               dgxdt(2,3,ilmn,ia,ispinor) =-scale*buffer_rc
               dgxdt(1,4,ilmn,ia,ispinor) = scale2*buffer_id
               dgxdt(2,4,ilmn,ia,ispinor) =-scale2*buffer_rd
               dgxdt(1,5,ilmn,ia,ispinor) = scale2*buffer_ie
               dgxdt(2,5,ilmn,ia,ispinor) =-scale2*buffer_re
               dgxdt(1,6,ilmn,ia,ispinor) = scale2*buffer_if
               dgxdt(2,6,ilmn,ia,ispinor) =-scale2*buffer_rf
               dgxdt(1,7,ilmn,ia,ispinor) =-scale3*buffer_i1
               dgxdt(2,7,ilmn,ia,ispinor) = scale3*buffer_r1
               dgxdt(1,8,ilmn,ia,ispinor) =-scale3*buffer_i2
               dgxdt(2,8,ilmn,ia,ispinor) = scale3*buffer_r2
               dgxdt(1,9,ilmn,ia,ispinor) =-scale3*buffer_i3
               dgxdt(2,9,ilmn,ia,ispinor) = scale3*buffer_r3
             end if
           else
             if (parity) then
               buffer_r1 = zero
               buffer_r2 = zero ; buffer_r3 = zero
               buffer_ra = zero ; buffer_rb = zero
               buffer_rc = zero ; buffer_rd = zero
               buffer_re = zero ; buffer_rf = zero
               do ipw=1,npw
                 aux_r2 = scalr(ipw)*ffnl(ipw,2,ilmn)
                 aux_r3 = scalr(ipw)*ffnl(ipw,3,ilmn)
                 aux_r4 = scalr(ipw)*ffnl(ipw,4,ilmn)
                 aux_i1 = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_r1 = buffer_r1 - aux_i1*kpg(ipw,1)
                 buffer_r2 = buffer_r2 - aux_i1*kpg(ipw,2)
                 buffer_r3 = buffer_r3 - aux_i1*kpg(ipw,3)
                 buffer_ra = buffer_ra + aux_r2*kpg(ipw,1)
                 buffer_rb = buffer_rb + aux_r3*kpg(ipw,2)
                 buffer_rc = buffer_rc + aux_r4*kpg(ipw,3)
                 buffer_rd = buffer_rd + aux_r4*kpg(ipw,2) + aux_r3*kpg(ipw,3)
                 buffer_re = buffer_re + aux_r4*kpg(ipw,1) + aux_r2*kpg(ipw,3)
                 buffer_rf = buffer_rf + aux_r3*kpg(ipw,1) + aux_r2*kpg(ipw,2)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_ra
               dgxdt(1,2,ilmn,ia,ispinor) =-scale*buffer_rb
               dgxdt(1,3,ilmn,ia,ispinor) =-scale*buffer_rc
               dgxdt(1,4,ilmn,ia,ispinor) =-scale2*buffer_rd
               dgxdt(1,5,ilmn,ia,ispinor) =-scale2*buffer_re
               dgxdt(1,6,ilmn,ia,ispinor) =-scale2*buffer_rf
               dgxdt(1,7,ilmn,ia,ispinor) = scale3*buffer_r1
               dgxdt(1,8,ilmn,ia,ispinor) = scale3*buffer_r2
               dgxdt(1,9,ilmn,ia,ispinor) = scale3*buffer_r3
             else
               buffer_i1 = zero
               buffer_i2 = zero ; buffer_i3 = zero
               buffer_ia = zero ; buffer_ib = zero
               buffer_ic = zero ; buffer_id = zero
               buffer_ie = zero ; buffer_if = zero
               do ipw=1,npw
                 aux_r1 = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i2 = scali(ipw)*ffnl(ipw,2,ilmn)
                 aux_i3 = scali(ipw)*ffnl(ipw,3,ilmn)
                 aux_i4 = scali(ipw)*ffnl(ipw,4,ilmn)
                 buffer_i1 = buffer_i1 + aux_r1*kpg(ipw,1)
                 buffer_i2 = buffer_i2 + aux_r1*kpg(ipw,2)
                 buffer_i3 = buffer_i3 + aux_r1*kpg(ipw,3)
                 buffer_ia = buffer_ia + aux_i2*kpg(ipw,1)
                 buffer_ib = buffer_ib + aux_i3*kpg(ipw,2)
                 buffer_ic = buffer_ic + aux_i4*kpg(ipw,3)
                 buffer_id = buffer_id + aux_i4*kpg(ipw,2) + aux_i3*kpg(ipw,3)
                 buffer_ie = buffer_ie + aux_i4*kpg(ipw,1) + aux_i2*kpg(ipw,3)
                 buffer_if = buffer_if + aux_i3*kpg(ipw,1) + aux_i2*kpg(ipw,2)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_ia
               dgxdt(1,2,ilmn,ia,ispinor) = scale*buffer_ib
               dgxdt(1,3,ilmn,ia,ispinor) = scale*buffer_ic
               dgxdt(1,4,ilmn,ia,ispinor) = scale2*buffer_id
               dgxdt(1,5,ilmn,ia,ispinor) = scale2*buffer_ie
               dgxdt(1,6,ilmn,ia,ispinor) = scale2*buffer_if
               dgxdt(1,7,ilmn,ia,ispinor) =-scale3*buffer_i1
               dgxdt(1,8,ilmn,ia,ispinor) =-scale3*buffer_i2
               dgxdt(1,9,ilmn,ia,ispinor) =-scale3*buffer_i3
             end if
           end if
         end do
       end if

!      ------------------------------------------
       if (signs==2) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_ra = zero ; buffer_ia = zero
             do ipw=1,npw
               buffer_ra = buffer_ra - scalr(ipw)*ffnl(ipw,2,ilmn)
               buffer_ia = buffer_ia - scali(ipw)*ffnl(ipw,2,ilmn)
             end do
             if (parity) then
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_ra
               dgxdt(2,1,ilmn,ia,ispinor) = scale*buffer_ia
             else
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_ia
               dgxdt(2,1,ilmn,ia,ispinor) = scale*buffer_ra
             end if
           else
             if (parity) then
               buffer_ra = zero
               do ipw=1,npw
                 buffer_ra = buffer_ra - scalr(ipw)*ffnl(ipw,2,ilmn)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_ra
             else
               buffer_ia = zero
               do ipw=1,npw
                 buffer_ia = buffer_ia - scali(ipw)*ffnl(ipw,2,ilmn)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_ia
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE 4 or 24:
!      Accumulate dGxdt --- derivative wrt atm pos. ---
!      Accumulate d2Gxdt --- 2nd derivative wrt 2 atm pos. ---
!      --------------------------------------------------------------------
     else if (choice_==4.or.choice_==24) then
       if (signs==1) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale;scale3=two_pi2*scale
           if (cplex==2) then
             buffer_r1 = zero
             buffer_r2 = zero ; buffer_r3 = zero
             buffer_ra = zero ; buffer_rb = zero
             buffer_rc = zero ; buffer_rd = zero
             buffer_re = zero ; buffer_rf = zero
             buffer_i1 = zero
             buffer_i2 = zero ; buffer_i3 = zero
             buffer_ia = zero ; buffer_ib = zero
             buffer_ic = zero ; buffer_id = zero
             buffer_ie = zero ; buffer_if = zero
             do ipw=1,npw
               aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
               aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
               buffer_r1 = buffer_r1 - aux_i*kpg(ipw,1)
               buffer_r2 = buffer_r2 - aux_i*kpg(ipw,2)
               buffer_r3 = buffer_r3 - aux_i*kpg(ipw,3)
               buffer_i1 = buffer_i1 + aux_r*kpg(ipw,1)
               buffer_i2 = buffer_i2 + aux_r*kpg(ipw,2)
               buffer_i3 = buffer_i3 + aux_r*kpg(ipw,3)
               buffer_ra = buffer_ra - aux_r*kpg(ipw,4)
               buffer_rb = buffer_rb - aux_r*kpg(ipw,5)
               buffer_rc = buffer_rc - aux_r*kpg(ipw,6)
               buffer_rd = buffer_rd - aux_r*kpg(ipw,7)
               buffer_re = buffer_re - aux_r*kpg(ipw,8)
               buffer_rf = buffer_rf - aux_r*kpg(ipw,9)
               buffer_ia = buffer_ia - aux_i*kpg(ipw,4)
               buffer_ib = buffer_ib - aux_i*kpg(ipw,5)
               buffer_ic = buffer_ic - aux_i*kpg(ipw,6)
               buffer_id = buffer_id - aux_i*kpg(ipw,7)
               buffer_ie = buffer_ie - aux_i*kpg(ipw,8)
               buffer_if = buffer_if - aux_i*kpg(ipw,9)
             end do
             if (parity) then
               dgxdt(1,1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(2,1,ilmn,ia,ispinor) = scale2*buffer_i1
               dgxdt(1,2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(2,2,ilmn,ia,ispinor) = scale2*buffer_i2
               dgxdt(1,3,ilmn,ia,ispinor) = scale2*buffer_r3
               dgxdt(2,3,ilmn,ia,ispinor) = scale2*buffer_i3
               d2gxdt(1,1,ilmn,ia,ispinor) = scale3*buffer_ra
               d2gxdt(2,1,ilmn,ia,ispinor) = scale3*buffer_ia
               d2gxdt(1,2,ilmn,ia,ispinor) = scale3*buffer_rb
               d2gxdt(2,2,ilmn,ia,ispinor) = scale3*buffer_ib
               d2gxdt(1,3,ilmn,ia,ispinor) = scale3*buffer_rc
               d2gxdt(2,3,ilmn,ia,ispinor) = scale3*buffer_ic
               d2gxdt(1,4,ilmn,ia,ispinor) = scale3*buffer_rd
               d2gxdt(2,4,ilmn,ia,ispinor) = scale3*buffer_id
               d2gxdt(1,5,ilmn,ia,ispinor) = scale3*buffer_re
               d2gxdt(2,5,ilmn,ia,ispinor) = scale3*buffer_ie
               d2gxdt(1,6,ilmn,ia,ispinor) = scale3*buffer_rf
               d2gxdt(2,6,ilmn,ia,ispinor) = scale3*buffer_if
             else
               dgxdt(1,1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(2,1,ilmn,ia,ispinor) = scale2*buffer_r1
               dgxdt(1,2,ilmn,ia,ispinor) =-scale2*buffer_i2
               dgxdt(2,2,ilmn,ia,ispinor) = scale2*buffer_r2
               dgxdt(1,3,ilmn,ia,ispinor) =-scale2*buffer_i3
               dgxdt(2,3,ilmn,ia,ispinor) = scale2*buffer_r3
               d2gxdt(1,1,ilmn,ia,ispinor) =-scale3*buffer_ia
               d2gxdt(2,1,ilmn,ia,ispinor) = scale3*buffer_ra
               d2gxdt(1,2,ilmn,ia,ispinor) =-scale3*buffer_ib
               d2gxdt(2,2,ilmn,ia,ispinor) = scale3*buffer_rb
               d2gxdt(1,3,ilmn,ia,ispinor) =-scale3*buffer_ic
               d2gxdt(2,3,ilmn,ia,ispinor) = scale3*buffer_rc
               d2gxdt(1,4,ilmn,ia,ispinor) =-scale3*buffer_id
               d2gxdt(2,4,ilmn,ia,ispinor) = scale3*buffer_rd
               d2gxdt(1,5,ilmn,ia,ispinor) =-scale3*buffer_ie
               d2gxdt(2,5,ilmn,ia,ispinor) = scale3*buffer_re
               d2gxdt(1,6,ilmn,ia,ispinor) =-scale3*buffer_if
               d2gxdt(2,6,ilmn,ia,ispinor) = scale3*buffer_rf
             end if
           else
             if (parity) then
               buffer_r1 = zero
               buffer_r2 = zero ; buffer_r3 = zero
               buffer_ra = zero ; buffer_rb = zero
               buffer_rc = zero ; buffer_rd = zero
               buffer_re = zero ; buffer_rf = zero
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_r1 = buffer_r1 - aux_i*kpg(ipw,1)
                 buffer_r2 = buffer_r2 - aux_i*kpg(ipw,2)
                 buffer_r3 = buffer_r3 - aux_i*kpg(ipw,3)
                 buffer_ra = buffer_ra - aux_r*kpg(ipw,4)
                 buffer_rb = buffer_rb - aux_r*kpg(ipw,5)
                 buffer_rc = buffer_rc - aux_r*kpg(ipw,6)
                 buffer_rd = buffer_rd - aux_r*kpg(ipw,7)
                 buffer_re = buffer_re - aux_r*kpg(ipw,8)
                 buffer_rf = buffer_rf - aux_r*kpg(ipw,9)
               end do
               dgxdt(1,1,ilmn,ia,ispinor)  = scale2*buffer_r1
               dgxdt(1,2,ilmn,ia,ispinor)  = scale2*buffer_r2
               dgxdt(1,3,ilmn,ia,ispinor)  = scale2*buffer_r3
               d2gxdt(1,1,ilmn,ia,ispinor) = scale3*buffer_ra
               d2gxdt(1,2,ilmn,ia,ispinor) = scale3*buffer_rb
               d2gxdt(1,3,ilmn,ia,ispinor) = scale3*buffer_rc
               d2gxdt(1,4,ilmn,ia,ispinor) = scale3*buffer_rd
               d2gxdt(1,5,ilmn,ia,ispinor) = scale3*buffer_re
               d2gxdt(1,6,ilmn,ia,ispinor) = scale3*buffer_rf
             else
               buffer_i1 = zero
               buffer_i2 = zero ; buffer_i3 = zero
               buffer_ia = zero ; buffer_ib = zero
               buffer_ic = zero ; buffer_id = zero
               buffer_ie = zero ; buffer_if = zero
               do ipw=1,npw
                 aux_r = scalr(ipw)*ffnl(ipw,1,ilmn)
                 aux_i = scali(ipw)*ffnl(ipw,1,ilmn)
                 buffer_i1 = buffer_i1 + aux_r*kpg(ipw,1)
                 buffer_i2 = buffer_i2 + aux_r*kpg(ipw,2)
                 buffer_i3 = buffer_i3 + aux_r*kpg(ipw,3)
                 buffer_ia = buffer_ia - aux_i*kpg(ipw,4)
                 buffer_ib = buffer_ib - aux_i*kpg(ipw,5)
                 buffer_ic = buffer_ic - aux_i*kpg(ipw,6)
                 buffer_id = buffer_id - aux_i*kpg(ipw,7)
                 buffer_ie = buffer_ie - aux_i*kpg(ipw,8)
                 buffer_if = buffer_if - aux_i*kpg(ipw,9)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) =-scale2*buffer_i1
               dgxdt(1,2,ilmn,ia,ispinor) =-scale2*buffer_i2
               dgxdt(1,3,ilmn,ia,ispinor) =-scale2*buffer_i3
               d2gxdt(1,1,ilmn,ia,ispinor) =-scale3*buffer_ia
               d2gxdt(1,2,ilmn,ia,ispinor) =-scale3*buffer_ib
               d2gxdt(1,3,ilmn,ia,ispinor) =-scale3*buffer_ic
               d2gxdt(1,4,ilmn,ia,ispinor) =-scale3*buffer_id
               d2gxdt(1,5,ilmn,ia,ispinor) =-scale3*buffer_ie
               d2gxdt(1,6,ilmn,ia,ispinor) =-scale3*buffer_if
             end if
           end if
         end do
       end if

!      --------------------------------------------------------------------
!      CHOICE 5, 51, 52, 53:
!      Accumulate dGxdt --- derivative wrt k ---
!      --------------------------------------------------------------------
     else if (choice_==5 .or. choice_==51 .or. choice_==52 .or. choice_==53) then
       if (signs==1) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_r2 = zero
             buffer_r3 = zero
             buffer_i1 = zero ; buffer_i2 = zero
             buffer_i3 = zero
             do ipw=1,npw
               buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,2,ilmn)
               buffer_r2 = buffer_r2 + scalr(ipw)*ffnl(ipw,3,ilmn)
               buffer_r3 = buffer_r3 + scalr(ipw)*ffnl(ipw,4,ilmn)
               buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,2,ilmn)
               buffer_i2 = buffer_i2 + scali(ipw)*ffnl(ipw,3,ilmn)
               buffer_i3 = buffer_i3 + scali(ipw)*ffnl(ipw,4,ilmn)
             end do
             if (parity) then
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,1,ilmn,ia,ispinor) = scale*buffer_i1
               dgxdt(1,2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(2,2,ilmn,ia,ispinor) = scale*buffer_i2
               dgxdt(1,3,ilmn,ia,ispinor) = scale*buffer_r3
               dgxdt(2,3,ilmn,ia,ispinor) = scale*buffer_i3
             else
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(1,2,ilmn,ia,ispinor) =-scale*buffer_i2
               dgxdt(2,2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(1,3,ilmn,ia,ispinor) =-scale*buffer_i3
               dgxdt(2,3,ilmn,ia,ispinor) = scale*buffer_r3
             end if
           else
             if (parity) then
               buffer_r1 = zero ; buffer_r2 = zero
               buffer_r3 = zero
               do ipw=1,npw
                 buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,2,ilmn)
                 buffer_r2 = buffer_r2 + scalr(ipw)*ffnl(ipw,3,ilmn)
                 buffer_r3 = buffer_r3 + scalr(ipw)*ffnl(ipw,4,ilmn)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(1,2,ilmn,ia,ispinor) = scale*buffer_r2
               dgxdt(1,3,ilmn,ia,ispinor) = scale*buffer_r3
             else
               buffer_i1 = zero ; buffer_i2 = zero
               buffer_i3 = zero
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,2,ilmn)
                 buffer_i2 = buffer_i2 + scali(ipw)*ffnl(ipw,3,ilmn)
                 buffer_i3 = buffer_i3 + scali(ipw)*ffnl(ipw,4,ilmn)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(1,2,ilmn,ia,ispinor) =-scale*buffer_i2
               dgxdt(1,3,ilmn,ia,ispinor) =-scale*buffer_i3
             end if
           end if
         end do
       end if

!      ------------------------------------------
       if (signs==2) then
!        here in the simple case of a single direction idir, the original call
!        to mkffnl would have produced ffnl with a single derivative in direction
!        idir in location ffnl(:,2,:). Case 53 though need multiple directions
!        and here ffnl will contain the derivative information in locations 2,3, and 4
!        corresponding to idir = 1, 2, and 3. Moreover, choice 53 needs the derivatives
!        in direction idir+1 and idir-1. The parameter vector
!        ffnl_dir_dat contains the necessary translations in locations 1,2 for idir=1;
!        3,4 for idir=2; and 5,6 for idir=3.
         if (choice_==5 .or. choice_==51 .or. choice_==52) ffnl_dir = 2 ! idir derivative
         if (choice_==53) ffnl_dir = ffnl_dir_dat(2*idir-1) ! idir+1 derivative
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           if (cplex==2) then
             buffer_r1 = zero ; buffer_i1 = zero
             do ipw=1,npw
               buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,ffnl_dir,ilmn)
               buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,ffnl_dir,ilmn)
             end do
             if (parity) then
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_r1
               dgxdt(2,1,ilmn,ia,ispinor) = scale*buffer_i1
             else
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_i1
               dgxdt(2,1,ilmn,ia,ispinor) = scale*buffer_r1
             end if
           else
             if (parity) then
               buffer_r1 = zero
               do ipw=1,npw
                 buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,ffnl_dir,ilmn)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) = scale*buffer_r1
             else
               buffer_i1 = zero
               do ipw=1,npw
                 buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,ffnl_dir,ilmn)
               end do
               dgxdt(1,1,ilmn,ia,ispinor) =-scale*buffer_i1
             end if
           end if
         end do
         if (choice_==53) then ! accumlate idir-1 derivative
           ffnl_dir = ffnl_dir_dat(2*idir)
           do ilmn=1,nlmn
             il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
             scale=wt;if (il>1) scale=-scale
             if (cplex==2) then
               buffer_r1 = zero ; buffer_i1 = zero
               do ipw=1,npw
                 buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,ffnl_dir,ilmn)
                 buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,ffnl_dir,ilmn)
               end do ! loop over npw
               if (parity) then
                 dgxdt(1,2,ilmn,ia,ispinor) = scale*buffer_r1
                 dgxdt(2,2,ilmn,ia,ispinor) = scale*buffer_i1
               else
                 dgxdt(1,2,ilmn,ia,ispinor) =-scale*buffer_i1
                 dgxdt(2,2,ilmn,ia,ispinor) = scale*buffer_r1
               end if ! end if on parity
             else
               if (parity) then
                 buffer_r1 = zero
                 do ipw=1,npw
                   buffer_r1 = buffer_r1 + scalr(ipw)*ffnl(ipw,ffnl_dir,ilmn)
                 end do ! loop over npw
                 dgxdt(1,2,ilmn,ia,ispinor) = scale*buffer_r1
               else
                 buffer_i1 = zero
                 do ipw=1,npw
                   buffer_i1 = buffer_i1 + scali(ipw)*ffnl(ipw,ffnl_dir,ilmn)
                 end do
                 dgxdt(1,2,ilmn,ia,ispinor) =-scale*buffer_i1
               end if !end if on parity
             end if ! end if on cplx
           end do ! end loop over nlmn
         end if ! end second part of choice=53
       end if

!      --------------------------------------------------------------------
!      CHOICE 6:
!      Accumulate dGxdt --- derivative wrt atm pos. ---
!      Accumulate dGxdt --- derivative wrt strain ---
!      Accumulate d2Gxdt --- 2nd derivative wrt 2 strains ---
!      Accumulate d2Gxdt --- 2nd derivative wrt strain and atm. pos ---
!      --------------------------------------------------------------------
     else if (choice_==6) then
       if (signs==1) then
         ABI_ALLOCATE(scalarr,(npw,10))
         ABI_ALLOCATE(scalari,(npw,10))
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4);parity=(mod(il,2)==0)
           scale=wt;if (il>1) scale=-scale
           scale2=two_pi*scale;scale3=half*scale
           scale4=quarter*scale;scale5=pi*scale

           do mu=1,10
             do ipw=1,npw
               scalarr(ipw,mu)=scalr(ipw)*ffnl(ipw,mu,ilmn)
               scalari(ipw,mu)=scali(ipw)*ffnl(ipw,mu,ilmn)
             end do
           end do

!          ===== Accumulate derivative of Gx wrt atm pos. =====
           ishift=6
           if (cplex==2) then
             do mu=1,3
               buffer_r = zero ; buffer_i = zero
               do ipw=1,npw
                 buffer_r = buffer_r - scalari(ipw,1)*kpg(ipw,mu)
                 buffer_i = buffer_i + scalarr(ipw,1)*kpg(ipw,mu)
               end do
               if (parity) then
                 dgxdt(1,ishift+mu,ilmn,ia,ispinor) = scale2*buffer_r
                 dgxdt(2,ishift+mu,ilmn,ia,ispinor) = scale2*buffer_i
               else
                 dgxdt(1,ishift+mu,ilmn,ia,ispinor) =-scale2*buffer_i
                 dgxdt(2,ishift+mu,ilmn,ia,ispinor) = scale2*buffer_r
               end if
             end do
           else
             if (parity) then
               do mu=1,3
                 buffer_r = zero
                 do ipw=1,npw
                   buffer_r = buffer_r - scalari(ipw,1)*kpg(ipw,mu)
                 end do
                 dgxdt(1,ishift+mu,ilmn,ia,ispinor) = scale2**buffer_r
               end do
             else
               do mu=1,3
                 buffer_i = zero
                 do ipw=1,npw
                   buffer_i = buffer_i + scalarr(ipw,1)*kpg(ipw,mu)
                 end do
                 dgxdt(1,ishift+mu,ilmn,ia,ispinor) =-scale2**buffer_i
               end do
             end if
           end if

!          ===== Accumulate derivative of Gx wrt strain =====
           if(cplex==2) then
             do mu=1,6
               mua=alpha(mu);mub=beta(mu)
               buffer_r = zero ; buffer_i = zero
               do ipw=1,npw
                 buffer_r = buffer_r + scalarr(ipw,1+mua)*kpg(ipw,mub) &
&                 + scalarr(ipw,1+mub)*kpg(ipw,mua)
                 buffer_i = buffer_i + scalari(ipw,1+mua)*kpg(ipw,mub) &
&                 + scalari(ipw,1+mub)*kpg(ipw,mua)
               end do
               if (parity) then
                 dgxdt(1,mu,ilmn,ia,ispinor) =-scale3*buffer_r
                 dgxdt(2,mu,ilmn,ia,ispinor) =-scale3*buffer_i
               else
                 dgxdt(1,mu,ilmn,ia,ispinor) = scale3*buffer_i
                 dgxdt(2,mu,ilmn,ia,ispinor) =-scale3*buffer_r
               end if
             end do
           else
             if (parity) then
               do mu=1,6
                 mua=alpha(mu);mub=beta(mu)
                 buffer_r = zero
                 do ipw=1,npw
                   buffer_r = buffer_r + scalarr(ipw,1+mua)*kpg(ipw,mub) &
&                   + scalarr(ipw,1+mub)*kpg(ipw,mua)
                 end do
                 dgxdt(1,mu,ilmn,ia,ispinor) =-scale3*buffer_r
               end do
             else
               do mu=1,6
                 mua=alpha(mu);mub=beta(mu)
                 buffer_i = zero
                 do ipw=1,npw
                   buffer_i = buffer_i + scalari(ipw,1+mua)*kpg(ipw,mub) &
&                   + scalari(ipw,1+mub)*kpg(ipw,mua)
                 end do
                 dgxdt(1,mu,ilmn,ia,ispinor) = scale3*buffer_i
               end do
             end if
           end if

!          ===== Accumulate 2nd derivative of Gx wrt two strains =====
           if(cplex==2) then
             do mub=1,6
               nub1=alpha(mub);nub2=beta(mub)
               do mua=1,6
                 mu=mua+6*(mub-1)
                 nua1=alpha(mua);nua2=beta(mua)
                 gama=gamma(nub1,nua1)
                 gamb=gamma(nub2,nua1)
                 gamc=gamma(nub1,nua2)
                 gamd=gamma(nub2,nua2)
                 buffer_r = zero ; buffer_i = zero
                 do ipw=1,npw
                   kpga=kpg(ipw,nub1)
                   kpgb=kpg(ipw,nub2)
                   kpgc=kpg(ipw,nua1)
                   kpgd=kpg(ipw,nua2)
                   buffer_r = buffer_r + scalarr(ipw,4+gama)*kpgb*kpgd+scalarr(ipw,4+gamb)*kpga*kpgd &
&                   + scalarr(ipw,4+gamc)*kpgb*kpgc+scalarr(ipw,4+gamd)*kpga*kpgc
                   buffer_i = buffer_i + scalari(ipw,4+gama)*kpgb*kpgd+scalari(ipw,4+gamb)*kpga*kpgd &
&                   + scalari(ipw,4+gamc)*kpgb*kpgc+scalari(ipw,4+gamd)*kpga*kpgc
                 end do
                 if (parity) then
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale4*buffer_r
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale4*buffer_i
                 else
                   d2gxdt(1,mu,ilmn,ia,ispinor) =-scale4*buffer_i
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale4*buffer_r
                 end if
               end do
             end do
           else
             if (parity) then
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,6
                   mu=mua+6*(mub-1)
                   nua1=alpha(mua);nua2=beta(mua)
                   gama=gamma(nub1,nua1)
                   gamb=gamma(nub2,nua1)
                   gamc=gamma(nub1,nua2)
                   gamd=gamma(nub2,nua2)
                   buffer_r = zero
                   do ipw=1,npw
                     kpga=kpg(ipw,nub1)
                     kpgb=kpg(ipw,nub2)
                     kpgc=kpg(ipw,nua1)
                     kpgd=kpg(ipw,nua2)
                     buffer_r = buffer_r + scalarr(ipw,4+gama)*kpgb*kpgd+scalarr(ipw,4+gamb)*kpga*kpgd &
&                     + scalarr(ipw,4+gamc)*kpgb*kpgc+scalarr(ipw,4+gamd)*kpga*kpgc
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale4*buffer_r
                 end do
               end do
             else
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,6
                   mu=mua+6*(mub-1)
                   nua1=alpha(mua);nua2=beta(mua)
                   gama=gamma(nub1,nua1)
                   gamb=gamma(nub2,nua1)
                   gamc=gamma(nub1,nua2)
                   gamd=gamma(nub2,nua2)
                   buffer_i = zero
                   do ipw=1,npw
                     kpga=kpg(ipw,nub1)
                     kpgb=kpg(ipw,nub2)
                     kpgc=kpg(ipw,nua1)
                     kpgd=kpg(ipw,nua2)
                     buffer_i = buffer_i + scalari(ipw,4+gama)*kpgb*kpgd+scalari(ipw,4+gamb)*kpga*kpgd &
&                     + scalari(ipw,4+gamc)*kpgb*kpgc+scalari(ipw,4+gamd)*kpga*kpgc
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) =-scale4*buffer_i
                 end do
               end do
             end if
           end if

!          ===== Accumulate 2nd derivative of Gx wrt strain and atm pos. =====
           ishift=36
           if(cplex==2) then
             do mub=1,6
               nub1=alpha(mub);nub2=beta(mub)
               do mua=1,3
                 mu=ishift+mua+3*(mub-1)
                 buffer_r = zero ; buffer_i = zero
                 do ipw=1,npw
                   buffer_r = buffer_r + kpg(ipw,mua)*(scalari(ipw,1+nub1)*kpg(ipw,nub2) + &
&                   scalari(ipw,1+nub2)*kpg(ipw,nub1))
                   buffer_i = buffer_i + kpg(ipw,mua)*(scalarr(ipw,1+nub1)*kpg(ipw,nub2) + &
&                   scalarr(ipw,1+nub2)*kpg(ipw,nub1))
                 end do
                 if (parity) then
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale5*buffer_r
                   d2gxdt(2,mu,ilmn,ia,ispinor) =-scale5*buffer_i
                 else
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale5*buffer_i
                   d2gxdt(2,mu,ilmn,ia,ispinor) = scale5*buffer_r
                 end if
               end do
             end do
           else
             if (parity) then
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,3
                   mu=ishift+mua+3*(mub-1)
                   buffer_r = zero
                   do ipw=1,npw
                     buffer_r = buffer_r + kpg(ipw,mua)*(scalari(ipw,1+nub1)*kpg(ipw,nub2) + &
&                     scalari(ipw,1+nub2)*kpg(ipw,nub1))
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale5*buffer_r
                 end do
               end do
             else
               do mub=1,6
                 nub1=alpha(mub);nub2=beta(mub)
                 do mua=1,3
                   mu=ishift+mua+3*(mub-1)
                   buffer_i = zero
                   do ipw=1,npw
                     buffer_i = buffer_i + kpg(ipw,mua)*(scalarr(ipw,1+nub1)*kpg(ipw,nub2) + &
&                     scalarr(ipw,1+nub2)*kpg(ipw,nub1))
                   end do
                   d2gxdt(1,mu,ilmn,ia,ispinor) = scale5*buffer_i
                 end do
               end do
             end if
           end if
         end do
         ABI_DEALLOCATE(scalarr)
         ABI_DEALLOCATE(scalari)
       end if

!      --------------------------------------------------------------------
!      END CHOICES
!      --------------------------------------------------------------------
     end if

!    End loop on atoms
   end do
!  $OMP END DO
!  $OMP END PARALLEL

!  End loop on spinorial components
 end do

!Deallocate temporary space
 ABI_DEALLOCATE(scali)
 ABI_DEALLOCATE(scalr)


!Has to reduce arrays in case of FFT parallelization
 if (mpi_enreg%paral_compil_fft==1) then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
   call timab(48,1,tsec)
   if (choice>=0) then
     call xsum_mpi(gx,spaceComm,ierr)
   end if
   if (choice_>1) then
     call xsum_mpi(dgxdt,spaceComm,ierr)
   end if
   if (choice_==4.or.choice_==24.or.choice_==6) then
     call xsum_mpi(d2gxdt,spaceComm,ierr)
   end if
   call timab(48,2,tsec)
   mpi_enreg%paral_level=old_paral_level
 end if

end subroutine opernla_ylm
!!***
