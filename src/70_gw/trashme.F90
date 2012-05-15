!{\src2tex{textfont=tt}}
!!****f* ABINIT/assemblychi0q0_sym
!! NAME
!! assemblychi0q0_sym
!!
!! FUNCTION
!! Update the independent particle susceptibility at q==0 for the contribution
!! of one pair of occupied-unoccupied band, for each frequency.
!! This routine takes advantage of the symmetries of the little group of the external q-point
!! to symmetrize the contribution arising from the input k-point located in the IBZ_q.
!! It computes:
!!
!!   $ \chi_0(G1,G2,io) = \chi_0(G1,G2,io)+\sum_S (rhotwg(G1)*rhotwg^\dagger(G2))*green_w(io) $
!!
!! where S is a symmetry in reciprocal space.
!! The matrix elements of the gradient operator and [V_{nl},r] are symmetrized as well.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ik_bz=Index of the k-point whose contribution has to be added to chi0.
!!  isym_kbz=Index of the symmetry such that k_bz = IS k_ibz
!!  itim_kbz=2 if time-reversal has to be used to obtain k_bz, 1 otherwise.
!!  npwepG0=Maximum number of G vectors
!!  nqlwl=Number of small q-points used to deal with the optical limit.
!!  qlwl(3,nqlwl)=reciprocal space coordinates of the q-point for long-wavelength limit treatment.
!!  rhotwg(npwepG0*nspinor**2)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
!!  rhotwx(3,nspinor**2)=Matrix element of the operator -i[H,r]/(e1-e2) in reciprocal lattice units.
!!  green_w(nomega)=Frequency dependent part of the Green function.
!!  Ltg_q<little_group_type>=Info on the little group associated to the external q-point.
!!    %timrev=2 it time-reversal is used, 1 otherwise.
!!    %nsym_sg=Number of space group symmetries.
!!    %wtksym(2,nsym,nkbz)=1 if the symmetry (with or without time-reversal) must be considered for this k-point.
!!  Gsph_epsG0<Gvectors_type> Information on the "enlarged" G-sphere used for chi0, it contains umklapp G0 vectors
!!    %ng=number of G vectors in the enlarged sphere, actually MUST be equal to the size of rhotwg.
!!    %rottbm1(ng,2,nsym)=index of (IR)^{-1} G where I is the identity or the inversion.
!!    %phmGt(ng,nsym)=phase factors associated to non-simmorphic operations.
!!  Cryst<Crystal_structure>=Structure defining the unit cell and its symmetries
!!    %nsym=Number of symmetries.
!!    %symrec(3,3,nsym)=Symmetry operation in reciprocal space (reduced coordinates)
!!  Ep<Epsilonm1_parameters>=Parameters of the chi0 calculation.
!!     %npwe=number of plane waves in chi0.
!!     %symchi=1 if symmetrization has to be performed.
!!     %nomega=number of frequencies in chi0.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0(npwe,npwe,nomega)= Updated independent-particle susceptibility matrix in reciprocal space at q==0.
!!  lwing(Ep%npwe*Ep%nI,Ep%nomega,3)=Lower wing (calculated only if nqlwl > 1 )
!!  uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)=Upper wing (calculated only if nqlwl > 1 )
!!
!! NOTES
!!
!!  1) Symmetrization of the oscilator matrix elements.
!!    If  Sq = q then  M_G( Sk,q)= e^{-i(q+G).\tau} M_{ S^-1G}  (k,q)
!!    If -Sq = q then  M_G(-Sk,q)= e^{-i(q+G).\tau} M_{-S^-1G}^*(k,q)
!!
!!    In the case of umklapps:
!!    If  Sq = q+G0 then  M_G( Sk,q)= e^{-i(q+G).\tau} M_{ S^-1(G-G0}   (k,q)
!!    If -Sq = q+G0 then  M_G(-Sk,q)= e^{-i(q+G).\tau} M_{-S^-1(G-G0)}^*(k,q)
!!
!!  In the equation below there is no need to take into account the phases due to q.t
!!  as they cancel each other in the scalar product ==> only phmGt(G,isym)=e^{-iG.\tau} is needed.
!!
!!  2) Symmetrization of the matrix elements of the position operator.
!!
!!    <Sk,b|\vec r| Sk,b'> = <k b| R\vec r + \tau|k b'>
!!
!!     where S is one of the symrec operation, R and \tau is the corresponding
!!     operation in real space. The term involving the fractional translation is zero provided that b /= b'.
!!
!! PARENTS
!!
!! CHILDREN
!!      matrginv,xgerc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine assemblychi0q0_sym(nqlwl,qlwl,ik_bz,isym_kbz,itim_kbz,gwcomp,nspinor,npwepG0,Ep,Cryst,Ltg_q,Gsph_epsG0,&
& chi0,rhotwx,rhotwg,green_w,green_enhigh_w,deltaf_b1b2,lwing,uwing)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

 use m_gwdefs,   only : GW_TOL_DOCC, czero_gw, epsilonm1_parameters
 use m_blas,     only : xgerc
 use m_geometry, only : vdotw
 use m_crystal,  only : crystal_structure
 use m_gsphere,  only : gvectors_type
 use m_bz_mesh,  only : little_group

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assemblychi0q0_sym'
 use interfaces_32_util
 use interfaces_70_gw, except_this_one => assemblychi0q0_sym
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,isym_kbz,itim_kbz
 integer,intent(in) :: npwepG0,nqlwl,nspinor,gwcomp
 real(dp),intent(in) :: deltaf_b1b2
 type(Little_group),intent(in) :: Ltg_q
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Crystal_structure),intent(in) :: Cryst
 type(Epsilonm1_parameters),intent(in) :: Ep
!arrays
 real(dp),intent(in) :: qlwl(3,nqlwl)
 complex(gwpc),intent(inout) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
 complex(dpc),intent(in) :: green_w(Ep%nomega),green_enhigh_w(Ep%nomega)
 complex(dpc),intent(inout) :: lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
 complex(dpc),intent(inout) :: uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)

!Local variables-------------------------------
!scalars
 integer :: itim,io,isym,igp,ig,iqlwl
 integer :: jj,ii,s_jj,pad_jj,pad_ii
 complex(gwpc) :: dd,mqg0,mqg0_sym,rhotwg0_bkp
 character(len=500) :: msg
!arrays
 integer,pointer :: Sm1G(:)
 real(dp) :: opinv(3,3),qrot(3)
 real(dp) :: b1(3),b2(3),b3(3)
 complex(gwpc),allocatable :: rhotwg_sym(:),rhotwg_I(:),rhotwg_J(:)
 complex(gwpc),allocatable :: rhotwg_sym_star(:),rhotwg_star(:)
 complex(gwpc),pointer :: phmGt(:)

!************************************************************************

 b1=two_pi*Gsph_epsG0%gprimd(:,1)
 b2=two_pi*Gsph_epsG0%gprimd(:,2)
 b3=two_pi*Gsph_epsG0%gprimd(:,3)

 SELECT CASE (Ep%symchi)

 CASE (0) ! Do not use symmetries.

   if (nspinor==1) then
    ! * Accumulate over the full BZ i.e
    !    chi0(G1,G2,io) = chi0(G1,G2,io) + (rhotwg(G1)*CONJG(rhotwg(G2)))*green_w(io)
    ! * The non-analytic term is symmetrized for this k-point in the BZ according to:
    !    rhotwg(1)= S^-1q * rhotwx_ibz
    !    rhotwg(1)=-S^-1q * CONJG(rhotwx_ibz) if time-reversal is used.
    opinv(:,:)=REAL(Cryst%symrec(:,:,isym_kbz),dp)
    call matrginv(opinv,3,3)
    qrot = (3-2*itim_kbz) * MATMUL(opinv,qlwl(:,1))

    rhotwg(1)=dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3) !TODO get rid of this
    if (itim_kbz==2) rhotwg(1)=CONJG(rhotwg(1))

    if (gwcomp==1) then ! Leave the head and wings uncorrected (does not matter much)
      if (ABS(deltaf_b1b2) < GW_TOL_DOCC) rhotwg(1)=czero_gw
      do igp=1,Ep%npwe
        chi0(1,igp,:) = chi0(1,igp,:) + rhotwg(1) *CONJG(rhotwg(igp))*green_enhigh_w(:)
      end do
      do ig=2,Ep%npwe
        chi0(ig,1,:)  = chi0(ig,1,:)  + rhotwg(ig)*CONJG(rhotwg(1))  *green_enhigh_w(:)
      end do
    end if

    ! Multiply elements G1,G2 of rhotwg_sym by green_w(io) and accumulate in chi0(G1,G2,io)
    do io=1,Ep%nomega
      dd=green_w(io)
      call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg,1,rhotwg,1,chi0(:,:,io),Ep%npwe)
    end do

    ! === Accumulate heads and wings for each small q ===
    ! * For better performance, this part is not done if nqlwl==1
    !   lwing and uwing will be filled in cchi0q0 after the MPI collective sum
    ! FIXME extrapolar method should be checked!!
    !do io=1,Ep%nomega
    ! lwing(:,io,1) =  chi0(:,1,io)
    ! uwing(:,io,1) =  chi0(1,:,io)
    !end do

    if (nqlwl>1.and..FALSE.) then
      rhotwg0_bkp = rhotwg(1) ! Save G=0 value of the first q
      ABI_ALLOCATE(rhotwg_star,(Ep%npwe))
      rhotwg_star = CONJG(rhotwg(1:Ep%npwe))

      do iqlwl=2,nqlwl
        qrot =  (3-2*itim_kbz) * MATMUL(opinv,qlwl(:,iqlwl))
        mqg0 = dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3) !TODO get rid of this
        if (itim_kbz==2) mqg0=CONJG(mqg0)
        rhotwg     (1) =mqg0
        rhotwg_star(1) =CONJG(mqg0)
        !
        ! here we might take advantage of Hermiticity along Im axis in RPA (see mkG0w)
        do io=1,Ep%nomega
          lwing(:,io,iqlwl) = lwing(:,io,iqlwl) + rhotwg     (1:Ep%npwe) * CONJG(mqg0) * green_w(io)
          uwing(:,io,iqlwl) = uwing(:,io,iqlwl) + rhotwg_star(1:Ep%npwe) *       mqg0  * green_w(io)
        end do
      end do ! iqlwl

      ABI_DEALLOCATE(rhotwg_star)
      rhotwg(1) = rhotwg0_bkp ! Reinstate previous value of rhotwg(1).
    end if ! nqlwl > 1

  else ! spinorial case
    ABI_ALLOCATE(rhotwg_I,(Ep%npwe))
    ABI_ALLOCATE(rhotwg_J,(Ep%npwe))

    ABI_CHECK(nqlwl==1,"nqlwl/=1 Not implemented")

    ! I can use symmetries to loop over the upper triangle but
    ! this makes using BLAS more difficult
    ! Important NOTE: treatment of q-->0 limit is correct only
    ! for i=j=0. Other components require additional terms.

    do jj=1,Ep%nJ
      s_jj=1 ; if (jj==4) s_jj=-1
      pad_jj=(jj-1)*Ep%npwe
      call mkrhotwg_sigma(jj,nspinor,Ep%npwe,rhotwg,rhotwg_J)

      rhotwg_J(1) = q0limit(jj,qlwl(:,1),nspinor,rhotwx,b1,b2,b3)
      !TODO RECHECK this
      if (itim_kbz==2) rhotwg_J(1)=-CONJG(rhotwg_J(1))

      do ii=1,Ep%nI
        pad_ii=(ii-1)*Ep%npwe

        if (ii/=jj) then
          call mkrhotwg_sigma(ii,nspinor,Ep%npwe,rhotwg,rhotwg_I)
          rhotwg_I(1) = q0limit(ii,qlwl(:,1),nspinor,rhotwx,b1,b2,b3)
          if (itim_kbz==2) rhotwg_I(1)=-CONJG(rhotwg_I(1))
        else
          rhotwg_I(:)=rhotwg_J(:)
        end if

        do io=1,Ep%nomega
          dd = s_jj*green_w(io)
          call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_I,1,rhotwg_J,1,&
&           chi0(pad_ii+1:pad_ii+Ep%npwe,pad_jj+1:pad_jj+Ep%npwe,io),Ep%npwe)
        end do

      end do !ii
    end do !jj

    ABI_DEALLOCATE(rhotwg_I)
    ABI_DEALLOCATE(rhotwg_J)
  end if

 CASE (1) ! Use symmetries to reconstruct the integrand.
   if (nspinor==1) then
     ABI_ALLOCATE(rhotwg_sym,(Ep%npwe))

     ! === Loop over symmetries of the space group and time-reversal ===
     do isym=1,Ltg_q%nsym_sg
       do itim=1,Ltg_q%timrev

         if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then
           ! === This operation belongs to the little group and has to be considered to reconstruct the BZ ===
           ! TODO this is a hot-spot, should add a test on the umklapp
           phmGt => Gsph_epsG0%phmGt  (1:Ep%npwe,isym) ! In the 2 lines below note the slicing (1:npwe)
           Sm1G  => Gsph_epsG0%rottbm1(1:Ep%npwe,itim,isym)

           opinv(:,:)=REAL(Cryst%symrec(:,:,isym),dp)
           call matrginv(opinv,3,3)
           qrot = (3-2*itim) * MATMUL(opinv,qlwl(:,1))

           SELECT CASE (itim)

           CASE (1)
             rhotwg_sym(1:Ep%npwe)=rhotwg(Sm1G(1:Ep%npwe))*phmGt(1:Ep%npwe)
             rhotwg_sym(1)=dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3)

           CASE (2)
             rhotwg_sym(1:Ep%npwe)=CONJG(rhotwg(Sm1G(1:Ep%npwe)))*phmGt(1:Ep%npwe)
             rhotwg_sym(1)=CONJG(dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3))

           CASE DEFAULT
             write(msg,'(a,i3)')'Wrong value of itim= ',itim
             MSG_BUG(msg)
           END SELECT

           if (gwcomp==1) then ! Leave the head and wings uncorrected (does not matter much)
             if (ABS(deltaf_b1b2) < GW_TOL_DOCC) rhotwg_sym(1)=czero_gw
             do igp=1,Ep%npwe
               chi0(1,igp,:) = chi0(1,igp,:) + rhotwg_sym(1) *CONJG(rhotwg_sym(igp))*green_enhigh_w(:)
             end do
             do ig=2,Ep%npwe
               chi0(ig,1,:)  = chi0(ig,1,:)  + rhotwg_sym(ig)*CONJG(rhotwg_sym(1))  *green_enhigh_w(:)
             end do
           end if

           ! Multiply elements G1,G2 of rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)
           do io=1,Ep%nomega
             dd=green_w(io)
             call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
           end do

           ! === Accumulate heads and wings for each small q ===
           ! * For better performance, this part is not done if nqlwl==1
           !   lwing and uwing will be filled in cchi0q0 after the MPI collective sum
           ! FIXME extrapolar method should be checked!!
           if (nqlwl>1.and..FALSE.) then
             ABI_ALLOCATE(rhotwg_sym_star,(Ep%npwe))
             rhotwg_sym_star = CONJG(rhotwg_sym)

             do iqlwl=2,nqlwl
               qrot = (3-2*itim) * MATMUL(opinv,qlwl(:,iqlwl))
               mqg0_sym = dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3)
               if (itim==2) mqg0_sym = CONJG(mqg0_sym)

               rhotwg_sym     (1) =       mqg0_sym
               rhotwg_sym_star(1) = CONJG(mqg0_sym)

               ! here we might take advantage of Hermiticity along Im axis in RPA (see mkG0w)
               do io=1,Ep%nomega
                 lwing(:,io,iqlwl) = lwing(:,io,iqlwl) + rhotwg_sym     (1:Ep%npwe) * CONJG(mqg0_sym) * green_w(io)
                 uwing(:,io,iqlwl) = uwing(:,io,iqlwl) + rhotwg_sym_star(1:Ep%npwe) *       mqg0_sym  * green_w(io)
               end do
             end do !iqlwl

             ABI_DEALLOCATE(rhotwg_sym_star)
           end if !nqlwl>1

         end if !wtksym
       end do !itim
     end do !isym

     ABI_DEALLOCATE(rhotwg_sym)

   else  !spinorial case
     write(msg,'(a,i3)')' symchi=1 with spinor not implemented '
     MSG_BUG(msg)
     ABI_CHECK(nqlwl==1,"nqlwl/=1 Not implemented")
   end if

 CASE DEFAULT
   write(msg,'(a,i3)')'Wrong value of symchi= ',Ep%symchi
   MSG_BUG(msg)
 END SELECT

end subroutine assemblychi0q0_sym
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/assemblychi0_sym
!! NAME
!! assemblychi0_sym
!!
!! FUNCTION
!! Update the independent particle susceptibility for the contribution
!! of one pair of occupied-unoccupied band, for each frequency.
!! If symchi=1 the expression is symmetrized taking into account the symmetries
!! of the little group associated to the external q-point.
!! Compute chi0(G1,G2,io)=chi0(G1,G2,io)+\sum_S \hat S (rhotwg(G1)*rhotwg*(G2))*green_w(io)
!! where S are the symmetries of the little group associated to the external q-point.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nspinor=Number of spinorial components.
!!  ik_bz=Index of the k-point in the BZ array whose contribution has to be symmetrized and added to cchi0
!!  npwepG0=Maximum number of G vectors taking into account possible umklapp G0, ie enlarged sphere G-G0
!!  rhotwg(npwe*nspinor**2)=Oscillator matrix elements for this k-point and the transition that has to be summed
!!  green_w(nomega)=frequency dependent part coming from the green function
!!  Gsph_epsG0<Gvectors_type> Information on the "enlarged" G-sphere used for chi0, it contains umklapp G0 vectors
!!    %ng=number of G vectors in the enlarged sphere, actually MUST be equal to the size of rhotwg
!!    %rottbm1(ng,2,nsym)=index of (IR)^{-1} G where I is the identity or the inversion
!!    %phmGt(ng,nsym)=phase factors associated to non-simmorphic operations
!!  Ltg_q<little_group_type>=Info on the little group associated to the external q-point.
!!    %timrev=2 it time-reversal is used, 1 otherwise
!!    %nsym_sg=Number of space group symmetries
!!    %wtksym(2,nsym,nkbz)=1 if the symmetry (with or without time-reversal) must be considered for this k-point
!!    %flag_umklp(timrev,nsym)= flag for umklapp processes
!!      if 1 that the particular operation (IS) requires a G_o to preserve Q, 0 otherwise
!!    %igmG0(npwepG0,timrev,nsym) index of G-G0 in the array gvec
!!  Ep<Epsilonm1_parameters>=Parameters related to the calculation of chi0/epsilon^-1
!!    %symchi
!!    %nomega=number of frequencies
!!    %npwe=number of plane waves for epsilon (input variable)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0(npwe,npwe,nomega)=independent-particle susceptibility matrix in reciprocal space
!!
!! PARENTS
!!      cchi0,cchi0q0_intraband,check_completeness
!!
!! CHILDREN
!!      matrginv,xgerc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine assemblychi0_sym(ik_bz,nspinor,Ep,Ltg_q,green_w,npwepG0,rhotwg,Gsph_epsG0,chi0)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

 use m_gwdefs,   only : epsilonm1_parameters
 use m_blas,     only : xgerc
 use m_gsphere,  only : gvectors_type
 use m_bz_mesh,  only : little_group

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assemblychi0_sym'
 use interfaces_70_gw, except_this_one => assemblychi0_sym
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,npwepG0,nspinor
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Little_group),intent(in) :: Ltg_q
 type(Epsilonm1_parameters),intent(in) :: Ep
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 complex(dpc),intent(in) :: green_w(Ep%nomega)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: itim,io,isym
 integer :: jj,ii,s_jj,pad_jj,pad_ii
 complex(gwpc) :: dd
 character(len=500) :: msg
!arrays
 integer,pointer :: gmG0(:)
 integer,allocatable :: Sm1_gmG0(:)
 complex(gwpc),allocatable :: rhotwg_sym(:),rhotwg_I(:),rhotwg_J(:)
 complex(gwpc),pointer :: phmGt(:)

! *************************************************************************

 SELECT CASE (Ep%symchi)

 CASE (0) ! Do not use symmetries
   if (nspinor==1) then
     do io=1,Ep%nomega
       dd=green_w(io)
       call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg,1,rhotwg,1,chi0(:,:,io),Ep%npwe)
     end do

   else ! spinorial case
     ABI_ALLOCATE(rhotwg_I,(Ep%npwe))
     ABI_ALLOCATE(rhotwg_J,(Ep%npwe))

     ! I can use symmetries to loop over the upper triangle but
     ! this makes using BLAS more difficult

     do jj=1,Ep%nJ
       s_jj=1 ; if (jj==4) s_jj=-1
       pad_jj=(jj-1)*Ep%npwe
       call mkrhotwg_sigma(jj,nspinor,Ep%npwe,rhotwg,rhotwg_J)
       !
       do ii=1,Ep%nI
         pad_ii=(ii-1)*Ep%npwe

         if (ii/=jj) then
          call mkrhotwg_sigma(ii,nspinor,Ep%npwe,rhotwg,rhotwg_I)
         else
          rhotwg_I(:)=rhotwg_J(:)
         end if

         do io=1,Ep%nomega
          dd = s_jj*green_w(io)
          call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_I,1,rhotwg_J,1,chi0(pad_ii+1:pad_ii+Ep%npwe,pad_jj+1:pad_jj+Ep%npwe,io),Ep%npwe)
         end do
         !
       end do !ii
     end do !jj

     ABI_DEALLOCATE(rhotwg_I)
     ABI_DEALLOCATE(rhotwg_J)
   end if

 CASE (1) ! Use symmetries to reconstruct the integrand in the BZ.
   !
   ! Notes on the symmetrization of the oscillator matrix elements
   !  If  Sq = q then  M_G^( Sk,q)= e^{-i(q+G).t} M_{ S^-1G}  (k,q)
   !  If -Sq = q then  M_G^(-Sk,q)= e^{-i(q+G).t} M_{-S^-1G}^*(k,q)
   !
   ! In case of an umklapp process
   !  If  Sq = q+G0 then  M_G( Sk,q)= e^{-i(q+G).t} M_{ S^-1(G-G0}   (k,q)
   !  If -Sq = q+G0 then  M_G(-Sk,q)= e^{-i(q+G).t} M_{-S^-1(G-G0)}^*(k,q)
   !
   ! Ltg_q%igmG0(ig,itim,isym) contains the index of G-G0 where ISq=q+G0
   ! Note that there is no need to take into account the phases due to q,
   ! They cancel in the scalar product ==> phmGt(G,isym)=e^{-iG\cdot t}
   !
   ! Mind the slicing of %rottbm1(npwepG0,timrev,nsym) and %phmGt(npwepG0,nsym) as
   ! these arrays, usually, do not conform to rho_twg_sym(npw) !
   !
   ABI_ALLOCATE(rhotwg_sym,(Ep%npwe))
   ABI_ALLOCATE(Sm1_gmG0  ,(Ep%npwe))
   !
   ! === Loop over symmetries of the space group and time-reversal ===
   do isym=1,Ltg_q%nsym_sg
     do itim=1,Ltg_q%timrev
       !
       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then
         ! === This operation belongs to the little group and has to be used to reconstruct the BZ ===
         ! * In the following 3 lines mind the slicing (1:npwe)
         ! TODO this is a hot-spot, should add a test on the umklapp
         !
         phmGt => Gsph_epsG0%phmGt(1:Ep%npwe,isym)
         gmG0  => Ltg_q%igmG0     (1:Ep%npwe,itim,isym)
         Sm1_gmG0(1:Ep%npwe)=Gsph_epsG0%rottbm1(gmG0(1:Ep%npwe),itim,isym)

         SELECT CASE (itim)
         CASE (1)
           rhotwg_sym(1:Ep%npwe)=rhotwg(Sm1_gmG0)*phmGt(1:Ep%npwe)
         CASE (2)
           rhotwg_sym(1:Ep%npwe)=CONJG(rhotwg(Sm1_gmG0))*phmGt(1:Ep%npwe)
         CASE DEFAULT
           write(msg,'(a,i3)')'Wrong itim= ',itim
           MSG_BUG(msg)
         END SELECT
         !
         ! Multiply rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)
         do io=1,Ep%nomega
           dd=green_w(io)
           call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
         end do
         !
       end if
       !
     end do
   end do

   ABI_DEALLOCATE(rhotwg_sym)
   ABI_DEALLOCATE(Sm1_gmG0)

 CASE DEFAULT
   write(msg,'(a,i3)')'Wrong symchi= ',Ep%symchi
   MSG_BUG(msg)
 END SELECT

end subroutine assemblychi0_sym
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/mkrhotwg_sigma
!! NAME
!! mkrhotwg_sigma
!!
!! FUNCTION
!!  Helper function used to calculate selected linear combination
!!  of the oscillator matrix elements in the case of noncollinear magnetism.
!!
!! INPUTS
!!  ii=Index selecting the particolar combination of spin components.
!!  npw=Number of plane-waves in the oscillators.
!!  nspinor=Number of spinorial components.
!!  rhotwg(npw*nspinor**2)=OScillator matrix elements.
!!
!! OUTPUT
!!  rhotwg_I(npw)=Required linear combination of the oscillator matrix elements.
!!
!! PARENTS
!!      accumulate_chi0_q0,trashme
!!
!! CHILDREN
!!      matrginv,xgerc
!!
!! SOURCE

subroutine mkrhotwg_sigma(ii,nspinor,npw,rhotwg,rhotwg_I)

 use m_profiling

 use defs_basis
 use m_errors

 use m_gwdefs, only : j_gw

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkrhotwg_sigma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,npw,nspinor
!arrays
 complex(gwpc),intent(in) :: rhotwg(npw*nspinor**2)
 complex(gwpc),intent(out) :: rhotwg_I(npw)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************

 SELECT CASE (ii)
 CASE (1) ! $ M_0 = M_{\up,\up} + M_{\down,\down} $
   rhotwg_I(:) = rhotwg(1:npw) + rhotwg(npw+1:2*npw)
 CASE (2) ! $ M_z = M_{\up,\up} - M_{\down,\down} $
   rhotwg_I(:) = rhotwg(1:npw) - rhotwg(npw+1:2*npw)
 CASE (3) ! $ M_x = M_{\up,\down} + M_{\down,\up} $
   rhotwg_I(:) = ( rhotwg(2*npw+1:3*npw) + rhotwg(3*npw+1:4*npw) )
 CASE (4) ! $ M_y = i * (M_{\up,\down} -M_{\down,\up}) $
   rhotwg_I(:) = (rhotwg(2*npw+1:3*npw) - rhotwg(3*npw+1:4*npw) )*j_gw
 CASE DEFAULT
   write(msg,'(a,i3)')'Wrong ii value= ',ii
   MSG_BUG(msg)
 END SELECT

end subroutine mkrhotwg_sigma
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/assemblychi0sfq0
!! NAME
!! assemblychi0sfq0
!!
!! FUNCTION
!! Update the spectral function of the independent particle susceptibility at q==0 for the contribution
!! of one pair of occupied-unoccupied band, for each frequency.
!! If symchi==1, the symmetries belonging to the little group of the external point q are used
!! to reconstrunct the contributions in the full Brillouin zone. In this case, the equation implented is:
!!
!!  $ chi0(G1,G2,io)=chi0(G1,G2,io)+\sum_S (rhotwg(G1)*rhotwg^\dagger(G2))* \delta(\omega -trans) $
!!
!! where S is a symmetry belonging to the little group of q.
!! The subroutine also performs the symmetrization of the matrix elements of the
!! gradient operator and of the commutator [V_{nl},r] with the position operator.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ikbz=Index in the BZ of the k-point whose contribution to chi0 has to be added,
!!   if we use symmetries, the contribution to chi0 by this k-point has to be symmetrized.
!!  isym_kbz=Index of the symmetry such as k_bz = IS k_ibz
!!  itim_kbz=2 if time-reversal has to be used to obtain k_bz, 1 otherwise.
!!  my_wl,my_wr=min and Max frequency index treated by this processor.
!!  npwe=Number of plane waves used to describe chi0.
!!  npwepG0=Maximum number of G vectors to account for umklapps.
!!  nomega=Number of frequencies in the imaginary part.
!!  nqlwl=Number of q-points used for the optical limit.
!!  qlwl(3,nqlwl)=Reciprocal space coordinates of the q-points for the long-wavelength limit treatment.
!!  rhotwg(npwepG0*nspinor**2)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
!!  rhotwx(3,nspinor**2)=Matrix elements of the gradient and of the commutator of the non-local operator with
!!    the position operator. The second term is present only if inclvkb=1,2.
!!  Gsph_epsG0<Gvectors_type> Information on the "enlarged" G-sphere used for chi0, it contains umklapp G0 vectors
!!    %ng=number of G vectors in the enlarged sphere. It MUST be equal to the size of rhotwg
!!    %rottbm1(ng,2,nsym)=index of (IR)^{-1} G where I is the identity or the inversion
!!    %phmGt(ng,nsym)=phase factors associated to non-symmorphic operations
!!  Ltg_q<little_group_type>=Info on the little group associated to the external q-point.
!!    %timrev=2 it time-reversal is used, 1 otherwise
!!    %nsym_sg=Number of space group symmetries
!!    %wtksym(2,nsym,nkbz)=1 if the symmetry (with or without time-reversal) must be considered for this k-point
!!    %flag_umklp(timrev,nsym)= flag for umklapp processes
!!     if 1 that the particular operation (IS) requires a G_o to preserve Q, 0 otherwise
!! Cryst<Crystal_structure>=Info on unit cell and it symmetries
!!    %nsym=Number of symmetry operations.
!!    %symrec(3,3,nsym)=Symmetry operations in reciprocal space (reduced coordinates).
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0sf(npwe,npwe,my_wl:my_wr)=Updated spectral function at q==0.
!!  lwing_sf(npwe,nomega,nqlwl)=Updated lower wing of the spectral function.
!!  uwing_sf(npwe,nomega,nqlwl)=Updated Upper wing of the spectral function.
!!
!! PARENTS
!!
!! CHILDREN
!!      matrginv,xgerc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine assemblychi0sfq0(nqlwl,qlwl,ikbz,isym_kbz,itim_kbz,nspinor,symchi,npwepG0,npwe,Cryst,Ltg_q,Gsph_epsG0,&
& factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx,rhotwg,nomegasf,chi0sf,lwing_sf,uwing_sf)

 use m_profiling

 use defs_basis
 use m_errors

 use m_gwdefs,   only : epsilonm1_parameters
 use m_blas,     only : xgerc
 use m_crystal,  only : crystal_structure
 use m_gsphere,  only : gvectors_type
 use m_bz_mesh,  only : little_group

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assemblychi0sfq0'
 use interfaces_32_util
 use interfaces_70_gw, except_this_one => assemblychi0sfq0
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikbz,my_wl,my_wr,nomegasf,npwe,npwepG0,nqlwl,nspinor
 integer,intent(in) :: isym_kbz,itim_kbz,symchi,iomegal,iomegar
 real(dp),intent(in) :: factocc,wl,wr
 type(Little_group),intent(in) :: Ltg_q
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 real(dp),intent(in) :: qlwl(3,nqlwl)
 complex(gwpc),intent(inout) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(in) :: rhotwx(3)
 complex(gwpc),intent(inout) :: chi0sf(npwe,npwe,my_wl:my_wr)
 complex(dpc),intent(inout) :: lwing_sf(npwe,my_wl:my_wr,3)
 complex(dpc),intent(inout) :: uwing_sf(npwe,my_wl:my_wr,3)

!Local variables-------------------------------
!scalars
 integer :: itim,isym,iqlwl
 complex(gwpc) :: num
 complex(gwpc) :: mqg0,mqg0_sym,rhotwg0_bkp
 character(len=500) :: msg
!arrays
 integer,pointer :: Sm1G(:)
 real(dp) :: opinv(3,3),qrot(3),b1(3),b2(3),b3(3)
 complex(gwpc),allocatable :: rhotwg_sym(:)
 complex(gwpc),allocatable :: rhotwg_sym_star(:),rhotwg_star(:)
 complex(gwpc),pointer :: phmGt(:)
!************************************************************************

 if (iomegal<my_wl .or. iomegar>my_wr) then
   write(msg,'(3a,2(a,i5,a,i5))')ch10,&
&    ' assemblychi0sfq0 : Indeces out of boundary ',ch10,&
&    '  my_wl = ',my_wl,' iomegal = ',iomegal,ch10,&
&    '  my_wr = ',my_wr,' iomegar = ',iomegar,ch10
   MSG_PERS_BUG(msg)
 end if

 b1(:)=two_pi*Gsph_epsG0%gprimd(:,1)
 b2(:)=two_pi*Gsph_epsG0%gprimd(:,2)
 b3(:)=two_pi*Gsph_epsG0%gprimd(:,3)

 SELECT CASE (symchi)

 CASE (0)
   !
   ! === Calculation without symmetries ===
   ! * rhotwg(1)= R^-1q*rhotwx_ibz
   ! * rhotwg(1)=-R^-1q*conjg(rhotwx_ibz) for inversion
   ! FIXME My equation reads  -iq* <cSk|\nabla|vSk> = -i \transpose S <ck_i|\nabla\|vk_i>
   if (nspinor==1) then
     opinv(:,:)=REAL(Cryst%symrec(:,:,isym_kbz),dp)
     call matrginv(opinv,3,3)
     qrot =  (3-2*itim_kbz) * MATMUL(opinv,qlwl(:,1))
     rhotwg(1)=dotproductqrc(qrot,rhotwx,b1,b2,b3)
     if (itim_kbz==2) rhotwg(1)=CONJG(rhotwg(1))

     if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
       num=-wl*factocc ! Num is single precision needed for cgerc check factocc
       call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegal),npwe)
     end if
     ! Last point, must accumulate left point but not the right one
     if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
       num=-wr*factocc
       call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegar),npwe)
     end if
     !
     ! === Accumulate heads and wings for each small q ===
     ! * For better performance, this part is not done if nqlwl==1
     !   lwing and uwing will be filled in cchi0q0 after the MPI collective sum
     !
     if (nqlwl>1.and..FALSE.) then
       rhotwg0_bkp = rhotwg(1) ! Save G=0 value of the first q
       ABI_ALLOCATE(rhotwg_star,(npwepG0))
       rhotwg_star = CONJG(rhotwg(1:npwepG0))

       do iqlwl=2,nqlwl
         qrot =  (3-2*itim_kbz) * MATMUL(opinv,qlwl(:,iqlwl))
         mqg0 = dotproductqrc(qrot,rhotwx,b1,b2,b3) !TODO get rid of this
         if (itim_kbz==2) mqg0=CONJG(mqg0)
         rhotwg     (1) =mqg0
         rhotwg_star(1) =CONJG(mqg0)
         !
         if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
           num=-wl*factocc ! Num is single precision needed for cgerc check factocc
           lwing_sf(:,iomegal,iqlwl) = lwing_sf(:,iomegal,iqlwl) + rhotwg     (1:npwepG0) * CONJG(mqg0) * num
           uwing_sf(:,iomegal,iqlwl) = uwing_sf(:,iomegal,iqlwl) + rhotwg_star(1:npwepG0) *       mqg0  * num
         end if
         !
         ! Last point, must accumulate left point but not the right one
         if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
           num=-wr*factocc
           lwing_sf(:,iomegar,iqlwl) = lwing_sf(:,iomegar,iqlwl) + rhotwg     (1:npwepG0) * CONJG(mqg0) * num
           uwing_sf(:,iomegar,iqlwl) = uwing_sf(:,iomegar,iqlwl) + rhotwg_star(1:npwepG0) *       mqg0  * num
         end if
       end do ! iqlwl

       ABI_DEALLOCATE(rhotwg_star)
       rhotwg(1) = rhotwg0_bkp ! Reinstate previous value of rhotwg(1).
     end if !nqlwl

   else ! spinorial case
     msg="Spectral method + nspinor==2 not implemented"
     MSG_BUG(msg)
   end if


 CASE (1)
   ! === Notes on the symmetrization of oscillator matrix elements ===
   ! If  Sq = q then  M_G( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1G}  (k,q)
   ! If -Sq = q then  M_G(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1G}^*(k,q)
   !
   ! In case of an umklapp process
   ! If  Sq = q+G_o then  M_G( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1(G-G_o}   (k,q)
   ! If -Sq = q+G_o then  M_G(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1(G-G-o)}^*(k,q)
   !
   ! rhotwg(1)= R^-1q*rhotwx_ibz
   ! rhotwg(1)=-R^-1q*conjg(rhotwx_ibz) for inversion
   !
   if (nspinor==1) then
     ABI_ALLOCATE(rhotwg_sym,(npwe))
     !
     ! === Loop over symmetries of the space group and time-reversal ===
     do isym=1,Ltg_q%nsym_sg
       do itim=1,Ltg_q%timrev

         if (Ltg_q%wtksym(itim,isym,ikbz)==1) then
           ! === This operation belongs to the little group and has to be considered to reconstruct the BZ ===
           ! TODO this is a hot-spot, should add a test on the umklapp
           !
           phmGt => Gsph_epsG0%phmGt(1:npwe,isym) ! In these 2 lines mind the slicing (1:npwe)
           Sm1G  => Gsph_epsG0%rottbm1(1:npwe,itim,isym)

           opinv(:,:)=REAL(Cryst%symrec(:,:,isym),dp)
           call matrginv(opinv,3,3)
           qrot = (3-2*itim) * MATMUL(opinv,qlwl(:,1))

           SELECT CASE (itim)

           CASE (1)
             rhotwg_sym(1:npwe)=rhotwg(Sm1G(1:npwe))*phmGt(1:npwe)
             rhotwg_sym(1)=dotproductqrc(qrot,rhotwx,b1,b2,b3)

           CASE (2)
             rhotwg_sym(1:npwe)=CONJG(rhotwg(Sm1G(1:npwe)))*phmGt(1:npwe)
             rhotwg_sym(1)=CONJG(dotproductqrc(qrot,rhotwx,b1,b2,b3))

           CASE DEFAULT
             write(msg,'(a,i4)')'Wrong value of itim= ',itim
             MSG_BUG(msg)
           END SELECT
           !
           ! === Multiply elements G,Gp of rhotwg_sym*num and accumulate in chi0sf(G,Gp,io) ===
           if (wl<huge(0.0_dp)*1.d-11) then
             num=-wl*factocc
             call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegal),npwe)
           end if
           !
           ! Last point, must accumulate left point but not the right one
           if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
             num=-wr*factocc
             call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegar),npwe)
           end if

           ! === Accumulate heads and wings for each small q ===
           ! * For better performance, this part is not done if nqlwl==1
           !   lwing and uwing will be filled in cchi0q0 after the MPI collective sum
           if (nqlwl>1.and..FALSE.) then
             ABI_ALLOCATE(rhotwg_sym_star,(npwe))
             rhotwg_sym_star = CONJG(rhotwg_sym(1:npwe))

             do iqlwl=2,nqlwl
               qrot =  (3-2*itim_kbz) * MATMUL(opinv,qlwl(:,iqlwl))
               mqg0_sym = dotproductqrc(qrot,rhotwx,b1,b2,b3) !TODO get rid of this
               if (itim_kbz==2) mqg0_sym=CONJG(mqg0_sym)
               rhotwg_sym     (1) =mqg0_sym
               rhotwg_sym_star(1) =CONJG(mqg0_sym)
               !
               if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
                 num=-wl*factocc ! Num is single precision needed for cgerc check factocc
                 lwing_sf(:,iomegal,iqlwl) = lwing_sf(:,iomegal,iqlwl) + rhotwg_sym_star(1:npwe) * CONJG(mqg0_sym) * num
                 uwing_sf(:,iomegal,iqlwl) = uwing_sf(:,iomegal,iqlwl) + rhotwg_sym_star(1:npwe) *       mqg0_sym  * num
               end if
               ! Last point, must accumulate left point but not the right one
               if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then
                 num=-wr*factocc
                 lwing_sf(:,iomegar,iqlwl) = lwing_sf(:,iomegar,iqlwl) + rhotwg_sym_star(1:npwe) * CONJG(mqg0_sym) * num
                 uwing_sf(:,iomegar,iqlwl) = uwing_sf(:,iomegar,iqlwl) + rhotwg_sym_star(1:npwe) *       mqg0_sym  * num
               end if
             end do ! iqlwl

             ABI_DEALLOCATE(rhotwg_sym_star)
           end if !nqlwl

         end if !wtksym
       end do !inv
     end do !isym
     ABI_DEALLOCATE(rhotwg_sym)

   else ! spinorial case
     msg="Spectral method + nspinor==2 not implemented"
     MSG_BUG(msg)
   end if

 CASE DEFAULT
   write(msg,'(a,i4)')'Wrong value of symchi= ',symchi
   MSG_BUG(msg)
 END SELECT

end subroutine assemblychi0sfq0
!!***
