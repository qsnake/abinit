!{\src2tex{textfont=tt}}
!!****f* ABINIT/accumulate_chi0_q0
!! NAME
!! accumulate_chi0_q0
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
!!  rhotwg(npwepG0*nspinor**2)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
!!  rhotwx(3,nspinor**2)=Matrix element of the operator $-i[H,r]/(e1-e2) = -i r$ in reciprocal lattice units.
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
!!  chi0_head(3,3,Ep%nomega)=Head.
!!  chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)=Lower wing.
!!  chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)=Upper wing.
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
!!      cchi0q0
!!
!! CHILDREN
!!      mkrhotwg_sigma,xgerc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine accumulate_chi0_q0(ik_bz,isym_kbz,itim_kbz,gwcomp,gw_eet,nspinor,npwepG0,Ep,Cryst,Ltg_q,Gsph_epsG0,&
& chi0,rhotwx,rhotwg,green_w,green_enhigh_w,deltaf_b1b2,chi0_head,chi0_lwing,chi0_uwing)

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
#define ABI_FUNC 'accumulate_chi0_q0'
 use interfaces_70_gw, except_this_one => accumulate_chi0_q0
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,isym_kbz,itim_kbz,npwepG0,nspinor,gwcomp,gw_eet
 real(dp),intent(in) :: deltaf_b1b2
 type(Little_group),intent(in) :: Ltg_q
 type(Gvectors_type),intent(in) :: Gsph_epsG0 
 type(Crystal_structure),intent(in) :: Cryst
 type(Epsilonm1_parameters),intent(in) :: Ep
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
 complex(dpc),intent(in) :: green_w(Ep%nomega),green_enhigh_w(Ep%nomega)
 complex(dpc),intent(inout) :: chi0_head(3,3,Ep%nomega)
 complex(dpc),intent(inout) :: chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
 complex(dpc),intent(inout) :: chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)

!Local variables-------------------------------
!scalars
 integer :: itim,io,isym,idir,jdir,jj,ii,s_jj,pad_jj,pad_ii
 complex(gwpc) :: dd
 character(len=500) :: msg
!arrays
 integer,pointer :: Sm1G(:) 
 complex(dpc) :: mir_kbz(3)
 complex(gwpc),allocatable :: rhotwg_sym(:),rhotwg_I(:),rhotwg_J(:) 
 complex(gwpc),pointer :: phmGt(:)

!************************************************************************

 ABI_UNUSED(deltaf_b1b2)

 SELECT CASE (Ep%symchi)
 
 CASE (0) ! Do not use symmetries.

   if (nspinor==1) then
    ! * Symmetrize rhotwg in the full BZ and accumulate over the full BZ i.e. 
    !     chi0(G1,G2,io) = chi0(G1,G2,io) + (rhotwg(G1)*CONJG(rhotwg(G2)))*green_w(io)
    ! * The non-analytic term is symmetrized for this k-point in the BZ according to:
    !    rhotwg(1)= S^-1q * rhotwx_ibz
    !    rhotwg(1)=-S^-1q * CONJG(rhotwx_ibz) if time-reversal is used.

    ! Multiply elements G1,G2 of rhotwg_sym by green_w(io) and accumulate in chi0(G1,G2,io)
    if (gw_eet==-1) then
      do io=1,Ep%nomega
        dd=green_w(io) 
        call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg,1,rhotwg,1,chi0(:,:,io),Ep%npwe)
      end do
    endif

    ! === Accumulate heads and wings for each small q ===
    ! FIXME extrapolar method should be checked!!

    ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
    mir_kbz =(3-2*itim_kbz) * MATMUL(Cryst%symrec(:,:,isym_kbz),rhotwx(:,1)) 
    if (itim_kbz==2) mir_kbz=CONJG(mir_kbz)

    ! here we might take advantage of Hermiticity along Im axis in RPA (see mkG0w)
    do idir=1,3
      do io=1,Ep%nomega                                                       
        chi0_uwing(:,io,idir) = chi0_uwing(:,io,idir) + green_w(io) * mir_kbz(idir)*CONJG(rhotwg)
        chi0_lwing(:,io,idir) = chi0_lwing(:,io,idir) + green_w(io) * rhotwg*CONJG(mir_kbz(idir)) 
        ! Add contribution due to extrapolar technique.
        !if (gwcomp==1.and.ABS(deltaf_b1b2) >= GW_TOL_DOCC) then
        if (gwcomp==1) then
          chi0_uwing(:,io,idir) = chi0_uwing(:,io,idir) + green_enhigh_w(io) * mir_kbz(idir)*CONJG(rhotwg)
          chi0_lwing(:,io,idir) = chi0_lwing(:,io,idir) + green_enhigh_w(io) * rhotwg*CONJG(mir_kbz(idir)) 
        end if
      end do
    end do
    !
    ! Accumulate the head.
    do io=1,Ep%nomega
      do jdir=1,3
        do idir=1,3
          chi0_head(idir,jdir,io) = chi0_head(idir,jdir,io) + green_w(io) * mir_kbz(idir)*CONJG(mir_kbz(jdir))
          ! Add contribution due to extrapolar technique.
          !if (gwcomp==1.and.ABS(deltaf_b1b2) >= GW_TOL_DOCC) then
          if (gwcomp==1) then
            chi0_head(idir,jdir,io) = chi0_head(idir,jdir,io) + green_enhigh_w(io) * mir_kbz(idir)*CONJG(mir_kbz(jdir))
          end if
        end do
      end do
    end do

  else ! spinorial case 
    ABI_ALLOCATE(rhotwg_I,(Ep%npwe))
    ABI_ALLOCATE(rhotwg_J,(Ep%npwe))

    MSG_ERROR("Check tensor")

    ! I can use symmetries to loop over the upper triangle but 
    ! this makes using BLAS more difficult
    ! Important NOTE: treatment of q-->0 limit is correct only
    ! for i=j=0. Other components require additional terms.

    do jj=1,Ep%nJ
      s_jj=1 ; if (jj==4) s_jj=-1
      pad_jj=(jj-1)*Ep%npwe
      call mkrhotwg_sigma(jj,nspinor,Ep%npwe,rhotwg,rhotwg_J)

      !rhotwg_J(1) = q0limit(jj,qlwl(:,1),nspinor,rhotwx,b1,b2,b3) 
      !TODO RECHECK this
      if (itim_kbz==2) rhotwg_J(1)=-CONJG(rhotwg_J(1))

      do ii=1,Ep%nI
        pad_ii=(ii-1)*Ep%npwe

        if (ii/=jj) then
          call mkrhotwg_sigma(ii,nspinor,Ep%npwe,rhotwg,rhotwg_I)
          !rhotwg_I(1) = q0limit(ii,qlwl(:,1),nspinor,rhotwx,b1,b2,b3) 
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
           ! This operation belongs to the little group and has to be considered to reconstruct the BZ.
           phmGt => Gsph_epsG0%phmGt  (1:Ep%npwe,isym) ! In the 2 lines below note the slicing (1:npwe)
           Sm1G  => Gsph_epsG0%rottbm1(1:Ep%npwe,itim,isym)

           SELECT CASE (itim)

           CASE (1)
             rhotwg_sym(1:Ep%npwe)=rhotwg(Sm1G(1:Ep%npwe))*phmGt(1:Ep%npwe) 

           CASE (2) 
             rhotwg_sym(1:Ep%npwe)=CONJG(rhotwg(Sm1G(1:Ep%npwe)))*phmGt(1:Ep%npwe) 

           CASE DEFAULT 
             write(msg,'(a,i3)')'Wrong value of itim= ',itim
             MSG_BUG(msg)
           END SELECT 

           ! Multiply elements G1,G2 of rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)
           if (gw_eet==-1) then
             do io=1,Ep%nomega
               dd=green_w(io) 
               call XGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
             end do
           endif

           ! === Accumulate heads and wings for each small q ===
           ! FIXME extrapolar method should be checked!!

           ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
           mir_kbz =(3-2*itim) * MATMUL(Cryst%symrec(:,:,isym),rhotwx(:,1)) 
           if (itim==2) mir_kbz=CONJG(mir_kbz)

           ! here we might take advantage of Hermiticity along Im axis in RPA (see mkG0w)
           do idir=1,3
             do io=1,Ep%nomega                                                       
               chi0_uwing(:,io,idir) = chi0_uwing(:,io,idir) + green_w(io) * mir_kbz(idir)*CONJG(rhotwg_sym)
               chi0_lwing(:,io,idir) = chi0_lwing(:,io,idir) + green_w(io) * rhotwg_sym*CONJG(mir_kbz(idir))
               ! Add contribution due to extrapolar technique.
               !if (gwcomp==1.and.ABS(deltaf_b1b2)>=GW_TOL_DOCC) then 
               if (gwcomp==1) then
                 chi0_uwing(:,io,idir) = chi0_uwing(:,io,idir) + green_enhigh_w(io) * mir_kbz(idir)*CONJG(rhotwg_sym)
                 chi0_lwing(:,io,idir) = chi0_lwing(:,io,idir) + green_enhigh_w(io) * rhotwg_sym*CONJG(mir_kbz(idir))
               end if
             end do
           end do
           !
           ! Accumulate the head.
           do io=1,Ep%nomega
             do jdir=1,3
               do idir=1,3
                  chi0_head(idir,jdir,io) = chi0_head(idir,jdir,io) +  green_w(io) * mir_kbz(idir)*CONJG(mir_kbz(jdir))
                  ! Add contribution due to extrapolar technique.
                  !if (gwcomp==1.and.ABS(deltaf_b1b2) >= GW_TOL_DOCC) then 
                  if (gwcomp==1) then
                    chi0_head(idir,jdir,io) = chi0_head(idir,jdir,io) + green_enhigh_w(io)*mir_kbz(idir)*CONJG(mir_kbz(jdir))
                  end if
               end do
             end do
           end do

         end if !wtksym
       end do !itim
     end do !isym
   
     ABI_DEALLOCATE(rhotwg_sym)

   else  !spinorial case
     MSG_ERROR('symchi=1 with spinor not implemented ')
     MSG_ERROR("Check tensor")
   end if

 CASE DEFAULT
   write(msg,'(a,i3)')'Wrong value of symchi= ',Ep%symchi
   MSG_BUG(msg)
 END SELECT

end subroutine accumulate_chi0_q0
!!***

!----------------------------------------------------------------------

!!****f* assemblychi0q0_sym/q0limit
!! NAME
!!  q0limit
!!
!! FUNCTION
!!  Return an appropriate linear combination of the spin-dependent oscilator matrix elements at G=0, taking
!!  into account the small q used for the treatment of the optical limit. Used when nspinor=2.
!!
!! INPUTS
!! ii=Index defining the required linear combination of the oscillator strengths.
!! nspinor=Number of spinorial components.
!! qlwl(3)=Reduced components of the small q around Gamma for the treatment of the long wave-length limit.
!! b1(3),b2(3),b3(3)=Lattice vectore of the reciprocal lattice.
!! rhotwx(3,nspinor**2)=Oscillator matrix elements at G=0, for each possible combination of the left and right spin component.
!!
!! OUTPUT
!!  q0limit=Linear combination, see doc below.
!!
!! TODO 
!!  This function should be "contained" to facilitate inlining but abilint crashes, dont know why!
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function q0limit(ii,qlwl,nspinor,rhotwx,b1,b2,b3)

 use defs_basis
 use m_errors

 use m_gwdefs, only : j_gw

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'q0limit'
 use interfaces_70_gw, except_this_one => q0limit
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,nspinor
 complex(gwpc) :: q0limit
!arrays
 real(dp),intent(in) :: qlwl(3)
 real(dp),intent(in) :: b1(3),b2(3),b3(3)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *********************************************************************

 SELECT CASE (ii)
 CASE (1) ! M_0(q-->0) = Lim M(up,up)+M(dwn,dwn). Exact, neglecting Vnl
   q0limit =  dotproductqrc(qlwl,rhotwx(:,1),b1,b2,b3) &
&            +dotproductqrc(qlwl,rhotwx(:,2),b1,b2,b3) 

 CASE (2) ! M_z(q-->0) = Lim M(up,up)-M(dwn,dwn). 
   ! WARNING off-diagonal elements of rV12 and rV12 are neglected
   q0limit =  dotproductqrc(qlwl,rhotwx(:,1),b1,b2,b3) &
&            -dotproductqrc(qlwl,rhotwx(:,2),b1,b2,b3) 

 CASE (3) ! M_x(q-->0) = M(up,dwn)+M(dwn,up). 
   ! Both diagonal elements of the form v12r-rv21 and similiar terms in 12 and 21 are neglected
   q0limit =  dotproductqrc(qlwl,rhotwx(:,3),b1,b2,b3) &
&            +dotproductqrc(qlwl,rhotwx(:,4),b1,b2,b3) 

 CASE (4)
   ! Both diagonal elements of the form v12r-rv21 and similiar terms in 12 and 21 are neglected
   q0limit =( dotproductqrc(qlwl,rhotwx(:,3),b1,b2,b3) &
&            -dotproductqrc(qlwl,rhotwx(:,4),b1,b2,b3) )*j_gw

 CASE DEFAULT
   write(msg,'(a,i3)')'Wrong value for ii= ',ii
   MSG_BUG(msg)
 END SELECT

end function q0limit
!!***
