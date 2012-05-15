!{\src2tex{textfont=tt}}
!!****f* ABINIT/accumulate_sfchi0_q0
!! NAME
!! accumulate_sfchi0_q0
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
!!  sf_chi0(npwe,npwe,my_wl:my_wr)=Updated spectral function at q==0.
!!  sf_lwing(npwe,my_wl:my_wr,3)=Updated lower wing of the spectral function.
!!  sf_uwing(npwe,mw_wl:my_wr,3)=Updated upper wing of the spectral function.
!!  sf_head(3,3,my_wl:my_wr)=Updated head of the spectral function.
!!
!! PARENTS
!!      cchi0q0
!!
!! CHILDREN
!!      xgerc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine accumulate_sfchi0_q0(ikbz,isym_kbz,itim_kbz,nspinor,symchi,npwepG0,npwe,Cryst,Ltg_q,Gsph_epsG0,&
& factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx,rhotwg,nomegasf,sf_chi0,sf_head,sf_lwing,sf_uwing)

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
#define ABI_FUNC 'accumulate_sfchi0_q0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikbz,my_wl,my_wr,nomegasf,npwe,npwepG0,nspinor
 integer,intent(in) :: isym_kbz,itim_kbz,symchi,iomegal,iomegar
 real(dp),intent(in) :: factocc,wl,wr 
 type(Little_group),intent(in) :: Ltg_q
 type(Gvectors_type),intent(in) :: Gsph_epsG0 
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 complex(gwpc),intent(inout) :: sf_chi0(npwe,npwe,my_wl:my_wr)
 complex(dpc),intent(inout) :: sf_head(3,3,my_wl:my_wr)
 complex(dpc),intent(inout) :: sf_lwing(npwe,my_wl:my_wr,3)
 complex(dpc),intent(inout) :: sf_uwing(npwe,my_wl:my_wr,3)

!Local variables-------------------------------
!scalars
 integer :: itim,isym,idir,jdir
 complex(gwpc) :: num
 character(len=500) :: msg
!arrays
 integer,pointer :: Sm1G(:) 
 complex(dpc) :: mir_kbz(3)
 complex(gwpc),pointer :: phmGt(:)
 complex(gwpc),allocatable :: rhotwg_sym(:)

!************************************************************************

 if (iomegal<my_wl .or. iomegar>my_wr) then 
   write(msg,'(3a,2(a,i0,a,i0,a))')ch10,&
&    ' accumulate_sfchi0_q0 : Indeces out of boundary ',ch10,&
&    '  my_wl = ',my_wl,' iomegal = ',iomegal,ch10,&
&    '  my_wr = ',my_wr,' iomegar = ',iomegar,ch10
   MSG_BUG(msg)
 end if 

 SELECT CASE (symchi)

 CASE (0)
   ! 
   ! === Calculation without symmetries ===
   ! * rhotwg(1)= R^-1q*rhotwx_ibz
   ! * rhotwg(1)=-R^-1q*conjg(rhotwx_ibz) for inversion
   if (nspinor==1) then

     if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
       num=-wl*factocc ! Num is single precision needed for cgerc check factocc
       call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,sf_chi0(:,:,iomegal),npwe)
     end if 
     ! Last point, must accumulate left point but not the right one
     if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then 
       num=-wr*factocc
       call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,sf_chi0(:,:,iomegar),npwe)
     end if 
     !
     ! === Accumulate heads and wings for each small q ===
     ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
     mir_kbz =(3-2*itim_kbz) * MATMUL(Cryst%symrec(:,:,isym_kbz),rhotwx(:,1)) 
     if (itim_kbz==2) mir_kbz=CONJG(mir_kbz)
     !
     ! ================================
     ! ==== Update heads and wings ====
     ! ================================
     do jdir=1,3
       !
       if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding    
         num=-wl*factocc ! Num is single precision needed for cgerc check factocc
         sf_uwing(:,iomegal,jdir) = sf_uwing(:,iomegal,jdir) + num * mir_kbz(jdir) * CONJG(rhotwg(1:npwepG0))
         sf_lwing(:,iomegal,jdir) = sf_lwing(:,iomegal,jdir) + num * rhotwg(1:npwepG0) * CONJG(mir_kbz(jdir))
         do idir=1,3
           sf_head(idir,jdir,iomegal) = sf_head(idir,jdir,iomegal) + num * mir_kbz(idir) * CONJG(mir_kbz(jdir))
         end do
       end if 
       !
       ! Last point, must accumulate left point but not the right one
       if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then 
         num=-wr*factocc
         sf_uwing(:,iomegar,jdir) = sf_uwing(:,iomegar,jdir) + num * mir_kbz(jdir) * CONJG(rhotwg(1:npwepG0))
         sf_lwing(:,iomegar,jdir) = sf_lwing(:,iomegar,jdir) + num * rhotwg(1:npwepG0) * CONJG(mir_kbz(jdir))
         do idir=1,3
           sf_head(idir,jdir,iomegar) = sf_head(idir,jdir,iomegar) + num * mir_kbz(idir) * CONJG(mir_kbz(jdir))
         end do
       end if 
     end do ! jdir

   else ! spinorial case
     MSG_BUG("Spectral method + nspinor==2 not implemented")
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
           !
           phmGt => Gsph_epsG0%phmGt(1:npwe,isym) ! In these 2 lines mind the slicing (1:npwe)
           Sm1G  => Gsph_epsG0%rottbm1(1:npwe,itim,isym)

           SELECT CASE (itim)

           CASE (1)
             rhotwg_sym(1:npwe)=rhotwg(Sm1G(1:npwe))*phmGt(1:npwe)

           CASE (2) 
             rhotwg_sym(1:npwe)=CONJG(rhotwg(Sm1G(1:npwe)))*phmGt(1:npwe)

           CASE DEFAULT
             write(msg,'(a,i0)')'Wrong value of itim= ',itim
             MSG_BUG(msg)
           END SELECT
           !
           ! === Multiply elements G,Gp of rhotwg_sym*num and accumulate in sf_chi0(G,Gp,io) ===
           if (wl<huge(0.0_dp)*1.d-11) then
             num=-wl*factocc
             call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,sf_chi0(:,:,iomegal),npwe)
           end if
           !
           ! Last point, must accumulate left point but not the right one
           if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then 
             num=-wr*factocc
             call XGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,sf_chi0(:,:,iomegar),npwe)
           end if 

           ! === Accumulate heads and wings for each small q ===
           ! Symmetrize <r> in full BZ: <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
           mir_kbz =(3-2*itim) * MATMUL(Cryst%symrec(:,:,isym),rhotwx(:,1)) 
           if (itim==2) mir_kbz=CONJG(mir_kbz)
 
           do jdir=1,3
             !
             if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding    
               num=-wl*factocc ! Num is single precision needed for cgerc check factocc
               sf_uwing(:,iomegal,jdir) = sf_uwing(:,iomegal,jdir) + num * mir_kbz(jdir) * CONJG(rhotwg_sym(1:npwe))
               sf_lwing(:,iomegal,jdir) = sf_lwing(:,iomegal,jdir) + num * rhotwg_sym(1:npwe) * CONJG(mir_kbz(jdir))
               do idir=1,3
                 sf_head(idir,jdir,iomegal) = sf_head(idir,jdir,iomegal) + num * mir_kbz(idir) * CONJG(mir_kbz(jdir))
               end do
             end if 
             !
             ! Last point, must accumulate left point but not the right one
             if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then 
               num=-wr*factocc
               sf_uwing(:,iomegar,jdir) = sf_uwing(:,iomegar,jdir) + num * mir_kbz(jdir) * CONJG(rhotwg_sym(1:npwe))
               sf_lwing(:,iomegar,jdir) = sf_lwing(:,iomegar,jdir) + num * rhotwg_sym(1:npwe) * CONJG(mir_kbz(jdir))
               do idir=1,3
                 sf_head(idir,jdir,iomegar) = sf_head(idir,jdir,iomegar) + num * mir_kbz(idir) * CONJG(mir_kbz(jdir))
               end do
             end if 
           end do ! jdir

         end if !wtksym
       end do !inv
     end do !isym
     ABI_DEALLOCATE(rhotwg_sym)

   else ! spinorial case
     MSG_BUG("Spectral method + nspinor==2 not implemented")
   end if

 CASE DEFAULT
   write(msg,'(a,i0)')'Wrong value of symchi= ',symchi
   MSG_BUG(msg)
 END SELECT

end subroutine accumulate_sfchi0_q0
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/assemblychi0sf
!! NAME
!! assemblychi0sf
!!
!! FUNCTION
!! Update the spectral function of the irreducible polarizability for the contribution
!! of one pair of occupied-unoccupied states, for each frequenciy.
!! If symchi==1, the symmetries of the little group of the external q-point are used 
!! to symmetrize the contribution in the full Brillouin zone. In this case, the routine computes: 
!! 
!!   $ chi0(G1,G2,io)=chi0(G1,G2,io)+\sum_S (rhotwg(G1)*rhotwg^\dagger(G2))*\delta(w - trans) $
!! 
!! where S are the symmetries of the little group of the external q-point.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ik_bz=Index of the k-point in the BZ whose contribution has to be added to the spectral function of chi0 
!!    If symchi=1, the contribution is symmetrized.
!!  my_wl,my_wr=min and Max frequency index treated by this processor.
!!  npwe=Number of plane waves used to describe chi0.
!!  npwepG0=Maximum number of G vectors taking into account umklapp vectors.
!!  nomegasf=Number of frequencies for the spectral function.
!!  nspinor=Number of spinorial components.
!!  symchi=1 if symmetries are used, 0 otherwise
!!  rhotwg(npwepG0*nspinor**2)=Oscillator matrix elements corresponding to an occupied-unoccupied pair of states.
!!  timrev=if 2, time reversal has to be used to obtain k_bz; 1 otherwise.
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
!!  factocc=occupation factor=f_occ*(ockp-occk) (see cchi0.F90)  
!!  wl,wr=Weights used to approximate the delta function.
!!    
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0sf(npwe,npwe,my_wl:my_wr)= updated spectral function.
!!
!! NOTES
!!  Umklapp processes are not yet implemented 
!! 
!! PARENTS
!!      cchi0
!!
!! CHILDREN
!!      xgerc
!!
!! SOURCE

subroutine assemblychi0sf(ik_bz,nspinor,symchi,Ltg_q,npwepG0,npwe,rhotwg,Gsph_epsG0,&
& factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,nomegasf,chi0sf)

 use m_profiling

 use defs_basis
 use m_errors

 use m_gwdefs,   only : epsilonm1_parameters
 use m_blas,     only : xgerc
 use m_gsphere,  only : gvectors_type
 use m_bz_mesh,  only : little_group

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'assemblychi0sf'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,iomegal,iomegar,my_wl,my_wr,nomegasf,npwe,npwepG0
 integer,intent(in) :: nspinor,symchi
 real(dp),intent(in) :: factocc,wl,wr
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Little_group),intent(in) :: Ltg_q
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(inout) :: chi0sf(npwe,npwe,my_wl:my_wr)

!Local variables-------------------------------
!scalars
 integer :: isym,itim
 complex(gwpc) :: num
 character(len=500) :: msg
!arrays
 integer,allocatable :: Sm1_gmG0(:)
 integer,pointer :: gmG0(:)
 complex(gwpc),allocatable :: rhotwg_sym(:)
 complex(gwpc),pointer :: phmGt(:)

! *************************************************************************

 if (iomegal < my_wl .or. iomegar > my_wr) then 
   write(msg,'(3a,2(a,i0,a,i0,a))')ch10,&
&    ' Indeces out of boundary ',ch10,&
&    '  my_wl = ',my_wl,' iomegal = ',iomegal,ch10,&
&    '  my_wr = ',my_wr,' iomegar = ',iomegar,ch10
   MSG_BUG(msg)
 end if 

 SELECT CASE (symchi)

 CASE (0)  ! Do not use symmetries ===
   if (wl<huge(0.0_dp)*1.d-11) then !FIXME this is awful
     num=-wl*factocc 
     call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegal),npwe)
   end if
   !
   ! Last point, must accumulate left point but not the right one
   if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then 
     num=-wr*factocc
     call XGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegar),npwe)
   end if 

 CASE (1)
   ! Use symmetries to reconstruct oscillator matrix elements
   ! Notes on the symmetrization of the oscillator maxtri elements:
   ! 
   ! If  Sq=q then  M_G^( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1G}  (k,q)
   ! If -Sq=q then  M_G^(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1G}^*(k,q)
   ! 
   ! In case of an umklapp process 
   ! If  Sq=q+G_o then  M_G( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1(G-G_o}   (k,q)
   ! If -Sq=q+G_o then  M_G(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1(G-G_o)}^*(k,q)
   !
   ! Ltg_q%igmG0(ig,itim,isym) contains the index of G-G0 where ISq=q+G0
   ! Note that there is no need to take into account the phases due to q, 
   ! They cancel in the scalar product ==> phmGt(G,isym)=e^{-iG\cdot t}
   ! 
   ! Mind the slicing of %rottbm1(npwepG0,timrev,nsym) and %phgt(npwepG0,nsym) as 
   ! these arrays, usually, do not conform to rho_twg_sym(npw) !
   ! 
   ABI_ALLOCATE(rhotwg_sym,(npwe))
   ABI_ALLOCATE(Sm1_gmG0,(npwe))
   !
   ! === Loop over symmetries of the space group and time-reversal ===
   do isym=1,Ltg_q%nsym_sg
     do itim=1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then 
         ! === This operation belongs to the little group and has to be used to reconstruct BZ ===
         ! TODO this is a hot-spot, should add a test on the umklapp
         !
         phmGt => Gsph_epsG0%phmGt(1:npwe,isym) ! In these 3 lines mind the slicing (1:npwe)
         gmG0 => Ltg_q%igmG0(1:npwe,itim,isym) 
         Sm1_gmG0(1:npwe)=Gsph_epsG0%rottbm1(gmG0(1:npwe),itim,isym)

         SELECT CASE (itim)
         CASE (1)
           !rhotwg_sym(1:npwe)=rhotwg(Sm1_gmG0)*Gsph_epsG0%phmGt(1:npwe,isym)
           rhotwg_sym(1:npwe)=rhotwg(Sm1_gmG0(1:npwe))*phmGt(1:npwe)
         CASE (2) 
           rhotwg_sym(1:npwe)=CONJG(rhotwg(Sm1_gmG0(1:npwe)))*phmGt(1:npwe)
         CASE DEFAULT
           write(msg,'(a,i0)')'Wrong value for itim= ',itim
           MSG_BUG(msg)
         END SELECT

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
       end if !wtksym 

     end do !inv
   end do !isym
   ABI_DEALLOCATE(rhotwg_sym)
   ABI_DEALLOCATE(Sm1_gmG0)

 CASE DEFAULT
   write(msg,'(a,i0)')'Wrong value for symchi= ',symchi
   MSG_BUG(msg)
 END SELECT

end subroutine assemblychi0sf
!!***
