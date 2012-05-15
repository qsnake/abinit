!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_delta_ppm
!! NAME
!! calc_delta_ppm
!!
!! FUNCTION
!! Calculation of the function \delta required for the EET
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calc_delta_ppm(Sigp,ptwsq,niter)

 use m_profiling

 use defs_basis
 use m_errors
 use m_gwdefs

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_delta_ppm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Sigma_parameters),intent(in) :: Sigp

 integer,intent(in) :: niter
!arrays
 complex(gwpc),intent(inout) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)

!Local variables-------------------------------
!scalars
 integer :: ig,igp
!arrays
 real(dp) :: delta_huge,denchk
 complex(gwpc) :: num,den

!*************************************************************************

 delta_huge = 1.0d8

 if(niter==1) then

   do igp=1,Sigp%npwc
     do ig=igp,Sigp%npwc
       num = ptwsq(ig,igp,2)+conjg(ptwsq(igp,ig,2))
       den = ptwsq(ig,igp,1)+conjg(ptwsq(igp,ig,1))
       denchk=abs(real(den))+abs(aimag(den))
       if (denchk>0.0) then
         ptwsq(ig,igp,2)=num/den
       else
         ptwsq(ig,igp,2)=cmplx(delta_huge,0.0)
       endif
     end do !ig
   end do !igp

   do igp=1,Sigp%npwc
     do ig=1,igp-1
       ptwsq(ig,igp,2)=conjg(ptwsq(igp,ig,2))
     enddo
   enddo

 elseif(niter==2) then

   do igp=1,Sigp%npwc
     do ig=igp,Sigp%npwc
       num = ptwsq(ig,igp,3)+conjg(ptwsq(igp,ig,3))
       den = ptwsq(ig,igp,2)+conjg(ptwsq(igp,ig,2))
       denchk=abs(real(den))+abs(aimag(den))
       if (denchk>0.0) then
         ptwsq(ig,igp,3)=num/den
       else
         ptwsq(ig,igp,3)=cmplx(delta_huge,0.0)
       endif
       num = ptwsq(ig,igp,2)+conjg(ptwsq(igp,ig,2))
       den = ptwsq(ig,igp,1)+conjg(ptwsq(igp,ig,1))
       denchk=abs(real(den))+abs(aimag(den))
       if (denchk>0.0) then
         ptwsq(ig,igp,2)=num/den
       else
         ptwsq(ig,igp,2)=cmplx(delta_huge,0.0)
       endif
     end do !ig
   end do !igp

   do igp=1,Sigp%npwc
     do ig=1,igp-1
       ptwsq(ig,igp,2)=conjg(ptwsq(igp,ig,2))
       ptwsq(ig,igp,3)=conjg(ptwsq(igp,ig,3))
     enddo
   enddo

 endif

end subroutine calc_delta_ppm
!!***

!!****f* ABINIT/calc_delta_ppm_sc
!! NAME
!! calc_delta_ppm_sc
!!
!! FUNCTION
!! Calculation of the function \delta required for OPTIMAL GW
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calc_delta_ppm_sc(Sigp,nomega,otq,omegame0k,omegame0lumo,npwc2, &
&                          qbzpg,ikbz,jkbz,ptwsq,delta,niter)

 use m_profiling

 use defs_basis
 use m_errors
 use m_gwdefs

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_delta_ppm_sc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Sigma_parameters),intent(in) :: Sigp

 integer,intent(in) :: npwc2,nomega,niter
!arrays
 real(dp),intent(in) :: omegame0k(nomega),omegame0lumo(nomega)
 complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
 complex(gwpc),intent(in) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
 complex(gwpc),intent(out) :: delta(Sigp%npwc,Sigp%npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ios,ikbz,jkbz
 real(dp) :: otw,omegakg
 complex(gwpc) :: daux(niter)
!arrays
 real(dp),intent(in) :: qbzpg(Sigp%npwc)
 real(dp) :: bq

 integer :: iter

 real(dp) :: delta_huge,denchk,test

 logical :: ggpnonzero

 complex(gwpc) :: dfrac(niter)
 complex(gwpc) :: num,den

!*************************************************************************

 delta_huge = 1.0d8

 if (niter==0) then

   do ig = 1, Sigp%npwc
     do igp = ig, Sigp%npwc
       delta(ig,igp,:) = 0.25*(qbzpg(ig)*qbzpg(ig)+qbzpg(igp)*qbzpg(igp))
     enddo
   enddo
   do ig = 1, Sigp%npwc
     do igp = 1, ig-1
       delta(ig,igp,:)=conjg(delta(igp,ig,:))
      enddo
   enddo

 elseif (niter==1) then

   do ig=1,Sigp%npwc
     do igp=1,Sigp%npwc
       ggpnonzero=(ig/=1.and.igp/=1)
       if (ikbz/=jkbz.or.ggpnonzero) then
         num = ptwsq(ig,igp,2)
         den = ptwsq(ig,igp,1)
         denchk=abs(real(den))+abs(aimag(den))
         if (denchk>0.0) then
           dfrac(1)=num/den
         else
           dfrac(1)=cmplx(delta_huge,0.0)
         endif
         bq=0.25*(qbzpg(ig)*qbzpg(ig)+qbzpg(igp)*qbzpg(igp))
         delta(ig,igp,1) = bq + dfrac(1)
       else
         delta(ig,igp,1) = (0.0,0.0)
       endif
     end do !igp
   end do !ig

   do ios = 2, nomega
     delta(:,:,ios) = delta(:,:,1)
   enddo

 elseif (niter==2) then

   do ig=1,Sigp%npwc
     do igp=1,Sigp%npwc
       ggpnonzero=(ig/=1.and.igp/=1)
       if (ikbz/=jkbz.or.ggpnonzero) then
         do iter = 1, 2
           num = ptwsq(ig,igp,iter+1)+conjg(ptwsq(igp,ig,iter+1))
           den = ptwsq(ig,igp,iter)+conjg(ptwsq(igp,ig,iter))
           denchk=abs(real(den))+abs(aimag(den))
           if (denchk>0.0) then
             dfrac(iter)=num/den
           else
             dfrac(iter)=cmplx(delta_huge,0.0)
           endif
         enddo
         bq=0.25*(qbzpg(ig)*qbzpg(ig)+qbzpg(igp)*qbzpg(igp))
         otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
         do ios=1,nomega
           omegakg=omegame0k(ios)-otw
           do iter = 1, 2
             daux(iter) = omegakg-bq-dfrac(iter)
           enddo
           delta(ig,igp,ios) = bq + dfrac(1)*daux(1)/daux(2)
         enddo
       else
         delta(ig,igp,:) = (0.0,0.0)
       endif
     end do !igp
   end do !ig

 endif

 do ig = 1, Sigp%npwc
   do igp = 1, Sigp%npwc
     do ios = 1, nomega
      ggpnonzero=(ig/=1.and.igp/=1)
      if (ikbz/=jkbz.or.ggpnonzero) then
        if (ig==igp) then
          delta(ig,igp,ios)=real(delta(ig,igp,ios))
          test=omegame0k(ios)-real(delta(ig,igp,ios))
          if (test>omegame0lumo(ios).and.niter>0) then
            delta(ig,igp,ios)=0.25*qbzpg(ig)**2+0.25*qbzpg(igp)**2
            test=omegame0k(ios)-real(delta(ig,igp,ios))
          endif
          if (test>omegame0lumo(ios)) then
            delta(ig,igp,ios)=omegame0k(ios)-omegame0lumo(ios)
          endif
        else
          test=omegame0k(ios)-real(delta(ig,igp,ios))
          if (test>omegame0lumo(ios)) then
            delta(ig,igp,ios)=omegame0k(ios)-omegame0lumo(ios)
          endif
        endif
      else
        delta(ig,igp,ios)=(0.0,0.0)
      endif
     enddo
   enddo
 enddo

end subroutine calc_delta_ppm_sc
!!***

!!****f* ABINIT/calc_sig_ppm_delta
!! NAME
!! calc_sig_ppm_delta
!!
!! FUNCTION
!! Calculation of the part of the matrix elements of the self-energy coming from the sum
!! over the valence states in OPTIMAL GW
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine calc_sig_ppm_delta(npwc,nomega,rhotwgp,botsq,otq,omegame0i,zcut,theta_mu_minus_e0i, &
&                             ket,npwx,npwc1,npwc2,omega4sd,e0,delta)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_sig_ppm_delta'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,npwc1,npwc2,npwx
 real(dp),intent(in) :: theta_mu_minus_e0i,zcut
!arrays
 real(dp),intent(in) :: omegame0i(nomega),e0
 complex(dpc),intent(in) :: omega4sd(nomega)
 complex(gwpc),intent(in) :: otq(npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(npwc,npwc1)
 complex(gwpc),intent(in) :: rhotwgp(npwx)
 complex(gwpc),intent(in) :: delta(npwc,npwc,nomega)
 complex(gwpc),intent(inout) :: ket(npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ios
 real(dp) :: den,omegame0i_io,otw,twofm1_zcut
 complex(gwpc) :: num,den2,rhotwgdp_igp
 logical :: fully_occupied,totally_empty

!*************************************************************************

  fully_occupied=(abs(theta_mu_minus_e0i-1.)<0.001)
  totally_empty=(abs(theta_mu_minus_e0i)<0.001)

  if (.not.(totally_empty)) then
   twofm1_zcut=zcut

   do ios=1,nomega
    omegame0i_io=omegame0i(ios)
    do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
      otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
      num = botsq(ig,igp)*rhotwgdp_igp
      den = omegame0i_io+otw

      if (den**2>zcut**2)then
       ket(ig,ios) = ket(ig,ios) + num/(den*otw) * theta_mu_minus_e0i
      else
       ket(ig,ios) = ket(ig,ios) + num*cmplx(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw)&
&                            *theta_mu_minus_e0i
      end if
     end do !ig
    end do !igp
   end do !ios

  end if !not totally empty


  if (.not.(totally_empty)) then
   twofm1_zcut=-zcut

   do ios=1,nomega
    omegame0i_io=DBLE(omega4sd(ios)) - e0
    do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
      otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
      num = botsq(ig,igp)*rhotwgdp_igp
      den2 = omegame0i_io-otw-delta(ig,igp,ios)

      if (real(conjg(den2)*den2)>zcut**2) then
       ket(ig,ios) = ket(ig,ios) - num/(den2*otw)*theta_mu_minus_e0i
      else
       ket(ig,ios) = ket(ig,ios) - num*(den2+cmplx(0.0,twofm1_zcut))/((den2**2+twofm1_zcut**2)*otw) &
&                           *theta_mu_minus_e0i
      end if
     end do !ig
    end do !igp
   end do !ios

  end if

  if (.not.(fully_occupied)) then
   twofm1_zcut=-zcut

   do ios=1,nomega
    omegame0i_io=DBLE(omega4sd(ios)) - e0
    do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
      otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
      num = botsq(ig,igp)*rhotwgdp_igp
      den2 = omegame0i_io-otw-delta(ig,igp,ios)

      if (real(conjg(den2)*den2)>zcut**2) then
       ket(ig,ios) = ket(ig,ios) - num/(den2*otw)*(1.-theta_mu_minus_e0i)
      else
       ket(ig,ios) = ket(ig,ios) - num*(den2+cmplx(0.0,twofm1_zcut))/((den2**2+twofm1_zcut**2)*otw) &
&                           *(1.-theta_mu_minus_e0i)
      end if
     end do !ig
    end do !igp
   end do !ios

   do ios=1,nomega
    omegame0i_io=omegame0i(ios)
    do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
      otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
      num = botsq(ig,igp)*rhotwgdp_igp

      den = omegame0i_io-otw
      if (den**2>zcut**2) then
       ket(ig,ios) = ket(ig,ios) + num/(den*otw)*(1.-theta_mu_minus_e0i)
      else
       ket(ig,ios) = ket(ig,ios) + num*cmplx(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw) &
&                           *(1.-theta_mu_minus_e0i)
      end if
     end do !ig
    end do !igp
   end do !ios

  end if

  ket(:,:)=ket(:,:)*0.5

end subroutine calc_sig_ppm_delta
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/calc_sig_ppm_delta_clos
!! NAME
!! calc_sig_ppm_delta_clos
!!
!! FUNCTION
!! Calculation of the part of the matrix elements of the self-energy coming from
!! the sum over all the states in OPTIMAL GW
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE
subroutine calc_sig_ppm_delta_clos(npwc,nomega,ikbz,jkbz,qbzpg,botsq,otq,omegame0k,omegame0lumo,zcut,ptwsq, &
&                                   ket,npwc1,npwc2,gw_eet_scale,niter)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_sig_ppm_delta_clos'
 use interfaces_70_gw, except_this_one => calc_sig_ppm_delta_clos
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,npwc1,npwc2,niter
 integer,intent(in) :: ikbz,jkbz
 real(dp),intent(in) :: zcut,gw_eet_scale
!arrays
 real(dp),intent(in) :: omegame0k(nomega),omegame0lumo(nomega)
 real(dp),intent(in) :: qbzpg(npwc)
 complex(gwpc),intent(in) :: otq(npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(npwc,npwc1)
 complex(gwpc),intent(in) :: ptwsq(npwc,npwc,niter+1)
 complex(dpc),intent(inout) :: ket(nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ios
 real(dp) :: omegame0i_io,otw,omegakg,twofm1_zcut
 complex(gwpc) :: den,num

 real(dp), allocatable :: qpgsq(:)
 real(dp) :: bq
 complex(gwpc) :: delta,deltaux

 logical :: ggpnonzero

!*************************************************************************

   twofm1_zcut=-zcut

   ABI_ALLOCATE(qpgsq,(npwc))

   do ig = 1, npwc
     qpgsq(ig)=half*qbzpg(ig)*qbzpg(ig)
   enddo

   do igp=1,npwc
     do ig=1,npwc
       otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
       ggpnonzero=(ig/=1.and.igp/=1)
       if (ikbz/=jkbz.or.ggpnonzero) then
         bq=half*(qpgsq(ig)+qpgsq(igp))
         if (niter==0) then
           delta = bq
         elseif(niter==1) then
           delta = bq+ptwsq(ig,igp,2)
         endif
         do ios=1,nomega
           if(niter==2) then
             omegakg=omegame0k(ios)-otw
             delta = bq + ptwsq(ig,igp,2)*(omegakg-bq-ptwsq(ig,igp,2))/(omegakg-bq-ptwsq(ig,igp,3))
           endif
           if (niter<2) deltaux=delta
           call check_delta_sigma(qpgsq(ig),qpgsq(igp),delta,omegame0k(ios),omegame0lumo(ios),ig,igp,gw_eet_scale,niter)
           omegame0i_io=omegame0k(ios)
           num = botsq(ig,igp)*ptwsq(ig,igp,1)
           den = omegame0i_io-otw-delta
           if (real(conjg(den)*den)>zcut**2) then
             ket(ios) = ket(ios) + 0.5*num/(den*otw)
           else
             ket(ios) = ket(ios) + 0.5*num*(den+cmplx(0.0,twofm1_zcut))/((den**2+twofm1_zcut**2)*otw)
           end if
           if (niter<2) delta=deltaux
         end do !ios
       else
         do ios=1,nomega
           omegame0i_io=omegame0k(ios)
           num = botsq(ig,igp)*ptwsq(ig,igp,1)
           den = omegame0i_io-otw
           if (real(conjg(den)*den)>zcut**2) then
             ket(ios) = ket(ios) + 0.5*num/(den*otw)
           else
             ket(ios) = ket(ios) + 0.5*num*(den+cmplx(0.0,twofm1_zcut))/((den**2+twofm1_zcut**2)*otw)
           end if
         end do !ios
       endif
     end do !ig
   end do !igp

 ABI_DEALLOCATE(qpgsq)

end subroutine calc_sig_ppm_delta_clos
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/gw_eet_sigma
!! NAME
!! gw_eet_sigma
!!
!! FUNCTION
!! Wrapper routine for the calculation of the matrix elements of Sigma using GW_OPTIMAL
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine gw_eet_sigma(Sigp,Sr,Dtset,Cryst,Wfs,Kmesh,Qmesh,Gsph_Max,Gsph_c,Psps,Vcp,QP_BSt,PPm, &
&                       isppol,iq_bz,ik_bz,jk_bz,ik_ibz,jk_ibz,itim_q,isym_q,iq_ibz,tabr_ki, &
&                       tabr_kj,spinrot_ki,spinrot_kj,ph_mkit,ph_mkjt,nfftot_gw,ngfft_gw, &
&                       use_padfft,igfftcg0,gw_gbound,gw_mgfft,ib1,ib2,nomega_tot,nomega_sigc, &
&                       fact_sp,nspinor,botsq,otq,sigcme_tmp,sigc,nbhomo,tim_fourdp,wtqp,wtqm, &
&                       MPI_enreg,extrapolar_distrb,can_symmetrize)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_errors
 use m_ppmodel

 use m_geometry,      only : normv, vdotw
 use m_vcoul,         only : vcoul_t
 use m_sigma_results, only : sigma_results

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_eet_sigma'
 use interfaces_70_gw, except_this_one => gw_eet_sigma
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh,Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(Gvectors_type),intent(in) :: Gsph_Max
 type(Gvectors_type),intent(in) :: Gsph_c
 type(Dataset_type),intent(in) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(Sigma_parameters),intent(in) :: Sigp
 type(Sigma_results),intent(in) :: Sr
 type(wfs_descriptor),intent(inout) :: Wfs
 type(Bandstructure_type),intent(in) :: QP_BSt
 type(PPmodel_type),intent(inout) :: PPm
 type(MPI_type),intent(inout) :: MPI_enreg

 integer,intent(in) :: nfftot_gw,ngfft_gw(18),tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: igfftcg0(Sigp%npwc)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)

 integer,intent(in) :: isppol
 integer,intent(in) :: itim_q,isym_q,iq_ibz
 integer,intent(in) :: nomega_tot,nomega_sigc
 integer,intent(in) :: nspinor
 integer,intent(in) :: tabr_ki(nfftot_gw),tabr_kj(nfftot_gw)
 integer,intent(in) :: iq_bz,ik_bz,jk_bz,ik_ibz,jk_ibz,ib1,ib2
 integer,intent(in) :: extrapolar_distrb(ib1:ib2,ib1:ib2,Kmesh%nbz,Wfs%nsppol)
 integer,intent(in) :: wtqp,wtqm
 integer,intent(out) :: nbhomo
 real(dp),intent(in) :: spinrot_ki(4),spinrot_kj(4)
 complex(dpc),intent(in) ::  ph_mkit,ph_mkjt
 complex(dpc),intent(inout) :: sigcme_tmp(nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nsppol*Sigp%nsig_ab)
 complex(dpc),intent(inout) :: sigc(2,nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nsppol*Sigp%nsig_ab)
 logical,intent(in) :: can_symmetrize(Wfs%nsppol)
 real(dp),intent(in) :: fact_sp
 complex(gwpc),intent(in) :: otq(Sigp%npwc,PPm%dm2_otq)
 complex(gwpc),intent(in) :: botsq(Sigp%npwc,PPm%dm2_botsq)
 real(dp),pointer :: qp_ene(:,:,:),qp_occ(:,:,:)

 real(dp),allocatable :: qplg(:,:)
 real(dp),allocatable :: kplqg(:)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:)
 complex(gwpc),allocatable :: fnlkr(:,:,:)
 complex(gwpc),allocatable :: fnlkpr(:,:,:)
 complex(gwpc),allocatable :: wfr1(:,:)
 complex(gwpc),allocatable :: mtwk(:,:)
 complex(gwpc),allocatable :: mtwkp(:,:)
 complex(gwpc),allocatable :: ptwsq(:,:,:)

 integer :: io,kb,jb
 integer :: niter,nptwg,iter,nbmax
 integer :: ig,igp,ib,ibv
 integer :: isym_kgw,isym_ki,iik,jik

 real(dp),allocatable :: omegame0k(:),omegame0lumo(:)

 real(dp),allocatable :: qbzpg(:) 

 complex(gwpc),allocatable :: otq_transp(:,:)
 complex(gwpc),allocatable :: botsq_conjg_transp(:,:)

 integer :: fnlloc(Cryst%ntypat,2),fnlmax(Cryst%ntypat)

 complex(dpc) :: sigctmp(nomega_sigc)
 complex(dpc) :: sigctmp2(nomega_sigc)

!************************************************************************

 qp_ene => QP_BSt%eig(:,:,:)
 qp_occ => QP_BSt%occ(:,:,:)

 !@Arjan:
 ! Gsph_c gamma-centered sphere for W and thus Sigma_x
 ! Gsph_Max gamma-centered sphere with gvec(3,npwvec) where Sigp%npwvec=MAX(Sigp%npwwfn,Sigp%npwx)
 ! it is used for the wavefunctions but it will be REMOVED when we switch to k-centered
 ! G-spheres for the wavefunctions

 ! min and Max band indeces for GW corrections (for this k-point)

 nbhomo=1
 do ib = 2, Sigp%nbnds 
   if (fact_sp*qp_occ(ib,jk_ibz,isppol)<GW_TOL_DOCC) exit
   nbhomo=nbhomo+1
 enddo
 nbmax=max(nbhomo,Dtset%gw_eet_nband)

 !allocate(val_idx(QP_BSt%nkpt,QP_BSt%nsppol))
 !val_idx = get_valence_idx(QP_BSt,tol3)
 !do isppol=1,nsppol
 ! nbhomo(isppol) = val_idx(1,isppol)
 ! ltest = ALL(val_idx(:,isppol))==val_idx(1,isppol),
 ! ABI_CHECK(ltest,"Optimal GW + metals not coded")
 !end do
 !deallocate(val_idx)

 niter = Dtset%gw_eet
 nptwg=1
 do iter = 1, niter
   nptwg=nptwg+mod(iter,2)
 enddo

 ABI_ALLOCATE(vc_sqrt_qbz,(Sigp%npwc))
 do ig=1,Sigp%npwc
   vc_sqrt_qbz(Gsph_c%rottb(ig,itim_q,isym_q))=Vcp%vc_sqrt(ig,iq_ibz)
 end do

 ABI_ALLOCATE(omegame0k,(nomega_tot))
 ABI_ALLOCATE(omegame0lumo,(nomega_tot))
 ABI_ALLOCATE(qplg,(Sigp%npwc,3))
 ABI_ALLOCATE(kplqg,(Sigp%npwc))
 ABI_ALLOCATE(qbzpg,(Sigp%npwc))

 isym_kgw = Kmesh%tabo(jk_bz)
 jik = (3-Kmesh%tabi(jk_bz))/2

 isym_ki = Kmesh%tabo(ik_bz)
 iik = (3-Kmesh%tabi(ik_bz))/2

 do ig=1,Sigp%npwc
   qplg(ig,:) = Qmesh%bz(:,iq_bz) + Gsph_c%gvec(:,ig)
   kplqg(ig) = -vdotw(Kmesh%bz(:,jk_bz),qplg(ig,:),Cryst%gmet,"G")
   qbzpg(ig) = normv(qplg(ig,:),Cryst%gmet,"G")
 end do

 if (niter>0.and.Dtset%gw_eet_inclvkb==1) then

   ABI_ALLOCATE(fnlkr,(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fnlkpr,(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(mtwk,(Wfs%nfftot*nspinor,nbmax))
   ABI_ALLOCATE(mtwkp,(Wfs%nfftot*nspinor,ib1:ib2))

   call gw_eet_sigma_vkb(Sigp,Cryst,Wfs,Kmesh,Psps,isppol,ik_ibz,jk_ibz,ib1,ib2,nspinor, &
&                        tim_fourdp,MPI_enreg,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)
 endif

 ABI_ALLOCATE(wfr1,(Wfs%nfftot*nspinor,nbmax))

 do ibv = 1, nbmax
   call wfd_get_ur(Wfs,ibv,ik_ibz,isppol,wfr1(:,ibv))
 enddo

 do jb = ib1,ib2

   do kb = ib1,ib2

     if (Sigp%gwcalctyp/=28.and.kb/=jb) CYCLE

     if (extrapolar_distrb(jb,kb,ik_bz,isppol)/=Wfs%my_rank) CYCLE

     ABI_ALLOCATE(ptwsq,(Sigp%npwc,Sigp%npwc,niter+1))

     if (Sigp%gwcalctyp==28) then
       do io=1,Sr%nomega4sd
         omegame0k(io)=real(Sr%omega4sd(kb,jk_ibz,io,isppol))-half*(qp_ene(kb,jk_ibz,isppol)+qp_ene(jb,jk_ibz,isppol))
         omegame0lumo(io)= real(Sr%omega4sd(kb,jk_ibz,io,isppol)) - qp_ene(nbmax+1,ik_ibz,isppol)
       end do
     else
       do io=1,Sr%nomega4sd
         omegame0k(io)  = real(Sr%omega4sd(kb,jk_ibz,io,isppol)) - qp_ene(kb,jk_ibz,isppol)
         omegame0lumo(io)= real(Sr%omega4sd(kb,jk_ibz,io,isppol)) - qp_ene(nbmax+1,ik_ibz,isppol)
       end do
     endif

     sigctmp=czero_gw

     if (Sigp%gwcalctyp==28) then

       call fft4eet_sig_sc(Sigp,Dtset,Cryst,Wfs,Gsph_c,Sr,nbhomo,nbmax,nomega_tot,isppol,nfftot_gw, &
&                   ngfft_gw,use_padfft,igfftcg0,gw_gbound,gw_mgfft,iik,tabr_ki,ph_mkit,spinrot_ki, &
&                   ik_ibz,jk_ibz,isym_kgw,jik,tabr_kj,ph_mkjt,spinrot_kj,Gsph_Max%rottbm1, &
&                   nspinor,tim_fourdp,MPI_enreg,wfr1,vc_sqrt_qbz,Vcp%i_sz,jb,kb,qplg,kplqg,niter, &
&                   ptwsq,ik_bz,jk_bz,PPm%dm2_botsq,PPm%dm2_otq,botsq,otq,sigctmp)

     elseif (niter==0.or.Dtset%gw_eet_inclvkb==0) then

       call fft4eet_sig(Sigp,Dtset,Cryst,Wfs,Gsph_c,Sr,nbhomo,nbmax,nomega_tot,isppol,nfftot_gw, &
&                   ngfft_gw,use_padfft,igfftcg0,gw_gbound,gw_mgfft,iik,tabr_ki,ph_mkit,spinrot_ki, &
&                   ik_ibz,jk_ibz,isym_kgw,jik,tabr_kj,ph_mkjt,spinrot_kj,Gsph_Max%rottbm1, &
&                   nspinor,tim_fourdp,MPI_enreg,wfr1,vc_sqrt_qbz,Vcp%i_sz,kb,qplg,kplqg,niter, &
&                   ptwsq,ik_bz,jk_bz,PPm%dm2_botsq,PPm%dm2_otq,botsq,otq,sigctmp)

     else

       call fft4eet_sig_kb(Sigp,Dtset,Cryst,Wfs,Kmesh,Gsph_c,Psps,Sr,nbhomo,nbmax,nomega_tot,isppol,nfftot_gw, &
&                   ngfft_gw,use_padfft,igfftcg0,gw_gbound,gw_mgfft,iik,tabr_ki,ph_mkit,spinrot_ki, &
&                   ik_ibz,jk_ibz,isym_kgw,jik,tabr_kj,ph_mkjt,spinrot_kj,Gsph_Max%rottbm1, &
&                   nspinor,tim_fourdp,MPI_enreg,fnlloc,fnlmax,fnlkr,mtwk,mtwkp(:,kb),wfr1, &
&                   vc_sqrt_qbz,Vcp%i_sz,kb,qplg,kplqg,niter,ptwsq,ik_bz,jk_bz,PPm%dm2_botsq,PPm%dm2_otq,botsq,otq, &
&                   sigctmp)

     endif

!     if (Sigp%gwcalctyp==28) then
!       call calc_delta_ppm_sc(Sigp,nomega_tot,otq,omegame0k,omegame0lumo,PPm%dm2_otq,qbzpg,ik_bz,jk_bz,ptwsq,delta, &
!&                           niter)
!     else
       call calc_delta_ppm(Sigp,ptwsq,niter)
!     endif

     do ig = 1, Sigp%npwc
       do igp = 1, Sigp%npwc
         if (ik_bz==jk_bz) then
           if (ig/=1.and.igp/=1) then
             ptwsq(ig,igp,1)= ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)*vc_sqrt_qbz(igp)
           else
             if (jb/=kb.or.kb<=nbmax.or.jb<=nbmax) then
               ptwsq(ig,igp,1)=(0.0,0.0)
             else
               if (ig==1.and.igp==1) then
                 ptwsq(ig,igp,1) = cmplx(Vcp%i_sz,0.0_gwp)
               elseif (ig==1.and.igp/=1) then
                 ptwsq(ig,igp,1) = cmplx(sqrt(Vcp%i_sz),0.0_gwp)*ptwsq(ig,igp,1)*vc_sqrt_qbz(igp)
               elseif (igp==1.and.ig/=1) then
                 ptwsq(ig,igp,1) = cmplx(sqrt(Vcp%i_sz),0.0_gwp)*ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)
               else
                 ptwsq(ig,igp,1)= ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)*vc_sqrt_qbz(igp)
               endif
             endif
           endif
         else
           ptwsq(ig,igp,1)= ptwsq(ig,igp,1)*vc_sqrt_qbz(ig)*vc_sqrt_qbz(igp)
         endif
       enddo
     enddo

     call calc_sig_ppm_delta_clos(Sigp%npwc,nomega_tot,ik_bz,jk_bz,qbzpg,botsq,otq,omegame0k,omegame0lumo,Sigp%zcut, &
&                                 ptwsq,sigctmp,PPm%dm2_botsq,PPm%dm2_otq,Dtset%gw_eet_scale,niter)

     if (Sigp%gwcalctyp==28) then
       ABI_ALLOCATE(botsq_conjg_transp,(PPm%dm2_botsq,Sigp%npwc))
       botsq_conjg_transp=TRANSPOSE(botsq) ! Keep these two lines separated, otherwise gfortran messes up
       botsq_conjg_transp=CONJG(botsq_conjg_transp)
       ABI_ALLOCATE(otq_transp,(PPm%dm2_otq,PPm%npwc))
       otq_transp=TRANSPOSE(otq)

       sigctmp2=czero_gw
       call calc_sig_ppm_delta_clos(Sigp%npwc,nomega_tot,ik_bz,jk_bz,qbzpg,botsq_conjg_transp,otq_transp, &
&                                   omegame0k,omegame0lumo,Sigp%zcut,ptwsq,sigctmp2,PPm%dm2_botsq,PPm%dm2_otq, &
&                                   Dtset%gw_eet_scale,niter)
       ABI_DEALLOCATE(botsq_conjg_transp)
       ABI_DEALLOCATE(otq_transp)
       sigctmp=half*(sigctmp+sigctmp2)
     endif
   
     if (can_symmetrize(isppol)) then
       sigcme_tmp(:,jb,kb,isppol)=sigcme_tmp(:,jb,kb,isppol) + &
         (wtqp+wtqm)*DBLE(sigctmp(:)) + (wtqp-wtqm)*j_gw*AIMAG(sigctmp(:))
       sigc(1,:,jb,kb,isppol)=sigc(1,:,jb,kb,isppol) + wtqp*      sigctmp(:)
       sigc(2,:,jb,kb,isppol)=sigc(2,:,jb,kb,isppol) + wtqm*CONJG(sigctmp(:))
     else
       sigcme_tmp(:,jb,kb,isppol)=sigcme_tmp(:,jb,kb,isppol)+sigctmp(:)
     endif

     ABI_DEALLOCATE(ptwsq)

   enddo
 enddo

 ABI_DEALLOCATE(wfr1)
 if (niter>0.and.Dtset%gw_eet_inclvkb==1) then
   ABI_DEALLOCATE(mtwk)
   ABI_DEALLOCATE(mtwkp)
   ABI_DEALLOCATE(fnlkr)
   ABI_DEALLOCATE(fnlkpr)
 endif
 ABI_DEALLOCATE(vc_sqrt_qbz)
 ABI_DEALLOCATE(omegame0k)
 ABI_DEALLOCATE(omegame0lumo)
 ABI_DEALLOCATE(qplg)
 ABI_DEALLOCATE(kplqg)
 ABI_DEALLOCATE(qbzpg)

end subroutine gw_eet_sigma
!!***
 
!----------------------------------------------------------------------

!!****f* ABINIT/gw_eet_chi0
!! NAME
!! gw_eet_chi0
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine gw_eet_chi0(Ep,Dtset,Cryst,Wfs,Kmesh,Gsph_epsG0,Gsph_wfn,Psps,Ltg_q,nbvw,qpoint, &
&                      nfftot_gw,ngfft_gw,use_padfft,igfftepsG0,gw_gbound,gw_mgfft,is, &
&                      ik_bz,ik_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                      ikmq_ibz,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg, &
&                      qp_energy,chi0,spin_fact,qp_occ,nspinor,tim_fourdp,bbp_ks_distrb,nbmax)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_errors

 use m_geometry,  only : normv,vdotw
 use m_vcoul,     only : vcoul_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_eet_chi0'
 use interfaces_70_gw, except_this_one => gw_eet_chi0
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Gvectors_type),intent(in) :: Gsph_epsG0,Gsph_wfn
 type(Dataset_type),intent(in) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfs_descriptor),intent(inout) :: Wfs
 type(Little_group),intent(in) :: Ltg_q

 integer,intent(in) :: nbvw,nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: is
 integer,intent(in) :: isym_k,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Ep%npwepG0)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: dim_rtwg
 integer,intent(in) :: bbp_ks_distrb(nbvw,Kmesh%nbz,Wfs%nsppol)
 integer,intent(out) :: nbmax
 real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4),spin_fact

 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)

 real(dp),intent(in) :: qpoint(3)

 integer,intent(in) :: ik_bz,ik_ibz,ikmq_ibz

 integer :: ib,ibv
 integer :: niter,nptwg,iter
 integer :: ig

 real(dp),allocatable :: qpgsq(:)
 real(dp),allocatable :: qplg(:,:)
 real(dp),allocatable :: kplqg(:)
 complex(gwpc),allocatable :: fnlkr(:,:,:)
 complex(gwpc),allocatable :: fnlkpr(:,:,:)
 complex(gwpc),allocatable :: wfr1(:,:)
 complex(gwpc),allocatable :: wfr2(:)
 complex(gwpc),allocatable :: mtwk(:,:)
 complex(gwpc),allocatable :: mtwkp(:,:)

 complex(gwpc),allocatable :: frhorho(:)
 complex(gwpc),allocatable :: frhoj(:,:)
 complex(gwpc),allocatable :: fjj(:)

 integer :: nbhomo(2)
 integer :: istat

 integer :: fnlloc(Cryst%ntypat,2),fnlmax(Cryst%ntypat)

!************************************************************************

 nbhomo(:)=1
 do ib = 2, Ep%nbnds
   if (spin_fact*qp_occ(ib,ikmq_ibz,is)<GW_TOL_DOCC) exit
   nbhomo(1)=nbhomo(1)+1
 enddo
 do ib = 2, Ep%nbnds
   if (spin_fact*qp_occ(ib,ikmq_ibz,is)<(one-GW_TOL_DOCC)) exit
   nbhomo(2)=nbhomo(2)+1
 enddo
 nbmax=max(nbhomo(1),Dtset%gw_eet_nband)

 niter = Dtset%gw_eet
 nptwg=1
 do iter = 1, niter
   nptwg=nptwg+mod(iter,2)
 enddo

 ABI_ALLOCATE(qpgsq,(Ep%npwe))
 ABI_ALLOCATE(qplg,(Ep%npwe,3))
 ABI_ALLOCATE(kplqg,(Ep%npwe))

 do ig=1,Ep%npwe
   qplg(ig,:) = qpoint(:)+ Gsph_epsG0%gvec(:,ig)
   kplqg(ig)=-vdotw(Kmesh%bz(:,ik_bz),qplg(ig,:),Cryst%gmet,"G")
   qpgsq(ig) = normv(qplg(ig,:),Cryst%gmet,"G")
 enddo
 qpgsq(:)=half*qpgsq(:)**2

 if (niter>0.and.Dtset%gw_eet_inclvkb==1) then
   ABI_ALLOCATE(fnlkpr,(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fnlkr,(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(mtwk,(Wfs%nfftot*nspinor,nbhomo(1)))
   ABI_ALLOCATE(mtwkp,(Wfs%nfftot*nspinor,nbmax))

   call gw_eet_chi0_vkb(Ep,Cryst,Wfs,Kmesh,Psps,is,ik_ibz,ikmq_ibz,nspinor,tim_fourdp, &
&                       nbhomo,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)
 endif

 ABI_ALLOCATE(wfr1,(Wfs%nfftot*nspinor,nbmax))

 do ibv = 1, nbmax
   call wfd_get_ur(Wfs,ibv,ikmq_ibz,is,wfr1(:,ibv))
 enddo

 do ibv = 1, nbhomo(1)

   if ((bbp_ks_distrb(ibv,ik_bz,is) /= Wfs%my_rank)) CYCLE

   ABI_ALLOCATE(frhorho,(Ep%npwe*(Ep%npwe+1)/2))
   istat = ABI_ALLOC_STAT
   if(istat/=0) stop 'out of memory: gw_eet_chi0: frhorho'
   if (niter>0) then
     ABI_ALLOCATE(frhoj,(Ep%npwe,Ep%npwe))
     istat = ABI_ALLOC_STAT
     if(istat/=0) stop 'out of memory: gw_eet_chi0: frhoj'
     ABI_ALLOCATE(fjj,(Ep%npwe*(Ep%npwe+1)/2*(niter-1)))
     istat = ABI_ALLOC_STAT
     if(istat/=0) stop 'out of memory: gw_eet_chi0: fjj'
   endif

   ABI_ALLOCATE(wfr2,(Wfs%nfftot*nspinor))
   call wfd_get_ur(Wfs,ibv,ik_ibz,is,wfr2)

   if (niter==0) then
     call fft4eet_0(ik_bz,Ep,Wfs,Kmesh,Gsph_epsG0,Ltg_q,nbhomo,nbmax,is,nfftot_gw,ngfft_gw, &
&                   use_padfft,igfftepsG0,gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,itim_k,tabr_k, &
&                   ph_mkt,spinrot_k,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg, &
&                   nspinor,tim_fourdp,wfr1,wfr2,ibv,frhorho,spin_fact,qp_occ,qp_energy,chi0,1)
   else
     if (Dtset%gw_eet_inclvkb==1) then
       call fft4eet_kb(ik_bz,Ep,Cryst,Wfs,Kmesh,Gsph_epsG0,Ltg_q,Psps,nbhomo,nbmax, &
&                      is,nfftot_gw,ngfft_gw,use_padfft,igfftepsG0, &
&                      gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                      itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,Gsph_wfn%rottbm1, &
&                      nspinor,tim_fourdp,fnlloc,fnlmax,fnlkpr,mtwk,mtwkp,wfr1,wfr2,ibv,qplg,kplqg, &
&                      niter,frhorho,frhoj,fjj,spin_fact,qp_occ,qp_energy,chi0,1)
     else
       call fft4eet(ik_bz,Ep,Cryst,Wfs,Kmesh,Gsph_epsG0,Ltg_q,nbhomo,nbmax, &
&                   is,nfftot_gw,ngfft_gw,use_padfft,igfftepsG0, &
&                   gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                   itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,Gsph_wfn%rottbm1, &
&                   nspinor,tim_fourdp,wfr1,wfr2,ibv,qplg,kplqg,niter,frhorho,frhoj,fjj, &
&                   spin_fact,qp_occ,qp_energy,chi0,1)
     endif
   endif

   ABI_DEALLOCATE(wfr2)

   if (niter==0) then
     frhorho=spin_fact*qp_occ(ibv,ik_ibz,is)*frhorho
     call calc_chi0_delta0(ik_bz,Dtset,Ep,Gsph_epsG0,Ltg_q,qp_energy(ibv,ik_ibz,is), &
&                          qp_energy(nbmax+1,ikmq_ibz,is),qpgsq,frhorho,chi0)
   else
     call calc_delta_chi0(Ep,frhorho,frhoj,fjj,niter)
     frhorho=spin_fact*qp_occ(ibv,ik_ibz,is)*frhorho
     call calc_chi0_delta_clos(ik_bz,Dtset,Ep,Gsph_epsG0,Ltg_q,niter,qp_energy(ibv,ik_ibz,is), &
&                              qp_energy(nbmax+1,ikmq_ibz,is),qpgsq,frhorho,frhoj,fjj,chi0)
     ABI_DEALLOCATE(frhoj)
     ABI_DEALLOCATE(fjj)
   endif

   ABI_DEALLOCATE(frhorho)

 enddo

 ABI_DEALLOCATE(wfr1)
 if (niter>0.and.Dtset%gw_eet_inclvkb==1) then
   ABI_DEALLOCATE(mtwk)
   ABI_DEALLOCATE(mtwkp)
   ABI_DEALLOCATE(fnlkr)
   ABI_DEALLOCATE(fnlkpr)
 endif

 ABI_DEALLOCATE(qpgsq)
 ABI_DEALLOCATE(qplg)
 ABI_DEALLOCATE(kplqg)

 nbmax=0

end subroutine gw_eet_chi0

!!***
 
!----------------------------------------------------------------------

!!****f* ABINIT/drho_tw_g
!! NAME                  
!! drho_tw_g
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine drho_tw_g(paral_kgb,nspinor,npwvec,nr,ngfft,map2sphere,use_padfft,igfftg0,gbound,&
&                    wfn1,i1,ktabr1,ktabp1,wfn2,dim_rtwg,rhotwg,tim_fourdp,MPI_enreg)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_gwdefs,    only : czero_gw

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'drho_tw_g'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: paral_kgb,i1,npwvec,nr,tim_fourdp,nspinor,dim_rtwg,map2sphere,use_padfft
 complex(dpc),intent(in) :: ktabp1
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 !integer,intent(in) :: gbound(2*mgfft+8,2)
 integer,intent(in) :: gbound(:,:)
 integer,intent(in) :: igfftg0(npwvec*map2sphere),ngfft(18)
 integer,intent(in) :: ktabr1(nr)
 complex(gwpc),intent(in) :: wfn1(nr*nspinor),wfn2(nr*nspinor)
 complex(gwpc),intent(out) :: rhotwg(npwvec*dim_rtwg)

!Local variables-------------------------------
!scalars
 integer :: ig,igfft
 integer :: nx,ny,nz,ldx,ldy,ldz,mgfft
!arrays
 complex(dpc),allocatable :: usk(:),uu(:)

! *************************************************************************

 SELECT CASE (nspinor)

 CASE (1) ! Collinear case.
  !
  ! Form rho-twiddle(r)=u_1^*(r,b1,kbz1) u_2(r,b2,kbz2), to account for
  ! symmetries:
  ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz)
  !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
  !
  ABI_ALLOCATE(uu,(nr))
  ABI_ALLOCATE(usk,(nr))

  uu  = wfn1(ktabr1)*ktabp1; if (i1==1) uu  = CONJG(uu)
!  usk = wfn2(ktabr2)*ktabp2; if (i2==2) usk = CONJG(usk)
  usk = wfn2
  uu  = uu * usk

  SELECT CASE (map2sphere)

  CASE (0) ! Need results on the full FFT box thus cannot use zero-padded FFT.

    call fourdp_c2c_ip(uu,-1,MPI_enreg,nr,ngfft,paral_kgb,tim_fourdp)
    rhotwg=uu
    !call fourdp_c2c_op(uu,rhotwg,-1,MPI_enreg,nr,ngfft,paral_kgb,tim_fourdp)

  CASE (1) ! Need results on the G-sphere. Call zero-padded FFT routines if required.

    if (use_padfft==1) then
      nx =ngfft(1); ny =ngfft(2); nz =ngfft(3); mgfft = MAXVAL(ngfft(1:3))
      ldx=nx      ; ldy=ny      ; ldz=nz
      call padded_fourwf_cplx(uu,ngfft,nx,ny,nz,ldx,ldy,ldz,mgfft,-1,gbound)
    else
      call fourdp_c2c_ip(uu,-1,MPI_enreg,nr,ngfft,paral_kgb,tim_fourdp)
    end if

    do ig=1,npwvec       ! Have to map FFT to G-sphere.
      igfft=igfftg0(ig)
      if (igfft/=0) then ! G-G0 belong to the FFT mesh.
        rhotwg(ig)=uu(igfft)
      else               ! Set this component to zero.
        rhotwg(ig)=czero_gw
      end if
    end do

  CASE DEFAULT
    MSG_BUG("Wrong map2sphere")
  END SELECT

  ABI_DEALLOCATE(uu)
  ABI_DEALLOCATE(usk)

  RETURN

 CASE DEFAULT
   MSG_BUG('Wrong nspinor')
 END SELECT

end subroutine drho_tw_g

!!***
 
!----------------------------------------------------------------------

!!****f* ABINIT/calc_dwfwfg
!! NAME                  
!! calc_dwfwfg
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine calc_dwfwfg(MPI_enreg,paral_kgb,tim_fourdp,ktabr_k,ktabi_k,nfftot,ngfft_gw,ph_mkt,wfr_jb,wfr_kb,wfg2_jk)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_dwfwfg'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ktabi_k,nfftot,paral_kgb,tim_fourdp
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: ktabr_k(nfftot),ngfft_gw(18)
 complex(dpc),intent(in) :: ph_mkt
 complex(gwpc),intent(in) :: wfr_jb(nfftot),wfr_kb(nfftot)
 complex(gwpc),intent(out) :: wfg2_jk(nfftot)

!Local variables-------------------------------
!arrays
 complex(dpc),allocatable :: wfr2_dpcplx(:)
#if ! defined HAVE_GW_DPC
 complex(dpc),allocatable :: wfg2_dpcplx(:)
#endif

! *************************************************************************

 ABI_ALLOCATE(wfr2_dpcplx,(nfftot))

 SELECT CASE (ktabi_k)

 CASE (1)
!   wfr2_dpcplx = CONJG(wfr_jb(ktabr_k)) * wfr_kb(ktabr_k)
   wfr2_dpcplx = CONJG(ph_mkt*wfr_jb(ktabr_k)) * wfr_kb

 CASE (2) ! Conjugate the product if time-reversal is used to reconstruct this k-point
!   wfr2_dpcplx = wfr_jb(ktabr_k) * CONJG(wfr_kb(ktabr_k))
   wfr2_dpcplx = ph_mkt*wfr_jb(ktabr_k) * wfr_kb

 CASE DEFAULT
   MSG_ERROR("Wrong ktabi_k")
 END SELECT

 ! Transform to Fourier space (result in wfg2_jk)
#if defined HAVE_GW_DPC
 call fourdp_c2c_op(wfr2_dpcplx,wfg2_jk,-1,MPI_enreg,nfftot,ngfft_gw,paral_kgb,tim_fourdp)
#else
 ABI_ALLOCATE(wfg2_dpcplx,(nfftot))
 call fourdp_c2c_op(wfr2_dpcplx,wfg2_dpcplx,-1,MPI_enreg,nfftot,ngfft_gw,paral_kgb,tim_fourdp)
 wfg2_jk=wfg2_dpcplx
 ABI_DEALLOCATE(wfg2_dpcplx)
#endif

 ABI_DEALLOCATE(wfr2_dpcplx)

end subroutine calc_dwfwfg

!!***
 
!----------------------------------------------------------------------

!!****f* ABINIT/calc_ddwfwfg
!! NAME                  
!! calc_ddwfwfg
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine calc_ddwfwfg(MPI_enreg,paral_kgb,tim_fourdp,ktabi_k,nfftot,ngfft_gw,wfr_jb,wfr_kb,wfg2_jk)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_ddwfwfg'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ktabi_k,nfftot,paral_kgb,tim_fourdp
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: ngfft_gw(18)
 complex(gwpc),intent(in) :: wfr_jb(nfftot),wfr_kb(nfftot)
 complex(gwpc),intent(out) :: wfg2_jk(nfftot)

!Local variables-------------------------------
!arrays
 complex(dpc),allocatable :: wfr2_dpcplx(:)
#if ! defined HAVE_GW_DPC
 complex(dpc),allocatable :: wfg2_dpcplx(:)
#endif

! *************************************************************************

 ! There is no need to take into account phases arising from non-symmorphic
 ! operations since the wavefunctions are evaluated at the same k-point.
 ABI_ALLOCATE(wfr2_dpcplx,(nfftot))

 SELECT CASE (ktabi_k)

 CASE (1)
!   wfr2_dpcplx = CONJG(wfr_jb(ktabr_k)) * wfr_kb(ktabr_k)
   wfr2_dpcplx = CONJG(wfr_jb) * wfr_kb

 CASE (2) ! Conjugate the product if time-reversal is used to reconstruct this k-point
!   wfr2_dpcplx = wfr_jb(ktabr_k) * CONJG(wfr_kb(ktabr_k))
   wfr2_dpcplx = CONJG(wfr_jb) * wfr_kb

 CASE DEFAULT
   MSG_ERROR("Wrong ktabi_k")
 END SELECT

 ! Transform to Fourier space (result in wfg2_jk)
#if defined HAVE_GW_DPC
 call fourdp_c2c_op(wfr2_dpcplx,wfg2_jk,-1,MPI_enreg,nfftot,ngfft_gw,paral_kgb,tim_fourdp)
#else
 ABI_ALLOCATE(wfg2_dpcplx,(nfftot))
 call fourdp_c2c_op(wfr2_dpcplx,wfg2_dpcplx,-1,MPI_enreg,nfftot,ngfft_gw,paral_kgb,tim_fourdp)
 wfg2_jk=wfg2_dpcplx
 ABI_DEALLOCATE(wfg2_dpcplx)
#endif

 ABI_DEALLOCATE(wfr2_dpcplx)

end subroutine calc_ddwfwfg

!!***
 
!----------------------------------------------------------------------

!!****f* ABINIT/fft4eet_0
!! NAME                  
!! fft4eet_0
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine fft4eet_0(ik_bz,Ep,Wfs,Kmesh,Gsph_epsG0,Ltg_q,nbhomo,nbmax,is,nfftot_gw,ngfft_gw, &
&                    use_padfft,igfftepsG0,gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,itim_k,tabr_k,ph_mkt, &
&                    spinrot_k,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,nspinor,tim_fourdp, &
&                    wfr1,wfr2,ibv,frhorho,spin_fact,qp_occ,qp_energy,chi0,igstart)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_errors

 use m_geometry,      only : vdotw
 use m_vcoul,         only : vcoul_t
 use m_oscillators,   only : rho_tw_g, calc_wfwfg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_0'
 use interfaces_14_hidewrite
 use interfaces_70_gw, except_this_one => fft4eet_0
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfs_descriptor),intent(inout) :: Wfs
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Little_group),intent(in) :: Ltg_q

 integer,intent(in) :: ik_bz,nbhomo(2),nbmax,is,ibv,igstart
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_ibz,ikmq_ibz
 integer,intent(in) :: itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: dim_rtwg
 integer,intent(in) :: igfftepsG0(Ep%npwepG0)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: spin_fact
 real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt
 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: wfr2(Wfs%nfftot*nspinor)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

 complex(gwpc),intent(out) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)

 integer,pointer :: gbound_k(:,:),kg_k(:,:),igfft_k(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: wfwfg(:)

 integer :: ibvp,ig,igp,istwf_k,npw_k
 integer :: ig4,ig4x,ig4y,ig4z,igaux
 integer :: gmgp(3)
 integer :: outofbox
 integer,save :: enough=0
 character(len=500) :: msg

!************************************************************************

 ABI_ALLOCATE(rhotwg,(Ep%npwepG0*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 istwf_k  =  Wfs%istwfk(ik_ibz)
 npw_k    =  Wfs%npwarr(ik_ibz)
 gbound_k => Wfs%Kdata(ik_ibz)%gbound
 kg_k     => Wfs%Kdata(ik_ibz)%kg_k    
 igfft_k  => Wfs%Kdata(ik_ibz)%igfft0

 call calc_wfwfg(Wfs%MPI_enreg,Wfs%paral_kgb,tim_fourdp,tabr_k,itim_k,nfftot_gw,ngfft_gw, &
&                wfr2,wfr2,wfwfg(:))

 do ibvp = 1, nbmax
   call rho_tw_g(Wfs%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,wfr2,itim_k,tabr_k,ph_mkt, &
&                spinrot_k,dim_rtwg,rhotwg(:,ibvp),tim_fourdp,Wfs%MPI_enreg)

   if (ibvp>nbhomo(2)) then
     call calc_corr_chi0(ik_bz,Ep,Kmesh,Gsph_epsG0,Ltg_q,rhotwg(:,ibvp),spin_fact,qp_occ,qp_energy, &
&                        ibv,ibvp,ik_ibz,ikmq_ibz,is,chi0)
   endif

 enddo

 frhorho(:)=czero_gw
 outofbox=0
 do ig=igstart,Ep%npwe
   do igp=ig,Ep%npwe
     gmgp(:)=kg_k(:,ig)-kg_k(:,igp)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)
     igaux=igp*(igp-1)/2+ig
     frhorho(igaux)=wfwfg(ig4)
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if

 do ibvp = 1, nbmax
   do ig=igstart,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       frhorho(igaux)=frhorho(igaux)-conjg(rhotwg(igp,ibvp))*rhotwg(ig,ibvp)
     enddo
   end do
 end do

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)

end subroutine fft4eet_0

!!***

!----------------------------------------------------------------------

!!****f* ABINIT/fft4eet
!! NAME                  
!! fft4eet
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine fft4eet(ik_bz,Ep,Cryst,Wfs,Kmesh,Gsph_epsG0,Ltg_q,nbhomo,nbmax, &
&                  is,nfftot_gw,ngfft_gw,use_padfft,igfftepsG0, &
&                  gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                  itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,grottbm1, &
&                  nspinor,tim_fourdp,wfr1,wfr2,ibv,qplg,kplqg, &
&                  niter,frhorho,frhoj,fjj,spin_fact,qp_occ,qp_energy,chi0,igstart)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_errors

 use m_geometry,      only : vdotw
 use m_vcoul,         only : vcoul_t
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_blas,          only : xgerc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet'
 use interfaces_14_hidewrite
 use interfaces_70_gw, except_this_one => fft4eet
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfs_descriptor),intent(inout) :: Wfs
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Little_group),intent(in) :: Ltg_q

 integer,intent(in) :: ik_bz,nbhomo(2),nbmax,is,ibv,niter,igstart
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_ibz,ikmq_ibz
 integer,intent(in) :: isym_k,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: dim_rtwg
 integer,intent(in) :: igfftepsG0(Ep%npwepG0)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: grottbm1(Ep%npwvec,2,Cryst%nsym)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: qplg(Ep%npwe,3),kplqg(Ep%npwe)
 real(dp),intent(in) :: spin_fact
 real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt
 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: wfr2(Wfs%nfftot*nspinor)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

 complex(gwpc),intent(out) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(out) :: frhoj(Ep%npwe,Ep%npwe)
 complex(gwpc),intent(out) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))

 integer,pointer :: gbound_k(:,:),kg_k(:,:),igfft_k(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: drhotwg(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)

 complex(gwpc),allocatable :: dwfr(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg(:,:)

 integer :: ibvp,ig,igp,igbz,igaux,istwf_k,npw_k
 integer :: i,j
 integer :: ig4,ig4x,ig4y,ig4z
 integer :: gmgp(3)
 integer :: outofbox
 integer,save :: enough=0
 character(len=500) :: msg

 complex(gwpc) :: faux,minusone
 complex(dpc) :: ph_Gt
 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux(3)

!************************************************************************

 ABI_ALLOCATE(rhotwg,(Ep%npwepG0*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))
 if (niter>0) then
   ABI_ALLOCATE(drhotwg,(Ep%npwepG0*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(Ep%npwwfn))
   ABI_ALLOCATE(dwfr,(Wfs%nfftot*nspinor,3))
 endif
 if (niter>1)  then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
 end if

 istwf_k  =  Wfs%istwfk(ik_ibz)
 npw_k    =  Wfs%npwarr(ik_ibz)
 gbound_k => Wfs%Kdata(ik_ibz)%gbound
 kg_k     => Wfs%Kdata(ik_ibz)%kg_k    
 igfft_k  => Wfs%Kdata(ik_ibz)%igfft0

 call calc_wfwfg(Wfs%MPI_enreg,Wfs%paral_kgb,tim_fourdp,tabr_k,itim_k,nfftot_gw,ngfft_gw, &
&                wfr2,wfr2,wfwfg(:))

 do ibvp = 1, nbmax
   call rho_tw_g(Wfs%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,wfr2,itim_k,tabr_k,ph_mkt, &
&                spinrot_k,dim_rtwg,rhotwg(:,ibvp),tim_fourdp,Wfs%MPI_enreg)

   if (ibvp>nbhomo(2)) then
     call calc_corr_chi0(ik_bz,Ep,Kmesh,Gsph_epsG0,Ltg_q,rhotwg(:,ibvp),spin_fact,qp_occ,qp_energy, &
&                        ibv,ibvp,ik_ibz,ikmq_ibz,is,chi0)
   endif

 enddo

 if (niter>0) then

   do i = 1, 3
     do ig = 1, npw_k
       igbz = grottbm1(ig,itim_k,isym_k)
       ph_Gt = EXP(-(0.,1.)*two_pi*DOT_PRODUCT(kg_k(:,ig),Cryst%tnons(:,isym_k)))
       if (itim_k==1) then
         gwfg(ig)=-kg_k(i,ig)*Wfs%Wave(ibv,ik_ibz,is)%ug(igbz)*ph_mkt*ph_Gt
       else
         gwfg(ig)=-kg_k(i,ig)*conjg(Wfs%Wave(ibv,ik_ibz,is)%ug(igbz)*ph_mkt)*ph_Gt
       endif
     enddo

     call fft_onewfn(Wfs%paral_kgb,istwf_k,nspinor,npw_k,Wfs%nfftot,Wfs%mgfft,Wfs%ngfft,&
&                    gwfg,dwfr(:,i),igfft_k,kg_k,gbound_k,tim_fourdp,Wfs%MPI_enreg)
   enddo

   do ibvp = 1, nbmax
     do i = 1, 3
       call drho_tw_g(Wfs%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,dwfr(:,i), &
&                     dim_rtwg,drhotwg(:,ibvp,i),tim_fourdp,Wfs%MPI_enreg)
     enddo
   enddo

   do i = 1, 3
     call calc_dwfwfg(Wfs%MPI_enreg,Wfs%paral_kgb,tim_fourdp,tabr_k,itim_k,nfftot_gw,ngfft_gw,ph_mkt, &
&                     wfr2,dwfr(:,i),dwfwfg(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(Wfs%MPI_enreg,Wfs%paral_kgb,tim_fourdp,itim_k,nfftot_gw,ngfft_gw, &
&                          dwfr(:,i),dwfr(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo

   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr)
   ABI_ALLOCATE(cauxg,(Ep%npwe,nbmax))

   do ibvp = 1, nbmax
     do ig=1,Ep%npwe
       drhaux(:)=cmplx(real(drhotwg(ig,ibvp,:)),aimag(drhotwg(ig,ibvp,:)))
       cauxg(ig,ibvp)=vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
     enddo
   enddo

 endif

 frhorho=czero_gw
 frhoj=czero_gw
 if (niter>1) fjj=czero_gw
 outofbox=0
 do ig=igstart,Ep%npwe
   do igp=igstart,Ep%npwe
     gmgp(:)=kg_k(:,ig)-kg_k(:,igp)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)
     if (ig<=igp) then
       igaux=igp*(igp-1)/2+ig
       frhorho(igaux)=wfwfg(ig4)
     endif
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       frhoj(ig,igp)=frhoj(ig,igp) + vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
       if (niter>1.and.ig<=igp) then
         igaux=igp*(igp-1)/2+ig
         paux(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux(i)=paux(i) + vdotw(qplg(igp,:),drhaux,Cryst%gmet,"G")
         enddo
         fjj(igaux)=fjj(igaux) + vdotw(qplg(ig,:),paux,Cryst%gmet,"G")
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if

 do ibvp = 1, nbmax
   do ig=igstart,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       frhorho(igaux)=frhorho(igaux)-conjg(rhotwg(igp,ibvp))*rhotwg(ig,ibvp)
     enddo
   end do
 end do

 if (niter>0) then
   minusone=(-1.,0.)
   do ibvp = 1, nbmax
     call XGERC(Ep%npwe,Ep%npwe,minusone,cauxg(:,ibvp),1,rhotwg(:,ibvp),1,frhoj,Ep%npwe)
   end do
   do ig=igstart,Ep%npwe
     do igp=igstart,Ep%npwe
       if (ig<=igp) then
         igaux=igp*(igp-1)/2+ig
         faux=frhorho(igaux)
       else
         igaux=ig*(ig-1)/2+igp
         faux=conjg(frhorho(igaux))
       endif
       frhoj(ig,igp)=frhoj(ig,igp)+faux*kplqg(ig)
     end do !igp
   end do !ig
 end if

 if (niter>1) then

   do ibvp = 1, nbmax
     do ig=igstart,Ep%npwe
       do igp=ig,Ep%npwe
         igaux=igp*(igp-1)/2+ig
         fjj(igaux)=fjj(igaux)-conjg(cauxg(igp,ibvp))*cauxg(ig,ibvp)
       enddo
     end do
   end do
   do ig=igstart,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       faux=frhorho(igaux)
       fjj(igaux)=fjj(igaux)+kplqg(igp)*frhoj(ig,igp)+ kplqg(ig)*conjg(frhoj(igp,ig)) - &
&                            kplqg(igp)*faux*kplqg(ig)
     end do !igp
   end do !ig

 endif

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg)
   ABI_DEALLOCATE(drhotwg)
   ABI_DEALLOCATE(dwfwfg)
 endif
 if (niter>1)  then
   ABI_DEALLOCATE(ddwfwfg)
 end if

end subroutine fft4eet

!!***
 
!----------------------------------------------------------------------

!!****f* ABINIT/fft4eet_kb
!! NAME                  
!! fft4eet_kb
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine fft4eet_kb(ik_bz,Ep,Cryst,Wfs,Kmesh,Gsph_epsG0,Ltg_q,Psps,nbhomo,nbmax, &
&                     is,nfftot_gw,ngfft_gw,use_padfft,igfftepsG0, &
&                     gw_gbound,gw_mgfft,ik_ibz,ikmq_ibz,isym_k,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                     itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,dim_rtwg,grottbm1, &
&                     nspinor,tim_fourdp,fnlloc,fnlmax,fnlkpr,mtwk,mtwkp,wfr1,wfr2,ibv,qplg,kplqg, &
&                     niter,frhorho,frhoj,fjj,spin_fact,qp_occ,qp_energy,chi0,igstart)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_errors

 use m_geometry,      only : vdotw
 use m_vcoul,         only : vcoul_t
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_blas,          only : xgerc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_kb'
 use interfaces_14_hidewrite
 use interfaces_70_gw, except_this_one => fft4eet_kb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfs_descriptor),intent(inout) :: Wfs
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(kb_potential) :: KBff_k_ibz
 type(Little_group),intent(in) :: Ltg_q

 integer,intent(in) :: ik_bz,nbhomo(2),nbmax,is,ibv,niter,igstart
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_ibz,ikmq_ibz
 integer,intent(in) :: isym_k,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: dim_rtwg
 integer,intent(in) :: igfftepsG0(Ep%npwepG0)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: grottbm1(Ep%npwvec,2,Cryst%nsym)
 integer,intent(in) :: fnlloc(Cryst%ntypat,2)
 integer,intent(in) :: fnlmax(Cryst%ntypat)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: spin_fact
 real(dp),intent(in) :: qplg(Ep%npwe,3),kplqg(Ep%npwe)
 real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt

 complex(gwpc),intent(in) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(in) :: mtwk(Wfs%nfftot*nspinor,nbhomo(1))
 complex(gwpc),intent(in) :: mtwkp(Wfs%nfftot*nspinor,nbmax)

 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: wfr2(Wfs%nfftot*nspinor)

 complex(gwpc),intent(out) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(out) :: frhoj(Ep%npwe,Ep%npwe)
 complex(gwpc),intent(out) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

 integer,pointer :: gbound_k(:,:),kg_k(:,:),igfft_k(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: drhotwg(:,:,:)
 complex(gwpc),allocatable :: fnltwg(:,:,:)
 complex(gwpc),allocatable :: fnltwg2(:,:)
 complex(gwpc),allocatable :: fnltwg3(:,:)
 complex(gwpc),allocatable :: kns(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)
 complex(gwpc),allocatable :: fnlwfg(:)
 complex(gwpc),allocatable :: fkdwfg(:,:)
 complex(gwpc),allocatable :: fdrhotwg(:,:,:,:)
 complex(gwpc),allocatable :: lnkp(:)

 complex(gwpc),allocatable :: dwfr(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg(:,:)
 complex(gwpc),allocatable :: ff(:,:,:,:)
 complex(gwpc),allocatable :: vzn(:,:,:)

 complex(gwpc),allocatable :: paux(:,:)

 integer :: istwf_k,npw_k
 integer :: ilm,iat,ilm2,iat2,ibvp,ig,igp,igbz,igaux
 integer :: i,j
 integer :: ig4,ig4x,ig4y,ig4z
 integer :: ig5,ig5x,ig5y,ig5z
 integer :: nlx
 integer :: gmgp(3),ngfft(3)
 integer :: outofbox
 integer,save :: enough=0
 character(len=500) :: msg

 complex(gwpc) :: faux,minusone
 complex(dpc) :: ph_Gt
 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux2(3)

!************************************************************************

 nlx = min(Psps%mpsang,4)

 ABI_ALLOCATE(rhotwg,(Ep%npwepG0*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 if (niter>0) then
   ABI_ALLOCATE(drhotwg,(Ep%npwepG0*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(Ep%npwwfn))
   ABI_ALLOCATE(dwfr,(Wfs%nfftot*nspinor,3))
   ABI_ALLOCATE(fnltwg,(Ep%npwepG0*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fnlwfg,(nfftot_gw*nspinor**2))
   ABI_ALLOCATE(fnltwg2,(Ep%npwepG0*nspinor**2,nbmax))
   ABI_ALLOCATE(fnltwg3,(Ep%npwepG0*nspinor**2,nbmax))
 endif

 if (niter>1) then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
   ABI_ALLOCATE(lnkp,(nfftot_gw*nspinor**2))
   ABI_ALLOCATE(kns,(Ep%npwepG0*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fdrhotwg,(Ep%npwepG0*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom,3))
   ABI_ALLOCATE(fkdwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(ff,(nlx*nlx,Cryst%natom,nlx*nlx,Cryst%natom))
   ABI_ALLOCATE(vzn,(Ep%npwepG0*nspinor**2,nlx*nlx,Cryst%natom))
 endif

 istwf_k  =  Wfs%istwfk(ik_ibz)
 npw_k    =  Wfs%npwarr(ik_ibz)
 gbound_k => Wfs%Kdata(ik_ibz)%gbound
 kg_k     => Wfs%Kdata(ik_ibz)%kg_k    
 igfft_k  => Wfs%Kdata(ik_ibz)%igfft0

 call calc_wfwfg(Wfs%MPI_enreg,Wfs%paral_kgb,tim_fourdp,tabr_k,itim_k,nfftot_gw,ngfft_gw, &
&                wfr2,wfr2,wfwfg(:))

 do ibvp = 1, nbmax
   call rho_tw_g(Wfs%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,wfr2,itim_k,tabr_k,ph_mkt, &
&                spinrot_k,dim_rtwg,rhotwg(:,ibvp),tim_fourdp,Wfs%MPI_enreg)

   if (ibvp>nbhomo(2)) then
     call calc_corr_chi0(ik_bz,Ep,Kmesh,Gsph_epsG0,Ltg_q,rhotwg(:,ibvp),spin_fact,qp_occ,qp_energy, &
&                        ibv,ibvp,ik_ibz,ikmq_ibz,is,chi0)
   endif
 enddo

 if (niter>0) then

   do i = 1, 3
     do ig = 1,npw_k
       igbz = grottbm1(ig,itim_k,isym_k)
       ph_Gt = EXP(-(0.,1.)*two_pi*DOT_PRODUCT(kg_k(:,ig),Cryst%tnons(:,isym_k)))
       if (itim_k==1) then
         gwfg(ig)=-kg_k(i,ig)*Wfs%Wave(ibv,ik_ibz,is)%ug(igbz)*ph_mkt*ph_Gt
       else
         gwfg(ig)=-kg_k(i,ig)*conjg(Wfs%Wave(ibv,ik_ibz,is)%ug(igbz)*ph_mkt)*ph_Gt
       endif
     enddo

     call fft_onewfn(Wfs%paral_kgb,istwf_k,nspinor,npw_k,Wfs%nfftot,Wfs%mgfft,Wfs%ngfft,&
&                      gwfg,dwfr(:,i),igfft_k,kg_k,gbound_k,tim_fourdp,Wfs%MPI_enreg)
   enddo

   do ibvp = 1, nbmax
     do i = 1, 3
       call drho_tw_g(Wfs%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,dwfr(:,i), &
&                     dim_rtwg,drhotwg(:,ibvp,i),tim_fourdp,Wfs%MPI_enreg)
     enddo
   enddo

   do i = 1, 3
     call calc_dwfwfg(Wfs%MPI_enreg,Wfs%paral_kgb,tim_fourdp,tabr_k,itim_k,nfftot_gw,ngfft_gw,ph_mkt, &
&                     wfr2,dwfr(:,i),dwfwfg(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(Wfs%MPI_enreg,Wfs%paral_kgb,tim_fourdp,itim_k,nfftot_gw,ngfft_gw, &
&                          dwfr(:,i),dwfr(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo

   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       call rho_tw_g(Wfs%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                    fnlkpr(:,ilm,iat),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,wfr2,itim_k,tabr_k,ph_mkt, &
&                    spinrot_k,dim_rtwg,fnltwg(:,ilm,iat),tim_fourdp,Wfs%MPI_enreg)
     enddo
   enddo

   call calc_wfwfg(Wfs%MPI_enreg,Wfs%paral_kgb,tim_fourdp,tabr_k,itim_k,nfftot_gw,ngfft_gw, &
&                  wfr2,mtwk(:,ibv),fnlwfg)

   do ibvp = 1, nbmax
     call rho_tw_g(Wfs%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                  wfr1(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,mtwk(:,ibv),itim_k,tabr_k,ph_mkt, &
&                  spinrot_k,dim_rtwg,fnltwg2(:,ibvp),tim_fourdp,Wfs%MPI_enreg)
     call rho_tw_g(Wfs%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                  mtwkp(:,ibvp),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,wfr2,itim_k,tabr_k,ph_mkt, &
&                  spinrot_k,dim_rtwg,fnltwg3(:,ibvp),tim_fourdp,Wfs%MPI_enreg)
  
   enddo

   if (niter>1) then

     do iat = 1, Cryst%natom
       do ilm = 1, Psps%mpsang*Psps%mpsang
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         do i = 1, 3
           call drho_tw_g(Wfs%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                         fnlkpr(:,ilm,iat),itim_kmq,tabr_kmq,ph_mkmqt, &
&                         dwfr(:,i),dim_rtwg,fdrhotwg(:,ilm,iat,i),tim_fourdp,Wfs%MPI_enreg)
         enddo
         call rho_tw_g(Wfs%paral_kgb,nspinor,Ep%npwepG0,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                       fnlkpr(:,ilm,iat),itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,mtwk(:,ibv),itim_k,tabr_k,ph_mkt, &
&                       spinrot_k,dim_rtwg,kns(:,ilm,iat),tim_fourdp,Wfs%MPI_enreg)
       enddo
     enddo

     do i = 1, 3
       call calc_dwfwfg(Wfs%MPI_enreg,Wfs%paral_kgb,tim_fourdp,tabr_k,itim_k,nfftot_gw,ngfft_gw,ph_mkt, &
&                       mtwk(:,ibv),dwfr(:,i),fkdwfg(:,i))
     enddo
     call calc_wfwfg(Wfs%MPI_enreg,Wfs%paral_kgb,tim_fourdp,tabr_k,itim_k,nfftot_gw,ngfft_gw, &
&                    mtwk(:,ibv),mtwk(:,ibv),lnkp)

   endif

   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr)
   ABI_ALLOCATE(cauxg,(Ep%npwe,nbmax))

   do ibvp = 1, nbmax
     do ig=1,Ep%npwe
       drhaux(:)=cmplx(real(drhotwg(ig,ibvp,:)),aimag(drhotwg(ig,ibvp,:)))
       cauxg(ig,ibvp)=vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
     enddo
   enddo
   do ig=1,Ep%npwe
     cauxg(ig,:)=cauxg(ig,:)-fnltwg2(ig,:)+fnltwg3(ig,:)
   enddo

 endif

 frhorho=czero_gw
 frhoj=czero_gw
 if (niter>1) fjj=czero_gw
 outofbox=0
 do ig=igstart,Ep%npwe
   do igp=igstart,Ep%npwe
     gmgp(:)=Gsph_epsG0%gvec(:,ig)-Gsph_epsG0%gvec(:,igp)
     ngfft(1)=ngfft_gw(1)
     ngfft(2)=ngfft_gw(2)
     ngfft(3)=ngfft_gw(3)
     if (ANY(gmgp(:)>ngfft(1:3)/2) .or. ANY(gmgp(:)<-(ngfft(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft(1))
     ig4y= modulo(gmgp(2),ngfft(2))
     ig4z= modulo(gmgp(3),ngfft(3))
     ig4= 1+ig4x+ig4y*ngfft(1)+ig4z*ngfft(1)*ngfft(2)
     ig5x= modulo(-gmgp(1),ngfft(1))
     ig5y= modulo(-gmgp(2),ngfft(2))
     ig5z= modulo(-gmgp(3),ngfft(3))
     ig5= 1+ig5x+ig5y*ngfft(1)+ig5z*ngfft(1)*ngfft(2)
     if (ig<=igp) then
       igaux=igp*(igp-1)/2+ig
       frhorho(igaux)=wfwfg(ig4)
     endif
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       frhoj(ig,igp)=frhoj(ig,igp) + vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")-fnlwfg(ig4)
       if (niter>1.and.ig<=igp) then
         igaux=igp*(igp-1)/2+ig
         paux2(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux2(i)=paux2(i) + vdotw(qplg(igp,:),drhaux,Cryst%gmet,"G")
         enddo
         drhaux(:)=paux2(:)-cmplx(real(fkdwfg(ig4,:)),aimag(fkdwfg(ig4,:)))
         fjj(igaux)=fjj(igaux) + vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
         drhaux(:)=cmplx(real(fkdwfg(ig5,:)),-aimag(fkdwfg(ig5,:)))
         fjj(igaux)=fjj(igaux) - vdotw(qplg(igp,:),drhaux,Cryst%gmet,"G")+lnkp(ig4)
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if

 do ibvp = 1, nbmax
   do ig=igstart,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       frhorho(igaux)=frhorho(igaux)-conjg(rhotwg(igp,ibvp))*rhotwg(ig,ibvp)
     enddo
   end do
 end do

 if (niter>0) then
   ABI_ALLOCATE(paux,(Ep%npwe,Ep%npwe))
   paux(:,:)=(0.0,0.0)
   minusone=(-1.,0.)
   do ibvp = 1, nbmax
     call XGERC(Ep%npwe,Ep%npwe,minusone,cauxg(:,ibvp),1,rhotwg(:,ibvp),1,frhoj,Ep%npwe)
   enddo
   do ig=igstart,Ep%npwe
     do igp=igstart,Ep%npwe
       if (ig>=igp) then
         do iat = 1, Cryst%natom
           do ilm = 1, nlx*nlx
             if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
             if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
             paux(ig,igp)=paux(ig,igp)+conjg(fnltwg(igp,ilm,iat))*fnltwg(ig,ilm,iat)
           enddo
         enddo
       endif
       if (ig<=igp) then
         igaux=igp*(igp-1)/2+ig
         faux=frhorho(igaux)
       else
         igaux=ig*(ig-1)/2+igp
         faux=conjg(frhorho(igaux))
       endif
       frhoj(ig,igp)=frhoj(ig,igp)+faux*kplqg(ig)
     end do !igp
   end do !ig
   do ig=igstart,Ep%npwe
     do igp=igstart,Ep%npwe
       if (ig>=igp) then
         frhoj(ig,igp)=frhoj(ig,igp)+paux(ig,igp)
       else
         frhoj(ig,igp)=frhoj(ig,igp)+conjg(paux(igp,ig))
       endif
     end do !igp
   end do !ig
   ABI_DEALLOCATE(paux)
 end if

 if (niter>1) then

  !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz
   call nullify_kb_potential(KBff_k_ibz  )

   call init_kb_potential(KBff_k_ibz ,Cryst,Psps,2,istwf_k,Ep%npwwfn,Kmesh%ibz(:,ik_ibz),kg_k)
   ABI_DEALLOCATE(KBff_k_ibz%fnld)

   ff(:,:,:,:)=czero_gw
   do iat = 1, Cryst%natom
     do ilm = 1, nlx*nlx
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       do iat2 = 1, Cryst%natom
         do ilm2 = 1, nlx*nlx
           if (ilm2>fnlmax(Cryst%typat(iat2))) CYCLE
           if (ilm2>=fnlloc(Cryst%typat(iat2),1).and.ilm2<=fnlloc(Cryst%typat(iat2),2)) CYCLE
           do ig = 1, Ep%npwwfn
             igbz = grottbm1(ig,itim_k,isym_k)
             if (itim_k==1) then
               ff(ilm,iat,ilm2,iat2)=ff(ilm,iat,ilm2,iat2)+ &
                conjg(KBff_k_ibz%fnl(igbz,ilm,iat))*KBff_k_ibz%fnl(igbz,ilm2,iat2)
             else
               ff(ilm,iat,ilm2,iat2)=ff(ilm,iat,ilm2,iat2)+ &
                KBff_k_ibz%fnl(igbz,ilm,iat)*conjg(KBff_k_ibz%fnl(igbz,ilm2,iat2))
             endif
           enddo
         enddo
       enddo
     enddo
   enddo

   call destroy_kb_potential(KBff_k_ibz)

   vzn(:,:,:)=czero_gw
   do ig=igstart,Ep%npwe
     do iat = 1, Cryst%natom
       do ilm = 1, nlx*nlx
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         vzn(ig,:,:) = vzn(ig,:,:) + half*ff(ilm,iat,:,:)*fnltwg(ig,ilm,iat)
       enddo
     enddo
     do iat = 1, Cryst%natom
       do ilm = 1, nlx*nlx
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         drhaux(:)=cmplx(real(fdrhotwg(ig,ilm,iat,:)),aimag(fdrhotwg(ig,ilm,iat,:)))
         vzn(ig,ilm,iat)=vzn(ig,ilm,iat)+vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")-kns(ig,ilm,iat)
       enddo
     enddo
   enddo
   do ibvp = 1, nbmax
     do ig=igstart,Ep%npwe
       do igp=ig,Ep%npwe
         igaux=igp*(igp-1)/2+ig
         fjj(igaux)=fjj(igaux)-conjg(cauxg(igp,ibvp))*cauxg(ig,ibvp)
       enddo
     enddo
   enddo
   do ig=igstart,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       faux=frhorho(igaux)
       fjj(igaux)=fjj(igaux)+kplqg(igp)*frhoj(ig,igp) + kplqg(ig)*conjg(frhoj(igp,ig))&
&                           -kplqg(igp)*faux*kplqg(ig)
       do iat = 1, Cryst%natom
         do ilm = 1, nlx*nlx
           if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
           if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
           fjj(igaux)=fjj(igaux)+conjg(fnltwg(igp,ilm,iat))*vzn(ig,ilm,iat)+ &
&                                      fnltwg(ig,ilm,iat)*conjg(vzn(igp,ilm,iat))
         enddo
       enddo
     end do !igp
   end do !ig

 end if

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg)
   ABI_DEALLOCATE(drhotwg)
   ABI_DEALLOCATE(dwfwfg)
   ABI_DEALLOCATE(fnltwg)
   ABI_DEALLOCATE(fnltwg2)
   ABI_DEALLOCATE(fnltwg3)
   ABI_DEALLOCATE(fnlwfg)
 endif
 if (niter>1) then
   ABI_DEALLOCATE(ddwfwfg)
   ABI_DEALLOCATE(lnkp)
   ABI_DEALLOCATE(kns)
   ABI_DEALLOCATE(fdrhotwg)
   ABI_DEALLOCATE(fkdwfg)
   ABI_DEALLOCATE(vzn)
   ABI_DEALLOCATE(ff)
 endif

end subroutine fft4eet_kb

!!***
 
!----------------------------------------------------------------------

!!****f* ABINIT/calc_eet_prep
!! NAME                  
!! calc_eet_prep
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine calc_eet_prep(Ep,Cryst,Wfs,Kmesh,Psps,is,nbhomo,nbmax,ik_ibz,ikmq_ibz,nspinor, &
&                        fnlloc,fnlmax,fnlkr,fnlkpr,mtwk,mtwkp)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_errors

 use m_oscillators,   only : rho_tw_g
 use m_geometry,      only : normv
 use m_vcoul,         only : vcoul_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_eet_prep'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfs_descriptor),intent(in) :: Wfs
 type(kb_potential) :: KBff_k_ibz,KBff_kmq_ibz

 integer,intent(in) :: is,nbhomo,nbmax
 integer,intent(in) :: nspinor
 integer,intent(in) :: ik_ibz,ikmq_ibz
 integer,intent(in) :: fnlloc(Cryst%ntypat,2)
 integer,intent(in) :: fnlmax(Cryst%ntypat)

 complex(gwpc),intent(in) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(in) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)

 complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbhomo)
 complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,nbmax)

 complex(gwpc),allocatable :: maux(:,:)

 integer :: ibv,ilm,iat,ig
 integer :: istwf_k,istwf_kmq,npw_k,npw_kmq
 integer,pointer :: kg_k(:,:),kg_kmq(:,:)

!************************************************************************

 ABI_UNUSED(Ep%npwe)

 !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz
 call nullify_kb_potential(KBff_k_ibz  )
 call nullify_kb_potential(KBff_kmq_ibz)

 npw_k   = Wfs%npwarr(ik_ibz)
 istwf_k = Wfs%istwfk(ik_ibz)
 kg_k    => Wfs%Kdata(ik_ibz)%kg_k

 npw_kmq   = Wfs%npwarr(ikmq_ibz)
 istwf_kmq = Wfs%istwfk(ikmq_ibz)
 kg_kmq    => Wfs%Kdata(ikmq_ibz)%kg_k

 ! MG: TODO Be careful here as the symmetry properties of fnl are not the same as the 
 !     ones of the wavefunctions. One should write a separate methods for the FFT of fnl.
 if (istwf_k>1.or.istwf_k>1) then
   MSG_ERROR("istwfk /= 1 not coded")
 end if

 call init_kb_potential(KBff_k_ibz  ,Cryst,Psps,2,istwf_k,  npw_k,  Kmesh%ibz(:,  ik_ibz),kg_k)
 call init_kb_potential(KBff_kmq_ibz,Cryst,Psps,2,istwf_kmq,npw_kmq,Kmesh%ibz(:,ikmq_ibz),kg_kmq)

 ABI_DEALLOCATE(KBff_k_ibz%fnld)
 ABI_DEALLOCATE(KBff_kmq_ibz%fnld)

 ABI_ALLOCATE(maux,(Psps%mpsang*Psps%mpsang,Cryst%natom))

 mtwk(:,:)=(0.0,0.0)
 mtwkp(:,:)=(0.0,0.0)
 do ibv = 1, nbhomo
   maux(:,:)=(0.0,0)
   do ig=1,npw_k
     maux(:,:) = maux(:,:) + Wfs%Wave(ibv,ik_ibz,is)%ug(ig) * KBff_k_ibz%fnl(ig,:,:)
   enddo
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       mtwk(:,ibv)=mtwk(:,ibv)+maux(ilm,iat)*fnlkr(:,ilm,iat)
     enddo
   enddo
 enddo

 do ibv = 1, nbmax
   maux(:,:)=(0.0,0)
   do ig=1,npw_kmq
     maux(:,:) = maux(:,:) + Wfs%Wave(ibv,ikmq_ibz,is)%ug(ig)* KBff_kmq_ibz%fnl(ig,:,:)
   enddo
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       mtwkp(:,ibv)=mtwkp(:,ibv)+maux(ilm,iat)*fnlkpr(:,ilm,iat)
     enddo
   enddo
 enddo

 ABI_DEALLOCATE(maux)

 call destroy_kb_potential(KBff_k_ibz)
 call destroy_kb_potential(KBff_kmq_ibz)

end subroutine calc_eet_prep

!!***
 
!----------------------------------------------------------------------

!!****f* ABINIT/fft4eet_sig
!! NAME                  
!! fft4eet_sig
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine fft4eet_sig(Sigp,Dtset,Cryst,Wfs,Gsph_c,Sr,nbhomo,nbmax,nomega,is,nfftot_gw,ngfft_gw, &
&                  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                  ik_ibz,ikmq_ibz,isym_kmq,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,grottbm1, &
&                  nspinor,tim_fourdp,MPI_enreg,wfr1,vc_sqrt_qbz,i_sz,kb,qplg,kplqg,niter, &
&                  ptwsq,ik_bz,ikmq_bz,npwc1,npwc2,botsq,otq,sigmac)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_errors

 use m_geometry,      only : vdotw
 use m_vcoul,         only : vcoul_t
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_sigma_results, only : sigma_results
 use m_blas,          only : xgerc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_sig'
 use interfaces_14_hidewrite
 use interfaces_70_gw, except_this_one => fft4eet_sig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(Sigma_parameters),intent(in) :: Sigp
 type(wfs_descriptor),intent(inout) :: Wfs
 type(Dataset_type),intent(in) :: Dtset
 type(Sigma_results),intent(in) :: Sr
 type(Gvectors_type),intent(in) :: Gsph_c
 type(MPI_type),intent(inout) :: MPI_enreg

 integer,intent(in) :: nbhomo,nbmax,is,kb,niter,npwc1,npwc2,nomega
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_bz,ikmq_bz,ik_ibz,ikmq_ibz
 integer,intent(in) :: isym_kmq,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Sigp%npwc)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: grottbm1(Sigp%npwvec,2,Cryst%nsym)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: qplg(Sigp%npwc,3),kplqg(Sigp%npwc)
 real(dp),intent(in) :: i_sz
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt

 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)

 complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
 complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)

 complex(gwpc),intent(out) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
 complex(dpc),intent(inout) :: sigmac(nomega)

 integer,pointer :: gbound_kmq(:,:),kg_kmq(:,:),igfft_kmq(:)
 complex(gwpc),allocatable :: wfr2(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: drhotwg(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)

 complex(gwpc),allocatable :: dwfr(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg(:,:)

 integer :: istwf_kmq,npw_kmq
 integer :: ibv,ig,igp,igbz
 integer :: i,j
 integer :: ig4,ig4x,ig4y,ig4z
 integer :: gmgp(3)
 integer :: outofbox
 integer,save :: enough=0
 character(len=500) :: msg

 complex(dpc) :: ph_Gt

 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux(3)

 complex(gwpc) :: minusone

!************************************************************************

 gbound_kmq => Wfs%Kdata(ikmq_ibz)%gbound
 igfft_kmq  => Wfs%Kdata(ikmq_ibz)%igfft0
 kg_kmq     => Wfs%Kdata(ikmq_ibz)%kg_k
 istwf_kmq  = Wfs%istwfk(ikmq_ibz)
 npw_kmq    = Wfs%npwarr(ikmq_ibz)

 ABI_ALLOCATE(wfr2,(Wfs%nfftot*nspinor))

 ABI_ALLOCATE(rhotwg,(Sigp%npwc*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 if (niter>0) then
   ABI_ALLOCATE(drhotwg,(Sigp%npwc*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(Sigp%npwwfn))
   ABI_ALLOCATE(dwfr,(Wfs%nfftot*nspinor,3))
 endif

 if (niter>1) then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
 endif

 call wfd_get_ur(Wfs,kb,ikmq_ibz,is,wfr2)
 call calc_wfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw, &
&                wfr2,wfr2,wfwfg(:))

 do ibv = 1, nbmax
   call rho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                spinrot_kmq,nspinor,rhotwg(:,ibv),tim_fourdp,MPI_enreg)

   if (ibv>nbhomo) then
     call calc_corr_sig(Sigp,Sr,nomega,nspinor,npwc1,npwc2,botsq,otq,rhotwg(:,ibv), &
&                       is,ibv,kb,ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,i_sz,vc_sqrt_qbz,sigmac)
   endif
 enddo

 if (niter>0) then

   do i = 1, 3
     do ig = 1, Sigp%npwwfn
       igbz = grottbm1(ig,itim_kmq,isym_kmq)
       ph_Gt = EXP(-(0.,1.)*two_pi*DOT_PRODUCT(kg_kmq(:,ig),Cryst%tnons(:,isym_kmq)))
       if (itim_kmq==1) then
         gwfg(ig)=-kg_kmq(i,ig)*Wfs%Wave(kb,ikmq_ibz,is)%ug(igbz)*ph_mkmqt*ph_Gt
       else
         gwfg(ig)=-kg_kmq(i,ig)*conjg(Wfs%Wave(kb,ikmq_ibz,is)%ug(igbz)*ph_mkmqt)*ph_Gt
       endif
     enddo

     call fft_onewfn(Dtset%paral_kgb,istwf_kmq,nspinor,npw_kmq,Wfs%nfftot,Wfs%mgfft,Wfs%ngfft,&
&                    gwfg,dwfr(:,i),igfft_kmq,kg_kmq,gbound_kmq,tim_fourdp,MPI_enreg)
   enddo

   do ibv = 1, nbmax
     do i = 1, 3
       call drho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr1(:,ibv),itim_k,tabr_k,ph_mkt,dwfr(:,i), &
&                     nspinor,drhotwg(:,ibv,i),tim_fourdp,MPI_enreg)
     enddo
   enddo

   do i = 1, 3
     call calc_dwfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt, &
&                     wfr2,dwfr(:,i),dwfwfg(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,itim_kmq,nfftot_gw,ngfft_gw, &
&                          dwfr(:,i),dwfr(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo

   ABI_DEALLOCATE(wfr2)
   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr)
   ABI_ALLOCATE(cauxg,(Sigp%npwc,nbmax))

   cauxg(:,:)=(0.0,0.0)
   do ibv = 1, nbmax
     do ig = 1, Sigp%npwc
       drhaux(:)=cmplx(real(drhotwg(ig,ibv,:)),aimag(drhotwg(ig,ibv,:)))
       cauxg(ig,ibv)=vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
     enddo
   enddo

 endif

 ptwsq(:,:,:)=(0.0,0.0)
 outofbox=0
 do igp=1,Sigp%npwc
   do ig=1,Sigp%npwc
     gmgp(:)=Gsph_c%gvec(:,igp)-Gsph_c%gvec(:,ig)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)
     if (igp>=ig) then
       ptwsq(ig,igp,1)=wfwfg(ig4)
     endif
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2) + vdotw(qplg(igp,:),drhaux,Cryst%gmet,"G")
       if (niter>1.and.igp>=ig) then
         paux(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux(i)=paux(i) + vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
         enddo
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3) + vdotw(qplg(igp,:),paux,Cryst%gmet,"G")
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if

 do ibv = 1, nbmax
   do igp=1,Sigp%npwc
     do ig=1,igp
       ptwsq(ig,igp,1)=ptwsq(ig,igp,1)-conjg(rhotwg(ig,ibv))*rhotwg(igp,ibv)
     enddo
   end do
 end do

 do igp = 1, Sigp%npwc
   do ig = igp+1, Sigp%npwc
     ptwsq(ig,igp,1)=conjg(ptwsq(igp,ig,1))
   enddo
 enddo

 if (niter>0) then
   minusone=(-1.,0.)
   do ibv = 1, nbmax
     call XGERC(Sigp%npwc,Sigp%npwc,minusone,conjg(rhotwg(:,ibv)),1,conjg(cauxg(:,ibv)),1,ptwsq(:,:,2),Sigp%npwc)
   enddo
   do igp=1,Sigp%npwc
     do ig=1,Sigp%npwc
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+ptwsq(ig,igp,1)*kplqg(igp)
     end do !igp
   end do !ig
 endif

 if (niter>1) then

   do ibv = 1, nbmax
     do igp=1,Sigp%npwc
       do ig=1,igp
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3)-conjg(cauxg(ig,ibv))*cauxg(igp,ibv)
       enddo
     enddo
   enddo
   do igp=1,Sigp%npwc
     do ig=1,igp
       ptwsq(ig,igp,3)=ptwsq(ig,igp,3)+kplqg(igp)*ptwsq(ig,igp,2)+kplqg(ig)*conjg(ptwsq(igp,ig,2))&
&                                     -kplqg(ig)*ptwsq(ig,igp,1)*kplqg(igp)
     end do !igp
   end do !ig

   do igp = 1, Sigp%npwc
     do ig = igp+1, Sigp%npwc
       ptwsq(ig,igp,3)=conjg(ptwsq(igp,ig,3))
     enddo
   enddo

 endif

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg)
   ABI_DEALLOCATE(drhotwg)
   ABI_DEALLOCATE(dwfwfg)
 endif
 if (niter>1) then
   ABI_DEALLOCATE(ddwfwfg)
 endif

end subroutine fft4eet_sig

!!***
!----------------------------------------------------------------------

!!****f* ABINIT/fft4eet_sig_sc
!! NAME                  
!! fft4eet_sig_sc
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine fft4eet_sig_sc(Sigp,Dtset,Cryst,Wfs,Gsph_c,Sr,nbhomo,nbmax,nomega,is,nfftot_gw,ngfft_gw, &
&                  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                  ik_ibz,ikmq_ibz,isym_kmq,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,grottbm1, &
&                  nspinor,tim_fourdp,MPI_enreg,wfr,vc_sqrt_qbz,i_sz,jb,kb,qplg,kplqg,niter, &
&                  ptwsq,ik_bz,ikmq_bz,npwc1,npwc2,botsq,otq,sigmac)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_errors

 use m_geometry,      only : vdotw
 use m_vcoul,         only : vcoul_t
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_sigma_results, only : sigma_results

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_sig_sc'
 use interfaces_14_hidewrite
 use interfaces_70_gw, except_this_one => fft4eet_sig_sc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(Sigma_parameters),intent(in) :: Sigp
 type(wfs_descriptor),intent(inout) :: Wfs
 type(Dataset_type),intent(in) :: Dtset
 type(Sigma_results),intent(in) :: Sr
 type(Gvectors_type),intent(in) :: Gsph_c
 type(MPI_type),intent(inout) :: MPI_enreg

 integer,intent(in) :: nbhomo,nbmax,is,jb,kb,niter,npwc1,npwc2,nomega
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_bz,ikmq_bz,ik_ibz,ikmq_ibz
                                                       
 integer,intent(in) :: isym_kmq,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Sigp%npwc)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: grottbm1(Sigp%npwvec,2,Cryst%nsym)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: qplg(Sigp%npwc,3),kplqg(Sigp%npwc)
 real(dp),intent(in) :: i_sz
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt

 complex(gwpc),intent(in) :: wfr(Wfs%nfftot*nspinor,nbmax)

 complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
 complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)

 complex(gwpc),intent(out) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
 complex(dpc),intent(inout) :: sigmac(nomega)

 integer,pointer :: gbound_kmq(:,:),kg_kmq(:,:),igfft_kmq(:)
 complex(gwpc),allocatable :: wfr_kb(:)
 complex(gwpc),allocatable :: wfr_jb(:)
 complex(gwpc),allocatable :: rhotwg_kb(:,:)
 complex(gwpc),allocatable :: rhotwg_jb(:,:)
 complex(gwpc),allocatable :: drhotwg_kb(:,:,:)
 complex(gwpc),allocatable :: drhotwg_jb(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: dwfwfg2(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)

 complex(gwpc),allocatable :: dwfr_kb(:,:)
 complex(gwpc),allocatable :: dwfr_jb(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg_kb(:,:)
 complex(gwpc),allocatable :: cauxg_jb(:,:)
 complex(gwpc),allocatable :: taux(:,:)
 complex(gwpc),allocatable :: taux2(:,:)

 integer :: istwf_kmq,npw_kmq
 integer :: ibv,ig,igp,igbz
 integer :: i,j
 integer :: ig4,ig4x,ig4y,ig4z
 integer :: ig5,ig5x,ig5y,ig5z
 integer :: gmgp(3)
 integer :: outofbox
 integer,save :: enough=0
 character(len=500) :: msg

 complex(dpc) :: ph_Gt

 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux(3)

!************************************************************************

 gbound_kmq => Wfs%Kdata(ikmq_ibz)%gbound
 igfft_kmq  => Wfs%Kdata(ikmq_ibz)%igfft0
 kg_kmq     => Wfs%Kdata(ikmq_ibz)%kg_k
 istwf_kmq  = Wfs%istwfk(ikmq_ibz)
 npw_kmq    = Wfs%npwarr(ikmq_ibz)

 ABI_ALLOCATE(wfr_jb,(Wfs%nfftot*nspinor))
 ABI_ALLOCATE(wfr_kb,(Wfs%nfftot*nspinor))

 ABI_ALLOCATE(rhotwg_kb,(Sigp%npwc*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))
 ABI_ALLOCATE(rhotwg_jb,(Sigp%npwc*nspinor**2,nbmax))

 if (niter>0) then
   ABI_ALLOCATE(drhotwg_kb,(Sigp%npwc*nspinor**2,nbmax,3))
   ABI_ALLOCATE(drhotwg_jb,(Sigp%npwc*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(dwfwfg2,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(Sigp%npwwfn))
   ABI_ALLOCATE(dwfr_kb,(Wfs%nfftot*nspinor,3))
   ABI_ALLOCATE(dwfr_jb,(Wfs%nfftot*nspinor,3))
   ABI_ALLOCATE(taux,(Sigp%npwc*nspinor**2,Sigp%npwc*nspinor**2))
   ABI_ALLOCATE(taux2,(Sigp%npwc*nspinor**2,Sigp%npwc*nspinor**2))
 endif

 if (niter>1) then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
 endif

 call wfd_get_ur(Wfs,kb,ikmq_ibz,is,wfr_kb)
 call wfd_get_ur(Wfs,jb,ikmq_ibz,is,wfr_jb)

 call calc_wfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw, &
&                wfr_jb,wfr_kb,wfwfg(:))

 do ibv = 1, nbmax
   call rho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr_kb,itim_kmq,tabr_kmq,ph_mkmqt, &
&                spinrot_kmq,nspinor,rhotwg_kb(:,ibv),tim_fourdp,MPI_enreg)

   call rho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr_jb,itim_kmq,tabr_kmq,ph_mkmqt, &
&                spinrot_kmq,nspinor,rhotwg_jb(:,ibv),tim_fourdp,MPI_enreg)
   if (ibv>nbhomo.and..false.) then
     call calc_corr_sig(Sigp,Sr,nomega,nspinor,npwc1,npwc2,botsq,otq,rhotwg_kb(:,ibv), &
&                       is,ibv,kb,ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,i_sz,vc_sqrt_qbz,sigmac)
   endif
 enddo

 if (niter>0) then

   do i = 1, 3
     do ig = 1, Sigp%npwwfn
       igbz = grottbm1(ig,itim_kmq,isym_kmq)
       ph_Gt = EXP(-(0.,1.)*two_pi*DOT_PRODUCT(kg_kmq(:,ig),Cryst%tnons(:,isym_kmq)))
       if (itim_kmq==1) then
         gwfg(ig)=-kg_kmq(i,ig)*Wfs%Wave(kb,ikmq_ibz,is)%ug(igbz)*ph_mkmqt*ph_Gt
       else
         gwfg(ig)=-kg_kmq(i,ig)*conjg(Wfs%Wave(kb,ikmq_ibz,is)%ug(igbz)*ph_mkmqt)*ph_Gt
       endif
     enddo

     call fft_onewfn(Dtset%paral_kgb,istwf_kmq,nspinor,npw_kmq,Wfs%nfftot,Wfs%mgfft,Wfs%ngfft,&
&                    gwfg,dwfr_kb(:,i),igfft_kmq,kg_kmq,gbound_kmq,tim_fourdp,MPI_enreg)

     do ig = 1, Sigp%npwwfn
       igbz = grottbm1(ig,itim_kmq,isym_kmq)
       ph_Gt = EXP(-(0.,1.)*two_pi*DOT_PRODUCT(kg_kmq(:,ig),Cryst%tnons(:,isym_kmq)))
       if (itim_kmq==1) then
         gwfg(ig)=-kg_kmq(i,ig)*Wfs%Wave(jb,ikmq_ibz,is)%ug(igbz)*ph_mkmqt*ph_Gt
       else
         gwfg(ig)=-kg_kmq(i,ig)*conjg(Wfs%Wave(jb,ikmq_ibz,is)%ug(igbz)*ph_mkmqt)*ph_Gt
       endif
     enddo

     call fft_onewfn(Dtset%paral_kgb,istwf_kmq,nspinor,npw_kmq,Wfs%nfftot,Wfs%mgfft,Wfs%ngfft,&
&                    gwfg,dwfr_jb(:,i),igfft_kmq,kg_kmq,gbound_kmq,tim_fourdp,MPI_enreg)

   enddo

   do ibv = 1, nbmax
     do i = 1, 3
       call drho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr(:,ibv),itim_k,tabr_k,ph_mkt,dwfr_kb(:,i), &
&                     nspinor,drhotwg_kb(:,ibv,i),tim_fourdp,MPI_enreg)
       call drho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr(:,ibv),itim_k,tabr_k,ph_mkt,dwfr_jb(:,i), &
&                     nspinor,drhotwg_jb(:,ibv,i),tim_fourdp,MPI_enreg)
     enddo

   enddo

   do i = 1, 3
     call calc_dwfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt, &
&                     wfr_jb,dwfr_kb(:,i),dwfwfg(:,i))
     call calc_dwfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt, &
&                     wfr_kb,dwfr_jb(:,i),dwfwfg2(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,itim_kmq,nfftot_gw,ngfft_gw, &
&                          dwfr_jb(:,i),dwfr_kb(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo

   ABI_DEALLOCATE(wfr_kb)
   ABI_DEALLOCATE(wfr_jb)
   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr_kb)
   ABI_DEALLOCATE(dwfr_jb)
   ABI_ALLOCATE(cauxg_kb,(Sigp%npwc,nbmax))
   ABI_ALLOCATE(cauxg_jb,(Sigp%npwc,nbmax))

   cauxg_kb(:,:)=(0.0,0.0)
   do ibv = 1, nbmax
     do ig=1,Sigp%npwc
       drhaux(:)=cmplx(real(drhotwg_kb(ig,ibv,:)),aimag(drhotwg_kb(ig,ibv,:)))
       cauxg_kb(ig,ibv)=vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
     enddo
   enddo

   cauxg_jb(:,:)=(0.0,0.0)
   do ibv = 1, nbmax
     do ig=1,Sigp%npwc
       drhaux(:)=cmplx(real(drhotwg_jb(ig,ibv,:)),aimag(drhotwg_jb(ig,ibv,:)))
       cauxg_jb(ig,ibv)=vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
     enddo
   enddo

 endif

 ptwsq(:,:,:)=(0.0,0.0)
 if (niter>0) then
   taux=czero_gw
   taux2=czero_gw
 endif
 outofbox=0
 do igp=1,Sigp%npwc
   do ig=1,Sigp%npwc
     gmgp(:)=Gsph_c%gvec(:,igp)-Gsph_c%gvec(:,ig)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)
     ig5x= modulo(-gmgp(1),ngfft_gw(1))
     ig5y= modulo(-gmgp(2),ngfft_gw(2))
     ig5z= modulo(-gmgp(3),ngfft_gw(3))
     ig5= 1+ig5x+ig5y*ngfft_gw(1)+ig5z*ngfft_gw(1)*ngfft_gw(2)
     ptwsq(ig,igp,1)=wfwfg(ig4)
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       taux(ig,igp)=taux(ig,igp) + vdotw(qplg(igp,:),drhaux,Cryst%gmet,"G")
       drhaux(:)=cmplx(real(dwfwfg2(ig5,:)),-aimag(dwfwfg2(ig5,:)))
       taux2(ig,igp)=taux2(ig,igp) + vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
       if (niter>1) then
         paux(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux(i)=paux(i) + vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
         enddo
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3) + vdotw(qplg(igp,:),paux,Cryst%gmet,"G")
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if

 do ibv = 1, nbmax
   do igp=1,Sigp%npwc
     do ig=1,Sigp%npwc
       ptwsq(ig,igp,1)=ptwsq(ig,igp,1)-conjg(rhotwg_jb(ig,ibv))*rhotwg_kb(igp,ibv)
     enddo
   end do
 end do

 if (niter>0) then
   do igp=1,Sigp%npwc
     do ig=1,Sigp%npwc
       do ibv = 1, nbmax
         taux(ig,igp)=taux(ig,igp)-conjg(rhotwg_jb(ig,ibv))*cauxg_kb(igp,ibv)
         taux2(ig,igp)=taux2(ig,igp)-rhotwg_kb(igp,ibv)*conjg(cauxg_jb(ig,ibv))
       enddo
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+half*(kplqg(ig)+kplqg(igp))*ptwsq(ig,igp,1)
     end do !igp
   end do !ig
   ptwsq(:,:,2)=ptwsq(:,:,2)+half*(taux(:,:)+taux2(:,:))
 endif

 if (niter>1) then

   do igp=1,Sigp%npwc
     do ig=1,Sigp%npwc
       do ibv = 1, nbmax
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3)-conjg(cauxg_jb(ig,ibv))*cauxg_kb(igp,ibv)
       enddo
       ptwsq(ig,igp,3)=ptwsq(ig,igp,3)+kplqg(ig)*taux(ig,igp)+kplqg(igp)*taux2(ig,igp)&
&                                     +kplqg(ig)*ptwsq(ig,igp,1)*kplqg(igp)
     end do !igp
   end do !ig

 endif

 ABI_DEALLOCATE(rhotwg_kb)
 ABI_DEALLOCATE(rhotwg_jb)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg_kb)
   ABI_DEALLOCATE(cauxg_jb)
   ABI_DEALLOCATE(drhotwg_kb)
   ABI_DEALLOCATE(drhotwg_jb)
   ABI_DEALLOCATE(dwfwfg)
   ABI_DEALLOCATE(dwfwfg2)
   ABI_DEALLOCATE(taux)
   ABI_DEALLOCATE(taux2)
 endif
 if (niter>1) then
   ABI_DEALLOCATE(ddwfwfg)
 endif

end subroutine fft4eet_sig_sc
!!***
!----------------------------------------------------------------------

!!****f* ABINIT/fft4eet_sig_kb
!! NAME                  
!! fft4eet_sig_kb
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine fft4eet_sig_kb(Sigp,Dtset,Cryst,Wfs,Kmesh,Gsph_c,Psps,Sr,nbhomo,nbmax,nomega,is,nfftot_gw,ngfft_gw, &
&                  use_padfft,igfftepsG0,gw_gbound,gw_mgfft,itim_k,tabr_k,ph_mkt,spinrot_k, &
&                  ik_ibz,ikmq_ibz,isym_kmq,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,grottbm1, &
&                  nspinor,tim_fourdp,MPI_enreg,fnlloc,fnlmax,fnlkr,mtwk,mtwkp,wfr1, &
&                  vc_sqrt_qbz,i_sz,kb,qplg,kplqg,niter,ptwsq,ik_bz,ikmq_bz, &
&                  npwc1,npwc2,botsq,otq,sigmac)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_errors

 use m_geometry,      only : vdotw
 use m_vcoul,         only : vcoul_t
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_sigma_results, only : sigma_results
 use m_blas,          only : xgerc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft4eet_sig_kb'
 use interfaces_14_hidewrite
 use interfaces_70_gw, except_this_one => fft4eet_sig_kb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Sigma_parameters),intent(in) :: Sigp
 type(wfs_descriptor),intent(inout) :: Wfs
 type(kb_potential) :: KBff_kmq_ibz
 type(Dataset_type),intent(in) :: Dtset
 type(Sigma_results),intent(in) :: Sr
 type(Gvectors_type),intent(in) :: Gsph_c
 type(MPI_type),intent(inout) :: MPI_enreg

 integer,intent(in) :: nbhomo,nbmax,is,kb,niter,npwc1,npwc2,nomega
 integer,intent(in) :: nfftot_gw,ngfft_gw(18),nspinor,tim_fourdp,use_padfft,gw_mgfft
 integer,intent(in) :: ik_bz,ikmq_bz,ik_ibz,ikmq_ibz
 integer,intent(in) :: isym_kmq,itim_k,itim_kmq
 integer,intent(in) :: tabr_k(nfftot_gw),tabr_kmq(nfftot_gw)
 integer,intent(in) :: igfftepsG0(Sigp%npwc)
 integer,intent(in) :: gw_gbound(2*gw_mgfft+8,2*use_padfft)
 integer,intent(in) :: grottbm1(Sigp%npwvec,2,Cryst%nsym)
 integer,intent(in) :: fnlloc(Cryst%ntypat,2)
 integer,intent(in) :: fnlmax(Cryst%ntypat)
 real(dp),intent(in) :: spinrot_k(4),spinrot_kmq(4)
 real(dp),intent(in) :: qplg(Sigp%npwc,3),kplqg(Sigp%npwc)
 real(dp),intent(in) :: i_sz
 complex(dpc),intent(in) :: ph_mkmqt,ph_mkt

 complex(gwpc),intent(in) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(in) :: mtwk(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(in) :: mtwkp(Wfs%nfftot*nspinor)

 complex(gwpc),intent(in) :: wfr1(Wfs%nfftot*nspinor,nbmax)

 complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
 complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)

 complex(gwpc),intent(out) :: ptwsq(Sigp%npwc,Sigp%npwc,niter+1)
 complex(dpc),intent(inout) :: sigmac(nomega)

 integer,pointer :: gbound_kmq(:,:),kg_kmq(:,:),igfft_kmq(:)
 complex(gwpc),allocatable :: wfr2(:)
 complex(gwpc),allocatable :: rhotwg(:,:)
 complex(gwpc),allocatable :: drhotwg(:,:,:)
 complex(gwpc),allocatable :: fnltwg(:,:,:)
 complex(gwpc),allocatable :: fnltwg2(:,:)
 complex(gwpc),allocatable :: fnltwg3(:,:)
 complex(gwpc),allocatable :: kns(:,:,:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: dwfwfg(:,:)
 complex(gwpc),allocatable :: ddwfwfg(:,:,:)
 complex(gwpc),allocatable :: fnlwfg(:)
 complex(gwpc),allocatable :: fkdwfg(:,:)
 complex(gwpc),allocatable :: fdrhotwg(:,:,:,:)
 complex(gwpc),allocatable :: lnkp(:)

 complex(gwpc),allocatable :: dwfr(:,:)
 complex(gwpc),allocatable :: gwfg(:)
 complex(gwpc),allocatable :: cauxg(:,:)
 complex(gwpc),allocatable :: ff(:,:,:,:)
 complex(gwpc),allocatable :: vzn(:,:,:)

 complex(gwpc),allocatable :: paux(:,:)

 integer :: istwf_kmq,npw_kmq
 integer :: ibv,ilm,iat,ilm2,iat2,ig,igp,igbz
 integer :: i,j
 integer :: ig4,ig4x,ig4y,ig4z
 integer :: ig5,ig5x,ig5y,ig5z
 integer :: nlx
 integer :: gmgp(3)
 integer :: outofbox
 integer,save :: enough=0
 character(len=500) :: msg

 complex(dpc) :: ph_Gt

 complex(dpc) :: drhaux(3)
 complex(dpc) :: paux2(3)

 complex(gwpc) :: minusone

!************************************************************************

 gbound_kmq => Wfs%Kdata(ikmq_ibz)%gbound
 igfft_kmq  => Wfs%Kdata(ikmq_ibz)%igfft0
 kg_kmq     => Wfs%Kdata(ikmq_ibz)%kg_k
 istwf_kmq  = Wfs%istwfk(ikmq_ibz)
 npw_kmq    = Wfs%npwarr(ikmq_ibz)

 nlx = min(Psps%mpsang,4)

 ABI_ALLOCATE(wfr2,(Wfs%nfftot*nspinor))
 ABI_ALLOCATE(rhotwg,(Sigp%npwc*nspinor**2,nbmax))
 ABI_ALLOCATE(wfwfg,(nfftot_gw*nspinor**2))

 if (niter>0) then
   ABI_ALLOCATE(drhotwg,(Sigp%npwc*nspinor**2,nbmax,3))
   ABI_ALLOCATE(dwfwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(gwfg,(Sigp%npwwfn))
   ABI_ALLOCATE(dwfr,(Wfs%nfftot*nspinor,3))
   ABI_ALLOCATE(fnltwg,(Sigp%npwc*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fnlwfg,(nfftot_gw*nspinor**2))
   ABI_ALLOCATE(fnltwg2,(Sigp%npwc*nspinor**2,nbmax))
   ABI_ALLOCATE(fnltwg3,(Sigp%npwc*nspinor**2,nbmax))
 endif

 if (niter>1) then
   ABI_ALLOCATE(ddwfwfg,(nfftot_gw*nspinor**2,3,3))
   ABI_ALLOCATE(lnkp,(nfftot_gw*nspinor**2))
   ABI_ALLOCATE(kns,(Sigp%npwc*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom))
   ABI_ALLOCATE(fdrhotwg,(Sigp%npwc*nspinor**2,Psps%mpsang*Psps%mpsang,Cryst%natom,3))
   ABI_ALLOCATE(fkdwfg,(nfftot_gw*nspinor**2,3))
   ABI_ALLOCATE(ff,(nlx*nlx,Cryst%natom,nlx*nlx,Cryst%natom))
   ABI_ALLOCATE(vzn,(Sigp%npwc*nspinor**2,nlx*nlx,Cryst%natom))
 endif

 call wfd_get_ur(Wfs,kb,ikmq_ibz,is,wfr2)
 call calc_wfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw, &
&                wfr2,wfr2,wfwfg(:))

 do ibv = 1, nbmax
   call rho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                wfr1(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                spinrot_kmq,nspinor,rhotwg(:,ibv),tim_fourdp,MPI_enreg)

   if (ibv>nbhomo) then
     call calc_corr_sig(Sigp,Sr,nomega,nspinor,npwc1,npwc2,botsq,otq,rhotwg(:,ibv), &
&                       is,ibv,kb,ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,i_sz,vc_sqrt_qbz,sigmac)
   endif
 enddo

 if (niter>0) then

   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       call rho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                    fnlkr(:,ilm,iat),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                    spinrot_kmq,nspinor,fnltwg(:,ilm,iat),tim_fourdp,MPI_enreg)
     enddo
   enddo

   call calc_wfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw, &
&                  wfr2,mtwkp,fnlwfg)

   do i = 1, 3
     do ig = 1, Sigp%npwwfn
       igbz = grottbm1(ig,itim_kmq,isym_kmq)
       ph_Gt = EXP(-(0.,1.)*two_pi*DOT_PRODUCT(kg_kmq(:,ig),Cryst%tnons(:,isym_kmq)))
       if (itim_kmq==1) then
         gwfg(ig)=-kg_kmq(i,ig)*Wfs%Wave(kb,ikmq_ibz,is)%ug(igbz)*ph_mkmqt*ph_Gt
       else
         gwfg(ig)=-kg_kmq(i,ig)*conjg(Wfs%Wave(kb,ikmq_ibz,is)%ug(igbz)*ph_mkmqt)*ph_Gt
       endif
     enddo

     call fft_onewfn(Dtset%paral_kgb,istwf_kmq,nspinor,npw_kmq,Wfs%nfftot,Wfs%mgfft,Wfs%ngfft,&
&                    gwfg,dwfr(:,i),igfft_kmq,kg_kmq,gbound_kmq,tim_fourdp,MPI_enreg)
   enddo

   do ibv = 1, nbmax
     do i = 1, 3
       call drho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                     wfr1(:,ibv),itim_k,tabr_k,ph_mkt,dwfr(:,i), &
&                     nspinor,drhotwg(:,ibv,i),tim_fourdp,MPI_enreg)
     enddo
     call rho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                  wfr1(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,mtwkp,itim_kmq,tabr_kmq,ph_mkmqt, &
&                  spinrot_kmq,nspinor,fnltwg2(:,ibv),tim_fourdp,MPI_enreg)
     call rho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound, &
&                  mtwk(:,ibv),itim_k,tabr_k,ph_mkt,spinrot_k,wfr2,itim_kmq,tabr_kmq,ph_mkmqt, &
&                  spinrot_kmq,nspinor,fnltwg3(:,ibv),tim_fourdp,MPI_enreg)
   enddo

   do i = 1, 3
     call calc_dwfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt, &
&                     wfr2,dwfr(:,i),dwfwfg(:,i))
     if (niter>1) then
       do j = 1, 3
         call calc_ddwfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,itim_kmq,nfftot_gw,ngfft_gw, &
&                          dwfr(:,i),dwfr(:,j),ddwfwfg(:,i,j))
       enddo
     endif
   enddo

   if (niter>1) then

     do iat = 1, Cryst%natom
       do ilm = 1, Psps%mpsang*Psps%mpsang
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         do i = 1, 3
           call drho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                         fnlkr(:,ilm,iat),itim_k,tabr_k,ph_mkt,dwfr(:,i), &
&                         nspinor,fdrhotwg(:,ilm,iat,i),tim_fourdp,MPI_enreg)
         enddo
         call rho_tw_g(Dtset%paral_kgb,nspinor,Sigp%npwc,nfftot_gw,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&                      fnlkr(:,ilm,iat),itim_k,tabr_k,ph_mkt,spinrot_k,mtwkp,itim_kmq,tabr_kmq,ph_mkmqt, &
&                      spinrot_kmq,nspinor,kns(:,ilm,iat),tim_fourdp,MPI_enreg)
       enddo
     enddo

     do i = 1, 3
       call calc_dwfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw,ph_mkmqt, &
&                       mtwkp,dwfr(:,i),fkdwfg(:,i))
     enddo
     call calc_wfwfg(MPI_enreg,Dtset%paral_kgb,tim_fourdp,tabr_kmq,itim_kmq,nfftot_gw,ngfft_gw, &
&                    mtwkp,mtwkp,lnkp)

   endif

   ABI_DEALLOCATE(wfr2)
   ABI_DEALLOCATE(gwfg)
   ABI_DEALLOCATE(dwfr)
   ABI_ALLOCATE(cauxg,(Sigp%npwc,nbmax))

   cauxg(:,:)=(0.0,0.0)
   do ibv = 1, nbmax
     do ig=1,Sigp%npwc
       drhaux(:)=cmplx(real(drhotwg(ig,ibv,:)),aimag(drhotwg(ig,ibv,:)))
       cauxg(ig,ibv)=vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
     enddo
   enddo
   do ig=1,Sigp%npwc
     cauxg(ig,:)=cauxg(ig,:)-fnltwg2(ig,:)+fnltwg3(ig,:)
   enddo

 endif

 ptwsq(:,:,:)=(0.0,0.0)
 outofbox=0
 do igp=1,Sigp%npwc
   do ig=1,Sigp%npwc
     gmgp(:)=Gsph_c%gvec(:,igp)-Gsph_c%gvec(:,ig)
     if (ANY(gmgp(:)>ngfft_gw(1:3)/2) .or. ANY(gmgp(:)<-(ngfft_gw(1:3)-1)/2)) then
       outofbox = outofbox+1; CYCLE
     end if
     ig4x= modulo(gmgp(1),ngfft_gw(1))
     ig4y= modulo(gmgp(2),ngfft_gw(2))
     ig4z= modulo(gmgp(3),ngfft_gw(3))
     ig4= 1+ig4x+ig4y*ngfft_gw(1)+ig4z*ngfft_gw(1)*ngfft_gw(2)

     ig5x= modulo(-gmgp(1),ngfft_gw(1))
     ig5y= modulo(-gmgp(2),ngfft_gw(2))
     ig5z= modulo(-gmgp(3),ngfft_gw(3))
     ig5= 1+ig5x+ig5y*ngfft_gw(1)+ig5z*ngfft_gw(1)*ngfft_gw(2)

     if (igp>=ig) then
       ptwsq(ig,igp,1)=wfwfg(ig4)
     endif
     if (niter>0) then
       drhaux(:)=cmplx(real(dwfwfg(ig4,:)),aimag(dwfwfg(ig4,:)))
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2) + vdotw(qplg(igp,:),drhaux,Cryst%gmet,"G")-fnlwfg(ig4)
       if (niter>1.and.igp>=ig) then
         paux2(:)=czero
         do i = 1, 3
           drhaux(:)=cmplx(real(ddwfwfg(ig4,:,i)),aimag(ddwfwfg(ig4,:,i)))
           paux2(i)=paux2(i) + vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")
         enddo
         drhaux(:)=paux2(:)-cmplx(real(fkdwfg(ig4,:)),aimag(fkdwfg(ig4,:)))
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3) + vdotw(qplg(igp,:),drhaux,Cryst%gmet,"G")
         drhaux(:)=cmplx(real(fkdwfg(ig5,:)),-aimag(fkdwfg(ig5,:)))
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3) - vdotw(qplg(ig,:),drhaux,Cryst%gmet,"G")+lnkp(ig4)
       endif
     endif
   enddo
 enddo

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=10) then
     write(msg,'(a,i5)')' Number of G1-G2 pairs outside the G-sphere for Wfns = ',outofbox
     MSG_WARNING(msg)
     if (enough==10) then
       write(msg,'(a)')' ========== Stop writing Warnings =========='
       call wrtout(std_out,msg,'COLL')
     end if
   end if
 end if

 do ibv = 1, nbmax
   do igp=1,Sigp%npwc
     do ig=1,igp
       ptwsq(ig,igp,1)=ptwsq(ig,igp,1)-conjg(rhotwg(ig,ibv))*rhotwg(igp,ibv)
     enddo
   end do !igp
 end do !ig

 do igp = 1, Sigp%npwc
   do ig = igp+1, Sigp%npwc
     ptwsq(ig,igp,1)=conjg(ptwsq(igp,ig,1))
   enddo
 enddo

 if (niter>0) then
   ABI_ALLOCATE(paux,(Sigp%npwc,Sigp%npwc))
   paux(:,:)=(0.0,0.0)
   minusone=(-1.,0.)
   do ibv = 1, nbmax
     call XGERC(Sigp%npwc,Sigp%npwc,minusone,conjg(rhotwg(:,ibv)),1,conjg(cauxg(:,ibv)),1,ptwsq(:,:,2),Sigp%npwc)
   enddo
   do igp=1,Sigp%npwc
     do ig=1,Sigp%npwc
       if (ig>=igp) then
         do iat = 1, Cryst%natom
           do ilm = 1, nlx*nlx
             if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
             if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
             paux(ig,igp)=paux(ig,igp)+conjg(fnltwg(ig,ilm,iat))*fnltwg(igp,ilm,iat)
           enddo
         enddo
       endif
       ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+ptwsq(ig,igp,1)*kplqg(igp)
     end do !igp
   end do !ig
   do ig=1,Sigp%npwc
     do igp=1,Sigp%npwc
       if (ig>=igp) then
         ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+paux(ig,igp)
       else
         ptwsq(ig,igp,2)=ptwsq(ig,igp,2)+conjg(paux(igp,ig))
       endif
     end do !igp
   end do !ig
   ABI_DEALLOCATE(paux)
 endif

 if (niter>1) then

   !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz
   call nullify_kb_potential(KBff_kmq_ibz)

   call init_kb_potential(KBff_kmq_ibz,Cryst,Psps,2,istwf_kmq,npw_kmq,Kmesh%ibz(:,ikmq_ibz),kg_kmq)
   ABI_DEALLOCATE(KBff_kmq_ibz%fnld)

   ff(:,:,:,:)=(0.0,0.0)
   do iat = 1, Cryst%natom
     do ilm = 1, nlx*nlx
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       do iat2 = 1, Cryst%natom
         do ilm2 = 1, nlx*nlx
           if (ilm2>fnlmax(Cryst%typat(iat2))) CYCLE
           if (ilm2>=fnlloc(Cryst%typat(iat2),1).and.ilm2<=fnlloc(Cryst%typat(iat2),2)) CYCLE
           do ig = 1, Sigp%npwwfn
             igbz = grottbm1(ig,itim_kmq,isym_kmq)
             if (itim_kmq==1) then
               ff(ilm,iat,ilm2,iat2)=ff(ilm,iat,ilm2,iat2)+ &
                conjg(KBff_kmq_ibz%fnl(igbz,ilm,iat))*KBff_kmq_ibz%fnl(igbz,ilm2,iat2)
             else
               ff(ilm,iat,ilm2,iat2)=ff(ilm,iat,ilm2,iat2)+ &
                KBff_kmq_ibz%fnl(igbz,ilm,iat)*conjg(KBff_kmq_ibz%fnl(igbz,ilm2,iat2))
             endif
           enddo
         enddo
       enddo
     enddo
   enddo

   call destroy_kb_potential(KBff_kmq_ibz)

   vzn(:,:,:)=(0.0,0.0)
   do igp=1,Sigp%npwc
     do iat = 1, Cryst%natom
       do ilm = 1, nlx*nlx
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         vzn(igp,:,:) = vzn(igp,:,:) + half*ff(ilm,iat,:,:)*fnltwg(igp,ilm,iat)
       enddo
     enddo
     do iat = 1, Cryst%natom
       do ilm = 1, nlx*nlx
         if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
         if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
         drhaux(:)=cmplx(real(fdrhotwg(igp,ilm,iat,:)),aimag(fdrhotwg(igp,ilm,iat,:)))
         vzn(igp,ilm,iat)=vzn(igp,ilm,iat)+vdotw(qplg(igp,:),drhaux,Cryst%gmet,"G")-kns(igp,ilm,iat)
       enddo
     enddo
     do ig=1,igp
       do ibv = 1, nbmax
         ptwsq(ig,igp,3)=ptwsq(ig,igp,3)-conjg(cauxg(ig,ibv))*cauxg(igp,ibv)
       enddo
       ptwsq(ig,igp,3)=ptwsq(ig,igp,3)+kplqg(igp)*ptwsq(ig,igp,2)+kplqg(ig)*conjg(ptwsq(igp,ig,2))&
&                                     -kplqg(ig)*ptwsq(ig,igp,1)*kplqg(igp)
       do iat = 1, Cryst%natom
         do ilm = 1, nlx*nlx
           if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
           if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
           ptwsq(ig,igp,3)=ptwsq(ig,igp,3)+conjg(fnltwg(ig,ilm,iat))*vzn(igp,ilm,iat)+ &
&                                                fnltwg(igp,ilm,iat)*conjg(vzn(ig,ilm,iat))
         enddo
       enddo
     end do !igp
   end do !ig

   do igp = 1, Sigp%npwc
     do ig = igp+1, Sigp%npwc
       ptwsq(ig,igp,3)=conjg(ptwsq(igp,ig,3))
     enddo
   enddo

 endif

 ABI_DEALLOCATE(rhotwg)
 ABI_DEALLOCATE(wfwfg)
 if (niter>0) then
   ABI_DEALLOCATE(cauxg)
   ABI_DEALLOCATE(drhotwg)
   ABI_DEALLOCATE(dwfwfg)
   ABI_DEALLOCATE(fnltwg)
   ABI_DEALLOCATE(fnltwg2)
   ABI_DEALLOCATE(fnltwg3)
   ABI_DEALLOCATE(fnlwfg)
 endif
 if (niter>1) then
   ABI_DEALLOCATE(ddwfwfg)
   ABI_DEALLOCATE(lnkp)
   ABI_DEALLOCATE(kns)
   ABI_DEALLOCATE(fdrhotwg)
   ABI_DEALLOCATE(fkdwfg)
   ABI_DEALLOCATE(vzn)
   ABI_DEALLOCATE(ff)
 endif

end subroutine fft4eet_sig_kb

!!***
 
!----------------------------------------------------------------------

!!****f* ABINIT/calc_eet_sig_prep
!! NAME                  
!! calc_eet_sig_prep
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine calc_eet_sig_prep(Sigp,Cryst,Wfs,Kmesh,Psps,is,nbmax,ib1,ib2,ik_ibz, &
&                            jk_ibz,nspinor,fnlloc,fnlmax,fnlkr,fnlkpr,mtwk,mtwkp)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_errors

 use m_oscillators,   only : rho_tw_g
 use m_geometry,      only : normv
 use m_vcoul,         only : vcoul_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_eet_sig_prep'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(in) :: Wfs
 type(Sigma_parameters),intent(in) :: Sigp
 type(kb_potential) :: KBff_ki,KBff_kj

 integer,intent(in) :: is,nbmax,ib1,ib2
 integer,intent(in) :: nspinor
 integer,intent(in) :: ik_ibz,jk_ibz
 integer,intent(in) :: fnlloc(Cryst%ntypat,2)
 integer,intent(in) :: fnlmax(Cryst%ntypat)

 complex(gwpc),intent(in) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(in) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)

 complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,ib1:ib2)

!Local variables-------------------------------
!scalars
 integer :: ibv,kb,ilm,iat,ig
 integer :: istwf_ki,istwf_kj,npw_ki,npw_kj
!arrays 
 integer,pointer :: kg_ki(:,:),kg_kj(:,:)
 complex(gwpc),allocatable :: maux(:,:)

!************************************************************************

 !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz
 call nullify_kb_potential(KBff_ki  )
 call nullify_kb_potential(KBff_kj)

 istwf_ki = Wfs%istwfk(ik_ibz)
 istwf_kj = Wfs%istwfk(jk_ibz)

 npw_ki = Wfs%npwarr(ik_ibz)
 npw_kj = Wfs%npwarr(jk_ibz)

 kg_ki => Wfs%Kdata(ik_ibz)%kg_k
 kg_kj => Wfs%Kdata(jk_ibz)%kg_k

 call init_kb_potential(KBff_ki,Cryst,Psps,2,istwf_ki,npw_ki,Kmesh%ibz(:,ik_ibz),kg_ki)
 call init_kb_potential(KBff_kj,Cryst,Psps,2,istwf_kj,npw_kj,Kmesh%ibz(:,jk_ibz),kg_kj)
 ABI_DEALLOCATE(KBff_ki%fnld)
 ABI_DEALLOCATE(KBff_kj%fnld)

 ABI_ALLOCATE(maux,(Psps%mpsang*Psps%mpsang,Cryst%natom))

 mtwk(:,:)=(0.0,0.0)
 do ibv = 1, nbmax
   maux(:,:)=(0.0,0)
   do ig = 1, npw_ki
     maux(:,:) = maux(:,:) + Wfs%Wave(ibv,ik_ibz,is)%ug(ig)*KBff_ki%fnl(ig,:,:)
   enddo
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       mtwk(:,ibv)=mtwk(:,ibv)+maux(ilm,iat)*fnlkr(:,ilm,iat)
     enddo
   enddo
 enddo

 mtwkp(:,:)=(0.0,0.0)
 do kb = ib1, ib2
   maux(:,:)=(0.0,0)
   do ig = 1, npw_kj
     maux(:,:) = maux(:,:) + Wfs%Wave(kb,jk_ibz,is)%ug(ig)*KBff_kj%fnl(ig,:,:)
   enddo
   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE
       mtwkp(:,kb)=mtwkp(:,kb)+maux(ilm,iat)*fnlkpr(:,ilm,iat)
     enddo
   enddo
 enddo

 ABI_DEALLOCATE(maux)

 call destroy_kb_potential(KBff_ki)
 call destroy_kb_potential(KBff_kj)

 RETURN 
 ABI_UNUSED(Sigp%npwc)

end subroutine calc_eet_sig_prep
!!***

!!****f* ABINIT/check_delta_sigma
!!
!! NAME
!! check_delta_sigma
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (AB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (to be filled)
!!
!! OUTPUT
!!  (to be filled)
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine check_delta_sigma(qpgsq,qpgpsq,delta,omegame0k,omegame0lumo,ig,igp,gw_eet_scale,niter)

 use m_profiling

 use defs_basis
 use m_errors
 use m_gwdefs
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_delta_sigma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 integer,intent(in) :: ig,igp,niter
 real(dp),intent(in) :: qpgsq,qpgpsq
 real(dp),intent(in) :: omegame0k,omegame0lumo
 real(dp),intent(in) :: gw_eet_scale
 complex(gwpc),intent(inout) :: delta

!Local variables-------------------------------
!scalars

 real(dp) :: test

!*************************************************************************

 if (gw_eet_scale>0.01) then
   delta=gw_eet_scale*delta
 endif

 if (ig==igp) then
   delta=real(delta)
   test=omegame0k-real(delta)
   if (test>omegame0lumo.and.niter>0) then
     delta=half*(qpgsq+qpgpsq)
     test=omegame0k-real(delta)
   endif
   if (test>omegame0lumo) then
     delta=omegame0k-omegame0lumo
   endif
 else
   test=omegame0k-real(delta)
   if (test>omegame0lumo) then
     delta=omegame0k-omegame0lumo
   endif
 endif

end subroutine check_delta_sigma
!!***

!!****f* ABINIT/gw_eet_sigma_vkb
!! NAME
!! gw_eet_sigma_vkb
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine gw_eet_sigma_vkb(Sigp,Cryst,Wfs,Kmesh,Psps,isppol,ik_ibz,jk_ibz,ib1,ib2,nspinor, &
&                           tim_fourdp,MPI_enreg,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_errors
 use m_ppmodel

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_eet_sigma_vkb'
 use interfaces_70_gw, except_this_one => gw_eet_sigma_vkb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Sigma_parameters),intent(in) :: Sigp
 type(wfs_descriptor),intent(inout) :: Wfs
 type(MPI_type),intent(inout) :: MPI_enreg

 type(kb_potential) :: KBff_ki,KBff_kj

 integer,intent(in) :: tim_fourdp

 integer,intent(in) :: isppol
 integer,intent(in) :: nspinor
 integer,intent(in) :: ik_ibz,jk_ibz,ib1,ib2
 integer,intent(in) :: nbmax

 integer,intent(out) :: fnlloc(Cryst%ntypat,2),fnlmax(Cryst%ntypat)

 complex(gwpc),intent(out) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(out) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbmax)
 complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,ib1:ib2)

 integer :: istwf_ki,istwf_kj,npw_ki,npw_kj
 integer :: iat,ilm,i

 integer,pointer :: gbound_ki(:,:),gbound_kj(:,:),igfft0_ki(:),igfft0_kj(:)

 character(len=fnlen) :: title
 integer :: lloc,lmax,mmax,pspcod,pspdat,pspxc
 integer :: ityp
 real(dp) :: r2well,zion,znucl

 integer,pointer :: kg_ki(:,:),kg_kj(:,:)

!************************************************************************

 fnlloc(:,:)=0
 fnlmax(:)=0
 do ityp = 1, Cryst%ntypat
   open (unit=tmp_unit,file=psps%filpsp(ityp),form='formatted',status='old')
   rewind (unit=tmp_unit)
   read (tmp_unit,'(a)') title
   read (tmp_unit,*) znucl,zion,pspdat
   read (tmp_unit,*) pspcod,pspxc,lmax,lloc,mmax,r2well
   do i = 1, lloc
     fnlloc(ityp,1) = fnlloc(ityp,1) + 2*(lloc-1)+1
   enddo
   fnlloc(ityp,1)=fnlloc(ityp,1)+1
   do i = 1, lloc+1
     fnlloc(ityp,2) = fnlloc(ityp,2) + 2*(lloc-1)+1
   enddo
   do i = 1, lmax+1
     fnlmax(ityp) = fnlmax(ityp) + 2*(lmax-1)+1
   enddo
 enddo

 istwf_ki  =  Wfs%istwfk(ik_ibz)
 npw_ki    =  Wfs%npwarr(ik_ibz)
 kg_ki     => Wfs%Kdata(ik_ibz)%kg_k
 igfft0_ki => Wfs%Kdata(ik_ibz)%igfft0
 gbound_ki => Wfs%Kdata(ik_ibz)%gbound
 
 istwf_kj  =  Wfs%istwfk(jk_ibz)
 npw_kj    =  Wfs%npwarr(jk_ibz)
 kg_kj     => Wfs%Kdata(jk_ibz)%kg_k
 igfft0_kj => Wfs%Kdata(jk_ibz)%igfft0
 gbound_kj => Wfs%Kdata(jk_ibz)%gbound

 !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz
 call nullify_kb_potential(KBff_ki)
 call nullify_kb_potential(KBff_kj)

 call init_kb_potential(KBff_ki,Cryst,Psps,2,istwf_ki,npw_ki,Kmesh%ibz(:,ik_ibz),kg_ki)
 call init_kb_potential(KBff_kj,Cryst,Psps,2,istwf_kj,npw_kj,Kmesh%ibz(:,jk_ibz),kg_kj)
 ABI_DEALLOCATE(KBff_ki%fnld)
 ABI_DEALLOCATE(KBff_kj%fnld)

 !MG: note that it is not necessary to rotate fnl in the full BZ to
 !calculate
 !<SK|V_nl|Sk>.
 !MG check this part in parallel. Keep in mind the difference between
 !Wfs_braket (bdgw states= and Wfs (full set distributed
 !across the nod

 ! MG: TODO Be careful here as the symmetry properties of fnl are not the same as the 
 !     ones of the wavefunctions. One should write a separate methods for the FFT of fnl.
 if (istwf_kj>1.or.istwf_ki>1) then
   MSG_ERROR("istwfk /= 1 not coded")
 end if

 do iat = 1, Cryst%natom
   do ilm = 1, Psps%mpsang*Psps%mpsang
     if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
     if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE

     call fft_onewfn(Wfs%paral_kgb,istwf_ki,nspinor,npw_ki,Wfs%nfftot,Wfs%mgfft,Wfs%ngfft,&
&                    conjg(KBff_ki%fnl(:,ilm,iat)),fnlkr(:,ilm,iat),igfft0_ki, &
&                    kg_ki,gbound_ki,tim_fourdp,MPI_enreg)


     call fft_onewfn(Wfs%paral_kgb,istwf_kj,nspinor,npw_kj,Wfs%nfftot,Wfs%mgfft,Wfs%ngfft,&
&                    conjg(KBff_kj%fnl(:,ilm,iat)),fnlkpr(:,ilm,iat),igfft0_kj, &
&                    kg_kj,gbound_kj,tim_fourdp,MPI_enreg)
   enddo
 enddo
 call destroy_kb_potential(KBff_ki)
 call destroy_kb_potential(KBff_kj)

 call calc_eet_sig_prep(Sigp,Cryst,Wfs,Kmesh,Psps,isppol,nbmax,ib1,ib2,ik_ibz,jk_ibz, &
&                       nspinor,fnlloc,fnlmax,fnlkr,fnlkpr,mtwk,mtwkp)

end subroutine gw_eet_sigma_vkb
!!***

!!****f* ABINIT/gw_eet_chi0_vkb
!! NAME
!! gw_eet_chi0_vkb
!!
!! FUNCTION
!! (to be provided)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_sig_ppm_eet
!!
!! CHILDREN
!!      calc_eet_prep,destroy_kb_potential,fft_onewfn,init_kb_potential
!!      nullify_kb_potential
!!
!! SOURCE

subroutine gw_eet_chi0_vkb(Ep,Cryst,Wfs,Kmesh,Psps,is,ik_ibz,ikmq_ibz,nspinor,tim_fourdp, &
&                          nbhomo,nbmax,fnlkr,fnlkpr,mtwk,mtwkp,fnlloc,fnlmax)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_crystal
 use m_bz_mesh
 use m_gsphere
 use m_wfs
 use m_commutator_vkbr
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_eet_chi0_vkb'
 use interfaces_70_gw, except_this_one => gw_eet_chi0_vkb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Pseudopotential_type),intent(in) :: Psps
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(wfs_descriptor),intent(inout) :: Wfs
 type(kb_potential) :: KBff_k_ibz,KBff_kmq_ibz

 integer,intent(in) :: nspinor,tim_fourdp
 integer,intent(in) :: is
 integer,intent(in) :: nbhomo(2), nbmax
 integer,intent(in) :: ik_ibz,ikmq_ibz

 integer,intent(out) :: fnlloc(Cryst%ntypat,2),fnlmax(Cryst%ntypat)

 complex(gwpc),intent(out) :: fnlkr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(out) :: fnlkpr(Wfs%nfftot*nspinor,Psps%mpsang*Psps%mpsang,Cryst%natom)
 complex(gwpc),intent(out) :: mtwk(Wfs%nfftot*nspinor,nbhomo(1))
 complex(gwpc),intent(out) :: mtwkp(Wfs%nfftot*nspinor,nbmax)

 integer :: npw_k,npw_kmq
 integer :: i,ilm,iat
 integer :: istwf_k,istwf_kmq

 character(len=fnlen) :: title
 integer :: lloc,lmax,mmax,pspcod,pspdat,pspxc
 integer :: ityp
 integer,pointer :: gbound_k(:,:),gbound_kmq(:,:),kg_k(:,:),kg_kmq(:,:)
 integer,pointer :: igfft0_k(:),igfft0_kmq(:)
 real(dp) :: r2well,zion,znucl

!************************************************************************

   fnlloc(:,:)=0
   fnlmax(:)=0
   do ityp = 1, Cryst%ntypat
     open (unit=tmp_unit,file=psps%filpsp(ityp),form='formatted',status='old')
     rewind (unit=tmp_unit)
     read (tmp_unit,'(a)') title
     read (tmp_unit,*) znucl,zion,pspdat
     read (tmp_unit,*) pspcod,pspxc,lmax,lloc,mmax,r2well
     do i = 1, lloc
       fnlloc(ityp,1) = fnlloc(ityp,1) + 2*(lloc-1)+1
     enddo
     fnlloc(ityp,1)=fnlloc(ityp,1)+1
     do i = 1, lloc+1
       fnlloc(ityp,2) = fnlloc(ityp,2) + 2*(lloc-1)+1
     enddo
     do i = 1, lmax+1
       fnlmax(ityp) = fnlmax(ityp) + 2*(lmax-1)+1
     enddo
   enddo

   !MG KB form factors: to be done outside. array dimensioned with Kmesh%nibz
   call nullify_kb_potential(KBff_k_ibz  )
   call nullify_kb_potential(KBff_kmq_ibz)

   npw_k    = Wfs%npwarr(ik_ibz)
   istwf_k  = Wfs%istwfk(ik_ibz)
   kg_k     => Wfs%Kdata(ik_ibz)%kg_k
   igfft0_k => Wfs%Kdata(ik_ibz)%igfft0
   gbound_k => Wfs%Kdata(ik_ibz)%gbound

   npw_kmq   = Wfs%npwarr(ikmq_ibz)
   istwf_kmq = Wfs%istwfk(ikmq_ibz)
   kg_kmq     => Wfs%Kdata(ikmq_ibz)%kg_k
   igfft0_kmq => Wfs%Kdata(ikmq_ibz)%igfft0
   gbound_kmq => Wfs%Kdata(ikmq_ibz)%gbound

   call init_kb_potential(KBff_k_ibz  ,Cryst,Psps,2,istwf_k  ,npw_k,  Kmesh%ibz(:,ik_ibz),kg_k)
   call init_kb_potential(KBff_kmq_ibz,Cryst,Psps,2,istwf_kmq,npw_kmq,Kmesh%ibz(:,ikmq_ibz),kg_kmq)

   ABI_DEALLOCATE(KBff_k_ibz%fnld)
   ABI_DEALLOCATE(KBff_kmq_ibz%fnld)

   ! MG: TODO Be careful here as the symmetry properties of fnl are not the same as the 
   !     ones of the wavefunctions. One should write a separate methods for the FFT of fnl.
   if (istwf_k>1.or.istwf_kmq>1) then
     MSG_ERROR("istwfk /= 1 not coded")
   end if

   do iat = 1, Cryst%natom
     do ilm = 1, Psps%mpsang*Psps%mpsang
       if (ilm>fnlmax(Cryst%typat(iat))) CYCLE
       if (ilm>=fnlloc(Cryst%typat(iat),1).and.ilm<=fnlloc(Cryst%typat(iat),2)) CYCLE

       call fft_onewfn(Wfs%paral_kgb,istwf_k,nspinor,npw_k,Wfs%nfftot,Wfs%mgfft,Wfs%ngfft,&
&                      conjg(KBff_k_ibz%fnl(:,ilm,iat)),fnlkr(:,ilm,iat),igfft0_k,&
&                      kg_k,gbound_k,tim_fourdp,Wfs%MPI_enreg)

       call fft_onewfn(Wfs%paral_kgb,istwf_kmq,nspinor,npw_kmq,Wfs%nfftot,Wfs%mgfft,Wfs%ngfft,&
&                      conjg(KBff_kmq_ibz%fnl(:,ilm,iat)),fnlkpr(:,ilm,iat),igfft0_kmq,&
&                      kg_kmq,gbound_kmq,tim_fourdp,Wfs%MPI_enreg)
     enddo
   enddo

   call destroy_kb_potential(KBff_k_ibz  )
   call destroy_kb_potential(KBff_kmq_ibz)

   call calc_eet_prep(Ep,Cryst,Wfs,Kmesh,Psps,is,nbhomo(1),nbmax,ik_ibz,ikmq_ibz,nspinor, &
&                     fnlloc,fnlmax,fnlkr,fnlkpr,mtwk,mtwkp)

end subroutine gw_eet_chi0_vkb

!!***
