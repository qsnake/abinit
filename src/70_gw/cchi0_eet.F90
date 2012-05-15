!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_delta_chi0
!!
!! NAME
!! calc_delta_chi0
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
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calc_delta_chi0(Ep,frhorho,frhoj,fjj,niter)

 use m_profiling

 use defs_basis
 use m_errors
 use m_gwdefs
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_delta_chi0'
!End of the abilint section

 implicit none

 type(Epsilonm1_parameters),intent(in) :: Ep

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter
!arrays
 complex(gwpc),intent(in) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(inout) :: frhoj(Ep%npwe,Ep%npwe)
 complex(gwpc),intent(inout) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))

!Local variables-------------------------------
!scalars

 integer :: ig,igp,igaux!,ios
! integer :: iter

! real(dp) :: bq
 real(dp) :: delta_huge,denchk

 complex(gwpc) :: num,den

!*************************************************************************

   delta_huge = 1.0d8

   if (niter==1) then

     do ig=1,Ep%npwe
       do igp=ig,Ep%npwe
         igaux=igp*(igp-1)/2+ig
         num = frhoj(ig,igp)+conjg(frhoj(igp,ig))
         den = 2*frhorho(igaux)
         denchk=abs(real(den))+abs(aimag(den))
         if (denchk>0.0) then
           frhoj(ig,igp)=num/den
         else
           frhoj(ig,igp)=cmplx(delta_huge,0.0)
         endif
       end do !igp
     end do !ig

     do ig=1,Ep%npwe
       do igp=1,ig-1
         frhoj(ig,igp)=conjg(frhoj(igp,ig))
       end do !igp
     end do !ig

   elseif (niter==2) then

     do ig=1,Ep%npwe
       do igp=ig,Ep%npwe
         igaux=igp*(igp-1)/2+ig
         num = 2*fjj(igaux)
         den = frhoj(ig,igp)+conjg(frhoj(igp,ig))
         denchk=abs(real(den))+abs(aimag(den))
         if (denchk>0.0) then
           fjj(igaux)=num/den
         else
           fjj(igaux)=cmplx(delta_huge,0.0)
         endif
       end do !igp
     end do !ig

     do ig=1,Ep%npwe
       do igp=ig,Ep%npwe
         igaux=igp*(igp-1)/2+ig
         num = frhoj(ig,igp)+conjg(frhoj(igp,ig))
         den = 2*frhorho(igaux)
         denchk=abs(real(den))+abs(aimag(den))
         if (denchk>0.0) then
           frhoj(ig,igp)=num/den
         else
           frhoj(ig,igp)=cmplx(delta_huge,0.0)
         endif
       end do !igp
     end do !ig

     do ig=1,Ep%npwe
       do igp=1,ig-1
         frhoj(ig,igp)=conjg(frhoj(igp,ig))
       end do !igp
     end do !ig

   endif

end subroutine calc_delta_chi0
!!***

!!****f* ABINIT/calc_chi0_delta_clos
!!
!! NAME
!! calc_chi0_delta_clos
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
!!
!! SOURCE

subroutine calc_chi0_delta_clos(ik_bz,Dtset,Ep,Gsph_epsG0,Ltg_q,niter,epsv,epslumo,qpgsq,frhorho,frhoj,fjj,chi0)

 use m_profiling

 use defs_abitypes
 use defs_basis
 use m_errors
 use m_bz_mesh,  only : little_group
 use m_gwdefs,   only : epsilonm1_parameters
 use m_gsphere,  only : gvectors_type
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_chi0_delta_clos'
 use interfaces_70_gw, except_this_one => calc_chi0_delta_clos
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars

 type(Dataset_type),intent(in) :: Dtset
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Little_group),intent(in) :: Ltg_q

 integer,intent(in) :: ik_bz
 integer,intent(in) :: niter
!arrays
 real(dp),intent(in) :: qpgsq(Ep%npwe)
 real(dp),intent(in) :: epsv,epslumo
 complex(gwpc),intent(in) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(in) :: frhoj(Ep%npwe,Ep%npwe)
 complex(gwpc),intent(in) :: fjj(Ep%npwe*(Ep%npwe+1)/2*(niter-1))
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,igaux,igaux2,ios,isym,itim
 character(len=500) :: msg

 integer,pointer :: gmG0(:)
 integer,allocatable :: Sm1_gmG0(:)
 complex(gwpc),pointer :: phmGt(:)

 real(dp) :: bq
 complex(gwpc) :: delta1,delta2
 complex(gwpc) :: faux
 complex(gwpc) :: frhorho_sym,frhoj_sym,fjj_sym

!*************************************************************************

 SELECT CASE (Ep%symchi)

 CASE (0)

   do ig=1,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       bq=half*(qpgsq(ig)+qpgsq(igp))
       faux=frhorho(igaux)

       SELECT CASE (niter)
       CASE (1)
         delta1=bq+frhoj(ig,igp)
         if (Dtset%gw_eet_scale>0.01) delta1=Dtset%gw_eet_scale*delta1
         delta2=delta1
         do ios=1,Ep%nomega
           if (abs(aimag(Ep%omega(ios)))<0.001) then
             call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
           else
             delta1=delta2
           endif
           chi0(ig,igp,ios) = chi0(ig,igp,ios) - faux/(Ep%omega(ios)+delta1) &
                                               + faux/(Ep%omega(ios)-delta1)
         end do !ios
       CASE (2)
         do ios=1,Ep%nomega
           delta1=bq+frhoj(ig,igp)*(Ep%omega(ios)+bq+frhoj(ig,igp))/(Ep%omega(ios)+bq+fjj(igaux))
           delta2=bq+frhoj(ig,igp)*(Ep%omega(ios)-bq-frhoj(ig,igp))/(Ep%omega(ios)-bq-fjj(igaux))
           if (Dtset%gw_eet_scale>0.01) then
             delta1=Dtset%gw_eet_scale*delta1
             delta2=Dtset%gw_eet_scale*delta2
           endif
           if (abs(aimag(Ep%omega(ios)))<0.001) then
             call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
             call check_delta(Ep,qpgsq,delta2,epsv,epslumo,ig,igp)
           endif
           chi0(ig,igp,ios) = chi0(ig,igp,ios) - faux/(Ep%omega(ios)+delta1) &
                                               + faux/(Ep%omega(ios)-delta2)
         end do !ios
       END SELECT

     end do !igp
   end do !ig

 CASE (1)

   ABI_ALLOCATE(Sm1_gmG0,(Ep%npwe))

   do isym=1,Ltg_q%nsym_sg
     do itim = 1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then

         phmGt => Gsph_epsG0%phmGt(1:Ep%npwe,isym)
         gmG0  => Ltg_q%igmG0(1:Ep%npwe,itim,isym)
         Sm1_gmG0(1:Ep%npwe)=Gsph_epsG0%rottbm1(gmG0(1:Ep%npwe),itim,isym)

         do ig = 1, Ep%npwe
           do igp = ig, Ep%npwe
             igaux=igp*(igp-1)/2+ig
             bq=half*(qpgsq(ig)+qpgsq(igp))
             if (Sm1_gmG0(ig)<=Sm1_gmG0(igp)) then
               igaux2=Sm1_gmG0(igp)*(Sm1_gmG0(igp)-1)/2+Sm1_gmG0(ig)
               frhorho_sym=frhorho(igaux2)
               if (niter==2) fjj_sym=fjj(igaux2)
             else
               igaux2=Sm1_gmG0(ig)*(Sm1_gmG0(ig)-1)/2+Sm1_gmG0(igp)
               frhorho_sym=conjg(frhorho(igaux2))
               if (niter==2) fjj_sym=conjg(fjj(igaux2))
             endif
             frhoj_sym=frhoj(Sm1_gmG0(ig),Sm1_gmG0(igp))
             if (itim==2) then
               frhorho_sym=conjg(frhorho_sym)
               frhoj_sym=conjg(frhoj_sym)
               if (niter==2) fjj_sym=conjg(fjj_sym)
             endif
             frhorho_sym=frhorho_sym*conjg(phmGt(igp))*phmGt(ig)

             SELECT CASE (niter)
             CASE (1)
               delta1=bq+frhoj_sym
               if (Dtset%gw_eet_scale>0.01) delta1=Dtset%gw_eet_scale*delta1
               delta2=delta1
               do ios=1,Ep%nomega
                 if (abs(aimag(Ep%omega(ios)))<0.001) then
                   call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
                 else
                   delta1=delta2
                 endif
                 chi0(ig,igp,ios) = chi0(ig,igp,ios) - frhorho_sym/(Ep%omega(ios)+delta1) &
                                                     + frhorho_sym/(Ep%omega(ios)-delta1)
               enddo
             CASE (2)
               do ios=1,Ep%nomega
                 delta1=bq+frhoj_sym*(Ep%omega(ios)+bq+frhoj_sym)/(Ep%omega(ios)+bq+fjj_sym)
                 delta2=bq+frhoj_sym*(Ep%omega(ios)-bq-frhoj_sym)/(Ep%omega(ios)-bq-fjj_sym)
                 if (Dtset%gw_eet_scale>0.01) then
                   delta1=Dtset%gw_eet_scale*delta1
                   delta2=Dtset%gw_eet_scale*delta2
                 endif
                 if (abs(aimag(Ep%omega(ios)))<0.001) then
                   call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
                   call check_delta(Ep,qpgsq,delta2,epsv,epslumo,ig,igp)
                 endif
                 chi0(ig,igp,ios) = chi0(ig,igp,ios) - frhorho_sym/(Ep%omega(ios)+delta1) & 
                                                     + frhorho_sym/(Ep%omega(ios)-delta2)
               enddo
             END SELECT

           end do !igp
         end do !ig

       endif
     enddo
   enddo

   ABI_DEALLOCATE(Sm1_gmG0)

 CASE DEFAULT
  write(msg,'(a,i3)')'Wrong symchi= ',Ep%symchi
  MSG_BUG(msg)
 END SELECT

end subroutine calc_chi0_delta_clos
!!***

!!****f* ABINIT/calc_chi0_delta0
!!
!! NAME
!! calc_chi0_delta0
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
!!
!! SOURCE

subroutine calc_chi0_delta0(ik_bz,Dtset,Ep,Gsph_epsG0,Ltg_q,epsv,epslumo,qpgsq,frhorho,chi0)

 use m_profiling

 use defs_abitypes
 use defs_basis
 use m_errors
 use m_bz_mesh,  only : little_group
 use m_gwdefs,   only : epsilonm1_parameters
 use m_gsphere,  only : gvectors_type
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_chi0_delta0'
 use interfaces_70_gw, except_this_one => calc_chi0_delta0
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars

 type(Dataset_type),intent(in) :: Dtset
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Little_group),intent(in) :: Ltg_q

 integer,intent(in) :: ik_bz
!arrays
 real(dp),intent(in) :: qpgsq(Ep%npwe)
 real(dp),intent(in) :: epsv,epslumo
 complex(gwpc),intent(in) :: frhorho(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ios,isym,itim,igaux,igaux2
 character(len=500) :: msg

 integer,pointer :: gmG0(:)
 integer,allocatable :: Sm1_gmG0(:)
 complex(gwpc),pointer :: phmGt(:)

 complex(gwpc) :: faux,delta1,delta2,frhorho_sym

!*************************************************************************

 SELECT CASE (Ep%symchi)

 CASE (0)

   do ig=1,Ep%npwe
     do igp=ig,Ep%npwe
       igaux=igp*(igp-1)/2+ig
       faux=frhorho(igaux)
       delta1=half*(qpgsq(ig)+qpgsq(igp))
       if (Dtset%gw_eet_scale>0.01) delta1=Dtset%gw_eet_scale*delta1
       delta2=delta1
       do ios=1,Ep%nomega
         if (abs(aimag(Ep%omega(ios)))<0.001) then
           call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
         else
           delta1=delta2
         endif
         chi0(ig,igp,ios) = chi0(ig,igp,ios) - faux/(Ep%omega(ios)+delta1) &
                                             + faux/(Ep%omega(ios)-delta1)
       end do !ios
     end do !igp
   end do !ig

 CASE (1)

   ABI_ALLOCATE(Sm1_gmG0,(Ep%npwe))

   do isym=1,Ltg_q%nsym_sg
     do itim = 1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then

         phmGt => Gsph_epsG0%phmGt(1:Ep%npwe,isym)
         gmG0  => Ltg_q%igmG0(1:Ep%npwe,itim,isym)
         Sm1_gmG0(1:Ep%npwe)=Gsph_epsG0%rottbm1(gmG0(1:Ep%npwe),itim,isym)

         do ig = 1, Ep%npwe
           do igp = ig, Ep%npwe
             igaux=igp*(igp-1)/2+ig
             if (Sm1_gmG0(ig)<=Sm1_gmG0(igp)) then
               igaux2=Sm1_gmG0(igp)*(Sm1_gmG0(igp)-1)/2+Sm1_gmG0(ig)
               faux=frhorho(igaux2)
             else
               igaux2=Sm1_gmG0(ig)*(Sm1_gmG0(ig)-1)/2+Sm1_gmG0(igp)
               faux=conjg(frhorho(igaux2))
             endif
             frhorho_sym=faux*conjg(phmGt(igp))*phmGt(ig)
             if (itim==2) frhorho_sym=conjg(frhorho_sym)
             delta1=half*(qpgsq(ig)+qpgsq(igp))
             if (Dtset%gw_eet_scale>0.01) delta1=Dtset%gw_eet_scale*delta1
             delta2=delta1
             do ios=1,Ep%nomega
               if (abs(aimag(Ep%omega(ios)))<0.001) then
                 call check_delta(Ep,qpgsq,delta1,epsv,epslumo,ig,igp)
               else
                 delta1=delta2
               endif
               chi0(ig,igp,ios) = chi0(ig,igp,ios) - frhorho_sym/(Ep%omega(ios)+delta1) &
                                                   + frhorho_sym/(Ep%omega(ios)-delta1)
             end do !ios
           end do !igp
         end do !ig

       endif
     enddo
   enddo

   ABI_DEALLOCATE(Sm1_gmG0)

 CASE DEFAULT
  write(msg,'(a,i3)')'Wrong symchi= ',Ep%symchi
  MSG_BUG(msg)
 END SELECT

end subroutine calc_chi0_delta0
!!***

!!****f* ABINIT/calc_corr_chi0
!!
!! NAME
!! calc_corr_chi0
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
!!
!! SOURCE

subroutine calc_corr_chi0(ik_bz,Ep,Kmesh,Gsph_epsG0,Ltg_q,rhotwg,spin_fact,qp_occ,qp_energy, &
&                         ibv,ibc,ik_ibz,ikmq_ibz,is,chi0)

 use m_profiling

 use defs_basis
 use m_bz_mesh
 use m_gwdefs
 use m_gsphere,  only : gvectors_type
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_corr_chi0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Little_group),intent(in) :: Ltg_q

 integer,intent(in) :: ik_bz,ibv,ibc,ik_ibz,ikmq_ibz,is
 real(dp),intent(in) :: spin_fact
 real(dp),intent(in) :: qp_occ(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 real(dp),intent(in) :: qp_energy(Ep%nbnds,Kmesh%nibz,Ep%nsppol)
 complex(gwpc),intent(in) :: rhotwg(Ep%npwepG0)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

 integer :: ig,igp,ios,isym,itim
 real(dp) :: deltae,deltaf
 complex(dpc) :: green_w

 character(len=500) :: msg
!arrays
 integer,pointer :: gmG0(:)
 integer,allocatable :: Sm1_gmG0(:)
 complex(gwpc),allocatable :: rhotwg_sym(:)
 complex(gwpc),pointer :: phmGt(:)

!************************************************************************

 deltaf=spin_fact*(qp_occ(ibc,ikmq_ibz,is)-qp_occ(ibv,ik_ibz,is))
 deltae=qp_energy(ibc,ikmq_ibz,is)-qp_energy(ibv,ik_ibz,is)

 SELECT CASE (Ep%symchi)

 CASE (0) ! Do not use symmetries

   do ios = 1, Ep%nomega
     if (ibc==ibv) then
       green_w = g0g0w(Ep%omega(ios),deltaf,deltae,Ep%zcut,GW_TOL_W0,1)
     else
       green_w = g0g0w(Ep%omega(ios),deltaf,deltae,Ep%zcut,GW_TOL_W0,2)
     endif
     if (ABS(REAL(Ep%omega(ios)))<0.00001) then
       do ig=1,Ep%npwe
         do igp=ig,Ep%npwe
           chi0(ig,igp,ios) = chi0(ig,igp,ios)+rhotwg(ig)*conjg(rhotwg(igp))*green_w
         end do !igp
       end do !ig
     else
       do ig=1,Ep%npwe
         do igp=1,Ep%npwe
           chi0(ig,igp,ios) = chi0(ig,igp,ios)+rhotwg(ig)*conjg(rhotwg(igp))*green_w
         end do !igp
       end do !ig
     endif
   end do !ios

 CASE (1) ! Use symmetries to reconstruct the integrand in the BZ.

   ABI_ALLOCATE(rhotwg_sym,(Ep%npwe))
   ABI_ALLOCATE(Sm1_gmG0  ,(Ep%npwe))

   do isym=1,Ltg_q%nsym_sg
     do itim=1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then

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

        do ios=1, Ep%nomega
          if (ibc==ibv) then
            green_w = g0g0w(Ep%omega(ios),deltaf,deltae,Ep%zcut,GW_TOL_W0,1)
          else
            green_w = g0g0w(Ep%omega(ios),deltaf,deltae,Ep%zcut,GW_TOL_W0,2)
          endif
          if (ABS(REAL(Ep%omega(ios)))<0.00001) then
            do ig=1,Ep%npwe
              do igp=ig,Ep%npwe
                chi0(ig,igp,ios) = chi0(ig,igp,ios)+rhotwg_sym(ig)*conjg(rhotwg_sym(igp))*green_w
              end do !igp
            end do !ig
          else
            do ig=1,Ep%npwe
              do igp=1,Ep%npwe
                chi0(ig,igp,ios) = chi0(ig,igp,ios)+rhotwg_sym(ig)*conjg(rhotwg_sym(igp))*green_w
              end do !igp
            end do !ig
          endif
        end do

       end if
     end do
   end do

   ABI_DEALLOCATE(rhotwg_sym)
   ABI_DEALLOCATE(Sm1_gmG0)

 CASE DEFAULT
   write(msg,'(a,i3)')'Wrong symchi= ',Ep%symchi
   MSG_BUG(msg)
 END SELECT

end subroutine calc_corr_chi0
!!***

!!****f* ABINIT/calc_corr_sig
!!
!! NAME
!! calc_corr_sig
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
!!
!! SOURCE

subroutine calc_corr_sig(Sigp,Sr,nomega,nspinor,npwc1,npwc2,botsq,otq,rhotwg,is,ibv,kb, &
&                        ik_bz,ikmq_bz,ik_ibz,ikmq_ibz,i_sz,vc_sqrt_qbz,sigmac)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes

 use m_sigma_results,   only : sigma_results

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_corr_sig'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Sigma_parameters),intent(in) :: Sigp
 type(Sigma_results),intent(in) :: Sr

 integer,intent(in) :: nomega,nspinor,npwc1,npwc2
 integer,intent(in) :: is,ibv,kb
 integer,intent(in) :: ik_bz,ikmq_bz,ik_ibz,ikmq_ibz
 real(dp),intent(in) :: i_sz
 complex(gwpc),intent(in) :: rhotwg(Sigp%npwc)
 complex(gwpc),intent(in) :: vc_sqrt_qbz(Sigp%npwc)
 complex(gwpc),intent(in) :: otq(Sigp%npwc,npwc2)
 complex(gwpc),intent(in) :: botsq(Sigp%npwc,npwc1)
 complex(dpc),intent(inout) :: sigmac(nomega)

 complex(gwpc),allocatable :: rtaux(:)

 integer :: ig,igp,ios
 real(dp) :: otw,twofm1_zcut
 real(dp) :: den,omegame0i
 complex(gwpc) :: num

!************************************************************************

 ABI_ALLOCATE(rtaux,(Sigp%npwc*nspinor**2))

 twofm1_zcut=-Sigp%zcut

 do ig = 1,Sigp%npwc
   rtaux(ig)=rhotwg(ig)*vc_sqrt_qbz(ig)
 enddo
 if (ik_bz==ikmq_bz) then
   rtaux(1)=czero_gw
   if (ibv==kb) then
     rtaux(1)=cmplx(sqrt(i_sz),0.0_gwp)
   endif
 endif
 do ios=1,nomega
   omegame0i = real(Sr%omega4sd(kb,ikmq_ibz,ios,is)) - Sr%e0(ibv,ik_ibz,is)
   do ig=1,Sigp%npwc
     do igp=1,Sigp%npwc
       otw = DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
       num = botsq(ig,igp)*conjg(rtaux(ig))*rtaux(igp)
       den = omegame0i-otw
       if (real(den*den)>Sigp%zcut**2) then
         sigmac(ios) = sigmac(ios) + 0.5*num/(den*otw)
       else
         sigmac(ios) = sigmac(ios) + 0.5*num*cmplx(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw)
       end if
     end do !igp
   end do !ig
 end do !ios

 ABI_DEALLOCATE(rtaux)

end subroutine calc_corr_sig
!!***

!!****f* ABINIT/calc_delta0
!!
!! NAME
!! calc_delta0
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
!!
!! CHILDREN    
!!
!! SOURCE      

subroutine calc_delta0(Dtset,Ep,qpgsq,delta)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_delta0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Dataset_type),intent(in) :: Dtset
 type(Epsilonm1_parameters),intent(in) :: Ep

 real(dp), intent(in) :: qpgsq(Ep%npwe)
 complex(gwpc),intent(out) :: delta(Ep%npwe,Ep%npwe,Ep%nomega)

 integer :: ig,igp,ios

!************************************************************************

 do ig = 1, Ep%npwe
   do igp = ig, Ep%npwe
     delta(ig,igp,1)=half*(qpgsq(ig)+qpgsq(igp))
     if (Dtset%gw_eet_scale>0.01) then
       delta(ig,igp,1)=Dtset%gw_eet_scale*delta(ig,igp,1)
     endif
   enddo
 enddo
 do ig = 1, Ep%npwe
   do igp = 1, ig-1
     delta(ig,igp,1)=delta(igp,ig,1)
   enddo
 enddo
 do ios = 2, Ep%nomega
   delta(:,:,ios)=delta(:,:,1)
 enddo

end subroutine calc_delta0
!!***

!!****f* ABINIT/calc_chi0_delta0_bis
!!
!! NAME
!! calc_chi0_delta0_bis
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
!!
!! CHILDREN
!!
!! SOURCE

subroutine calc_chi0_delta0_bis(ik_bz,Ep,Gsph_epsG0,Ltg_q,paux,delta,chi0)

 use m_profiling

 use defs_basis
 use m_errors
 use m_bz_mesh,  only : little_group
 use m_gwdefs,   only : epsilonm1_parameters
 use m_gsphere,  only : gvectors_type
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_chi0_delta0_bis'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Epsilonm1_parameters),intent(in) :: Ep
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Little_group),intent(in) :: Ltg_q

 integer,intent(in) :: ik_bz

 complex(gwpc),intent(in) :: paux(Ep%npwe*(Ep%npwe+1)/2)
 complex(gwpc),intent(in) :: delta(Ep%npwe,Ep%npwe,Ep%nomega)
 complex(gwpc),intent(out) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

 integer :: ig,igp,ios,isym,itim,igaux,igaux2
 character(len=500) :: msg

 integer,pointer :: gmG0(:)
 integer,allocatable :: Sm1_gmG0(:)
 complex(gwpc),allocatable :: paux_sym(:)
 complex(gwpc),pointer :: phmGt(:)

 complex(gwpc) :: delta_sym
 complex(gwpc) :: faux

!************************************************************************

 SELECT CASE (Ep%symchi)

 CASE (0)

   do ios = 1, Ep%nomega
     if (ABS(REAL(Ep%omega(ios)))<0.00001) then
       if (abs(aimag(Ep%omega(ios)))>=0.001) then
         do ig = 1, Ep%npwe
           do igp = ig, Ep%npwe
             igaux=igp*(igp-1)/2+ig
             faux=paux(igaux)
             chi0(ig,igp,ios) = chi0(ig,igp,ios) - faux/(Ep%omega(ios)+delta(ig,igp,ios)) &
                                                 + faux/(Ep%omega(ios)-delta(ig,igp,ios))
           enddo
         enddo
       endif
     else
       if (abs(aimag(Ep%omega(ios)))>=0.001) then
         do ig = 1, Ep%npwe
           do igp = 1, Ep%npwe
             if (ig<=igp) then
               igaux=igp*(igp-1)/2+ig
               faux=paux(igaux)
             else
               igaux=ig*(ig-1)/2+igp
               faux=conjg(paux(igaux))
             endif
             chi0(ig,igp,ios) = chi0(ig,igp,ios) - faux/(Ep%omega(ios)+delta(ig,igp,ios)) &
                                                 + faux/(Ep%omega(ios)-delta(ig,igp,ios))
           enddo
         enddo
       endif
     endif
   enddo

 CASE (1)

   ABI_ALLOCATE(paux_sym,(Ep%npwe*(Ep%npwe+1)/2))
   ABI_ALLOCATE(Sm1_gmG0,(Ep%npwe))

   do isym=1,Ltg_q%nsym_sg
     do itim = 1,Ltg_q%timrev

       if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then

         phmGt => Gsph_epsG0%phmGt(1:Ep%npwe,isym)
         gmG0  => Ltg_q%igmG0(1:Ep%npwe,itim,isym)
         Sm1_gmG0(1:Ep%npwe)=Gsph_epsG0%rottbm1(gmG0(1:Ep%npwe),itim,isym)

         SELECT CASE (itim)
         CASE (1)
           do ig = 1, Ep%npwe
             do igp = ig, Ep%npwe
               igaux=igp*(igp-1)/2+ig
               if (Sm1_gmG0(ig)<=Sm1_gmG0(igp)) then
                 igaux2=Sm1_gmG0(igp)*(Sm1_gmG0(igp)-1)/2+Sm1_gmG0(ig)
                 faux=paux(igaux2)
               else
                 igaux2=Sm1_gmG0(ig)*(Sm1_gmG0(ig)-1)/2+Sm1_gmG0(igp)
                 faux=conjg(paux(igaux2))
               endif
               paux_sym(igaux)=faux*conjg(phmGt(igp))*phmGt(ig)
             enddo
           enddo
         CASE (2)
           do ig = 1, Ep%npwe
             do igp = ig, Ep%npwe
               igaux=igp*(igp-1)/2+ig
               if (Sm1_gmG0(ig)<=Sm1_gmG0(igp)) then
                 igaux2=Sm1_gmG0(igp)*(Sm1_gmG0(igp)-1)/2+Sm1_gmG0(ig)
                 faux=conjg(paux(igaux2))
               else
                 igaux2=Sm1_gmG0(ig)*(Sm1_gmG0(ig)-1)/2+Sm1_gmG0(igp)
                 faux=paux(igaux2)
               endif
               paux_sym(igaux)=faux*conjg(phmGt(igp))*phmGt(ig)
             enddo
           enddo
         CASE DEFAULT
           write(msg,'(a,i3)')'Wrong itim= ',itim
           MSG_BUG(msg)
         END SELECT

         do ios = 1, Ep%nomega
           if (ABS(REAL(Ep%omega(ios)))<0.00001) then
             if (abs(aimag(Ep%omega(ios)))>=0.001) then
               do ig = 1, Ep%npwe
                 do igp = ig, Ep%npwe
                   igaux=igp*(igp-1)/2+ig
                   faux=paux_sym(igaux)
                   delta_sym=delta(Sm1_gmG0(ig),Sm1_gmG0(igp),ios)
                   chi0(ig,igp,ios) = chi0(ig,igp,ios) - faux/(Ep%omega(ios)+delta_sym) &
                                                       + faux/(Ep%omega(ios)-delta_sym)
                 enddo
               enddo
             endif
           else
             if (abs(aimag(Ep%omega(ios)))>=0.001) then
               do ig = 1, Ep%npwe
                 do igp = 1, Ep%npwe
                   if (ig<=igp) then
                     igaux=igp*(igp-1)/2+ig
                     faux=paux_sym(igaux)
                   else
                     igaux=ig*(ig-1)/2+igp
                     faux=conjg(paux_sym(igaux))
                   endif
                   delta_sym=delta(Sm1_gmG0(ig),Sm1_gmG0(igp),ios)
                   chi0(ig,igp,ios) = chi0(ig,igp,ios) - faux/(Ep%omega(ios)+delta_sym) &
                                                       + faux/(Ep%omega(ios)-delta_sym)
                 enddo
               enddo
             endif
           endif
         enddo

       endif
     enddo
   enddo

   ABI_DEALLOCATE(paux_sym)
   ABI_DEALLOCATE(Sm1_gmG0)


 CASE DEFAULT
  write(msg,'(a,i3)')'Wrong symchi= ',Ep%symchi
  MSG_BUG(msg)
 END SELECT

end subroutine calc_chi0_delta0_bis
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/check_delta
!!
!! NAME
!! check_delta
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
!!      cchi0_eet
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine check_delta(Ep,qpgsq,delta,epsv,epslumo,ig,igp)

 use m_profiling

 use defs_basis
 use m_errors
 use m_gwdefs
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_delta'
!End of the abilint section

 implicit none

 type(Epsilonm1_parameters),intent(in) :: Ep

!Arguments ------------------------------------
!scalars
!arrays
 integer,intent(in) :: ig,igp
 real(dp),intent(in) :: qpgsq(Ep%npwe)
 real(dp),intent(in) :: epsv,epslumo
 complex(gwpc),intent(inout) :: delta

!Local variables-------------------------------
!scalars

 real(dp) :: test

!*************************************************************************

 if (ig==igp) then
   delta=real(delta)
   test=delta+epsv
   if (test<epslumo) then
     delta=half*(qpgsq(ig)+qpgsq(igp))
     test=delta+epsv
   endif
   if (test<epslumo) then
     delta=epslumo-epsv
   endif
 else
   test=delta+epsv
   if (test<epslumo) then
     delta=epslumo-epsv
   endif
 endif

end subroutine check_delta
!!***
