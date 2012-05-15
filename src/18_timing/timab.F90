!{\src2tex{textfont=tt}}
!!****f* ABINIT/timab
!! NAME
!!  timab
!!
!! FUNCTION
!!  Timing subroutine.  Calls machine-dependent "timein" which
!!  returns elapsed cpu and wall clock times in sec.
!!
!!  Depending on value of "option" routine will:
!!  (0) zero all accumulators
!!  (1) start with new incremental time slice for accumulator n
!!        using explicit call to timein (or PAPI)
!!  (2) stop time slice; add time to accumulator n also increase by one the counter for this accumulator
!!  (3) start with new incremental time slice for accumulator n
!!        using stored values for cpu, wall, and PAPI infos ( ! do not use for stop )
!!  (4) report time slice for accumlator n (not full time accumlated)
!!  (5) option to suppress timing (nn should be 0) or reenable it (nn /=0)
!!
!!  If, on first entry, subroutine is not being initialized, it
!!  will automatically initialize as well as rezero accumulator n.
!!  However, initialization SHOULD be done explicitly by the user
!!  so that it can be done near the top of his/her main routine.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nn=index of accumulator (distinguish what is being timed); NOT used if option=0
!!  option=see comment above
!!
!! OUTPUT
!!  on option=4 :
!!    tottim(2,nn)=accumulated time for accumulator nn; otherwise
!!     tottim is a dummy variable.
!!    option gives the number of times that the
!!     accumulator has been incremented
!!
!! PARENTS
!!      abinit,acfd_dyson,acfd_intexact,afterscfloop,atm2fft,back_wf,bestwfs
!!      bethe_salpeter,calc_sigc_me,calc_sigx_me,calcdensph,cchi0,cgwf,cgwf3
!!      check_completeness,cohsex_me,corrmetalwf1,density_rec,dielmt,dielmt2
!!      dieltcel,dotprod_g,dotprod_v,dotprod_vn,dotprodm_v,dotprodm_vn,driver
!!      dyfnl3,dyfro3,dyxc13,eig2tot,eltfrhar3,eltfrkin3,eltfrloc3,eltfrnl3
!!      eltfrxc3,energy,entropyrec,etotfor,exc_build_block,exc_build_ham
!!      fermisolverec,fftw3_fourdp,filnam_comm,filterpot,first_rec,forces
!!      forstr,forstrnps,forw_wf,fourdp,fourwf,fxphas,getgh1c,getghc,getgsc
!!      getngrec,gran_potrec,green_kernel,gstate,gstateimg,hartre,hartre1
!!      haydock,initwf,initylmg,inkpts,invars2,inwffil,inwffil3,kpgio,kpgsph
!!      ladielmt,lavnl,leave_test,lobpcgccIIwf,lobpcgwf,loop3dte,loper3
!!      m_ab6_invars_f90,m_ab6_mixing,m_dyson_solver,m_eigen,m_lobpcg,m_wfutils
!!      matrixelmt_g,mean_fftr,meanvalue_g,mkcore,mkffnl,mklocl_realspace
!!      mklocl_recipspace,mkresi,mkrho,mkrho3,mkvxc3,mkvxcstr3,newkpt,newocc
!!      newrho,newvtr,newvtr3,nhatgrid,nlenergyrec,nonlinear,nonlop,nstdy3
!!      nstpaw3,nstwf3,odamix,opernla_ylm,optics_paw,optics_paw_core
!!      optics_vloc,outkss,outscfcv,outwf,pareigocc,partial_dos_fractions_paw
!!      pawdenpot,pawdij,pawinit,pawmknhat,pawmkrho,pawmkrhoij,pawnstd2e
!!      pawpolev,pawxc,pawxc3,pawxc3_gga,pawxcm,pawxcm3,prctfvw1,prctfvw2
!!      precon,precon2,prep_bandfft_tabs,prep_fourwf,prep_getghc,prep_nonlop
!!      projbd,pspheads_comm,pspini,pw_orthon,recursion,recursion_nl,redgr
!!      respfn,rhofermi3,rhohxc,rhoij_utils,rhotov,rhotov3,rwwf,scfcv,scfcv3
!!      screening,setsym,setvtr,sigma,sqnorm_g,sqnorm_v,sqnormm_v,status,stress
!!      strhar,subdiago,suscep,suscep_dyn,suscep_kxc_dyn,suscep_stat,susk
!!      susk_dyn,susk_dyn_pgg,susk_kxc_dyn,suskmm,suskmm_dyn,suskmm_kxc_dyn
!!      symrhg,symsgcube,tddft,timana,vn_nl_rec,vtorho,vtorho3,vtorhorec
!!      vtorhotf,vtowfk,vtowfk3,wfconv,wfkfermi3,wfsinp,xcden,xcpot,zprecon3
!!
!! CHILDREN
!!      leave_new,papif_flops,papif_perror,timein,wrtout
!!
!! SOURCE
!!


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine timab(nn,option,tottim)

 use m_profiling

 use defs_basis
 use defs_time

#ifdef HAVE_FC_ISO_C_BINDING
 use iso_c_binding
#else
 use m_iso_c_binding
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'timab'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing, except_this_one => timab
!End of the abilint section

 implicit none

#if defined HAVE_TIMER_PAPI
#include "f90papi.h"
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn,option
!arrays
 real(dp),intent(out) :: tottim(2)

!Local variables-------------------------------
!scalars
! real(dp) :: cpu,wall
 character(len=500) :: message
! integer(C_LONG_LONG) :: flops1
! real(C_FLOAT) :: real_time, proc_time
#if defined HAVE_TIMER_PAPI
 integer(C_INT) :: retval 
 real(C_FLOAT) :: mflops1
 character(len=PAPI_MAX_STR_LEN) :: papi_errstr
#endif
! *************************************************************************

!DEBUG
!write(std_out,*)' timab : enter with  nn, option, timopt',nn,option,timopt
!if(entry==5)stop
!ENDDEBUG

 if (option==5) timopt=nn

!If timopt was set to zero by a call with option=5, suppress
!all action of this routine (might as well return at this point !)
 if(timopt/=0 .and. option/=5)then
!  
!  Check that nn lies in sensible bounds
   if (nn<1.or.nn>mtim) then
     write(message, '(a,a,a,a,i6,a,i8,a)' ) ch10,&
&     ' timab: BUG -',ch10,&
&     '  dim mtim=',mtim,' but input nn=',nn,'.'
     call wrtout(std_out,message,'PERS')
     call leave_new('PERS')
   end if
#ifdef HAVE_TIMER_PAPI
!  ie for all active options for time  
!  if papi analysis has been selected or for initializing papi
!  even if papi has not been choose 
   if (option/=3.and.(papiopt==1.or.initpapiopt==1)) then 
     initpapiopt=0
     call PAPIf_flops(real_time, proc_time, flops1, mflops1, retval)
     if (retval.NE.PAPI_OK) then
       write(std_out,*) 'Problem to initialize papi high level inteface'
       call papif_perror(retval,papi_errstr,retval)
       write(std_out,*) 'Error code', papi_errstr
     end if ! DEBUG
!    write(std_out,*) 'flops  = ', flops1, 'mflops= ',  mflops1
     if (flops1 >= HUGE(flops1)) then  
       write(std_out,*) 'Mflops analysis : Number of floating point instruction Overflow '
       flops(:)=-1            
     end if
   end if
#endif
   
   select case (option)
     case (0)  ! Zero out all accumulators of time and init timers
       acctim(:,:)=0.0d0
       tzero(:,:)=0.0d0
       ncount(:)=0
       flops(:)=0
       papi_acctim(:,:)=0. 
       papi_accflops(:)=0. 
       papi_tzero(:,:)=0. 

     case (1)  ! Initialize timab for nn
       call timein(cpu,wall)
       tzero(1,nn)=cpu
       tzero(2,nn)=wall
#ifdef HAVE_TIMER_PAPI
       flops(nn)=flops1       ! Initialize megaflops for nn
       papi_tzero(1,nn) = proc_time
       papi_tzero(2,nn) = real_time
#endif

     case (2)  ! Accumulate time for nn (also keep the values of cpu, wall, proc_time, real_time, flops1)
       call timein(cpu,wall)
       acctim(1,nn)=acctim(1,nn)+cpu -tzero(1,nn)
       acctim(2,nn)=acctim(2,nn)+wall-tzero(2,nn)
       ncount(nn)=ncount(nn)+1
#ifdef HAVE_TIMER_PAPI
!      accumulate time and flops for nn Difference entre 2 calls a Papif_flops 
       papi_acctim(1,nn)=papi_acctim(1,nn)+ proc_time - papi_tzero(1,nn)
       papi_acctim(2,nn)=papi_acctim(2,nn)+ real_time - papi_tzero(2,nn)
       papi_accflops(nn)=papi_accflops(nn)+ flops1- flops(nn) 
#endif

     case (3) ! Use previously obtained values to initialize timab for nn
       tzero(1,nn)=cpu
       tzero(2,nn)=wall
#ifdef HAVE_TIMER_PAPI
       flops(nn)=flops1
       papi_tzero(1,nn) = proc_time
       papi_tzero(2,nn) = real_time
#endif

     case (4) ! Return elapsed time for nn (do not accumulate)
       call timein(cpu,wall)
       tottim(1)=cpu-tzero(1,nn)
       tottim(2)=wall-tzero(2,nn)
#ifdef HAVE_TIMER_PAPI
!      return ellapsed floating point operationfor nn (do not accumulate)
       papi_tottim(1,nn)= proc_time - papi_tzero(1,nn)
       papi_tottim(2,nn)= real_time - papi_tzero(2,nn)
       papi_totflops(nn)= flops1- flops(nn) 
#endif

       case default
       write(message, '(a,a,a,a,i10,a)' ) ch10,&
&       ' timab: BUG -',ch10,&
&       '  Input option not valid, =',option,'.'
       call wrtout(std_out,message,'PERS')
       call leave_new('PERS')
   end select 
 end if

!DEBUG
!write(std_out,*)' timab : exit '
!ENDDEBUG

end subroutine timab
!!***
