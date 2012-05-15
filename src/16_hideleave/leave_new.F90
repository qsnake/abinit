!{\src2tex{textfont=tt}}
!!****f* ABINIT/leave_new
!! NAME
!!  leave_new
!!
!! FUNCTION
!!  Routine for clean exit of f90 code, taking into account possible
!!  parallelization.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, NCJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mode_paral=
!!   'COLL' if all procs are calling the routine with the same message to be
!!     written once only or
!!   'PERS' if the procs are calling the routine with different mesgs
!!     each to be written, or if one proc is calling the routine
!!  print_config=(optional, default=true)
!!       if true print out several informations before leaving
!!
!! OUTPUT
!!  (only writing, then stop)
!!
!! NOTES
!!  By default, it uses "call exit(1)", that is not completely portable.
!!
!! PARENTS
!!      abi_etsf_init,abinit,acfd_dyson,acfd_intexact,afterscfloop,appdig
!!      append_cml2,append_xyz,asrif9,asrprs,berryphase,berryphase_new,besjm
!!      bestwfs,bfactor,bldgrp,bonds_lgth_angles,bound,bound_new,calc_cs
!!      calc_efg,calc_fc,canat9,ccfft,cgwf,cgwf3,chkdilatmx,chkdpr,chkexi,chki8
!!      chkilwf,chkin9,chkinp,chkint_prt,chkneu,chknm8,chkorthsy,chkpawovlp
!!      chkprimit,chkr8,chkrp9,chkvars,constrf,contract_dp_ge_val
!!      contract_int_ge_val,contract_int_le_val,contract_int_list
!!      cprj_utils_mpi,cut3d,cvxclda,datafordmft,ddkten,deloc2xcart,der_int
!!      diel9,dielmt2,dieltcel,distrb2,dmft_solve,drivexc,dyson,dyson_sc
!!      eig2tot,eigen_meandege,elph2_fanddw,elpolariz,eltxccore,energy,etot3
!!      ewald3,ewald9,extraprho,extrapwf,fappnd,fermisolverec,fftw,fillcell
!!      filterpot,find_getdtset,finddistrproc,findmin,fixsym,forstr,ftgkk
!!      fxphas,gath3,gensymshub,gensymshub4,gensymspgr,get_g_tiny,get_veloc_tr
!!      getcut,getfreqsus,getgh1c,getghc,getlambda,getshell,gstate,gtblk9
!!      hartre,hartre1,hartrestr,hdr_check,hdr_io_etsf,hdr_io_netcdf,hermit
!!      hessinit,hubbard_one,importcml,importxyz,inarray,ingeo,ingeobld,init8
!!      init_occ_ent,initberry,initberry3,initmpi_respfn,initmv,initrhoij
!!      initro,initylmg,inkpts,inpgkk,inprep8,inpspheads,inqpt,inread,instrng
!!      int2char,int2char4,int_ang,intagm,integrate_gamma_alt
!!      integrate_gamma_tr,inupper,invars0,invars1,invars1m,invars2,invars9
!!      inwffil,inwffil3,ioarr,ioddb8_in,iofn1,iofn2,irrzg,isfile,jellium
!!      klocal,kpgstr,kxc_alda,kxc_eok,ladielmt,lavnl,linemin,listkk,lobpcgIIwf
!!      lobpcgccIIwf,local_ks_green,lwf,m_ab6_invars_f90,m_energy,m_errors
!!      m_green,m_header,m_libxc_functionals,m_matlu,m_matrix,m_oper,m_paw_dmft
!!      m_rec,m_results_respfn,m_self,m_timer,m_xredistribute,mat_mlms2jmj
!!      mat_slm2ylm,matcginv,matcginv_dpc,mati3inv,matrginv,matrixelmt_g
!!      mblktyp1,mblktyp5,mean_fftr,meanvalue_g,memana,metcon,metstr,mka2f_tr
!!      mkcor3,mkcore,mkdenpos,mkeuler,mkffnl,mkfilename,mkkpg,mklocl
!!      mklocl_realspace,mklocl_recipspace,mklocl_wavelets,mknormpath,mkvxc3
!!      mkvxcstr3,mlwfovlp,mlwfovlp_proj,mlwfovlp_projpaw,mlwfovlp_pw
!!      mlwfovlp_radial,mlwfovlp_setup,mlwfovlp_ylmfac,mlwfovlp_ylmfar,moddiel
!!      mpi_enreg_tools,mrgddb,mrggkk,nanal9,newkpt,newrho,newsp,newvtr,newvtr3
!!      nmsq_gam,nmsq_pure_gkk,nonlinear,nonlop,nonlop_pl,normev,nres2vres
!!      occeig,odamix,old_setmesh,opernl2,opernl3,opernl4b,optics_paw
!!      optics_paw_core,out1dm,outg2f,outphdos,outqmc,outscfcv,outwant,outwf
!!      outxfhist,overlap_g,parsefile,partial_dos_fractions_paw,pawdij
!!      pawfrnhat_recipspace,pawinit,pawpuxinit,pawtwdij,pawxcsph
!!      pawxcsphpositron,ph1d3d,pl_deriv,plm_coeff,plm_d2theta,plm_dphi
!!      plm_dtheta,prcref,prcref_PMA,prctfvw1,prctfvw2,pred_isothermal
!!      prep_bandfft_tabs,prep_kpgio,print_ierr,print_psps,prteigrs,prtfatbands
!!      prtocc,prtph3,prtrhomxmn,prtspgroup,prttagm,psddb8
!!      psichi_renormalization,psolver_hartree,psolver_kernel,psolver_rhohxc
!!      psp10nl,psp11nl,psp1cc,psp1in,psp1nl,psp2in,psp3nl,psp4cc,psp5in,psp5nl
!!      psp6in,psp7in,psp7nl,psp8cc,psp8in,pspatm,psxml2ab,ramansus
!!      randomcellpos,rdm,rdnpw,read_blok8,readeig,relaxpol,respfn,rhofermi3
!!      rhohxc,rhohxcpositron,rhotov3,rotmat,rrho,rsiaf9,rwwan,scalewf_nonlop
!!      scfcv,scfcv3,scphon_ft_fcart,scprqt,setnoccmmp,setrhoijpbe0,setup1
!!      setup_G_rotation_old,setup_positron,sg_ctrig,sg_fft,sg_fftpx
!!      sg_fftrisc_2,sg_fftx,sg_ffty,smallprim,smatrix,smatrix_paw
!!      smatrix_pawinit,smeared_delta,spectral_function,sphereboundary
!!      sphericaldens,spin_current,status,subdiago,suscep_stat,suskmm
!!      suskmm_dyn,suskmm_kxc_dyn,sym_cprj_kn,sym_gkk,symdet,symdm9,symfind
!!      symg,symkchk,symlatt,symptgroup,symrelrot,symrhg,symspgr,tddft
!!      testkgrid,tetrahedron,thm9,thmeig,timab,time_accu,timein,uderiv,ujdet
!!      uniformrandom,vlocalstr,vso_realspace_nonlop,vtorho,vtorho3,vtorhorec
!!      vtowfk,vtowfk3,wfconv,wffile,wffopen,wfkfermi3,wfsinp,wght9,wvl_memory
!!      wvl_mkrho,wvl_newvtr,wvl_nl_gradient,wvl_projectors_free
!!      wvl_projectors_set,wvl_rwwf,wvl_setboxgeometry,wvl_setngfft
!!      wvl_tail_corrections,wvl_vtorho,wvl_wfs_free,wvl_wfs_set
!!      wvl_wfsinp_disk,wvl_wfsinp_reformat,wvl_wfsinp_scratch,xcacfd,xcden
!!      xchcth,xchelu,xcpbe,xcpot,xcpzca,xcspol,xctetr,xcwign,xcxalp
!!      xfpack_f2vout,xfpack_vin2x,xfpack_x2vin,xredxcart
!!
!! CHILDREN
!!      dump_config,leave_myproc,print_kinds,wrtout,xmpi_abort,xmpi_end
!!      xmpi_show_info
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine leave_new(mode_paral,print_config)

 use m_profiling

 use defs_basis
 use m_build_info
 use m_xmpi
#if defined FC_NAG
 use f90_unix
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'leave_new'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave, except_this_one => leave_new
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=4),intent(in) :: mode_paral
 logical,intent(in),optional :: print_config

!Local variables-------------------------------
 integer :: nproc
 logical :: print_config_
 character(len=500) :: msg

! **********************************************************************

 write(msg,'(2a)') ch10,' leave_new : decision taken to exit ...'
 call wrtout(std_out,msg,'PERS')

 print_config_=.true.;if (present(print_config)) print_config_=print_config

!Determine nproc in COMM_WORLD
 nproc = xcomm_size(xmpi_world)

!Dump configuration before exiting
 if (print_config_) then
   call print_kinds()
   call xmpi_show_info()
   call dump_config()
 end if

 if( (mode_paral=='COLL') .or. (nproc==1) ) then

   if (abinit_comm_leave==xmpi_world) then
     call xmpi_end()
     call leave_myproc()
   else
     call leave_myproc(option=1)
     call xmpi_abort()
   end if

 else if(mode_paral=='PERS') then  ! Caveat: Do not use MPI collective calls!

   call wrtout(std_out,' leave_new : calling XMPI_ABORT...','PERS')
   call leave_myproc(option=1)
   call xmpi_abort()

 end if

end subroutine leave_new
!!***
