!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrtout
!! NAME
!!  wrtout
!!
!! FUNCTION
!!  Organizes the sequential or parallel version of the write intrinsic
!!  Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
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
!!  msg=(character(len=*)) message to be written
!!  unit=unit number for writing. The named constant dev_null defined in defs_basis can be used to avoid any printing.
!!  [mode_paral]= --optional argument--
!!   'COLL' if all procs are calling the routine with the same message to be written once only. Default.
!!   'PERS' if the procs are calling the routine with different messages each to be written,
!!          or if one proc is calling the routine
!!   "INIT" to change the rank of the master node that prints the message if "COLL" is used.
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      abi_etsf_electrons_put,abi_etsf_geo_put,abi_etsf_init,abinit,acfd_dyson
!!      acfd_intexact,afterscfloop,anaddb,append_cml2,append_xyz,asria_calc
!!      asria_corr,asrif9,asrprs,atm2fft,atomden,berryphase,berryphase_new
!!      besjm,bestwfs,bethe_salpeter,bfactor,bldgrp,bonds_lgth_angles,bound
!!      bound_new,calc_cs,calc_efg,calc_fc,calc_optical_mels
!!      calc_rpa_functional,calc_sig_ppm_eet,calc_sigc_me,calc_sigx_me
!!      calc_vhxc_me,calcdensph,canat9,ccfft,cchi0,cchi0q0,cchi0q0_intraband
!!      cgwf,cgwf3,check_completeness,chiscwrt,chkdilatmx,chkdpr,chkexi,chki8
!!      chkilwf,chkin9,chkinp,chkint_prt,chkneu,chknm8,chkorthsy,chkpawovlp
!!      chkph3,chkprimit,chkr8,chkrp9,chkvars,chneu9,clnup1,clnup2,cohsex_me
!!      constrf,contract_dp_ge_val,contract_int_ge_val,contract_int_le_val
!!      contract_int_list,cprj_utils_mpi,crho,cut3d,cvxclda,d3output
!!      datafordmft,ddkten,debug_tools,defs_scalapack,deloc2xcart,denfgr
!!      der_int,diel9,dielmt,dielmt2,dieltcel,distrb2,dmft_solve,dos_hdr_write
!!      driver,drivexc,dsksta,dtsetcopy,dyson,dyson_sc,echo_xc_name,eig1fixed
!!      eig2tot,eigen_meandege,elast9,electrooptic,eliashberg_1d,elph2_fanddw
!!      elphon,elpolariz,eltxccore,energy,entropyrec,ep_setupqpt,etot3
!!      evdw_wannier,ewald,ewald3,ewald4,ewald9,exc_build_block,exc_build_ham
!!      exc_den,exc_diago,exc_iterative_diago,exc_plot,exc_spectra,extraprho
!!      extrapwf,fconv,fermi,fermi_green,fermisolverec,fftprof,fftw,fftwfn
!!      fillcell,filterpot,find_getdtset,finddistrproc,findmin,findminscf
!!      first_rec,fixsym,forstr,forstrnps,fred2fdeloc,fsumrule,ftgkk,fxphas
!!      gath3,gaus_dos,gensymshub,gensymshub4,gensymspgr,get_all_gkq
!!      get_fs_bands,get_full_gsphere,get_g_tiny,get_veloc_tr,getcut
!!      getdim_nloc,getfreqsus,getgh1c,getghc,getkgrid,getlambda,getmpw,getnel
!!      getng,getshell,getspinrot,gran_potrec,green_kernel,gstate,gstateimg
!!      gtblk9,gtdyn9,gw_driver,gw_tools,gwcompleteness,hartre,hartre1
!!      hartrestr,haydock,haydock_psherm,hdr_check,hdr_io_etsf,hdr_io_netcdf
!!      hdr_vs_dtset,herald,hermit,hessinit,hubbard_one,importcml,importxyz
!!      impurity_solve,inarray,incomprs,ingeo,ingeobld,init8,init_occ_ent
!!      initberry,initberry3,initmpi_respfn,initmv,initrhoij,initro,initwf
!!      initylmg,inkpts,inpgkk,inprep8,inpspheads,inqpt,instr9,instrng,insy3
!!      int2char,int2char4,int_ang,intagm,integrate_gamma,integrate_gamma_alt
!!      integrate_gamma_tr,inupper,invars0,invars1,invars1m,invars2,invars2m
!!      invars9,inwffil,inwffil3,ioarr,ioddb8_in,ioddb8_out,iofn1,iofn2
!!      ioniondist,irred_perts,irrzg,isfile,jellium,klocal,kpgio,kpgstr
!!      kramerskronig,ks_ddiago,kss2wfk,kxc_alda,kxc_eok,ladielmt,lattice,lavnl
!!      ldau_self,leave_new,leave_test,linemin,listkk,lobpcgIIwf,lobpcgccIIwf
!!      lobpcgwf,local_ks_green,loop3dte,loper3,lwf,m_ab6_invars_f90,m_abilasi
!!      m_atom,m_bands_sym,m_bs_defs,m_bse_io,m_bz_mesh,m_cppopts_dumper
!!      m_crystal,m_crystal_io,m_dyson_solver,m_ebands,m_energy,m_errors
!!      m_fft_mesh,m_fft_prof,m_fftw3,m_gamma,m_geometry,m_gpu_detect,m_green
!!      m_gsphere,m_hamiltonian,m_header,m_hidecudarec,m_hu,m_initcuda,m_io_gkk
!!      m_io_kss,m_io_screening,m_iterators,m_libxc_functionals,m_matlu
!!      m_matrix,m_melemts,m_numeric_tools,m_oper,m_optim_dumper,m_paw_dmft
!!      m_paw_pwij,m_paw_slater,m_paw_toolbox,m_phdos,m_pimd,m_ppmodel
!!      m_pretty_rec,m_ptgroups,m_qparticles,m_radmesh,m_rec,m_results_respfn
!!      m_screen,m_screening,m_self,m_shexc,m_shirley,m_sigma_results,m_timer
!!      m_vcoul,m_wfs,m_xredistribute,mag_out,mat_mlms2jmj,mat_slm2ylm,matcginv
!!      matcginv_dpc,mati3inv,matrginv,matrixelmt_g,mblktyp1,mblktyp5,mean_fftr
!!      meanvalue_g,memana,memkss,memorf,memory,metcon,metric,metstr,mka2f
!!      mka2fQgrid,mka2f_tr,mkcor3,mkcore,mkdenpos,mkeuler,mkffnl,mkfilename
!!      mkfskgrid,mkifc9,mkkpg,mklocl,mklocl_realspace,mklocl_recipspace
!!      mklocl_wavelets,mknesting,mknormpath,mkph_linwid,mkphbs,mkqptequiv
!!      mkrho,mkrho3,mkvxc3,mkvxcstr3,mlwfovlp,mlwfovlp_proj,mlwfovlp_projpaw
!!      mlwfovlp_pw,mlwfovlp_qp,mlwfovlp_radial,mlwfovlp_seedname
!!      mlwfovlp_setup,mlwfovlp_ylmfac,mlwfovlp_ylmfar,moddiel,mover
!!      mpi_enreg_tools,mrgddb,mrggkk,mrgscr,multipoles_fftr,mv_3dte
!!      my_calc_wfwfg,nanal9,newfermie1,newkpt,newocc,newrho,newsp,newton
!!      newvtr,newvtr3,nlenergyrec,nmsq_gam,nmsq_pure_gkk,nonlinear,nonlop
!!      nonlop_pl,normev,normsq_gkq,nres2vres,nselt3,nstdy3,nstpaw3,occeig
!!      occred,odamix,old_setmesh,opernl2,opernl3,opernl4b,optics_paw
!!      optics_paw_core,orthonormalize,out1dm,outelph,outg2f,outgkk,outkss
!!      outphdos,outqmc,outscfcv,outvars,outwant,outwf,outxfhist,overlap_g
!!      pareigocc,parsefile,partial_dos_fractions_paw,paw_mknewh0,paw_qpscgw
!!      pawdenpot,pawdensities,pawdij,pawfrnhat_recipspace,pawinit,pawlsylm
!!      pawmkaewf,pawmkrhoij,pawprt,pawpupot,pawpuxinit,pawtwdij,pawuenergy
!!      pawuj_det,pawuj_red,pawuj_utils,pawxcpositron,pawxcsph,pawxcsphpositron
!!      pawxenergy,pawxpot,pclock,ph1d3d,phfrq3,piezo9,pimd_nosehoover_nvt
!!      pl_deriv,plm_coeff,plm_d2theta,plm_dphi,plm_dtheta,polcart,poslifetime
!!      prcref,prcref_PMA,prctfvw1,prctfvw2,precon,precon2,pred_delocint
!!      pred_isokinetic,pred_isothermal,pred_langevin,pred_nose,pred_verlet
!!      predictimg,prep_bandfft_tabs,prep_kpgio,print_ierr,print_ij,print_psps
!!      prmat,projbd,prt_cml2,prteigrs,prtene,prtene3,prtfatbands,prtimg,prtocc
!!      prtph3,prtrhomxmn,prtspgroup,prttagm,prtvsound,prtxf,prtxfase,prtxvf
!!      psddb8,psichi_renormalization,psolver_hartree,psolver_kernel
!!      psolver_rhohxc,psp10in,psp10nl,psp1cc,psp1in,psp1nl,psp2in,psp2lo
!!      psp3in,psp3nl,psp4cc,psp5in,psp5nl,psp6in,psp7in,psp7nl,psp8cc,psp8in
!!      psp9in,pspatm,pspini,pspnl_hgh_rec,pspnl_operat_rec,psxml2ab,ramansus
!!      randac,randomcellpos,rdddb9,rdm,rdnpw,read_blok8,read_gkk,readeig
!!      relaxpol,remove_inversion,respfn,rhofermi3,rhohxc,rhohxcpositron
!!      rhotov3,rotmat,rrho,rsiaf9,rwwan,scalewf_nonlop,scfcge,scfcv,scfcv3
!!      scfeig,scfopt,scphon,scphon_build_qsym_map,scphon_dynmat_to_freq2
!!      scphon_free_energy,scphon_ft_fcart,scphon_supercell_vectors_init,scprqt
!!      screening,setmqgrid,setnoccmmp,setrhoijpbe0,setshells,setsymrhoij
!!      setup1,setup2,setup_G_rotation_old,setup_bse,setup_positron,setup_qmesh
!!      setup_screening,setup_sigma,setvtr,sg_ctrig,sg_fft,sg_fftpx
!!      sg_fftrisc_2,sg_fftx,sg_ffty,shellstruct,sigma,smallprim,smatrix
!!      smatrix_paw,smatrix_pawinit,smeared_delta,smpbz,spectral
!!      spectral_function,sphereboundary,sphericaldens,spin_current,status
!!      stress,subdiago,sumrule,suscep,suscep_dyn,suscep_kxc_dyn,suscep_stat
!!      suskmm,suskmm_dyn,suskmm_kxc_dyn,sym_cprj_kn,sym_gkk,symanal,symatm
!!      symaxes,symbrav,symcharac,symdet,symdij,symdm9,symfind,symg,symkchk
!!      symkpt,symlatt,symmetrize_afm_chi0,symmultsg,symph3,symplanes
!!      symptgroup,symq3,symrelrot,symrhg,symrhoij,symspgr,tddft,testkgrid
!!      tetrahedron,thm9,thmeig,timab,timana,time_accu,timein,uderiv,ujdet
!!      vlocalstr,vso_realspace_nonlop,vtorho,vtorho3,vtorhorec,vtorhotf,vtowfk
!!      vtowfk3,wfconv,wfd_mkrho,wfd_pawrhoij,wffile,wffopen,wfkfermi3,wfsinp
!!      wght9,wrt_moldyn_netcdf,wrtloctens,wvl_memory,wvl_mkrho,wvl_newvtr
!!      wvl_nl_gradient,wvl_projectors_free,wvl_projectors_set,wvl_rwwf
!!      wvl_setboxgeometry,wvl_setngfft,wvl_tail_corrections,wvl_vtorho
!!      wvl_wfs_free,wvl_wfs_set,wvl_wfsinp_disk,wvl_wfsinp_reformat
!!      wvl_wfsinp_scratch,xc_kernel,xc_kernel_ADA,xcacfd,xcden,xchcth,xchelu
!!      xcpbe,xcpot,xcpzca,xcspol,xctetr,xcwign,xcxalp,xfpack_f2vout
!!      xfpack_vin2x,xfpack_x2vin,xredxcart,zorthonormalize,zprecon3
!!
!! CHILDREN
!!      wrtout_myproc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wrtout(unit,msg,mode_paral)

 use m_profiling

 use defs_basis
 use m_xmpi,    only : xmpi_world, xcomm_rank, xcomm_size

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrtout'
 use interfaces_14_hidewrite, except_this_one => wrtout
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unit
 character(len=*),optional,intent(in) :: mode_paral
 character(len=*),intent(in) :: msg

!Local variables-------------------------------
 integer :: comm
 integer,save :: master=0
 integer :: me,nproc,rtnpos
 character(len=7) :: tag
 character(len=len(msg)) :: my_msg
 character(len=len(msg)+50) :: string
 character(len=500) :: my_mode_paral

!******************************************************************

!Be careful with the coding  of the parallel case ...

!MG: Be careful**2, me and nproc are defined in MPI_COMM_WORLD.
!One should pass the MPI communicator

 if ((unit == std_out).and.(.not.do_write_log)) RETURN
 if (unit == dev_null) RETURN

 my_mode_paral = "COLL"; if (PRESENT(mode_paral)) my_mode_paral = mode_paral

!Communicator is xmpi_world by default, except for the parallelization over images
 if (abinit_comm_output/=-1) then
   comm=abinit_comm_output
 else
   comm=xmpi_world
 end if

!Determine who I am in COMM_WORLD
 nproc = xcomm_size(comm)
 me    = xcomm_rank(comm)

!msg is not changed therefore we can pass literal strings as well.
 my_msg = msg

 if( (my_mode_paral=='COLL') .or. (nproc==1) ) then

   if (me==master) then
     call wrtout_myproc(unit, my_msg)
   end if

 else if (my_mode_paral=='PERS') then

   if(me<10) then
     write(tag,'("-P-000",i1)') me
   elseif(me<100) then
     write(tag,'("-P-00",i2)') me
   elseif(me<1000) then
     write(tag,'("-P-0",i3)') me
   elseif(me<10000) then
     write(tag,'("-P-",i4)') me
   else
     tag=' ######'
   end if

   rtnpos=index(my_msg,ch10)
   do while(rtnpos/=0)
     write(string,'(3a)') tag, ' ', my_msg(1:rtnpos-1)
     write(unit,'(A)') trim(string)
     my_msg=my_msg(rtnpos+1:len(my_msg))
     rtnpos=index(my_msg,ch10)
   end do
   write(string, "(3a)") tag, ' ', my_msg
   write(unit,'(A)') trim(string)

 else if (my_mode_paral=='INIT') then

   master=unit

 else
   write(string,'(7a)')ch10,&
&   '  wrtout: ERROR -',ch10,&
&   '  Unknown write mode: ',my_mode_paral,ch10,&
&   '  Continuing anyway ...'
   write(unit, '(A)' ) trim(string)
 end if

end subroutine wrtout
!!***
