!{\src2tex{textfont=tt}}
!!****f* ABINIT/dtsetcopy
!! NAME
!! dtsetcopy
!!
!! FUNCTION
!! Copy all values of dataset dtin to dataset dtout. Pointers of dtout are
!! allocated if required. Use dtsetFree() to free a dataset after use.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MF, GZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtin <type(dataset_type)>=all input variables in this dataset
!!
!! OUTPUT
!!  dtout <type(dataset_type)>
!!
!! PARENTS
!!      afterscfloop,chkinp,cvxclda,driver,kxc_alda,m_io_kss,xc_kernel
!!      xc_kernel_ADA
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dtsetCopy(dtout, dtin)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_copy
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtsetCopy'
 use interfaces_53_abiutil, except_this_one => dtsetCopy
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtin
 type(dataset_type),intent(out) :: dtout

!Local variables-------------------------------
!scalars
 integer :: chosen_size1,chosen_size2,chosen_size3
 character(len=12) :: name

! *************************************************************************

 DBG_ENTER("COLL")

!BEGIN VARIABLES FOR @Bethe-Salpeter
 dtout%bs_algorithm     = dtin%bs_algorithm
 dtout%bs_haydock_niter = dtin%bs_haydock_niter
 dtout%bs_exchange_term = dtin%bs_exchange_term
 dtout%bs_coulomb_term  = dtin%bs_coulomb_term
 dtout%bs_calctype      = dtin%bs_calctype
 dtout%bs_coupling      = dtin%bs_coupling

 dtout%bs_haydock_tol   = dtin%bs_haydock_tol
 dtout%bs_hayd_term     = dtin%bs_hayd_term

 dtout%bs_loband        = dtin%bs_loband

 dtout%bs_eh_cutoff(:) = dtin%bs_eh_cutoff(:)
 dtout%bs_freq_mesh(:) = dtin%bs_freq_mesh(:)
!END VARIABLES FOR @Bethe-Salpeter.

!Copy integers from dtin to dtout
 dtout%accesswff          = dtin%accesswff
 dtout%awtr               = dtin%awtr
 dtout%bandpp             = dtin%bandpp
 dtout%bdeigrf            = dtin%bdeigrf
 dtout%berryopt           = dtin%berryopt
 dtout%berrystep          = dtin%berrystep
 dtout%brvltt             = dtin%brvltt
 dtout%bs_nstates         = dtin%bs_nstates
 dtout%cd_custom_imfrqs   = dtin%cd_custom_imfrqs
 dtout%cd_use_tangrid     = dtin%cd_use_tangrid
 dtout%cd_full_grid       = dtin%cd_full_grid
 dtout%chkexit            = dtin%chkexit
 dtout%chkgwcomp          = dtin%chkgwcomp
 dtout%chkprim            = dtin%chkprim
 dtout%chksymbreak        = dtin%chksymbreak
 dtout%delayperm          = dtin%delayperm
 dtout%diismemory         = dtin%diismemory
 dtout%dmatpuopt          = dtin%dmatpuopt
 dtout%dmatudiag          = dtin%dmatudiag
 dtout%dmft_dc            = dtin%dmft_dc
 dtout%dmft_iter          = dtin%dmft_iter
 dtout%dmft_mxsf          = dtin%dmft_mxsf
 dtout%dmft_nwlo          = dtin%dmft_nwlo
 dtout%dmft_nwli          = dtin%dmft_nwli
 dtout%dmft_rslf          = dtin%dmft_rslf
 dtout%dmft_solv          = dtin%dmft_solv
 dtout%dmftbandi          = dtin%dmftbandi
 dtout%dmftbandf          = dtin%dmftbandf
 dtout%dmftcheck          = dtin%dmftcheck
 dtout%enunit             = dtin%enunit
 dtout%exchn2n3d          = dtin%exchn2n3d
 dtout%pawfatbnd          = dtin%pawfatbnd
 dtout%fermie_nest        = dtin%fermie_nest
 dtout%fftgw              = dtin%fftgw
 dtout%fft_opt_lob        = dtin%fft_opt_lob
 dtout%freqremin          = dtin%freqremin
 dtout%freqremax          = dtin%freqremax
 dtout%freqspmin          = dtin%freqspmin
 dtout%freqspmax          = dtin%freqspmax
 dtout%frzfermi           = dtin%frzfermi
 dtout%getbseig           = dtin%getbseig
 dtout%getbsreso          = dtin%getbsreso
 dtout%getbscoup          = dtin%getbscoup
 dtout%getcell            = dtin%getcell
 dtout%getddk             = dtin%getddk
 dtout%getden             = dtin%getden
 dtout%getgam_eig2nkq     = dtin%getgam_eig2nkq
 dtout%gethaydock         = dtin%gethaydock
 dtout%getkss             = dtin%getkss
 dtout%getocc             = dtin%getocc
 dtout%getpawden          = dtin%getpawden
 dtout%getqps             = dtin%getqps
 dtout%getscr             = dtin%getscr
 dtout%getsuscep          = dtin%getsuscep
 dtout%getvel             = dtin%getvel
 dtout%getwfk             = dtin%getwfk
 dtout%getwfq             = dtin%getwfq
 dtout%getxcart           = dtin%getxcart
 dtout%getxred            = dtin%getxred
 dtout%get1den            = dtin%get1den
 dtout%get1wf             = dtin%get1wf
 dtout%goprecon           = dtin%goprecon
 dtout%gwcalctyp          = dtin%gwcalctyp
 dtout%gwcomp             = dtin%gwcomp
 dtout%gwencomp           = dtin%gwencomp
 dtout%gwmem              = dtin%gwmem
 dtout%gwpara             = dtin%gwpara
 dtout%gwgamma            = dtin%gwgamma
 dtout%gwrpacorr          = dtin%gwrpacorr
 dtout%gw_custom_freqsp   = dtin%gw_custom_freqsp
 dtout%gw_nqlwl           = dtin%gw_nqlwl
 dtout%gw_eet_nband       = dtin%gw_eet_nband
 dtout%gw_nstep           = dtin%gw_nstep
 dtout%gw_eet             = dtin%gw_eet
 dtout%gw_eet_inclvkb     = dtin%gw_eet_inclvkb
 dtout%gw_sctype          = dtin%gw_sctype
 dtout%gw_sigxcore        = dtin%gw_sigxcore
 dtout%gw_toldfeig        = dtin%gw_toldfeig
 dtout%gw_use_pole_scr    = dtin%gw_use_pole_scr
 dtout%gw_reconst_scr     = dtin%gw_reconst_scr
 dtout%gw_npoles          = dtin%gw_npoles
 dtout%iboxcut            = dtin%iboxcut
 dtout%icoulomb           = dtin%icoulomb
 dtout%icutcoul           = dtin%icutcoul
 dtout%idyson             = dtin%idyson
 dtout%ieig2rf            = dtin%ieig2rf
 dtout%iextrapwf          = dtin%iextrapwf
 dtout%ikhxc              = dtin%ikhxc
 dtout%imgmov             = dtin%imgmov
 dtout%inclvkb            = dtin%inclvkb
 dtout%intexact           = dtin%intexact
 dtout%intxc              = dtin%intxc
 dtout%ionmov             = dtin%ionmov
 dtout%iprcch             = dtin%iprcch
 dtout%iprcel             = dtin%iprcel
 dtout%iprctfvw           = dtin%iprctfvw
 dtout%iprcfc             = dtin%iprcfc
 dtout%irandom            = dtin%irandom
 dtout%irdbseig           = dtin%irdbseig
 dtout%irdbsreso          = dtin%irdbsreso
 dtout%irdbscoup          = dtin%irdbscoup
 dtout%irdddk             = dtin%irdddk
 dtout%irdden             = dtin%irdden
 dtout%irdhaydock         = dtin%irdhaydock
 dtout%irdkss             = dtin%irdkss
 dtout%irdpawden          = dtin%irdpawden
 dtout%irdqps             = dtin%irdqps
 dtout%irdscr             = dtin%irdscr
 dtout%irdsuscep          = dtin%irdsuscep
 dtout%irdwfk             = dtin%irdwfk
 dtout%irdwfq             = dtin%irdwfq
 dtout%ird1den            = dtin%ird1den
 dtout%ird1wf             = dtin%ird1wf
 dtout%iscf               = dtin%iscf
 dtout%isecur             = dtin%isecur
 dtout%istatimg           = dtin%istatimg
 dtout%istatr             = dtin%istatr
 dtout%istatshft          = dtin%istatshft
 dtout%ixc                = dtin%ixc
 dtout%ixcpositron        = dtin%ixcpositron
 dtout%jdtset             = dtin%jdtset
 dtout%jellslab           = dtin%jellslab
 dtout%kptopt             = dtin%kptopt
 dtout%kssform            = dtin%kssform
 dtout%ldgapp             = dtin%ldgapp
 dtout%localrdwf          = dtin%localrdwf
 dtout%maxnsym            = dtin%maxnsym
 dtout%mband              = dtin%mband
 dtout%mffmem             = dtin%mffmem
 dtout%mgfft              = dtin%mgfft
 dtout%mgfftdg            = dtin%mgfftdg
 dtout%mkmem              = dtin%mkmem
 dtout%mkqmem             = dtin%mkqmem
 dtout%mk1mem             = dtin%mk1mem
 dtout%mpw                = dtin%mpw
 dtout%mqgrid             = dtin%mqgrid
 dtout%mqgriddg           = dtin%mqgriddg
 dtout%natom              = dtin%natom
 dtout%natrd              = dtin%natrd
 dtout%natsph             = dtin%natsph
 dtout%natpawu            = dtin%natpawu
 dtout%natvshift          = dtin%natvshift
 dtout%nbandsus           = dtin%nbandsus
 dtout%nbdblock           = dtin%nbdblock
 dtout%nbdbuf             = dtin%nbdbuf
 dtout%nberry             = dtin%nberry
 dtout%nbandkss           = dtin%nbandkss
 dtout%nconeq             = dtin%nconeq
 dtout%nctime             = dtin%nctime
 dtout%ndtset             = dtin%ndtset
 dtout%ndyson             = dtin%ndyson
 dtout%ndynimage          = dtin%ndynimage
 dtout%nfft               = dtin%nfft
 dtout%nfftdg             = dtin%nfftdg
 dtout%nfreqim            = dtin%nfreqim
 dtout%nfreqre            = dtin%nfreqre
 dtout%nfreqsp            = dtin%nfreqsp
 dtout%nfreqsus           = dtin%nfreqsus
 dtout%ngroup_rf          = dtin%ngroup_rf
 dtout%nimage             = dtin%nimage
 dtout%nkptgw             = dtin%nkptgw
 dtout%nkpt               = dtin%nkpt
 dtout%nline              = dtin%nline
 dtout%nnsclo             = dtin%nnsclo
 dtout%nomegasf           = dtin%nomegasf
 dtout%nomegasi           = dtin%nomegasi
 dtout%nomegasrd          = dtin%nomegasrd
 dtout%npband             = dtin%npband
 dtout%npfft              = dtin%npfft
 dtout%npimage            = dtin%npimage
 dtout%npkpt              = dtin%npkpt
 dtout%npspinor           = dtin%npspinor
 dtout%npsp               = dtin%npsp
 dtout%npspalch           = dtin%npspalch
 dtout%npulayit           = dtin%npulayit
 dtout%npweps             = dtin%npweps
 dtout%npwkss             = dtin%npwkss
 dtout%npwsigx            = dtin%npwsigx
 dtout%npwwfn             = dtin%npwwfn
 dtout%nqpt               = dtin%nqpt
 dtout%nqptdm             = dtin%nqptdm
 dtout%nscforder          = dtin%nscforder
 dtout%nsheps             = dtin%nsheps
 dtout%nshiftk            = dtin%nshiftk
 dtout%nshsigx            = dtin%nshsigx
 dtout%nshwfn             = dtin%nshwfn
 dtout%nspden             = dtin%nspden
 dtout%nspinor            = dtin%nspinor
 dtout%nsppol             = dtin%nsppol
 dtout%nstep              = dtin%nstep
 dtout%nsym               = dtin%nsym
 dtout%ntime              = dtin%ntime
 dtout%ntimimage          = dtin%ntimimage
 dtout%ntypalch           = dtin%ntypalch
 dtout%ntypat             = dtin%ntypat
 dtout%ntyppure           = dtin%ntyppure
 dtout%nwfshist           = dtin%nwfshist
 dtout%occopt             = dtin%occopt
 dtout%optcell            = dtin%optcell
 dtout%optdriver          = dtin%optdriver
 dtout%optforces          = dtin%optforces
 dtout%optfreqsus         = dtin%optfreqsus
 dtout%optnlxccc          = dtin%optnlxccc
 dtout%optstress          = dtin%optstress
 dtout%ortalg             = dtin%ortalg
 dtout%paral_kgb          = dtin%paral_kgb
 dtout%paral_rf           = dtin%paral_rf
 dtout%pawcpxocc          = dtin%pawcpxocc
 dtout%pawcross           = dtin%pawcross
 dtout%pawlcutd           = dtin%pawlcutd
 dtout%pawlmix            = dtin%pawlmix
 dtout%pawmixdg           = dtin%pawmixdg
 dtout%pawnhatxc          = dtin%pawnhatxc
 dtout%pawnphi            = dtin%pawnphi
 dtout%pawntheta          = dtin%pawntheta
 dtout%pawnzlm            = dtin%pawnzlm
 dtout%pawoptmix          = dtin%pawoptmix
 dtout%pawprtden          = dtin%pawprtden
 dtout%pawprtdos          = dtin%pawprtdos
 dtout%pawprtvol          = dtin%pawprtvol
 dtout%pawprtwf           = dtin%pawprtwf
 dtout%pawprt_k           = dtin%pawprt_k
 dtout%pawprt_b           = dtin%pawprt_b
 dtout%pawspnorb          = dtin%pawspnorb
 dtout%pawstgylm          = dtin%pawstgylm
 dtout%pawsushat          = dtin%pawsushat
 dtout%pawusecp           = dtin%pawusecp
 dtout%pawujat            = dtin%pawujat
 dtout%macro_uj           = dtin%macro_uj
 dtout%pawujrad           = dtin%pawujrad
 dtout%pawujv             = dtin%pawujv
 dtout%pawxcdev           = dtin%pawxcdev
 dtout%pitransform        = dtin%pitransform
 dtout%positron           = dtin%positron
 dtout%posnstep           = dtin%posnstep
 dtout%ppmodel            = dtin%ppmodel
 dtout%prepanl            = dtin%prepanl
 dtout%prepgkk            = dtin%prepgkk
 dtout%prepscphon         = dtin%prepscphon
 dtout%prtbbb             = dtin%prtbbb
 dtout%prtbltztrp         = dtin%prtbltztrp
 dtout%prtcif             = dtin%prtcif
 dtout%prtcml             = dtin%prtcml
 dtout%prtcs              = dtin%prtcs
 dtout%prtden             = dtin%prtden
 dtout%prtdensph          = dtin%prtdensph
 dtout%prtdipole          = dtin%prtdipole
 dtout%prtdos             = dtin%prtdos
 dtout%prtdosm            = dtin%prtdosm
 dtout%prtefg             = dtin%prtefg
 dtout%prteig             = dtin%prteig
 dtout%prtelf             = dtin%prtelf
 dtout%prtfc              = dtin%prtfc
 dtout%prtfsurf           = dtin%prtfsurf
 dtout%prtgden            = dtin%prtgden
 dtout%prtgeo             = dtin%prtgeo
 dtout%prtgkk             = dtin%prtgkk
 dtout%prtkden            = dtin%prtkden
 dtout%prtkpt             = dtin%prtkpt
 dtout%prtlden            = dtin%prtlden
 dtout%prtnabla           = dtin%prtnabla
 dtout%prtnest            = dtin%prtnest
 dtout%prtposcar          = dtin%prtposcar
 dtout%prtpot             = dtin%prtpot
 dtout%prtspcur           = dtin%prtspcur
 dtout%prtstm             = dtin%prtstm
 dtout%prtvha             = dtin%prtvha
 dtout%prtvhxc            = dtin%prtvhxc
 dtout%prtvol             = dtin%prtvol
 dtout%prtvolimg          = dtin%prtvolimg
 dtout%prtvxc             = dtin%prtvxc
 dtout%prtwant            = dtin%prtwant
 dtout%prtwf              = dtin%prtwf
 dtout%prtxml             = dtin%prtxml
 dtout%prt1dm             = dtin%prt1dm
 dtout%ptgroupma          = dtin%ptgroupma
 dtout%random_atpos       = dtin%random_atpos
 dtout%recgratio          = dtin%recgratio
 dtout%recnpath           = dtin%recnpath
 dtout%recnrec            = dtin%recnrec
 dtout%recptrott          = dtin%recptrott
 dtout%rectesteg          = dtin%rectesteg
 dtout%rcut               = dtin%rcut
 dtout%rdmnb              = dtin%rdmnb
 dtout%restartxf          = dtin%restartxf
 dtout%rfasr              = dtin%rfasr
 dtout%rfddk              = dtin%rfddk
 dtout%rfelfd             = dtin%rfelfd
 dtout%rfmeth             = dtin%rfmeth
 dtout%rfphon             = dtin%rfphon
 dtout%rfstrs             = dtin%rfstrs
 dtout%rfuser             = dtin%rfuser
 dtout%rf1elfd            = dtin%rf1elfd
 dtout%rf1phon            = dtin%rf1phon
 dtout%rf2elfd            = dtin%rf2elfd
 dtout%rf2phon            = dtin%rf2phon
 dtout%rf3elfd            = dtin%rf3elfd
 dtout%rf3phon            = dtin%rf3phon
 dtout%rhoqpmix           = dtin%rhoqpmix
 dtout%signperm           = dtin%signperm
 dtout%slabwsrad          = dtin%slabwsrad
 dtout%slabzbeg           = dtin%slabzbeg
 dtout%slabzend           = dtin%slabzend
 dtout%smdelta            = dtin%smdelta
 dtout%spgaxor            = dtin%spgaxor
 dtout%spgorig            = dtin%spgorig
 dtout%spgroup            = dtin%spgroup
 dtout%spmeth             = dtin%spmeth
 dtout%suskxcrs           = dtin%suskxcrs
 dtout%symchi             = dtin%symchi
 dtout%symmorphi          = dtin%symmorphi
 dtout%symsigma           = dtin%symsigma
 dtout%td_mexcit          = dtin%td_mexcit
 dtout%tfkinfunc          = dtin%tfkinfunc
 dtout%timopt             = dtin%timopt
 dtout%use_gpu_cuda       = dtin%use_gpu_cuda
 dtout%use_slk            = dtin%use_slk
 dtout%useexexch          = dtin%useexexch
 dtout%usedmatpu          = dtin%usedmatpu
 dtout%usedmft            = dtin%usedmft
 dtout%usekden            = dtin%usekden
 dtout%usepaw             = dtin%usepaw
 dtout%usepawu            = dtin%usepawu
 dtout%userec             = dtin%userec
 dtout%useria             = dtin%useria
 dtout%userib             = dtin%userib
 dtout%useric             = dtin%useric
 dtout%userid             = dtin%userid
 dtout%userie             = dtin%userie
 dtout%usewvl             = dtin%usewvl
 dtout%usexcnhat          = dtin%usexcnhat
 dtout%useylm             = dtin%useylm
 dtout%vacnum             = dtin%vacnum
 dtout%vdw_nfrag          = dtin%vdw_nfrag
 dtout%vdw_xc             = dtin%vdw_xc
 dtout%wfoptalg           = dtin%wfoptalg
 dtout%w90iniprj          = dtin%w90iniprj
 dtout%w90prtunk          = dtin%w90prtunk
 dtout%xclevel            = dtin%xclevel
 dtout%xc_denpos          = dtin%xc_denpos

!Copy allocated integer arrays from dtin to dtout
 dtout%bdberry(:)         = dtin%bdberry(:)
 dtout%cd_subset_freq(:)  = dtin%cd_subset_freq(:)
 dtout%kptrlatt(:,:)      = dtin%kptrlatt(:,:)
 dtout%ngfft(:)           = dtin%ngfft(:)
 dtout%ngfftdg(:)         = dtin%ngfftdg(:)
 dtout%nloalg(:)          = dtin%nloalg(:)
 dtout%qprtrb(:)          = dtin%qprtrb(:)
 dtout%rfatpol(:)         = dtin%rfatpol(:)
 dtout%rfdir(:)           = dtin%rfdir(:)
 dtout%rf1atpol(:)        = dtin%rf1atpol(:)
 dtout%rf1dir(:)          = dtin%rf1dir(:)
 dtout%rf2atpol(:)        = dtin%rf2atpol(:)
 dtout%rf2dir(:)          = dtin%rf2dir(:)
 dtout%rf3atpol(:)        = dtin%rf3atpol(:)
 dtout%rf3dir(:)          = dtin%rf3dir(:)
 dtout%scphon_supercell(:)= dtin%scphon_supercell(:)
 dtout%supercell(:)       = dtin%supercell(:)
 dtout%vdw_supercell(:)   = dtin%vdw_supercell(:)
 dtout%vdw_typfrag(:)     = dtin%vdw_typfrag(:)
!Copy reals from dtin to dtout
 dtout%boxcutmin          = dtin%boxcutmin
 dtout%bxctmindg          = dtin%bxctmindg
 dtout%cd_halfway_freq    = dtin%cd_halfway_freq
 dtout%cd_max_freq        = dtin%cd_max_freq
 dtout%charge             = dtin%charge
 dtout%cpus               = dtin%cpus
 dtout%diecut             = dtin%diecut
 dtout%diegap             = dtin%diegap
 dtout%dielam             = dtin%dielam
 dtout%dielng             = dtin%dielng
 dtout%diemac             = dtin%diemac
 dtout%diemix             = dtin%diemix
 dtout%diemixmag          = dtin%diemixmag
 dtout%dilatmx            = dtin%dilatmx
 dtout%dosdeltae          = dtin%dosdeltae
 dtout%dtion              = dtin%dtion
 dtout%ecut               = dtin%ecut
 dtout%ecuteps            = dtin%ecuteps
 dtout%ecutsigx           = dtin%ecutsigx
 dtout%ecutsm             = dtin%ecutsm
 dtout%ecutwfn            = dtin%ecutwfn
 dtout%effmass            = dtin%effmass
 dtout%elph2_imagden      = dtin%elph2_imagden
 dtout%eshift             = dtin%eshift
 dtout%esmear             = dtin%esmear
 dtout%exchmix            = dtin%exchmix
 dtout%fband              = dtin%fband
 dtout%fixmom             = dtin%fixmom
 dtout%freqsusin          = dtin%freqsusin
 dtout%freqsuslo          = dtin%freqsuslo
 dtout%friction           = dtin%friction
 dtout%fxcartfactor       = dtin%fxcartfactor
 dtout%gw_eet_scale       = dtin%gw_eet_scale
 dtout%kptnrm             = dtin%kptnrm
 dtout%kptrlen            = dtin%kptrlen
 dtout%bmass              = dtin%bmass
 dtout%nnos               = dtin%nnos
 dtout%mdwall             = dtin%mdwall
 dtout%nelect             = dtin%nelect
 dtout%noseinert          = dtin%noseinert
 dtout%omegasimax         = dtin%omegasimax
 dtout%omegasrdmax        = dtin%omegasrdmax
 dtout%pawecutdg          = dtin%pawecutdg
 dtout%pawovlp            = dtin%pawovlp
 dtout%posocc             = dtin%posocc
 dtout%postoldfe          = dtin%postoldfe
 dtout%postoldff          = dtin%postoldff
 dtout%ppmfrq             = dtin%ppmfrq
 dtout%recrcut            = dtin%recrcut
 dtout%recefermi          = dtin%recefermi
 dtout%rectolden          = dtin%rectolden
 dtout%scphon_temp        = dtin%scphon_temp
 dtout%sciss              = dtin%sciss
 dtout%soenergy           = dtin%soenergy
 dtout%spbroad            = dtin%spbroad
 dtout%spnorbscl          = dtin%spnorbscl
 dtout%stmbias            = dtin%stmbias
 dtout%strfact            = dtin%strfact
 dtout%strprecon          = dtin%strprecon
 dtout%tl_radius          = dtin%tl_radius
 dtout%tl_nprccg          = dtin%tl_nprccg
 dtout%td_maxene          = dtin%td_maxene
 dtout%toldfe             = dtin%toldfe
 dtout%toldff             = dtin%toldff
 dtout%tolimg             = dtin%tolimg
 dtout%tolrff             = dtin%tolrff
 dtout%tolmxf             = dtin%tolmxf
 dtout%tolsym             = dtin%tolsym
 dtout%tolvrs             = dtin%tolvrs
 dtout%tolwfr             = dtin%tolwfr
 dtout%tphysel            = dtin%tphysel
 dtout%tsmear             = dtin%tsmear
 dtout%userra             = dtin%userra
 dtout%userrb             = dtin%userrb
 dtout%userrc             = dtin%userrc
 dtout%userrd             = dtin%userrd
 dtout%userre             = dtin%userre
 dtout%vacwidth           = dtin%vacwidth
 dtout%vis                = dtin%vis
 dtout%wvl_hgrid          = dtin%wvl_hgrid
 dtout%wvl_crmult         = dtin%wvl_crmult
 dtout%wvl_frmult         = dtin%wvl_frmult
 dtout%wvl_nprccg         = dtin%wvl_nprccg
 dtout%zcut               = dtin%zcut

!Copy allocated real arrays from dtin to dtout
 dtout%boxcenter(:)       = dtin%boxcenter(:)
 dtout%bfield(:)          = dtin%bfield(:)
 dtout%efield(:)          = dtin%efield(:)
 dtout%genafm(:)          = dtin%genafm(:)
 dtout%goprecprm(:)       = dtin%goprecprm(:)
 dtout%mdtemp(:)          = dtin%mdtemp(:)
 dtout%qptn(:)            = dtin%qptn(:)
 dtout%strtarget(:)       = dtin%strtarget(:)
 dtout%vcutgeo(:)         = dtin%vcutgeo(:)
 dtout%vprtrb(:)          = dtin%vprtrb(:)
 dtout%zeemanfield(:)     = dtin%zeemanfield(:)


!This list of variables are allocated using the maximum size
!over all datasets, so we just read the allocated size and reproduce
!it.
!atvshift(mxnatvshift,mxnsppol,mxnatom))
!bdgw(2,mxnkptgw))
!dynimage(mxnimage))
!kpt(3,mxnkpt))
!kptgw(3,mxnkptgw))
!kptns(3,mxnkpt))
!iatsph(mxnatsph))
!istwfk(mxnkpt))
!nband(mxnkpt*mxnsppol))
!occ_orig(mxmband_upper*mxnkpt*mxnsppol))
!qmass(mxnnos)
!qptdm(3,mxnqptdm))
!symafm(mxnsym))
!symrel(3,3,mxnsym))
!tnons(3,mxnsym))
!wtatcon(3,mxnatom,mxnconeq))
!wtk(mxnkpt))
!iatfix(3,mxnatom))
!spinat(3,mxnatom))
!typat(mxnatom))
!vel_orig(3,mxnatom,mxnimage))
!xred_orig(3,mxnatom,mxnimage))
!algalch(mxntypat))
!amu(mxntypat))
!dmatpawu(2*mxlpawu+1,*mxlpawu+1,mxnsppol*mxnspinor,mxnatpawu)
!gw_qlwl(3,gw_nqlwl)
!densty(mxntypat,4))
!jpawu(mxntypat))
!lpawu(mxntypat))
!mixalch(npsp,mxntypat))
!so_psp(npsp)
!upawu(mxntypat))
!ziontypat(mxntypat))
!znucl(npsp))

!#define DTSET_AREA_NAME_STR
!#define DTSET_AREA_NAME
!#define DTSET_AREA_SIZE1
!if (associated(dtin%DTSET_AREA_NAME)) then
!write(name,*) DTSET_AREA_NAME_STR
!call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
!allocate (dtout%DTSET_AREA_NAME(chosen_size1))
!dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
!end if
!#undef DTSET_AREA_NAME
!#undef DTSET_AREA_SIZE1
!#undef DTSET_AREA_NAME_STR
!
!#define DTSET_AREA_NAME_STR
!#define DTSET_AREA_NAME
!#define DTSET_AREA_SIZE1
!#define DTSET_AREA_SIZE2
!if (associated(dtin%DTSET_AREA_NAME)) then
!write(name,*) DTSET_AREA_NAME_STR
!call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
!call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
!allocate (dtout%DTSET_AREA_NAME(chosen_size1,chosen_size2))
!dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
!end if
!#undef DTSET_AREA_NAME
!#undef DTSET_AREA_SIZE1
!#undef DTSET_AREA_SIZE2
!#undef DTSET_AREA_NAME_STR


!!#define DEV_HAVE_NEWDTCOPY

!Use m_copy for performing the copy of the pointee.
!Rationale: it is shorter and much more readable!
!Ok, it does not perform any check on the size but this task should not be done here!
!The only purpose of this routine is performing a complete copy of the input structure.
!One should never change the basic dimensions of dtset during the execution.
!If the developer has to change the content of dtset before calling an abinit routine.
!that means that the interface of the routine should be changed not dtset.
!for example one might pass a smaller datatype with all the required information instead of dtset.
#ifdef DEV_HAVE_NEWDTCOPY

#include "copy_dt.finc"
#else

#define DTSET_AREA_NAME_STR "acell_orig"
#define DTSET_AREA_NAME  acell_orig
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  nimage
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,3,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%acell_orig) allocate(dtout%acell_orig(3,dtin%nimage))
!dtout%acell_orig(:,:)         = dtin%acell_orig(:,:)



#define DTSET_AREA_NAME_STR "algalch"
#define DTSET_AREA_NAME  algalch
#define DTSET_AREA_SIZE1  ntypalch
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%algalch)) allocate(dtout%algalch(dtin%ntypalch))
!dtout%algalch(:)         = dtin%algalch(:)


#define DTSET_AREA_NAME_STR "atvshift"
#define DTSET_AREA_NAME  atvshift
#define DTSET_AREA_SIZE1  natvshift
#define DTSET_AREA_SIZE2  nsppol
#define DTSET_AREA_SIZE3  natom
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   call tells_sizes(chosen_size3,name,3,dtin%DTSET_AREA_SIZE3,size(dtin%DTSET_AREA_NAME,3))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2,chosen_size3))
   dtout%DTSET_AREA_NAME(:,:,:)=dtin%DTSET_AREA_NAME(:,:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_SIZE3
#undef DTSET_AREA_NAME_STR
!allocate(dtout%atvshift(dtin%natvshift, dtin%nsppol, dtin%natom))
!dtout%atvshift(:,:,:)     = dtin%atvshift(:,:,:)



#define DTSET_AREA_NAME_STR "bdgw"
#define DTSET_AREA_NAME  bdgw
#define DTSET_AREA_SIZE1  2
#define DTSET_AREA_SIZE2  nkptgw
#define DTSET_AREA_SIZE3  nsppol
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,2,                    size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   call tells_sizes(chosen_size3,name,3,dtin%DTSET_AREA_SIZE3,size(dtin%DTSET_AREA_NAME,3))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2,chosen_size3))
   dtout%DTSET_AREA_NAME=dtin%DTSET_AREA_NAME(1:chosen_size1,1:chosen_size2,1:chosen_size3)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_SIZE3
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%bdgw)) allocate(dtout%bdgw(2, dtin%nkptgw,dtin%nsppol))
!dtout%bdgw          = dtin%bdgw

#define DTSET_AREA_NAME_STR "dynimage"
#define DTSET_AREA_NAME  dynimage
#define DTSET_AREA_SIZE1  nimage
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%dynimage)) allocate(dtout%dynimage(dtin%nimage))
!dtout%dynimage(:,:)        = dtin%dynimage(:,:)

#define DTSET_AREA_NAME_STR "cd_imfrqs"
#define DTSET_AREA_NAME  cd_imfrqs
#define DTSET_AREA_SIZE1  cd_custom_imfrqs
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%cd_imfrqs)) allocate(dtout%cd_imfrqs(dtin%cd_custom_imfrqs))
!dtout%cd_imfrqs(:)        = dtin%cd_imfrqs(:)

#define DTSET_AREA_NAME_STR "gw_freqsp"
#define DTSET_AREA_NAME  gw_freqsp
#define DTSET_AREA_SIZE1  gw_custom_freqsp
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%gw_freqsp)) allocate(dtout%gw_freqsp(dtin%gw_custom_freqsp))
!dtout%gw_freqsp(:)        = dtin%gw_freqsp(:)

#define DTSET_AREA_NAME_STR "iatfix"
#define DTSET_AREA_NAME  iatfix
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  natom
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%iatfix)) allocate(dtout%iatfix(3, dtin%natom))
!dtout%iatfix(:,:)        = dtin%iatfix(:,:)

#define DTSET_AREA_NAME_STR "iatsph"
#define DTSET_AREA_NAME  iatsph
#define DTSET_AREA_SIZE1  natsph
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%iatsph)) allocate(dtout%iatsph(dtin%natsph))
!dtout%iatsph(:)          = dtin%iatsph(:)

#define DTSET_AREA_NAME_STR "istwfk"
#define DTSET_AREA_NAME  istwfk
#define DTSET_AREA_SIZE1  nkpt
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%istwfk)) allocate(dtout%istwfk(dtin%nkpt))
!dtout%istwfk(:)          = dtin%istwfk(:)

#define DTSET_AREA_NAME_STR "kberry"
#define DTSET_AREA_NAME  kberry
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  nberry
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%kberry)) allocate(dtout%kberry(3, dtin%nberry))
!dtout%kberry(:,:)        = dtin%kberry(:,:)

#define DTSET_AREA_NAME_STR "lpawu"
#define DTSET_AREA_NAME  lpawu
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!if (associated(dtin%lpawu)) allocate(dtout%lpawu(dtin%ntypat))
!dtout%lpawu(:)        = dtin%lpawu(:)

#define DTSET_AREA_NAME_STR "lexexch"
#define DTSET_AREA_NAME  lexexch
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR

#define DTSET_AREA_NAME_STR "nband"
#define DTSET_AREA_NAME  nband
#define DTSET_AREA_SIZE1  nsppol
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%nband(dtin%nkpt * dtin%nsppol))
!dtout%nband(:)           = dtin%nband(:)

#define DTSET_AREA_NAME_STR "prtatlist"
#define DTSET_AREA_NAME  prtatlist
#define DTSET_AREA_SIZE1  natom
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR

#define DTSET_AREA_NAME_STR "qmass"
#define DTSET_AREA_NAME  qmass
#define DTSET_AREA_SIZE1  nnos
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%qmass(dtin%nnos))
!dtout%qmass(:)=dtin%qmass(:)

#define DTSET_AREA_NAME_STR "rprim_orig"
#define DTSET_AREA_NAME  rprim_orig
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  3
#define DTSET_AREA_SIZE3  nimage
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   call tells_sizes(chosen_size3,name,3,dtin%DTSET_AREA_SIZE3,size(dtin%DTSET_AREA_NAME,3))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2,chosen_size3))
   dtout%DTSET_AREA_NAME(:,:,:)=dtin%DTSET_AREA_NAME(:,:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_SIZE3
#undef DTSET_AREA_NAME_STR
!allocate(dtout%rprim_orig(3, 3, dtin%nimage))
!dtout%rprim_orig(:,:,:)      = dtin%rprim_orig(:,:,:)


#define DTSET_AREA_NAME_STR "rprimd_orig"
#define DTSET_AREA_NAME  rprimd_orig
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  3
#define DTSET_AREA_SIZE3  nimage
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   call tells_sizes(chosen_size3,name,3,dtin%DTSET_AREA_SIZE3,size(dtin%DTSET_AREA_NAME,3))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2,chosen_size3))
   dtout%DTSET_AREA_NAME(:,:,:)=dtin%DTSET_AREA_NAME(:,:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_SIZE3
#undef DTSET_AREA_NAME_STR
!allocate(dtout%rprimd_orig(3, 3, dtin%nimage))
!dtout%rprimd_orig(:,:,:)      = dtin%rprimd_orig(:,:,:)

#define DTSET_AREA_NAME_STR "so_psp"
#define DTSET_AREA_NAME  so_psp
#define DTSET_AREA_SIZE1  npsp
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%so_psp(dtin%npsp))
!dtout%so_psp(:)        = dtin%so_psp(:)

#define DTSET_AREA_NAME_STR "symafm"
#define DTSET_AREA_NAME  symafm
#define DTSET_AREA_SIZE1  nsym
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%symafm(dtin%nsym))
!dtout%symafm(:)          = dtin%symafm(:)

#define DTSET_AREA_NAME_STR "symrel"
#define DTSET_AREA_NAME  symrel
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  3
#define DTSET_AREA_SIZE3  nsym
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   call tells_sizes(chosen_size3,name,3,dtin%DTSET_AREA_SIZE3,size(dtin%DTSET_AREA_NAME,3))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2,chosen_size3))
   dtout%DTSET_AREA_NAME(:,:,:)=dtin%DTSET_AREA_NAME(:,:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_SIZE3
#undef DTSET_AREA_NAME_STR
!allocate(dtout%symrel(3, 3, dtin%nsym))
!dtout%symrel(:,:,:)      = dtin%symrel(:,:,:)

#define DTSET_AREA_NAME_STR "typat"
#define DTSET_AREA_NAME  typat
#define DTSET_AREA_SIZE1  natom
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%typat(dtin%natom))
!dtout%typat(:)           = dtin%typat(:)

!Allocate and copy real pointers
!same syntax as the one for integer pointers

#define DTSET_AREA_NAME_STR "corecs"
#define DTSET_AREA_NAME  corecs
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR

#define DTSET_AREA_NAME_STR "ptcharge"
#define DTSET_AREA_NAME  ptcharge
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR

#define DTSET_AREA_NAME_STR "quadmom"
#define DTSET_AREA_NAME  quadmom
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR

#define DTSET_AREA_NAME_STR "amu"
#define DTSET_AREA_NAME  amu
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%amu(dtin%ntypat))
!dtout%amu(:)             = dtin%amu(:)

#define DTSET_AREA_NAME_STR "densty"
#define DTSET_AREA_NAME  densty
#define DTSET_AREA_SIZE1  ntypat
#define DTSET_AREA_SIZE2  4
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
!  DEBUG
!  write(std_out,*)' dtsetcopy : allocated densty=',associated(dtout%densty)
!  ENDDEBUG
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%densty(dtin%ntypat, 4))
!dtout%densty(:,:)        = dtin%densty(:,:)

#define DTSET_AREA_NAME_STR "dmatpawu"
#define DTSET_AREA_NAME  dmatpawu
#define DTSET_AREA_SIZE4  natpawu
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,size(dtin%DTSET_AREA_NAME,1),size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,3,size(dtin%DTSET_AREA_NAME,3),size(dtin%DTSET_AREA_NAME,3))
   call tells_sizes(chosen_size3,name,4,dtin%DTSET_AREA_SIZE4,size(dtin%DTSET_AREA_NAME,4))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size1,chosen_size2,chosen_size3))
   dtout%DTSET_AREA_NAME(:,:,:,:)=dtin%DTSET_AREA_NAME(:,:,:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE4
#undef DTSET_AREA_NAME_STR
!allocate(dtout%dmatpawu(?2*lpawu+1?,?2*lpawu+1?,?nsppol*nspinor?,dtin%natpawu))
!dtout%dmatpawu(:,:)        = dtin%dmatpawu(:,:)

#define DTSET_AREA_NAME_STR "gw_qlwl"
#define DTSET_AREA_NAME  gw_qlwl
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  gw_nqlwl
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%gw_qlwl(3, dtin%gw_nqlwl))
!dtout%gw_qlwl(:,:)           = dtin%gw_qlwl(:,:)

#define DTSET_AREA_NAME_STR "jpawu"
#define DTSET_AREA_NAME  jpawu
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%jpawu(dtin%ntypat))
!dtout%jpawu(:)            = dtin%jpawu(:)

#define DTSET_AREA_NAME_STR "kpt"
#define DTSET_AREA_NAME  kpt
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  nkpt
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%kpt(3, dtin%nkpt))
!dtout%kpt(:,:)           = dtin%kpt(:,:)

#define DTSET_AREA_NAME_STR "kptgw"
#define DTSET_AREA_NAME  kptgw
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  nkptgw
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%kptgw(3, dtin%nkptgw))
!dtout%kptgw(:,:)         = dtin%kptgw(:,:)

#define DTSET_AREA_NAME_STR "kptns"
#define DTSET_AREA_NAME  kptns
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  nkpt
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%kptns(3, dtin%nkpt))
!dtout%kptns(:,:)         = dtin%kptns(:,:)

#define DTSET_AREA_NAME_STR "mixalch"
#define DTSET_AREA_NAME  mixalch
#define DTSET_AREA_SIZE1  npspalch
#define DTSET_AREA_SIZE2  ntypalch
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%mixalch(dtin%npspalch, dtin%ntypalch))
!dtout%mixalch(:,:)       = dtin%mixalch(:,:)

#define DTSET_AREA_NAME_STR "occ_orig"
#define DTSET_AREA_NAME  occ_orig
#define DTSET_AREA_SIZE1  (dtin%mband*dtin%nkpt*dtin%nsppol)
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%occ_orig(dtin%mband * dtin%nkpt * dtin%nsppol))
!dtout%occ_orig(:)        = dtin%occ_orig(:)

#define DTSET_AREA_NAME_STR "pimass"
#define DTSET_AREA_NAME  pimass
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%pimass(dtin%ntypat))
!dtout%pimass(:)          = dtin%pimass(:)

#define DTSET_AREA_NAME_STR "qptdm"
#define DTSET_AREA_NAME  qptdm
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  nqptdm
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%qptdm(3, dtin%nqptdm))
!dtout%qptdm(:,:)         = dtin%qptdm(:,:)

#define DTSET_AREA_NAME_STR "ratsph"
#define DTSET_AREA_NAME  ratsph
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR

#define DTSET_AREA_NAME_STR "shiftk"
#define DTSET_AREA_NAME  shiftk
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  nshiftk
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%shiftk(3, dtin%nshiftk))
!dtout%shiftk(:,:)        = dtin%shiftk(:,:)

#define DTSET_AREA_NAME_STR "spinat"
#define DTSET_AREA_NAME  spinat
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  natom
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%spinat(3, dtin%natom))
!dtout%spinat(:,:)        = dtin%spinat(:,:)

#define DTSET_AREA_NAME_STR "tnons"
#define DTSET_AREA_NAME  tnons
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  nsym
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2))
   dtout%DTSET_AREA_NAME(:,:)=dtin%DTSET_AREA_NAME(:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%tnons(3, dtin%nsym))
!dtout%tnons(:,:)         = dtin%tnons(:,:)

#define DTSET_AREA_NAME_STR "upawu"
#define DTSET_AREA_NAME  upawu
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%upawu(dtin%ntypat))
!dtout%upawu(:)            = dtin%upawu(:)

#define DTSET_AREA_NAME_STR "vel_orig"
#define DTSET_AREA_NAME  vel_orig
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  natom
#define DTSET_AREA_SIZE3  nimage
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   call tells_sizes(chosen_size3,name,3,dtin%DTSET_AREA_SIZE3,size(dtin%DTSET_AREA_NAME,3))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2,chosen_size3))
   dtout%DTSET_AREA_NAME(:,:,:)=dtin%DTSET_AREA_NAME(:,:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_SIZE3
#undef DTSET_AREA_NAME_STR
!allocate(dtout%vel_orig(3, dtin%natom, dtin%nimage))
!dtout%vel_orig(:,:,:)      = dtin%vel_orig(:,:,:)

#define DTSET_AREA_NAME_STR "wtatcon"
#define DTSET_AREA_NAME  wtatcon
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  natom
#define DTSET_AREA_SIZE3  nconeq
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   call tells_sizes(chosen_size3,name,3,dtin%DTSET_AREA_SIZE3,size(dtin%DTSET_AREA_NAME,3))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2,chosen_size3))
   dtout%DTSET_AREA_NAME(:,:,:)=dtin%DTSET_AREA_NAME(:,:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_SIZE3
#undef DTSET_AREA_NAME_STR
!allocate(dtout%wtatcon(3, dtin%natom, dtin%nconeq))
!dtout%wtatcon(:,:,:)     = dtin%wtatcon(:,:,:)

#define DTSET_AREA_NAME_STR "wtk"
#define DTSET_AREA_NAME  wtk
#define DTSET_AREA_SIZE1  nkpt
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%wtk(dtin%nkpt))
!dtout%wtk(:)             = dtin%wtk(:)

#define DTSET_AREA_NAME_STR "xred_orig"
#define DTSET_AREA_NAME  xred_orig
#define DTSET_AREA_SIZE1  3
#define DTSET_AREA_SIZE2  natom
#define DTSET_AREA_SIZE3  nimage
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   call tells_sizes(chosen_size2,name,2,dtin%DTSET_AREA_SIZE2,size(dtin%DTSET_AREA_NAME,2))
   call tells_sizes(chosen_size3,name,3,dtin%DTSET_AREA_SIZE3,size(dtin%DTSET_AREA_NAME,3))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1,chosen_size2,chosen_size3))
   dtout%DTSET_AREA_NAME(:,:,:)=dtin%DTSET_AREA_NAME(:,:,:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_SIZE2
#undef DTSET_AREA_NAME_STR
!allocate(dtout%xred_orig(3, dtin%natom, dtin%nimage))
!dtout%xred_orig(:,:,:)     = dtin%xred_orig(:,:,:)

#define DTSET_AREA_NAME_STR "ziontypat"
#define DTSET_AREA_NAME  ziontypat
#define DTSET_AREA_SIZE1  ntypat
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%ziontypat(dtin%ntypat))
!dtout%ziontypat(:)       = dtin%ziontypat(:)

#define DTSET_AREA_NAME_STR "znucl"
#define DTSET_AREA_NAME  znucl
#define DTSET_AREA_SIZE1  npsp
 if (associated(dtin%DTSET_AREA_NAME)) then
   write(name,*) DTSET_AREA_NAME_STR
   call tells_sizes(chosen_size1,name,1,dtin%DTSET_AREA_SIZE1,size(dtin%DTSET_AREA_NAME,1))
   ABI_ALLOCATE(dtout%DTSET_AREA_NAME,(chosen_size1))
   dtout%DTSET_AREA_NAME(:)=dtin%DTSET_AREA_NAME(:)
 end if
#undef DTSET_AREA_NAME
#undef DTSET_AREA_SIZE1
#undef DTSET_AREA_NAME_STR
!allocate(dtout%znucl(dtin%npsp))
!dtout%znucl(:)           = dtin%znucl(:)

#endif
!def DEV_HAVE_NEWDTCOPY

 DBG_EXIT("COLL")

end subroutine dtsetCopy
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/tells_sizes
!! NAME
!! tells_sizes
!!
!! FUNCTION
!! This routine compares the default size of an array with
!! the one that it actually has.
!! It returns the actual size as the chosen size for the allocation
!! of a similar array.
!! It gives a warning if the actual size and the default size are not equal.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  name         : name of the array
!!  index        : index associated with a dimension of the array
!!  default_size : default size of this array
!!  actual_size  : actual size of this array
!!
!! OUTPUT
!!  chosen_size : chosen size (=actual size) for the array that will be allocated
!!
!! PARENTS
!!      dtsetcopy
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine  tells_sizes(chosen_size,name,index,default_size,actual_size)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tells_sizes'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: actual_size,default_size,index
 integer,intent(out) :: chosen_size
 character(len=12),intent(in) :: name

!Local variables-------------------------------
!scalars
 logical :: warn
 character(len=500) :: message

!*********************************************************************

 warn = .FALSE.
 if (default_size .ne. actual_size) then
   warn = .TRUE.
 end if
 chosen_size=actual_size
 if (chosen_size < 0) then
   write(message, '(5a)' ) 'dtsetcopy : copying area ',name,&
&   ' whose size (',chosen_size,') appear to be uncorrectly defined'
   call wrtout(std_out,message,'COLL')
   chosen_size = 0
 end if
 if (warn) then
   write(message, '(3a,i6,a,i6,a,i6,a)')&
&   'dtsetcopy : copying area ',name, 'the actual size (',actual_size,&
&   ') of the index (',index,')  differs from its standard size (',default_size,')'
   call wrtout(std_out,message,'COLL')
 end if

end subroutine tells_sizes
!!***
