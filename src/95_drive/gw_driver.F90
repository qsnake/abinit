!{\src2tex{textfont=tt}}
!!****f* ABINIT/gw_driver
!! NAME
!! gw_driver
!!
!! FUNCTION
!!  Driver routine for self-consistent GW calculations (G0W0, G0W, GW0, GW) 
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! codvsn=code version
!! Dtfil<type(datafiles_type)>=variables related to files
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Pawang<type(pawang_type)>=paw angular mesh and related data
!! Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data
!! Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data
!! Psps<type(pseudopotential_type)>=variables related to pseudopotentials
!! rprim(3,3)=dimensionless real space primitive translations
!! xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  QP corrections are reportedd in the main abinit output files. 
!!  W and QP corrections are also saved in external files at each iteration.
!!
!! NOTES
!!  The particular approximation used for W and Sigma, the algorithms used 
!!  to evaluate these two quantities as well as the kind of self-consistency
!!  ("true" GW or model QPSCGW calculations) are specifies via the standard
!!  abinit routines. Here we only call screening and sigma inside a loop until
!!  the QP corrections are convergend withing the accuracy specified by the user.
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gw_driver(idtset,jdtset_,ndtset,acell,codvsn,filnam,Dtfil,Dtset,&
&  Pawang,Pawrad,Pawtab,Psps,rprim,xred)

 use m_profiling

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gw_driver'
 use interfaces_14_hidewrite
 use interfaces_27_toolbox_oop
 use interfaces_53_abiutil
 use interfaces_95_drive, except_this_one => gw_driver
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idtset,ndtset
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),intent(in) :: filnam(5)
 type(Datafiles_type),intent(inout) :: Dtfil
 type(Dataset_type),intent(inout) :: Dtset
 type(Pawang_type),intent(inout) :: Pawang
 type(Pseudopotential_type),intent(inout) :: Psps
!arrays
 integer,intent(in) :: jdtset_(0:ndtset)
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,Dtset%natom)
 type(Pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: in_screening=1,in_sigma=2
 integer :: istep,read_qps1,read_scr1,read_sus1,gw_sctype,ierr,gw_nstep
 logical :: converged,did_screening
 character(len=4) :: stringfile
 character(len=9) :: stringvar
 character(len=10) :: tag_it,tag_itm1
 character(len=500) :: msg
 character(len=fnlen) :: filchi0_it1,filqps_it1,filscr_it1,base_prefix

!*************************************************************************

 ierr=0  ! Check input consistency. 
 if ( ALL( (/Dtset%npweps,Dtset%nsheps/) <tol6 ).and. Dtset%ecuteps<tol6 ) then
   msg = ' One of the three variables ecuteps, npweps, or nsheps must be non-null.'
   MSG_WARNING(msg)
   ierr=ierr+1
 end if

 if ( ALL( (/Dtset%npwsigx,Dtset%nshsigx/) <tol6 ) .and. Dtset%ecutsigx<tol6 ) then
   msg = ' One of the three variables ecutsigx, npwsigx, or nshsigx must be non-null.'
   MSG_WARNING(msg)
   ierr=ierr+1
 end if

 if (ierr/=0) then
   msg = " Cannot continue, check input file"
   MSG_ERROR(msg)
 end if
!
!NOTES:
!
!* In screening one writes or reads the following files containing data that might change due to self-consistency:
!
!OUTPUT:
!1) fnameabo_scr --> SCR file
!2) fnameabo_sus --> SUS file
!if nomegaer>2:
!3) fnameabo_eelf     -->  EELF
!4) fnameabo_em1_lf   -->  Absorption spectra with local fields effects
!5) fnameabo_em1_nlf  -->  Absorption spectra without local fields effects
!
!INPUT:
!1) fnameabi_qps        --> QP band structure used to update W.
!
!* In sigma one writes or reads the following files containing data that might change due to self-consistency:
!
!INPUT:
!1) fnameabi_qps                  --> QP band structure used to update W.
!2) fnameabi_scr or fnameabi_sus  --> W
!
!OUTPUT (Only files important for self-consistency or for post-processing are reported)
!1) fnameabo_qps                  --> New QP band structure file
!2) fnameabo_qp_den               --> New QP density
!3) fnameabo_qp_pawden            --> New full QP density
!4) fnameabo_gw, fnameabo_sig,    --> QP corrections, Sigma(w), Spectral function.
!fnameabo_sgr, fnameabo_sgm

!Type of self-consistency and number of iterations.
 gw_sctype = Dtset%gw_sctype
 gw_nstep  = Dtset%gw_nstep; if (gw_sctype==GWSC_one_shot) gw_nstep=1

 base_prefix = Dtfil%filnam_ds(4)

 do istep=1,gw_nstep 
   
   call int2char(istep  ,tag_it  )
   call int2char(istep-1,tag_itm1)

!  Rename input SCR, SUSC and QPS files.
   Dtfil%fnameabi_scr = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_itm1)//"_SCR"  
   Dtfil%fnameabi_sus = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_itm1)//"_SUS"  

   read_sus1=0 
   read_scr1=0
   read_qps1=0 

   if (istep==1) then ! Allow the user to read a previous (QPS,SCR,SUS) file but _only_ at the first step.
!    
!    Build the _QPS file name filqps_it1 according to getqps and irdqps. A default is available if getqps is 0
     stringfile='_QPS'; stringvar='qps'
     call mkfilename(filnam,filqps_it1,dtset%getqps,idtset,dtset%irdqps,jdtset_,ndtset,stringfile,stringvar,read_qps1)

     if (read_qps1/=0) then ! Change the default.
       Dtfil%fnameabi_qps=filqps_it1
     else                   ! Save the IT1 name for GW0 runs
!      filqps_it1=TRIM(filnam_ds(3))//'_QPS'
       filqps_it1=Dtfil%fnameabi_qps
     end if
!    
!    Build the _SCR file name filscr_it1 according to getscr and irdscr. A default is available if getscr is 0
     stringfile='_SCR' ; stringvar='scr'                                        
     call mkfilename(filnam,filscr_it1,dtset%getscr,idtset,dtset%irdscr,jdtset_,ndtset,stringfile,stringvar,read_scr1)

     if (read_scr1/=0) then  ! Change the default.
       Dtfil%fnameabi_scr=filscr_it1
     else                    ! Save the IT1 name for GW0 runs
!      filscr_it1=TRIM(filnam_ds(3))//'_SCR'
       filscr_it1=Dtfil%fnameabi_scr
     end if
!    
!    Build the _SUS file name filchi0_it1 according to getsuscep and irdsuscep. A default is available if getsuscep is 0
     stringfile='_SUS' ; stringvar='sus'
     call mkfilename(filnam,filchi0_it1,dtset%getsuscep,idtset,dtset%irdsuscep,jdtset_,ndtset,stringfile,stringvar,read_sus1)

     if (read_sus1/=0) then ! Change the default.
       Dtfil%fnameabi_sus=filchi0_it1
     else                   ! Save the IT1 name for GW0 runs
!      filchi0_it1=TRIM(filnam_ds(3))//'_SUS'
       filchi0_it1 =Dtfil%fnameabi_sus
     end if
   end if

   did_screening=.FALSE.
   if (make_screening()) then

!    Change the prefix used for generic output files of screening.
     Dtfil%filnam_ds(4) = TRIM(base_prefix)//"_GWSCR_IT"//TRIM(tag_it)

!    Rename the output files produced by screening.
     Dtfil%fnameabo_scr     = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_SCR"  
     Dtfil%fnameabo_sus     = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_SUS"  
     Dtfil%fnameabo_eelf    = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_EELF"  
     Dtfil%fnameabo_em1_lf  = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_EM1_LF"  
     Dtfil%fnameabo_em1_nlf = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_EM1_NLF"  

     call show_gwfiles(in_screening,header="Calling screening",unit=std_out)
     call show_gwfiles(in_screening,header="Calling screening",unit=ab_out)

     call screening(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim)

!    Hack: when screening returns, %ecuteps, %npweps and %nsheps 
!    have been set to a non-null value and setshells complains in setup_sigma.
     Dtset%ecuteps=zero; Dtset%nsheps =zero
     did_screening = .TRUE.
   end if

!  Connect screening to sigma.
   Dtfil%fnameabi_scr =Dtfil%fnameabo_scr
   Dtfil%fnameabi_sus =Dtfil%fnameabo_sus
!  We are reading an existing file or we skipped screening ; or GW0 starting from previous file.
   if ( (istep==1 .and. .not.did_screening) .or.          & 
&   gw_sctype==GWSC_only_G .and. .not.did_screening) & 
&   then 
     Dtfil%fnameabi_scr =filscr_it1
     Dtfil%fnameabi_sus =filchi0_it1
   end if

   if (gw_sctype==GWSC_only_W) then ! sigma reads a fixed QPS file (either KS states or QPS from a previous GWSC run)
     MSG_WARNING("GWSC_only_W not tested")
     Dtfil%fnameabi_qps = filqps_it1 
   end if

!  Change the prefix used for generic output files of sigma.
   Dtfil%filnam_ds(4) = TRIM(base_prefix)//"_GWSIG_IT"//TRIM(tag_it)
   
!  Rename the output files produced by sigma
   Dtfil%fnameabo_qps    = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_QPS"  
   Dtfil%fnameabo_qp_den = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_QP_DEN"  
   Dtfil%fnameabo_qp_pawden = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_QP_PAWDEN"  
   Dtfil%fnameabo_qp_dos = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_QP_DOS"  
   Dtfil%fnameabo_qp_eig = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_QP_DB.nc"   ! TODO change name
   Dtfil%fnameabo_gw     = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_GW"         ! TODO change name  
   Dtfil%fnameabo_sig    = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_SIG"  
   Dtfil%fnameabo_sgr    = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_SGR"  
   Dtfil%fnameabo_sgm    = TRIM(base_prefix)//"_GW_IT"//TRIM(tag_it)//"_SGM"  

   call show_gwfiles(in_sigma,"Calling sigma",unit=std_out)
   call show_gwfiles(in_sigma,"Calling sigma",unit=ab_out)

   call sigma(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim,xred,converged)

   Dtfil%fnameabi_qps = Dtfil%fnameabo_qps ! output QPS --> input QPS for next iteration.
   
   if (converged.and.gw_sctype/=GWSC_one_shot) then 
     write(msg,'(a,i3,a)')" GW iterations converged after ",istep," iterations."
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")
     EXIT 
   end if
 end do

 if (.not.converged.and.gw_sctype/=GWSC_one_shot) then 
   write(msg,'(a,i3,a)')" WARNING: GW self-consistent iterations not converged after ",istep-1," steps."
   call wrtout(std_out,msg,"COLL")
   call wrtout(ab_out,msg,"COLL")
 end if

!Reinstate old value.
 Dtfil%filnam_ds(4) = base_prefix
 RETURN

 CONTAINS  !===================================================================
!!***

!!****f* gw_driver/make_screening
!! NAME
!!  make_screening
!!
!! FUNCTION
!!  Return true if screening has to be called for this iteration.
!!
!! SOURCE

logical function make_screening() 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_screening'
!End of the abilint section

 implicit none

!Local variables-------------------------------
!scalars
 logical :: have_W1file 

! *************************************************************************

   have_W1file = (read_sus1/=0.or.read_scr1/=0) ! Do we have an initial (SCR|SUSC) file?

   if (gw_sctype==GWSC_one_shot    ) make_screening = (  .not.have_W1file .and. istep==1) 
   if (gw_sctype==GWSC_only_G      ) make_screening = (  .not.have_W1file .and. istep==1) 
   if (gw_sctype==GWSC_only_W      ) make_screening = ( (.not.have_W1file .and. istep==1) .or. (istep>1) ) 
   if (gw_sctype==GWSC_both_G_and_W) make_screening = ( (.not.have_W1file .and. istep==1) .or. (istep>1) ) 

 end function make_screening
!!***

!----------------------------------------------------------------------

!!****f* gw_driver/show_gwfiles
!! NAME
!!  show_gwfiles
!!
!! FUNCTION
!!  printout of the files used in screening or sigma for a particular iteration.
!!
!! PARENTS
!!      gw_driver
!!
!! CHILDREN
!!
!! SOURCE

   subroutine show_gwfiles(where_iam,header,unit,mode_paral) 

 use m_profiling

   use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'show_gwfiles'
!End of the abilint section

   implicit none

!  Arguments ------------------------------------
!  scalars
   integer,intent(in) :: where_iam
   integer,optional,intent(in) :: unit
   character(len=*),optional,intent(in) :: header
   character(len=4),optional,intent(in) :: mode_paral 

!  Local variables-------------------------------
   integer :: my_unt
   character(len=4) :: my_mode
   character(len=500) :: msg      
!  ********************************************************************* 

   my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
   my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

#define LAZY_P(msg) call wrtout(my_unt,msg,my_mode)

   LAZY_P(" ")
   if (PRESENT(header)) then
     LAZY_P(header)
   end if

   msg = "None"
   if (gw_sctype==GWSC_one_shot    ) msg = " One-shot GW method " 
   if (gw_sctype==GWSC_only_G      ) msg = " Self-consistent GW method (G-only) " 
   if (gw_sctype==GWSC_only_W      ) msg = " Self-consistent GW method (W-only) " 
   if (gw_sctype==GWSC_both_G_and_W) msg = " Self-consistent GW method (both G and W) " 
   LAZY_P(msg)

   if (msg == "None") then
     write(msg,'(a,i0)')" Unknown value for gw_sctype= ",gw_sctype
     MSG_ERROR(msg)
   end if

   write(msg,'(a,i3,a,f6.3,a)')" Iteration number: ",istep, "; gw_toldfeig: ",Dtset%gw_toldfeig*Ha_eV," [eV]. "
   LAZY_P(msg)

   if (where_iam == in_screening) then
     LAZY_P(" Input QPS  file     : "//TRIM(Dtfil%fnameabi_qps) )

     LAZY_P(" Output SCR  file    : "//TRIM(Dtfil%fnameabo_scr) )
     LAZY_P(" Output SUSC file    : "//TRIM(Dtfil%fnameabo_sus) )
     LAZY_P(" Output EELF file    : "//TRIM(Dtfil%fnameabo_eelf)    )
     LAZY_P(" Output EM1_LF file  : "//TRIM(Dtfil%fnameabo_em1_lf)  ) 
     LAZY_P(" Output EM1_NLF file : "//TRIM(Dtfil%fnameabo_em1_nlf) )

   else if (where_iam == in_sigma)  then
     LAZY_P(" Input QPS  file  : "//TRIM(Dtfil%fnameabi_qps) )
     LAZY_P(" Input SCR  file  : "//TRIM(Dtfil%fnameabi_scr) )
     LAZY_P(" Input SUSC file  : "//TRIM(Dtfil%fnameabi_sus) )

     LAZY_P(" Output QPS  file : "//TRIM(Dtfil%fnameabo_qps) )
!    LAZY_P(" Output SCR  file : "//TRIM(Dtfil%fnameabo_scr) )
!    LAZY_P(" Output SUSC file : "//TRIM(Dtfil%fnameabo_sus) )

   else 
     write(msg,'(a,i0)')" Unknown value for where_iam= ",where_iam
     MSG_ERROR(msg)
   end if

   LAZY_P(" Prefix for output files: "//TRIM(Dtfil%filnam_ds(4)) )
   LAZY_P(" ")

 end subroutine show_gwfiles

end subroutine gw_driver
!!***
