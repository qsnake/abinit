!{\src2tex{textfont=tt}}
!!****f* ABINIT/scprqt
!! NAME
!! scprqt
!!
!! FUNCTION
!! Conducts printing inside the scfcv.F90 routine, according to the value of choice.
!! Also checks the convergence with respect to the different criteria.
!! Eventually send a signal to quit the SCF cycle.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG,AF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  choice= if 1 => called at the initialisation of scfcv.f
!!          if 2 => called during the loop in scfcv.f
!!          if 3 => called at the end of scfcv.f
!!  cpus=cpu time limit in seconds
!!  deltae=change in energy between the previous and present SCF cycle
!!  diffor=maximum absolute change in component of fcart between present
!!          and previous SCF cycle.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | chkexit= if non-zero, check whether the user wishes to exit
!!   | enunit=parameter determining units of output energies
!!   | ionmov=governs the movement of atoms (see help file)
!!   | kptopt=option for the generation of k points
!!   | mband=maximum number of bands
!!   | natom=number of atoms in cell.
!!   | nnsclo_now=number of non-self-consistent loops for the current vtrial
!!   |  (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | occopt=option for occupancies
!!   | prtxml=1 if values have to be stored in an XML file.
!!   | prteig=
!!   | prtstm=print STM input variable
!!   | prtvol= control print volume
!!   | usedmatpu=LDA+U: number of SCF steps keeping occ. matrix fixed
!!   | usepawu=0 if no LDA+U; 1 if LDA+U
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  etotal=total energy (hartree)
!!  favg(3)=average of forces (ha/bohr)
!!  fcart(3,natom)=cartesian forces (hartree/bohr)
!!  fermie=fermi energy (Hartree)
!!  fname_eig=filename for printing of the eigenenergies
!!  character(len=fnlen) :: filnam1=character strings giving input file name
!!  initGS= 1 if one GS SCF cycle has already be done
!!  iscf=( <= 0 =>non-SCF), >0 => SCF)
!!   iscf =1 => determination of the largest eigenvalue of the SCF cycle
!!   iscf =2 => SCF cycle, simple mixing
!!   iscf =3 => SCF cycle, anderson mixing
!!   iscf =5 => SCF cycle, CG based on estimations of gradients of the energy
!!   iscf =6 => SCF cycle, CG based on true minimization of the energy
!!   iscf =-3, although non-SCF, the energy is computed, so print it here.
!!  istep=number of the SCF iteration (needed if choice=2)
!!  kpt(3,nkpt)=reduced coordinates of k points.
!!  maxfor=maximum absolute value of fcart
!!  moved_atm_inside: if==1, the atoms are allowed to move.
!!  mpi_enreg=informations about MPI parallelization
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  nkpt=number of k points
!!  nstep=number of steps expected in iterations.
!!  occ(mband*nkpt*nsppol)=occupation number for each band at each k point.
!!  optres=0 if the residual (res2) is a POTENTIAL residual
!!         1 if the residual (res2) is a DENSITY residual
!!  prtfor=1 only if forces have to be printed (0 otherwise)
!!  prtxml=1 if XML file has to be output
!!  res2=square of the density/potential residual
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!         in Wavelets mode, it is used as the maximum value for the gradient norm.
!!  response= if 0, GS case, if 1, RF case.
!!  tollist(12)=tolerance list. Presently, the following are defined :
!!    tollist(1)=tolmxf ; tollist(2)=tolwfr ; tollist(3)=toldff
!!    tollist(4)=toldfe ; tollist(5)=toleig ; tollist(6)=tolvrs
!!    tollist(7)=tolrff
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vxcavg=mean of the vxc potential
!!  wtk(nkpt)=weight assigned to each k point.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  quit= 0 if the SCF cycle is not finished ; 1 otherwise.
!!
!! PARENTS
!!      afterscfloop,scfcv,scfcv3
!!
!! CHILDREN
!!      chkexi,leave_new,prteigrs,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scprqt(choice,cpus,deltae,diffor,dtset,&
&  eigen,etotal,favg,fcart,fermie,fname_eig,filnam1,initGS,&
&  iscf,istep,kpt,maxfor,moved_atm_inside,mpi_enreg,&
&  nband,nkpt,nstep,occ,optres,&
&  prtfor,prtxml,quit,res2,resid,residm,response,tollist,usepaw,&
&  vxcavg,wtk,xred,&
&  electronpositron) ! optional argument)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_electronpositron, only : electronpositron_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scprqt'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_59_io_mpi
 use interfaces_67_common, except_this_one => scprqt
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,initGS,iscf,istep,moved_atm_inside,nkpt,nstep
 integer,intent(in) :: optres,prtfor,prtxml,response,usepaw
 integer,intent(out) :: quit
 real(dp),intent(in) :: cpus,deltae,diffor,etotal,fermie,maxfor,res2,residm
 real(dp),intent(in) :: vxcavg
 character(len=fnlen),intent(in) :: fname_eig,filnam1
 type(electronpositron_type),pointer,optional :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: nband(nkpt*dtset%nsppol)
 real(dp),intent(in) :: eigen(dtset%mband*nkpt*dtset%nsppol),favg(3)
 real(dp),intent(in) :: fcart(3,dtset%natom),kpt(3,nkpt)
 real(dp),intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
 real(dp),intent(in) :: resid(dtset%mband*nkpt*dtset%nsppol),tollist(12)
 real(dp),intent(in) :: wtk(nkpt),xred(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer,save :: toldfe_ok,toldff_ok,tolrff_ok,ttoldfe,ttoldff,ttolrff,ttolvrs
 integer,save :: ttolwfr
 integer :: iatom,iband,iexit,ikpt,isppol,nband_index,nband_k,openexit,option
 integer :: tmagnet
 real(dp),save :: toldfe,toldff,tolrff,tolvrs,tolwfr
 real(dp) :: diff_e,diff_f,magnet,rhodn,rhoup
 logical :: noquit
 character(len=500) :: message
 character(len=8) :: colname
!arrays
 real(dp) :: f_tmp(3)

! *********************************************************************

!DEBUG
!write(std_out,*)' scprqt : enter, iscf,nstep= ',iscf,nstep
!stop
!ENDDEBUG

 quit=0

 tmagnet=0
 if(response==0.and.(iscf>0.or.iscf==-3).and.dtset%nsppol==2.and.dtset%occopt>2)tmagnet=1

 if(choice==1)then

!  Examine tolerance criteria
   tolwfr=tollist(2)
   toldff=tollist(3)
   toldfe=tollist(4)
   tolvrs=tollist(6)
   tolrff=tollist(7)
   ttolwfr=0 ; ttoldff=0 ; ttoldfe=0 ; ttolvrs=0; ttolrff=0;
   if(abs(tolwfr)>tiny(0.0_dp))ttolwfr=1
   if(abs(toldff)>tiny(0.0_dp))ttoldff=1
   if(abs(tolrff)>tiny(0.0_dp))ttolrff=1
   if(abs(toldfe)>tiny(0.0_dp))ttoldfe=1
   if(abs(tolvrs)>tiny(0.0_dp))ttolvrs=1
!  If non-scf calculations, tolwfr must be defined
!  FIXME: MJV 26/3/2010: this should be just after initialization, not here.
   if(ttolwfr /= 1 .and. (iscf<=0 .and. iscf/=-3) )then
     write(message, '(a,a,a,a,a,a,es14.6,a,a)' ) ch10,&
&     ' scprqt : ERROR -',ch10,&
&     '   when iscf <=0 and /= -3, tolwfr must be strictly',ch10,&
&     '   positive, while it is ',tolwfr,ch10,&
&     '   Action : change tolwfr in your input file and resubmit the job.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
!  toldff only allowed when prtfor==1
!  FIXME: this test should be done on input, not during calculation
   if((ttoldff == 1 .or. ttolrff == 1) .and. prtfor==0 )then
     write(message, '(a,a,a,a)' ) ch10,&
&     ' scprqt : BUG -',ch10,&
&     '   toldff only allowed when prtfor=1 !'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
!  If SCF calculations, one and only one of these can differ from zero
!  FIXME: this test should be done on input, not during calculation
   if(ttolwfr+ttoldff+ttoldfe+ttolvrs+ttolrff /= 1 .and. (iscf>0 .or. iscf==-3))then
     write(message, &
&     '(9a,es14.6,a,es14.6,a,es14.6,a,es14.6,a,a,es14.6,a,a,a)' ) ch10,&
&     ' scprqt : ERROR -',ch10,&
&     '  For the SCF case, one and only one of the input tolerance criteria ',ch10,&
&     '  tolwfr, toldff, tolrff, toldfe or tolvrs ',&
&     'must differ from zero, while they are',ch10,&
&     '  tolwfr=',tolwfr,', toldff=',toldff,', tolrff=',tolrff,', toldfe=',toldfe,ch10,&
&     '  and tolvrs=',tolvrs,' .',ch10,&
&     '  Action : change your input file and resubmit the job.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   if (dtset%usewvl == 1) then
     write(colname, "(A)") "grdnorm "
   else
     write(colname, "(A)") "residm  "
   end if
   if (nstep>0 .and. (iscf>0 .or.iscf==-3) .and. dtset%prtstm==0) then
     if(tmagnet==1)then
       if (prtfor==0) then
         if (optres==0) then
           write(message, '(a,a,a,a)' ) ch10,&
&           '     iter   Etot(hartree)     deltaE(h) ',colname,&
&           '  vres2    magn'
         else
           write(message, '(a,a,a,a)' ) ch10,&
&           '     iter   Etot(hartree)     deltaE(h) ',colname,&
&           '  nres2    magn'
         end if
       else
         if (optres==0) then
           write(message, '(a,a,a,a)' ) ch10,&
&           '     iter   Etot(hartree)     deltaE(h) ',colname,&
&           '  vres2   diffor   maxfor   magn'
         else
           write(message, '(a,a,a,a)' ) ch10,&
&           '     iter   Etot(hartree)     deltaE(h) ',colname,&
&           '  nres2   diffor   maxfor   magn'
         end if
       end if
     else
       if(response==0)then
         if (prtfor==0) then
           if (optres==0) then
             write(message, '(a,a,a,a)' ) ch10,&
&             '     iter   Etot(hartree)      deltaE(h)  ', colname, '   vres2'
           else
             write(message, '(a,a,a,a)' ) ch10,&
&             '     iter   Etot(hartree)      deltaE(h)  ', colname, '   nres2'
           end if
         else
           if (optres==0) then
             write(message, '(a,a,a,a)' ) ch10,&
&             '     iter   Etot(hartree)      deltaE(h)  ',colname,&
&             '   vres2    diffor    maxfor '
           else
             write(message, '(a,a,a,a)' ) ch10,&
&             '     iter   Etot(hartree)      deltaE(h)  ',colname,&
&             '   nres2    diffor    maxfor '
           end if
         end if
       else
         if (optres==0) then
           write(message, '(a,a,a,a)' ) ch10,&
&           '     iter   2DEtotal(Ha)        deltaE(Ha) ', colname, '  vres2'
         else
           write(message, '(a,a,a,a)' ) ch10,&
&           '     iter   2DEtotal(Ha)        deltaE(Ha) ', colname, '  nres2'
         end if
       end if
     end if
     call wrtout(ab_out,message,'COLL')
!    call wrtout(std_out,message,'COLL')
   end if

 else if(choice==2)then

!  Conduct printing

!  If extra output follows, then put a blank line into the output here
   if (dtset%prtvol>=10) then
     write(message, '(a)' )' '
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if

!  Calculate up and down charge and magnetization
   if(tmagnet==1) then
     rhoup = zero
     rhodn = zero
     nband_index = 1
     do isppol=1,dtset%nsppol
       do ikpt=1,nkpt
         nband_k=nband(ikpt+(isppol-1)*nkpt)
         do iband=1,nband_k
           if(isppol==1) rhoup = rhoup + wtk(ikpt)*occ(nband_index)
           if(isppol==2) rhodn = rhodn + wtk(ikpt)*occ(nband_index)
           nband_index = nband_index + 1
         end do
       end do
     end do
     magnet = abs(rhoup - rhodn)
   end if

   if (prtxml == 1) then
     write(ab_xml_out, "(A)", advance = "NO") '      <scfcvStep'
     write(message, "(es22.10)") etotal
     write(ab_xml_out, "(A,A,A)", advance = "NO") ' eTotal="', trim(message) ,'"'
     write(message, "(es20.8)") deltae
     write(ab_xml_out, "(A,A,A)", advance = "NO") ' deltaETotal="', trim(message) ,'"'
     write(message, "(es20.8)") residm
     write(ab_xml_out, "(A,A,A)", advance = "NO") ' maxResid="', trim(message) ,'"'
     write(message, "(es20.8)") res2
     if (optres == 0) then
       write(ab_xml_out, "(A,A,A)", advance = "NO") ' potResid="', trim(message) ,'"'
     else
       write(ab_xml_out, "(A,A,A)", advance = "NO") ' denResid="', trim(message) ,'"'
     end if
     if (tmagnet == 1) then
       write(message, "(es20.8)") magnet
       write(ab_xml_out, "(A,A,A)", advance = "NO") ' magn="', trim(message) ,'"'
     end if
     if (prtfor == 1) then
       write(message, "(es20.8)") diffor
       write(ab_xml_out, "(A,A,A)", advance = "NO") ' deltaForces="', trim(message) ,'"'
       write(message, "(es20.8)") maxfor
       write(ab_xml_out, "(A,A,A)", advance = "NO") ' maxForces="', trim(message) ,'"'
     end if
     write(ab_xml_out, "(A)") " />"
   end if

!  Print total (free) energy (hartree) and other convergence measures
   if(dtset%prtstm==0)then
     if(tmagnet==1)then
       if (prtfor==0) then
         if(magnet<10)then
           write(message, '(a,i3,1p,g22.14,3es9.2,0p,f6.3)' ) &
&           ' ETOT',istep,etotal,deltae,residm,res2,magnet
         else if(magnet<100)then
           write(message, '(a,i3,1p,g22.14,3es9.2,0p,f6.2)' ) &
&           ' ETOT',istep,etotal,deltae,residm,res2,magnet
         else
           write(message, '(a,i3,1p,g22.14,3es9.2,0p,f6.1)' ) &
&           ' ETOT',istep,etotal,deltae,residm,res2,magnet
         end if
       else
         if(magnet<10)then
           write(message, '(a,i3,1p,g22.14,3es9.2,es8.1,es9.2,0p,f6.3)' ) &
&           ' ETOT',istep,etotal,deltae,residm,res2,diffor,maxfor,magnet
         else if(magnet<100)then
           write(message, '(a,i3,1p,g22.14,3es9.2,es8.1,es9.2,0p,f6.2)' ) &
&           ' ETOT',istep,etotal,deltae,residm,res2,diffor,maxfor,magnet
         else
           write(message, '(a,i3,1p,g22.14,3es9.2,es8.1,es9.2,0p,f6.1)' ) &
&           ' ETOT',istep,etotal,deltae,residm,res2,diffor,maxfor,magnet
         end if
       end if
     else
       if(response==0)then
         if (prtfor==0) then
           write(message, '(a,i3,1p,g22.14,3es10.3)' ) &
&           ' ETOT',istep,etotal,deltae,residm,res2
         else
           write(message, '(a,i3,1p,g22.14,5es10.3)' ) &
&           ' ETOT',istep,etotal,deltae,residm,res2,diffor,maxfor
         end if
       else
         if(istep==1)then
           write(message, '(a,i3,1p,g22.14,1x,es10.3,es10.3,es10.3)' ) &
&           '-ETOT',istep,etotal,deltae,residm,res2
         else
           write(message, '(a,i3,1p,g22.14,1x,es10.3,es10.3,es10.3)' ) &
&           ' ETOT',istep,etotal,deltae,residm,res2
         end if
       end if
     end if
     call wrtout(ab_out,message,'COLL')

!    BEGIN TF_CHANGES
     if(mpi_enreg%paral_compil_respfn==1 .and. mpi_enreg%me_respfn==0) then
       call wrtout(std_out,  message,'PERS')
     elseif(mpi_enreg%paral_compil_respfn==0) then
       call wrtout(std_out,  message,'COLL')
     end if
!    END TF_CHANGES

   end if ! dtset%prtstm==0

!  Print positions/forces every step if dtset%prtvol>=10 and iscf>0 or -3 and GS case
   if (dtset%prtvol>=10.and.(iscf>0.or.iscf==-3).and.response==0.and.dtset%prtstm==0) then

     write(message, '(a)' ) ' '
     call wrtout(ab_out,message,'COLL')

!    Print up and down charge and magnetization
     if(tmagnet==1) then
       write(message,'(a,f11.6,a,f11.6,a,f10.6)')&
&       ' #electrons spin up=',rhoup,&
&       ', spin down=',rhodn,', magnetisation=',magnet
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if

!    Moreover, print atomic positions if dtset%ionmov==4, and moved_atm_inside==1
     if (dtset%ionmov==4 .and. moved_atm_inside==1)then
       write(message, '(a)' ) ' reduced coordinates :'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
       do iatom=1,dtset%natom
         write(message, '(i5,1x,3es21.11)' ) iatom,xred(:,iatom)
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end do
     end if

!    Slightly change favg for printing reasons
     if (prtfor>0) then
       f_tmp(:)=favg(:)
       if(abs(favg(1))<1.0d-13)f_tmp(1)=zero
       if(abs(favg(2))<1.0d-13)f_tmp(2)=zero
       if(abs(favg(3))<1.0d-13)f_tmp(3)=zero
       write(message, '(a,3es10.2)' ) &
&       ' cartesian forces (ha/bohr); non-corrected avg=',f_tmp(:)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
       do iatom=1,dtset%natom
         f_tmp(:)=fcart(:,iatom)
         if(abs(fcart(1,iatom))<1.0d-13)f_tmp(1)=zero
         if(abs(fcart(2,iatom))<1.0d-13)f_tmp(2)=zero
         if(abs(fcart(3,iatom))<1.0d-13)f_tmp(3)=zero
         write(message, '(i5,1x,3es21.11)' ) iatom,f_tmp(:)
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end do
     end if

   end if

!  Print eigenvalues every step if dtset%prtvol>=10 and GS case
   if (dtset%prtvol>=10 .and. response==0 .and. dtset%tfkinfunc==0) then
     option=1
     call prteigrs(eigen,dtset%enunit,fermie,fname_eig,ab_out,iscf,kpt,dtset%kptopt,dtset%mband,&
&     nband,nkpt,dtset%nnsclo,dtset%nsppol,occ,dtset%occopt,option,dtset%prteig,dtset%prtvol,resid,tolwfr,vxcavg,wtk)
     call prteigrs(eigen,dtset%enunit,fermie,fname_eig,  06,iscf,kpt,dtset%kptopt,dtset%mband,&
&     nband,nkpt,dtset%nnsclo,dtset%nsppol,occ,dtset%occopt,option,dtset%prteig,dtset%prtvol,resid,tolwfr,vxcavg,wtk)
   end if

   if(response==0)then
     write(message, '(a,1p,e15.7,a)'  ) ' scprqt: <Vxc>=',vxcavg,' hartree'
     call wrtout(std_out,message,'COLL')
   end if

!  Check whether exiting was required by the user.
   openexit=1 ; if(dtset%chkexit==0) openexit=0
   call chkexi(cpus,filnam1,iexit,ab_out,mpi_enreg,openexit)
   if (iexit/=0) quit=1

!  In special cases, do not quit even if convergence is reached
   noquit=((istep<nstep).and.(usepaw==1).and.(dtset%usepawu>0).and.&
&   (dtset%usedmatpu/=0).and.(istep<=abs(dtset%usedmatpu)).and.&
&   (dtset%usedmatpu<0.or.initGS==0))

!  Additional stuff for electron/positron
   if (present(electronpositron)) then
     if (associated(electronpositron)) then
       if (electronpositron%istep_scf==1) then
         toldff_ok=0;tolrff_ok=0;toldfe_ok=0
       end if
     end if
   end if

!  Stopping criteria in the SCF case
   if(iscf>1 .or. iscf==-3) then
!    Here treat the tolwfr criterion : if maximum residual is less than
!    input tolwfr, stop steps (exit loop here)
     if( ttolwfr==1 .and. residm<tolwfr .and. (.not.noquit)) then
       if (dtset%usewvl == 0) then
         write(message, '(a,a,i5,a,1p,e10.2,a,e10.2,a,a)' )ch10, &
&         ' At SCF step',istep,'   max residual=',residm,&
&         ' < tolwfr=',tolwfr,' =>converged.'
       else
         write(message, '(a,a,i5,a,1p,e10.2,a,e10.2,a,a)' )ch10, &
&         ' At SCF step',istep,'   max grdnorm=',residm,&
&         ' < tolwfr=',tolwfr,' =>converged.'
       end if
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       quit=1
     end if
!    Here treat the toldff criterion : if maximum change of fcart is less than
!    input toldff twice consecutively, stop steps (exit loop here)
     if( ttoldff==1 ) then
       if( istep==1 )then
         toldff_ok=0
       else if (diffor<toldff) then
         toldff_ok=toldff_ok+1
       else
         toldff_ok=0
       end if
       if(toldff_ok==2 .and. (.not.noquit))then
         write(message, '(a,a,i5,a,a,a,es11.3,a,es11.3)' ) ch10, &
&         ' At SCF step',istep,', forces are converged : ',ch10,&
&         '  for the second time, max diff in force=',diffor,' < toldff=',toldff
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         quit=1
       end if
     end if
!    Here treat the tolrff criterion : if maximum change of fcart is less than
!    input tolrff times fcart itself twice consecutively, stop steps (exit loop here)
     if( ttolrff==1 ) then
       if( istep==1 )then
         tolrff_ok=0
!        27/7/2009: added test for absolute value of maxfor, otherwise if it is 0 this never exits the scf loop.
       else if (diffor<tolrff*maxfor .or. (maxfor < tol16 .and. diffor < tol16)) then
         tolrff_ok=tolrff_ok+1
       else
         tolrff_ok=0
       end if
       if(tolrff_ok==2 .and. (.not.noquit))then
         write(message, '(a,a,i5,a,a,a,es11.3,a,es11.3,a)' ) ch10, &
&         ' At SCF step',istep,', forces are sufficiently converged : ',ch10,&
&         '  for the second time, max diff in force=',diffor,&
&         ' is less than < tolrff=',tolrff, ' times max force'
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         quit=1
       end if
     end if
!    Here treat the toldfe criterion : if the change of energy is less than
!    input toldfe twice consecutively, stop steps (exit loop here)
     if( ttoldfe==1 ) then
       if( istep==1 )then
         toldfe_ok=0
       else if (abs(deltae)<toldfe) then
         toldfe_ok=toldfe_ok+1
       else
         toldfe_ok=0
       end if
       if(toldfe_ok==2 .and. (.not.noquit))then
         write(message, '(a,a,i5,a,a,a,es11.3,a,es11.3)' ) ch10, &
&         ' At SCF step',istep,', etot is converged : ',ch10,&
&         '  for the second time, diff in etot=',abs(deltae),' < toldfe=',toldfe
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')
         quit=1
       end if
     end if
!    Here treat the tolvrs criterion : if density/potential residual (squared)
!    is less than input tolvrs, stop steps (exit loop here)
     if( ttolvrs==1 .and. res2<tolvrs .and. (.not.noquit)) then
       if (optres==0) then
         write(message, '(a,a,i5,a,1p,e10.2,a,e10.2,a)' ) ch10,&
&         ' At SCF step',istep,'       vres2   =',res2,&
&         ' < tolvrs=',tolvrs,' =>converged.'
       else
         write(message, '(a,a,i5,a,1p,e10.2,a,e10.2,a)' ) ch10,&
&         ' At SCF step',istep,'       nres2   =',res2,&
&         ' < tolvrs=',tolvrs,' =>converged.'
       end if
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       quit=1
     end if

     if (quit==1.and.noquit) then
       write(message, '(a,a,a)' ) ch10, &
&       ' SCF cycle will continue as it is in an initialization stage',&
&       ' (occ. matrix was kept constant)...'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if

   end if


 else if(choice==3)then

!  If wavefunction convergence was not reached (for nstep>0),
!  print a warning
   if(nstep>0)then
     if ( (ttolwfr==1 .and. residm > tolwfr) .or.  &
&     (ttoldff==1 .and. diffor > toldff) .or.  &
&     (ttolrff==1 .and. diffor > tolrff*maxfor) .or.  &
&     (ttoldfe==1 .and. abs(deltae) > toldfe) .or.  &
&     (ttolvrs==1 .and. res2  > tolvrs)        )     then
       if(iscf>=1 .or. iscf==-3)then
         write(message, '(a,a,a,a,i5,a)' ) ch10,&
&         ' scprqt:  WARNING -',ch10,&
&         '  nstep=',nstep,' was not enough SCF cycles to converge;'
       else
         write(message, '(a,a,a,a,i5,a)' ) ch10,&
&         ' scprqt:  WARNING -',ch10,&
&         '  nstep=',nstep,' was not enough non-SCF iterations to converge;'
       end if
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       if(ttolwfr==1)then
         if (dtset%usewvl == 0) then
           write(message, '(a,es11.3,a,es11.3,a)' ) &
&           '  maximum residual=',residm,' exceeds tolwfr=',tolwfr,ch10
         else
           write(message, '(a,es11.3,a,es11.3,a)' ) &
&           '  maximum grdnorm=',residm,' exceeds tolwfr=',tolwfr,ch10
         end if
       else if(ttoldff==1)then
         write(message, '(a,es11.3,a,es11.3,a)' ) &
&         '  maximum force difference=',diffor,' exceeds toldff=',toldff,ch10
       else if(ttolrff==1)then
         write(message, '(a,es11.3,a,es11.3,a)' ) &
&         '  maximum force difference=',diffor,' exceeds tolrff*maxfor=',&
&         tolrff*maxfor,ch10
       else if(ttoldfe==1)then
         write(message, '(a,es11.3,a,es11.3,a)' ) &
&         '  maximum energy difference=',abs(deltae),&
&         ' exceeds toldfe=',toldfe,ch10
       else if(ttolvrs==1)then
         if (optres==0) then
           write(message, '(a,es11.3,a,es11.3,a)' ) &
&           '  potential residual=',res2,' exceeds tolvrs=',tolvrs,ch10
         else
           write(message, '(a,es11.3,a,es11.3,a)' ) &
&           '  density residual=',res2,' exceeds tolvrs=',tolvrs,ch10
         end if
       end if
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       if (prtxml == 1) then
         write(ab_xml_out, "(A)", advance = "NO") '      <status cvState="Failed"'
       end if
     else    ! Convergence is OK
       if (prtxml == 1) then
         write(ab_xml_out, "(A)", advance = "NO") '      <status cvState="Ok"'
       end if
     end if ! test for convergence reached or not

     if (prtxml == 1) then
       if (ttoldfe == 1) then
         write(ab_xml_out, "(A)") ' stop-criterion="toldfe" />'
       else if (ttoldff == 1) then
         write(ab_xml_out, "(A)") ' stop-criterion="toldff" />'
       else if (ttolrff == 1) then
         write(ab_xml_out, "(A)") ' stop-criterion="tolrff" />'
       else if (ttolvrs == 1) then
         write(ab_xml_out, "(A)") ' stop-criterion="tolvrs" />'
       else if (ttolwfr == 1) then
         write(ab_xml_out, "(A)") ' stop-criterion="tolwfr" />'
       else
         write(ab_xml_out, "(A)") ' />'
       end if
     end if

   end if ! nstep == 0 : no output

 else

   write(message, '(a,a,a,a,i3,a)' ) ch10,&
&   ' scprqt : BUG -',ch10,&
&   '  choice=',choice,' is not an allowed value.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')

 end if

!Additional stuff for the two-component DFT SCF cycle (electrons+positron)
 if (present(electronpositron)) then
   if (associated(electronpositron)) then
     electronpositron%scf_converged=(quit==1)
     if (dtset%positron<0) then
       diff_e=abs(etotal-electronpositron%etotal_prev)
       diff_f=abs(maxfor-electronpositron%maxfor_prev)
     end if
     if (choice==1) then
       ttoldff=0;ttoldfe=0
       if(abs(dtset%postoldff)>tiny(0.0_dp))ttoldff=1
       if(abs(dtset%postoldfe)>tiny(0.0_dp))ttoldfe=1
       if (dtset%positron<0.and.ttoldff+ttoldfe/=1.and.iscf>0) then
         write(message, '(4a)' ) ch10,&
&         ' scprqt : ERROR -',ch10,&
&         '  one and only one of toldff or toldfe must differ from zero !'
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')
       end if
     end if
     if (choice==2) then
       if (dtset%positron<0.and.istep<=nstep) then
         if (electronpositron%scf_converged) then
           if (electronpositron%istep==electronpositron%nstep) then
             quit=1
           else if ((.not.noquit).and.&
&             (diff_e<electronpositron%postoldfe.or.diff_f<electronpositron%postoldff).and.&
&             (mod(electronpositron%calctype,2)==0.or.(dtset%positron>-20.and.dtset%positron/=-2))) then
             if (diff_e<electronpositron%postoldfe) then
               write(message, '(2a,i5,5a,es11.3,a,es11.3)' ) ch10, &
&               ' At SCF step',istep,', the difference between',ch10,&
&               ' etotal from electronic calculation and etotal from positronic calculation',ch10,&
&               ' is converged :  diff(etot_el-etot_pos)=',diff_e,' < postoldfe=',electronpositron%postoldfe
             else
               write(message, '(2a,i5,5a,es11.3,a,es11.3)' ) ch10, &
&               ' At SCF step',istep,', the difference between',ch10,&
&               ' max. force from electronic calculation and max. force from positronic calculation',ch10,&
&               ' is converged :  diff(maxfor_el-maxfor_pos)=',diff_f,' < postoldff=',electronpositron%postoldff
             end if
             call wrtout(ab_out,message,'COLL')
             call wrtout(std_out,message,'COLL')
           else
             quit=0
           end if
         end if
       end if
     end if
     if (choice==3) then
       if (dtset%positron<0.and.nstep>0)then
         if (diff_e>=electronpositron%postoldfe.and.abs(dtset%postoldfe)>tiny(0.0_dp)) then
           write(message, '(4a,i5,5a,es11.3,a,es11.3)' ) ch10,&
&           ' scprqt:  WARNING -',ch10,&
&           '  posnstep=',dtset%posnstep,' was not enough SCF cycles to converge difference between',ch10,&
&           '  etotal from electronic calculation and etotal from positronic calculation;',ch10,&
&           '  diff=',diff_e,' exceeds postoldfe=',electronpositron%postoldfe
           call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')
         end if
         if (diff_f>=electronpositron%postoldff.and.abs(dtset%postoldff)>tiny(0.0_dp)) then
           write(message, '(4a,i5,5a,es11.3,a,es11.3)' ) ch10,&
&           ' scprqt:  WARNING -',ch10,&
&           '  posnstep=',dtset%posnstep,' was not enough SCF cycles to converge difference between',ch10,&
&           '  max. force from electronic calculation and max. force from positronic calculation;',ch10,&
&           '  diff=',diff_e,' exceeds postoldff=',electronpositron%postoldff
           call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')
         end if
       end if
     end if
   end if
 end if

end subroutine scprqt
!!***
