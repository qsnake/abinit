!{\src2tex{textfont=tt}}
!!****p* ABINIT/aim
!! NAME
!! aim
!!
!! FUNCTION
!! Main routine for Bader Atom-In-Molecule analysis.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! WARNING
!! ABINIT rules are not yet followed in the present routine.
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,adini,defad,drvaim,dump_config,herald,inpar,timein
!!      xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program aim

 use defs_basis
 use defs_aimprom
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_build_info

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'aim'
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_63_bader
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
 integer :: fin,ii,ios,iunt,ivst
! Allow for maximum of 100 fc files
 integer,parameter :: natm=100
 integer :: lenstr
 real(dp) :: tcpu,tcpui,twall,twalli
 character(len=fnlen) :: dnfile,fcfile,hname,infile,ofile
 character(len=strlen) :: instr
 character(len=24) :: codename
 type(aim_dataset_type) :: aim_dtset

!******************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

!Initialize MPI
 call xmpi_init()

 call timein(tcpui,twalli)

 unt=21 ! WARNING : this number is used to define other unit numbers, in init.f
 unt0=9
 untout=14

!unto=16
 unto=6  ! XG020629 use standard IO unit => easier testing

 untc=11
 unts=12
 untad=19
 untd=17
 untl=18
 unta=15
 untp=13
 untg=20

 read(*,'(a)') infile
 read(*,'(a)') dnfile
 read(*,'(a)') ofile

 open(unt0,file=infile,status='old',form='formatted',iostat=ivst)
 open(untad,file=dnfile,status='old',form='unformatted',iostat=ivst)

 if (ivst/=0) stop 'err opening input file'
 do ii=1,natm
   iunt=unt+ii
   read(*,'(a)',iostat=ios) fcfile
   if (ios /=0) exit
   open(iunt,file=fcfile,status='old',form='formatted',iostat=ivst)
   if (ivst/=0) stop 'err opening fc-file'
 end do

 fin=len_trim(ofile)
 hname(1:fin)=ofile(1:fin)
 hname(fin+1:fin+1)='.'

 hname(fin+2:fin+4)='out'
 open(untout,file=hname(1:fin+4),status='unknown',form='formatted')
 hname(fin+2:fin+4)='log'

!XG020629 Standard IO is already open
!open(unto,file=hname(1:fin+4),status='unknown',form='formatted')

 codename='AIM   '//repeat(' ',18)
 call herald(codename,abinit_version,untout)
 call herald(codename,abinit_version,unto)
 call dump_config()


!READING OF THE INPUT FILE

!Setting the input variables to their default values
 call defad(aim_dtset)

!Reading of the input file -> one string called instr
 call inpar(instr,lenstr)

!Analysis of the input string, setting of input variables in aim_dtset
 call adini(aim_dtset,instr,lenstr)

!OPENING OF THE OUTPUT FILES

 if (aim_dtset%isurf/=0) hname(fin+2:fin+5)='surf'
 if (aim_dtset%isurf==1) then
   open(unts,file=hname(1:fin+5),status='unknown',form='formatted')
 elseif (aim_dtset%isurf==-1) then
   open(unts,file=hname(1:fin+5),status='old',action='read',form='formatted')
 end if
 if (aim_dtset%crit/=0) hname(fin+2:fin+5)='crit'
 if (aim_dtset%crit>0) then
   open(untc,file=hname(1:fin+5),status='unknown',form='formatted')
 elseif (aim_dtset%crit==-1) then
   open(untc,file=hname(1:fin+5),status='old',action='read',form='formatted')
 end if

 if (aim_dtset%denout==1) then
   hname(fin+2:fin+4)='dn1'
   open(untd,file=hname(1:fin+4),status='unknown',form='formatted')
 elseif (aim_dtset%denout==2) then
   hname(fin+2:fin+4)='dn2'
   open(untd,file=hname(1:fin+4),status='unknown',form='formatted')
 elseif (aim_dtset%denout==3) then
   hname(fin+2:fin+4)='dn3'
   open(untd,file=hname(1:fin+4),status='unknown',form='formatted')
 elseif (aim_dtset%denout==-1) then
   hname(fin+2:fin+4)='dna'
   open(untd,file=hname(1:fin+4),status='unknown',form='unformatted')
 end if

 if (aim_dtset%lapout==1) then
   hname(fin+2:fin+4)='lp1'
   open(untl,file=hname(1:fin+4),status='unknown',form='formatted')
 elseif (aim_dtset%lapout==2) then
   hname(fin+2:fin+4)='lp2'
   open(untl,file=hname(1:fin+4),status='unknown',form='formatted')
 elseif (aim_dtset%lapout==3) then
   hname(fin+2:fin+4)='lp3'
   open(untl,file=hname(1:fin+4),status='unknown',form='formatted')
 elseif (aim_dtset%lapout==-1) then
   hname(fin+2:fin+4)='lpa'
   open(untl,file=hname(1:fin+4),status='unknown',form='unformatted')
 end if
 if (aim_dtset%gpsurf==1) then
   hname(fin+2:fin+3)='gp'
   open(untg,file=hname(1:fin+3),status='unknown',form='formatted')
 end if

 if (aim_dtset%plden==1) then
   hname(fin+2:fin+4)='pld'
   open(untp,file=hname(1:fin+5),status='unknown',form='formatted')
 end if

!MAIN DRIVER OF THE ANALYSIS

 write(std_out,'(a,a,a,a,i4)' )char(10),&
& ' aim : read density file ',trim(dnfile),' from unit number ',untad

 call drvaim(aim_dtset)

!THE TOTAL TIME NEEDED

 call timein(tcpu,twall)

 write(untout,*)
 write(untout,*) "TIME ANALYSIS"
 write(untout,*) "============"
 write(untout,'(/," Time needed (seconds) - total, CP analyse, SURF determination:",/,/,"-         ",3F16.8)') &
& tcpu-tcpui,ttcp,ttsrf
 write(unto,'(/," Time needed (seconds) - total, CP analyse, SURF determination:",/,/,"-         ",3F16.8)')&
& tcpu-tcpui,ttcp,ttsrf
 write(untout,'(a,a,f11.3,a,f11.3,a,a,a,a)') char(10),&
& '+Total cpu time',tcpu-tcpui,&
& '  and wall time',twall-twalli,' sec',char(10),char(10),&
& ' aim : the run completed succesfully.'

!CLOSING OF THE FILES

 close(unt0)
!XG 020629 close(unto)
 close(untout)
 if (aim_dtset%isurf/=0) close(unts)
 if (aim_dtset%crit/=0) close(untc)
 if (aim_dtset%denout/=0) close(untd)
 if (aim_dtset%lapout/=0) close(untl)
 if (aim_dtset%gpsurf/=0) close(untg)
 if (aim_dtset%plden/=0) close(untp)

!Eventual cleaning of MPI run
 call xmpi_end()

 end program aim
!!***
