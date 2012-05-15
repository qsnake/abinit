#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


        subroutine ctrig(n,trig,after,before,now,isign,ic)

 use m_profiling
        use defs_basis
       use defs_fftdata
!       RESTRICTIONS on USAGE
!  Copyright (C) 2002-2007 Stefan Goedecker, CEA Grenoble
!  This file is distributed under the terms of the
!  GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ctrig'
!End of the abilint section

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer after,before,n,isign,ic
        dimension now(mdata),after(mdata),before(mdata),trig(2,n)
!Local variables-------------------------------


        do 111,i=1,ndata
        if (n.eq.ifftdata(1,i)) then
        ic=0
        do 11,j=1,(mdata-1)
        itt=ifftdata(1+j,i)
        if (itt.gt.1) then
        ic=ic+1
        now(j)=ifftdata(1+j,i)
        else
        goto 1000
        end if
11        continue
        goto 1000
        end if
111        continue
        write(std_out,*) 'VALUE OF',n,'NOT ALLOWED FOR FFT, ALLOWED VALUES ARE:'
37        format(15(i5))
        write(std_out,37) (ifftdata(1,i),i=1,ndata)
        stop
1000        continue

        after(1)=1
        before(ic)=1
        do 22,i=2,ic
        after(i)=after(i-1)*now(i-1)
22        before(ic-i+1)=before(ic-i+2)*now(ic-i+2)

!12        format(6(i3))
!        write(std_out,12) (after(i),i=1,ic)
!        write(std_out,12) (now(i),i=1,ic)
!        write(std_out,12) (before(i),i=1,ic)

        angle=isign*two_pi/n
        if (mod(n,2).eq.0) then
        nh=n/2
        trig(1,1)=one
        trig(2,1)=zero
        trig(1,nh+1)=-one
        trig(2,nh+1)=zero
        do 40,i=1,nh-1
        trigc=cos(i*angle)
        trigs=sin(i*angle)
        trig(1,i+1)=trigc
        trig(2,i+1)=trigs
        trig(1,n-i+1)=trigc
        trig(2,n-i+1)=-trigs
40      continue
        else
        nh=(n-1)/2
        trig(1,1)=one
        trig(2,1)=zero
        do 20,i=1,nh
        trigc=cos(i*angle)
        trigs=sin(i*angle)
        trig(1,i+1)=trigc
        trig(2,i+1)=trigs
        trig(1,n-i+1)=trigc
        trig(2,n-i+1)=-trigs
20      continue
        end if


        return
end subroutine ctrig
