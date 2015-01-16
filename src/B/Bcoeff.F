* Bcoeff.F
* invoke the two-point tensor coefficients
* this file is part of LoopTools
* last modified 11 Sep 12 th

#include "externals.h"
#include "types.h"

#define npoint 2
#include "defs.h"


	subroutine Bcoeff(res, para)
	implicit none
	ComplexType res(*)
	RealType para(1,*)

#include "lt.h"

	ComplexType cmp(Nbb)
	RealType p, m1, m2
	integer ier(Nbb), ierall, i
	logical dump

	m1 = M(1)
	m2 = M(2)
	p = P(1)

	serial = serial + 1
	dump = ibits(debugkey, DebugB, 1) .ne. 0 .and.
     &    serial .ge. debugfrom .and. serial .le. debugto
	if( dump ) call XDumpPara(2, para, "Bcoeff")

	do i = 1, Nbb
	  res(i) = 0
	enddo

	if( abs(p) + abs(m1) + abs(m2) .lt. eps ) goto 9

	if( lambda .lt. 0 ) then
	  if( lambda .ne. -1 ) goto 9
	  res(bb0) = 1
	  res(bb1) = -.5D0
	  res(bb00) = -(p - 3*(m1 + m2))/12D0
	  res(bb11) = 1/3D0
	  res(bb001) = (p - 2*m1 - 4*m2)/24D0
	  res(bb111) = -.25D0
	  if( m1*m2 .eq. 0 .and. abs(p - m1 - m2) .lt. acc )
     &      res(dbb0) = -.5D0/p
	  if( m2 .eq. 0 .and. abs(p - m1) .lt. acc )
     &      res(dbb1) = .5D0/p
	  res(dbb00) = -1/12D0
	  goto 9
	endif

	goto (1, 2, 3) ibits(versionkey, KeyBget, 2)

	call BcoeffAD(res, para)
	goto 9

1	call BcoeffFF(res, para, ier)
	ierall = 0
	do i = 1, Nbb
	  ierall = max(ierall, ier(i))
	enddo
	if( ierall .gt. warndigits ) then
	  call BcoeffAD(cmp, para)
	  call Bcheck(cmp, res, ier, para)
	endif
	goto 9

2	call BcoeffAD(res, para)
	call BcoeffFF(cmp, para, ier)
	call Bcheck(res, cmp, ier, para)
	goto 9

3	call BcoeffFF(res, para, ier)
	call BcoeffAD(cmp, para)
	call Bcheck(cmp, res, ier, para)

9	if( dump ) call DumpCoeff(2, res)
	end

************************************************************************

	subroutine Bcheck(Ba, Bb, ier, para)
	implicit none
	ComplexType Ba(*), Bb(*)
	integer ier(*)
	DVAR para(1,*)

#include "lt.h"

	integer i
	logical ini

	character*5 name(Nbb)
	data name /"bb0", "bb1", "bb00", "bb11", "bb001", "bb111",
     &    "dbb0", "dbb1", "dbb00", "dbb11"/

	ini = .TRUE.
	do i = 1, Nbb
	  if( abs(Ba(i) - Bb(i)) .gt. maxdev*abs(Ba(i)) ) then
	    if( ini ) then
	      call DumpPara(2, para, "Discrepancy in Bget")
	      ini = .FALSE.
	    endif
	    print *, name(i), " a =", Ba(i)
	    print *, name(i), " b =", Bb(i)
	    if( ier(i) .gt. errdigits ) Bb(i) = Ba(i)
	  endif
	enddo
	end

