* C0.F
* the scalar three-point function
* this file is part of LoopTools
* last modified 25 Nov 99 th

#include "defs.h"


	double complex function C0(p1, p2, p1p2, m1, m2, m3)
	implicit none
	double precision p1, p2, p1p2, m1, m2, m3

#if CHECK
	external C0_check
#endif

	double precision xpi(6)
	integer ier

	xpi(1) = m1
	xpi(2) = m2
	xpi(3) = m3
	xpi(4) = p1
	xpi(5) = p2
	xpi(6) = p1p2
	ier = 0
	call ffxc0(C0, xpi, ier)
	if(ier .gt. MAX_DIGITS) then
	  ier = 0
	  call ffxc0r(C0, xpi, ier)
	  if(ier .gt. MAX_DIGITS)
     +      print *, "C0 lost ", ier, " digits (m1 = ",
     +        dsqrt(m1), ")"
	endif

#if CHECK
	call check_C(C0, C0_check, "C0", p1, p2, p1p2, m1, m2, m3)
#endif
	end

************************************************************************

	subroutine c0sub(result, p1, p2, p1p2, m1, m2, m3)
	implicit none
	double complex result
	double precision p1, p2, p1p2, m1, m2, m3

	double complex C0
	external C0

	result = C0(p1, p2, p1p2, m1, m2, m3)
	end

