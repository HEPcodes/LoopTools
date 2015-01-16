* B1.F
* the two-point tensor coefficient of p_mu
* this file is part of LoopTools
* last modified 25 Nov 99 th

#include "defs.h"


	double complex function B1(p, m1, m2)
	implicit none
	double precision p, m1, m2
#include "../ff/ff.h"

	double complex a0i(2), b0p, b1p
	common /bsave/ a0i, b0p, b1p

#if CHECK
	external B1_check
#endif

	integer ier

	ldot = .TRUE.
	ier = 0
	call ffxb0(b0p, p, m1, m2, ier)
	call ffxa0(a0i(1), m1, ier)
	call ffxa0(a0i(2), m2, ier)
	call ffxb1(B1, b0p, a0i, p, m1, m2, fpij2, ier)
	if(ier .gt. MAX_DIGITS)
     +    print *, "B1 lost ", ier, " digits (m1 = ",
     +      dsqrt(m1), ")"
	ldot = .FALSE.

#if CHECK
	call check_B(B1, B1_check, "B1", p, m1, m2)
#endif
	end

************************************************************************

	subroutine b1sub(result, p, m1, m2)
	implicit none
	double complex result
	double precision p, m1, m2

	double complex B1
	external B1

	result = B1(p, m1, m2)
	end

