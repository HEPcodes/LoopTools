* B11.F
* the two-point tensor coefficient of p_mu p_nu
* this file is part of LoopTools
* last modified 25 Nov 99 th

#include "defs.h"


	double complex function B11(p, m1, m2)
	implicit none
	double precision p, m1, m2
#include "../ff/ff.h"

	double complex a0i(2), b0p, b1p
	common /bsave/ a0i, b0p, b1p

	double complex B1
	external B1

#if CHECK
	external B11_check
#endif

	double complex b2(2)
	integer ier

	b1p = B1(p, m1, m2)
	ier = 0
	call ffxb2p(b2, b1p, b0p, a0i, p, m1, m2, fpij2, ier)
	if(ier .gt. MAX_DIGITS)
     +    print *, "B11 lost ", ier, " digits (m1 = ",
     +      dsqrt(m1), ")"
	B11 = b2(1)

#if CHECK
	call check_B(B11, B11_check, "B11", p, m1, m2)
#endif
	end

************************************************************************

	subroutine b11sub(result, p, m1, m2)
	implicit none
	double complex result
	double precision p, m1, m2

	double complex B11
	external B11

	result = B11(p, m1, m2)
	end

