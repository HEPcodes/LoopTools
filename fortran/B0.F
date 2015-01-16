* B0.F
* the scalar two-point function
* this file is part of LoopTools
* last modified 25 Nov 99 th

#include "defs.h"


	double complex function B0(p, m1, m2)
	implicit none
	double precision p, m1, m2

#if CHECK
	external B0_check
#endif

	integer ier

	ier = 0
	call ffxb0(B0, p, m1, m2, ier)
	if(ier .gt. MAX_DIGITS)
     +    print *, "B0 lost ", ier, " digits (m1 = ",
     +      dsqrt(m1), ")"

#if CHECK
	call check_B(B0, B0_check, "B0", p, m1, m2)
#endif
	end

************************************************************************

	subroutine b0sub(result, p, m1, m2)
	implicit none
	double complex result
	double precision p, m1, m2

	double complex B0
	external B0

	result = B0(p, m1, m2)
	end

