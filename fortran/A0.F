* A0.F
* the scalar one-point function
* this file is part of LoopTools
* last modified 25 Nov 99 th

#include "defs.h"


	double complex function A0(m)
	implicit none
	double precision m

#if CHECK
	external A0_check
#endif

	integer ier

	ier = 0
	call ffxa0(A0, m, ier)
	if(ier .gt. MAX_DIGITS)
     +     print *, "A0 lost ", ier, " digits (m = ",
     +       dsqrt(m), ")"

#if CHECK
	call check_A(A0, A0_check, "A0", m)
#endif
	end

************************************************************************

	subroutine a0sub(result, m)
	implicit none
	double complex result
	double precision m

	double complex A0
	external A0

	result = A0(m)
	end

