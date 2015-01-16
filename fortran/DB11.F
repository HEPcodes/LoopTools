* DB11.F
* the derivative of B11 with respect to p^2
* this file is part of LoopTools
* last modified 25 Nov 99 th

#include "defs.h"


	double complex function DB11(p, m1, m2)
	implicit none
	double precision p, m1, m2

#if CHECK
	external DB11_check
#endif

	integer ier

	ier = 0
	call ffxdb11(DB11, p, m1, m2, ier)

#if CHECK
	call check_B(DB11, DB11_check, "DB11", p, m1, m2)
#endif
	end

************************************************************************

	subroutine db11sub(result, p, m1, m2)
	implicit none
	double complex result
	double precision p, m1, m2

	double complex DB11
	external DB11

	result = DB11(p, m1, m2)
	end

