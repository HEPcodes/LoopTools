* DB1.F
* the derivative of B1 with respect to p^2
* this file is part of LoopTools
* last modified 25 Nov 99 th

#include "defs.h"


	double complex function DB1(p, m1, m2)
	implicit none
	double precision p, m1, m2

#if CHECK
	external DB1_check
#endif

	integer ier

	ier = 0
	call ffxdb1(DB1, p, m1, m2, ier)

#if CHECK
	call check_B(DB1, DB1_check, "DB1", p, m1, m2)
#endif
	end

************************************************************************

	subroutine db1sub(result, p, m1, m2)
	implicit none
	double complex result
	double precision p, m1, m2

	double complex DB1
	external DB1

	result = DB1(p, m1, m2)
	end

