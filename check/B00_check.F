* B00_check.F
* the two-point tensor coefficient of g_{mu nu}
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools;
* it is used for double-checking the results of FF
* last modified 8 Feb 00 th

#include "defs.h"


	double complex function B00_check(p, m1, m2)
	implicit none
	double precision p, m1, m2

	double complex A0_check, B0_check, B1_check
	external A0_check, B0_check, B1_check

	if(abs(p) .gt. CALACC*(abs(p) + m1 + m2)) then
	  B00_check = ((p + m1 - m2)*B1_check(p, m1, m2) +
     +      2*m1*B0_check(p, m1, m2) + A0_check(m2) +
     +      m1 + m2 - p/3D0)/6D0
	else
	  B00_check = (2*(m2*B0_check(0D0, m1, m2) +
     +      A0_check(m1)) + m1 + m2)/8D0
	endif
	end

