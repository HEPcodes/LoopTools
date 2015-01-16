* DB00_check.F
* the derivative of B00 with respect to p^2
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools;
* it is used for double-checking the results of FF
* last modified 8 Feb 00 th

#include "defs.h"


	double complex function DB00_check(p, m1, m2)
	implicit none
	double precision p, m1, m2

	double complex B1_check, DB0_check, DB1_check
	external B1_check, DB0_check, DB1_check

	DB00_check = 1/6D0*(2*m1*DB0_check(p, m1, m2) +
     +    B1_check(p, m1, m2) +
     +    (p + m1 - m2)*DB1_check(p, m1, m2)) - 1/18D0
	end

