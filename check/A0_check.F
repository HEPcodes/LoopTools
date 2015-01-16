* A0_check.F
* the scalar one-point function
* this function is adapted from Ansgar Denner's bcanew.f
* to the conventions of LoopTools;
* it is used for double-checking the results of FF
* last modified 4 Mar 00 th

#include "defs.h"


	double complex function A0_check(m)
	implicit none
	double precision m

#include "../ff/ff.h"

	if(m .eq. 0) then
	  A0_check = 0
	else
	  A0_check = m*(1 - log(m/mudim2) + divergence)
	endif
	end

