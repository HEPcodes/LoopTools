* A00.F
* the scalar one-point function
* this file is part of LoopTools
* written by M. Rauch
* last modified 15 Apr 11 th

#include "externals.h"
#include "types.h"
#include "defs.h"


	ComplexType function XA00(m)
	implicit none
	DVAR m

#include "lt.h"

	ComplexType XA0
	external XA0

	XA00 = .25D0*m*XA0(m)
	if( lambda .ge. 0 ) XA00 = XA00 + .125D0*m**2
	end

************************************************************************
* adapter code for C++

	subroutine XA00sub(res, m)
	implicit none
	ComplexType res
	DVAR m

	ComplexType XA00
	external XA00

	res = XA00(m)
	end

