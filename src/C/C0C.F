* C0C.F
* the scalar three-point function for complex parameters
* this file is part of LoopTools
* last modified 13 Sep 12 th

#include "externals.h"
#include "types.h"

#define npoint 3
#include "defs.h"


	ComplexType function C0C(p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	ComplexType para(1,Pcc)

	ComplexType C0
	external C0

	P(1) = p1
	P(2) = p2
	P(3) = p1p2
	if( abs(Im(P(1))) + abs(Im(P(2))) + abs(Im(P(3))) .gt. 0 )
     &    print *, "C0C: Complex momenta not implemented"

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
	M(3) = m3
	if( abs(M(3)) .lt. minmass ) M(3) = 0
	if( abs(Im(M(1))) + abs(Im(M(2))) + abs(Im(M(3))) .eq. 0 ) then
	  C0C = C0(p1, p2, p1p2, m1, m2, m3)
          return
        endif

	call C0Cpara(C0C, para)
	end

************************************************************************
* subroutine version for C++

	subroutine c0subc(res, p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType res, p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	ComplexType para(1,Pcc)

	P(1) = p1
	P(2) = p2
	P(3) = p1p2
	if( abs(Im(P(1))) + abs(Im(P(2))) + abs(Im(P(3))) .gt. 0 )
     &    print *, "C0C: Complex momenta not implemented"

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
	M(3) = m3
	if( abs(M(3)) .lt. minmass ) M(3) = 0
	if( abs(Im(M(1))) + abs(Im(M(2))) + abs(Im(M(3))) .eq. 0 ) then
	  call c0sub(res, p1, p2, p1p2, m1, m2, m3)
          return
        endif

	call C0Cpara(res, para)
	end

************************************************************************

	subroutine C0Cpara(res, para)
	implicit none
	ComplexType res, para(1,*)

	integer ier

#include "lt.h"

	ier = 0
	call ffcc0(res, para, ier)

	if( ier .gt. warndigits ) then
	  ier = 0
	  call ffcc0r(res, para, ier)
	  if( ier .gt. warndigits ) then
	    print *, "C0C claims ", ier, " lost digits"
	    print *, "  p1   =", P(1)
	    print *, "  p2   =", P(2)
	    print *, "  p1p2 =", P(3)
	    print *, "  m1   =", M(1)
	    print *, "  m2   =", M(2)
	    print *, "  m3   =", M(3)
	  endif
	endif
	end

