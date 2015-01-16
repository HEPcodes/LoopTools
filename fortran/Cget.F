* Cget.F
* the three-point tensor coefficients
* this file is part of LoopTools
* last modified 18 Oct 01 th

#include "defs.h"


	integer function Cget(p1, p2, p1p2, m1, m2, m3)
	implicit none
	double precision p1, p2, p1p2, m1, m2, m3

	integer cachelookup
	external cachelookup, Cijk

	double precision para(6)

	double precision Ccache(2)
	common /Cbase/ Ccache

* The following data statement initializes two *pointers* to NULL
* (see cache.c). This is sneaky but ok because IEEE says 0D0 is
* represented as 8 0-bytes.
	data Ccache /0, 0/

	para(1) = p1
	para(2) = p2
	para(3) = p1p2
	para(4) = m1
	para(5) = m2
	para(6) = m3

	Cget = cachelookup(para, Ccache, Cijk, 6, 13)
	end

************************************************************************

	subroutine Cijk(P, C)
	implicit none
	double precision P(6)
	double complex C(13)

	double complex a0i(2), b0p, b1p
	common /bsave/ a0i, b0p, b1p

	double complex B11, C0
	external B11, C0

	XREAL M11, M12, M22, det2, f1, f2
	double complex b1123, b023, b123
	double complex b1113, b113, b1112, b112
	double complex s1, s2

#ifndef LANGUAGE_FORTRAN_90
	double complex lc1, lc2
	lc1() = dcmplx(M11*dble(s1) + M12*dble(s2),
     +    M11*dimag(s1) + M12*dimag(s2))
	lc2() = dcmplx(M12*dble(s1) + M22*dble(s2),
     +    M12*dimag(s1) + M22*dimag(s2))
#endif

	M12 = .5D0*(QEXT(P(2)) - QEXT(P(3)) - P(1))
	det2 = 2*(QEXT(P(1))*P(3) - M12*M12)
	M12 = M12/det2
	M11 = P(3)/det2
	M22 = P(1)/det2
	f1 = QEXT(P(4)) - QEXT(P(5)) + P(1)
	f2 = QEXT(P(4)) - QEXT(P(6)) + P(3)

	C(cc0) = C0(P(1), P(2), P(3), P(4), P(5), P(6))
	b1123 = B11(P(2), P(6), P(5))
	b023 = b0p
	b123 = b1p
	b1113 = B11(P(3), P(4), P(6))
	b113 = b1p
	s1 = b0p - b023 - f1*C(cc0)
	b1112 = B11(P(1), P(4), P(5))
	b112 = b1p
	s2 = b0p - b023 - f2*C(cc0)
	C(cc1) = lc1()
	C(cc2) = lc2()

	C(cc00) = .5D0*(P(4)*C(cc0) +
     +    .5D0*(f1*C(cc1) + f2*C(cc2) + b023)) + .25D0
	s1 = -(f1*C(cc1) + b123) - 2*C(cc00)
	s2 = -(f2*C(cc1) + b123 - b112)
	C(cc11) = lc1()
	C(cc12) = lc2()
	b023 = b023 + b123
	s1 = b023 + b113 - f1*C(cc2)
	s2 = b023 - f2*C(cc2) - 2*C(cc00)
	C(cc22) = lc2()

	C(cc001) = (P(4)*C(cc1) + .5D0*(f1*C(cc11) +
     +    f2*C(cc12) + b123))/3D0 - 1/18D0
	C(cc002) = (P(4)*C(cc2) + .5D0*(f1*C(cc12) +
     +    f2*C(cc22) - b023))/3D0 - 1/18D0
	s1 = -(b1123 + f1*C(cc11)) - 4*C(cc001)
	s2 = -(b1123 + f2*C(cc11) - b1112)
	C(cc111) = lc1()
	C(cc112) = lc2()
	b1123 = b1123 + b023 + b123
	s1 = -(b1123 + f1*C(cc22) - b1113)
	s2 = -(b1123 + f2*C(cc22)) - 4*C(cc002)
	C(cc122) = lc1()
	C(cc222) = lc2()

#ifdef LANGUAGE_FORTRAN_90
	contains

	double complex function lc1
	lc1 = dcmplx(M11*dble(s1) + M12*dble(s2),
     +    M11*dimag(s1) + M12*dimag(s2))
	end function lc1

	double complex function lc2
	lc2 = dcmplx(M12*dble(s1) + M22*dble(s2),
     +    M12*dimag(s1) + M22*dimag(s2))
	end function lc2
#endif
	end

