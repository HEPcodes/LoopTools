*	tools.F
*		this file is part of LoopTools;
*		it provides the scalar and tensor
*		one-loop functions in Fortran
*		last modified 19 Jul 99 th


#ifndef __TOOLS_F__
#define __TOOLS_F__

* the following user #defines (defined _before_ including tools.F)
* will be recognized:
*   MUDIM         -- to check UV finiteness (independence of MUDIM)
*   DELTA         -- to check UV finiteness (independence of DELTA)
*   LAMBDA        -- to check IR finiteness (independence of LAMBDA)
*   FFCHECK       -- to double-check the FF results
*   MAX_DEVIATION -- the maximum deviation allowed when using FFCHECK,
*                    i.e. |ffresult - check|/|ffresult| < MAX_DEVIATION

#ifndef MUDIM
#define MUDIM 1D0
#endif

#ifndef LAMBDA
#define LAMBDA 1D0
#endif

#ifndef DELTA
#define DELTA 0D0
#endif

#ifndef MAX_DEVIATION
#define MAX_DEVIATION 1D-12
#endif

#define MAX_DIGITS 9
#define MINMASS 1D-10

#include "defs.h"

#define IN_TOOLS_F

************************************************************************
	double complex function A0(m)
	implicit none
	double precision m

	double precision mudim, divergence
	common /cutoff/ mudim, divergence

#ifdef FFCHECK
	external A0_check
#endif

	integer ier

	ier = 0
	call ffxa0(A0, divergence, mudim, m, ier)
	if(ier .gt. MAX_DIGITS)
     +     print *, "A0 lost ", ier, " digits (m = ",
     +       dsqrt(m), ")"

#ifdef FFCHECK
	call check_A(A0, A0_check, "A0", m)
#endif
	end

************************************************************************
	double complex function B0(p, m1, m2)
	implicit none
	double precision p, m1, m2

#ifdef FFCHECK
	external B0_check
#endif

	double precision mudim, divergence
	common /cutoff/ mudim, divergence

	integer ier

	ier = 0
	call ffxb0(B0, divergence, mudim, p, m1, m2, ier)
	if(ier .gt. MAX_DIGITS)
     +    print *, "B0 lost ", ier, " digits (m1 = ",
     +      dsqrt(m1), ")"

#ifdef FFCHECK
	call check_B(B0, B0_check, "B0", p, m1, m2)
#endif
	end

************************************************************************
	double complex function DB0(p, m1, m2)
	implicit none
	double precision p, m1, m2

	double complex db0p, pdb0p
	common /deriv/ db0p, pdb0p

#ifdef FFCHECK
	external DB0_check
#endif

	double precision dm
	integer ier

	dm = (sqrt(m1) - sqrt(m2))**2
	if(p .eq. dm .and. p .ne. 0D0 .and.
     +      m1 .ne. 0D0 .and. m2 .ne. 0D0) then
	  pdb0p = (m2 - m1)/2D0/dm*log(m2/m1) - 2D0
	  db0p = pdb0p/dm
	else
	  ier = 0
	  call ffxdb0(db0p, pdb0p, p, m1, m2, ier)
	  if(ier .gt. MAX_DIGITS)
     +      print *, "DB0 lost ", ier, " digits (m1 = ",
     +        dsqrt(m1), ")"
	endif
	DB0 = db0p

#ifdef FFCHECK
	call check_B(DB0, DB0_check, "DB0", p, m1, m2)
#endif
	end

************************************************************************
	double complex function C0(p1, p2, p1p2, m1, m2, m3)
	implicit none
	double precision p1, p2, p1p2, m1, m2, m3

#ifdef FFCHECK
	external C0_check
#endif

	double precision xpi(6)
	integer ier

	xpi(1) = m1
	xpi(2) = m2
	xpi(3) = m3
	xpi(4) = p1
	xpi(5) = p2
	xpi(6) = p1p2
	ier = 0
	call ffxc0(C0, xpi, ier)
	if(ier .gt. MAX_DIGITS) then
	  ier = 0
	  call ffxc0r(C0, xpi, ier)
	  if(ier .gt. MAX_DIGITS)
     +      print *, "C0 lost ", ier, " digits (m1 = ",
     +        dsqrt(m1), ")"
	endif

#ifdef FFCHECK
	call check_C(C0, C0_check, "C0", p1, p2, p1p2, m1, m2, m3)
#endif
	end

************************************************************************
	double complex function D0(p1, p2, p3, p4,
     +    p1p2, p2p3, m1, m2, m3, m4)
	implicit none
	double precision p1, p2, p3, p4, p1p2, p2p3
	double precision m1, m2, m3, m4

#ifdef FFCHECK
	external D0_check
#endif

	double precision xpi(13)
	integer ier

	xpi(1) = m1
	xpi(2) = m2
	xpi(3) = m3
	xpi(4) = m4
	xpi(5) = p1 
	xpi(6) = p2
	xpi(7) = p3
	xpi(8) = p4
	xpi(9) = p1p2
	xpi(10) = p2p3
	xpi(11) = 0D0
	xpi(12) = 0D0
	xpi(13) = 0D0

	if(p1p2 .lt. 0D0 .and. p2p3 .lt. 0D0 .and.
     +     p1p2 .gt. p2p3 .and.
     +     m1 .eq. m2 .and. m1 .eq. m3 .and. m1 .eq. m4) then
	  xpi(9) = p2p3
	  xpi(10) = p1p2
	endif
	if(m1 + m2 + m3 + m4 .lt. MINMASS) then
	  call ffxd0m0(D0, xpi, ier)
	else
	  ier = 0
	  call ffxd0(D0, xpi, ier)
	  if(ier .gt. MAX_DIGITS) then
	    ier = 0
	    call ffxd0r(D0, xpi, ier)
	    if(ier .gt. MAX_DIGITS)
     +        print *, "D0 lost ", ier, " digits (m1 = ",
     +          dsqrt(m1), ")"
	  endif
	endif

#ifdef FFCHECK
	call check_D(D0, D0_check, "D0",
     +    p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)
#endif
	end

************************************************************************
	double complex function B1(p, m1, m2)
	implicit none
	double precision p, m1, m2
#include "ff.h"

	double precision mudim, divergence
	common /cutoff/ mudim, divergence

	double complex a0i(2), b0p, b1p
	common /bpave/ a0i, b0p, b1p

#ifdef FFCHECK
	external B1_check
#endif

	integer ier

	ldot = .TRUE.
	ier = 0
	call ffxb0(b0p, divergence, mudim, p, m1, m2, ier)
	call ffxa0(a0i(1), divergence, mudim, m1, ier)
	call ffxa0(a0i(2), divergence, mudim, m2, ier)
	call ffxb1(B1, b0p, a0i, p, m1, m2, fpij2, ier)
	if(ier .gt. MAX_DIGITS)
     +    print *, "B1 lost ", ier, " digits (m1 = ",
     +      dsqrt(m1), ")"
	ldot = .FALSE.

#ifdef FFCHECK
	call check_B(B1, B1_check, "B1", p, m1, m2)
#endif
	end

************************************************************************
	double complex function DB1(p, m1, m2)
	implicit none
	double precision p, m1, m2

#ifdef FFCHECK
	external DB1_check
#endif

	integer ier

	ier = 0
	call ffxdb1(DB1, p, m1, m2, ier)

#ifdef FFCHECK
	call check_B(DB1, DB1_check, "DB1", p, m1, m2)
#endif
	end

************************************************************************
	double complex function B00(p, m1, m2)
	implicit none
	double precision p, m1, m2
#include "ff.h"

	double precision mudim, divergence
	common /cutoff/ mudim, divergence

	double complex a0i(2), b0p, b1p
	common /bpave/ a0i, b0p, b1p

	double complex B1
	external B1

#ifdef FFCHECK
	external B00_check
#endif

	double complex b2(2)
	integer ier

	b1p = B1(p, m1, m2)
	ier = 0
	call ffxb2p(b2, b1p, b0p, a0i, p, m1, m2, fpij2, ier)
	if(ier .gt. MAX_DIGITS)
     +    print *,"B00 lost ", ier, " digits (m1 = ",
     +      dsqrt(m1), ")"
	B00 = b2(2)

#ifdef FFCHECK
	call check_B(B00, B00_check, "B00", p, m1, m2)
#endif
	end

************************************************************************
	double complex function DB00(p, m1, m2)
	implicit none
	double precision p, m1, m2

	double complex B1, DB0, DB1
	external B1, DB0, DB1

#ifdef FFCHECK
	external DB00_check
#endif

	DB00 =
     +    1D0/6D0*(2D0*m1*DB0(p, m1, m2) + B1(p, m1, m2) +
     +      (p + m1 - m2)*DB1(p, m1, m2)) -
     +    1D0/18D0

#ifdef FFCHECK
	call check_B(DB00, DB00_check, "DB00", p, m1, m2)
#endif
	end

************************************************************************
	double complex function B11(p, m1, m2)
	implicit none
	double precision p, m1, m2
#include "ff.h"

	double precision mudim, divergence
	common /cutoff/ mudim, divergence

	double complex a0i(2), b0p, b1p
	common /bpave/ a0i, b0p, b1p

	double complex B1
	external B1

#ifdef FFCHECK
	external B11_check
#endif

	double complex b2(2)
	integer ier

	b1p = B1(p, m1, m2)
	ier = 0
	call ffxb2p(b2, b1p, b0p, a0i, p, m1, m2, fpij2, ier)
	if(ier .gt. MAX_DIGITS)
     +    print *, "B11 lost ", ier, " digits (m1 = ",
     +      dsqrt(m1), ")"
	B11 = b2(1)

#ifdef FFCHECK
	call check_B(B11, B11_check, "B11", p, m1, m2)
#endif
	end

************************************************************************
	double complex function DB11(p, m1, m2)
	implicit none
	double precision p, m1, m2

#ifdef FFCHECK
	external DB11_check
#endif

	integer ier

	ier = 0
	call ffxdb11(DB11, p, m1, m2, ier)

#ifdef FFCHECK
	call check_B(DB11, DB11_check, "DB11", p, m1, m2)
#endif
	end

************************************************************************
	subroutine Cijk(P, C)
	implicit none
	double precision P(6)
	double complex C(13)

	double complex a0i(2), b0p, b1p
	common /bpave/ a0i, b0p, b1p

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

	M12 = .5D0*(XPREC(P(2)) - XPREC(P(3)) - XPREC(P(1)))
	det2 = 2D0*(XPREC(P(1))*XPREC(P(3)) - M12*M12)
	M12 = M12/det2
	M11 = XPREC(P(3))/det2
	M22 = XPREC(P(1))/det2
	f1 = XPREC(P(4)) - XPREC(P(5)) + XPREC(P(1))
	f2 = XPREC(P(4)) - XPREC(P(6)) + XPREC(P(3))

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
	s1 = -(f1*C(cc1) + b123) - 2D0*C(cc00)
	s2 = -(f2*C(cc1) + b123 - b112)
	C(cc11) = lc1()
	C(cc12) = lc2()
	b023 = b023 + b123
	s1 = b023 + b113 - f1*C(cc2)
	s2 = b023 - f2*C(cc2) - 2D0*C(cc00)
	C(cc22) = lc2()

	C(cc001) = (P(4)*C(cc1) + .5D0*(f1*C(cc11) +
     +    f2*C(cc12) + b123))/3D0 - 1D0/18D0
	C(cc002) = (P(4)*C(cc2) + .5D0*(f1*C(cc12) +
     +    f2*C(cc22) - b023))/3D0 - 1D0/18D0
	s1 = -(b1123 + f1*C(cc11)) - 4D0*C(cc001)
	s2 = -(b1123 + f2*C(cc11) - b1112)
	C(cc111) = lc1()
	C(cc112) = lc2()
	b1123 = b1123 + b023 + b123
	s1 = -(b1123 + f1*C(cc22) - b1113)
	s2 = -(b1123 + f2*C(cc22)) - 4D0*C(cc002)
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

************************************************************************
	integer function Cget(p1, p2, p1p2, m1, m2, m3)
	implicit none
	double precision p1, p2, p1p2, m1, m2, m3
#include "tools.h"

	integer i, j

	Cpara(1, Cptr) = p1
	Cpara(2, Cptr) = p2
	Cpara(3, Cptr) = p1p2
	Cpara(4, Cptr) = m1
	Cpara(5, Cptr) = m2
	Cpara(6, Cptr) = m3

	do 20 i = 1, Cptr - 1
	  do 30 j = 1, 6
	    if(Cpara(j, Cptr) .ne. Cpara(j, i)) goto 20
30	  continue
	  Cget = i
	  return
20	continue

	Cget = Cptr
	call Cijk(Cpara(1, Cptr), Cval(1, Cptr))
	Cptr = Cptr + 1
	if(Cptr .gt. CSTORE) then
	  print *, "Cget: Not enough memory. Increase CSTORE."
	  Cptr = CSTORE
	endif
	end

************************************************************************
	double complex function C0i(id,
     +    p1, p2, p1p2, m1, m2, m3)
	implicit none
	integer id
	double precision p1, p2, p1p2, m1, m2, m3
#include "tools.h"

	integer Cget
	external Cget

	C0i = Cval(id, Cget(p1, p2, p1p2, m1, m2, m3))
	end

************************************************************************
	subroutine Dijkl(P, D)
	implicit none
	double precision P(10)
	double complex D(46)
#include "tools.h"

	integer Cget
	double complex D0
	external Cget, D0

	XREAL f1, f2, f3, det3
	XREAL M11, M12, M13, M22, M23, M33
	double complex s1, s2, s3
	double complex c1, c2, c3, c4, c5
	integer C234, C134, C124, C123

#ifndef LANGUAGE_FORTRAN_90
	double complex lc1, lc2, lc3
	lc1() = dcmplx(
     +    M11*dble(s1) + M12*dble(s2) + M13*dble(s3),
     +    M11*dimag(s1) + M12*dimag(s2) + M13*dimag(s3))
	lc2() = dcmplx(
     +    M12*dble(s1) + M22*dble(s2) + M23*dble(s3),
     +    M12*dimag(s1) + M22*dimag(s2) + M23*dimag(s3))
	lc3() = dcmplx(
     +    M13*dble(s1) + M23*dble(s2) + M33*dble(s3),
     +    M13*dimag(s1) + M23*dimag(s2) + M33*dimag(s3))
#endif

	f1 = .5D0*(XPREC(P(1)) + XPREC(P(5)) - XPREC(P(2)))
	f2 = .5D0*(XPREC(P(1)) + XPREC(P(4)) - XPREC(P(6)))
	f3 = .5D0*(XPREC(P(5)) + XPREC(P(4)) - XPREC(P(3)))
	M22 = XPREC(P(1))*XPREC(P(4)) - f2*f2
	M23 = f1*f2 - XPREC(P(1))*f3
	det3 = 2D0*(XPREC(P(5))*M22 - XPREC(P(4))*f1*f1 +
     +    f3*(f1*f2 + M23))
	M11 = (XPREC(P(4))*XPREC(P(5)) - f3*f3)/det3
	M12 = (f2*f3 - XPREC(P(4))*f1)/det3
	M13 = (f1*f3 - XPREC(P(5))*f2)/det3
	M22 = M22/det3
	M23 = M23/det3
	M33 = (XPREC(P(1))*XPREC(P(5)) - f1*f1)/det3

	f1 = XPREC(P(1)) + XPREC(P(7)) - XPREC(P(8))
	f2 = XPREC(P(5)) + XPREC(P(7)) - XPREC(P(9))
	f3 = XPREC(P(4)) + XPREC(P(7)) - XPREC(P(10))

	C234 = Cget(P(2), P(3), P(6), P(8), P(9), P(10))
	C134 = Cget(P(5), P(3), P(4), P(7), P(9), P(10))
	C124 = Cget(P(1), P(6), P(4), P(7), P(8), P(10))
	C123 = Cget(P(1), P(2), P(5), P(7), P(8), P(9))

	c1 = Cval(cc0, C234) + Cval(cc1, C234) + Cval(cc2, C234)
	c2 = Cval(cc12, C234) + Cval(cc22, C234) + Cval(cc2, C234)
	c3 = 2D0*(Cval(cc1, C234) + c2) + Cval(cc11, C234) -
     +    Cval(cc22, C234) + Cval(cc0, C234)
	c4 = Cval(cc1, C234) + Cval(cc11, C234) + Cval(cc12,C234)
	c5 = Cval(cc111, C234) + Cval(cc112, C234) + Cval(cc11, C234)

	D(dd0) = D0(P(1), P(2), P(3), P(4), P(5), P(6),
     +    P(7), P(8), P(9), P(10))
	s1 = Cval(cc0, C134) - Cval(cc0, C234) - f1*D(dd0)
	s2 = Cval(cc0, C124) - Cval(cc0, C234) - f2*D(dd0)
	s3 = Cval(cc0, C123) - Cval(cc0, C234) - f3*D(dd0)
	D(dd1) = lc1()
	D(dd2) = lc2()
	D(dd3) = lc3()

	D(dd00) = P(7)*D(dd0) + .5D0*(Cval(cc0,C234) +
     +    D(dd1)*f1 + D(dd2)*f2 + D(dd3)*f3)
	s1 = c1 - f1*D(dd1) - 2D0*D(dd00)
	s2 = Cval(cc1, C124) + c1 - f2*D(dd1)
	s3 = Cval(cc1, C123) + c1 - f3*D(dd1)
	D(dd11) = lc1()
	s1 = Cval(cc1, C134) - Cval(cc1, C234) - f1*D(dd2)
	s2 = -(Cval(cc1, C234) + f2*D(dd2)) - 2D0*D(dd00)
	s3 = Cval(cc2, C123) - Cval(cc1, C234) - f3*D(dd2)
	D(dd12) = lc1()
	D(dd22) = lc2()
	s1 = Cval(cc2, C134) - Cval(cc2, C234) - f1*D(dd3)
	s2 = Cval(cc2, C124) - Cval(cc2, C234) - f2*D(dd3)
	s3 = -(Cval(cc2, C234) + f3*D(dd3)) - 2D0*D(dd00)
	D(dd13) = lc1()
	D(dd23) = lc2()
	D(dd33) = lc3()

	s1 = Cval(cc00, C134) - Cval(cc00, C234) - f1*D(dd00)
	s2 = Cval(cc00, C124) - Cval(cc00, C234) - f2*D(dd00)
	s3 = Cval(cc00, C123) - Cval(cc00, C234) - f3*D(dd00)
	D(dd001) = lc1()
	D(dd002) = lc2()
	D(dd003) = lc3()
	s1 = -f1*D(dd11) - c3 - 4D0*D(dd001)
	s2 = Cval(cc11,C124) - f2*D(dd11) - c3
	s3 = Cval(cc11,C123) - f3*D(dd11) - c3
	D(dd111) = lc1()
	D(dd112) = lc2()
	D(dd113) = lc3()
	s1 = Cval(cc11, C134) - Cval(cc11, C234) - f1*D(dd22)
	s2 = -(Cval(cc11, C234) + f2*D(dd22)) - 4D0*D(dd002)
	s3 = Cval(cc22, C123) - Cval(cc11, C234) - f3*D(dd22)
	D(dd122) = lc1()
	D(dd222) = lc2()
	D(dd223) = lc3()
	s1 = Cval(cc22, C134) - Cval(cc22, C234) - f1*D(dd33)
	s2 = Cval(cc22, C124) - Cval(cc22, C234) - f2*D(dd33)
	s3 = -(Cval(cc22, C234) + f3*D(dd33)) - 4D0*D(dd003)
	D(dd133) = lc1()
	D(dd233) = lc2()
	D(dd333) = lc3()
	s1 = c2 - f1*D(dd13) - 2D0*D(dd003)
	s2 = Cval(cc12,C124) + c2 - f2*D(dd13)
	s3 = c2 - f3*D(dd13) - 2D0*D(dd001)
	D(dd123) = lc2()

	D(dd0000) = (P(7)*D(dd00) + .5D0*(Cval(cc00, C234) +
     +    f1*D(dd001) + f2*D(dd002) + f3*D(dd003)))/3D0 +
     +    1D0/36D0
	D(dd0011) = (P(7)*D(dd11) + .5D0*(c3 + f1*D(dd111) +
     +    f2*D(dd112) + f3*D(dd113)))/3D0
	D(dd0012) = (P(7)*D(dd12) - .5D0*(c4 - f1*D(dd112) -
     +    f2*D(dd122) - f3*D(dd123)))/3D0
	D(dd0013) = (P(7)*D(dd13) - .5D0*(c2 - f1*D(dd113) -
     +    f2*D(dd123) - f3*D(dd133)))/3D0
	D(dd0022) = (P(7)*D(dd22) + .5D0*(Cval(cc11, C234)+
     +    f1*D(dd122) + f2*D(dd222) + f3*D(dd223)))/3D0
	D(dd0023) = (P(7)*D(dd23) + .5D0*(Cval(cc12, C234)+
     +    f1*D(dd123) + f2*D(dd223) + f3*D(dd233)))/3D0
	D(dd0033) = (P(7)*D(dd33) + .5D0*(Cval(cc22, C234)+
     +    f1*D(dd133) + f2*D(dd233) + f3*D(dd333)))/3D0

	c3 = c2 + Cval(cc112, C234) + Cval(cc122, C234)
	c2 = c2 + c3 - Cval(cc2, C234) + Cval(cc122, C234) +
     +    Cval(cc222, C234)
	c4 = c4 + c1 + c2 + c3 + c4 + c5 + Cval(cc12, C234)

	s1 = c4 - f1*D(dd111) - 6D0*D(dd0011)
	s2 = Cval(cc111, C124) + c4 - f2*D(dd111)
	s3 = Cval(cc111, C123) + c4 - f3*D(dd111)
	D(dd1111) = lc1()
	D(dd1112) = lc2()
	D(dd1113) = lc3()
	s1 = -(c2 + f1*D(dd113)) - 4D0*D(dd0013)
	s2 = Cval(cc112, C124) - c2 - f2*D(dd113)
	s3 = -(c2 + f3*D(dd113)) - 2D0*D(dd0011)
	D(dd1123) = lc2()
	D(dd1133) = lc3()
	s1 = c5 - f1*D(dd122) - 2D0*D(dd0022)
	s2 = c5 - f2*D(dd122) - 4D0*D(dd0012)
	s3 = Cval(cc122, C123) + c5 - f3*D(dd122)
	D(dd1122) = lc1()
	D(dd1223) = lc3()
	s1 = Cval(cc111, C134) - Cval(cc111, C234) - f1*D(dd222)
	s2 = -(Cval(cc111, C234) + f2*D(dd222)) - 6D0*D(dd0022)
	s3 = Cval(cc222, C123) - Cval(cc111, C234) - f3*D(dd222)
	D(dd1222) = lc1()
	D(dd2222) = lc2()
	D(dd2223) = lc3()
	s1 = Cval(cc122, C134) - Cval(cc122, C234) - f1*D(dd233)
	s2 = -(Cval(cc122, C234) + f2*D(dd233)) - 2D0*D(dd0033)
	s3 = -(Cval(cc122, C234) + f3*D(dd233)) - 4D0*D(dd0023)
	D(dd1233) = lc1()
	D(dd2233) = lc2()
	s1 = Cval(cc222, C134) - Cval(cc222, C234) - f1*D(dd333)
	s2 = Cval(cc222, C124) - Cval(cc222, C234) - f2*D(dd333)
	s3 = -(Cval(cc222, C234) + f3*D(dd333)) - 6D0*D(dd0033)
	D(dd1333) = lc1()
	D(dd2333) = lc2()
	D(dd3333) = lc3()

#ifdef LANGUAGE_FORTRAN_90
	contains

	double complex function lc1()
	lc1 = dcmplx(
     +    M11*dble(s1) + M12*dble(s2) + M13*dble(s3),
     +    M11*dimag(s1) + M12*dimag(s2) + M13*dimag(s3))
	end function lc1

	double complex function lc2()
	lc2 = dcmplx(
     +    M12*dble(s1) + M22*dble(s2) + M23*dble(s3),
     +    M12*dimag(s1) + M22*dimag(s2) + M23*dimag(s3))
	end function lc2

	double complex function lc3()
	lc3 = dcmplx(
     +    M13*dble(s1) + M23*dble(s2) + M33*dble(s3),
     +    M13*dimag(s1) + M23*dimag(s2) + M33*dimag(s3))
	end function lc3
#endif
	end

************************************************************************
	integer function Dget(p1, p2, p3, p4, p1p2, p2p3,
     +    m1, m2, m3, m4)
	implicit none
	double precision p1, p2, p3, p4, p1p2, p2p3
	double precision m1, m2, m3, m4
#include "tools.h"

	integer i, j

	Dpara(1, Dptr) = p1
	Dpara(2, Dptr) = p2
	Dpara(3, Dptr) = p3
	Dpara(4, Dptr) = p4
	Dpara(5, Dptr) = p1p2
	Dpara(6, Dptr) = p2p3
	Dpara(7, Dptr) = m1
	Dpara(8, Dptr) = m2
	Dpara(9, Dptr) = m3
	Dpara(10, Dptr) = m4

	do 40 i = 1, Dptr - 1
	  do 50 j = 1, 10
	    if(Dpara(j, Dptr) .ne. Dpara(j, i)) goto 40
50	  continue
	  Dget = i
	  return
40	continue

	Dget = Dptr
	call Dijkl(Dpara(1, Dptr), Dval(1, Dptr))
	Dptr = Dptr + 1
	if(Dptr .gt. DSTORE) then
	  print *, "Dget: Not enough memory. Increase DSTORE."
	  Dptr = DSTORE
	endif
	end

************************************************************************
	double complex function D0i(id, p1, p2, p3, p4,
     +    p1p2, p2p3, m1, m2, m3, m4)
	implicit none
	integer id
	double precision p1, p2, p3, p4, p1p2, p2p3
	double precision m1, m2, m3, m4
#include "tools.h"

	integer Dget
	external Dget

	D0i = Dval(id,
     +    Dget(p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4))
	end

************************************************************************
	subroutine bcaini
	implicit none
#include "ff.h"
#include "tools.h"

	double precision mudim, divergence
	common /cutoff/ mudim, divergence

	double precision delta
	common /ffcut/ delta

#ifdef FFCHECK
	double precision maxdev
	common /ffcheck/ maxdev

	maxdev = MAX_DEVIATION
#endif

	call ffini
	ltest = .FALSE.
	lwarn = .FALSE.
	mudim = MUDIM
	divergence = DELTA
	delta = LAMBDA
	Cptr = 1
	Dptr = 1
	end

#undef IN_TOOLS_F

#endif
