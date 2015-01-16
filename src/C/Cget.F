* Cget.F
* the three-point tensor coefficients
* this file is part of LoopTools
* improvements by M. Rauch
* last modified 13 Sep 12 th

#include "externals.h"
#include "types.h"

#define npoint 3
#include "defs.h"


	memindex function XCget(p1, p2, p1p2, m1, m2, m3)
	implicit none
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	memindex cacheindex
	external cacheindex, XCcoeff

#ifdef COMPLEXPARA
	memindex Cget
	external Cget
#endif

	DVAR para(1,Pcc)

	P(1) = p1
	P(2) = p2
	P(3) = p1p2
#ifdef COMPLEXPARA
	if( abs(Im(P(1))) + abs(Im(P(2))) + abs(Im(P(3))) .gt. 0 )
     &    print *, "CgetC: Complex momenta not implemented"
#endif

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
	M(3) = m3
	if( abs(M(3)) .lt. minmass ) M(3) = 0
#ifdef COMPLEXPARA
	if( abs(Im(M(1))) + abs(Im(M(2))) + abs(Im(M(3))) .eq. 0 ) then
	  XCget = Cget(p1, p2, p1p2, m1, m2, m3) - offsetC
	  return
	endif
#endif

	XCget = cacheindex(para, Cval(1,0), XCcoeff, RC*Pcc, Ncc)
	end

************************************************************************

	subroutine XCput(res, p1, p2, p1p2, m1, m2, m3)
	implicit none
	ComplexType res(*)
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	external XCcoeff

	DVAR para(1,Pcc)

	P(1) = p1
	P(2) = p2
	P(3) = p1p2
#ifdef COMPLEXPARA
	if( abs(Im(P(1))) + abs(Im(P(2))) + abs(Im(P(3))) .gt. 0 )
     &    print *, "CgetC: Complex momenta not implemented"
#endif

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
	M(3) = m3
	if( abs(M(3)) .lt. minmass ) M(3) = 0
#ifdef COMPLEXPARA
	if( abs(Im(M(1))) + abs(Im(M(2))) + abs(Im(M(3))) .eq. 0 ) then
	  call Cput(res, p1, p2, p1p2, m1, m2, m3)
	  return
	endif
#endif

	call cachecopy(res, para, Cval(1,0), XCcoeff, RC*Pcc, Ncc)
	end

************************************************************************

	ComplexType function XC0i(i, p1, p2, p1p2, m1, m2, m3)
	implicit none
	integer i
	DVAR p1, p2, p1p2, m1, m2, m3

#include "lt.h"

	memindex XCget
	external XCget

	memindex b

	b = XCget(p1, p2, p1p2, m1, m2, m3)
	XC0i = Cval(i,b)
	end

************************************************************************

	subroutine XCcoeff(C, para)
	implicit none
	ComplexType C(*)
	DVAR para(1,*)

#include "lt.h"

	memindex XBget
	external XBget

	DVAR p1, p2, p1p2, m1, m2, m3
	DVAR f1, f2
	QVAR G(2,2)
	ComplexType bsum, b1sum, b00sum, b11sum, in(2)
	integer finite
	memindex B12, B23, B13
	logical dump

#ifdef SOLVE_EIGEN
	QVAR Ginv(2,2)
#define SOLVE_SETUP XInverse(2, G,2, Ginv,2)
#define SOLVE(b) XSolve(2, G,2, Ginv,2, b)
#else
	integer perm(2)
#define IN(i) in(perm(i))
#define SOLVE_SETUP XDecomp(2, G,2, perm)
#define SOLVE(b) XSolve(2, G,2, b)
#endif

	m1 = M(1)
	m2 = M(2)
	m3 = M(3)
	p1 = P(1)
	p2 = P(2)
	p1p2 = P(3)

	finite = 1
	if( lambda .lt. 0 ) finite = 0

	B12 = XBget(p1, m1, m2)
	B23 = XBget(p2, m2, m3)
	B13 = XBget(p1p2, m1, m3)

	serial = serial + 1
	dump = ibits(debugkey, DebugC, 1) .ne. 0 .and.
     &    serial .ge. debugfrom .and. serial .le. debugto
	if( dump ) call XDumpPara(3, para, "Ccoeff")

	f1 = m2
	  f1 = f1 - m1
	  f1 = f1 - p1
	f2 = m3
	  f2 = f2 - m1
	  f2 = f2 - p1p2

	G(1,1) = 2*p1
	G(2,2) = 2*p1p2
	G(1,2) = p1
	  G(1,2) = G(1,2) + p1p2
	  G(1,2) = G(1,2) - p2
	G(2,1) = G(1,2)
	call SOLVE_SETUP

	bsum = Bval(bb0,B23) + Bval(bb1,B23)
	b1sum = Bval(bb1,B23) + Bval(bb11,B23)
	b00sum = Bval(bb00,B23) + Bval(bb001,B23)
	b11sum = Bval(bb11,B23) + Bval(bb111,B23)

	call XC0para(C(cc0), para)

	IN(1) = f1*C(cc0) - Bval(bb0,B23) + Bval(bb0,B13)
	IN(2) = f2*C(cc0) - Bval(bb0,B23) + Bval(bb0,B12)
	call SOLVE(in)
	C(cc1) = in(1)
	C(cc2) = in(2)

	C(cc00) = .5D0*(m1*C(cc0) -
     &    .5D0*(f1*C(cc1) + f2*C(cc2) - Bval(bb0,B23) - finite))

	IN(1) = f1*C(cc1) + bsum - 2*C(cc00)
	IN(2) = f2*C(cc1) + bsum + Bval(bb1,B12)
	call SOLVE(in)
	C(cc11) = in(1)
	C(cc12) = in(2)

	IN(1) = f1*C(cc2) - Bval(bb1,B23) + Bval(bb1,B13)
	IN(2) = f2*C(cc2) - Bval(bb1,B23) - 2*C(cc00)
	call SOLVE(in)
	C(cc12) = .5D0*(C(cc12) + in(1))
	C(cc22) = in(2)

	C(cc001) = 1/3D0*(m1*C(cc1) -
     &    .5D0*(f1*C(cc11) + f2*C(cc12) + bsum + finite/3D0))
	C(cc002) = 1/3D0*(m1*C(cc2) -
     &    .5D0*(f1*C(cc12) + f2*C(cc22) - Bval(bb1,B23) + finite/3D0))

	bsum = bsum + b1sum
	IN(1) = f1*C(cc11) - bsum - 4*C(cc001)
	IN(2) = f2*C(cc11) - bsum + Bval(bb11,B12)
	call SOLVE(in)
	C(cc111) = in(1)
	C(cc112) = in(2)

	IN(1) = f1*C(cc22) - Bval(bb11,B23) + Bval(bb11,B13)
	IN(2) = f2*C(cc22) - Bval(bb11,B23) - 4*C(cc002)
	call SOLVE(in)
	C(cc122) = in(1)
	C(cc222) = in(2)

	C(cc0000) = 1/4D0*(m1*C(cc00) -
     &    .5D0*(f1*C(cc001) + f2*C(cc002) - Bval(bb00,B23) -
     &      finite*(m1 + m2 + m3 - .25D0*(p1 + p2 + p1p2))/6D0))

	IN(1) = f1*C(cc001) + b00sum - 2*C(cc0000)
	IN(2) = f2*C(cc001) + b00sum + Bval(bb001,B12)
	call SOLVE(in)
	C(cc0011) = in(1)
	C(cc0012) = in(2)

	IN(1) = f1*C(cc002) - Bval(bb001,B23) + Bval(bb001,B13)
	IN(2) = f2*C(cc002) - Bval(bb001,B23) - 2*C(cc0000)
	call SOLVE(in)
	C(cc0012) = .5D0*(C(cc0012) + in(1))
	C(cc0022) = in(2)

	bsum = bsum + b1sum + b11sum
	IN(1) = f1*C(cc111) + bsum - 6*C(cc0011)
	IN(2) = f2*C(cc111) + bsum + Bval(bb111,B12)
	call SOLVE(in)
	C(cc1111) = in(1)
	C(cc1112) = in(2)

	IN(1) = f1*C(cc222) - Bval(bb111,B23) + Bval(bb111,B13)
	IN(2) = f2*C(cc222) - Bval(bb111,B23) - 6*C(cc0022)
	call SOLVE(in)
	C(cc1222) = in(1)
	C(cc2222) = in(2)

	IN(1) = f1*C(cc122) + b11sum - 2*C(cc0022)
	IN(2) = f2*C(cc122) + b11sum - 4*C(cc0012)
	call SOLVE(in)
	C(cc1122) = in(1)
	C(cc1222) = .5D0*(C(cc1222) + in(2))

	if( dump ) call XDumpCoeff(3, C)
	end

