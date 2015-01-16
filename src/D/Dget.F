* Dget.F
* retrieve the four-point tensor coefficients
* this file is part of LoopTools
* improvements by M. Rauch
* last modified 26 Jan 14 th

#include "externals.h"
#include "types.h"

#define npoint 4
#include "defs.h"


	memindex function XDget(p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	DVAR p1, p2, p3, p4, p1p2, p2p3
	DVAR m1, m2, m3, m4

#include "lt.h"

	memindex cacheindex
	external cacheindex, XDcoeff

#ifdef COMPLEXPARA
	memindex Dget
	external Dget
#endif

	DVAR para(1,Pdd+3)

	P(1) = p1
	P(2) = p2
	P(3) = p3
	P(4) = p4
	P(5) = p1p2
	P(6) = p2p3
#ifdef COMPLEXPARA
	if( abs(Im(P(1))) + abs(Im(P(2))) +
     &      abs(Im(P(3))) + abs(Im(P(4))) +
     &      abs(Im(P(5))) + abs(Im(P(6))) .gt. 0 )
     &    print *, "DgetC: Complex momenta not implemented"
#endif

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
	M(3) = m3
	if( abs(M(3)) .lt. minmass ) M(3) = 0
	M(4) = m4
	if( abs(M(4)) .lt. minmass ) M(4) = 0
#ifdef COMPLEXPARA
	if( abs(Im(M(1))) + abs(Im(M(2))) +
     &      abs(Im(M(3))) + abs(Im(M(4))) .eq. 0 ) then
	  XDget = Dget(p1, p2, p3, p4, p1p2, p2p3,
     &      m1, m2, m3, m4) - offsetC
	  return
	endif
#endif

	XDget = cacheindex(para, Dval(1,0), XDcoeff, RC*Pdd, Ndd)
	end

************************************************************************

	subroutine XDput(res, p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	ComplexType res(*)
	DVAR p1, p2, p3, p4, p1p2, p2p3
	DVAR m1, m2, m3, m4

#include "lt.h"

	external XDcoeff

	DVAR para(1,Pdd+3)

	P(1) = p1
	P(2) = p2
	P(3) = p3
	P(4) = p4
	P(5) = p1p2
	P(6) = p2p3
#ifdef COMPLEXPARA
	if( abs(Im(P(1))) + abs(Im(P(2))) +
     &      abs(Im(P(3))) + abs(Im(P(4))) +
     &      abs(Im(P(5))) + abs(Im(P(6))) .gt. 0 )
     &    print *, "DputC: Complex momenta not implemented"
#endif

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
	M(3) = m3
	if( abs(M(3)) .lt. minmass ) M(3) = 0
	M(4) = m4
	if( abs(M(4)) .lt. minmass ) M(4) = 0
#ifdef COMPLEXPARA
	if( abs(Im(M(1))) + abs(Im(M(2))) +
     &      abs(Im(M(3))) + abs(Im(M(4))) .eq. 0 ) then
	  call Dput(res, p1, p2, p3, p4, p1p2, p2p3,
     &      m1, m2, m3, m4)
	  return
	endif
#endif

	call cachecopy(res, para, Dval(1,0), XDcoeff, RC*Pdd, Ndd)
	end

************************************************************************

	ComplexType function XD0i(i, p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	integer i
	DVAR p1, p2, p3, p4, p1p2, p2p3
	DVAR m1, m2, m3, m4

#include "lt.h"

	memindex XDget
	external XDget

	memindex b

	b = XDget(p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)
	XD0i = Dval(i+epsi,b)
	end

************************************************************************

	ComplexType function XD0(p1, p2, p3, p4, p1p2, p2p3,
     &    m1, m2, m3, m4)
	implicit none
	DVAR p1, p2, p3, p4, p1p2, p2p3
	DVAR m1, m2, m3, m4

#include "lt.h"

	ComplexType XD0i
	external XD0i

	XD0 = XD0i(dd0, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)
	end

************************************************************************

	subroutine XDcoeff(D, para)
	implicit none
	ComplexType D(*)
	DVAR para(1,*)

#include "lt.h"

	memindex XCget
	external XCget

	DVAR p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4
	DVAR f1, f2, f3
	QVAR G(3,3)
	ComplexType c0sum(0:2), c1sum(0:2), c2sum(0:2), csum(0:2)
	ComplexType c00sum(0:2), c11sum(0:2), c12sum(0:2), c22sum(0:2)
	ComplexType in(0:2,3)
	memindex C234, C134, C124, C123
	logical dump

#ifdef SOLVE_EIGEN
	QVAR Ginv(3,3)
#define IN(i) in(:,i)
#define SOLVE_SETUP XInverse(3, G,3, Ginv,3)
#define SOLVE(b) XSolve(3, G,3, Ginv,3, b)
#else
	integer perm(3)
#define IN(i) in(:,perm(i))
#define SOLVE_SETUP XDecomp(3, G,3, perm)
#define SOLVE(b) XSolve(3, G,3, b)
#endif
#define OUT(i) in(:,i)

	integer finite(0:2)
	data finite /1, 0, 0/

	m1 = M(1)
	m2 = M(2)
	m3 = M(3)
	m4 = M(4)
	p1 = P(1)
	p2 = P(2)
	p3 = P(3)
	p4 = P(4)
	p1p2 = P(5)
	p2p3 = P(6)

	C234 = XCget(p2, p3, p2p3, m2, m3, m4)
	C134 = XCget(p1p2, p3, p4, m1, m3, m4)
	C124 = XCget(p1, p2p3, p4, m1, m2, m4)
	C123 = XCget(p1, p2, p1p2, m1, m2, m3)

	serial = serial + 1
	dump = ibits(debugkey, DebugD, 1) .ne. 0 .and.
     &    serial .ge. debugfrom .and. serial .le. debugto
	if( dump ) call XDumpPara(4, para, "Dcoeff")

	f1 = m2
	  f1 = f1 - m1
	  f1 = f1 - p1
	f2 = m3
	  f2 = f2 - m1
	  f2 = f2 - p1p2
	f3 = m4
	  f3 = f3 - m1
	  f3 = f3 - p4

	G(1,1) = 2*p1
	G(2,2) = 2*p1p2
	G(3,3) = 2*p4
	G(1,2) = p1
	  G(1,2) = G(1,2) + p1p2
	  G(1,2) = G(1,2) - p2
	G(2,1) = G(1,2)
	G(1,3) = p1
	  G(1,3) = G(1,3) + p4
	  G(1,3) = G(1,3) - p2p3
	G(3,1) = G(1,3)
	G(2,3) = p1p2
	  G(2,3) = G(2,3) - p3
	  G(2,3) = G(2,3) + p4
	G(3,2) = G(2,3)
	call SOLVE_SETUP

	c0sum = Cval(CC0,C234) + Cval(CC1,C234) + Cval(CC2,C234)
	c1sum = Cval(CC1,C234) + Cval(CC11,C234) + Cval(CC12,C234)
	c2sum = Cval(CC2,C234) + Cval(CC12,C234) + Cval(CC22,C234)
	csum = c0sum + c1sum + c2sum

	c00sum = Cval(CC00,C234) +
     &    Cval(CC001,C234) + Cval(CC002,C234)
	c11sum = Cval(CC11,C234) +
     &    Cval(CC111,C234) + Cval(CC112,C234)
	c12sum = Cval(CC12,C234) +
     &    Cval(CC112,C234) + Cval(CC122,C234)
	c22sum = Cval(CC22,C234) +
     &    Cval(CC122,C234) + Cval(CC222,C234)

	call XD0para(D(DD0), para)

	IN(1) = f1*D(DD0) - Cval(CC0,C234) + Cval(CC0,C134)
	IN(2) = f2*D(DD0) - Cval(CC0,C234) + Cval(CC0,C124)
	IN(3) = f3*D(DD0) - Cval(CC0,C234) + Cval(CC0,C123)
	call SOLVE(in)
	D(DD1) = OUT(1)
	D(DD2) = OUT(2)
	D(DD3) = OUT(3)

	D(DD00) = m1*D(DD0) - .5D0*
     &    (D(DD1)*f1 + D(DD2)*f2 + D(DD3)*f3 - Cval(CC0,C234))

	IN(1) = f1*D(DD1) + c0sum - 2*D(DD00)
	IN(2) = f2*D(DD1) + c0sum + Cval(CC1,C124)
	IN(3) = f3*D(DD1) + c0sum + Cval(CC1,C123)
	call SOLVE(in)
	D(DD11) = OUT(1)
	D(DD12) = OUT(2)
	D(DD13) = OUT(3)

	IN(1) = f1*D(DD2) - Cval(CC1,C234) + Cval(CC1,C134)
	IN(2) = f2*D(DD2) - Cval(CC1,C234) - 2*D(DD00)
	IN(3) = f3*D(DD2) - Cval(CC1,C234) + Cval(CC2,C123)
	call SOLVE(in)
	D(DD12) = .5D0*(D(DD12) + OUT(1))
	D(DD22) = OUT(2)
	D(DD23) = OUT(3)

	IN(1) = f1*D(DD3) - Cval(CC2,C234) + Cval(CC2,C134)
	IN(2) = f2*D(DD3) - Cval(CC2,C234) + Cval(CC2,C124)
	IN(3) = f3*D(DD3) - Cval(CC2,C234) - 2*D(DD00)
	call SOLVE(in)
	D(DD13) = .5D0*(D(DD13) + OUT(1))
	D(DD23) = .5D0*(D(DD23) + OUT(2))
	D(DD33) = OUT(3)

	IN(1) = f1*D(DD00) - Cval(CC00,C234) + Cval(CC00,C134)
	IN(2) = f2*D(DD00) - Cval(CC00,C234) + Cval(CC00,C124)
	IN(3) = f3*D(DD00) - Cval(CC00,C234) + Cval(CC00,C123)
	call SOLVE(in)
	D(DD001) = OUT(1)
	D(DD002) = OUT(2)
	D(DD003) = OUT(3)

	IN(1) = f1*D(DD11) - csum - 4*D(DD001)
	IN(2) = f2*D(DD11) - csum + Cval(CC11,C124)
	IN(3) = f3*D(DD11) - csum + Cval(CC11,C123)
	call SOLVE(in)
	D(DD111) = OUT(1)
	D(DD112) = OUT(2)
	D(DD113) = OUT(3)

	IN(1) = f1*D(DD22) - Cval(CC11,C234) + Cval(CC11,C134)
	IN(2) = f2*D(DD22) - Cval(CC11,C234) - 4*D(DD002)
	IN(3) = f3*D(DD22) - Cval(CC11,C234) + Cval(CC22,C123)
	call SOLVE(in)
	D(DD122) = OUT(1)
	D(DD222) = OUT(2)
	D(DD223) = OUT(3)

	IN(1) = f1*D(DD33) - Cval(CC22,C234) + Cval(CC22,C134)
	IN(2) = f2*D(DD33) - Cval(CC22,C234) + Cval(CC22,C124)
	IN(3) = f3*D(DD33) - Cval(CC22,C234) - 4*D(DD003)
	call SOLVE(in)
	D(DD133) = OUT(1)
	D(DD233) = OUT(2)
	D(DD333) = OUT(3)

	IN(1) = f1*D(DD13) + c2sum - 2*D(DD003)
	IN(2) = f2*D(DD13) + c2sum + Cval(CC12,C124)
	IN(3) = f3*D(DD13) + c2sum - 2*D(DD001)
	call SOLVE(in)
	D(DD113) = .5D0*(D(DD113) + OUT(1))
	D(DD123) = OUT(2)
	D(DD133) = .5D0*(D(DD133) + OUT(3))

	D(DD0000) = 1/3D0*(m1*D(DD00) -
     &    .5D0*(f1*D(DD001) + f2*D(DD002) + f3*D(DD003) -
     &      Cval(CC00,C234) - finite/6D0))
	D(DD0011) = 1/3D0*(m1*D(DD11) -
     &    .5D0*(f1*D(DD111) + f2*D(DD112) + f3*D(DD113) - csum))
	D(DD0012) = 1/3D0*(m1*D(DD12) -
     &    .5D0*(f1*D(DD112) + f2*D(DD122) + f3*D(DD123) + c1sum))
	D(DD0013) = 1/3D0*(m1*D(DD13) -
     &    .5D0*(f1*D(DD113) + f2*D(DD123) + f3*D(DD133) + c2sum))
	D(DD0022) = 1/3D0*(m1*D(DD22) -
     &    .5D0*(f1*D(DD122) + f2*D(DD222) + f3*D(DD223) -
     &      Cval(CC11,C234)))
	D(DD0023) = 1/3D0*(m1*D(DD23) -
     &    .5D0*(f1*D(DD123) + f2*D(DD223) + f3*D(DD233) -
     &      Cval(CC12,C234)))
	D(DD0033) = 1/3D0*(m1*D(DD33) -
     &    .5D0*(f1*D(DD133) + f2*D(DD233) + f3*D(DD333) -
     &      Cval(CC22,C234)))

	c1sum = c1sum + c11sum + c12sum
	c2sum = c2sum + c12sum + c22sum
	csum = csum + c1sum + c2sum

	IN(1) = f1*D(DD111) + csum - 6*D(DD0011)
	IN(2) = f2*D(DD111) + csum + Cval(CC111,C124)
	IN(3) = f3*D(DD111) + csum + Cval(CC111,C123)
	call SOLVE(in)
	D(DD1111) = OUT(1)
	D(DD1112) = OUT(2)
	D(DD1113) = OUT(3)

	IN(1) = f1*D(DD113) - c2sum - 4*D(DD0013)
	IN(2) = f2*D(DD113) - c2sum + Cval(CC112,C124)
	IN(3) = f3*D(DD113) - c2sum - 2*D(DD0011)
	call SOLVE(in)
	D(DD1113) = .5D0*(D(DD1113) + OUT(1))
	D(DD1123) = OUT(2)
	D(DD1133) = OUT(3)

	IN(1) = f1*D(DD122) + c11sum - 2*D(DD0022)
	IN(2) = f2*D(DD122) + c11sum - 4*D(DD0012)
	IN(3) = f3*D(DD122) + c11sum + Cval(CC122,C123)
	call SOLVE(in)
	D(DD1122) = OUT(1)
	D(DD1222) = OUT(2)
	D(DD1223) = OUT(3)

	IN(1) = f1*D(DD222) - Cval(CC111,C234) + Cval(CC111,C134)
	IN(2) = f2*D(DD222) - Cval(CC111,C234) - 6*D(DD0022)
	IN(3) = f3*D(DD222) - Cval(CC111,C234) + Cval(CC222,C123)
	call SOLVE(in)
	D(DD1222) = .5D0*(D(DD1222) + OUT(1))
	D(DD2222) = OUT(2)
	D(DD2223) = OUT(3)

	IN(1) = f1*D(DD233) - Cval(CC122,C234) + Cval(CC122,C134)
	IN(2) = f2*D(DD233) - Cval(CC122,C234) - 2*D(DD0033)
	IN(3) = f3*D(DD233) - Cval(CC122,C234) - 4*D(DD0023)
	call SOLVE(in)
	D(DD1233) = OUT(1)
	D(DD2233) = OUT(2)
	D(DD2333) = OUT(3)

	IN(1) = f1*D(DD333) - Cval(CC222,C234) + Cval(CC222,C134)
	IN(2) = f2*D(DD333) - Cval(CC222,C234) + Cval(CC222,C124)
	IN(3) = f3*D(DD333) - Cval(CC222,C234) - 6*D(DD0033)
	call SOLVE(in)
	D(DD1333) = OUT(1)
	D(DD2333) = .5D0*(D(DD2333) + OUT(2))
	D(DD3333) = OUT(3)

	c00sum = c00sum +
     &    Cval(CC001,C234) + Cval(CC0011,C234) + Cval(CC0012,C234) +
     &    Cval(CC002,C234) + Cval(CC0012,C234) + Cval(CC0022,C234)
	c11sum = c11sum +
     &    Cval(CC111,C234) + Cval(CC1111,C234) + Cval(CC1112,C234) +
     &    Cval(CC112,C234) + Cval(CC1112,C234) + Cval(CC1122,C234)
	c12sum = c12sum +
     &    Cval(CC112,C234) + Cval(CC1112,C234) + Cval(CC1122,C234) +
     &    Cval(CC122,C234) + Cval(CC1122,C234) + Cval(CC1222,C234)
	c22sum = c22sum +
     &    Cval(CC122,C234) + Cval(CC1122,C234) + Cval(CC1222,C234) +
     &    Cval(CC222,C234) + Cval(CC1222,C234) + Cval(CC2222,C234)
	c1sum = c1sum + c11sum + c12sum
	c2sum = c2sum + c12sum + c22sum
	csum = csum + c1sum + c2sum

	IN(1) = f1*D(DD0000) - Cval(CC0000,C234) + Cval(CC0000,C134)
	IN(2) = f2*D(DD0000) - Cval(CC0000,C234) + Cval(CC0000,C124)
	IN(3) = f3*D(DD0000) - Cval(CC0000,C234) + Cval(CC0000,C123)
	call SOLVE(in)
	D(DD00001) = OUT(1)
	D(DD00002) = OUT(2)
	D(DD00003) = OUT(3)

	IN(1) = f1*D(DD0011) - c00sum - 4*D(DD00001)
	IN(2) = f2*D(DD0011) - c00sum + Cval(CC0011,C124)
	IN(3) = f3*D(DD0011) - c00sum + Cval(CC0011,C123)
	call SOLVE(in)
	D(DD00111) = OUT(1)
	D(DD00112) = OUT(2)
	D(DD00113) = OUT(3)

	IN(1) = f1*D(DD0022) - Cval(CC0011,C234) + Cval(CC0011,C134)
	IN(2) = f2*D(DD0022) - Cval(CC0011,C234) - 4*D(DD00002)
	IN(3) = f3*D(DD0022) - Cval(CC0011,C234) + Cval(CC0022,C123)
	call SOLVE(in)
	D(DD00122) = OUT(1)
	D(DD00222) = OUT(2)
	D(DD00223) = OUT(3)

	IN(1) = f1*D(DD0033) - Cval(CC0022,C234) + Cval(CC0022,C134)
	IN(2) = f2*D(DD0033) - Cval(CC0022,C234) + Cval(CC0022,C124)
	IN(3) = f3*D(DD0033) - Cval(CC0022,C234) - 4*D(DD00003)
	call SOLVE(in)
	D(DD00133) = OUT(1)
	D(DD00233) = OUT(2)
	D(DD00333) = OUT(3)

	IN(1) = f1*D(DD0023) - Cval(CC0012,C234) + Cval(CC0012,C134)
	IN(2) = f2*D(DD0023) - Cval(CC0012,C234) - 2*D(DD00003)
	IN(3) = f3*D(DD0023) - Cval(CC0012,C234) - 2*D(DD00002)
	call SOLVE(in)
	D(DD00123) = OUT(1)
	D(DD00223) = .5D0*(D(DD00223) + OUT(2))
	D(DD00233) = .5D0*(D(DD00233) + OUT(3))

	IN(1) = f1*D(DD1111) - csum - 8*D(DD00111)
	IN(2) = f2*D(DD1111) - csum + Cval(CC1111,C124)
	IN(3) = f3*D(DD1111) - csum + Cval(CC1111,C123)
	call SOLVE(in)
	D(DD11111) = OUT(1)
	D(DD11112) = OUT(2)
	D(DD11113) = OUT(3)

	IN(1) = f1*D(DD2222) - Cval(CC1111,C234) + Cval(CC1111,C134)
	IN(2) = f2*D(DD2222) - Cval(CC1111,C234) - 8*D(DD00222)
	IN(3) = f3*D(DD2222) - Cval(CC1111,C234) + Cval(CC2222,C123)
	call SOLVE(in)
	D(DD12222) = OUT(1)
	D(DD22222) = OUT(2)
	D(DD22223) = OUT(3)

	IN(1) = f1*D(DD3333) - Cval(CC2222,C234) + Cval(CC2222,C134)
	IN(2) = f2*D(DD3333) - Cval(CC2222,C234) + Cval(CC2222,C124)
	IN(3) = f3*D(DD3333) - Cval(CC2222,C234) - 8*D(DD00333)
	call SOLVE(in)
	D(DD13333) = OUT(1)
	D(DD23333) = OUT(2)
	D(DD33333) = OUT(3)

	IN(1) = f1*D(DD1122) - c11sum - 4*D(DD00122)
	IN(2) = f2*D(DD1122) - c11sum - 4*D(DD00112)
	IN(3) = f3*D(DD1122) - c11sum + Cval(CC1122,C123)
	call SOLVE(in)
	D(DD11122) = OUT(1)
	D(DD11222) = OUT(2)
	D(DD11223) = OUT(3)

	IN(1) = f1*D(DD1133) - c22sum - 4*D(DD00133)
	IN(2) = f2*D(DD1133) - c22sum + Cval(CC1122,C124)
	IN(3) = f3*D(DD1133) - c22sum - 4*D(DD00113)
	call SOLVE(in)
	D(DD11133) = OUT(1)
	D(DD11233) = OUT(2)
	D(DD11333) = OUT(3)

	IN(1) = f1*D(DD2233) - Cval(CC1122,C234) + Cval(CC1122,C134)
	IN(2) = f2*D(DD2233) - Cval(CC1122,C234) - 4*D(DD00233)
	IN(3) = f3*D(DD2233) - Cval(CC1122,C234) - 4*D(DD00223)
	call SOLVE(in)
	D(DD12233) = OUT(1)
	D(DD22233) = OUT(2)
	D(DD22333) = OUT(3)

	IN(1) = f1*D(DD1123) - c12sum - 4*D(DD00123)
	IN(2) = f2*D(DD1123) - c12sum - 2*D(DD00113)
	IN(3) = f3*D(DD1123) - c12sum - 2*D(DD00112)
	call SOLVE(in)
	D(DD11123) = OUT(1)
	D(DD11223) = .5D0*(D(DD11223) + OUT(2))
	D(DD11233) = .5D0*(D(DD11233) + OUT(3))

	IN(1) = f1*D(DD2223) - Cval(CC1112,C234) + Cval(CC1112,C134)
	IN(2) = f2*D(DD2223) - Cval(CC1112,C234) - 6*D(DD00223)
	IN(3) = f3*D(DD2223) - Cval(CC1112,C234) - 2*D(DD00222)
	call SOLVE(in)
	D(DD12223) = OUT(1)
	D(DD22223) = .5D0*(D(DD22223) + OUT(2))
	D(DD22233) = .5D0*(D(DD22233) + OUT(3))

	IN(1) = f1*D(DD2333) - Cval(CC1222,C234) + Cval(CC1222,C134)
	IN(2) = f2*D(DD2333) - Cval(CC1222,C234) - 2*D(DD00333)
	IN(3) = f3*D(DD2333) - Cval(CC1222,C234) - 6*D(DD00233)
	call SOLVE(in)
	D(DD12333) = OUT(1)
	D(DD22333) = .5D0*(D(DD22333) + OUT(2))
	D(DD23333) = .5D0*(D(DD23333) + OUT(3))

	if( dump ) call XDumpCoeff(4, D)
	end

