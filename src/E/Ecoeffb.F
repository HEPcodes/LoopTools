* Ecoeffb.F
* the five-point tensor coefficients via Passarino-Veltman decomposition
* this file is part of LoopTools
* written by M. Rauch
* last modified 9 Sep 12 th

#include "externals.h"
#include "types.h"

#define npoint 5
#include "defs.h"


	subroutine XEcoeffb(E, para)
	implicit none
	ComplexType E(*)
	DVAR para(1,*)

#include "lt.h"

	memindex XDget
	external XDget

	DVAR p1, p2, p3, p4, p5
	DVAR p1p2, p2p3, p3p4, p4p5, p5p1
	DVAR m1, m2, m3, m4, m5

	DVAR f1, f2, f3, f4
	ComplexType di, d0sum
	ComplexType d1i, d1sum, d2i, d2sum, dii
	ComplexType d1ii, d2ii, diii
	ComplexType d00sum, d22sum, d33sum
	ComplexType in(4)
	memindex Di1, Di2, Di3, Di4, Di5
	logical dump

	QVAR G(4,4), Ginv(4,4)
	common /XInvGramE/ Ginv

#ifdef SOLVE_EIGEN
#define SOLVE_SETUP XInverse(4, G,4, Ginv,4)
#define SOLVE(b) XSolve(4, G,4, Ginv,4, b)
#else
	integer perm(4)
#define IN(i) in(perm(i))
#define SOLVE_SETUP XInverse(4, G,4, Ginv,4, perm)
#define SOLVE(b) XSolve(4, G,4, b)
#endif

	m1 = M(1)
	m2 = M(2)
	m3 = M(3)
	m4 = M(4)
	m5 = M(5)
	p1 = P(1)
	p2 = P(2)
	p3 = P(3)
	p4 = P(4)
	p5 = P(5)
	p1p2 = P(6)
	p2p3 = P(7)
	p3p4 = P(8)
	p4p5 = P(9)
	p5p1 = P(10)

	Di1 = XDget(p2, p3, p4, p5p1, p2p3, p3p4, m2, m3, m4, m5)
	Di2 = XDget(p1p2, p3, p4, p5, p4p5, p3p4, m1, m3, m4, m5)
	Di3 = XDget(p1, p2p3, p4, p5, p4p5, p5p1, m1, m2, m4, m5)
	Di4 = XDget(p1, p2, p3p4, p5, p1p2, p5p1, m1, m2, m3, m5)
	Di5 = XDget(p1, p2, p3, p4p5, p1p2, p2p3, m1, m2, m3, m4)

	serial = serial + 1
	dump = ibits(debugkey, DebugE, 1) .ne. 0 .and.
     &    serial .ge. debugfrom .and. serial .le. debugto
	if( dump ) call XDumpPara(5, para)

	f1 = m2
	  f1 = f1 - m1
	  f1 = f1 - p1
	f2 = m3
	  f2 = f2 - m1
	  f2 = f2 - p1p2
	f3 = m4
	  f3 = f3 - m1
	  f3 = f3 - p4p5
	f4 = m5
	  f4 = f4 - m1
	  f4 = f4 - p5

* build up G and calculate matrix decomposition and inverse Y
	G(1,1) = 2*p1
	G(2,2) = 2*p1p2
	G(3,3) = 2*p4p5
	G(4,4) = 2*p5
	G(1,2) = p1
	  G(1,2) = G(1,2) + p1p2
	  G(1,2) = G(1,2) - p2
	G(2,1) = G(1,2)
	G(1,3) = p1
	  G(1,3) = G(1,3) - p2p3
	  G(1,3) = G(1,3) + p4p5
	G(3,1) = G(1,3)
	G(1,4) = p1
	  G(1,4) = G(1,4) - p5p1
	  G(1,4) = G(1,4) + p5
	G(4,1) = G(1,4)
	G(2,3) = p1p2
	  G(2,3) = G(2,3) - p3
	  G(2,3) = G(2,3) + p4p5
	G(3,2) = G(2,3)
	G(2,4) = p1p2
	  G(2,4) = G(2,4) - p3p4
	  G(2,4) = G(2,4) + p5
	G(4,2) = G(2,4)
	G(3,4) = p5
	  G(3,4) = G(3,4) + p4p5
	  G(3,4) = G(3,4) - p4
	G(4,3) = G(3,4)
	call SOLVE_SETUP

	di = Dval(dd1,Di1) + Dval(dd2,Di1) + Dval(dd3,Di1)
	d0sum = Dval(dd0,Di1) + di

	d1i = Dval(dd11,Di1) + Dval(dd12,Di1) + Dval(dd13,Di1)
	d1sum = Dval(dd1,Di1) + d1i
	d2i = Dval(dd12,Di1) + Dval(dd22,Di1) + Dval(dd23,Di1)
	d2sum = Dval(dd2,Di1) + d2i
	dii = d1i + d2i +
     &    Dval(dd13,Di1) + Dval(dd23,Di1) + Dval(dd33,Di1)

	d1ii = Dval(dd111,Di1) + Dval(dd122,Di1) +
     &    Dval(dd133,Di1) + 2*(Dval(dd112,Di1) +
     &    Dval(dd113,Di1) + Dval(dd123,Di1))
	d2ii = Dval(dd112,Di1) + Dval(dd222,Di1) +
     &    Dval(dd233,Di1) + 2*(Dval(dd122,Di1) +
     &    Dval(dd123,Di1) + Dval(dd223,Di1))
	diii = d1ii + d2ii +
     &    Dval(dd113,Di1) + Dval(dd223,Di1) +
     &    Dval(dd333,Di1) + 2*(Dval(dd123,Di1) +
     &    Dval(dd133,Di1) + Dval(dd233,Di1))

	d00sum = Dval(dd00,Di1) +
     &    Dval(dd001,Di1) + Dval(dd002,Di1) + Dval(dd003,Di1)
	d22sum = Dval(dd22,Di1) +
     &    Dval(dd122,Di1) + Dval(dd222,Di1) + Dval(dd223,Di1)
	d33sum = Dval(dd33,Di1) +
     &    Dval(dd133,Di1) + Dval(dd233,Di1) + Dval(dd333,Di1)

	call XE0para(E(ee0), para)

	IN(1) = f1*E(ee0) - Dval(dd0,Di1) + Dval(dd0,Di2)
	IN(2) = f2*E(ee0) - Dval(dd0,Di1) + Dval(dd0,Di3)
	IN(3) = f3*E(ee0) - Dval(dd0,Di1) + Dval(dd0,Di4)
	IN(4) = f4*E(ee0) - Dval(dd0,Di1) + Dval(dd0,Di5)
	call SOLVE(in)
	E(ee1) = in(1)
	E(ee2) = in(2)
	E(ee3) = in(3)
	E(ee4) = in(4)

	E(ee00) = 0

	IN(1) = f1*E(ee1) + d0sum
	IN(2) = f2*E(ee1) + d0sum + Dval(dd1,Di3)
	IN(3) = f3*E(ee1) + d0sum + Dval(dd1,Di4)
	IN(4) = f4*E(ee1) + d0sum + Dval(dd1,Di5)
	call SOLVE(in)
	E(ee11) = in(1)
	E(ee12) = in(2)
	E(ee13) = in(3)
	E(ee14) = in(4)

	IN(1) = f1*E(ee2) - Dval(dd1,Di1) + Dval(dd1,Di2)
	IN(2) = f2*E(ee2) - Dval(dd1,Di1)
	IN(3) = f3*E(ee2) - Dval(dd1,Di1) + Dval(dd2,Di4)
	IN(4) = f4*E(ee2) - Dval(dd1,Di1) + Dval(dd2,Di5)
	call SOLVE(in)
	E(ee12) = .5D0*(E(ee12) + in(1))
	E(ee22) = in(2)
	E(ee23) = in(3)
	E(ee24) = in(4)

	IN(1) = f1*E(ee3) - Dval(dd2,Di1) + Dval(dd2,Di2)
	IN(2) = f2*E(ee3) - Dval(dd2,Di1) + Dval(dd2,Di3)
	IN(3) = f3*E(ee3) - Dval(dd2,Di1)
	IN(4) = f4*E(ee3) - Dval(dd2,Di1) + Dval(dd3,Di5)
	call SOLVE(in)
	E(ee13) = .5D0*(E(ee13) + in(1))
	E(ee23) = .5D0*(E(ee23) + in(2))
	E(ee33) = in(3)
	E(ee34) = in(4)

	IN(1) = f1*E(ee4) - Dval(dd3,Di1) + Dval(dd3,Di2)
	IN(2) = f2*E(ee4) - Dval(dd3,Di1) + Dval(dd3,Di3)
	IN(3) = f3*E(ee4) - Dval(dd3,Di1) + Dval(dd3,Di4)
	IN(4) = f4*E(ee4) - Dval(dd3,Di1)
	call SOLVE(in)
	E(ee14) = .5D0*(E(ee14) + in(1))
	E(ee24) = .5D0*(E(ee24) + in(2))
	E(ee34) = .5D0*(E(ee34) + in(3))
	E(ee44) = in(4)

	E(ee001) = 0
	E(ee002) = 0
	E(ee003) = 0
	E(ee004) = 0

	d0sum = d0sum + di + dii
	IN(1) = f1*E(ee11) - d0sum -
     &    2*Ginv(1,1)*(Dval(dd00,Di1) - Dval(dd00,Di2))
	IN(2) = f2*E(ee11) - d0sum + Dval(dd11,Di3) -
     &    2*Ginv(1,1)*(Dval(dd00,Di1) - Dval(dd00,Di3))
	IN(3) = f3*E(ee11) - d0sum + Dval(dd11,Di4) -
     &    2*Ginv(1,1)*(Dval(dd00,Di1) - Dval(dd00,Di4))
	IN(4) = f4*E(ee11) - d0sum + Dval(dd11,Di5) -
     &    2*Ginv(1,1)*(Dval(dd00,Di1) - Dval(dd00,Di5))
	call SOLVE(in)
	E(ee111) = in(1)
	E(ee112) = in(2)
	E(ee113) = in(3)
	E(ee114) = in(4)

	IN(1) = f1*E(ee22) - Dval(dd11,Di1) + Dval(dd11,Di2) -
     &    2*Ginv(2,2)*(Dval(dd00,Di1) - Dval(dd00,Di2))
	IN(2) = f2*E(ee22) - Dval(dd11,Di1) -
     &    2*Ginv(2,2)*(Dval(dd00,Di1) - Dval(dd00,Di3))
	IN(3) = f3*E(ee22) - Dval(dd11,Di1) + Dval(dd22,Di4) -
     &    2*Ginv(2,2)*(Dval(dd00,Di1) - Dval(dd00,Di4))
	IN(4) = f4*E(ee22) - Dval(dd11,Di1) + Dval(dd22,Di5) -
     &    2*Ginv(2,2)*(Dval(dd00,Di1) - Dval(dd00,Di5))
	call SOLVE(in)
	E(ee122) = in(1)
	E(ee222) = in(2)
	E(ee223) = in(3)
	E(ee224) = in(4)

	IN(1) = f1*E(ee33) - Dval(dd22,Di1) + Dval(dd22,Di2) -
     &    2*Ginv(3,3)*(Dval(dd00,Di1) - Dval(dd00,Di2))
	IN(2) = f2*E(ee33) - Dval(dd22,Di1) + Dval(dd22,Di3) -
     &    2*Ginv(3,3)*(Dval(dd00,Di1) - Dval(dd00,Di3))
	IN(3) = f3*E(ee33) - Dval(dd22,Di1) -
     &    2*Ginv(3,3)*(Dval(dd00,Di1) - Dval(dd00,Di4))
	IN(4) = f4*E(ee33) - Dval(dd22,Di1) + Dval(dd33,Di5) -
     &    2*Ginv(3,3)*(Dval(dd00,Di1) - Dval(dd00,Di5))
	call SOLVE(in)
	E(ee133) = in(1)
	E(ee233) = in(2)
	E(ee333) = in(3)
	E(ee334) = in(4)

	IN(1) = f1*E(ee44) - Dval(dd33,Di1) + Dval(dd33,Di2) -
     &    2*Ginv(4,4)*(Dval(dd00,Di1) - Dval(dd00,Di2))
	IN(2) = f2*E(ee44) - Dval(dd33,Di1) + Dval(dd33,Di3) -
     &    2*Ginv(4,4)*(Dval(dd00,Di1) - Dval(dd00,Di3))
	IN(3) = f3*E(ee44) - Dval(dd33,Di1) + Dval(dd33,Di4) -
     &    2*Ginv(4,4)*(Dval(dd00,Di1) - Dval(dd00,Di4))
	IN(4) = f4*E(ee44) - Dval(dd33,Di1) -
     &    2*Ginv(4,4)*(Dval(dd00,Di1) - Dval(dd00,Di5))
	call SOLVE(in)
	E(ee144) = in(1)
	E(ee244) = in(2)
	E(ee344) = in(3)
	E(ee444) = in(4)

	IN(1) = f1*E(ee12) + d1sum -
     &    2*Ginv(1,2)*(Dval(dd00,Di1) - Dval(dd00,Di2))
	IN(2) = f2*E(ee12) + d1sum -
     &    2*Ginv(1,2)*(Dval(dd00,Di1) - Dval(dd00,Di3))
	IN(3) = f3*E(ee12) + d1sum + Dval(dd12,Di4) -
     &    2*Ginv(1,2)*(Dval(dd00,Di1) - Dval(dd00,Di4))
	IN(4) = f4*E(ee12) + d1sum + Dval(dd12,Di5) -
     &    2*Ginv(1,2)*(Dval(dd00,Di1) - Dval(dd00,Di5))
	call SOLVE(in)
	E(ee112) = .5D0*(E(ee112) + in(1))
	E(ee122) = .5D0*(E(ee122) + in(2))
	E(ee123) = in(3)
	E(ee124) = in(4)

	IN(1) = f1*E(ee34) - Dval(dd23,Di1) + Dval(dd23,Di2) -
     &    2*Ginv(3,4)*(Dval(dd00,Di1) - Dval(dd00,Di2))
	IN(2) = f2*E(ee34) - Dval(dd23,Di1) + Dval(dd23,Di3) -
     &    2*Ginv(3,4)*(Dval(dd00,Di1) - Dval(dd00,Di3))
	IN(3) = f3*E(ee34) - Dval(dd23,Di1) -
     &    2*Ginv(3,4)*(Dval(dd00,Di1) - Dval(dd00,Di4))
	IN(4) = f4*E(ee34) - Dval(dd23,Di1) -
     &    2*Ginv(3,4)*(Dval(dd00,Di1) - Dval(dd00,Di5))
	call SOLVE(in)
	E(ee134) = in(1)
	E(ee234) = in(2)
	E(ee334) = .5D0*(E(ee334) + in(3))
	E(ee344) = .5D0*(E(ee344) + in(4))

	E(ee0000) = 0
	E(ee0011) = 0
	E(ee0012) = 0
	E(ee0013) = 0
	E(ee0014) = 0
	E(ee0022) = 0
	E(ee0023) = 0
	E(ee0024) = 0
	E(ee0033) = 0
	E(ee0034) = 0
	E(ee0044) = 0

	d0sum = d0sum + di + 2*dii + diii
	IN(1) = f1*E(ee111) + d0sum +
     &    6*Ginv(1,1)*d00sum
	IN(2) = f2*E(ee111) + d0sum + Dval(dd111,Di3) +
     &    6*Ginv(1,1)*(d00sum + Dval(dd001,Di3))
	IN(3) = f3*E(ee111) + d0sum + Dval(dd111,Di4) +
     &    6*Ginv(1,1)*(d00sum + Dval(dd001,Di4))
	IN(4) = f4*E(ee111) + d0sum + Dval(dd111,Di5) +
     &    6*Ginv(1,1)*(d00sum + Dval(dd001,Di5))
	call SOLVE(in)
	E(ee1111) = in(1)
	E(ee1112) = in(2)
	E(ee1113) = in(3)
	E(ee1114) = in(4)

	IN(1) = f1*E(ee222) - Dval(dd111,Di1) + Dval(dd111,Di2) -
     &    6*Ginv(2,2)*(Dval(dd001,Di1) - Dval(dd001,Di2))
	IN(2) = f2*E(ee222) - Dval(dd111,Di1) -
     &    6*Ginv(2,2)*Dval(dd001,Di1)
	IN(3) = f3*E(ee222) - Dval(dd111,Di1) + Dval(dd222,Di4) -
     &    6*Ginv(2,2)*(Dval(dd001,Di1) - Dval(dd002,Di4))
	IN(4) = f4*E(ee222) - Dval(dd111,Di1) + Dval(dd222,Di5) -
     &    6*Ginv(2,2)*(Dval(dd001,Di1) - Dval(dd002,Di5))
	call SOLVE(in)
	E(ee1222) = in(1)
	E(ee2222) = in(2)
	E(ee2223) = in(3)
	E(ee2224) = in(4)

	IN(1) = f1*E(ee333) - Dval(dd222,Di1) + Dval(dd222,Di2) -
     &    6*Ginv(3,3)*(Dval(dd002,Di1) - Dval(dd002,Di2))
	IN(2) = f2*E(ee333) - Dval(dd222,Di1) + Dval(dd222,Di3) -
     &    6*Ginv(3,3)*(Dval(dd002,Di1) - Dval(dd002,Di3))
	IN(3) = f3*E(ee333) - Dval(dd222,Di1) -
     &    6*Ginv(3,3)*Dval(dd002,Di1)
	IN(4) = f4*E(ee333) - Dval(dd222,Di1) + Dval(dd333,Di5) -
     &    6*Ginv(3,3)*(Dval(dd002,Di1) - Dval(dd003,Di5))
	call SOLVE(in)
	E(ee1333) = in(1)
	E(ee2333) = in(2)
	E(ee3333) = in(3)
	E(ee3334) = in(4)

	IN(1) = f1*E(ee444) - Dval(dd333,Di1) + Dval(dd333,Di2) -
     &    6*Ginv(4,4)*(Dval(dd003,Di1) - Dval(dd003,Di2))
	IN(2) = f2*E(ee444) - Dval(dd333,Di1) + Dval(dd333,Di3) -
     &    6*Ginv(4,4)*(Dval(dd003,Di1) - Dval(dd003,Di3))
	IN(3) = f3*E(ee444) - Dval(dd333,Di1) + Dval(dd333,Di4) -
     &    6*Ginv(4,4)*(Dval(dd003,Di1) - Dval(dd003,Di4))
	IN(4) = f4*E(ee444) - Dval(dd333,Di1) -
     &    6*Ginv(4,4)*Dval(dd003,Di1)
	call SOLVE(in)
	E(ee1444) = in(1)
	E(ee2444) = in(2)
	E(ee3444) = in(3)
	E(ee4444) = in(4)

	d1sum = d1sum + d1i + d1ii
	IN(1) = f1*E(ee112) - d1sum -
     &    2*Ginv(1,1)*(Dval(dd001,Di1) - Dval(dd001,Di2)) +
     &    4*Ginv(1,2)*d00sum
	IN(2) = f2*E(ee112) - d1sum -
     &    2*Ginv(1,1)*Dval(dd001,Di1) +
     &    4*Ginv(1,2)*(d00sum + Dval(dd001,Di3))
	IN(3) = f3*E(ee112) - d1sum + Dval(dd112,Di4) -
     &    2*Ginv(1,1)*(Dval(dd001,Di1) - Dval(dd002,Di4)) +
     &    4*Ginv(1,2)*(d00sum + Dval(dd001,Di4))
	IN(4) = f4*E(ee112) - d1sum + Dval(dd112,Di5) -
     &    2*Ginv(1,1)*(Dval(dd001,Di1) - Dval(dd002,Di5)) +
     &    4*Ginv(1,2)*(d00sum + Dval(dd001,Di5))
	call SOLVE(in)
	E(ee1112) = .5D0*(E(ee1112) + in(1))
	E(ee1122) = in(2)
	E(ee1123) = in(3)
	E(ee1124) = in(4)

	IN(1) = f1*E(ee223) - Dval(dd112,Di1) + Dval(dd112,Di2) -
     &    2*Ginv(2,2)*(Dval(dd002,Di1) - Dval(dd002,Di2)) -
     &    4*Ginv(2,3)*(Dval(dd001,Di1) - Dval(dd001,Di2))
	IN(2) = f2*E(ee223) - Dval(dd112,Di1) -
     &    2*Ginv(2,2)*(Dval(dd002,Di1) - Dval(dd002,Di3)) -
     &    4*Ginv(2,3)*Dval(dd001,Di1)
	IN(3) = f3*E(ee223) - Dval(dd112,Di1) -
     &    2*Ginv(2,2)*Dval(dd002,Di1) -
     &    4*Ginv(2,3)*(Dval(dd001,Di1) - Dval(dd002,Di4))
	IN(4) = f4*E(ee223) - Dval(dd112,Di1) + Dval(dd223,Di5) -
     &    2*Ginv(2,2)*(Dval(dd002,Di1) - Dval(dd003,Di5)) -
     &    4*Ginv(2,3)*(Dval(dd001,Di1) - Dval(dd002,Di5))
	call SOLVE(in)
	E(ee1223) = in(1)
	E(ee2223) = .5D0*(E(ee2223) + in(2))
	E(ee2233) = in(3)
	E(ee2234) = in(4)

	IN(1) = f1*E(ee334) - Dval(dd223,Di1) + Dval(dd223,Di2) -
     &    2*Ginv(3,3)*(Dval(dd003,Di1) - Dval(dd003,Di2)) -
     &    4*Ginv(3,4)*(Dval(dd002,Di1) - Dval(dd002,Di2))
	IN(2) = f2*E(ee334) - Dval(dd223,Di1) + Dval(dd223,Di3) -
     &    2*Ginv(3,3)*(Dval(dd003,Di1) - Dval(dd003,Di3)) -
     &    4*Ginv(3,4)*(Dval(dd002,Di1) - Dval(dd002,Di3))
	IN(3) = f3*E(ee334) - Dval(dd223,Di1) -
     &    2*Ginv(3,3)*(Dval(dd003,Di1) - Dval(dd003,Di4)) -
     &    4*Ginv(3,4)*Dval(dd002,Di1)
	IN(4) = f4*E(ee334) - Dval(dd223,Di1) -
     &    2*Ginv(3,3)*Dval(dd003,Di1) -
     &    4*Ginv(3,4)*(Dval(dd002,Di1) - Dval(dd003,Di5))
	call SOLVE(in)
	E(ee1334) = in(1)
	E(ee2334) = in(2)
	E(ee3334) = .5D0*(E(ee3334) + in(3))
	E(ee3344) = in(4)

	IN(1) = f1*E(ee144) + d33sum -
     &    4*Ginv(1,4)*(Dval(dd003,Di1) - Dval(dd003,Di2)) +
     &    2*Ginv(4,4)*d00sum
	IN(2) = f2*E(ee144) + d33sum + Dval(dd133,Di3) -
     &    4*Ginv(1,4)*(Dval(dd003,Di1) - Dval(dd003,Di3)) +
     &    2*Ginv(4,4)*(d00sum + Dval(dd001,Di3))
	IN(3) = f3*E(ee144) + d33sum + Dval(dd133,Di4) -
     &    4*Ginv(1,4)*(Dval(dd003,Di1) - Dval(dd003,Di4)) +
     &    2*Ginv(4,4)*(d00sum + Dval(dd001,Di4))
	IN(4) = f4*E(ee144) + d33sum -
     &    4*Ginv(1,4)*Dval(dd003,Di1) +
     &    2*Ginv(4,4)*(d00sum + Dval(dd001,Di5))
	call SOLVE(in)
	E(ee1144) = in(1)
	E(ee1244) = in(2)
	E(ee1344) = in(3)
	E(ee1444) = .5D0*(E(ee1444) + in(4))

	d2sum = d2sum + d2i + d2ii
	IN(1) = f1*E(ee113) - d2sum -
     &    2*Ginv(1,1)*(Dval(dd002,Di1) - Dval(dd002,Di2)) +
     &    4*Ginv(1,3)*d00sum
	IN(2) = f2*E(ee113) - d2sum + Dval(dd112,Di3) -
     &    2*Ginv(1,1)*(Dval(dd002,Di1) - Dval(dd002,Di3)) +
     &    4*Ginv(1,3)*(d00sum + Dval(dd001,Di3))
	IN(3) = f3*E(ee113) - d2sum -
     &    2*Ginv(1,1)*Dval(dd002,Di1) +
     &    4*Ginv(1,3)*(d00sum + Dval(dd001,Di4))
	IN(4) = f4*E(ee113) - d2sum + Dval(dd113,Di5) -
     &    2*Ginv(1,1)*(Dval(dd002,Di1) - Dval(dd003,Di5)) +
     &    4*Ginv(1,3)*(d00sum + Dval(dd001,Di5))
	call SOLVE(in)
	E(ee1113) = .5D0*(E(ee1113) + in(1))
	E(ee1123) = .5D0*(E(ee1123) + in(2))
	E(ee1133) = in(3)
	E(ee1134) = in(4)

	IN(1) = f1*E(ee224) - Dval(dd113,Di1) + Dval(dd113,Di2) -
     &    2*Ginv(2,2)*(Dval(dd003,Di1) - Dval(dd003,Di2)) -
     &    4*Ginv(2,4)*(Dval(dd001,Di1) - Dval(dd001,Di2))
	IN(2) = f2*E(ee224) - Dval(dd113,Di1) -
     &    2*Ginv(2,2)*(Dval(dd003,Di1) - Dval(dd003,Di3)) -
     &    4*Ginv(2,4)*Dval(dd001,Di1)
	IN(3) = f3*E(ee224) - Dval(dd113,Di1) + Dval(dd223,Di4) -
     &    2*Ginv(2,2)*(Dval(dd003,Di1) - Dval(dd003,Di4)) -
     &    4*Ginv(2,4)*(Dval(dd001,Di1) - Dval(dd002,Di4))
	IN(4) = f4*E(ee224) - Dval(dd113,Di1) -
     &    2*Ginv(2,2)*Dval(dd003,Di1) -
     &    4*Ginv(2,4)*(Dval(dd001,Di1) - Dval(dd002,Di5))
	call SOLVE(in)
	E(ee1224) = in(1)
	E(ee2224) = .5D0*(E(ee2224) + in(2))
	E(ee2234) = E(ee2234) + in(3)
	E(ee2244) = in(4)

	IN(1) = f1*E(ee234) - Dval(dd123,Di1) + Dval(dd123,Di2) -
     &    2*Ginv(2,3)*(Dval(dd003,Di1) - Dval(dd003,Di2)) -
     &    2*Ginv(3,4)*(Dval(dd001,Di1) - Dval(dd001,Di2)) -
     &    2*Ginv(2,4)*(Dval(dd002,Di1) - Dval(dd002,Di2))
	IN(2) = f2*E(ee234) - Dval(dd123,Di1) -
     &    2*Ginv(2,3)*(Dval(dd003,Di1) - Dval(dd003,Di3)) -
     &    2*Ginv(3,4)*Dval(dd001,Di1) -
     &    2*Ginv(2,4)*(Dval(dd002,Di1) - Dval(dd002,Di3))
	IN(3) = f3*E(ee234) - Dval(dd123,Di1) -
     &    2*Ginv(2,3)*(Dval(dd003,Di1) - Dval(dd003,Di4)) -
     &    2*Ginv(3,4)*(Dval(dd001,Di1) - Dval(dd002,Di4)) -
     &    2*Ginv(2,4)*Dval(dd002,Di1)
	IN(4) = f4*E(ee234) - Dval(dd123,Di1) -
     &    2*Ginv(2,3)*Dval(dd003,Di1) -
     &    2*Ginv(3,4)*(Dval(dd001,Di1) - Dval(dd002,Di5)) -
     &    2*Ginv(2,4)*(Dval(dd002,Di1) - Dval(dd003,Di5))
	call SOLVE(in)
	E(ee1234) = in(1)
	E(ee2234) = 1/3D0*(E(ee2234) + in(2))
	E(ee2334) = .5D0*(E(ee2334) + in(3))
	E(ee2344) = in(4)

	IN(1) = f1*E(ee133) + d22sum -
     &    4*Ginv(1,3)*(Dval(dd002,Di1) - Dval(dd002,Di2)) +
     &    2*Ginv(3,3)*d00sum
	IN(2) = f2*E(ee133) + d22sum + Dval(dd122,Di3) -
     &    4*Ginv(1,3)*(Dval(dd002,Di1) - Dval(dd002,Di3)) +
     &    2*Ginv(3,3)*(d00sum + Dval(dd001,Di3))
	IN(3) = f3*E(ee133) + d22sum -
     &    4*Ginv(1,3)*Dval(dd002,Di1) +
     &    2*Ginv(3,3)*(d00sum + Dval(dd001,Di4))
	IN(4) = f4*E(ee133) + d22sum + Dval(dd133,Di5) -
     &    4*Ginv(1,3)*(Dval(dd002,Di1) - Dval(dd003,Di5)) +
     &    2*Ginv(3,3)*(d00sum + Dval(dd001,Di5))
	call SOLVE(in)
	E(ee1133) = .5D0*(E(ee1133) + in(1))
	E(ee1233) = in(2)
	E(ee1333) = .5D0*(E(ee1333) + in(3))
	E(ee1334) = .5D0*(E(ee1334) + in(4))

	if( dump ) call XDumpCoeff(5, E)
	end

