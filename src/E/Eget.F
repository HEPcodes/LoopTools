* Eget.F
* retrieve the five-point tensor coefficients
* this file is part of LoopTools
* written by M. Rauch
* last modified 13 Sep 12 th

#include "externals.h"
#include "types.h"

#define npoint 5
#include "defs.h"


	memindex function XEget(p1, p2, p3, p4, p5, 
     &    p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5)
	implicit none
	DVAR p1, p2, p3, p4, p5
	DVAR p1p2, p2p3, p3p4, p4p5, p5p1
	DVAR m1, m2, m3, m4, m5

#include "lt.h"

	memindex cacheindex
	external cacheindex, XEcoeff

#ifdef COMPLEXPARA
	memindex Eget
	external Eget
#endif

	DVAR para(1,Pee)

	P(1) = p1
	P(2) = p2
	P(3) = p3
	P(4) = p4
	P(5) = p5
	P(6) = p1p2
	P(7) = p2p3
	P(8) = p3p4
	P(9) = p4p5
	P(10) = p5p1
#ifdef COMPLEXPARA
	if( abs(Im(P(1))) + abs(Im(P(2))) +
     &      abs(Im(P(3))) + abs(Im(P(4))) +
     &      abs(Im(P(5))) + abs(Im(P(6))) +
     &      abs(Im(P(7))) + abs(Im(P(8))) +
     &      abs(Im(P(9))) + abs(Im(P(10))) .gt. 0 )
     &    print *, "EgetC: Complex momenta not implemented"
#endif

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
	M(3) = m3
	if( abs(M(3)) .lt. minmass ) M(3) = 0
	M(4) = m4
	if( abs(M(4)) .lt. minmass ) M(4) = 0
	M(5) = m5
	if( abs(M(5)) .lt. minmass ) M(5) = 0
#ifdef COMPLEXPARA
	if( abs(Im(M(1))) + abs(Im(M(2))) +
     &      abs(Im(M(3))) + abs(Im(M(4))) + abs(Im(M(5))) .eq. 0 ) then
	  XEget = Eget(p1, p2, p3, p4, p1p2, p2p3,
     &      m1, m2, m3, m4) - offsetC
	  return
	endif
#endif

	XEget = cacheindex(para, Eval(1,0), XEcoeff, RC*Pee, Nee)
	end

************************************************************************

	subroutine XEput(res, p1, p2, p3, p4, p5, 
     &    p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5)
	implicit none
	ComplexType res(*)
	DVAR p1, p2, p3, p4, p5
	DVAR p1p2, p2p3, p3p4, p4p5, p5p1
	DVAR m1, m2, m3, m4, m5

#include "lt.h"

	external XEcoeff

	DVAR para(1,Pee)

	P(1) = p1
	P(2) = p2
	P(3) = p3
	P(4) = p4
	P(5) = p5
	P(6) = p1p2
	P(7) = p2p3
	P(8) = p3p4
	P(9) = p4p5
	P(10) = p5p1
#ifdef COMPLEXPARA
	if( abs(Im(P(1))) + abs(Im(P(2))) +
     &      abs(Im(P(3))) + abs(Im(P(4))) +
     &      abs(Im(P(5))) + abs(Im(P(6))) +
     &      abs(Im(P(7))) + abs(Im(P(8))) +
     &      abs(Im(P(9))) + abs(Im(P(10))) .gt. 0 )
     &    print *, "EgetC: Complex momenta not implemented"
#endif

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
	M(3) = m3
	if( abs(M(3)) .lt. minmass ) M(3) = 0
	M(4) = m4
	if( abs(M(4)) .lt. minmass ) M(4) = 0
	M(5) = m5
	if( abs(M(5)) .lt. minmass ) M(5) = 0
#ifdef COMPLEXPARA
	if( abs(Im(M(1))) + abs(Im(M(2))) +
     &      abs(Im(M(3))) + abs(Im(M(4))) + abs(Im(M(5))) .eq. 0 ) then
	  call Eput(res, p1, p2, p3, p4, p1p2, p2p3,
     &      m1, m2, m3, m4)
	  return
	endif
#endif

	call cachecopy(res, para, Eval(1,0), XEcoeff, RC*Pee, Nee)
	end

************************************************************************

	ComplexType function XE0i(i, p1, p2, p3, p4, p5,
     &    p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5)
	implicit none
	integer i
	DVAR p1, p2, p3, p4, p5
	DVAR p1p2, p2p3, p3p4, p4p5, p5p1
	DVAR m1, m2, m3, m4, m5

#include "lt.h"

	memindex XEget
	external XEget

	memindex b

	b = XEget(p1, p2, p3, p4, p5,
     &    p1p2, p2p3, p3p4, p4p5, p5p1, m1, m2, m3, m4, m5)
	XE0i = Eval(i,b)
	end

************************************************************************

	subroutine XEcoeff(res, para)
	implicit none
	ComplexType res(*)
	DVAR para(1,*)

#include "lt.h"

	ComplexType cmp(Nee)

#ifdef COMPLEXPARA
	goto (1, 2, 3) ibits(versionkey, KeyEgetC, 2)
#else
	goto (1, 2, 3) ibits(versionkey, KeyEget, 2)
#endif

	call XEcoeffa(res, para)
	return

1	call XEcoeffb(res, para)
	return

2	call XEcheck(res, cmp, para)
	return

3	call XEcheck(cmp, res, para)
	end

************************************************************************

	subroutine XEcheck(Ea, Eb, para)
	implicit none
	ComplexType Ea(*), Eb(*)
	DVAR para(1,*)

#include "lt.h"

	ComplexType dE(Nee)
	integer i
	logical ini

	QVAR Ginv(4,4)
	common /XInvGramE/ Ginv

	character*8 coeffname(Nee,2:5)
	common /ltcoeffnames/ coeffname

	call XEcoeffa(Ea, para)
	call XEcoeffb(Eb, para)

	dE(ee0) = 0
	dE(ee1) = 0
	dE(ee2) = 0
	dE(ee3) = 0
	dE(ee4) = 0
	dE(ee11) = -2*Ea(ee00)*Ginv(1,1)
	dE(ee12) = -2*Ea(ee00)*Ginv(1,2)
	dE(ee13) = -2*Ea(ee00)*Ginv(1,3)
	dE(ee14) = -2*Ea(ee00)*Ginv(1,4)
	dE(ee22) = -2*Ea(ee00)*Ginv(2,2)
	dE(ee23) = -2*Ea(ee00)*Ginv(2,3)
	dE(ee24) = -2*Ea(ee00)*Ginv(2,4)
	dE(ee33) = -2*Ea(ee00)*Ginv(3,3)
	dE(ee34) = -2*Ea(ee00)*Ginv(3,4)
	dE(ee44) = -2*Ea(ee00)*Ginv(4,4)
	dE(ee00) = Ea(ee00)
	dE(ee111) = -6*Ea(ee001)*Ginv(1,1)
	dE(ee112) = -2*(Ea(ee002)*Ginv(1,1) +
     &    Ea(ee001)*(Ginv(1,2) + Ginv(2,1)))
	dE(ee113) = -2*(Ea(ee003)*Ginv(1,1) +
     &    Ea(ee001)*(Ginv(1,3) + Ginv(3,1)))
	dE(ee114) = -2*(Ea(ee004)*Ginv(1,1) +
     &    Ea(ee001)*(Ginv(1,4) + Ginv(4,1)))
	dE(ee122) = -2*(Ea(ee002)*(Ginv(1,2) + Ginv(2,1)) +
     &    Ea(ee001)*Ginv(2,2))
        dE(ee123) = -2*(Ea(ee003)*Ginv(1,2) + Ea(ee001)*Ginv(2,3) + 
     &    Ea(ee002)*Ginv(3,1))
        dE(ee124) = -2*(Ea(ee004)*Ginv(1,2) + Ea(ee001)*Ginv(2,4) + 
     &    Ea(ee002)*Ginv(4,1))
	dE(ee133) = -2*(Ea(ee003)*(Ginv(1,3) + Ginv(3,1)) +
     &    Ea(ee001)*Ginv(3,3))
        dE(ee134) = -2*(Ea(ee004)*Ginv(1,3) + Ea(ee001)*Ginv(3,4) + 
     &    Ea(ee003)*Ginv(4,1))
	dE(ee144) = -2*(Ea(ee004)*(Ginv(1,4) + Ginv(4,1)) +
     &    Ea(ee001)*Ginv(4,4))
	dE(ee222) = -6*Ea(ee002)*Ginv(2,2)
	dE(ee223) = -2*(Ea(ee003)*Ginv(2,2) +
     &    Ea(ee002)*(Ginv(2,3) + Ginv(3,2)))
	dE(ee224) = -2*(Ea(ee004)*Ginv(2,2) +
     &    Ea(ee002)*(Ginv(2,4) + Ginv(4,2)))
	dE(ee233) = -2*(Ea(ee003)*(Ginv(2,3) + Ginv(3,2)) +
     &    Ea(ee002)*Ginv(3,3))
        dE(ee234) = -2*(Ea(ee004)*Ginv(2,3) + Ea(ee002)*Ginv(3,4) + 
     &    Ea(ee003)*Ginv(4,2))
	dE(ee244) = -2*(Ea(ee004)*(Ginv(2,4) + Ginv(4,2)) +
     &    Ea(ee002)*Ginv(4,4))
	dE(ee333) = -6*Ea(ee003)*Ginv(3,3)
	dE(ee334) = -2*(Ea(ee004)*Ginv(3,3) +
     &    Ea(ee003)*(Ginv(3,4) + Ginv(4,3)))
	dE(ee344) = -2*(Ea(ee004)*(Ginv(3,4) + Ginv(4,3)) +
     &    Ea(ee003)*Ginv(4,4))
	dE(ee444) = -6*Ea(ee004)*Ginv(4,4)
	dE(ee001) = Ea(ee001)
	dE(ee002) = Ea(ee002)
	dE(ee003) = Ea(ee003)
	dE(ee004) = Ea(ee004)
	dE(ee1111) = -12*Ginv(1,1)*(Ea(ee0011) + Ea(ee0000)*Ginv(1,1))
        dE(ee1112) = -6*(Ea(ee0012)*Ginv(1,1) + 
     &     (Ea(ee0011) + 2*Ea(ee0000)*Ginv(1,1))*Ginv(1,2))
        dE(ee1113) = -6*(Ea(ee0013)*Ginv(1,1) + 
     &     (Ea(ee0011) + 2*Ea(ee0000)*Ginv(1,1))*Ginv(1,3))
        dE(ee1114) = -6*(Ea(ee0014)*Ginv(1,1) + 
     &     (Ea(ee0011) + 2*Ea(ee0000)*Ginv(1,1))*Ginv(1,4))
        dE(ee1122) = -2*(Ea(ee0022)*Ginv(1,1) + 
     &     4*Ginv(1,2)*(Ea(ee0012) + Ea(ee0000)*Ginv(1,2)) + 
     &     (Ea(ee0011) + 2*Ea(ee0000)*Ginv(1,1))*Ginv(2,2))
        dE(ee1123) = -2*(Ea(ee0023)*Ginv(1,1) + 
     &     2*Ea(ee0013)*Ginv(1,2) + 
     &     2*(Ea(ee0012) + 2*Ea(ee0000)*Ginv(1,2))*Ginv(1,3) + 
     &     (Ea(ee0011) + 2*Ea(ee0000)*Ginv(1,1))*Ginv(2,3))
        dE(ee1124) = -2*(Ea(ee0024)*Ginv(1,1) + 
     &     2*Ea(ee0014)*Ginv(1,2) + 
     &     2*(Ea(ee0012) + 2*Ea(ee0000)*Ginv(1,2))*Ginv(1,4) + 
     &     (Ea(ee0011) + 2*Ea(ee0000)*Ginv(1,1))*Ginv(2,4))
        dE(ee1133) = -2*(Ea(ee0033)*Ginv(1,1) + 
     &     4*Ginv(1,3)*(Ea(ee0013) + Ea(ee0000)*Ginv(1,3)) + 
     &     (Ea(ee0011) + 2*Ea(ee0000)*Ginv(1,1))*Ginv(3,3))
        dE(ee1134) = -2*(Ea(ee0034)*Ginv(1,1) + 
     &     2*Ea(ee0014)*Ginv(1,3) + 
     &     2*(Ea(ee0013) + 2*Ea(ee0000)*Ginv(1,3))*Ginv(1,4) + 
     &     (Ea(ee0011) + 2*Ea(ee0000)*Ginv(1,1))*Ginv(3,4))
        dE(ee1144) = -2*(Ea(ee0044)*Ginv(1,1) + 
     &     4*Ginv(1,4)*(Ea(ee0014) + Ea(ee0000)*Ginv(1,4)) + 
     &     (Ea(ee0011) + 2*Ea(ee0000)*Ginv(1,1))*Ginv(4,4))
        dE(ee1222) = -6*(Ea(ee0022)*Ginv(1,2) + 
     &     (Ea(ee0012) + 2*Ea(ee0000)*Ginv(1,2))*Ginv(2,2))
        dE(ee1223) = -2*(2*Ea(ee0023)*Ginv(1,2) + 
     &     Ea(ee0022)*Ginv(1,3) + 
     &     (Ea(ee0013) + 2*Ea(ee0000)*Ginv(1,3))*Ginv(2,2) + 
     &     2*(Ea(ee0012) + 2*Ea(ee0000)*Ginv(1,2))*Ginv(2,3))
        dE(ee1224) = -2*(2*Ea(ee0024)*Ginv(1,2) + 
     &     Ea(ee0022)*Ginv(1,4) + 
     &     (Ea(ee0014) + 2*Ea(ee0000)*Ginv(1,4))*Ginv(2,2) + 
     &     2*(Ea(ee0012) + 2*Ea(ee0000)*Ginv(1,2))*Ginv(2,4))
        dE(ee1233) = -2*(Ea(ee0033)*Ginv(1,2) + 
     &     2*Ea(ee0023)*Ginv(1,3) + 
     &     2*(Ea(ee0013) + 2*Ea(ee0000)*Ginv(1,3))*Ginv(2,3) + 
     &     (Ea(ee0012) + 2*Ea(ee0000)*Ginv(1,2))*Ginv(3,3))
        dE(ee1234) = -2*(Ea(ee0023)*Ginv(1,4) + 
     &     (Ea(ee0014) + 2*Ea(ee0000)*Ginv(1,4))*Ginv(2,3) + 
     &     Ea(ee0013)*Ginv(2,4) + 
     &     Ginv(1,3)*(Ea(ee0024) + 2*Ea(ee0000)*Ginv(2,4)) + 
     &     Ea(ee0012)*Ginv(3,4) + 
     &     Ginv(1,2)*(Ea(ee0034) + 2*Ea(ee0000)*Ginv(3,4)))
        dE(ee1244) = -2*(Ea(ee0044)*Ginv(1,2) + 
     &     2*Ea(ee0024)*Ginv(1,4) + 
     &     2*(Ea(ee0014) + 2*Ea(ee0000)*Ginv(1,4))*Ginv(2,4) + 
     &     (Ea(ee0012) + 2*Ea(ee0000)*Ginv(1,2))*Ginv(4,4))
        dE(ee1333) = -6*(Ea(ee0033)*Ginv(1,3) + 
     &     (Ea(ee0013) + 2*Ea(ee0000)*Ginv(1,3))*Ginv(3,3))
        dE(ee1334) = -2*(2*Ea(ee0034)*Ginv(1,3) + 
     &     Ea(ee0033)*Ginv(1,4) + 
     &     (Ea(ee0014) + 2*Ea(ee0000)*Ginv(1,4))*Ginv(3,3) + 
     &     2*(Ea(ee0013) + 2*Ea(ee0000)*Ginv(1,3))*Ginv(3,4))
        dE(ee1344) = -2*(Ea(ee0044)*Ginv(1,3) + 
     &     2*Ea(ee0034)*Ginv(1,4) + 
     &     2*(Ea(ee0014) + 2*Ea(ee0000)*Ginv(1,4))*Ginv(3,4) + 
     &     (Ea(ee0013) + 2*Ea(ee0000)*Ginv(1,3))*Ginv(4,4))
        dE(ee1444) = -6*(Ea(ee0044)*Ginv(1,4) + 
     &     (Ea(ee0014) + 2*Ea(ee0000)*Ginv(1,4))*Ginv(4,4))
	dE(ee2222) = -12*Ginv(2,2)*(Ea(ee0022) + Ea(ee0000)*Ginv(2,2))
        dE(ee2223) = -6*(Ea(ee0023)*Ginv(2,2) + 
     &     (Ea(ee0022) + 2*Ea(ee0000)*Ginv(2,2))*Ginv(2,3))
        dE(ee2224) = -6*(Ea(ee0024)*Ginv(2,2) + 
     &     (Ea(ee0022) + 2*Ea(ee0000)*Ginv(2,2))*Ginv(2,4))
        dE(ee2233) = -2*(Ea(ee0033)*Ginv(2,2) + 
     &     4*Ginv(2,3)*(Ea(ee0023) + Ea(ee0000)*Ginv(2,3)) + 
     &     (Ea(ee0022) + 2*Ea(ee0000)*Ginv(2,2))*Ginv(3,3))
        dE(ee2234) = -2*(Ea(ee0034)*Ginv(2,2) + 
     &     2*Ea(ee0024)*Ginv(2,3) + 
     &     2*(Ea(ee0023) + 2*Ea(ee0000)*Ginv(2,3))*Ginv(2,4) + 
     &     (Ea(ee0022) + 2*Ea(ee0000)*Ginv(2,2))*Ginv(3,4))
        dE(ee2244) = -2*(Ea(ee0044)*Ginv(2,2) + 
     &     4*Ginv(2,4)*(Ea(ee0024) + Ea(ee0000)*Ginv(2,4)) + 
     &     (Ea(ee0022) + 2*Ea(ee0000)*Ginv(2,2))*Ginv(4,4))
        dE(ee2333) = -6*(Ea(ee0033)*Ginv(2,3) + 
     &     (Ea(ee0023) + 2*Ea(ee0000)*Ginv(2,3))*Ginv(3,3))
        dE(ee2334) = -2*(2*Ea(ee0034)*Ginv(2,3) + 
     &     Ea(ee0033)*Ginv(2,4) + 
     &     (Ea(ee0024) + 2*Ea(ee0000)*Ginv(2,4))*Ginv(3,3) + 
     &     2*(Ea(ee0023) + 2*Ea(ee0000)*Ginv(2,3))*Ginv(3,4))
        dE(ee2344) = -2*(Ea(ee0044)*Ginv(2,3) + 
     &     2*Ea(ee0034)*Ginv(2,4) + 
     &     2*(Ea(ee0024) + 2*Ea(ee0000)*Ginv(2,4))*Ginv(3,4) + 
     &     (Ea(ee0023) + 2*Ea(ee0000)*Ginv(2,3))*Ginv(4,4))
        dE(ee2444) = -6*(Ea(ee0044)*Ginv(2,4) + 
     &     (Ea(ee0024) + 2*Ea(ee0000)*Ginv(2,4))*Ginv(4,4))
	dE(ee3333) = -12*Ginv(3,3)*(Ea(ee0033) + Ea(ee0000)*Ginv(3,3))
        dE(ee3334) = -6*(Ea(ee0034)*Ginv(3,3) + 
     &     (Ea(ee0033) + 2*Ea(ee0000)*Ginv(3,3))*Ginv(3,4))
        dE(ee3344) = -2*(Ea(ee0044)*Ginv(3,3) + 
     &     4*Ginv(3,4)*(Ea(ee0034) + Ea(ee0000)*Ginv(3,4)) + 
     &     (Ea(ee0033) + 2*Ea(ee0000)*Ginv(3,3))*Ginv(4,4))
        dE(ee3444) = -6*(Ea(ee0044)*Ginv(3,4) + 
     &     (Ea(ee0034) + 2*Ea(ee0000)*Ginv(3,4))*Ginv(4,4))
	dE(ee4444) = -12*Ginv(4,4)*(Ea(ee0044) + Ea(ee0000)*Ginv(4,4))
	dE(ee0000) = Ea(ee0000)
	dE(ee0011) = Ea(ee0011)
	dE(ee0012) = Ea(ee0012)
	dE(ee0013) = Ea(ee0013)
	dE(ee0014) = Ea(ee0014)
	dE(ee0022) = Ea(ee0022)
	dE(ee0023) = Ea(ee0023)
	dE(ee0024) = Ea(ee0024)
	dE(ee0033) = Ea(ee0033)
	dE(ee0034) = Ea(ee0034)
	dE(ee0044) = Ea(ee0044)

	ini = .TRUE.
	do i = 1, Nee
	  if( abs(Ea(i) - Eb(i) - dE(i)) .gt.
     &        .5D0*maxdev*abs(Ea(i) + Eb(i)) ) then
	    if( ini ) then
#ifdef COMPLEXPARA
	      print *, "Discrepancy in CEget:"
#else
	      print *, "Discrepancy in Eget:"
#endif
	      call XDumpPara(5, para, " ")
	      ini = .FALSE.
	    endif
	    print *, coeffname(i,5), " a =", Ea(i)
	    print *, coeffname(i,5), " b =", Eb(i) + dE(i)
	  endif
	enddo
	end

