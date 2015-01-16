* Ecoeffa.F
* the five-point tensor coefficients
* this file is part of LoopTools
* written by M. Rauch
* last modified 7 Oct 12 th

#include "externals.h"
#include "types.h"

#define npoint 5
#include "defs.h"


	subroutine XEcoeffa(E, para)
	implicit none
	ComplexType E(*)
	DVAR para(1,*)

#include "lt.h"

	memindex XDget
	external XDget

	DVAR p1, p2, p3, p4, p5
	DVAR p1p2, p2p3, p3p4, p4p5, p5p1
	DVAR m1, m2, m3, m4, m5
	QVAR Y(5,5), Yi(5,5), Z(4,4), Zij(3,3)
	QVAR eta(5), zeta(4,4), detY, detZ
	RealType del, del4
	integer i, j, k, l, finite
	memindex D2345, D1345, D1245, D1235, D1234
	logical dump

	QVAR Yflat(25), Yiflat(25)
	equivalence (Y, Yflat)
	equivalence (Yi, Yiflat)

	ComplexType help1(0:2), help2(0:2), help3(0:2), help4(0:2)
	ComplexType dabbr41(0:2), dabbr48(0:2), dabbr65(0:2)
	ComplexType dabbr60(0:2), dabbr55(0:2), dabbr50(0:2)
	ComplexType dabbr49(0:2), dabbr45(0:2), dabbr42(0:2)
	ComplexType dabbr10(0:2), dabbr84(0:2), dabbr91(0:2)
	ComplexType dabbr81(0:2), dabbr52(0:2), dabbr88(0:2)
	ComplexType dabbr77(0:2), dabbr90(0:2), dabbr46(0:2)
	ComplexType dabbr87(0:2), dabbr74(0:2), dabbr80(0:2)
	ComplexType dabbr71(0:2), dabbr83(0:2), dabbr76(0:2)
	ComplexType dabbr79(0:2), dabbr43(0:2), dabbr70(0:2)
	ComplexType dabbr73(0:2), dabbr35(0:2), dabbr37(0:2)
	ComplexType dabbr39(0:2), dabbr27(0:2), dabbr31(0:2)
	ComplexType dabbr24(0:2), dabbr20(0:2), dabbr16(0:2)
	ComplexType dabbr30(0:2), dabbr13(0:2), dabbr23(0:2)
	ComplexType dabbr26(0:2), dabbr5(0:2), dabbr7(0:2)
	ComplexType dabbr9(0:2), dabbr33(0:2), dabbr92(0:2)
	ComplexType dabbr89(0:2), dabbr86(0:2), dabbr85(0:2)
	ComplexType dabbr82(0:2), dabbr78(0:2), dabbr75(0:2)
	ComplexType dabbr72(0:2), dabbr69(0:2), dabbr68(0:2)
	ComplexType dabbr38(0:2), dabbr36(0:2), dabbr34(0:2)
	ComplexType dabbr32(0:2), dabbr29(0:2), dabbr28(0:2)
	ComplexType dabbr25(0:2), dabbr22(0:2), dabbr21(0:2)
	ComplexType dabbr3(0:2), dabbr8(0:2), dabbr6(0:2)
	ComplexType dabbr4(0:2), dabbr1(0:2), dabbr61(0:2)
	ComplexType dabbr66(0:2), dabbr57(0:2), dabbr47(0:2)
	ComplexType dabbr51(0:2), dabbr64(0:2), dabbr44(0:2)
	ComplexType dabbr56(0:2), dabbr59(0:2), dabbr14(0:2)
	ComplexType dabbr17(0:2), dabbr19(0:2), dabbr67(0:2)
	ComplexType dabbr63(0:2), dabbr62(0:2), dabbr58(0:2)
	ComplexType dabbr54(0:2), dabbr53(0:2), dabbr11(0:2)
	ComplexType dabbr18(0:2), dabbr15(0:2), dabbr12(0:2)
	ComplexType dabbr2(0:2), dabbr40(0:2)

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

	finite = 0
	if( lambda .ge. 0 ) then
	  finite = 1
	  del = (delta + log(mudim))/24D0
	else if( lambda .eq. -1 ) then
	  del = 1/24D0
	endif
	del4 = .25D0*del

	D2345 = XDget(p2, p3, p4, p5p1, p2p3, p3p4, m2, m3, m4, m5)
	D1345 = XDget(p1p2, p3, p4, p5, p4p5, p3p4, m1, m3, m4, m5)
	D1245 = XDget(p1, p2p3, p4, p5, p4p5, p5p1, m1, m2, m4, m5)
	D1235 = XDget(p1, p2, p3p4, p5, p1p2, p5p1, m1, m2, m3, m5)
	D1234 = XDget(p1, p2, p3, p4p5, p1p2, p2p3, m1, m2, m3, m4)

	serial = serial + 1
	dump = ibits(debugkey, DebugE, 1) .ne. 0 .and.
     &    serial .ge. debugfrom .and. serial .le. debugto
	if( dump ) call XDumpPara(5, para, "Ecoeffa")

	Y(1,1) = 2*m1
	Y(2,2) = 2*m2
	Y(3,3) = 2*m3
	Y(4,4) = 2*m4
	Y(5,5) = 2*m5
	Y(1,2) = m1
	  Y(1,2) = Y(1,2) + m2
	  Y(1,2) = Y(1,2) - p1
	Y(2,1) = Y(1,2)
	Y(1,3) = m1
	  Y(1,3) = Y(1,3) + m3
	  Y(1,3) = Y(1,3) - p1p2
	Y(3,1) = Y(1,3)
	Y(1,4) = m1
	  Y(1,4) = Y(1,4) + m4
	  Y(1,4) = Y(1,4) - p4p5
	Y(4,1) = Y(1,4)
	Y(1,5) = m1
	  Y(1,5) = Y(1,5) + m5
	  Y(1,5) = Y(1,5) - p5
	Y(5,1) = Y(1,5)
	Y(2,3) = m2
	  Y(2,3) = Y(2,3) + m3
	  Y(2,3) = Y(2,3) - p2
	Y(3,2) = Y(2,3)
	Y(2,4) = m2
	  Y(2,4) = Y(2,4) + m4
	  Y(2,4) = Y(2,4) - p2p3
	Y(4,2) = Y(2,4)
	Y(2,5) = m2
	  Y(2,5) = Y(2,5) + m5
	  Y(2,5) = Y(2,5) - p5p1
	Y(5,2) = Y(2,5)
	Y(3,4) = m3
	  Y(3,4) = Y(3,4) + m4
	  Y(3,4) = Y(3,4) - p3
	Y(4,3) = Y(3,4)
	Y(3,5) = m3
	  Y(3,5) = Y(3,5) + m5
	  Y(3,5) = Y(3,5) - p3p4
	Y(5,3) = Y(3,5)
	Y(4,5) = m4
	  Y(4,5) = Y(4,5) + m5
	  Y(4,5) = Y(4,5) - p4
	Y(5,4) = Y(4,5)

* calculate the Y(i), their determinants, and eta(i)
	do i = 1, 5
	  do j = 1, 25
	    Yiflat(j) = Yflat(j)
	  enddo
	  do j = 1, 5
	    Yi(j,i) = 1
	  enddo
	  call XDet(5, Yi,5, eta(i))
	enddo

* Y is no longer needed, now calculate its determinant and
* add the missing factor 1/detY to eta
	call XDet(5, Y,5, detY)
	do i = 1, 5
	  eta(i) = eta(i)/detY
	enddo

	Z(1,1) = 2*p1
	Z(2,2) = 2*p1p2
	Z(3,3) = 2*p4p5
	Z(4,4) = 2*p5
	Z(1,2) = p1
	  Z(1,2) = Z(1,2) + p1p2
	  Z(1,2) = Z(1,2) - p2
	Z(2,1) = Z(1,2)
	Z(1,3) = p1
	  Z(1,3) = Z(1,3) - p2p3
	  Z(1,3) = Z(1,3) + p4p5
	Z(3,1) = Z(1,3)
	Z(1,4) = p1
	  Z(1,4) = Z(1,4) - p5p1
	  Z(1,4) = Z(1,4) + p5
	Z(4,1) = Z(1,4)
	Z(2,3) = p1p2
	  Z(2,3) = Z(2,3) - p3
	  Z(2,3) = Z(2,3) + p4p5
	Z(3,2) = Z(2,3)
	Z(2,4) = p1p2
	  Z(2,4) = Z(2,4) - p3p4
	  Z(2,4) = Z(2,4) + p5
	Z(4,2) = Z(2,4)
	Z(3,4) = p5
	  Z(3,4) = Z(3,4) + p4p5
	  Z(3,4) = Z(3,4) - p4
	Z(4,3) = Z(3,4)

* calculate the zeta(i,j)
	do i = 1, 4
	  do j = i, 4
* generate the submatrix Z_ij
	    do l = 1, 3
	      do k = 1, 3
	        Zij(k,l) = Z(k + ibits(not(k - i), 3, 1),
     &                       l + ibits(not(l - j), 3, 1))
	      enddo
	    enddo
	    call XDet(3, Zij,3, detZ)
	    zeta(i,j) = Sgn(i + j)*detZ/detY
	    zeta(j,i) = zeta(i,j)
	  enddo
	enddo
	call XDet(4, Z,4, detZ)

	help1 = Dval(DD003,D2345) + 2*Dval(DD0033,D2345) + 
     &   Dval(DD00333,D2345)
	help2 = Dval(DD002,D2345) + 2*Dval(DD0023,D2345) + 
     &   Dval(DD00233,D2345)
	help3 = Dval(DD001,D2345) + 2*Dval(DD0013,D2345) + 
     &   Dval(DD00133,D2345)
	help4 = Dval(DD00,D2345) + 2*Dval(DD003,D2345) + 
     &   Dval(DD0033,D2345)
	dabbr41 = Dval(DD00223,D2345) + Dval(DD00233,D2345)
	dabbr48 = Dval(DD00113,D2345) + Dval(DD00133,D2345)
	dabbr65 = Dval(DD0023,D2345) + Dval(DD00123,D2345)
	dabbr60 = Dval(DD0013,D2345) + Dval(DD00123,D2345)
	dabbr55 = Dval(DD0012,D2345) + Dval(DD00123,D2345)
	dabbr50 = Dval(DD00112,D2345) + Dval(DD00122,D2345)
	dabbr49 = Dval(DD2233,D2345) + Dval(DD2333,D2345)
	dabbr45 = Dval(DD2223,D2345) + Dval(DD2233,D2345)
	dabbr42 = Dval(DD1223,D2345) + Dval(DD1233,D2345)
	dabbr10 = Dval(DD223,D2345) + Dval(DD233,D2345)
	dabbr84 = Dval(DD0033,D2345) + Dval(DD00233,D2345) + 
     &   Dval(DD00333,D2345)
	dabbr91 = Dval(DD0033,D2345) + Dval(DD00133,D2345) + 
     &   Dval(DD00333,D2345)
	dabbr81 = dabbr41 + Dval(DD0023,D2345)
	dabbr52 = Dval(DD0033,D2345) + Dval(DD00133,D2345) + 
     &   Dval(DD00233,D2345)
	dabbr88 = dabbr65 + Dval(DD00233,D2345)
	dabbr77 = Dval(DD0022,D2345) + Dval(DD00222,D2345) + 
     &   Dval(DD00223,D2345)
	dabbr90 = dabbr65 + Dval(DD00223,D2345)
	dabbr46 = Dval(DD0022,D2345) + Dval(DD00122,D2345) + 
     &   Dval(DD00223,D2345)
	dabbr87 = Dval(DD0022,D2345) + Dval(DD00122,D2345) + 
     &   Dval(DD00222,D2345)
	dabbr74 = dabbr60 + Dval(DD00133,D2345)
	dabbr80 = dabbr48 + Dval(DD0013,D2345)
	dabbr71 = dabbr55 + Dval(DD00122,D2345)
	dabbr83 = dabbr60 + Dval(DD00113,D2345)
	dabbr76 = dabbr55 + Dval(DD00112,D2345)
	dabbr79 = dabbr50 + Dval(DD0012,D2345)
	dabbr43 = Dval(DD0011,D2345) + Dval(DD00112,D2345) + 
     &   Dval(DD00113,D2345)
	dabbr70 = Dval(DD0011,D2345) + Dval(DD00111,D2345) + 
     &   Dval(DD00113,D2345)
	dabbr73 = Dval(DD0011,D2345) + Dval(DD00111,D2345) + 
     &   Dval(DD00112,D2345)
	dabbr35 = (Dval(DD0000,D2345) - del) + 
     &   (Dval(DD00002,D2345) + del4) + 
     &   (Dval(DD00003,D2345) + del4)
	dabbr37 = (Dval(DD0000,D2345) - del) + 
     &   (Dval(DD00001,D2345) + del4) + 
     &   (Dval(DD00003,D2345) + del4)
	dabbr39 = (Dval(DD0000,D2345) - del) + 
     &   (Dval(DD00001,D2345) + del4) + 
     &   (Dval(DD00002,D2345) + del4)
	dabbr27 = Dval(DD003,D2345) + Dval(DD0023,D2345) + 
     &   Dval(DD0033,D2345)
	dabbr31 = Dval(DD003,D2345) + Dval(DD0013,D2345) + 
     &   Dval(DD0033,D2345)
	dabbr24 = Dval(DD002,D2345) + Dval(DD0022,D2345) + 
     &   Dval(DD0023,D2345)
	dabbr20 = Dval(DD003,D2345) + Dval(DD0013,D2345) + 
     &   Dval(DD0023,D2345)
	dabbr16 = Dval(DD002,D2345) + Dval(DD0012,D2345) + 
     &   Dval(DD0023,D2345)
	dabbr30 = Dval(DD002,D2345) + Dval(DD0012,D2345) + 
     &   Dval(DD0022,D2345)
	dabbr13 = Dval(DD001,D2345) + Dval(DD0012,D2345) + 
     &   Dval(DD0013,D2345)
	dabbr23 = Dval(DD001,D2345) + Dval(DD0011,D2345) + 
     &   Dval(DD0013,D2345)
	dabbr26 = Dval(DD001,D2345) + Dval(DD0011,D2345) + 
     &   Dval(DD0012,D2345)
	dabbr5 = Dval(DD00,D2345) + Dval(DD002,D2345) + 
     &   Dval(DD003,D2345)
	dabbr7 = Dval(DD00,D2345) + Dval(DD001,D2345) + 
     &   Dval(DD003,D2345)
	dabbr9 = Dval(DD00,D2345) + Dval(DD001,D2345) + 
     &   Dval(DD002,D2345)
	dabbr33 = dabbr35 + (Dval(DD00001,D2345) + del4)
	dabbr92 = Dval(DD333,D2345) + Dval(DD1333,D2345) + 
     &   Dval(DD2333,D2345) + Dval(DD3333,D2345)
	dabbr89 = dabbr49 + Dval(DD233,D2345) + Dval(DD1233,D2345)
	dabbr86 = dabbr45 + Dval(DD223,D2345) + Dval(DD1223,D2345)
	dabbr85 = Dval(DD222,D2345) + Dval(DD1222,D2345) + 
     &   Dval(DD2222,D2345) + Dval(DD2223,D2345)
	dabbr82 = Dval(DD133,D2345) + Dval(DD1133,D2345) + 
     &   Dval(DD1233,D2345) + Dval(DD1333,D2345)
	dabbr78 = dabbr42 + Dval(DD123,D2345) + Dval(DD1123,D2345)
	dabbr75 = Dval(DD122,D2345) + Dval(DD1122,D2345) + 
     &   Dval(DD1222,D2345) + Dval(DD1223,D2345)
	dabbr72 = Dval(DD113,D2345) + Dval(DD1113,D2345) + 
     &   Dval(DD1123,D2345) + Dval(DD1133,D2345)
	dabbr69 = Dval(DD112,D2345) + Dval(DD1112,D2345) + 
     &   Dval(DD1122,D2345) + Dval(DD1123,D2345)
	dabbr68 = Dval(DD111,D2345) + Dval(DD1111,D2345) + 
     &   Dval(DD1112,D2345) + Dval(DD1113,D2345)
	dabbr38 = dabbr27 + Dval(DD0013,D2345)
	dabbr36 = dabbr24 + Dval(DD0012,D2345)
	dabbr34 = dabbr13 + Dval(DD0011,D2345)
	dabbr32 = Dval(DD33,D2345) + Dval(DD133,D2345) + 
     &   Dval(DD233,D2345) + Dval(DD333,D2345)
	dabbr29 = dabbr10 + Dval(DD23,D2345) + Dval(DD123,D2345)
	dabbr28 = Dval(DD22,D2345) + Dval(DD122,D2345) + 
     &   Dval(DD222,D2345) + Dval(DD223,D2345)
	dabbr25 = Dval(DD13,D2345) + Dval(DD113,D2345) + 
     &   Dval(DD123,D2345) + Dval(DD133,D2345)
	dabbr22 = Dval(DD12,D2345) + Dval(DD112,D2345) + 
     &   Dval(DD122,D2345) + Dval(DD123,D2345)
	dabbr21 = Dval(DD11,D2345) + Dval(DD111,D2345) + 
     &   Dval(DD112,D2345) + Dval(DD113,D2345)
	dabbr3 = dabbr5 + Dval(DD001,D2345)
	dabbr8 = Dval(DD3,D2345) + Dval(DD13,D2345) + 
     &   Dval(DD23,D2345) + Dval(DD33,D2345)
	dabbr6 = Dval(DD2,D2345) + Dval(DD12,D2345) + 
     &   Dval(DD22,D2345) + Dval(DD23,D2345)
	dabbr4 = Dval(DD1,D2345) + Dval(DD11,D2345) + 
     &   Dval(DD12,D2345) + Dval(DD13,D2345)
	dabbr1 = Dval(DD0,D2345) + Dval(DD1,D2345) + 
     &   Dval(DD2,D2345) + Dval(DD3,D2345)
	dabbr61 = help1 + 2*Dval(DD0023,D2345) + 
     &   Dval(DD00223,D2345) + 2*Dval(DD00233,D2345)
	dabbr66 = help1 + 2*Dval(DD0013,D2345) + 
     &   Dval(DD00113,D2345) + 2*Dval(DD00133,D2345)
	dabbr57 = help2 + 2*Dval(DD0022,D2345) + 
     &   Dval(DD00222,D2345) + 2*Dval(DD00223,D2345)
	dabbr47 = help2 + 2*Dval(DD0012,D2345) + 
     &   Dval(DD00112,D2345) + 2*Dval(DD00123,D2345)
	dabbr51 = Dval(DD003,D2345) + 2*Dval(DD0013,D2345) + 
     &   2*Dval(DD0023,D2345) + Dval(DD00113,D2345) + 
     &   2*Dval(DD00123,D2345) + Dval(DD00223,D2345)
	dabbr64 = Dval(DD002,D2345) + 2*Dval(DD0012,D2345) + 
     &   2*Dval(DD0022,D2345) + Dval(DD00112,D2345) + 
     &   2*Dval(DD00122,D2345) + Dval(DD00222,D2345)
	dabbr44 = help3 + 2*Dval(DD0012,D2345) + 
     &   Dval(DD00122,D2345) + 2*Dval(DD00123,D2345)
	dabbr56 = help3 + 2*Dval(DD0011,D2345) + 
     &   Dval(DD00111,D2345) + 2*Dval(DD00113,D2345)
	dabbr59 = Dval(DD001,D2345) + 2*Dval(DD0011,D2345) + 
     &   2*Dval(DD0012,D2345) + Dval(DD00111,D2345) + 
     &   2*Dval(DD00112,D2345) + Dval(DD00122,D2345)
	dabbr14 = help4 + 2*Dval(DD002,D2345) + 
     &   Dval(DD0022,D2345) + 2*Dval(DD0023,D2345)
	dabbr17 = help4 + 2*Dval(DD001,D2345) + 
     &   Dval(DD0011,D2345) + 2*Dval(DD0013,D2345)
	dabbr19 = Dval(DD00,D2345) + 2*Dval(DD001,D2345) + 
     &   2*Dval(DD002,D2345) + Dval(DD0011,D2345) + 
     &   2*Dval(DD0012,D2345) + Dval(DD0022,D2345)
	dabbr67 = Dval(DD33,D2345) + 2*Dval(DD133,D2345) + 
     &   2*Dval(DD233,D2345) + 2*Dval(DD333,D2345) + 
     &   Dval(DD1133,D2345) + 2*Dval(DD1233,D2345) + 
     &   2*Dval(DD1333,D2345) + Dval(DD2233,D2345) + 
     &   2*Dval(DD2333,D2345) + Dval(DD3333,D2345)
	dabbr63 = Dval(DD23,D2345) + 2*Dval(DD123,D2345) + 
     &   2*Dval(DD223,D2345) + 2*Dval(DD233,D2345) + 
     &   Dval(DD1123,D2345) + 2*Dval(DD1223,D2345) + 
     &   2*Dval(DD1233,D2345) + Dval(DD2223,D2345) + 
     &   2*Dval(DD2233,D2345) + Dval(DD2333,D2345)
	dabbr62 = Dval(DD22,D2345) + 2*Dval(DD122,D2345) + 
     &   2*Dval(DD222,D2345) + 2*Dval(DD223,D2345) + 
     &   Dval(DD1122,D2345) + 2*Dval(DD1222,D2345) + 
     &   2*Dval(DD1223,D2345) + Dval(DD2222,D2345) + 
     &   2*Dval(DD2223,D2345) + Dval(DD2233,D2345)
	dabbr58 = Dval(DD13,D2345) + 2*Dval(DD113,D2345) + 
     &   2*Dval(DD123,D2345) + 2*Dval(DD133,D2345) + 
     &   Dval(DD1113,D2345) + 2*Dval(DD1123,D2345) + 
     &   2*Dval(DD1133,D2345) + Dval(DD1223,D2345) + 
     &   2*Dval(DD1233,D2345) + Dval(DD1333,D2345)
	dabbr54 = Dval(DD12,D2345) + 2*Dval(DD112,D2345) + 
     &   2*Dval(DD122,D2345) + 2*Dval(DD123,D2345) + 
     &   Dval(DD1112,D2345) + 2*Dval(DD1122,D2345) + 
     &   2*Dval(DD1123,D2345) + Dval(DD1222,D2345) + 
     &   2*Dval(DD1223,D2345) + Dval(DD1233,D2345)
	dabbr53 = Dval(DD11,D2345) + 2*Dval(DD111,D2345) + 
     &   2*Dval(DD112,D2345) + 2*Dval(DD113,D2345) + 
     &   Dval(DD1111,D2345) + 2*Dval(DD1112,D2345) + 
     &   2*Dval(DD1113,D2345) + Dval(DD1122,D2345) + 
     &   2*Dval(DD1123,D2345) + Dval(DD1133,D2345)
	dabbr11 = dabbr14 + 2*Dval(DD001,D2345) + 
     &   Dval(DD0011,D2345) + 2*Dval(DD0012,D2345) + 
     &   2*Dval(DD0013,D2345)
	dabbr18 = Dval(DD3,D2345) + 2*Dval(DD13,D2345) + 
     &   2*Dval(DD23,D2345) + 2*Dval(DD33,D2345) + 
     &   Dval(DD113,D2345) + 2*Dval(DD123,D2345) + 
     &   2*Dval(DD133,D2345) + Dval(DD223,D2345) + 
     &   2*Dval(DD233,D2345) + Dval(DD333,D2345)
	dabbr15 = Dval(DD2,D2345) + 2*Dval(DD12,D2345) + 
     &   2*Dval(DD22,D2345) + 2*Dval(DD23,D2345) + 
     &   Dval(DD112,D2345) + 2*Dval(DD122,D2345) + 
     &   2*Dval(DD123,D2345) + Dval(DD222,D2345) + 
     &   2*Dval(DD223,D2345) + Dval(DD233,D2345)
	dabbr12 = Dval(DD1,D2345) + 2*Dval(DD11,D2345) + 
     &   2*Dval(DD12,D2345) + 2*Dval(DD13,D2345) + 
     &   Dval(DD111,D2345) + 2*Dval(DD112,D2345) + 
     &   2*Dval(DD113,D2345) + Dval(DD122,D2345) + 
     &   2*Dval(DD123,D2345) + Dval(DD133,D2345)
	dabbr2 = Dval(DD0,D2345) + 2*Dval(DD1,D2345) + 
     &   2*Dval(DD2,D2345) + 2*Dval(DD3,D2345) + 
     &   Dval(DD11,D2345) + 2*Dval(DD12,D2345) + 
     &   2*Dval(DD13,D2345) + Dval(DD22,D2345) + 
     &   2*Dval(DD23,D2345) + Dval(DD33,D2345)
	dabbr40 = Dval(DD0,D2345) + 4*Dval(DD1,D2345) + 
     &   4*Dval(DD2,D2345) + 4*Dval(DD3,D2345) + 
     &   6*Dval(DD11,D2345) + 12*Dval(DD12,D2345) + 
     &   12*Dval(DD13,D2345) + 6*Dval(DD22,D2345) + 
     &   12*Dval(DD23,D2345) + 6*Dval(DD33,D2345) + 
     &   4*Dval(DD111,D2345) + 12*Dval(DD112,D2345) + 
     &   12*Dval(DD113,D2345) + 12*Dval(DD122,D2345) + 
     &   24*Dval(DD123,D2345) + 12*Dval(DD133,D2345) + 
     &   4*Dval(DD222,D2345) + 12*Dval(DD223,D2345) + 
     &   12*Dval(DD233,D2345) + 4*Dval(DD333,D2345) + 
     &   Dval(DD1111,D2345) + 4*Dval(DD1112,D2345) + 
     &   4*Dval(DD1113,D2345) + 6*Dval(DD1122,D2345) + 
     &   12*Dval(DD1123,D2345) + 6*Dval(DD1133,D2345) + 
     &   4*Dval(DD1222,D2345) + 12*Dval(DD1223,D2345) + 
     &   12*Dval(DD1233,D2345) + 4*Dval(DD1333,D2345) + 
     &   Dval(DD2222,D2345) + 4*Dval(DD2223,D2345) + 
     &   6*Dval(DD2233,D2345) + 4*Dval(DD2333,D2345) + 
     &   Dval(DD3333,D2345)

	call XE0para(E(EE0), para)

	E(EE1) = dabbr1*eta(1) - eta(3)*Dval(DD1,D1245) - 
     &   eta(4)*Dval(DD1,D1235) - eta(5)*Dval(DD1,D1234) + 
     &   2*(zeta(1,1)*Dval(DD00,D1345) + 
     &      zeta(1,2)*Dval(DD00,D1245) + 
     &      zeta(1,3)*Dval(DD00,D1235) + 
     &      zeta(1,4)*Dval(DD00,D1234) - 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,3) + zeta(1,4))*
     &       Dval(DD00,D2345))
	E(EE2) = -(eta(2)*Dval(DD1,D1345)) - 
     &   eta(4)*Dval(DD2,D1235) - eta(5)*Dval(DD2,D1234) - 
     &   eta(1)*Dval(DD1,D2345) + 
     &   2*zeta(1,2)*Dval(DD00,D1345) + 
     &   2*zeta(2,2)*Dval(DD00,D1245) + 
     &   2*zeta(2,3)*Dval(DD00,D1235) + 
     &   2*zeta(2,4)*Dval(DD00,D1234) - 
     &   2*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD00,D2345)
	E(EE3) = -(eta(2)*Dval(DD2,D1345)) - 
     &   eta(3)*Dval(DD2,D1245) - eta(5)*Dval(DD3,D1234) - 
     &   eta(1)*Dval(DD2,D2345) + 
     &   2*zeta(1,3)*Dval(DD00,D1345) + 
     &   2*zeta(2,3)*Dval(DD00,D1245) + 
     &   2*zeta(3,3)*Dval(DD00,D1235) + 
     &   2*zeta(3,4)*Dval(DD00,D1234) - 
     &   2*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD00,D2345)
	E(EE4) = -(eta(2)*Dval(DD3,D1345)) - 
     &   eta(3)*Dval(DD3,D1245) - eta(4)*Dval(DD3,D1235) - 
     &   eta(1)*Dval(DD3,D2345) + 
     &   2*zeta(1,4)*Dval(DD00,D1345) + 
     &   2*zeta(2,4)*Dval(DD00,D1245) + 
     &   2*zeta(3,4)*Dval(DD00,D1235) + 
     &   2*zeta(4,4)*Dval(DD00,D1234) - 
     &   2*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD00,D2345)

	E(EE00) = -(eta(2)*Dval(DD00,D1345)) - 
     &   eta(3)*Dval(DD00,D1245) - eta(4)*Dval(DD00,D1235) - 
     &   eta(5)*Dval(DD00,D1234) - eta(1)*Dval(DD00,D2345)
	E(EE11) = -(dabbr2*eta(1)) - eta(3)*Dval(DD11,D1245) - 
     &   eta(4)*Dval(DD11,D1235) - eta(5)*Dval(DD11,D1234) + 
     &   4*(dabbr3*(zeta(1,1) + zeta(1,2) + zeta(1,3) + 
     &         zeta(1,4)) + zeta(1,2)*Dval(DD001,D1245) + 
     &      zeta(1,3)*Dval(DD001,D1235) + 
     &      zeta(1,4)*Dval(DD001,D1234))
	E(EE12) = dabbr4*eta(1) - eta(4)*Dval(DD12,D1235) - 
     &   eta(5)*Dval(DD12,D1234) + 
     &   2*(dabbr5*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + zeta(1,1)*Dval(DD001,D1345) + 
     &      zeta(2,2)*Dval(DD001,D1245) + 
     &      zeta(2,3)*Dval(DD001,D1235) + 
     &      zeta(2,4)*Dval(DD001,D1234) + 
     &      zeta(1,3)*Dval(DD002,D1235) + 
     &      zeta(1,4)*Dval(DD002,D1234) - 
     &      (zeta(1,1) + zeta(1,3) + zeta(1,4) - zeta(2,2) - 
     &         zeta(2,3) - zeta(2,4))*Dval(DD001,D2345))
	E(EE13) = dabbr6*eta(1) - eta(3)*Dval(DD12,D1245) - 
     &   eta(5)*Dval(DD13,D1234) + 
     &   2*(dabbr7*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + zeta(2,3)*Dval(DD001,D1245) + 
     &      zeta(3,3)*Dval(DD001,D1235) + 
     &      zeta(3,4)*Dval(DD001,D1234) + 
     &      zeta(1,1)*Dval(DD002,D1345) + 
     &      zeta(1,2)*Dval(DD002,D1245) + 
     &      zeta(1,4)*Dval(DD003,D1234) - 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,4) - zeta(2,3) - 
     &         zeta(3,3) - zeta(3,4))*Dval(DD002,D2345))
	E(EE14) = dabbr8*eta(1) - eta(3)*Dval(DD13,D1245) - 
     &   eta(4)*Dval(DD13,D1235) + 
     &   2*(dabbr9*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + zeta(2,4)*Dval(DD001,D1245) + 
     &      zeta(3,4)*Dval(DD001,D1235) + 
     &      zeta(4,4)*Dval(DD001,D1234) + 
     &      zeta(1,1)*Dval(DD003,D1345) + 
     &      zeta(1,2)*Dval(DD003,D1245) + 
     &      zeta(1,3)*Dval(DD003,D1235) - 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,3) - zeta(2,4) - 
     &         zeta(3,4) - zeta(4,4))*Dval(DD003,D2345))
	E(EE22) = -(eta(2)*Dval(DD11,D1345)) - 
     &   eta(4)*Dval(DD22,D1235) - eta(5)*Dval(DD22,D1234) - 
     &   eta(1)*Dval(DD11,D2345) + 
     &   4*zeta(1,2)*Dval(DD001,D1345) + 
     &   4*zeta(2,3)*Dval(DD002,D1235) + 
     &   4*zeta(2,4)*Dval(DD002,D1234) - 
     &   4*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD001,D2345)
	E(EE23) = -(eta(2)*Dval(DD12,D1345)) - 
     &   eta(5)*Dval(DD23,D1234) - eta(1)*Dval(DD12,D2345) + 
     &   2*(zeta(1,3)*Dval(DD001,D1345) + 
     &      zeta(1,2)*Dval(DD002,D1345) + 
     &      zeta(2,2)*Dval(DD002,D1245) + 
     &      zeta(3,3)*Dval(DD002,D1235) + 
     &      zeta(3,4)*Dval(DD002,D1234) + 
     &      zeta(2,4)*Dval(DD003,D1234) - 
     &      (zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD001,D2345) - 
     &      (zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       Dval(DD002,D2345))
	E(EE24) = -(eta(2)*Dval(DD13,D1345)) - 
     &   eta(4)*Dval(DD23,D1235) - eta(1)*Dval(DD13,D2345) + 
     &   2*(zeta(1,4)*Dval(DD001,D1345) + 
     &      zeta(3,4)*Dval(DD002,D1235) + 
     &      zeta(4,4)*Dval(DD002,D1234) + 
     &      zeta(1,2)*Dval(DD003,D1345) + 
     &      zeta(2,2)*Dval(DD003,D1245) + 
     &      zeta(2,3)*Dval(DD003,D1235) - 
     &      (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD001,D2345) - 
     &      (zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       Dval(DD003,D2345))
	E(EE33) = -(eta(2)*Dval(DD22,D1345)) - 
     &   eta(3)*Dval(DD22,D1245) - eta(5)*Dval(DD33,D1234) - 
     &   eta(1)*Dval(DD22,D2345) + 
     &   4*zeta(1,3)*Dval(DD002,D1345) + 
     &   4*zeta(2,3)*Dval(DD002,D1245) + 
     &   4*zeta(3,4)*Dval(DD003,D1234) - 
     &   4*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD002,D2345)
	E(EE34) = -(eta(2)*Dval(DD23,D1345)) - 
     &   eta(3)*Dval(DD23,D1245) - eta(1)*Dval(DD23,D2345) + 
     &   2*(zeta(1,4)*Dval(DD002,D1345) + 
     &      zeta(2,4)*Dval(DD002,D1245) + 
     &      zeta(1,3)*Dval(DD003,D1345) + 
     &      zeta(2,3)*Dval(DD003,D1245) + 
     &      zeta(3,3)*Dval(DD003,D1235) + 
     &      zeta(4,4)*Dval(DD003,D1234) - 
     &      (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD002,D2345) - 
     &      (zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD003,D2345))
	E(EE44) = -(eta(2)*Dval(DD33,D1345)) - 
     &   eta(3)*Dval(DD33,D1245) - eta(4)*Dval(DD33,D1235) - 
     &   eta(1)*Dval(DD33,D2345) + 
     &   4*zeta(1,4)*Dval(DD003,D1345) + 
     &   4*zeta(2,4)*Dval(DD003,D1245) + 
     &   4*zeta(3,4)*Dval(DD003,D1235) - 
     &   4*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD003,D2345)

	E(EE001) = dabbr3*eta(1) - eta(3)*Dval(DD001,D1245) - 
     &   eta(4)*Dval(DD001,D1235) - eta(5)*Dval(DD001,D1234) + 
     &   2*(zeta(1,1)*(Dval(DD0000,D1345) - del) + 
     &      zeta(1,2)*(Dval(DD0000,D1245) - del) + 
     &      zeta(1,3)*(Dval(DD0000,D1235) - del) + 
     &      zeta(1,4)*(Dval(DD0000,D1234) - del) - 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,3) + zeta(1,4))*
     &       (Dval(DD0000,D2345) - del))
	E(EE002) = -(eta(2)*Dval(DD001,D1345)) - 
     &   eta(4)*Dval(DD002,D1235) - eta(5)*Dval(DD002,D1234) - 
     &   eta(1)*Dval(DD001,D2345) + 
     &   2*zeta(1,2)*(Dval(DD0000,D1345) - del) + 
     &   2*zeta(2,2)*(Dval(DD0000,D1245) - del) + 
     &   2*zeta(2,3)*(Dval(DD0000,D1235) - del) + 
     &   2*zeta(2,4)*(Dval(DD0000,D1234) - del) - 
     &   2*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    (Dval(DD0000,D2345) - del)
	E(EE003) = -(eta(2)*Dval(DD002,D1345)) - 
     &   eta(3)*Dval(DD002,D1245) - eta(5)*Dval(DD003,D1234) - 
     &   eta(1)*Dval(DD002,D2345) + 
     &   2*zeta(1,3)*(Dval(DD0000,D1345) - del) + 
     &   2*zeta(2,3)*(Dval(DD0000,D1245) - del) + 
     &   2*zeta(3,3)*(Dval(DD0000,D1235) - del) + 
     &   2*zeta(3,4)*(Dval(DD0000,D1234) - del) - 
     &   2*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    (Dval(DD0000,D2345) - del)
	E(EE004) = -(eta(2)*Dval(DD003,D1345)) - 
     &   eta(3)*Dval(DD003,D1245) - eta(4)*Dval(DD003,D1235) - 
     &   eta(1)*Dval(DD003,D2345) + 
     &   2*zeta(1,4)*(Dval(DD0000,D1345) - del) + 
     &   2*zeta(2,4)*(Dval(DD0000,D1245) - del) + 
     &   2*zeta(3,4)*(Dval(DD0000,D1235) - del) + 
     &   2*zeta(4,4)*(Dval(DD0000,D1234) - del) - 
     &   2*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    (Dval(DD0000,D2345) - del)
	E(EE111) = -6*dabbr11*
     &    (zeta(1,1) + zeta(1,2) + zeta(1,3) + zeta(1,4)) - 
     &   eta(3)*Dval(DD111,D1245) - eta(4)*Dval(DD111,D1235) - 
     &   eta(5)*Dval(DD111,D1234) + 
     &   eta(1)*(3*dabbr10 + Dval(DD0,D2345) + 
     &      3*Dval(DD1,D2345) + 3*Dval(DD2,D2345) + 
     &      3*Dval(DD3,D2345) + 3*Dval(DD11,D2345) + 
     &      6*Dval(DD12,D2345) + 6*Dval(DD13,D2345) + 
     &      3*Dval(DD22,D2345) + 6*Dval(DD23,D2345) + 
     &      3*Dval(DD33,D2345) + Dval(DD111,D2345) + 
     &      3*Dval(DD112,D2345) + 3*Dval(DD113,D2345) + 
     &      3*Dval(DD122,D2345) + 6*Dval(DD123,D2345) + 
     &      3*Dval(DD133,D2345) + Dval(DD222,D2345) + 
     &      Dval(DD333,D2345)) + 
     &   6*zeta(1,2)*Dval(DD0011,D1245) + 
     &   6*zeta(1,3)*Dval(DD0011,D1235) + 
     &   6*zeta(1,4)*Dval(DD0011,D1234)
	E(EE112) = -(dabbr12*eta(1)) - 
     &   2*dabbr14*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &      zeta(2,4)) - eta(4)*Dval(DD112,D1235) - 
     &   eta(5)*Dval(DD112,D1234) + 
     &   2*(2*dabbr13*(zeta(1,1) + zeta(1,3) + zeta(1,4) - 
     &         zeta(2,2) - zeta(2,3) - zeta(2,4)) + 
     &      zeta(2,2)*Dval(DD0011,D1245) + 
     &      zeta(2,3)*Dval(DD0011,D1235) + 
     &      zeta(2,4)*Dval(DD0011,D1234) + 
     &      2*zeta(1,3)*Dval(DD0012,D1235) + 
     &      2*zeta(1,4)*Dval(DD0012,D1234) + 
     &      (2*zeta(1,1) + zeta(1,2) + 2*zeta(1,3) + 2*zeta(1,4) - 
     &         zeta(2,2) - zeta(2,3) - zeta(2,4))*
     &       Dval(DD0011,D2345))
	E(EE113) = -(dabbr15*eta(1)) - 
     &   2*dabbr17*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &      zeta(3,4)) - eta(3)*Dval(DD112,D1245) - 
     &   eta(5)*Dval(DD113,D1234) + 
     &   2*(2*dabbr16*(zeta(1,1) + zeta(1,2) + zeta(1,4) - 
     &         zeta(2,3) - zeta(3,3) - zeta(3,4)) + 
     &      zeta(2,3)*Dval(DD0011,D1245) + 
     &      zeta(3,3)*Dval(DD0011,D1235) + 
     &      zeta(3,4)*Dval(DD0011,D1234) + 
     &      2*zeta(1,2)*Dval(DD0012,D1245) + 
     &      2*zeta(1,4)*Dval(DD0013,D1234) + 
     &      (2*zeta(1,1) + 2*zeta(1,2) + zeta(1,3) + 2*zeta(1,4) - 
     &         zeta(2,3) - zeta(3,3) - zeta(3,4))*
     &       Dval(DD0022,D2345))
	E(EE114) = -(dabbr18*eta(1)) - 
     &   2*dabbr19*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &      zeta(4,4)) - eta(3)*Dval(DD113,D1245) - 
     &   eta(4)*Dval(DD113,D1235) + 
     &   2*(2*dabbr20*(zeta(1,1) + zeta(1,2) + zeta(1,3) - 
     &         zeta(2,4) - zeta(3,4) - zeta(4,4)) + 
     &      zeta(2,4)*Dval(DD0011,D1245) + 
     &      zeta(3,4)*Dval(DD0011,D1235) + 
     &      zeta(4,4)*Dval(DD0011,D1234) + 
     &      2*zeta(1,2)*Dval(DD0013,D1245) + 
     &      2*zeta(1,3)*Dval(DD0013,D1235) + 
     &      (2*zeta(1,1) + 2*zeta(1,2) + 2*zeta(1,3) + zeta(1,4) - 
     &         zeta(2,4) - zeta(3,4) - zeta(4,4))*
     &       Dval(DD0033,D2345))
	E(EE122) = dabbr21*eta(1) - eta(4)*Dval(DD122,D1235) - 
     &   eta(5)*Dval(DD122,D1234) + 
     &   2*(2*dabbr13*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + zeta(1,1)*Dval(DD0011,D1345) + 
     &      2*zeta(2,3)*Dval(DD0012,D1235) + 
     &      2*zeta(2,4)*Dval(DD0012,D1234) + 
     &      zeta(1,3)*Dval(DD0022,D1235) + 
     &      zeta(1,4)*Dval(DD0022,D1234) - 
     &      (zeta(1,1) - zeta(1,2) + zeta(1,3) + zeta(1,4) - 
     &         2*(zeta(2,2) + zeta(2,3) + zeta(2,4)))*
     &       Dval(DD0011,D2345))
	E(EE123) = dabbr22*eta(1) - eta(5)*Dval(DD123,D1234) + 
     &   2*(dabbr24*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + 
     &      dabbr23*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + zeta(1,1)*Dval(DD0012,D1345) + 
     &      zeta(2,2)*Dval(DD0012,D1245) + 
     &      zeta(3,3)*Dval(DD0012,D1235) + 
     &      zeta(3,4)*Dval(DD0012,D1234) + 
     &      zeta(2,4)*Dval(DD0013,D1234) + 
     &      zeta(1,4)*Dval(DD0023,D1234) + 
     &      (-zeta(1,1) - zeta(1,4) + zeta(2,2) + 2*zeta(2,3) + 
     &         zeta(2,4) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD0012,D2345))
	E(EE124) = dabbr25*eta(1) - eta(4)*Dval(DD123,D1235) + 
     &   2*(dabbr27*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + 
     &      dabbr26*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + zeta(3,4)*Dval(DD0012,D1235) + 
     &      zeta(4,4)*Dval(DD0012,D1234) + 
     &      zeta(1,1)*Dval(DD0013,D1345) + 
     &      zeta(2,2)*Dval(DD0013,D1245) + 
     &      zeta(2,3)*Dval(DD0013,D1235) + 
     &      zeta(1,3)*Dval(DD0023,D1235) + 
     &      (-zeta(1,1) - zeta(1,3) + zeta(2,2) + zeta(2,3) + 
     &         2*zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD0013,D2345))
	E(EE133) = dabbr28*eta(1) - eta(3)*Dval(DD122,D1245) - 
     &   eta(5)*Dval(DD133,D1234) + 
     &   2*(2*dabbr16*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + 2*zeta(2,3)*Dval(DD0012,D1245) + 
     &      2*zeta(3,4)*Dval(DD0013,D1234) + 
     &      zeta(1,1)*Dval(DD0022,D1345) + 
     &      zeta(1,2)*Dval(DD0022,D1245) + 
     &      zeta(1,4)*Dval(DD0033,D1234) - 
     &      (zeta(1,1) + zeta(1,2) - zeta(1,3) + zeta(1,4) - 
     &         2*(zeta(2,3) + zeta(3,3) + zeta(3,4)))*
     &       Dval(DD0022,D2345))
	E(EE134) = dabbr29*eta(1) - eta(3)*Dval(DD123,D1245) + 
     &   2*(dabbr31*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + 
     &      dabbr30*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + zeta(2,4)*Dval(DD0012,D1245) + 
     &      zeta(2,3)*Dval(DD0013,D1245) + 
     &      zeta(3,3)*Dval(DD0013,D1235) + 
     &      zeta(4,4)*Dval(DD0013,D1234) + 
     &      zeta(1,1)*Dval(DD0023,D1345) + 
     &      zeta(1,2)*Dval(DD0023,D1245) + 
     &      (-zeta(1,1) - zeta(1,2) + zeta(2,3) + zeta(2,4) + 
     &         zeta(3,3) + 2*zeta(3,4) + zeta(4,4))*
     &       Dval(DD0023,D2345))
	E(EE144) = dabbr32*eta(1) - eta(3)*Dval(DD133,D1245) - 
     &   eta(4)*Dval(DD133,D1235) + 
     &   2*(2*dabbr20*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + 2*zeta(2,4)*Dval(DD0013,D1245) + 
     &      2*zeta(3,4)*Dval(DD0013,D1235) + 
     &      zeta(1,1)*Dval(DD0033,D1345) + 
     &      zeta(1,2)*Dval(DD0033,D1245) + 
     &      zeta(1,3)*Dval(DD0033,D1235) - 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,3) - zeta(1,4) - 
     &         2*(zeta(2,4) + zeta(3,4) + zeta(4,4)))*
     &       Dval(DD0033,D2345))
	E(EE222) = -(eta(2)*Dval(DD111,D1345)) - 
     &   eta(4)*Dval(DD222,D1235) - eta(5)*Dval(DD222,D1234) - 
     &   eta(1)*Dval(DD111,D2345) + 
     &   6*zeta(1,2)*Dval(DD0011,D1345) + 
     &   6*zeta(2,3)*Dval(DD0022,D1235) + 
     &   6*zeta(2,4)*Dval(DD0022,D1234) - 
     &   6*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD0011,D2345)
	E(EE223) = -(eta(2)*Dval(DD112,D1345)) - 
     &   eta(5)*Dval(DD223,D1234) - eta(1)*Dval(DD112,D2345) + 
     &   2*zeta(1,3)*Dval(DD0011,D1345) + 
     &   4*zeta(1,2)*Dval(DD0012,D1345) + 
     &   2*zeta(3,3)*Dval(DD0022,D1235) + 
     &   2*zeta(3,4)*Dval(DD0022,D1234) + 
     &   4*zeta(2,4)*Dval(DD0023,D1234) - 
     &   2*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD0011,D2345) - 
     &   4*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD0012,D2345)
	E(EE224) = -(eta(2)*Dval(DD113,D1345)) - 
     &   eta(4)*Dval(DD223,D1235) - eta(1)*Dval(DD113,D2345) + 
     &   2*zeta(1,4)*Dval(DD0011,D1345) + 
     &   4*zeta(1,2)*Dval(DD0013,D1345) + 
     &   2*zeta(3,4)*Dval(DD0022,D1235) + 
     &   2*zeta(4,4)*Dval(DD0022,D1234) + 
     &   4*zeta(2,3)*Dval(DD0023,D1235) - 
     &   2*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD0011,D2345) - 
     &   4*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD0013,D2345)
	E(EE233) = -(eta(2)*Dval(DD122,D1345)) - 
     &   eta(5)*Dval(DD233,D1234) - eta(1)*Dval(DD122,D2345) + 
     &   4*zeta(1,3)*Dval(DD0012,D1345) + 
     &   2*zeta(1,2)*Dval(DD0022,D1345) + 
     &   2*zeta(2,2)*Dval(DD0022,D1245) + 
     &   4*zeta(3,4)*Dval(DD0023,D1234) + 
     &   2*zeta(2,4)*Dval(DD0033,D1234) - 
     &   4*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD0012,D2345) - 
     &   2*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD0022,D2345)
	E(EE234) = -(eta(2)*Dval(DD123,D1345)) - 
     &   eta(1)*Dval(DD123,D2345) + 
     &   2*(zeta(1,4)*Dval(DD0012,D1345) + 
     &      zeta(1,3)*Dval(DD0013,D1345) + 
     &      zeta(1,2)*Dval(DD0023,D1345) + 
     &      zeta(2,2)*Dval(DD0023,D1245) + 
     &      zeta(3,3)*Dval(DD0023,D1235) + 
     &      zeta(4,4)*Dval(DD0023,D1234) - 
     &      (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD0012,D2345) - 
     &      (zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD0013,D2345) - 
     &      (zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       Dval(DD0023,D2345))
	E(EE244) = -(eta(2)*Dval(DD133,D1345)) - 
     &   eta(4)*Dval(DD233,D1235) - eta(1)*Dval(DD133,D2345) + 
     &   4*zeta(1,4)*Dval(DD0013,D1345) + 
     &   4*zeta(3,4)*Dval(DD0023,D1235) + 
     &   2*zeta(1,2)*Dval(DD0033,D1345) + 
     &   2*zeta(2,2)*Dval(DD0033,D1245) + 
     &   2*zeta(2,3)*Dval(DD0033,D1235) - 
     &   4*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD0013,D2345) - 
     &   2*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD0033,D2345)
	E(EE333) = -(eta(2)*Dval(DD222,D1345)) - 
     &   eta(3)*Dval(DD222,D1245) - eta(5)*Dval(DD333,D1234) - 
     &   eta(1)*Dval(DD222,D2345) + 
     &   6*zeta(1,3)*Dval(DD0022,D1345) + 
     &   6*zeta(2,3)*Dval(DD0022,D1245) + 
     &   6*zeta(3,4)*Dval(DD0033,D1234) - 
     &   6*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD0022,D2345)
	E(EE334) = -(eta(2)*Dval(DD223,D1345)) - 
     &   eta(3)*Dval(DD223,D1245) - eta(1)*Dval(DD223,D2345) + 
     &   2*zeta(1,4)*Dval(DD0022,D1345) + 
     &   2*zeta(2,4)*Dval(DD0022,D1245) + 
     &   4*zeta(1,3)*Dval(DD0023,D1345) + 
     &   4*zeta(2,3)*Dval(DD0023,D1245) + 
     &   2*zeta(4,4)*Dval(DD0033,D1234) - 
     &   2*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD0022,D2345) - 
     &   4*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD0023,D2345)
	E(EE344) = -(eta(2)*Dval(DD233,D1345)) - 
     &   eta(3)*Dval(DD233,D1245) - eta(1)*Dval(DD233,D2345) + 
     &   4*zeta(1,4)*Dval(DD0023,D1345) + 
     &   4*zeta(2,4)*Dval(DD0023,D1245) + 
     &   2*zeta(1,3)*Dval(DD0033,D1345) + 
     &   2*zeta(2,3)*Dval(DD0033,D1245) + 
     &   2*zeta(3,3)*Dval(DD0033,D1235) - 
     &   4*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD0023,D2345) - 
     &   2*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD0033,D2345)
	E(EE444) = -(eta(2)*Dval(DD333,D1345)) - 
     &   eta(3)*Dval(DD333,D1245) - eta(4)*Dval(DD333,D1235) - 
     &   eta(1)*Dval(DD333,D2345) + 
     &   6*zeta(1,4)*Dval(DD0033,D1345) + 
     &   6*zeta(2,4)*Dval(DD0033,D1245) + 
     &   6*zeta(3,4)*Dval(DD0033,D1235) - 
     &   6*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD0033,D2345)

	E(EE0000) = (-finite*(detZ/detY) - 
     &     48*(eta(2)*(Dval(DD0000,D1345) - del) + 
     &        eta(3)*(Dval(DD0000,D1245) - del) + 
     &        eta(4)*(Dval(DD0000,D1235) - del) + 
     &        eta(5)*(Dval(DD0000,D1234) - del) + 
     &        eta(1)*(Dval(DD0000,D2345) - del)))/48D0
	E(EE0011) = -(dabbr11*eta(1)) - 
     &   eta(3)*Dval(DD0011,D1245) - 
     &   eta(4)*Dval(DD0011,D1235) - 
     &   eta(5)*Dval(DD0011,D1234) + 
     &   4*(dabbr33*(zeta(1,1) + zeta(1,2) + zeta(1,3) + 
     &         zeta(1,4)) + 
     &      zeta(1,2)*(Dval(DD00001,D1245) + del4) + 
     &      zeta(1,3)*(Dval(DD00001,D1235) + del4) + 
     &      zeta(1,4)*(Dval(DD00001,D1234) + del4))
	E(EE0012) = dabbr34*eta(1) - 
     &   eta(4)*Dval(DD0012,D1235) - 
     &   eta(5)*Dval(DD0012,D1234) + 
     &   2*(dabbr35*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + 
     &      zeta(1,1)*(Dval(DD00001,D1345) + del4) + 
     &      zeta(2,2)*(Dval(DD00001,D1245) + del4) + 
     &      zeta(2,3)*(Dval(DD00001,D1235) + del4) + 
     &      zeta(2,4)*(Dval(DD00001,D1234) + del4) + 
     &      zeta(1,3)*(Dval(DD00002,D1235) + del4) + 
     &      zeta(1,4)*(Dval(DD00002,D1234) + del4) - 
     &      (zeta(1,1) + zeta(1,3) + zeta(1,4) - zeta(2,2) - 
     &         zeta(2,3) - zeta(2,4))*
     &       (Dval(DD00001,D2345) + del4))
	E(EE0013) = dabbr36*eta(1) - 
     &   eta(3)*Dval(DD0012,D1245) - 
     &   eta(5)*Dval(DD0013,D1234) + 
     &   2*(dabbr37*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + 
     &      zeta(2,3)*(Dval(DD00001,D1245) + del4) + 
     &      zeta(3,3)*(Dval(DD00001,D1235) + del4) + 
     &      zeta(3,4)*(Dval(DD00001,D1234) + del4) + 
     &      zeta(1,1)*(Dval(DD00002,D1345) + del4) + 
     &      zeta(1,2)*(Dval(DD00002,D1245) + del4) + 
     &      zeta(1,4)*(Dval(DD00003,D1234) + del4) - 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,4) - zeta(2,3) - 
     &         zeta(3,3) - zeta(3,4))*
     &       (Dval(DD00002,D2345) + del4))
	E(EE0014) = dabbr38*eta(1) - 
     &   eta(3)*Dval(DD0013,D1245) - 
     &   eta(4)*Dval(DD0013,D1235) + 
     &   2*(dabbr39*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + 
     &      zeta(2,4)*(Dval(DD00001,D1245) + del4) + 
     &      zeta(3,4)*(Dval(DD00001,D1235) + del4) + 
     &      zeta(4,4)*(Dval(DD00001,D1234) + del4) + 
     &      zeta(1,1)*(Dval(DD00003,D1345) + del4) + 
     &      zeta(1,2)*(Dval(DD00003,D1245) + del4) + 
     &      zeta(1,3)*(Dval(DD00003,D1235) + del4) - 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,3) - zeta(2,4) - 
     &         zeta(3,4) - zeta(4,4))*
     &       (Dval(DD00003,D2345) + del4))
	E(EE0022) = -(eta(2)*Dval(DD0011,D1345)) - 
     &   eta(4)*Dval(DD0022,D1235) - 
     &   eta(5)*Dval(DD0022,D1234) - 
     &   eta(1)*Dval(DD0011,D2345) + 
     &   4*zeta(1,2)*(Dval(DD00001,D1345) + del4) + 
     &   4*zeta(2,3)*(Dval(DD00002,D1235) + del4) + 
     &   4*zeta(2,4)*(Dval(DD00002,D1234) + del4) - 
     &   4*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    (Dval(DD00001,D2345) + del4)
	E(EE0023) = -(eta(2)*Dval(DD0012,D1345)) - 
     &   eta(5)*Dval(DD0023,D1234) - 
     &   eta(1)*Dval(DD0012,D2345) + 
     &   2*(zeta(1,3)*(Dval(DD00001,D1345) + del4) + 
     &      zeta(1,2)*(Dval(DD00002,D1345) + del4) + 
     &      zeta(2,2)*(Dval(DD00002,D1245) + del4) + 
     &      zeta(3,3)*(Dval(DD00002,D1235) + del4) + 
     &      zeta(3,4)*(Dval(DD00002,D1234) + del4) + 
     &      zeta(2,4)*(Dval(DD00003,D1234) + del4) - 
     &      (zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       (Dval(DD00001,D2345) + del4) - 
     &      (zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       (Dval(DD00002,D2345) + del4))
	E(EE0024) = -(eta(2)*Dval(DD0013,D1345)) - 
     &   eta(4)*Dval(DD0023,D1235) - 
     &   eta(1)*Dval(DD0013,D2345) + 
     &   2*(zeta(1,4)*(Dval(DD00001,D1345) + del4) + 
     &      zeta(3,4)*(Dval(DD00002,D1235) + del4) + 
     &      zeta(4,4)*(Dval(DD00002,D1234) + del4) + 
     &      zeta(1,2)*(Dval(DD00003,D1345) + del4) + 
     &      zeta(2,2)*(Dval(DD00003,D1245) + del4) + 
     &      zeta(2,3)*(Dval(DD00003,D1235) + del4) - 
     &      (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       (Dval(DD00001,D2345) + del4) - 
     &      (zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       (Dval(DD00003,D2345) + del4))
	E(EE0033) = -(eta(2)*Dval(DD0022,D1345)) - 
     &   eta(3)*Dval(DD0022,D1245) - 
     &   eta(5)*Dval(DD0033,D1234) - 
     &   eta(1)*Dval(DD0022,D2345) + 
     &   4*zeta(1,3)*(Dval(DD00002,D1345) + del4) + 
     &   4*zeta(2,3)*(Dval(DD00002,D1245) + del4) + 
     &   4*zeta(3,4)*(Dval(DD00003,D1234) + del4) - 
     &   4*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    (Dval(DD00002,D2345) + del4)
	E(EE0034) = -(eta(2)*Dval(DD0023,D1345)) - 
     &   eta(3)*Dval(DD0023,D1245) - 
     &   eta(1)*Dval(DD0023,D2345) + 
     &   2*(zeta(1,4)*(Dval(DD00002,D1345) + del4) + 
     &      zeta(2,4)*(Dval(DD00002,D1245) + del4) + 
     &      zeta(1,3)*(Dval(DD00003,D1345) + del4) + 
     &      zeta(2,3)*(Dval(DD00003,D1245) + del4) + 
     &      zeta(3,3)*(Dval(DD00003,D1235) + del4) + 
     &      zeta(4,4)*(Dval(DD00003,D1234) + del4) - 
     &      (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       (Dval(DD00002,D2345) + del4) - 
     &      (zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       (Dval(DD00003,D2345) + del4))
	E(EE0044) = -(eta(2)*Dval(DD0033,D1345)) - 
     &   eta(3)*Dval(DD0033,D1245) - 
     &   eta(4)*Dval(DD0033,D1235) - 
     &   eta(1)*Dval(DD0033,D2345) + 
     &   4*zeta(1,4)*(Dval(DD00003,D1345) + del4) + 
     &   4*zeta(2,4)*(Dval(DD00003,D1245) + del4) + 
     &   4*zeta(3,4)*(Dval(DD00003,D1235) + del4) - 
     &   4*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    (Dval(DD00003,D2345) + del4)
	E(EE1111) = -(dabbr40*eta(1)) - 
     &   eta(3)*Dval(DD1111,D1245) - 
     &   eta(4)*Dval(DD1111,D1235) - 
     &   eta(5)*Dval(DD1111,D1234) + 
     &   8*zeta(1,2)*Dval(DD00111,D1245) + 
     &   8*zeta(1,3)*Dval(DD00111,D1235) + 
     &   8*zeta(1,4)*Dval(DD00111,D1234) + 
     &   8*(zeta(1,1) + zeta(1,2) + zeta(1,3) + zeta(1,4))*
     &    (3*dabbr41 + Dval(DD00,D2345) + 
     &      3*Dval(DD001,D2345) + 3*Dval(DD002,D2345) + 
     &      3*Dval(DD003,D2345) + 3*Dval(DD0011,D2345) + 
     &      6*Dval(DD0012,D2345) + 6*Dval(DD0013,D2345) + 
     &      3*Dval(DD0022,D2345) + 6*Dval(DD0023,D2345) + 
     &      3*Dval(DD0033,D2345) + Dval(DD00111,D2345) + 
     &      3*Dval(DD00112,D2345) + 3*Dval(DD00113,D2345) + 
     &      3*Dval(DD00122,D2345) + 6*Dval(DD00123,D2345) + 
     &      3*Dval(DD00133,D2345) + Dval(DD00222,D2345) + 
     &      Dval(DD00333,D2345))
	E(EE1112) = -6*dabbr44*
     &    (zeta(1,1) + zeta(1,3) + zeta(1,4) - zeta(2,2) - 
     &      zeta(2,3) - zeta(2,4)) - 
     &   6*dabbr43*(2*zeta(1,1) + zeta(1,2) + 2*zeta(1,3) + 
     &      2*zeta(1,4) - zeta(2,2) - zeta(2,3) - zeta(2,4)) - 
     &   eta(4)*Dval(DD1112,D1235) - 
     &   eta(5)*Dval(DD1112,D1234) + 
     &   eta(1)*(3*dabbr42 + Dval(DD1,D2345) + 
     &      3*Dval(DD11,D2345) + 3*Dval(DD12,D2345) + 
     &      3*Dval(DD13,D2345) + 3*Dval(DD111,D2345) + 
     &      6*Dval(DD112,D2345) + 6*Dval(DD113,D2345) + 
     &      3*Dval(DD122,D2345) + 6*Dval(DD123,D2345) + 
     &      3*Dval(DD133,D2345) + Dval(DD1111,D2345) + 
     &      3*Dval(DD1112,D2345) + 3*Dval(DD1113,D2345) + 
     &      3*Dval(DD1122,D2345) + 6*Dval(DD1123,D2345) + 
     &      3*Dval(DD1133,D2345) + Dval(DD1222,D2345) + 
     &      Dval(DD1333,D2345)) + 
     &   2*(zeta(2,2)*Dval(DD00111,D1245) + 
     &      zeta(2,3)*Dval(DD00111,D1235) + 
     &      zeta(2,4)*Dval(DD00111,D1234) + 
     &      3*zeta(1,3)*Dval(DD00112,D1235) + 
     &      3*zeta(1,4)*Dval(DD00112,D1234) + 
     &      (-3*zeta(1,1) - 2*zeta(1,2) - 3*zeta(1,3) - 
     &         3*zeta(1,4) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       Dval(DD00111,D2345)) + 
     &   2*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    (3*dabbr41 + Dval(DD00,D2345) + 
     &      3*Dval(DD002,D2345) + 3*Dval(DD003,D2345) + 
     &      3*Dval(DD0022,D2345) + 6*Dval(DD0023,D2345) + 
     &      3*Dval(DD0033,D2345) + Dval(DD00222,D2345) + 
     &      Dval(DD00333,D2345))
	E(EE1113) = -6*dabbr47*
     &    (zeta(1,1) + zeta(1,2) + zeta(1,4) - zeta(2,3) - 
     &      zeta(3,3) - zeta(3,4)) - 
     &   6*dabbr46*(2*zeta(1,1) + 2*zeta(1,2) + zeta(1,3) + 
     &      2*zeta(1,4) - zeta(2,3) - zeta(3,3) - zeta(3,4)) - 
     &   eta(3)*Dval(DD1112,D1245) - 
     &   eta(5)*Dval(DD1113,D1234) + 
     &   eta(1)*(3*dabbr45 + Dval(DD2,D2345) + 
     &      3*Dval(DD12,D2345) + 3*Dval(DD22,D2345) + 
     &      3*Dval(DD23,D2345) + 3*Dval(DD112,D2345) + 
     &      6*Dval(DD122,D2345) + 6*Dval(DD123,D2345) + 
     &      3*Dval(DD222,D2345) + 6*Dval(DD223,D2345) + 
     &      3*Dval(DD233,D2345) + Dval(DD1112,D2345) + 
     &      3*Dval(DD1122,D2345) + 3*Dval(DD1123,D2345) + 
     &      3*Dval(DD1222,D2345) + 6*Dval(DD1223,D2345) + 
     &      3*Dval(DD1233,D2345) + Dval(DD2222,D2345) + 
     &      Dval(DD2333,D2345)) + 
     &   2*(zeta(2,3)*Dval(DD00111,D1245) + 
     &      zeta(3,3)*Dval(DD00111,D1235) + 
     &      zeta(3,4)*Dval(DD00111,D1234) + 
     &      3*zeta(1,2)*Dval(DD00112,D1245) + 
     &      3*zeta(1,4)*Dval(DD00113,D1234) + 
     &      (-3*zeta(1,1) - 3*zeta(1,2) - 2*zeta(1,3) - 
     &         3*zeta(1,4) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD00222,D2345)) + 
     &   2*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    (3*dabbr48 + Dval(DD00,D2345) + 
     &      3*Dval(DD001,D2345) + 3*Dval(DD003,D2345) + 
     &      3*Dval(DD0011,D2345) + 6*Dval(DD0013,D2345) + 
     &      3*Dval(DD0033,D2345) + Dval(DD00111,D2345) + 
     &      Dval(DD00333,D2345))
	E(EE1114) = -6*dabbr51*
     &    (zeta(1,1) + zeta(1,2) + zeta(1,3) - zeta(2,4) - 
     &      zeta(3,4) - zeta(4,4)) - 
     &   6*dabbr52*(2*zeta(1,1) + 2*zeta(1,2) + 2*zeta(1,3) + 
     &      zeta(1,4) - zeta(2,4) - zeta(3,4) - zeta(4,4)) - 
     &   eta(3)*Dval(DD1113,D1245) - 
     &   eta(4)*Dval(DD1113,D1235) + 
     &   eta(1)*(3*dabbr49 + Dval(DD3,D2345) + 
     &      3*Dval(DD13,D2345) + 3*Dval(DD23,D2345) + 
     &      3*Dval(DD33,D2345) + 3*Dval(DD113,D2345) + 
     &      6*Dval(DD123,D2345) + 6*Dval(DD133,D2345) + 
     &      3*Dval(DD223,D2345) + 6*Dval(DD233,D2345) + 
     &      3*Dval(DD333,D2345) + Dval(DD1113,D2345) + 
     &      3*Dval(DD1123,D2345) + 3*Dval(DD1133,D2345) + 
     &      3*Dval(DD1223,D2345) + 6*Dval(DD1233,D2345) + 
     &      3*Dval(DD1333,D2345) + Dval(DD2223,D2345) + 
     &      Dval(DD3333,D2345)) + 
     &   2*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    (3*dabbr50 + Dval(DD00,D2345) + 
     &      3*Dval(DD001,D2345) + 3*Dval(DD002,D2345) + 
     &      3*Dval(DD0011,D2345) + 6*Dval(DD0012,D2345) + 
     &      3*Dval(DD0022,D2345) + Dval(DD00111,D2345) + 
     &      Dval(DD00222,D2345)) + 
     &   2*(zeta(2,4)*Dval(DD00111,D1245) + 
     &      zeta(3,4)*Dval(DD00111,D1235) + 
     &      zeta(4,4)*Dval(DD00111,D1234) + 
     &      3*zeta(1,2)*Dval(DD00113,D1245) + 
     &      3*zeta(1,3)*Dval(DD00113,D1235) + 
     &      (-3*zeta(1,1) - 3*zeta(1,2) - 3*zeta(1,3) - 
     &         2*zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD00333,D2345))
	E(EE1122) = -(dabbr53*eta(1)) - 
     &   4*dabbr44*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &      zeta(2,4)) - eta(4)*Dval(DD1122,D1235) - 
     &   eta(5)*Dval(DD1122,D1234) + 
     &   4*(-2*dabbr43*(zeta(2,2) + zeta(2,3) + zeta(2,4)) + 
     &      zeta(1,1)*Dval(DD0011,D2345) - 
     &      zeta(1,2)*Dval(DD0011,D2345) + 
     &      zeta(1,3)*Dval(DD0011,D2345) + 
     &      zeta(1,4)*Dval(DD0011,D2345) + 
     &      zeta(2,3)*Dval(DD00112,D1235) + 
     &      zeta(2,4)*Dval(DD00112,D1234) + 
     &      zeta(1,3)*Dval(DD00122,D1235) + 
     &      zeta(1,4)*Dval(DD00122,D1234) + 
     &      (zeta(1,1) + zeta(1,3) + zeta(1,4) - zeta(2,2) - 
     &         zeta(2,3) - zeta(2,4))*Dval(DD00111,D2345) + 
     &      zeta(1,1)*Dval(DD00112,D2345) - 
     &      zeta(1,2)*Dval(DD00112,D2345) + 
     &      zeta(1,3)*Dval(DD00112,D2345) + 
     &      zeta(1,4)*Dval(DD00112,D2345) + 
     &      (zeta(1,1) - zeta(1,2) + zeta(1,3) + zeta(1,4))*
     &       Dval(DD00113,D2345))
	E(EE1123) = -(dabbr54*eta(1)) - 
     &   2*dabbr57*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &      zeta(2,4)) - 2*dabbr56*
     &    (zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4)) - 
     &   eta(5)*Dval(DD1123,D1234) + 
     &   2*(2*dabbr55*(zeta(1,1) + zeta(1,4) - zeta(2,2) - 
     &         2*zeta(2,3) - zeta(2,4) - zeta(3,3) - zeta(3,4)) + 
     &      zeta(2,2)*Dval(DD00112,D1245) + 
     &      zeta(3,3)*Dval(DD00112,D1235) + 
     &      zeta(3,4)*Dval(DD00112,D1234) + 
     &      zeta(2,4)*Dval(DD00113,D1234) + 
     &      2*zeta(1,4)*Dval(DD00123,D1234) + 
     &      (2*zeta(1,1) + zeta(1,2) + 2*zeta(1,4) - zeta(2,2) - 
     &         3*zeta(2,3) - zeta(2,4) - 2*zeta(3,3) - 2*zeta(3,4))
     &        *Dval(DD00112,D2345) + 
     &      (2*zeta(1,1) + zeta(1,3) + 2*zeta(1,4) - 2*zeta(2,2) - 
     &         3*zeta(2,3) - 2*zeta(2,4) - zeta(3,3) - zeta(3,4))*
     &       Dval(DD00122,D2345))
	E(EE1124) = -(dabbr58*eta(1)) - 
     &   2*dabbr61*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &      zeta(2,4)) - 2*dabbr59*
     &    (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4)) - 
     &   eta(4)*Dval(DD1123,D1235) + 
     &   2*(2*dabbr60*(zeta(1,1) + zeta(1,3) - zeta(2,2) - 
     &         zeta(2,3) - 2*zeta(2,4) - zeta(3,4) - zeta(4,4)) + 
     &      zeta(3,4)*Dval(DD00112,D1235) + 
     &      zeta(4,4)*Dval(DD00112,D1234) + 
     &      zeta(2,2)*Dval(DD00113,D1245) + 
     &      zeta(2,3)*Dval(DD00113,D1235) + 
     &      2*zeta(1,3)*Dval(DD00123,D1235) + 
     &      (2*zeta(1,1) + zeta(1,2) + 2*zeta(1,3) - zeta(2,2) - 
     &         zeta(2,3) - 3*zeta(2,4) - 2*zeta(3,4) - 2*zeta(4,4))
     &        *Dval(DD00113,D2345) + 
     &      (2*zeta(1,1) + 2*zeta(1,3) + zeta(1,4) - 2*zeta(2,2) - 
     &         2*zeta(2,3) - 3*zeta(2,4) - zeta(3,4) - zeta(4,4))*
     &       Dval(DD00133,D2345))
	E(EE1133) = -(dabbr62*eta(1)) - 
     &   4*dabbr47*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &      zeta(3,4)) - eta(3)*Dval(DD1122,D1245) - 
     &   eta(5)*Dval(DD1133,D1234) + 
     &   4*(-2*dabbr46*(zeta(2,3) + zeta(3,3) + zeta(3,4)) + 
     &      zeta(1,1)*Dval(DD0022,D2345) + 
     &      zeta(1,2)*Dval(DD0022,D2345) - 
     &      zeta(1,3)*Dval(DD0022,D2345) + 
     &      zeta(1,4)*Dval(DD0022,D2345) + 
     &      zeta(2,3)*Dval(DD00112,D1245) + 
     &      zeta(3,4)*Dval(DD00113,D1234) + 
     &      zeta(1,2)*Dval(DD00122,D1245) + 
     &      zeta(1,4)*Dval(DD00133,D1234) + 
     &      zeta(1,1)*Dval(DD00122,D2345) + 
     &      zeta(1,2)*Dval(DD00122,D2345) - 
     &      zeta(1,3)*Dval(DD00122,D2345) + 
     &      zeta(1,4)*Dval(DD00122,D2345) + 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,4) - zeta(2,3) - 
     &         zeta(3,3) - zeta(3,4))*Dval(DD00222,D2345) + 
     &      (zeta(1,1) + zeta(1,2) - zeta(1,3) + zeta(1,4))*
     &       Dval(DD00223,D2345))
	E(EE1134) = -(dabbr63*eta(1)) - 
     &   2*dabbr66*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &      zeta(3,4)) - 2*dabbr64*
     &    (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4)) - 
     &   eta(3)*Dval(DD1123,D1245) + 
     &   2*(2*dabbr65*(zeta(1,1) + zeta(1,2) - zeta(2,3) - 
     &         zeta(2,4) - zeta(3,3) - 2*zeta(3,4) - zeta(4,4)) + 
     &      zeta(2,4)*Dval(DD00112,D1245) + 
     &      zeta(2,3)*Dval(DD00113,D1245) + 
     &      zeta(3,3)*Dval(DD00113,D1235) + 
     &      zeta(4,4)*Dval(DD00113,D1234) + 
     &      2*zeta(1,2)*Dval(DD00123,D1245) + 
     &      (2*zeta(1,1) + 2*zeta(1,2) + zeta(1,3) - zeta(2,3) - 
     &         2*zeta(2,4) - zeta(3,3) - 3*zeta(3,4) - 2*zeta(4,4))
     &        *Dval(DD00223,D2345) + 
     &      (2*zeta(1,1) + 2*zeta(1,2) + zeta(1,4) - 2*zeta(2,3) - 
     &         zeta(2,4) - 2*zeta(3,3) - 3*zeta(3,4) - zeta(4,4))*
     &       Dval(DD00233,D2345))
	E(EE1144) = -(dabbr67*eta(1)) - 
     &   4*dabbr51*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &      zeta(4,4)) - eta(3)*Dval(DD1133,D1245) - 
     &   eta(4)*Dval(DD1133,D1235) + 
     &   4*(-2*dabbr52*(zeta(2,4) + zeta(3,4) + zeta(4,4)) + 
     &      zeta(1,1)*Dval(DD0033,D2345) + 
     &      zeta(1,2)*Dval(DD0033,D2345) + 
     &      zeta(1,3)*Dval(DD0033,D2345) - 
     &      zeta(1,4)*Dval(DD0033,D2345) + 
     &      zeta(2,4)*Dval(DD00113,D1245) + 
     &      zeta(3,4)*Dval(DD00113,D1235) + 
     &      zeta(1,2)*Dval(DD00133,D1245) + 
     &      zeta(1,3)*Dval(DD00133,D1235) + 
     &      zeta(1,1)*Dval(DD00133,D2345) + 
     &      zeta(1,2)*Dval(DD00133,D2345) + 
     &      zeta(1,3)*Dval(DD00133,D2345) - 
     &      zeta(1,4)*Dval(DD00133,D2345) + 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,3) - zeta(1,4))*
     &       Dval(DD00233,D2345) + 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,3) - zeta(2,4) - 
     &         zeta(3,4) - zeta(4,4))*Dval(DD00333,D2345))
	E(EE1222) = dabbr68*eta(1) - 
     &   eta(4)*Dval(DD1222,D1235) - 
     &   eta(5)*Dval(DD1222,D1234) + 
     &   2*(3*dabbr43*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + zeta(1,1)*Dval(DD00111,D1345) + 
     &      3*zeta(2,3)*Dval(DD00122,D1235) + 
     &      3*zeta(2,4)*Dval(DD00122,D1234) + 
     &      zeta(1,3)*Dval(DD00222,D1235) + 
     &      zeta(1,4)*Dval(DD00222,D1234) - 
     &      (zeta(1,1) - 2*zeta(1,2) + zeta(1,3) + zeta(1,4) - 
     &         3*(zeta(2,2) + zeta(2,3) + zeta(2,4)))*
     &       Dval(DD00111,D2345))
	E(EE1223) = dabbr69*eta(1) - 
     &   eta(5)*Dval(DD1223,D1234) + 
     &   2*(2*dabbr71*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + 
     &      dabbr70*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + zeta(1,1)*Dval(DD00112,D1345) + 
     &      zeta(3,3)*Dval(DD00122,D1235) + 
     &      zeta(3,4)*Dval(DD00122,D1234) + 
     &      2*zeta(2,4)*Dval(DD00123,D1234) + 
     &      zeta(1,4)*Dval(DD00223,D1234) + 
     &      (-zeta(1,1) + zeta(1,2) - zeta(1,4) + 2*zeta(2,2) + 
     &         3*zeta(2,3) + 2*zeta(2,4) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD00112,D2345))
	E(EE1224) = dabbr72*eta(1) - 
     &   eta(4)*Dval(DD1223,D1235) + 
     &   2*(2*dabbr74*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + 
     &      dabbr73*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + zeta(1,1)*Dval(DD00113,D1345) + 
     &      zeta(3,4)*Dval(DD00122,D1235) + 
     &      zeta(4,4)*Dval(DD00122,D1234) + 
     &      2*zeta(2,3)*Dval(DD00123,D1235) + 
     &      zeta(1,3)*Dval(DD00223,D1235) + 
     &      (-zeta(1,1) + zeta(1,2) - zeta(1,3) + 2*zeta(2,2) + 
     &         2*zeta(2,3) + 3*zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD00113,D2345))
	E(EE1233) = dabbr75*eta(1) - 
     &   eta(5)*Dval(DD1233,D1234) + 
     &   2*(dabbr77*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + 
     &      2*dabbr76*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + zeta(1,1)*Dval(DD00122,D1345) + 
     &      zeta(2,2)*Dval(DD00122,D1245) + 
     &      2*zeta(3,4)*Dval(DD00123,D1234) + 
     &      zeta(2,4)*Dval(DD00133,D1234) + 
     &      zeta(1,4)*Dval(DD00233,D1234) + 
     &      (-zeta(1,1) + zeta(1,3) - zeta(1,4) + zeta(2,2) + 
     &         3*zeta(2,3) + zeta(2,4) + 2*(zeta(3,3) + zeta(3,4)))
     &        *Dval(DD00122,D2345))
	E(EE1234) = dabbr78*eta(1) + 
     &   2*(dabbr81*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + 
     &      dabbr80*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + 
     &      dabbr79*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + zeta(1,1)*Dval(DD00123,D1345) + 
     &      zeta(2,2)*Dval(DD00123,D1245) + 
     &      zeta(3,3)*Dval(DD00123,D1235) + 
     &      zeta(4,4)*Dval(DD00123,D1234) + 
     &      (-zeta(1,1) + zeta(2,2) + 2*zeta(2,3) + 2*zeta(2,4) + 
     &         zeta(3,3) + 2*zeta(3,4) + zeta(4,4))*
     &       Dval(DD00123,D2345))
	E(EE1244) = dabbr82*eta(1) - 
     &   eta(4)*Dval(DD1233,D1235) + 
     &   2*(dabbr84*(zeta(1,2) + zeta(2,2) + zeta(2,3) + 
     &         zeta(2,4)) + 
     &      2*dabbr83*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + 2*zeta(3,4)*Dval(DD00123,D1235) + 
     &      zeta(1,1)*Dval(DD00133,D1345) + 
     &      zeta(2,2)*Dval(DD00133,D1245) + 
     &      zeta(2,3)*Dval(DD00133,D1235) + 
     &      zeta(1,3)*Dval(DD00233,D1235) + 
     &      (-zeta(1,1) - zeta(1,3) + zeta(1,4) + zeta(2,2) + 
     &         zeta(2,3) + 3*zeta(2,4) + 2*(zeta(3,4) + zeta(4,4)))
     &        *Dval(DD00133,D2345))
	E(EE1333) = dabbr85*eta(1) - 
     &   eta(3)*Dval(DD1222,D1245) - 
     &   eta(5)*Dval(DD1333,D1234) + 
     &   2*(3*dabbr46*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + 3*zeta(2,3)*Dval(DD00122,D1245) + 
     &      3*zeta(3,4)*Dval(DD00133,D1234) + 
     &      zeta(1,1)*Dval(DD00222,D1345) + 
     &      zeta(1,2)*Dval(DD00222,D1245) + 
     &      zeta(1,4)*Dval(DD00333,D1234) - 
     &      (zeta(1,1) + zeta(1,2) - 2*zeta(1,3) + zeta(1,4) - 
     &         3*(zeta(2,3) + zeta(3,3) + zeta(3,4)))*
     &       Dval(DD00222,D2345))
	E(EE1334) = dabbr86*eta(1) - 
     &   eta(3)*Dval(DD1223,D1245) + 
     &   2*(2*dabbr88*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + 
     &      dabbr87*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + zeta(2,4)*Dval(DD00122,D1245) + 
     &      2*zeta(2,3)*Dval(DD00123,D1245) + 
     &      zeta(4,4)*Dval(DD00133,D1234) + 
     &      zeta(1,1)*Dval(DD00223,D1345) + 
     &      zeta(1,2)*Dval(DD00223,D1245) + 
     &      (-zeta(1,1) - zeta(1,2) + zeta(1,3) + 2*zeta(2,3) + 
     &         zeta(2,4) + 2*zeta(3,3) + 3*zeta(3,4) + zeta(4,4))*
     &       Dval(DD00223,D2345))
	E(EE1344) = dabbr89*eta(1) - 
     &   eta(3)*Dval(DD1233,D1245) + 
     &   2*(dabbr91*(zeta(1,3) + zeta(2,3) + zeta(3,3) + 
     &         zeta(3,4)) + 
     &      2*dabbr90*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + 2*zeta(2,4)*Dval(DD00123,D1245) + 
     &      zeta(2,3)*Dval(DD00133,D1245) + 
     &      zeta(3,3)*Dval(DD00133,D1235) + 
     &      zeta(1,1)*Dval(DD00233,D1345) + 
     &      zeta(1,2)*Dval(DD00233,D1245) + 
     &      (-zeta(1,1) - zeta(1,2) + zeta(1,4) + zeta(2,3) + 
     &         2*zeta(2,4) + zeta(3,3) + 3*zeta(3,4) + 2*zeta(4,4))
     &        *Dval(DD00233,D2345))
	E(EE1444) = dabbr92*eta(1) - 
     &   eta(3)*Dval(DD1333,D1245) - 
     &   eta(4)*Dval(DD1333,D1235) + 
     &   2*(3*dabbr52*(zeta(1,4) + zeta(2,4) + zeta(3,4) + 
     &         zeta(4,4)) + 3*zeta(2,4)*Dval(DD00133,D1245) + 
     &      3*zeta(3,4)*Dval(DD00133,D1235) + 
     &      zeta(1,1)*Dval(DD00333,D1345) + 
     &      zeta(1,2)*Dval(DD00333,D1245) + 
     &      zeta(1,3)*Dval(DD00333,D1235) - 
     &      (zeta(1,1) + zeta(1,2) + zeta(1,3) - 2*zeta(1,4) - 
     &         3*(zeta(2,4) + zeta(3,4) + zeta(4,4)))*
     &       Dval(DD00333,D2345))
	E(EE2222) = -(eta(2)*Dval(DD1111,D1345)) - 
     &   eta(4)*Dval(DD2222,D1235) - 
     &   eta(5)*Dval(DD2222,D1234) - 
     &   eta(1)*Dval(DD1111,D2345) + 
     &   8*zeta(1,2)*Dval(DD00111,D1345) + 
     &   8*zeta(2,3)*Dval(DD00222,D1235) + 
     &   8*zeta(2,4)*Dval(DD00222,D1234) - 
     &   8*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD00111,D2345)
	E(EE2223) = -(eta(2)*Dval(DD1112,D1345)) - 
     &   eta(5)*Dval(DD2223,D1234) - 
     &   eta(1)*Dval(DD1112,D2345) + 
     &   2*zeta(1,3)*Dval(DD00111,D1345) + 
     &   6*zeta(1,2)*Dval(DD00112,D1345) + 
     &   2*zeta(3,3)*Dval(DD00222,D1235) + 
     &   2*zeta(3,4)*Dval(DD00222,D1234) + 
     &   6*zeta(2,4)*Dval(DD00223,D1234) - 
     &   2*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD00111,D2345) - 
     &   6*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD00112,D2345)
	E(EE2224) = -(eta(2)*Dval(DD1113,D1345)) - 
     &   eta(4)*Dval(DD2223,D1235) - 
     &   eta(1)*Dval(DD1113,D2345) + 
     &   2*zeta(1,4)*Dval(DD00111,D1345) + 
     &   6*zeta(1,2)*Dval(DD00113,D1345) + 
     &   2*zeta(3,4)*Dval(DD00222,D1235) + 
     &   2*zeta(4,4)*Dval(DD00222,D1234) + 
     &   6*zeta(2,3)*Dval(DD00223,D1235) - 
     &   2*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD00111,D2345) - 
     &   6*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD00113,D2345)
	E(EE2233) = -(eta(2)*Dval(DD1122,D1345)) - 
     &   eta(5)*Dval(DD2233,D1234) - 
     &   eta(1)*Dval(DD1122,D2345) + 
     &   4*(zeta(1,3)*Dval(DD00112,D1345) + 
     &      zeta(1,2)*Dval(DD00122,D1345) + 
     &      zeta(3,4)*Dval(DD00223,D1234) + 
     &      zeta(2,4)*Dval(DD00233,D1234) - 
     &      (zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD00112,D2345) - 
     &      (zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       Dval(DD00122,D2345))
	E(EE2234) = -(eta(2)*Dval(DD1123,D1345)) - 
     &   eta(1)*Dval(DD1123,D2345) + 
     &   2*(zeta(1,4)*Dval(DD00112,D1345) + 
     &      zeta(1,3)*Dval(DD00113,D1345) + 
     &      2*zeta(1,2)*Dval(DD00123,D1345) + 
     &      zeta(3,3)*Dval(DD00223,D1235) + 
     &      zeta(4,4)*Dval(DD00223,D1234) - 
     &      (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD00112,D2345) - 
     &      (zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD00113,D2345) - 
     &      2*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       Dval(DD00123,D2345))
	E(EE2244) = -(eta(2)*Dval(DD1133,D1345)) - 
     &   eta(4)*Dval(DD2233,D1235) - 
     &   eta(1)*Dval(DD1133,D2345) + 
     &   4*(zeta(1,4)*Dval(DD00113,D1345) + 
     &      zeta(1,2)*Dval(DD00133,D1345) + 
     &      zeta(3,4)*Dval(DD00223,D1235) + 
     &      zeta(2,3)*Dval(DD00233,D1235) - 
     &      (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD00113,D2345) - 
     &      (zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       Dval(DD00133,D2345))
	E(EE2333) = -(eta(2)*Dval(DD1222,D1345)) - 
     &   eta(5)*Dval(DD2333,D1234) - 
     &   eta(1)*Dval(DD1222,D2345) + 
     &   6*zeta(1,3)*Dval(DD00122,D1345) + 
     &   2*zeta(1,2)*Dval(DD00222,D1345) + 
     &   2*zeta(2,2)*Dval(DD00222,D1245) + 
     &   6*zeta(3,4)*Dval(DD00233,D1234) + 
     &   2*zeta(2,4)*Dval(DD00333,D1234) - 
     &   6*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD00122,D2345) - 
     &   2*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD00222,D2345)
	E(EE2334) = -(eta(2)*Dval(DD1223,D1345)) - 
     &   eta(1)*Dval(DD1223,D2345) + 
     &   2*(zeta(1,4)*Dval(DD00122,D1345) + 
     &      2*zeta(1,3)*Dval(DD00123,D1345) + 
     &      zeta(1,2)*Dval(DD00223,D1345) + 
     &      zeta(2,2)*Dval(DD00223,D1245) + 
     &      zeta(4,4)*Dval(DD00233,D1234) - 
     &      (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD00122,D2345) - 
     &      2*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD00123,D2345) - 
     &      (zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       Dval(DD00223,D2345))
	E(EE2344) = -(eta(2)*Dval(DD1233,D1345)) - 
     &   eta(1)*Dval(DD1233,D2345) + 
     &   2*(2*zeta(1,4)*Dval(DD00123,D1345) + 
     &      zeta(1,3)*Dval(DD00133,D1345) + 
     &      zeta(1,2)*Dval(DD00233,D1345) + 
     &      zeta(2,2)*Dval(DD00233,D1245) + 
     &      zeta(3,3)*Dval(DD00233,D1235) - 
     &      2*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD00123,D2345) - 
     &      (zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD00133,D2345) - 
     &      (zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &       Dval(DD00233,D2345))
	E(EE2444) = -(eta(2)*Dval(DD1333,D1345)) - 
     &   eta(4)*Dval(DD2333,D1235) - 
     &   eta(1)*Dval(DD1333,D2345) + 
     &   6*zeta(1,4)*Dval(DD00133,D1345) + 
     &   6*zeta(3,4)*Dval(DD00233,D1235) + 
     &   2*zeta(1,2)*Dval(DD00333,D1345) + 
     &   2*zeta(2,2)*Dval(DD00333,D1245) + 
     &   2*zeta(2,3)*Dval(DD00333,D1235) - 
     &   6*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD00133,D2345) - 
     &   2*(zeta(1,2) + zeta(2,2) + zeta(2,3) + zeta(2,4))*
     &    Dval(DD00333,D2345)
	E(EE3333) = -(eta(2)*Dval(DD2222,D1345)) - 
     &   eta(3)*Dval(DD2222,D1245) - 
     &   eta(5)*Dval(DD3333,D1234) - 
     &   eta(1)*Dval(DD2222,D2345) + 
     &   8*zeta(1,3)*Dval(DD00222,D1345) + 
     &   8*zeta(2,3)*Dval(DD00222,D1245) + 
     &   8*zeta(3,4)*Dval(DD00333,D1234) - 
     &   8*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD00222,D2345)
	E(EE3334) = -(eta(2)*Dval(DD2223,D1345)) - 
     &   eta(3)*Dval(DD2223,D1245) - 
     &   eta(1)*Dval(DD2223,D2345) + 
     &   2*zeta(1,4)*Dval(DD00222,D1345) + 
     &   2*zeta(2,4)*Dval(DD00222,D1245) + 
     &   6*zeta(1,3)*Dval(DD00223,D1345) + 
     &   6*zeta(2,3)*Dval(DD00223,D1245) + 
     &   2*zeta(4,4)*Dval(DD00333,D1234) - 
     &   2*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD00222,D2345) - 
     &   6*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD00223,D2345)
	E(EE3344) = -(eta(2)*Dval(DD2233,D1345)) - 
     &   eta(3)*Dval(DD2233,D1245) - 
     &   eta(1)*Dval(DD2233,D2345) + 
     &   4*(zeta(1,4)*Dval(DD00223,D1345) + 
     &      zeta(2,4)*Dval(DD00223,D1245) + 
     &      zeta(1,3)*Dval(DD00233,D1345) + 
     &      zeta(2,3)*Dval(DD00233,D1245) - 
     &      (zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &       Dval(DD00223,D2345) - 
     &      (zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &       Dval(DD00233,D2345))
	E(EE3444) = -(eta(2)*Dval(DD2333,D1345)) - 
     &   eta(3)*Dval(DD2333,D1245) - 
     &   eta(1)*Dval(DD2333,D2345) + 
     &   6*zeta(1,4)*Dval(DD00233,D1345) + 
     &   6*zeta(2,4)*Dval(DD00233,D1245) + 
     &   2*zeta(1,3)*Dval(DD00333,D1345) + 
     &   2*zeta(2,3)*Dval(DD00333,D1245) + 
     &   2*zeta(3,3)*Dval(DD00333,D1235) - 
     &   6*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD00233,D2345) - 
     &   2*(zeta(1,3) + zeta(2,3) + zeta(3,3) + zeta(3,4))*
     &    Dval(DD00333,D2345)
	E(EE4444) = -(eta(2)*Dval(DD3333,D1345)) - 
     &   eta(3)*Dval(DD3333,D1245) - 
     &   eta(4)*Dval(DD3333,D1235) - 
     &   eta(1)*Dval(DD3333,D2345) + 
     &   8*zeta(1,4)*Dval(DD00333,D1345) + 
     &   8*zeta(2,4)*Dval(DD00333,D1245) + 
     &   8*zeta(3,4)*Dval(DD00333,D1235) - 
     &   8*(zeta(1,4) + zeta(2,4) + zeta(3,4) + zeta(4,4))*
     &    Dval(DD00333,D2345)

	if( dump ) call XDumpCoeff(5, E)
	end

