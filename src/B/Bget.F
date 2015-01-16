* Bget.F
* retrieve the two-point tensor coefficients
* this file is part of LoopTools
* last modified 13 Sep 12 th

#include "externals.h"
#include "types.h"

#define npoint 2
#include "defs.h"


	memindex function XBget(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex cacheindex
	external cacheindex, XBcoeff

#ifdef COMPLEXPARA
	memindex Bget
	external Bget
#endif

	DVAR para(1,Pbb)

	P(1) = p
#ifdef COMPLEXPARA
	if( abs(Im(P(1))) .gt. 0 )
     &    print *, "BgetC: Complex momenta not implemented"
#endif

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
#ifdef COMPLEXPARA
	if( abs(Im(M(1))) + abs(Im(M(2))) .eq. 0 ) then
	  XBget = Bget(p, m1, m2) - offsetC
	  return
	endif
#endif

	XBget = cacheindex(para, Bval(1,0), XBcoeff, RC*Pbb, Nbb)
	end

************************************************************************

	subroutine XBput(res, p, m1, m2)
	implicit none
	ComplexType res(*)
	DVAR p, m1, m2

#include "lt.h"

	external XBcoeff

	DVAR para(1,Pbb)

	P(1) = p
#ifdef COMPLEXPARA
	if( abs(Im(P(1))) .gt. 0 )
     &    print *, "BputC: Complex momenta not implemented"
#endif

	M(1) = m1
	if( abs(M(1)) .lt. minmass ) M(1) = 0
	M(2) = m2
	if( abs(M(2)) .lt. minmass ) M(2) = 0
#ifdef COMPLEXPARA
	if( abs(Im(M(1))) + abs(Im(M(2))) .eq. 0 ) then
	  call Bput(res, p, m1, m2)
	  return
	endif
#endif

	call cachecopy(res, para, Bval(1,0), XBcoeff, RC*Pbb, Nbb)
	end

************************************************************************

	ComplexType function XB0i(i, p, m1, m2)
	implicit none
	integer i
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	memindex b

	b = XBget(p, m1, m2)
	XB0i = Bval(i,b)
	end

************************************************************************

	ComplexType function XB0(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	XB0 = Bval(bb0,XBget(p, m1, m2))
	end

************************************************************************

	ComplexType function XB1(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	XB1 = Bval(bb1,XBget(p, m1, m2))
	end

************************************************************************

	ComplexType function XB00(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	XB00 = Bval(bb00,XBget(p, m1, m2))
	end

************************************************************************

	ComplexType function XB11(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	XB11 = Bval(bb11,XBget(p, m1, m2))
	end

************************************************************************

	ComplexType function XB001(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	XB001 = Bval(bb001,XBget(p, m1, m2))
	end

************************************************************************

	ComplexType function XB111(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	XB111 = Bval(bb111,XBget(p, m1, m2))
	end

************************************************************************

	ComplexType function XDB0(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	XDB0 = Bval(dbb0,XBget(p, m1, m2))
	end

************************************************************************

	ComplexType function XDB1(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	XDB1 = Bval(dbb1,XBget(p, m1, m2))
	end

************************************************************************

	ComplexType function XDB00(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	XDB00 = Bval(dbb00,XBget(p, m1, m2))
	end

************************************************************************

	ComplexType function XDB11(p, m1, m2)
	implicit none
	DVAR p, m1, m2

#include "lt.h"

	memindex XBget
	external XBget

	XDB11 = Bval(dbb11,XBget(p, m1, m2))
	end

