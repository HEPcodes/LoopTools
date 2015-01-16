* Dump.F
* dumps the parameters and coefficients on screen
* this file is part of LoopTools
* last modified 10 Sep 12 th

#include "externals.h"
#include "types.h"

#include "defs.h"


	subroutine XDumpPara(npoint, para, origin)
	implicit none
	integer npoint
	DVAR para(1,*)
	character*(*) origin

#include "lt.h"

	integer i

	character*6 paraname(Pee,2:5)
	common /ltparanames/ paraname

	integer npara(2:5)
	data npara /Pbb, Pcc, Pdd, Pee/

#ifdef COMPLEXPARA
	if( len(origin) .gt. 1 ) print *, origin, "C", serial
#else
	if( len(origin) .gt. 1 ) print *, origin, serial
#endif

	do i = npoint + 1, npara(npoint)
	  print *, "  ", paraname(i,npoint), "=", para(1,i)
	enddo
	do i = 1, npoint
	  print *, "  ", paraname(i,npoint), "=", para(1,i)
	enddo

	call flush(6)
	end

************************************************************************

	subroutine XDumpCoeff(npoint, coeff)
	implicit none
	integer npoint
	ComplexType coeff(*)

#include "lt.h"

	integer i

	character*8 coeffname(Nee,2:5)
	common /ltcoeffnames/ coeffname

	integer ncoeff(2:5)
	data ncoeff /Nbb, Ncc, Ndd, Nee/

	do i = 1, ncoeff(npoint)
	  print *, coeffname(i,npoint), "=", coeff(i)
	enddo

	print *, "===================================================="

	call flush(6)
	end

