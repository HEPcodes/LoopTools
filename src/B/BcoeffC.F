* BcoeffC.F
* invoke the two-point tensor coefficients
* this file is part of LoopTools
* last modified 19 Sep 12 th

#include "externals.h"
#include "types.h"

#define npoint 2
#include "defs.h"


	subroutine BcoeffC(B, para)
	implicit none
	ComplexType B(*), para(1,*)

#include "lt.h"

	integer ier(Nbb), i
	logical ini

	character*5 name(Nbb)
	data name /"bb0", "bb1", "bb11", "bb00", "bb001", "bb111",
     &    "dbb0", "dbb1", "dbb00", "dbb11"/

	if( lambda .lt. 0 ) then
	  do i = 1, Nbb
	    B(i) = 0
	  enddo
	  if( lambda .eq. -1 ) then
	    B(bb0) = 1
	    B(bb1) = -.5D0
	    B(bb00) = -(P(1) - 3*(M(1) + M(2)))/12D0
	    B(bb11) = 1/3D0
	    B(bb001) = -(P(1) - 2*M(1) - 4*M(2))/24D0
	    B(bb111) = -.25D0
	    B(dbb00) = -1/12D0
	  endif
	  return
	endif

	call BcoeffFFC(B, para, ier)

	ini = .TRUE.
	do i = 1, Nbb
	  if( ier(i) .gt. warndigits ) then
	    if( ini ) then
	      print *, "Loss of digits in BgetC for:"
	      call DumpParaC(2, para, " ")
	      ini = .FALSE.
	    endif
	    print *, name(i), " claims ", ier(i), "lost digits"
	  endif
	enddo
	end

