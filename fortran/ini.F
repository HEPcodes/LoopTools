* ini.F
* routines for initializing and setting some parameters
* this file is part of LoopTools
* last modified 22 Feb 00 th


	subroutine setmudim(newmudim2)
	implicit none
	double precision newmudim2

#include "../ff/ff.h"

	double complex Ccache(1)
	common /Cbase/ Ccache

	double complex Dcache(1)
	common /Dbase/ Dcache

	mudim2 = newmudim2

	call setcachelast(Ccache, 0)
	call setcachelast(Dcache, 0)
	end

************************************************************************

	double precision function getmudim()
	implicit none

#include "../ff/ff.h"

	getmudim = mudim2
	end

************************************************************************

	subroutine setdelta(newdivergence)
	implicit none
	double precision newdivergence

#include "../ff/ff.h"

	double complex Ccache(1)
	common /Cbase/ Ccache

	double complex Dcache(1)
	common /Dbase/ Dcache

	divergence = newdivergence

	call setcachelast(Ccache, 0)
	call setcachelast(Dcache, 0)
	end

************************************************************************

	double precision function getdelta()
	implicit none

#include "../ff/ff.h"

	getdelta = divergence
	end

************************************************************************

	subroutine setlambda(newlambda2)
	implicit none
	double precision newlambda2

#include "../ff/ff.h"

	double complex Ccache(1)
	common /Cbase/ Ccache

	double complex Dcache(1)
	common /Dbase/ Dcache

	lambda2 = newlambda2

	call setcachelast(Ccache, 0)
	call setcachelast(Dcache, 0)
	end

************************************************************************

	double precision function getlambda()
	implicit none

#include "../ff/ff.h"

	getlambda = lambda2
	end

