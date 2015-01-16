	double complex Cval(13, CSTORE)
	double precision Cpara(6, CSTORE)
	common /cpave/ Cval, Cpara

	double complex Dval(46, DSTORE)
	double precision Dpara(10, DSTORE)
	common /dpave/ Dval, Dpara

	integer Cptr, Dptr
	common /cdptrs/ Cptr, Dptr

#ifndef IN_TOOLS_F
	double complex A0, B0, DB0, C0, D0
	double complex B1, DB1, B00, DB00, B11, DB11
	double complex C0i, D0i
	integer Cget, Dget
	external A0, B0, DB0, C0, D0
	external B1, DB1, B00, DB00, B11, DB11
	external C0i, D0i
	external Cget, Dget
#endif

