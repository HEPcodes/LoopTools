:Begin:
:Function:	fa0
:Pattern:	A0[m_?doubleQ]
:Arguments:	{N[m]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fb0
:Pattern:	B0[p_?doubleQ, m1_?doubleQ, m2_?doubleQ]
:Arguments:	{N[p], N[m1], N[m2]}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fdb0
:Pattern:	DB0[p_?doubleQ, m1_?doubleQ, m2_?doubleQ]
:Arguments:	{N[p], N[m1], N[m2]}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fb1
:Pattern:	B1[p_?doubleQ, m1_?doubleQ, m2_?doubleQ]
:Arguments:	{N[p], N[m1], N[m2]}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fdb1
:Pattern:	DB1[p_?doubleQ, m1_?doubleQ, m2_?doubleQ]
:Arguments:	{N[p], N[m1], N[m2]}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fb00
:Pattern:	B00[p_?doubleQ, m1_?doubleQ, m2_?doubleQ]
:Arguments:	{N[p], N[m1], N[m2]}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fdb00
:Pattern:	DB00[p_?doubleQ, m1_?doubleQ, m2_?doubleQ]
:Arguments:	{N[p], N[m1], N[m2]}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fb11
:Pattern:	B11[p_?doubleQ, m1_?doubleQ, m2_?doubleQ]
:Arguments:	{N[p], N[m1], N[m2]}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fdb11
:Pattern:	DB11[p_?doubleQ, m1_?doubleQ, m2_?doubleQ]
:Arguments:	{N[p], N[m1], N[m2]}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fc0
:Pattern:	C0[p1_?doubleQ, p2_?doubleQ, p1p2_?doubleQ,
		  m1_?doubleQ, m2_?doubleQ, m3_?doubleQ]
:Arguments:	{N[p1], N[p2], N[p1p2], N[m1], N[m2], N[m3]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fc0i
:Pattern:	C0i[id_Symbol,
		  p1_?doubleQ, p2_?doubleQ, p1p2_?doubleQ,
		  m1_?doubleQ, m2_?doubleQ, m3_?doubleQ]
:Arguments:	{id, N[p1], N[p2], N[p1p2], N[m1], N[m2], N[m3]}
:ArgumentTypes:	{Symbol, Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fcget
:Pattern:	Cget[p1_?doubleQ, p2_?doubleQ, p1p2_?doubleQ,
		  m1_?doubleQ, m2_?doubleQ, m3_?doubleQ]
:Arguments:	{N[p1], N[p2], N[p1p2], N[m1], N[m2], N[m3]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fd0
:Pattern:	D0[p1_?doubleQ, p2_?doubleQ, p3_?doubleQ, p4_?doubleQ,
		  p1p2_?doubleQ, p2p3_?doubleQ,
		  m1_?doubleQ, m2_?doubleQ, m3_?doubleQ, m4_?doubleQ]
:Arguments:	{N[p1], N[p2], N[p3], N[p4], N[p1p2], N[p2p3],
		  N[m1], N[m2], N[m3], N[m4]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fd0i
:Pattern:	D0i[id_Symbol,
		  p1_?doubleQ, p2_?doubleQ, p3_?doubleQ, p4_?doubleQ,
		  p1p2_?doubleQ, p2p3_?doubleQ,
		  m1_?doubleQ, m2_?doubleQ, m3_?doubleQ, m4_?doubleQ]
:Arguments:	{id, N[p1], N[p2], N[p3], N[p4], N[p1p2], N[p2p3],
		  N[m1], N[m2], N[m3], N[m4]}
:ArgumentTypes:	{Symbol, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fdget
:Pattern:	Dget[p1_?doubleQ, p2_?doubleQ, p3_?doubleQ, p4_?doubleQ,
		  p1p2_?doubleQ, p2p3_?doubleQ,
		  m1_?doubleQ, m2_?doubleQ, m3_?doubleQ, m4_?doubleQ]
:Arguments:	{N[p1], N[p2], N[p3], N[p4], N[p1p2], N[p2p3],
		  N[m1], N[m2], N[m3], N[m4]}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fsetmudim
:Pattern:	SetMudim[newmudim_?doubleQ]
:Arguments:	{N[newmudim]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fgetmudim
:Pattern:	GetMudim[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Real
:End:

:Begin:
:Function:	fsetlambda
:Pattern:	SetLambda[newlambda_?doubleQ]
:Arguments:	{N[newlambda]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fgetlambda
:Pattern:	GetLambda[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Real
:End:

:Begin:
:Function:	fsetdelta
:Pattern:	SetDelta[newdelta_?doubleQ]
:Arguments:	{N[newdelta]}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fgetdelta
:Pattern:	GetDelta[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Real
:End:

:Begin:
:Function:	fgetcachelast
:Pattern:	GetCacheLast[buffer_Symbol]
:Arguments:	{buffer}
:ArgumentTypes:	{Symbol}
:ReturnType:	Integer
:End:

:Begin:
:Function:	fsetcachelast
:Pattern:	SetCacheLast[buffer_Symbol, offset_Integer]
:Arguments:	{buffer, offset}
:ArgumentTypes:	{Symbol, Integer}
:ReturnType:	Manual
:End:


:Evaluate:	doubleQ = IntegerQ[Floor[#]]&

:Evaluate:	A0[0] = 0

:Evaluate:	LoopTools::undef = "Identifier `1` unknown."

:Evaluate:	Derivative[1, 0, 0][B0] = DB0

:Evaluate:	Derivative[1, 0, 0][B1] = DB1

:Evaluate:	Derivative[1, 0, 0][B00] = DB00

:Evaluate:	Derivative[1, 0, 0][B11] = DB11

:Evaluate:	PaVe[i__Integer, p_List, {m1_, m2_, m3_, m4_}] := D0i[
		  StringJoin["dd", ToString/@ Sort[{i}]]//ToExpression,
		  Sequence@@ p, m1, m2, m3, m4]

:Evaluate:	PaVe[i__Integer, p_List, {m1_, m2_, m3_}] := C0i[
		  StringJoin["cc", ToString/@ Sort[{i}]]//ToExpression,
		  Sequence@@ p, m1, m2, m3]


/*
LoopTools.tm
provides the scalar and tensor one-loop functions in Mathematica
this file is part of LoopTools
last modified 7 Jun 00 th
*/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "mathlink.h"
#include "../include/ltproto.h"

#define fgetmudim getmudim
#define fgetlambda getlambda
#define fgetdelta getdelta

#ifndef MLCONST
#define MLCONST
#endif


int savestdout;

char *clookup[] = {
  "cc0", "cc1", "cc2", "cc00", "cc11", "cc12", "cc22",
  "cc001", "cc002", "cc111", "cc112", "cc122", "cc222" };

char *dlookup[] = {
  "dd0", "dd1", "dd2", "dd3", "dd00", "dd11", "dd12", "dd13", "dd22",
  "dd23", "dd33", "dd001", "dd002", "dd003", "dd111", "dd112",
  "dd113", "dd122", "dd123", "dd133", "dd222", "dd223", "dd233",
  "dd333", "dd0000", "dd0011", "dd0012", "dd0013", "dd0022",
  "dd0023", "dd0033", "dd1111", "dd1112", "dd1113", "dd1122",
  "dd1123", "dd1133", "dd1222", "dd1223", "dd1233", "dd1333",
  "dd2222", "dd2223", "dd2233", "dd2333", "dd3333" };


void return_complex(dcomplex c)
{
  dup2(savestdout, 1);
  fflush(stderr);
  if(c.i == 0.) MLPutReal(stdlink, c.r);
  else {
    MLPutFunction(stdlink, "Complex", 2);
    MLPutReal(stdlink, c.r);
    MLPutReal(stdlink, c.i);
    MLEndPacket(stdlink);
  }
}


void report_error(MLCONST char *id)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", 2);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "LoopTools");
  MLPutString(stdlink, "undef");
  MLPutSymbol(stdlink, id);
  MLPutSymbol(stdlink, "$Failed");
  MLEndPacket(stdlink);
}


void fa0(double m)
{
  dcomplex result;

  dup2(2, 1);
  A0(&result, &m);
  return_complex(result);
}


void fb0(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  B0(&result, &p, &m1, &m2);
  return_complex(result);
}


void fdb0(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  DB0(&result, &p, &m1, &m2);
  return_complex(result);
}


void fb1(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  B1(&result, &p, &m1, &m2);
  return_complex(result);
}


void fdb1(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  DB1(&result, &p, &m1, &m2);
  return_complex(result);
}


void fb00(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  B00(&result, &p, &m1, &m2);
  return_complex(result);
}


void fdb00(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  DB00(&result, &p, &m1, &m2);
  return_complex(result);
}


void fb11(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  B11(&result, &p, &m1, &m2);
  return_complex(result);
}


void fdb11(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  DB11(&result, &p, &m1, &m2);
  return_complex(result);
}


void fc0(double p1, double p2, double p1p2,
  double m1, double m2, double m3)
{
  dcomplex result;

  dup2(2, 1);
  C0(&result, &p1, &p2, &p1p2, &m1, &m2, &m3);
  return_complex(result);
}


void fd0(double p1, double p2, double p3, double p4,
  double p1p2, double p2p3,
  double m1, double m2, double m3, double m4)
{
  dcomplex result;

  dup2(2, 1);
  D0(&result, &p1, &p2, &p3, &p4, &p1p2, &p2p3, &m1, &m2, &m3, &m4);
  return_complex(result);
}


void fc0i(MLCONST char *id,
  double p1, double p2, double p1p2,
  double m1, double m2, double m3)
{
  dcomplex result;
  int id_ = 1;
  char **lp;

  for(lp = clookup; lp < &clookup[13]; ++lp, ++id_)
    if(strcmp(id, *lp) == 0) goto found;
  report_error(id);
  return;

found:
  dup2(2, 1);
  C0i(&result, &id_, &p1, &p2, &p1p2, &m1, &m2, &m3);
  return_complex(result);
}


void fcget(double p1, double p2, double p1p2,
  double m1, double m2, double m3)
{
  int offset;
  dcomplex *cp;
  char **lp;

  dup2(2, 1);
  offset = Cget(&p1, &p2, &p1p2, &m1, &m2, &m3);
  dup2(savestdout, 1);
  fflush(stderr);

  MLPutFunction(stdlink, "List", 13);
  cp = &Ccache[offset];
  for(lp = clookup; lp < &clookup[13]; ++lp, ++cp) {
    MLPutFunction(stdlink, "Rule", 2);
    MLPutSymbol(stdlink, *lp);
    if(cp->i == 0.) MLPutReal(stdlink, cp->r);
    else {
      MLPutFunction(stdlink, "Complex", 2);
      MLPutReal(stdlink, cp->r);
      MLPutReal(stdlink, cp->i);
    }
  }
  MLEndPacket(stdlink);
}


void fd0i(MLCONST char *id,
  double p1, double p2, double p3, double p4,
  double p1p2, double p2p3,
  double m1, double m2, double m3, double m4)
{
  dcomplex result;
  int id_ = 1;
  char **lp;

  for(lp = dlookup; lp < &dlookup[46]; ++lp, ++id_)
    if(strcmp(id, *lp) == 0) goto found;
  report_error(id);
  return;

found:
  dup2(2, 1);
  D0i(&result, &id_, &p1, &p2, &p3, &p4, &p1p2, &p2p3,
    &m1, &m2, &m3, &m4);
  return_complex(result);
}


void fdget(double p1, double p2, double p3, double p4,
  double p1p2, double p2p3,
  double m1, double m2, double m3, double m4)
{
  int offset;
  dcomplex *dp;
  char **lp;

  dup2(2, 1);
  offset = Dget(&p1, &p2, &p3, &p4, &p1p2, &p2p3,
    &m1, &m2, &m3, &m4);
  dup2(savestdout, 1);
  fflush(stderr);

  MLPutFunction(stdlink, "List", 46);
  dp = &Dcache[offset];
  for(lp = dlookup; lp < &dlookup[46]; ++lp, ++dp) {
    MLPutFunction(stdlink, "Rule", 2);
    MLPutSymbol(stdlink, *lp);
    if(dp->i == 0.) MLPutReal(stdlink, dp->r);
    else {
      MLPutFunction(stdlink, "Complex", 2);
      MLPutReal(stdlink, dp->r);
      MLPutReal(stdlink, dp->i);
    }
  }
  MLEndPacket(stdlink);
}


void fsetmudim(double newmudim2)
{
  setmudim(&newmudim2);
  MLPutSymbol(stdlink, "Null");
  MLEndPacket(stdlink);
}


void fsetdelta(double newdivergence)
{
  setdelta(&newdivergence);
  MLPutSymbol(stdlink, "Null");
  MLEndPacket(stdlink);
}


void fsetlambda(double newlambda2)
{
  setlambda(&newlambda2);
  MLPutSymbol(stdlink, "Null");
  MLEndPacket(stdlink);
}


int fgetcachelast(MLCONST char *buffer)
{
  dcomplex *buf;

  if(strcmp(buffer, "Ccache") == 0) buf = Ccache;
  else if(strcmp(buffer, "Dcache") == 0) buf = Dcache;
  else {
    report_error(buffer);
    return;
  }
  return getcachelast(buf);
}


void fsetcachelast(MLCONST char *buffer, int offset)
{
  dcomplex *buf;

  if(strcmp(buffer, "Ccache") == 0) buf = Ccache;
  else if(strcmp(buffer, "Dcache") == 0) buf = Dcache;
  else {
    report_error(buffer);
    return;
  }
  setcachelast(buf, &offset);
  MLPutSymbol(stdlink, "Null");
  MLEndPacket(stdlink);
}


int main(int argc, char **argv)
{
  int ret;

  savestdout = dup(1);
  dup2(2, 1);
  ffini();
  dup2(savestdout, 1);
  ret = MLMain(argc, argv);
  dup2(2, 1);
  ffexi();
  dup2(savestdout, 1);
  fflush(stderr);
  return ret;
}

