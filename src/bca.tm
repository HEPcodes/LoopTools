:Begin:
:Function:	fa0
:Pattern:	A0[m_Real]
:Arguments:	{m}
:ArgumentTypes:	{Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fb0
:Pattern:	B0[p_Real, m1_Real, m2_Real]
:Arguments:	{p, m1, m2}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fdb0
:Pattern:	DB0[p_Real, m1_Real, m2_Real]
:Arguments:	{p, m1, m2}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fc0
:Pattern:	C0[p1_Real, p2_Real, p1p2_Real, m1_Real, m2_Real, m3_Real]
:Arguments:	{p1, p2, p1p2, m1, m2, m3}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fd0
:Pattern:	D0[p1_Real, p2_Real, p3_Real, p4_Real,
		  p1p2_Real, p2p3_Real,
		  m1_Real, m2_Real, m3_Real, m4_Real]
:Arguments:	{p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fb1
:Pattern:	B1[p_Real, m1_Real, m2_Real]
:Arguments:	{p, m1, m2}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fdb1
:Pattern:	DB1[p_Real, m1_Real, m2_Real]
:Arguments:	{p, m1, m2}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fb00
:Pattern:	B00[p_Real, m1_Real, m2_Real]
:Arguments:	{p, m1, m2}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fdb00
:Pattern:	DB00[p_Real, m1_Real, m2_Real]
:Arguments:	{p, m1, m2}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fb11
:Pattern:	B11[p_Real, m1_Real, m2_Real]
:Arguments:	{p, m1, m2}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fdb11
:Pattern:	DB11[p_Real, m1_Real, m2_Real]
:Arguments:	{p, m1, m2}
:ArgumentTypes:	{Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fc0i
:Pattern:	C0i[id_Symbol, p1_Real, p2_Real, p1p2_Real,
		  m1_Real, m2_Real, m3_Real]
:Arguments:	{id, p1, p2, p1p2, m1, m2, m3}
:ArgumentTypes:	{Symbol, Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fcget
:Pattern:	Cget[p1_Real, p2_Real, p1p2_Real,
		  m1_Real, m2_Real, m3_Real]
:Arguments:	{p1, p2, p1p2, m1, m2, m3}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fd0i
:Pattern:	D0i[id_Symbol, p1_Real, p2_Real, p3_Real, p4_Real,
		  p1p2_Real, p2p3_Real,
		  m1_Real, m2_Real, m3_Real, m4_Real]
:Arguments:	{id, p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4}
:ArgumentTypes:	{Symbol, Real, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	fdget
:Pattern:	Dget[p1_Real, p2_Real, p3_Real, p4_Real,
		  p1p2_Real, p2p3_Real,
		  m1_Real, m2_Real, m3_Real, m4_Real]
:Arguments:	{p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4}
:ArgumentTypes:	{Real, Real, Real, Real, Real, Real, Real,
		  Real, Real, Real}
:ReturnType:	Manual
:End:

:Begin:
:Function:	set_mudim
:Pattern:	Mudim[mud_Real]
:Arguments:	{mud}
:ArgumentTypes:	{Real}
:ReturnType:	Real
:End:

:Begin:
:Function:	set_lambda
:Pattern:	Lambda[lambda_Real]
:Arguments:	{lambda}
:ArgumentTypes:	{Real}
:ReturnType:	Real
:End:

:Begin:
:Function:	set_delta
:Pattern:	Delta[del_Real]
:Arguments:	{del}
:ArgumentTypes:	{Real}
:ReturnType:	Real
:End:

:Begin:
:Function:	flushcache
:Pattern:	FlushCache[]
:Arguments:	{}
:ArgumentTypes:	{}
:ReturnType:	Manual
:End:

:Evaluate:	Scan[(#[0] := #[0.];
		  #[x:_Integer | _Rational] := #[x//N])&,
		  {A0, Mudim, Lambda, Delta}]
:Evaluate:	A0[0] = 0
:Evaluate:	Scan[(#[a___, 0, b___] := #[a, 0., b];
		  #[a___, x:_Integer | _Rational, b___] := #[a, x//N, b])&,
		  {B0, DB0, B1, DB1, B00, DB00, B11, DB11,
		   C0, D0, Cget, Dget}]
:Evaluate:	Scan[(#[a__, 0, b___] := #[a, 0., b];
		  #[a__, x:_Integer | _Rational, b___] := #[a, x//N, b])&,
		  {C0i, D0i}]
:Evaluate:	C0i::undef = D0i::undef = "Identifier `1` unknown"
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
	bca.tm
		this file is part of LoopTools;
		it provides the scalar and tensor
		one-loop functions in Mathematica
		last modified 15 Jul 99 th
*/


#include "mathlink.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "ffc.h"

#define MAXDIGITS 9
#define MUDIM 1.
#define LAMBDA 1.
#define DELTA 0.
#define MINMASS (1.e-5)

#ifndef MLCONST
#define MLCONST
#endif

#ifdef NO_UNDERSCORE
#define b1_ b1
#define db1_ db1
#define b00_ b00
#define db00_ db00
#define b11_ b11
#define db11_ db11
#define c0i_ c0i
#define cget_ cget
#define d0i_ d0i
#define dget_ dget
#define cpave_ cpave
#define dpave_ dpave
#define cdptrs_ cdptrs
#define cutoff_ cutoff
#endif

/* prototypes & common blocks from tensor.F */

int b1_(dcomplex *cb1,
  const double *p, const double *m1, const double *m2);
int b00_(dcomplex *cb00,
  const double *p, const double *m1, const double *m2);
int db00_(dcomplex *cdb00,
  const double *p, const double *m1, const double *m2);
int b11_(dcomplex *cb11,
  const double *p, const double *m1, const double *m2);
int db11_(dcomplex *cdb11,
  const double *p, const double *m1, const double *m2);
int c0i_(dcomplex *cc0i, const int *id,
  const double *p1, const double *p2, const double *p1p2,
  const double *m1, const double *m2, const double *m3);
int cget_(
  const double *p1, const double *p2, const double *p1p2,
  const double *m1, const double *m2, const double *m3);
int d0i_(dcomplex *cd0i, const int *id,
  const double *p1, const double *p2, const double *p3,
  const double *p4, const double *p1p2, const double *p2p3,
  const double *m1, const double *m2, const double *m3,
  const double *m4);
int dget_(
  const double *p1, const double *p2, const double *p3,
  const double *p4, const double *p1p2, const double *p2p3,
  const double *m1, const double *m2, const double *m3,
  const double *m4);

extern struct {
  dcomplex cval[13*CSTORE];
  double cpara[6*CSTORE];
} cpave_;

extern struct {
  dcomplex dval[46*DSTORE];
  double dpara[10*DSTORE];
} dpave_;

extern struct {
  int cptr, dptr;
} cdptrs_;

extern struct {
  double mudim, divergence;
} cutoff_;

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

void report_error(MLCONST char *msg, MLCONST char *id)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", 2);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, msg);
  MLPutString(stdlink, "undef");
  MLPutSymbol(stdlink, id);
  MLPutSymbol(stdlink, "$Failed");
  MLEndPacket(stdlink);
}

void fa0(double m)
{
  dcomplex result;
  int ier = 0;

  dup2(2, 1);
  ffxa0_(&result, &cutoff_.divergence, &cutoff_.mudim, &m, &ier);
  if(ier > MAXDIGITS)
    fprintf(stderr, "A0: lost %d digits (m = %lg)\n", ier, sqrt(m));
  return_complex(result);
}

void fb0(double p, double m1, double m2)
{
  dcomplex result;
  int ier = 0;

  dup2(2, 1);
  ffxb0_(&result, &cutoff_.divergence, &cutoff_.mudim, &p, &m1, &m2, &ier);
  if(ier > MAXDIGITS)
    fprintf(stderr, "B0: lost %d digits (m = %lg)\n", ier, sqrt(m1));
  return_complex(result);
}

void fdb0(double p, double m1, double m2)
{
  dcomplex result, pdb0;
  double dm;
  int ier = 0;

  dup2(2, 1);
  if(p != 0. && m1 != 0. && m2 != 0. &&
     p == (dm = sqrt(m1) - sqrt(m2), dm *= dm)) {
    result.r = ((m2 - m1)/2./dm*log(m2/m1) - 2.)/dm;
    result.i = 0.;
  }
  else {
    ffxdb0_(&result, &pdb0, &p, &m1, &m2, &ier);
    if(ier > MAXDIGITS)
      fprintf(stderr, "DB0: lost %d digits (m = %lg)\n", ier, sqrt(m1));
  }
  return_complex(result);
}

void fc0(double p1, double p2, double p1p2,
  double m1, double m2, double m3)
{
  dcomplex result;
  double xpi[6] = {m1, m2, m3, p1, p2, p1p2};
  int ier = 0;

  dup2(2, 1);
  ffxc0_(&result, xpi, &ier);
  if(ier > MAXDIGITS) {
    ier = 0;
    ffxc0r_(&result, xpi, &ier); 
    if(ier > MAXDIGITS)
      fprintf(stderr, "C0: lost %d digits (m = %lg)\n",
        ier, sqrt(m1));
  }
  return_complex(result);
}

void fd0(double p1, double p2, double p3, double p4,
  double p1p2, double p2p3,
  double m1, double m2, double m3, double m4)
{
  dcomplex result;
  double xpi[13] = {
    m1, m2, m3, m4,
    p1, p2, p3, p4, p1p2, p2p3,
    0., 0., 0.
  };
  int ier = 0;

  if(p1p2 < 0. && p2p3 < 0. && p1p2 > p2p3 &&
    m1 == m2 && m1 == m3 && m1 == m4) {  
    xpi[8] = p2p3;
    xpi[9] = p1p2;            /* numerically more stable */
  }
  dup2(2, 1);
  if(m1 + m2 + m3 + m4 <= MINMASS)
    ffxd0m0_(&result, xpi, &ier);
  else {
    ffxd0_(&result, xpi, &ier);
    if(ier > MAXDIGITS) {
      ier = 0;  
      ffxd0r_(&result, xpi, &ier);
      if(ier > MAXDIGITS)
        fprintf(stderr, "D0: lost %d digits (m = %lg)\n",
          ier, sqrt(m1));
    }
  }
  return_complex(result);
}

void fb1(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  b1_(&result, &p, &m1, &m2);
  return_complex(result);
}

void fdb1(double p, double m1, double m2)
{
  dcomplex result;
  int ier = 0;

  dup2(2, 1);
  ffxdb1_(&result, &p, &m1, &m2, &ier);
  return_complex(result);
}

void fb00(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  b00_(&result, &p, &m1, &m2);
  return_complex(result);
}

void fdb00(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  db00_(&result, &p, &m1, &m2);
  return_complex(result);
}

void fb11(double p, double m1, double m2)
{
  dcomplex result;

  dup2(2, 1);
  b11_(&result, &p, &m1, &m2);
  return_complex(result);
}

void fdb11(double p, double m1, double m2)
{
  dcomplex result;
  int ier = 0;

  dup2(2, 1);
  ffxdb11_(&result, &p, &m1, &m2, &ier);
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
  report_error("C0i", id);
  return;
found:
  dup2(2, 1);
  c0i_(&result, &id_, &p1, &p2, &p1p2, &m1, &m2, &m3);
  return_complex(result);
}

void fcget(double p1, double p2, double p1p2,
  double m1, double m2, double m3)
{
  int theC;
  dcomplex *cp;
  char **lp;

  dup2(2, 1);
  theC = cget_(&p1, &p2, &p1p2, &m1, &m2, &m3);
  dup2(savestdout, 1);
  fflush(stderr);
  MLPutFunction(stdlink, "List", 13);
  cp = &cpave_.cval[13*(theC - 1)];
  for(theC = 13, lp = clookup; theC; ++cp, --theC) {
    MLPutFunction(stdlink, "Rule", 2);
    MLPutSymbol(stdlink, *lp++);
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
  report_error("D0i", id);
  return;
found:
  dup2(2, 1);
  d0i_(&result, &id_, &p1, &p2, &p3, &p4, &p1p2, &p2p3,
    &m1, &m2, &m3, &m4);
  return_complex(result);
}

void fdget(double p1, double p2, double p3, double p4,
  double p1p2, double p2p3,
  double m1, double m2, double m3, double m4)
{
  int theD;
  dcomplex *dp;
  char **lp;

  dup2(2, 1);
  theD = dget_(&p1, &p2, &p3, &p4, &p1p2, &p2p3,
    &m1, &m2, &m3, &m4);
  dup2(savestdout, 1);
  fflush(stderr);
  MLPutFunction(stdlink, "List", 46);
  dp = &dpave_.dval[46*(theD - 1)];
  for(theD = 46, lp = dlookup; theD; ++dp, --theD) {
    MLPutFunction(stdlink, "Rule", 2);
    MLPutSymbol(stdlink, *lp++);
    if(dp->i == 0.) MLPutReal(stdlink, dp->r);
    else {
      MLPutFunction(stdlink, "Complex", 2);
      MLPutReal(stdlink, dp->r);
      MLPutReal(stdlink, dp->i);
    }
  }
  MLEndPacket(stdlink);
}

double set_mudim(double mud)
{
  double ret = cutoff_.mudim;

  cdptrs_.cptr = 1;
  cdptrs_.dptr = 1;
  cutoff_.mudim = mud;
  return ret;
}

double set_lambda(double lambda)
{
  double ret = ffcut_.delta;

  cdptrs_.cptr = 1;
  cdptrs_.dptr = 1;
  ffcut_.delta = lambda;
  return ret;
}

double set_delta(double del)
{
  double ret = cutoff_.divergence;

  cdptrs_.cptr = 1;
  cdptrs_.dptr = 1;
  cutoff_.divergence = del;
  return ret;
}

void flushcache(void)
{
  cdptrs_.cptr = 1;
  cdptrs_.dptr = 1;
  MLPutSymbol(stdlink, "Null");
  MLEndPacket(stdlink);
}

int main(int argc, char *argv[])
{
  int ret;

  savestdout = dup(1);
  dup2(2, 1);
  ffini_();
  dup2(savestdout, 1);
  ffflag_.ltest = 0;
  ffflag_.lwarn = 0;
  ffcut_.delta = LAMBDA;
  cutoff_.mudim = MUDIM;
  cutoff_.divergence = DELTA;
  cdptrs_.cptr = 1;
  cdptrs_.dptr = 1;
  ret = MLMain(argc, argv);
  dup2(2, 1);
  ffexi_();
  dup2(savestdout, 1);
  fflush(stderr);
  return ret;
}
