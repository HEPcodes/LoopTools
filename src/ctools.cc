//	ctools.cc
//		this file is part of LoopTools;
//		it proves the scalar and tensor
// 		one-loop functions for C++
//		last modified 15 Jul 99 th


#define Quad doubledouble

#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "ctools.h"
#include "doubledouble.h"
#include "qcomplex.cc"

#define MAXDIGITS 9
#define MINMASS (1.e-5)
#define DELTA_DEFAULT 0.
#define MUDIM_DEFAULT 1.
#define LAMBDA_DEFAULT 1.

void Cijk(const double *, double_complex *);
void Dijkl(const double *, double_complex *);

double mudim_, delta_;
dcomplex a0i_[2], b0p_, b1p_, db0p_, pdb0p_;
Cpack *cstart, *cbackp = NULL, *clast = NULL;
Dpack *dstart, *dbackp = NULL, *dlast = NULL;
int cini = 0, dini = 0;


double_complex A0(double m)
{
  dcomplex result;
  int ier = 0;

  ffxa0_(&result, &delta_, &mudim_, &m, &ier);
  if(ier > MAXDIGITS) 
    printf("A0: lost %d digits (m = %lg)\n", ier, sqrt(m));
  return double_complex(result.r, result.i);
}

double_complex B0(double p, double m1, double m2)
{
  dcomplex result;
  int ier = 0;

  ffxb0_(&result, &delta_, &mudim_, &p, &m1, &m2, &ier);
  if(ier > MAXDIGITS) 
    printf("B0: lost %d digits (m = %lg)\n", ier, sqrt(m1));
  return double_complex(result.r, result.i);
}

double_complex DB0(double p, double m1, double m2)
{
  double dm;
  int ier = 0;

  if(p != 0. && m1 != 0. && m2 != 0. &&
     p == (dm = sqrt(m1) - sqrt(m2), dm *= dm)) {
    pdb0p_.r = (m2 - m1)/2./dm*log(m2/m1) - 2.;
    db0p_.r = pdb0p_.r/dm;
    db0p_.i = pdb0p_.i = 0.;
  }
  else {
    ffxdb0_(&db0p_, &pdb0p_, &p, &m1, &m2, &ier);
    if(ier > MAXDIGITS)
      printf("DB0: lost %d digits (m = %lg)\n", ier, sqrt(m1));
  }
  return double_complex(db0p_.r, db0p_.i);
}

double_complex C0(double p1, double p2, double p1p2,
  double m1, double m2, double m3)
{
  dcomplex result;
  double xpi[6] = {m1, m2, m3, p1, p2, p1p2};
  int ier = 0;

  ffxc0_(&result, xpi, &ier);
  if(ier > MAXDIGITS) {
    ier = 0;
    ffxc0r_(&result, xpi, &ier);
    if(ier > MAXDIGITS)
      printf("C0: lost %d digits (m = %lg)\n", ier, sqrt(m1)); 
  }
  return double_complex(result.r, result.i);
}

double_complex D0(double p1, double p2, double p3, double p4,
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
    xpi[9] = p1p2;		// numerically more stable
  }
  if(m1 + m2 + m3 + m4 <= MINMASS) ffxd0m0_(&result, xpi, &ier);
  else {
    ffxd0_(&result, xpi, &ier);
    if(ier > MAXDIGITS) {
      ier = 0;
      ffxd0r_(&result, xpi, &ier);
      if(ier > MAXDIGITS)
        printf("D0: lost %d digits (m=%lg)\n", ier, sqrt(m1));
    }
  }
  return double_complex(result.r, result.i);
}

double_complex B1(double p, double m1, double m2)
{
  int ier = 0;

  ffflag_.ldot = 1;
  ffxb0_(&b0p_, &delta_, &mudim_, &p, &m1, &m2, &ier);
  ffxa0_(&a0i_[0], &delta_, &mudim_, &m1, &ier);
  ffxa0_(&a0i_[1], &delta_, &mudim_, &m2, &ier);
  ffxb1_(&b1p_, &b0p_, a0i_,
    &p, &m1, &m2, ffdot_.fpij2, &ier);
  if(ier > MAXDIGITS)
    printf("B1: lost %d digits (m = %lg)\n", ier, sqrt(m1));
  ffflag_.ldot = 0;
  return double_complex(b1p_.r, b1p_.i);
}

double_complex DB1(double p, double m1, double m2)
{
  dcomplex result;
  int ier = 0;

  ffxdb1_(&result, &p, &m1, &m2, &ier);
  return double_complex(result.r, result.i);
}

double_complex B00(double p, double m1, double m2)
{
  dcomplex result[2];
  int ier = 0;

  B1(p, m1, m2);
  ffxb2p_(result, &b1p_, &b0p_, a0i_, &p, &m1, &m2, ffdot_.fpij2, &ier);
  if(ier > MAXDIGITS)
    printf("B00: lost %d digits (m = %lg)\n", ier, sqrt(m1));
  return double_complex(result[1].r, result[1].i);
}

double_complex DB00(double p, double m1, double m2)
{
  return 1./6.*(2.*m1*DB0(p, m1, m2) + B1(p, m1, m2) +
    (p + m1 - m2)*DB1(p, m1, m2)) - 1./18.;
}

double_complex B11(double p, double m1, double m2)
{
  dcomplex result[2];
  int ier = 0;

  B1(p, m1, m2);
  ffxb2p_(result, &b1p_, &b0p_, a0i_, &p, &m1, &m2, ffdot_.fpij2, &ier);
  if(ier > MAXDIGITS)
    printf("B11: lost %d digits (m = %lg)\n", ier, sqrt(m1));
  return double_complex(result[0].r, result[0].i);
}

double_complex DB11(double p, double m1, double m2)
{
  dcomplex result;
  int ier = 0;

  ffxdb11_(&result, &p, &m1, &m2, &ier);
  return double_complex(result.r, result.i);
}

void resetheap()
{
  if(clast) {
    clast->next = cbackp;
    clast = NULL;
  }
  if(dlast) {
    dlast->next = dbackp;
    dlast = NULL;
  }
  cbackp = cstart;
  dbackp = dstart;
  cini = dini = 0;
}

Cpack *Cget(double p1, double p2, double p1p2,
  double m1, double m2, double m3)
{
  Cpack *cp, *op;
  double P[] = {p1, p2, p1p2, m1, m2, m3};
  int i;

  if(cini) {
    for(cp = cstart; cp; cp = cp->next) {
      for(i = 5; cp->P[i] == P[i]; )
        if(!i--) return cp;
      op = cp;
    }
  }
  if(cbackp) cbackp = (clast = cp = cbackp)->next;
  else clast = NULL, cp = new Cpack;
  (cini ? op->next : cstart) = cp;
  cp->next = NULL;
  memcpy(cp->P, P, sizeof(P));
  Cijk(P, cp->C);
  cini = 1;
  return cp;
}

double_complex C0i(Ccount i, double p1, double p2, double p1p2,
  double m1, double m2, double m3)
{
  return Cget(p1, p2, p1p2, m1, m2, m3)->C[i];
}

Dpack *Dget(double p1, double p2, double p3, double p4,
  double p1p2, double p2p3,
  double m1, double m2, double m3, double m4)
{
  Dpack *dp, *op;
  double P[] = {p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4};
  int i;

  if(dini) {
    for(dp = dstart; dp; dp = dp->next) {
      for(i = 9; dp->P[i] == P[i]; )
        if(!i--) return dp;
      op = dp;
    }
  }
  if(dbackp) dbackp = (dlast = dp = dbackp)->next;
  else dlast = NULL, dp = new Dpack;
  (dini ? op->next : dstart) = dp;
  dp->next = NULL;
  memcpy(dp->P, P, sizeof(P));
  Dijkl(P, dp->D);
  dini = 1;
  return dp;
}

double_complex D0i(Dcount i,
  double p1, double p2, double p3, double p4,
  double p1p2, double p2p3,
  double m1, double m2, double m3, double m4)
{
  return Dget(p1, p2, p3, p4, p1p2, p2p3, m1, m2, m3, m4)->D[i];
}

void Cijk(const double P[], double_complex C[])
{
  Quad M12 = .5*(Quad(P[1]) - Quad(P[2]) - Quad(P[0]));
  Quad det2 = 2.*(P[2]*Quad(P[0]) - sqr(M12));
  M12 /= det2;
  Quad M11 = Quad(P[2])/det2;
  Quad M22 = Quad(P[0])/det2;

#define LC1() q2d(M11*s1 + M12*s2)
#define LC2() q2d(M12*s1 + M22*s2)

  Quad f1 = Quad(P[3]) - Quad(P[4]) + Quad(P[0]);
  Quad f2 = Quad(P[3]) - Quad(P[5]) + Quad(P[2]);

  C[cc0] = C0(P[0], P[1], P[2], P[3], P[4], P[5]);
  quad_complex b1123 = B11(P[1], P[5], P[4]);
  quad_complex b023 = quad_complex(b0p_.r, b0p_.i);
  quad_complex b123 = quad_complex(b1p_.r, b1p_.i);
  quad_complex b1113 = B11(P[2], P[3], P[5]);
  quad_complex b113 = quad_complex(b1p_.r, b1p_.i);
  quad_complex s1 = quad_complex(b0p_.r, b0p_.i) - b023 - f1*C[cc0];
  quad_complex b1112 = B11(P[0], P[3], P[4]);
  quad_complex b112 = quad_complex(b1p_.r, b1p_.i);
  quad_complex s2 = quad_complex(b0p_.r, b0p_.i) - b023 - f2*C[cc0];
  C[cc1] = LC1();
  C[cc2] = LC2();

  C[cc00] = q2d(.5*(Quad(P[3])*C[cc0] + .5*(f1*C[cc1]
    + f2*C[cc2] + b023)) + .25);
  s1 = -(f1*C[cc1] + b123) - 2.*C[cc00];
  s2 = -(f2*C[cc1] + b123 - b112);
  C[cc11] = LC1();
  C[cc12] = LC2();
  b023 += b123;
  s1 = b023 + b113 - f1*C[cc2];
  s2 = b023 - f2*C[cc2] - 2.*C[cc00];
  C[cc22] = LC2();

  C[cc001] = q2d((Quad(P[3])*C[cc1] + .5*(f1*C[cc11]
    + f2*C[cc12] + b123))*(1./3.) - 1./18.);
  C[cc002] = q2d((Quad(P[3])*C[cc2] + .5*(f1*C[cc12]
    + f2*C[cc22] - b023))*(1./3.) - 1./18.);
  s1 = -(b1123 + f1*C[cc11]) - 4.*C[cc001];
  s2 = -(b1123 + f2*C[cc11] - b1112);
  C[cc111] = LC1();
  C[cc112] = LC2();
  b1123 += b023 + b123;
  s1 = -(b1123 + f1*C[cc22] - b1113);
  s2 = -(b1123 + f2*C[cc22]) - 4.*C[cc002];
  C[cc122] = LC1();
  C[cc222] = LC2();
}

void Dijkl(const double P[], double_complex D[])
{
  double_complex *C234 = Cget(P[1], P[2], P[5], P[7], P[8], P[9])->C;
  double_complex *C134 = Cget(P[4], P[2], P[3], P[6], P[8], P[9])->C;
  double_complex *C124 = Cget(P[0], P[5], P[3], P[6], P[7], P[9])->C;
  double_complex *C123 = Cget(P[0], P[1], P[4], P[6], P[7], P[8])->C;
  quad_complex c1 = quad_complex(C234[cc0]) + quad_complex(C234[cc1])
    + quad_complex(C234[cc2]);
  quad_complex c2 = quad_complex(C234[cc12]) + quad_complex(C234[cc22])
    + quad_complex(C234[cc2]);
  quad_complex c3 = 2.*(C234[cc1] + c2) + quad_complex(C234[cc11])
    - quad_complex(C234[cc22]) + quad_complex(C234[cc0]);
  quad_complex c4 = quad_complex(C234[cc1]) + quad_complex(C234[cc11])
    + quad_complex(C234[cc12]);
  quad_complex c5 = quad_complex(C234[cc111]) + quad_complex(C234[cc112])
    + quad_complex(C234[cc11]);

  Quad f1 = .5*(Quad(P[0]) + Quad(P[4]) - Quad(P[1]));
  Quad f2 = .5*(Quad(P[0]) + Quad(P[3]) - Quad(P[5]));
  Quad f3 = .5*(Quad(P[4]) + Quad(P[3]) - Quad(P[2]));
  Quad M22 = P[3]*Quad(P[0]) - sqr(f2);
  Quad M23 = f1*f2 - P[0]*f3;
  Quad det3 = 2.*(P[4]*M22 - P[3]*sqr(f1) + f3*(f1*f2 + M23));
  Quad M11 = (P[4]*Quad(P[3]) - sqr(f3))/det3;
  Quad M12 = (f2*f3 - P[3]*f1)/det3;
  Quad M13 = (f1*f3 - P[4]*f2)/det3;
  M22 /= det3;
  M23 /= det3;
  Quad M33 = (P[4]*Quad(P[0]) - sqr(f1))/det3;

#define LD1() q2d(M11*s1 + M12*s2 + M13*s3)
#define LD2() q2d(M12*s1 + M22*s2 + M23*s3)
#define LD3() q2d(M13*s1 + M23*s2 + M33*s3)

  f1 = Quad(P[0]) + Quad(P[6]) - Quad(P[7]);
  f2 = Quad(P[4]) + Quad(P[6]) - Quad(P[8]);
  f3 = Quad(P[3]) + Quad(P[6]) - Quad(P[9]);

  D[dd0] = D0(P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], P[8], P[9]);
  quad_complex s1 = quad_complex(C134[cc0]) - C234[cc0] - f1*D[dd0];
  quad_complex s2 = quad_complex(C124[cc0]) - C234[cc0] - f2*D[dd0];
  quad_complex s3 = quad_complex(C123[cc0]) - C234[cc0] - f3*D[dd0];
  D[dd1] = LD1();
  D[dd2] = LD2();
  D[dd3] = LD3();

  D[dd00] = q2d(Quad(P[6])*D[dd0] + .5*(C234[cc0]
    + f1*D[dd1] + f2*D[dd2] + f3*D[dd3]));
  s1 = c1 - f1*D[dd1] - 2.*D[dd00];
  s2 = C124[cc1] + c1 - f2*D[dd1];
  s3 = C123[cc1] + c1 - f3*D[dd1];
  D[dd11] = LD1();
  s1 = quad_complex(C134[cc1]) - C234[cc1] - f1*D[dd2];
  s2 = -(C234[cc1] + f2*D[dd2]) - 2.*D[dd00];
  s3 = quad_complex(C123[cc2]) - C234[cc1] - f3*D[dd2];
  D[dd12] = LD1();
  D[dd22] = LD2();
  s1 = quad_complex(C134[cc2]) - C234[cc2] - f1*D[dd3];
  s2 = quad_complex(C124[cc2]) - C234[cc2] - f2*D[dd3];
  s3 = -(C234[cc2] + f3*D[dd3]) - 2.*D[dd00];
  D[dd13] = LD1();
  D[dd23] = LD2();
  D[dd33] = LD3();

  s1 = quad_complex(C134[cc00]) - C234[cc00] - f1*D[dd00];
  s2 = quad_complex(C124[cc00]) - C234[cc00] - f2*D[dd00];
  s3 = quad_complex(C123[cc00]) - C234[cc00] - f3*D[dd00];
  D[dd001] = LD1();
  D[dd002] = LD2();
  D[dd003] = LD3();
  s1 = -f1*D[dd11] - c3 - 4.*D[dd001];
  s2 = C124[cc11] - f2*D[dd11] - c3;
  s3 = C123[cc11] - f3*D[dd11] - c3;
  D[dd111] = LD1();
  D[dd112] = LD2();
  D[dd113] = LD3();
  s1 = quad_complex(C134[cc11]) - C234[cc11] - f1*D[dd22];
  s2 = -(C234[cc11] + f2*D[dd22]) - 4.*D[dd002];
  s3 = quad_complex(C123[cc22]) - C234[cc11] - f3*D[dd22];
  D[dd122] = LD1();
  D[dd222] = LD2();
  D[dd223] = LD3();
  s1 = quad_complex(C134[cc22]) - C234[cc22] - f1*D[dd33];
  s2 = quad_complex(C124[cc22]) - C234[cc22] - f2*D[dd33];
  s3 = -(C234[cc22] + f3*D[dd33]) - 4.*D[dd003];
  D[dd133] = LD1();
  D[dd233] = LD2();
  D[dd333] = LD3();
  s1 = c2 - f1*D[dd13] - 2.*D[dd003];
  s2 = C124[cc12] + c2 - f2*D[dd13];
  s3 = c2 - f3*D[dd13] - 2.*D[dd001];
  D[dd123] = LD2();

  D[dd0000] = q2d((Quad(P[6])*D[dd00] + .5*(C234[cc00] + f1*D[dd001]
    + f2*D[dd002] + f3*D[dd003]))*(1./3.) + 1./36.);
  D[dd0011] = q2d((Quad(P[6])*D[dd11] + .5*(c3 + f1*D[dd111]
    + f2*D[dd112] + f3*D[dd113]))*(1./3.));
  D[dd0012] = q2d((Quad(P[6])*D[dd12] - .5*(c4 - f1*D[dd112]
    - f2*D[dd122] - f3*D[dd123]))*(1./3.));
  D[dd0013] = q2d((Quad(P[6])*D[dd13] - .5*(c2 - f1*D[dd113]
    - f2*D[dd123] - f3*D[dd133]))*(1./3.));
  D[dd0022] = q2d((Quad(P[6])*D[dd22] + .5*(C234[cc11] + f1*D[dd122]
    + f2*D[dd222] + f3*D[dd223]))*(1./3.));
  D[dd0023] = q2d((Quad(P[6])*D[dd23] + .5*(C234[cc12] + f1*D[dd123]
    + f2*D[dd223] + f3*D[dd233]))*(1./3.));
  D[dd0033] = q2d((Quad(P[6])*D[dd33] + .5*(C234[cc22] + f1*D[dd133]
    + f2*D[dd233] + f3*D[dd333]))*(1./3.));

  c3 = c2 + quad_complex(C234[cc112]) + quad_complex(C234[cc122]);
  c2 += c3 - quad_complex(C234[cc2]) + quad_complex(C234[cc122])
    + quad_complex(C234[cc222]);
  c4 += c1 + c2 + c3 + c4 + c5 + C234[cc12];

  s1 = c4 - f1*D[dd111] - 6.*D[dd0011];
  s2 = C124[cc111] + c4 - f2*D[dd111];
  s3 = C123[cc111] + c4 - f3*D[dd111];
  D[dd1111] = LD1();
  D[dd1112] = LD2();
  D[dd1113] = LD3();
  s1 = -(c2 + f1*D[dd113]) - 4.*D[dd0013];
  s2 = C124[cc112] - c2 - f2*D[dd113];
  s3 = -(c2 + f3*D[dd113]) - 2.*D[dd0011];
  D[dd1123] = LD2();
  D[dd1133] = LD3();
  s1 = c5 - f1*D[dd122] - 2.*D[dd0022];
  s2 = c5 - f2*D[dd122] - 4.*D[dd0012];
  s3 = C123[cc122] + c5 - f3*D[dd122];
  D[dd1122] = LD1();
  D[dd1223] = LD3();
  s1 = quad_complex(C134[cc111]) - C234[cc111] - f1*D[dd222];
  s2 = -(C234[cc111] + f2*D[dd222]) - 6.*D[dd0022];
  s3 = quad_complex(C123[cc222]) - C234[cc111] - f3*D[dd222];
  D[dd1222] = LD1();
  D[dd2222] = LD2();
  D[dd2223] = LD3();
  s1 = quad_complex(C134[cc122]) - C234[cc122] - f1*D[dd233];
  s2 = -(C234[cc122] + f2*D[dd233]) - 2.*D[dd0033];
  s3 = -(C234[cc122] + f3*D[dd233]) - 4.*D[dd0023];
  D[dd1233] = LD1();
  D[dd2233] = LD2();
  s1 = quad_complex(C134[cc222]) - C234[cc222] - f1*D[dd333];
  s2 = quad_complex(C124[cc222]) - C234[cc222] - f2*D[dd333];
  s3 = -(C234[cc222] + f3*D[dd333]) - 6.*D[dd0033];
  D[dd1333] = LD1();
  D[dd2333] = LD2();
  D[dd3333] = LD3();
}

double set_lambda(double lambda)
{
  double ret = ffcut_.delta;
  ffcut_.delta = lambda;
  resetheap();
  return ret;
}

double set_mudim(double mud)
{
  double ret = mudim_;
  mudim_ = mud;
  resetheap();
  return ret;
}

double set_delta(double del)
{
  double ret = delta_;
  delta_ = del;
  resetheap();
  return ret;
}

void bcaini()
{
  ffini_();
  ffflag_.ltest = 0;
  ffflag_.lwarn = 0;
  ffcut_.delta = LAMBDA_DEFAULT;
  mudim_ = MUDIM_DEFAULT;
  delta_ = DELTA_DEFAULT;
}

