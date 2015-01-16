#ifndef __CTOOLS_H__
#define __CTOOLS_H__

#include <complex.h>	// this works for gcc 2.7.2
#include "ffc.h"

enum Ccount {
  cc0,
  cc1, cc2,
  cc00, cc11, cc12, cc22,
  cc001, cc002, cc111, cc112, cc122, cc222,
  clength
};

enum Dcount {
  dd0,
  dd1, dd2, dd3,
  dd00, dd11, dd12, dd13, dd22, dd23, dd33,
  dd001, dd002, dd003, dd111, dd112, dd113, dd122, dd123,
    dd133, dd222, dd223, dd233, dd333,
  dd0000, dd0011, dd0012, dd0013, dd0022, dd0023, dd0033,
    dd1111, dd1112, dd1113, dd1122, dd1123, dd1133,
    dd1222, dd1223, dd1233, dd1333, dd2222, dd2223,
    dd2233, dd2333, dd3333,
  dlength
};

typedef struct cpack {
  struct cpack *next;
  double P[6];
  double_complex C[clength];
} Cpack;

typedef struct dpack {
  struct dpack *next;
  double P[10];
  double_complex D[dlength];
} Dpack;

double_complex A0(double);
double_complex B0(double, double, double);
double_complex DB0(double, double, double);
double_complex B1(double, double, double);
double_complex DB1(double, double, double);
double_complex B00(double, double, double);
double_complex DB00(double, double, double);
double_complex B11(double, double, double);
double_complex DB11(double, double, double);
double_complex C0(double, double, double,
  double, double, double);
double_complex D0(double, double, double, double,
  double, double, double, double, double, double);
double_complex C0i(Ccount, double, double, double,
  double, double, double);
double_complex D0i(Dcount, double, double, double, double,
  double, double, double, double, double, double);

Cpack *Cget(double, double, double,
  double, double, double);
Dpack *Dget(double, double, double, double,
  double, double,
  double, double, double, double);

void resetheap(void);

double set_lambda(double);
double set_mudim(double);
double set_delta(double);
inline double get_lambda() {
  return ffcut_.delta;
}
inline double get_mudim() {
  extern double mudim_;
  return mudim_;
}
inline double get_delta() {
  extern double delta_;
  return delta_;
}

void bcaini();
#define bcaexi ffexi_ 

#endif

