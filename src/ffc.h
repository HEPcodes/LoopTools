#ifndef __FFC_H__
#define __FFC_H__

/* this line is needed for configure */

/* some prototypes & common blocks from the FF package */

#ifdef NO_UNDERSCORE
#define ffini_ ffini
#define ffexi_ ffexi
#define ffxa0_ ffxa0
#define ffxb0_ ffxb0
#define ffxdb0_ ffxdb0
#define ffxb1_ ffxb1
#define ffxdb1_ ffxdb1
#define ffxdb11_ ffxdb11
#define ffxb2p_ ffxb2p
#define ffxc0_ ffxc0
#define ffxc0r_ ffxc0r
#define ffxd0_ ffxd0
#define ffxd0r_ ffxd0r
#define ffxd0m0_ ffxd0m0
#define ffflag_ ffflag
#define ffcut_ ffcut
#define ffdot_ ffdot
#endif

typedef struct {
  double r, i;
} dcomplex;

#ifdef __cplusplus
extern "C" {
#endif

int ffini_();

int ffexi_();

int ffxa0_(dcomplex *cb0,
  const double *d0, const double *xmu,
  const double *xm, int *ier);

int ffxb0_(dcomplex *cb0,
  const double *d0, const double *xmu,
  const double *xk, const double *xma,
  const double *xmb, int *ier);

int ffxdb0_(dcomplex *cdb0, dcomplex *cdb0p,
  const double *xp, const double *xma,
  const double *xmb, int *ier);

int ffxb1_(dcomplex *cb1,
  const dcomplex *cb0, const dcomplex *ca0i,
  const double *xp, const double *xma,
  const double *xmb, const double *piDpj, int *ier);

int ffxdb1_(dcomplex *cdb1,
  const double *p1, const double *m1,
  const double *m2, int *ier);

int ffxdb11_(dcomplex *cdb11,
  const double *p1, const double *m1,
  const double *m2, int *ier);

int ffxb2p_(dcomplex *cb2, const dcomplex *cb1,
  const dcomplex *cb0, const dcomplex *ca0i,
  const double *xp, const double *xma,
  const double *xmb, const double *piDpj, int *ier);

int ffxc0_(dcomplex *cc0, const double *xpi, int *ier);

int ffxc0r_(dcomplex *cc0, const double *xpi, int *ier);

int ffxd0_(dcomplex *cd0, const double *xpi, int *ier);

int ffxd0r_(dcomplex *cd0, const double *xpi, int *ier);

int ffxd0m0_(dcomplex *cd0, const double *xpi, int *ier);

extern struct {
  int lwrite, ltest, l4also, ldc3c4, lmem, lwarn, ldot;
  int nevent, ner, id, idsub, nwidth, nschem;
  int onshel;
  int idot;
} ffflag_;

extern struct {
  double delta;
} ffcut_;

extern struct {
  double fpij2[9], fpij3[36], fpij4[100],   
    fpij5[225], fpij6[441];
} ffdot_;

#ifdef __cplusplus
}
#endif

#endif

