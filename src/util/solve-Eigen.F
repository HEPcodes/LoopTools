* solve-Eigen.F
* computation of the inverse and solution of a linear system
* by diagonalizing the matrix with the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of LoopTools
* last modified 29 Jan 08 th

#include "defs.h"


* A matrix is considered diagonal if the sum of the squares
* of the off-diagonal elements is less than EPS.

#define EPS 2D0**(-102)


************************************************************************
** XEigen diagonalizes a general complex n-by-n matrix.
** Input: n, A = n-by-n matrix
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^-1.

	subroutine XEigen(n, Ain,ldAin, d, U,ldU)
	implicit none
	integer n, ldAin, ldU
	QVAR Ain(ldAin,*)
	double complex U(ldU,*), d(*)

	integer p, q, j, sweep
	QREAL red, off, thresh
	double complex delta, t, invc, sx, sy, tx, ty
	double complex x, y
	double complex ev(2,MAXDIM), A(MAXDIM,MAXDIM)

	QREAL sq
	QVAR c
	sq(c) = QPREC(c*QCC(c))

	do p = 1, n
	  do q = 1, n
	    A(q,p) = Ain(q,p)
	  enddo
	enddo

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = A(p,p)
	  d(p) = ev(2,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .01D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + sq(A(p,q)) + sq(A(q,p))
	    enddo
	  enddo
	  if( off .lt. EPS ) return

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      off = sq(A(p,q)) + sq(A(q,p))
	      if( sweep .gt. 4 .and. off .lt.
     &              EPS*max(sq(ev(2,p)), sq(ev(2,q))) ) then
	        A(p,q) = 0
	        A(q,p) = 0
	      else
	        if( off .gt. thresh ) then
	          delta = A(p,q)*A(q,p)
	          x = .5D0*(ev(2,p) - ev(2,q))
	          y = sqrt(x**2 + delta)
	          t = x - y
	          x = x + y
	          if( abs(t) .lt. abs(x) ) t = x

	          t = 1/t
	          delta = delta*t
	          ev(1,p) = ev(1,p) + delta
	          ev(2,p) = d(p) + ev(1,p)
	          ev(1,q) = ev(1,q) - delta
	          ev(2,q) = d(q) + ev(1,q)

	          invc = sqrt(delta*t + 1)
	          x = t/invc
	          t = t/(invc + 1)
	          sx = x*A(p,q)
	          ty = t*A(p,q)
	          sy = x*A(q,p)
	          tx = t*A(q,p)

	          do j = 1, n
	            x = A(j,p)
	            y = A(j,q)
	            A(j,p) = x + sy*(y - ty*x)
	            A(j,q) = y - sx*(x + tx*y)
	            x = A(p,j)
	            y = A(q,j)
	            A(p,j) = x + sx*(y - tx*x)
	            A(q,j) = y - sy*(x + ty*y)
	          enddo

	          A(p,q) = 0
	          A(q,p) = 0

	          do j = 1, n
	            x = U(p,j)
	            y = U(q,j)
	            U(p,j) = x + sx*(y - tx*x)
	            U(q,j) = y - sy*(x + ty*y)
	          enddo
	        endif
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in XEigen"
	end

************************************************************************

	subroutine XDet(n, A,ldA, det)
	implicit none
	integer n, ldA
	QVAR A(ldA,*), det

	double complex d(MAXDIM), U(MAXDIM,MAXDIM), cdet
	integer p

	call XEigen(n, A,ldA, d, U,MAXDIM)

	cdet = 1
	do p = 1, n
	  cdet = cdet*d(p)
	enddo
	det = cdet
	end

************************************************************************
** XInverse forms the (pseudo)inverse of a symmetric n-by-n matrix.
** Input: n, A = n-by-n matrix, symmetric
** (only the upper triangle of A needs to be filled).
** Output: Ainv = (pseudo)inverse of A

	subroutine XInverse(n, A,ldA, Ainv,ldAinv)
	implicit none
	integer n, ldA, ldAinv
	QVAR A(ldA,*), Ainv(ldAinv,*)

	integer p, q, j
	double complex U(MAXDIM,MAXDIM), d(MAXDIM)
	QVAR t

	call XEigen(n, A,ldA, d, U,MAXDIM)

* form (pseudo)inverse U^+ ev^-1 U

	do p = 1, n
	  do q = 1, n
	    t = 0
	    do j = 1, n
	      if( abs(d(j)) .gt. EPS )
     &          t = t + DCONJG(U(j,q))*U(j,p)/d(j)
	    enddo
	    Ainv(q,p) = t
	  enddo
	enddo
	end

************************************************************************

	subroutine XSolve(n, A,ldA, Ainv,ldAinv, b)
	implicit none
	integer n, ldA, ldAinv
	QVAR A(ldA,*), Ainv(ldAinv,*)
	double complex b(*)

	integer i, j
	double complex x(MAXDIM)

#if defined(QUAD) && !defined(COMPLEXPARA)

	QVAR tr, ti

    	do i = 1, n
	  tr = 0
	  ti = 0
	  do j = 1, n
	    tr = tr + Ainv(i,j)*DBLE(b(j))
	    ti = ti + Ainv(i,j)*DIMAG(b(j))
	  enddo
	  x(i) = DCMPLX(tr, ti)
	enddo

#else

	double complex delta(MAXDIM), t

	do i = 1, n
	  t = 0
	  do j = 1, n
	    t = t + Ainv(i,j)*b(j)
	  enddo
	  x(i) = t
	enddo

* improve on x
	do i = 1, n
	  t = 0
	  do j = 1, n
	    t = t + A(i,j)*x(j)
	  enddo
	  delta(i) = b(i) - t
	enddo
	do i = 1, n
	  t = 0
	  do j = 1, n
	    t = t + Ainv(i,j)*delta(j)
	  enddo
	  x(i) = x(i) + t
	enddo

#endif

	do i = 1, n
	  b(i) = x(i)
	enddo
	end

