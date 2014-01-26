package stats

import (
	"math"
)

// Optimize finds the minimum point of a univariate-input function
// within a specified interval.
// This is a direct translation from the Fortran code at
//   www.netlib.no/notlib/fmm/fmin.f
// The method is due to Brent.
// The function optimize in R's stats package uses the same algorithm.
func Optimize(f func(float64) float64, ax, bx, tol float64) float64 {

	/*
	       an approximation  x  to the point where  f  attains a minimum  on
	   the interval  (ax,bx)  is determined.

	   INPUT..

	   ax    left endpoint of initial interval
	   bx    right endpoint of initial interval
	   f     function which evaluates  f(x, info)  for any  x
	         in the interval  (ax,bx)
	   tol   desired length of the interval of uncertainty of the final
	         result ( >= 0.)

	   OUTPUT..

	   fmin  abcissa approximating the point where  f  attains a minimum

	       The method used is a combination of  golden  section  search  and
	   successive parabolic interpolation.  convergence is never much slower
	   than  that  for  a  Fibonacci search.  If  f  has a continuous second
	   derivative which is positive at the minimum (which is not  at  ax  or
	   bx),  then  convergence  is  superlinear, and usually of the order of
	   about  1.324....
	       The function  f  is never evaluated at two points closer together
	   than  eps*abs(fmin)+(tol/3), where eps is  approximately  the  square
	   root  of  the  relative  machine  precision.   if   f   is a unimodal
	   function and the computed values of   f   are  always  unimodal  when
	   separated  by  at least  eps*abs(x)+(tol/3), then  fmin  approximates
	   the abcissa of the global minimum of  f  on the interval  ax,bx  with
	   an error less than  3*eps*abs(fmin)+tol.  if   f   is  not  unimodal,
	   then fmin may approximate a local, but perhaps non-global, minimum to
	   the same accuracy.
	       This function subprogram is a slightly modified  version  of  the
	   Algol  60 procedure  localmin  given in Richard Brent, Algorithms for
	   Minimization without Derivatives, Prentice-Hall, Inc. (1973).
	*/

	/*  c is the squared inverse of the golden ratio */
	c := (3. - math.Sqrt(5.)) * .5
	// TODO: write down this const value directly.

	/* Local variables */
	var a, b, d, e, p, q, r, u, v, w, x float64
	var t2, fu, fv, fw, fx, xm, eps, tol1, tol3 float64

	/*  eps is approximately the square root of the relative machine precision. */
	eps = 1.e-16
	tol1 = eps + 1. /* the smallest 1.000... > 1 */
	eps = 1.e-8

	a = ax
	b = bx
	v = a + c*(b-a)
	w = v
	x = v

	d = 0. /* -Wall */
	e = 0.
	fx = f(x)
	fv = fx
	fw = fx
	tol3 = tol / 3.

	for { // main loop
		xm = (a + b) * .5
		tol1 = eps*math.Abs(x) + tol3
		t2 = tol1 * 2.

		/* check stopping criterion */

		if math.Abs(x-xm) <= t2-(b-a)*.5 {
			break
		}
		p = 0.
		q = 0.
		r = 0.
		if math.Abs(e) > tol1 { /* fit parabola */

			r = (x - w) * (fx - fv)
			q = (x - v) * (fx - fw)
			p = (x-v)*q - (x-w)*r
			q = (q - r) * 2.
			if q > 0. {
				p = -p
			} else {
				q = -q
			}
			r = e
			e = d
		}

		if math.Abs(p) >= math.Abs(q*.5*r) ||
			p <= q*(a-x) ||
			p >= q*(b-x) { /* a golden-section step */

			if x < xm {
				e = b - x
			} else {
				e = a - x
			}
			d = c * e
		} else { /* a parabolic-interpolation step */

			d = p / q
			u = x + d

			/* f must not be evaluated too close to ax or bx */

			if u-a < t2 || b-u < t2 {
				d = tol1
				if x >= xm {
					d = -d
				}
			}
		}

		/* f must not be evaluated too close to x */

		if math.Abs(d) >= tol1 {
			u = x + d
		} else if d > 0. {
			u = x + tol1
		} else {
			u = x - tol1
		}

		fu = f(u)

		/*  update  a, b, v, w, and x */

		if fu <= fx {
			if u < x {
				b = x
			} else {
				a = x
			}
			v = w
			w = x
			x = u
			fv = fw
			fw = fx
			fx = fu
		} else {
			if u < x {
				a = u
			} else {
				b = u
			}
			if fu <= fw || w == x {
				v = w
				fv = fw
				w = u
				fw = fu
			} else if fu <= fv || v == x || v == w {
				v = u
				fv = fu
			}
		}
	} /* end of main loop */

	return x
}
