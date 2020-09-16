/* Original: stl2.f -- translated by f2c (version 19970219). */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void linear_approx(double *x, double *y, double *e, int *m, double *u, 
		   double *v, double *w, int *k, int *ip) {

  /* System generated locals */
  int i__1;
  double r__1, r__2;

  /* Local variables */
  int indc, keep, itch, init;
  double smin;
  int ipiv;
  double smax, xeye, yeye, test, svmn=0, svmx=0, temp1, temp2;
  int i__, j, l, n, idiot;
  double epsln, slope=0, xinit, yinit;
  int kp=0;
  double dx;
  int it, igraze;
  double sgn, svx=0, svy=0;

  /* From D.G. Wilson, 1976, ACM Transactions on Mathematical Software,2,388 */

  /* Piecewise linear approximations of fewest line segments within
     given tolerances.  x, y, e and m contain input data.  u,v,k and
     possibly w contain output.  ip is a parameter determining the
     operation of the program.  x and y are input data arrays of m
     elements x(i),y(i) contains the ith data point.  e may be a
     single tolerance or a table of tolerances depending on the value
     of ip.  If e is an array, then e(i) is the tolerance associated
     with x(i),y(i) and e must contain m nonnegative elements.  u and
     v are output arrays of k+1 elements.  u is a partition of the
     interval (x(1),x(n)) with u(1)=x(1) and u(k+1)=x(n).  v(i) is an
     ordinate to be associated with u(i) in the approximation.  (If a
     continuous approximation is requested, then v(i) is 'the'
     ordinate to be associated with u(i).)  If a continuous
     approximation is requested, then w is not used.  In this case the
     ith approximating segment is the straight line from u(i),v(i) to
     u(i+1),v(i+1).  If a continuous approximation is not requested,
     then w is a k-element output array.  In this case the ith
     approximating segment is the straight line from u(i),w(i) to
     u(i+1),v(i+1), and v(1) is set equal to w(1).  k is the number of
     segments in the piece- wise linear approximation generated.  In
     case of an error return, k will be set to zero.

     The control parameter ip is the product of three indicators i1,i2
     and i3.  i1 indicates whether or not e is an array of tolerances.
     i1 = -1 indicates e is an array i1 = +1 indicates e is a single
     number.

     i2 indicates whether or not the approximation is to be restricted to
     the 'tolerance band' about the data.  
     i2 = 1 indicates no band restriction 
     i2 = 2 indicates apply this restriction 
     (the 'tolerance band' is a piecewise linear band centered at the
     data whose width is determined by the tolerances at the data
     points.)

     i3 indicates whether or not the approximation must be continuous. 
     i3 = 1 indicates continuity not required 
     i3 = 3 indicates continuity is required 

     call stl2 (x,y,e,m,x,y,e,m,ip) will not cause problems provided
     that either a continuous approximation is equested, or e is a
     sufficiently large array.

     The program performs the following data checks.  Are the x-values
     in increasing order?  Are the tolerance(s) nonnegative?  Is the
     number of data points greater than one?  If any check fails, the
     program returns with k set equal to 0.  In this case no further
     processing is attempted.  */

    /* Parameter adjustments */
    --w;
    --v;
    --u;
    --e;
    --y;
    --x;

    /* Function Body */
    n = *m;
    itch = *ip;
    j = 1;
/* ERROR CHECKS */
    if (n <= 1) {
      printf("linear_approx: n (%d) is <= 1\n",n);
      goto L400;
    }
    if (e[1] < (double)0.) {
      printf("linear_approx: tolerance[0] (%e) is < 0\n",e[1]);
      goto L400;
    }
    i__1 = n;
    for (l = 2; l <= i__1; ++l) {
	if (x[l - 1] >= x[l]) {
	  printf("linear_approx: x[%d] (%e) > x[%d] (%e)\n",l-2,x[l-1],
		   l-1,x[l]);
	  goto L400;
	}
	if (itch >= 0) {
	    goto L10;
	}
	if (e[l] < (double)0.) {
	  printf("linear_approx: tolerance[%d] (%e) is < 0 for (%e)\n",
		 l-1,e[l],y[l]);
	    goto L400;
	}
L10:
	;
    }
/* INITIALIZATION FOR ENTIRE PROGRAM */
    epsln = e[1];
    sgn = (double)1.;
    keep = 1;
    i__ = 1;
    u[1] = x[1];
    j = 2;
    init = 1;
    indc = 0;
    goto L30;
/* INITIALIZATION FOR EACH SEGMENT */
L20:
    ++j;
    init = i__;
    indc = 0;
    if (abs(itch) < 3) {
	keep = i__;
    }
    if ((i__1 = abs(itch) - 4, abs(i__1)) != 2) {
	goto L30;
    }
/* RESTRICTED TO TOLERANCE BAND */
    xeye = u[j - 1];
    yeye = v[j - 1];
    temp1 = epsln;
    if (itch < 0) {
	temp1 += (sgn * e[i__ - 1] - epsln) * (x[i__] - u[j - 1]) / (x[i__] - 
		x[i__ - 1]);
    }
    yinit = yeye - temp1 - temp1;
    goto L40;
L30:
/* NOT RESTRICTED TO TOLERANCE BAND */
    xeye = x[i__];
    yeye = y[i__] + epsln;
    yinit = y[i__] - epsln;
    if (abs(itch) == 1 || i__ == 1) {
	goto L40;
    }
    temp1 = epsln;
    if (itch < 0) {
	temp1 = sgn * e[i__ + 1];
    }
    smin = (y[i__ + 1] - yeye - temp1) / (x[i__ + 1] - xeye);
    if (itch < 0) {
	temp1 = sgn * e[i__ - 1];
    }
    smax = (yeye - y[i__ - 1] + temp1) / (xeye - x[i__ - 1]);
    if (keep == i__ - 1) {
	goto L50;
    }
    it = i__ - 2;
    xinit = xeye;
    ipiv = i__;
    igraze = i__;
    ++i__;
    goto L150;
L40:
    if (xeye >= x[i__]) {
	++i__;
    }
    if (itch < 0) {
	epsln = sgn * e[i__];
    }
    dx = x[i__] - xeye;
    smax = (y[i__] + epsln - yeye) / dx;
    smin = (y[i__] - epsln - yeye) / dx;
L50:
    xinit = xeye;
    ipiv = i__;
    igraze = i__;
/* DETERMINATION OF INDIVIDUAL SEGMENT */
L60:
    if (i__ == n) {
	goto L260;
    }
    ++i__;
L70:
/* TEST FOR NEW *MAX* SLOPE */
    dx = x[i__] - xeye;
    if (itch < 0) {
	epsln = sgn * e[i__];
    }
    temp1 = (y[i__] + epsln - yeye) / dx;
    test = temp1 - smax;
    if (sgn <= (double)0.) {
	test = -test;
    }
    if (test < (double)0.) {
	goto L80;
    } else if (test == 0) {
	goto L90;
    } else {
	goto L100;
    }
L80:
/* TEST FOR END OF CANDIDATE SEGMENT */
    test = temp1 - smin;
    if (sgn <= (double)0.) {
	test = -test;
    }
    if (test < (double)0.) {
	goto L210;
    }
    smax = temp1;
L90:
/* TEST FOR NEW *MIN* SLOPE */
    ipiv = i__;
L100:
    temp2 = (y[i__] - epsln - yeye) / dx;
    test = temp2 - smax;
    if (sgn <= (double)0.) {
	test = -test;
    }
    if (test < (double)0.) {
	goto L110;
    } else if (test == 0) {
	goto L120;
    } else {
	goto L140;
    }
L110:
    test = smin - temp2;
    if (sgn <= (double)0.) {
	test = -test;
    }
    if (test < (double)0.) {
	goto L120;
    } else if (test == 0) {
	goto L130;
    } else {
	goto L60;
    }
L120:
    smin = temp2;
L130:
    igraze = i__;
    goto L60;
/* CHECK FOR PIVOT AT NEW EYE POINT */
L140:
    if (xeye == x[ipiv]) {
	goto L220;
    }
    if (itch < 0) {
	epsln = sgn * e[ipiv];
    }
    indc = 1;
    svx = xeye;
    svy = yeye;
    svmn = smin;
    svmx = smax;
    xeye = x[ipiv];
    yeye = y[ipiv] + epsln;
    smin = smax;
    smax = (yinit - yeye) / (xinit - xeye);
    if (keep >= ipiv) {
	goto L170;
    }
    it = ipiv - 1;
L150:
    temp2 = yeye + epsln;
    i__1 = it;
    for (l = keep; l <= i__1; ++l) {
	if (itch < 0) {
	    temp2 = yeye + sgn * e[l];
	}
	temp1 = (y[l] - temp2) / (x[l] - xeye);
	test = temp1 - smax;
	if (sgn <= (double)0.) {
	    test = -test;
	}
	if (test < (double)0.) {
	    smax = temp1;
	}
/* L160: */
    }
L170:
    if (ipiv >= i__ - 1) {
	goto L70;
    }
    it = i__ - 2;
    temp2 = yeye - epsln;
    idiot = ipiv;
    i__1 = it;
    for (l = idiot; l <= i__1; ++l) {
	dx = x[l + 1] - xeye;
	if (itch < 0) {
	    temp2 = yeye - sgn * e[l + 1];
	}
	temp1 = (y[l + 1] - temp2) / dx;
	test = temp1 - smax;
	if (sgn <= (double)0.) {
	    test = -test;
	}
	if (test < (double)0.) {
	    goto L180;
	} else if (test == 0) {
	    goto L190;
	} else {
	    goto L200;
	}
L180:
	smax = temp1;
L190:
	ipiv = l + 1;
L200:
	;
    }
    goto L70;
/* END OF CURRENT SEGMENT */
L210:
    temp2 = smin;
    if (i__ == n) {
	goto L240;
    }
    keep = igraze;
    goto L250;
L220:
    temp2 = smax;
    if (i__ == n) {
	goto L230;
    }
    sgn = -sgn;
    epsln = -epsln;
    keep = ipiv;
    goto L250;
L230:
    if (indc == 0 || xeye != x[n - 1]) {
	goto L240;
    }
    xeye = svx;
    yeye = svy;
    smin = svmn;
    smax = svmx;
L240:
    u[j] = x[n - 1];
    yinit = y[n - 1];
    goto L270;
L250:
    if ((i__1 = abs(itch) - 4, abs(i__1)) != 2) {
	goto L300;
    }
/* DETERMINE KNOT ON EDGE OF TOLERANCE BAND */
    temp1 = (double)0.;
    if (itch < 0) {
	temp1 = epsln - sgn * e[i__ - 1];
    }
    temp1 = (y[i__] - y[i__ - 1] + temp1) / (x[i__] - x[i__ - 1]);
    u[j] = (y[i__] + epsln - yeye - temp1 * x[i__] + temp2 * xeye) / (temp2 - 
	    temp1);
    goto L310;
L260:
    u[j] = x[n];
    yinit = y[n];
L270:
/* CONTINUITY CHECK FOR LAST SEGMENT */
    if (abs(itch) >= 3 || init == 1) {
	goto L290;
    }
    it = init - 1;
    svmx = smax + sgn;
    temp2 = yeye + epsln;
    i__1 = it;
    for (l = kp; l <= i__1; ++l) {
	if (itch < 0) {
	    temp2 = yeye + sgn * e[l];
	}
	temp1 = (y[l] - temp2) / (x[l] - xeye);
	test = temp1 - svmx;
	if (sgn <= (double)0.) {
	    test = -test;
	}
	if (test < (double)0.) {
	    svmx = temp1;
	}
/* L280: */
    }
    if ((r__1 = svmx - smax + svmx - smin, fabs(r__1)) <= (r__2 = smax - smin,
	     fabs(r__2))) {
	smax = svmx;
    }
L290:
/* NEARNESS CHECK FOR LAST SEGMENT */
    temp2 = smax;
    temp1 = yeye + smax * (u[j] - xeye);
    test = yinit - temp1;
    if (sgn < (double)0.) {
	test = -test;
    }
    if (test > (double)0.) {
	goto L310;
    }
    temp2 = smin;
    temp1 = yeye + smin * (u[j] - xeye);
    test = yinit - temp1;
    if (sgn < (double)0.) {
	test = -test;
    }
    if (test < (double)0.) {
	goto L310;
    }
    temp2 = (yinit - yeye) / (u[j] - xeye);
    v[j] = yinit;
    goto L320;
L300:
    if (abs(itch) >= 3) {
	goto L330;
    }
    u[j] = (x[i__] + x[i__ - 1]) * (double).5;
L310:
    v[j] = yeye + temp2 * (u[j] - xeye);
L320:
    if (xeye != xinit) {
	goto L330;
    }
    if (abs(itch) == 2) {
	goto L360;
    }
    if (abs(itch) != 6) {
	goto L330;
    }
    if (j <= 2) {
	goto L380;
    }
    goto L390;
L330:
/* RECOMPUTATION OF KNOT FOR CONTINUITY */
    if (j <= 2) {
	goto L370;
    }
    if (slope == temp2) {
	goto L360;
    }
    yinit = v[j - 2];
    if (abs(itch) < 3) {
	yinit = w[j - 2];
    }
    temp1 = (xeye * temp2 - u[j - 2] * slope + yinit - yeye) / (temp2 - slope)
	    ;
    if (abs(itch) >= 3) {
	goto L350;
    }
    if (temp1 > xinit) {
	goto L360;
    }
    test = fabs(epsln);
    idiot = init - kp;
    i__1 = idiot;
    for (l = 1; l <= i__1; ++l) {
	it = init - l;
	if (temp1 >= x[it]) {
	    goto L350;
	}
	dx = y[it] - yeye - temp2 * (x[it] - xeye);
	if (itch < 0) {
	    test = e[it];
	}
	if (fabs(dx) > test) {
	    goto L360;
	}
/* L340: */
    }
L350:
    u[j - 1] = temp1;
    v[j - 1] = yeye + temp2 * (u[j - 1] - xeye);
    if (abs(itch) < 3) {
	w[j - 1] = v[j - 1];
    }
    goto L390;
L360:
    w[j - 1] = yeye + temp2 * (u[j - 1] - xeye);
    goto L390;
L370:
    if (abs(itch) < 3) {
	goto L360;
    }
L380:
    v[1] = yeye + temp2 * (u[1] - xeye);
L390:
    slope = temp2;
    kp = keep;
    if (i__ < n) {
	goto L20;
    }
    if (x[n] == u[j]) {
	goto L400;
    }
    if (abs(itch) < 3) {
	w[j] = v[j];
    }
    ++j;
    u[j] = x[n];
    v[j] = y[n];
L400:
    if (j >= 2 && abs(itch) < 3) {
	v[1] = w[1];
    }
    *k = j - 1;
}

