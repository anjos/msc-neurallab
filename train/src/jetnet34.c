/* jetnet34.f -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "jetnet.h"

/* Initialized data */

struct {
    integer e_1[40];
    real e_2[40];
    integer e_3[20];
    real e_4[20];
    integer fill_5[2001];
    } jndat1_ = { {3, 10, 1, 0, 0, 6, 0, 0, 100, 16, 8, 1, 0, 0, 0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 10, 0, 0, 10, 10, 0, 0, 0, 0},
	     {.001f, .5f, 1.f, .1f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.f, 1.f, 
	    1.f, 0.f, 1e-6f, .9f, .9f, 1.f, 0.f, 1.f, 1.75f, 1e3f, 0.f, .1f, 
	    .05f, .001f, 2.f, 1e-4f, 1e-6f, 1.2f, .5f, 50.f, 1e-6f, 0.f, 0.f, 
	    0.f, 0.f, 0.f, 0.f, 0.f}, {1, 0, 2, 1, 0, 6, 0, 0, 0, 10, 10, 1, 0,
	    0, 0, 0, 0, 0, 0, 0}, {.001f, 0.f, .01f, .5f, 0.f, 0.f, 0.f, 0.f, 
	    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f} };

struct {
    real e_1[10];
    integer e_2[10];
    real e_3[30];
    } jndat2_ = { {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}, {0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0}, {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
	    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 
	    0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f} };

struct {
    integer e_1[2];
    real e_2;
    integer e_3;
    real e_4[11];
    integer e_5[4];
    real e_6;
    } jnint4_ = { {0, 0}, 0.f, 0, {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
		  0.f, 0.f}, {0, 1, 0, 0}, 0.f }; 

struct {
    integer fill_1[45];
    real e_2[12];
    integer e_3;
    } jnint2_ = { {0}, {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 
	    0.f}, 0 };

struct {
    integer fill_1[15];
    integer e_2;
    } jmint1_ = { {0}, 10 };

struct {
    integer e_1;
    real e_2;
    } jngaus_ = { 0, 0.f };

struct {
    integer e_1[5];
    integer fill_2[100];
    } jndatr_ = { {19780503, 0, 0, 97, 33} };


/* Table of constant values */

static integer c__0 = 0;
static integer c__2 = 2;
static real c_b15 = 1.f;
static integer c__15 = 15;
static integer c__14 = 14;
static integer c__30 = 30;
static integer c__1 = 1;
static integer c__40 = 40;
static integer c__10 = 10;
static integer c__2000 = 2000;
static integer c_b168 = 150000;
static integer c__1000 = 1000;
static integer c__28 = 28;
static integer c__29 = 29;
static integer c__26 = 26;
static integer c__24 = 24;
static integer c__27 = 27;
static integer c__5 = 5;
static integer c__12 = 12;
static integer c__16 = 16;
static integer c__21 = 21;
static integer c__17 = 17;
static integer c__13 = 13;
static doublereal c_b469 = 10.;
static integer c__20 = 20;
static integer c__18 = 18;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;
static integer c__25 = 25;
static integer c__11 = 11;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__31 = 31;
static integer c__32 = 32;
static integer c__19 = 19;
static integer c__23 = 23;
static integer c__22 = 22;
static real c_b852 = .5f;
static integer c__9 = 9;
static integer c_n30 = -30;
static integer c__121 = 121;
static real c_b1164 = 2.f;
static real c_b1168 = 0.f;
static real c_b1194 = 85.2f;

/* **********************************************************************C */
/*                                                                      C */
/*                          J E T N E T - 3.4                           C */
/*                                                                      C */
/*           A Neural Network program for jet discrimination            C */
/*         and other High Energy Physics triggering situations          C */
/*                                                                      C */
/*                    Latest date of change 95.07.06                    C */
/*                                                                      C */
/*                              Authors :                               C */
/*                                                                      C */
/*                   Leif Lonnblad, Carsten Peterson,                   C */
/*                 Hong Pi and Thorsteinn Rognvaldsson                  C */
/*                                                                      C */
/*                   Department of Theoretical Physics                  C */
/*                  University of Lund, Solvegatan 14A,                 C */
/*                           S-223 62 Lund                              C */
/*                               Sweden                                 C */
/*                                                                      C */
/*                        tel  int+46-46109073                          C */
/*                        fax  int+46-46104438                          C */
/*                                                                      C */
/*                     BITNET/EARN THEPLL@SELDC52                       C */
/*                                 THEPCAP@SELDC52                      C */
/*                                 THEPHP@SELDC52                       C */
/*                                 THEPDR@SELDC52                       C */
/*                                                                      C */
/*                     internet    leif@thep.lu.se                      C */
/*                                 carsten@thep.lu.se                   C */
/*                                 pihong@thep.lu.se                    C */
/*                                 denni@thep.lu.se                     C */
/*                                                                      C */
/* Copyright 1991,1992,1993,1994,1995 L. Lonnblad & Th. Rognvaldsson    C */
/*                                                                      C */
/*          Please report any errors to: <denni@thep.lu.se>             C */
/*                                                                      C */
/* **********************************************************************C */
/* **********************************************************************C */
/*                                                                      C */
/* An updated version of the program is obtainable through anonymous    C */
/* ftp from thep.lu.se in directory /pub/LundPrograms/Jetnet/           C */
/*                                                                      C */
/* **********************************************************************C */
/* **********************************************************************C */
/*  A description of the models and the program can be found in:        C */
/*                                                                      C */
/*  (i) Lonnblad et. al., "Self-organizing Networks for Extracting      C */
/*   Jet Features", Computer Physics Communications, vol. 67,           C */
/*   pp. 193-209, 1991.                                                 C */
/*                                                                      C */
/*  (ii) Lonnblad et. al., "Pattern recognition in High Energy Physics  C */
/*  with Artificial Neural Networks - JETNET 2.0", Computer Physics     C */
/*  Communications, nr. 70, pp. 167-182, 1992.                          C */
/*                                                                      C */
/*  (iii) Lonnblad et. al. "JETNET 3.0 - A Versatile Artificial         C */
/*  Neural Network Package", Computer Physics Communications, vol. 81,  C */
/*  pp. 185-220, 1994.                                                  C */
/*                                                                      C */
/* **********************************************************************C */
/* **********************************************************************C */
/*         Order of appearance of subroutines and functions:            C */
/*                                                                      C */
/* 1) Feed-forward network (JN):                                        C */
/*   ERRJN, GAUSJN, GJN, GPJN, JNCHOP, JNCOGR, JNCGBE, JNDELT, JNDUMP,  C */
/*   JNERR, JNFEED, JNHEAD, JNHEIG, JNHESS, JNINDX, JNINIT, JNLINS,     C */
/*   JNREAD, JNROLD, JNSATM, JNSCGR, JNSEFI, JNSEPA, JNSTAT, JNTEST,    C */
/*   JNTRAL, JNTRED, JNTQLI                                             C */
/*                                                                      C */
/* 2) Self-organizing network (JM):                                     C */
/*   GJM, JMDUMP, JMERR, JMFEED, JMINDX, JMINIT, JMINWE, JMNBHD,        C */
/*   JMNORM, JMREAD, JMSEPA, JMSTAT, JMTEST, JMTRAL, JMWARN             C */
/*                                                                      C */
/* The block-data subroutine JNDATA is placed at the end, together      C */
/* with a test-deck in the subroutine JNTDEC, and a random number       C */
/* generator called RJN.                                                C */
/* **********************************************************************C */
/* **********************************************************************C */
/* **********************************************************************C */
/* PART ONE: FEED-FORWARD NETWORK                                       C */
/* **********************************************************************C */
doublereal errjn_(integer *idum)
{
    /* System generated locals */
    integer i__1;
    real ret_val, r__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__;
    extern integer jnindx_(integer *, integer *, integer *);
    static real outval, err;

/* ...JetNet function calculate ERRor. */
/* ...Returns the error function. */
/* ...The error measure is selected by MSTJN(4). */
    err = 0.f;
    if (jndat1_1.mstjn[3] == 0) {
/* ...Summed square error: */
	i__1 = jnint2_1.m[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    outval = jnint1_1.o[jnindx_(&jnint2_1.nl, &i__, &c__0) - 1];
/* Computing 2nd power */
	    r__1 = jndat1_1.out[i__ - 1] - outval;
	    err += r__1 * r__1 * .5f;
	    jndat1_1.out[i__ - 1] = outval;
/* L100: */
	}
    } else if (jndat1_1.mstjn[3] == 1) {
/* ...Cross Entropy error: */
	i__1 = jnint2_1.m[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
/* ...It is assumed that OUTVAL=0.0 or OUTVAL=1.0 never happens. */
	    outval = jnint1_1.o[jnindx_(&jnint2_1.nl, &i__, &c__0) - 1];
	    err -= jndat1_1.out[i__ - 1] * log(outval) + (1.f - jndat1_1.out[
		    i__ - 1]) * log(1.f - outval);
	    jndat1_1.out[i__ - 1] = outval;
/* L110: */
	}
    } else if (jndat1_1.mstjn[3] >= 2) {
/* ...Kullback error: */
	i__1 = jnint2_1.m[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
/* ...It is assumed that OUTVAL=0.0 never happens. */
	    outval = jnint1_1.o[jnindx_(&jnint2_1.nl, &i__, &c__0) - 1];
	    if (jndat1_1.out[i__ - 1] > 0.f) {
		err += jndat1_1.out[i__ - 1] * log(jndat1_1.out[i__ - 1] / 
			outval);
	    }
	    jndat1_1.out[i__ - 1] = outval;
/* L120: */
	}
    } else if (jndat1_1.mstjn[3] == -1) {
/* ...Log-squared error: */
	i__1 = jnint2_1.m[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
/* ...It is assumed that |OUT(I)-OUTVAL|=1.0 never happens. */
	    outval = jnint1_1.o[jnindx_(&jnint2_1.nl, &i__, &c__0) - 1];
/* Computing 2nd power */
	    r__1 = jndat1_1.out[i__ - 1] - outval;
	    err -= log(1.f - r__1 * r__1) * .5f;
	    jndat1_1.out[i__ - 1] = outval;
/* L130: */
	}
    }
    ret_val = err;
    return ret_val;
/* **** END OF ERRJN ***************************************************** */
} /* errjn_ */

/* *********************************************************************** */
doublereal gausjn_(integer *idum)
{
    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);

    /* Local variables */
    static real r__, v1, v2, fac;
    extern doublereal rjn_(integer *);

/* ...JetNet function GAUSsian random number. */
/* ...Generates Gaussian distributed random numbers with */
/* ...standard deviation 1.0 and mean 0.0. Polar method. */
    if (jngaus_1.iset == 0) {
L100:
	v1 = rjn_(idum) * 2.f - 1.f;
	v2 = rjn_(idum) * 2.f - 1.f;
/* Computing 2nd power */
	r__1 = v1;
/* Computing 2nd power */
	r__2 = v2;
	r__ = r__1 * r__1 + r__2 * r__2;
	if (r__ >= 1.f || r__ <= 1e-20f) {
	    goto L100;
	}
/* ...Box-Muller transformation: */
	fac = sqrt(log(r__) * -2.f / r__);
	ret_val = v1 * fac;
	jngaus_1.gasdev = v2 * fac;
	jngaus_1.iset = 1;
    } else {
	ret_val = jngaus_1.gasdev;
	jngaus_1.iset = 0;
    }
    return ret_val;
/* **** END OF GAUSJN **************************************************** */
} /* gausjn_ */

/* *********************************************************************** */
doublereal gjn_(integer *ind, real *x, integer *n)
{
    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double tanh(doublereal), exp(doublereal);
    integer pow_ii(integer *, integer *);
    double r_sign(real *, real *);

    /* Local variables */
    static real g;
    extern /* Subroutine */ int jnerr_(integer *);
    static integer ng, ns;
    static real ss, arg;

/* ...JetNet function G */
/* ...Gives sigmoid function N with argument X */
/* ...The derivative GPrime is also calculated and stored in GPJN. */
    if (*n == 1) {
/* ...        1 -> g(x)=1/(1+exp(-2x)) */
	arg = tanh(*x);
	ret_val = (arg + 1.f) * .5f;
/* Computing 2nd power */
	r__1 = arg;
	jnsigm_1.gpjn[*ind - 1] = (1.f - r__1 * r__1) * .5f + jndat1_1.parjn[
		22];
/* Computing 2nd power */
	r__1 = arg;
	jnsigm_1.gppjn[*ind - 1] = arg * (r__1 * r__1 - 1.f);
    } else if (*n == 2) {
/* ...        2 -> g(x)=tanh(x) */
	arg = tanh(*x);
	ret_val = arg;
/* Computing 2nd power */
	r__1 = arg;
	jnsigm_1.gpjn[*ind - 1] = 1.f - r__1 * r__1 + jndat1_1.parjn[22];
/* Computing 2nd power */
	r__1 = arg;
	jnsigm_1.gppjn[*ind - 1] = arg * 2.f * (r__1 * r__1 - 1.f);
    } else if (*n == 3) {
/* ...        3 -> g(x)=exp(x) (only used internally for Potts-nodes) */
/* Computing MAX */
	r__1 = -50.f, r__2 = dmin(*x,50.f);
	ret_val = exp((dmax(r__1,r__2)));
	jnsigm_1.gpjn[*ind - 1] = 1.f;
	jnsigm_1.gppjn[*ind - 1] = 0.f;
    } else if (*n == 4) {
/* ...        4 -> g(x)=x */
	ret_val = *x;
	jnsigm_1.gpjn[*ind - 1] = 1.f;
	jnsigm_1.gppjn[*ind - 1] = 0.f;
    } else if (*n == 5) {
/* ...        5 -> g(x)=1/(1+exp(-2x)) (only used internally for */
/* ...             entropy error) */
	ret_val = (tanh(*x) + 1.f) * .5f;
	jnsigm_1.gpjn[*ind - 1] = 2.f;
	jnsigm_1.gppjn[*ind - 1] = 0.f;
    } else if (*n == -1) {
/* ...        same as above, but with fixed precision */
	arg = tanh(*x);
	i__1 = abs(jndat1_1.mstjn[27]);
	ns = pow_ii(&c__2, &i__1);
	ss = 1.f / (ns - 1);
	g = (arg + 1.f) * .5f;
	ng = (integer) (g / ss + .5f);
	ret_val = (real) ng * ss;
/* Computing 2nd power */
	r__1 = arg;
	jnsigm_1.gpjn[*ind - 1] = (1.f - r__1 * r__1) * .5f + jndat1_1.parjn[
		22];
/* Computing 2nd power */
	r__1 = arg;
	jnsigm_1.gppjn[*ind - 1] = arg * (r__1 * r__1 - 1.f);
    } else if (*n == -2) {
	arg = tanh(*x);
	i__1 = abs(jndat1_1.mstjn[27]) - 1;
	ns = pow_ii(&c__2, &i__1);
	if (ns == 1) {
	    ret_val = r_sign(&c_b15, &arg);
	} else {
	    ss = 1.f / (ns - 1);
	    g = arg;
	    ng = (integer) (dabs(g) / ss + .5f);
	    r__1 = (real) ng * ss;
	    ret_val = r_sign(&r__1, &g);
	}
/* Computing 2nd power */
	r__1 = arg;
	jnsigm_1.gpjn[*ind - 1] = 1.f - r__1 * r__1 + jndat1_1.parjn[22];
/* Computing 2nd power */
	r__1 = arg;
	jnsigm_1.gppjn[*ind - 1] = arg * 2.f * (r__1 * r__1 - 1.f);
    } else {
	jndat1_1.mstjn[2] = *n;
	jnerr_(&c__15);
    }
    return ret_val;
/* **** END OF GJN ******************************************************* */
} /* gjn_ */

/* *********************************************************************** */
/* Subroutine */ int jncgbe_(real *betak, integer *iop)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static real beta, dgdw, fact1, fact2, oldg2;
    static integer i__;
    static real y, f1, f2, f3, f4;
    static integer il;

/* ...JetNet subroutine Conjugate Gradient BEta_k. */
/* ...If IOP=0, it only saves current value of (DW,DT)*(DW,DT) to be */
/* ...used later and returns BETAK=0.0. If IOP=1, it also calculates */
/* ...the beta_k value used to generate next search direction. Which */
/* ...formula to use is determined by MSTJN(5). */
/* ...Note: The vector (DW,DT) equals the negative gradient of E. */
    beta = 1.f;
    *betak = 0.f;
    oldg2 = jnint4_1.g2;
    jnint4_1.g2 = 0.f;
    if (jndat1_1.mstjn[4] == 4 || jndat1_1.mstjn[4] == 10) {
/* ...Polak-Ribiere's formula: */
	for (il = jnint2_1.nl; il >= 1; --il) {
	    if (jndat2_1.tinv[il - 1] == 0.f) {
		beta *= jndat1_1.parjn[2];
	    } else {
		beta *= (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	    }
	    i__1 = jnint2_1.mm0[il];
	    for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.dw[i__ - 1] /= (real) jndat1_1.mstjn[1];
/* Computing 2nd power */
		r__1 = jnint1_1.dw[i__ - 1] * beta;
		jnint4_1.g2 += (real) jnint1_1.nself[i__ - 1] * (r__1 * r__1);
/* Computing 2nd power */
		r__1 = beta;
		*betak += (real) jnint1_1.nself[i__ - 1] * jnint1_1.dw[i__ - 
			1] * jnint1_1.odw[i__ - 1] * (r__1 * r__1);
		jnint1_1.odw[i__ - 1] = jnint1_1.dw[i__ - 1];
/* L110: */
	    }
	    i__1 = jnint2_1.mv0[il];
	    for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.dt[i__ - 1] /= (real) jndat1_1.mstjn[1];
/* Computing 2nd power */
		r__1 = jnint1_1.dt[i__ - 1] * beta;
		jnint4_1.g2 += (real) jnint1_1.ntself[i__ - 1] * (r__1 * r__1)
			;
/* Computing 2nd power */
		r__1 = beta;
		*betak += (real) jnint1_1.ntself[i__ - 1] * jnint1_1.dt[i__ - 
			1] * jnint1_1.odt[i__ - 1] * (r__1 * r__1);
		jnint1_1.odt[i__ - 1] = jnint1_1.dt[i__ - 1];
/* L120: */
	    }
/* L100: */
	}
	if (*iop == 0) {
	    *betak = 0.f;
	} else {
	    if (dabs(oldg2) > 1e-8f) {
		*betak = (jnint4_1.g2 - *betak) / oldg2;
	    } else {
		*betak = 0.f;
	    }
	}
    } else if (jndat1_1.mstjn[4] == 5 || jndat1_1.mstjn[4] == 11) {
/* ...Hestenes-Stiefel's formula: */
	dgdw = 0.f;
	for (il = jnint2_1.nl; il >= 1; --il) {
	    if (jndat2_1.tinv[il - 1] == 0.f) {
		beta *= jndat1_1.parjn[2];
	    } else {
		beta *= (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	    }
	    i__1 = jnint2_1.mm0[il];
	    for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.dw[i__ - 1] /= (real) jndat1_1.mstjn[1];
/* Computing 2nd power */
		r__1 = jnint1_1.dw[i__ - 1] * beta;
		jnint4_1.g2 += (real) jnint1_1.nself[i__ - 1] * (r__1 * r__1);
/* Computing 2nd power */
		r__1 = beta;
		*betak += (real) jnint1_1.nself[i__ - 1] * jnint1_1.dw[i__ - 
			1] * jnint1_1.odw[i__ - 1] * (r__1 * r__1);
/* Computing 2nd power */
		r__1 = beta;
		dgdw += jnint1_1.g[i__ - 1] * (jnint1_1.odw[i__ - 1] - 
			jnint1_1.dw[i__ - 1]) * (r__1 * r__1);
		jnint1_1.odw[i__ - 1] = jnint1_1.dw[i__ - 1];
/* L210: */
	    }
	    i__1 = jnint2_1.mv0[il];
	    for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.dt[i__ - 1] /= (real) jndat1_1.mstjn[1];
/* Computing 2nd power */
		r__1 = jnint1_1.dt[i__ - 1] * beta;
		jnint4_1.g2 += (real) jnint1_1.ntself[i__ - 1] * (r__1 * r__1)
			;
/* Computing 2nd power */
		r__1 = beta;
		*betak += (real) jnint1_1.ntself[i__ - 1] * jnint1_1.dt[i__ - 
			1] * jnint1_1.odt[i__ - 1] * (r__1 * r__1);
/* Computing 2nd power */
		r__1 = beta;
		dgdw += jnint1_1.g[jnint2_1.mm0[jnint2_1.nl] + i__ - 1] * (
			jnint1_1.odt[i__ - 1] - jnint1_1.dt[i__ - 1]) * (r__1 
			* r__1);
		jnint1_1.odt[i__ - 1] = jnint1_1.dt[i__ - 1];
/* L220: */
	    }
/* L200: */
	}
	if (*iop == 0) {
	    *betak = 0.f;
	} else {
	    if (dabs(dgdw) > 1e-8f) {
		*betak = (jnint4_1.g2 - *betak) / dgdw;
	    } else {
		*betak = 0.f;
	    }
	}
    } else if (jndat1_1.mstjn[4] == 6 || jndat1_1.mstjn[4] == 12) {
/* ...Fletcher-Reeves' formula: */
	for (il = jnint2_1.nl; il >= 1; --il) {
	    if (jndat2_1.tinv[il - 1] == 0.f) {
		beta *= jndat1_1.parjn[2];
	    } else {
		beta *= (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	    }
	    i__1 = jnint2_1.mm0[il];
	    for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.dw[i__ - 1] /= (real) jndat1_1.mstjn[1];
/* Computing 2nd power */
		r__1 = jnint1_1.dw[i__ - 1] * beta;
		jnint4_1.g2 += (real) jnint1_1.nself[i__ - 1] * (r__1 * r__1);
		jnint1_1.odw[i__ - 1] = jnint1_1.dw[i__ - 1];
/* L310: */
	    }
	    i__1 = jnint2_1.mv0[il];
	    for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.dt[i__ - 1] /= (real) jndat1_1.mstjn[1];
/* Computing 2nd power */
		r__1 = jnint1_1.dt[i__ - 1] * beta;
		jnint4_1.g2 += (real) jnint1_1.ntself[i__ - 1] * (r__1 * r__1)
			;
		jnint1_1.odt[i__ - 1] = jnint1_1.dt[i__ - 1];
/* L320: */
	    }
/* L300: */
	}
	if (*iop == 0) {
	    *betak = 0.f;
	} else {
	    if (dabs(oldg2) > 1e-8f) {
		*betak = jnint4_1.g2 / oldg2;
	    } else {
		*betak = 0.f;
	    }
	}
    } else if (jndat1_1.mstjn[4] == 7 || jndat1_1.mstjn[4] == 13) {
/* ...Shanno's formula: */
	f1 = 0.f;
	f2 = 0.f;
	f3 = 0.f;
	f4 = 0.f;
	fact1 = 0.f;
	fact2 = 0.f;
	for (il = jnint2_1.nl; il >= 1; --il) {
	    if (jndat2_1.tinv[il - 1] == 0.f) {
		beta *= jndat1_1.parjn[2];
	    } else {
		beta *= (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	    }
	    i__1 = jnint2_1.mm0[il];
	    for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.dw[i__ - 1] /= (real) jndat1_1.mstjn[1];
/* Computing 2nd power */
		r__1 = jnint1_1.dw[i__ - 1] * beta;
		jnint4_1.g2 += (real) jnint1_1.nself[i__ - 1] * (r__1 * r__1);
		f1 -= jnint1_1.g[i__ - 1] * jnint1_1.odw[i__ - 1] * beta;
		f2 += jnint1_1.g[i__ - 1] * (jnint1_1.odw[i__ - 1] - 
			jnint1_1.dw[i__ - 1]) * beta;
/* Computing 2nd power */
		r__1 = beta;
		f3 += jnint1_1.odw[i__ - 1] * (jnint1_1.dw[i__ - 1] - 
			jnint1_1.odw[i__ - 1]) * (r__1 * r__1);
/* Computing 2nd power */
		r__1 = (jnint1_1.dw[i__ - 1] - jnint1_1.odw[i__ - 1]) * beta;
		f4 += r__1 * r__1;
		y = jnint1_1.odw[i__ - 1] - jnint1_1.dw[i__ - 1];
		jnint1_1.odw[i__ - 1] = jnint1_1.dw[i__ - 1];
		jnint1_1.dw[i__ - 1] = y;
/* L410: */
	    }
	    i__1 = jnint2_1.mv0[il];
	    for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.dt[i__ - 1] /= (real) jndat1_1.mstjn[1];
/* Computing 2nd power */
		r__1 = jnint1_1.dt[i__ - 1] * beta;
		jnint4_1.g2 += (real) jnint1_1.ntself[i__ - 1] * (r__1 * r__1)
			;
		f1 -= jnint1_1.g[jnint2_1.mm0[jnint2_1.nl] + i__ - 1] * 
			jnint1_1.odt[i__ - 1] * beta;
		f2 += jnint1_1.g[jnint2_1.mm0[jnint2_1.nl] + i__ - 1] * (
			jnint1_1.odt[i__ - 1] - jnint1_1.dt[i__ - 1]) * beta;
/* Computing 2nd power */
		r__1 = beta;
		f3 += jnint1_1.odt[i__ - 1] * (jnint1_1.dt[i__ - 1] - 
			jnint1_1.odt[i__ - 1]) * (r__1 * r__1);
/* Computing 2nd power */
		r__1 = (jnint1_1.dt[i__ - 1] - jnint1_1.odt[i__ - 1]) * beta;
		f4 += r__1 * r__1;
		y = jnint1_1.odt[i__ - 1] - jnint1_1.dt[i__ - 1];
		jnint1_1.odt[i__ - 1] = jnint1_1.dt[i__ - 1];
		jnint1_1.dt[i__ - 1] = y;
/* L420: */
	    }
	    f1 *= jnint4_1.stepln[0];
	    f2 *= jnint4_1.stepln[0];
	    if (dabs(f2) > 1e-8f && dabs(f2) < 1e8f) {
		fact1 = f1 / f2;
		fact2 = (f4 / f2 + 1.f) * fact1 - f3 / f2;
	    } else {
		fact1 = 0.f;
		fact2 = 0.f;
	    }
/* L400: */
	}
	if (*iop == 0) {
	    *betak = 0.f;
	    i__1 = jnint2_1.mm0[jnint2_1.nl];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jnint1_1.dw[i__ - 1] = jnint1_1.odw[i__ - 1];
/* L430: */
	    }
	    i__1 = jnint2_1.mv0[jnint2_1.nl];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jnint1_1.dt[i__ - 1] = jnint1_1.odt[i__ - 1];
/* L440: */
	    }
	} else {
	    *betak = -fact2 * jnint4_1.stepln[0];
	    i__1 = jnint2_1.mm0[jnint2_1.nl];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jnint1_1.dw[i__ - 1] = jnint1_1.odw[i__ - 1] + fact1 * 
			jnint1_1.dw[i__ - 1];
/* L450: */
	    }
	    i__1 = jnint2_1.mv0[jnint2_1.nl];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jnint1_1.dt[i__ - 1] = jnint1_1.odt[i__ - 1] + fact1 * 
			jnint1_1.dt[i__ - 1];
/* L460: */
	    }
	}
    }
    return 0;
/* **** END OF JNCGBE **************************************************** */
} /* jncgbe_ */

/* *********************************************************************** */
/* Subroutine */ int jnchop_(integer *ichp)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2, r__3;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);
    double r_sign(real *, real *);

    /* Local variables */
    static real tmax, wmax;
    static integer i__;
    extern /* Subroutine */ int jnerr_(integer *);
    static integer il;
    static real ss;
    static integer nst, nsw;

/* ...JetNet subroutine CHOP weights */
/* ...Switches on (ICHP>0) or off (ICHP<0) fixed precision weights */
/* ...thresholds and sigmoid functions. For IHCP >= 0 the weights and */
/* ...thresholds are chopped to the fixed precision. */
    if (*ichp > 0) {
	jnint2_1.icpon = 1;
	i__1 = jnint2_1.nl - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if ((i__2 = jnint2_1.ng[i__ - 1], abs(i__2)) != 1 && (i__3 = 
		    jnint2_1.ng[i__ - 1], abs(i__3)) != 2) {
		jnerr_(&c__14);
	    }
	    if (jndat1_1.mstjn[27] != 0) {
		jnint2_1.ng[i__ - 1] = -(i__2 = jnint2_1.ng[i__ - 1], abs(
			i__2));
	    } else {
		jnint2_1.ng[i__ - 1] = (i__2 = jnint2_1.ng[i__ - 1], abs(i__2)
			);
	    }
/* L100: */
	}
	if (jndat1_1.mstjn[27] > 0) {
	    if ((i__1 = jnint2_1.ng[jnint2_1.nl - 1], abs(i__1)) != 1 && (
		    i__2 = jnint2_1.ng[jnint2_1.nl - 1], abs(i__2)) != 2) {
		jnerr_(&c__14);
	    }
	    jnint2_1.ng[jnint2_1.nl - 1] = -(i__1 = jnint2_1.ng[jnint2_1.nl - 
		    1], abs(i__1));
	} else {
	    jnint2_1.ng[jnint2_1.nl - 1] = (i__1 = jnint2_1.ng[jnint2_1.nl - 
		    1], abs(i__1));
	}
    }
    if (*ichp < 0) {
	jnint2_1.icpon = 0;
	i__1 = jnint2_1.nl;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint2_1.ng[i__ - 1] = (i__2 = jnint2_1.ng[i__ - 1], abs(i__2));
/* L110: */
	}
    }
    if (*ichp >= 0) {
	i__1 = jndat1_1.mstjn[29] - 1;
	nsw = pow_ii(&c__2, &i__1);
	i__1 = jndat1_1.mstjn[28] - 1;
	nst = pow_ii(&c__2, &i__1);
	i__1 = jnint2_1.nl;
	for (il = 1; il <= i__1; ++il) {
	    if (jndat1_1.mstjn[29] >= 1) {
		wmax = 0.f;
		i__2 = jnint2_1.mm0[il];
		for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		    r__2 = wmax, r__3 = (r__1 = jnint1_1.w[i__ - 1], dabs(
			    r__1));
		    wmax = dmax(r__2,r__3);
/* L210: */
		}
		ss = 1.f;
		if (nsw > 1) {
		    ss = wmax / (nsw - 1);
		}
		i__2 = jnint2_1.mm0[il];
		for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__2; ++i__) {
		    if (nsw > 1) {
			r__2 = (real) ((integer) ((r__1 = jnint1_1.w[i__ - 1],
				 dabs(r__1)) / ss + .5f)) * ss;
			jnint1_1.w[i__ - 1] = r_sign(&r__2, &jnint1_1.w[i__ - 
				1]);
		    } else {
			jnint1_1.w[i__ - 1] = r_sign(&wmax, &jnint1_1.w[i__ - 
				1]);
		    }
/* L220: */
		}
	    }
	    if (jndat1_1.mstjn[28] >= 1) {
		tmax = 0.f;
		i__2 = jnint2_1.mv0[il];
		for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		    r__2 = tmax, r__3 = (r__1 = jnint1_1.t[i__ - 1], dabs(
			    r__1));
		    tmax = dmax(r__2,r__3);
/* L230: */
		}
		ss = 1.f;
		if (nst > 1) {
		    ss = tmax / (nst - 1);
		}
		i__2 = jnint2_1.mv0[il];
		for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__2; ++i__) {
		    if (nst > 1) {
			r__2 = (real) ((integer) ((r__1 = jnint1_1.t[i__ - 1],
				 dabs(r__1)) / ss + .5f)) * ss;
			jnint1_1.t[i__ - 1] = r_sign(&r__2, &jnint1_1.t[i__ - 
				1]);
		    } else {
			jnint1_1.t[i__ - 1] = r_sign(&tmax, &jnint1_1.t[i__ - 
				1]);
		    }
/* L240: */
		}
	    }
/* L200: */
	}
    }
    return 0;
/* **** END OF JNCHOP **************************************************** */
} /* jnchop_ */

/* *********************************************************************** */
/* Subroutine */ int jncogr_(void)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static real beta;
    static integer i__;
    static real betak;
    static integer il;
    extern /* Subroutine */ int jncgbe_(real *, integer *), jnlins_(void);

/* ...JetNet subroutine COnjugate GRadient */
/* ...Performs Conjugate Gradient updating. */
    if (jnint4_1.ilinon == 0) {
/* ...Calc. new conjugate search direction + calculate gradient: */
	i__1 = jnint4_1.nc % (jnint2_1.mm0[jnint2_1.nl] + jnint2_1.mv0[
		jnint2_1.nl]);
	jncgbe_(&betak, &i__1);
	++jnint4_1.nc;
	jnint4_1.nsc = 0;
	jnint4_1.derrln = 0.f;
	beta = 1.f;
	jnint4_1.stepmn = 0.f;
	for (il = jnint2_1.nl; il >= 1; --il) {
/* ...set effective beta in layer IL: */
	    if (jndat2_1.tinv[il - 1] == 0.f) {
		beta *= jndat1_1.parjn[2];
	    } else {
		beta *= (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	    }
	    i__1 = jnint2_1.mm0[il];
	    for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.g[i__ - 1] = betak * jnint1_1.g[i__ - 1] + 
			jnint1_1.dw[i__ - 1] * (real) jnint1_1.nself[i__ - 1] 
			* beta;
		jnint4_1.derrln -= jnint1_1.odw[i__ - 1] * (real) 
			jnint1_1.nself[i__ - 1] * beta * jnint1_1.g[i__ - 1];
/* L110: */
	    }
	    i__1 = jnint2_1.mv0[il];
	    for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.g[i__ + jnint2_1.mm0[jnint2_1.nl] - 1] = betak * 
			jnint1_1.g[i__ + jnint2_1.mm0[jnint2_1.nl] - 1] + 
			jnint1_1.dt[i__ - 1] * (real) jnint1_1.ntself[i__ - 1]
			 * beta;
		jnint4_1.derrln -= jnint1_1.odt[i__ - 1] * (real) 
			jnint1_1.ntself[i__ - 1] * beta * jnint1_1.g[i__ + 
			jnint2_1.mm0[jnint2_1.nl] - 1];
/* L120: */
	    }
/* L100: */
	}
	jnint4_1.ilinon = 1;
	jnint4_1.nit = 0;
	jnlins_();
    } else {
/* ...Do line search */
	jnlins_();
	if (jnint4_1.ilinon == 0) {
/* ...Zero (DW,DT) */
	    i__1 = jnint2_1.mm0[jnint2_1.nl];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jnint1_1.dw[i__ - 1] = 0.f;
/* L130: */
	    }
	    i__1 = jnint2_1.mv0[jnint2_1.nl];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jnint1_1.dt[i__ - 1] = 0.f;
/* L140: */
	    }
	}
    }
    return 0;
/* **** END OF JNCOGR **************************************************** */
} /* jncogr_ */

/* *********************************************************************** */
/* Subroutine */ int jndelt_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1;

    /* Local variables */
    static real diff;
    static integer i__, j, ihprf, nextl, ih, mi, il, mj, in, ix, iy, jy, iw, 
	    jx;
    static real sumrft, sumrfw;
    static integer mij, inx, iny;
    static real sum;

/* ...JetNet subroutine DELTa weights */
/* ...Calculates the change in weights and thresholds to minimize the */
/* ...cost function according to gradient descent */
/* ...(Learning rate and inverse temperature are multiplied in JNTRAL). */
/* ...calculate the deltas at nodes in output layer */
    if (jndat1_1.mstjn[3] == -1) {
	i__1 = jnint2_1.m[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    mi = jnint2_1.mv0[jnint2_1.nl - 1] + i__;
	    diff = jndat1_1.out[i__ - 1] - jnint1_1.o[mi - 1];
/* Computing 2nd power */
	    r__1 = diff;
	    jnint1_1.d__[mi - 1] = diff * jnsigm_1.gpjn[mi - 1] / (1.f - r__1 
		    * r__1);
/* L101: */
	}
    } else {
	i__1 = jnint2_1.m[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    mi = jnint2_1.mv0[jnint2_1.nl - 1] + i__;
	    jnint1_1.d__[mi - 1] = (jndat1_1.out[i__ - 1] - jnint1_1.o[mi - 1]
		    ) * jnsigm_1.gpjn[mi - 1];
/* L100: */
	}
    }
/* ...calculate deltas in following layers */
    for (il = jnint2_1.nl - 1; il >= 1; --il) {
/* ...calculate the deltas at nodes in layer IL */
	i__1 = jnint2_1.m[il];
	for (j = 1; j <= i__1; ++j) {
	    mj = jnint2_1.mv0[il - 1] + j;
	    mij = jnint2_1.mm0[il] + (j - 1) * jnint2_1.m[il + 1];
	    sum = 0.f;
	    i__2 = jnint2_1.mv0[il] + jnint2_1.m[il + 1];
	    for (i__ = jnint2_1.mv0[il] + 1; i__ <= i__2; ++i__) {
		++mij;
		sum += jnint1_1.d__[i__ - 1] * jnint1_1.w[mij - 1];
/* L220: */
	    }
	    jnint1_1.d__[mj - 1] = sum * jnsigm_1.gpjn[mj - 1];
/* L210: */
	}
/* L200: */
    }
/* ...calculate deltas at all weights between first and second layer */
    nextl = 2;
    if (jnint3_1.nxin == 0) {
/* ...normal first layer */
	i__1 = jnint2_1.m[1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    mij = i__ - jnint2_1.m[1];
	    i__2 = jnint2_1.m[0];
	    for (j = 1; j <= i__2; ++j) {
		mij += jnint2_1.m[1];
		jnint1_1.dw[mij - 1] += jnint1_1.d__[i__ - 1] * jndat1_1.oin[
			j - 1];
/* L310: */
	    }
	    jnint1_1.dt[i__ - 1] += jnint1_1.d__[i__ - 1];
/* L300: */
	}
    } else {
/* ...receptive fields in first layer */
	i__1 = jnint3_1.nhprf;
	for (ihprf = 1; ihprf <= i__1; ++ihprf) {
	    sumrft = 0.f;
	    i__2 = jnint3_1.nyhrf;
	    for (iy = 1; iy <= i__2; ++iy) {
		ih = iy - jnint3_1.nyhrf + (ihprf - 1) * jnint3_1.nhrf;
		i__3 = jnint3_1.nxhrf;
		for (ix = 1; ix <= i__3; ++ix) {
		    ih += jnint3_1.nyhrf;
		    i__4 = jnint3_1.nyrf;
		    for (jy = 1; jy <= i__4; ++jy) {
			iw = jy - jnint3_1.nyrf + (ihprf - 1) * jnint3_1.nrfw;
			i__5 = jnint3_1.nxrf;
			for (jx = 1; jx <= i__5; ++jx) {
			    iw += jnint3_1.nyrf;
			    inx = ix + jx - 1;
			    if (inx > abs(jnint3_1.nxin)) {
				inx -= abs(jnint3_1.nxin);
			    }
			    iny = iy + jy - 1;
			    if (iny > abs(jnint3_1.nyin)) {
				iny -= abs(jnint3_1.nyin);
			    }
			    in = (inx - 1) * abs(jnint3_1.nyin) + iny;
			    jnint1_1.dw[iw - 1] += jnint1_1.d__[ih - 1] * 
				    jndat1_1.oin[in - 1] / (real) 
				    jnint3_1.nhrf;
/* L360: */
			}
/* L350: */
		    }
		    i__5 = jnint2_1.m[0];
		    for (in = (i__4 = jnint3_1.nxin * jnint3_1.nyin, abs(i__4)
			    ) + 1; in <= i__5; ++in) {
			iw = jnint3_1.nxrf * jnint3_1.nyrf + in - (i__4 = 
				jnint3_1.nxin * jnint3_1.nyin, abs(i__4)) + (
				ihprf - 1) * jnint3_1.nrfw;
			jnint1_1.dw[iw - 1] += jnint1_1.d__[ih - 1] * 
				jndat1_1.oin[in - 1] / (real) jnint3_1.nhrf;
/* L370: */
		    }
		    sumrft += jnint1_1.d__[ih - 1];
/* L340: */
		}
/* L330: */
	    }
	    sumrft /= (real) (jnint3_1.nxhrf * jnint3_1.nyhrf);
	    i__2 = jnint3_1.nxhrf * jnint3_1.nyhrf;
	    for (ih = 1; ih <= i__2; ++ih) {
		jnint1_1.dt[ih + (ihprf - 1) * jnint3_1.nhrf - 1] += sumrft;
/* L380: */
	    }
/* L320: */
	}
	i__1 = jnint2_1.m[1];
	for (ih = jnint3_1.nhrf * jnint3_1.nhprf + 1; ih <= i__1; ++ih) {
	    iw = jnint3_1.nrfw * jnint3_1.nhprf + ih - jnint2_1.m[1];
	    i__2 = jnint2_1.m[0];
	    for (in = 1; in <= i__2; ++in) {
		iw += jnint2_1.m[1];
		jnint1_1.dw[iw - 1] += jnint1_1.d__[ih - 1] * jndat1_1.oin[in 
			- 1];
/* L400: */
	    }
	    jnint1_1.dt[ih - 1] += jnint1_1.d__[ih - 1];
/* L390: */
	}
	if (jndat1_1.mstjn[26] < 0) {
	    i__1 = jnint2_1.m[2];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		mi = jnint2_1.mv0[1] + i__;
		i__2 = jnint3_1.nhprf;
		for (ihprf = 1; ihprf <= i__2; ++ihprf) {
		    sumrfw = 0.f;
		    i__3 = ihprf * jnint3_1.nhrf;
		    for (j = (ihprf - 1) * jnint3_1.nhrf + 1; j <= i__3; ++j) 
			    {
			sumrfw += jnint1_1.d__[mi - 1] * jnint1_1.o[j - 1];
/* L520: */
		    }
		    mij = jnint2_1.mm0[1] + ((ihprf - 1) * jnint3_1.nhrf - 1) 
			    * jnint2_1.m[2] + i__;
		    sumrfw /= (real) jnint3_1.nhrf;
		    i__3 = jnint3_1.nhrf;
		    for (j = 1; j <= i__3; ++j) {
			mij += jnint2_1.m[2];
			jnint1_1.dw[mij - 1] += sumrfw;
/* L530: */
		    }
/* L510: */
		}
		mij = jnint2_1.mm0[1] + i__ + (jnint3_1.nhrf * jnint3_1.nhprf 
			- 1) * jnint2_1.m[2];
		i__2 = jnint2_1.m[1];
		for (j = jnint3_1.nhrf * jnint3_1.nhprf + 1; j <= i__2; ++j) {
		    mij += jnint2_1.m[2];
		    jnint1_1.dw[mij - 1] += jnint1_1.d__[mi - 1] * jnint1_1.o[
			    j - 1];
/* L540: */
		}
		jnint1_1.dt[mi - 1] += jnint1_1.d__[mi - 1];
/* L500: */
	    }
	    nextl = 3;
	}
    }
/* ...calculate deltas at all weights between following layers */
    i__1 = jnint2_1.nl;
    for (il = nextl; il <= i__1; ++il) {
	i__2 = jnint2_1.m[il];
	for (i__ = 1; i__ <= i__2; ++i__) {
	    mij = jnint2_1.mm0[il - 1] + i__ - jnint2_1.m[il];
	    mi = jnint2_1.mv0[il - 1] + i__;
	    i__3 = jnint2_1.mv0[il - 2] + jnint2_1.m[il - 1];
	    for (j = jnint2_1.mv0[il - 2] + 1; j <= i__3; ++j) {
		mij += jnint2_1.m[il];
		jnint1_1.dw[mij - 1] += jnint1_1.d__[mi - 1] * jnint1_1.o[j - 
			1];
/* L430: */
	    }
	    jnint1_1.dt[mi - 1] += jnint1_1.d__[mi - 1];
/* L420: */
	}
/* L410: */
    }
    return 0;
/* **** END OF JNDELT **************************************************** */
} /* jndelt_ */

/* *********************************************************************** */
/* Subroutine */ int jndump_(integer *nf)
{
    /* Format strings */
    static char fmt_600[] = "(26x,\002Dump of weights generated by\002)";
    static char fmt_610[] = "(21x,\002Values of weights between layer\002,"
	    "i2,\002 (rows) and\002,i2,\002 (columns)\002)";
    static char fmt_640[] = "(6(e12.4,a1))";
    static char fmt_650[] = "(21x,\002Values of weights in receptive field"
	    "s\002)";
    static char fmt_660[] = "(21x,\002Values of other weights between input "
	    "layer (rows)\002,\002 and layer  1 (columns)\002)";
    static char fmt_630[] = "(30x,\002Thresholds in layer\002,i2)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsue(cilist *), do_uio(integer *, char *, ftnlen), e_wsue(void),
	     s_wsfe(cilist *), e_wsfe(void), s_wsle(cilist *), e_wsle(void), 
	    do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static char c__[1*2];
    static integer i__, j, ihprf, jf, il;
    extern /* Subroutine */ int jnhead_(void);
    static integer iw, nfsave;
    extern integer jnindx_(integer *, integer *, integer *);
    extern /* Subroutine */ int jnstat_(integer *);

    /* Fortran I/O blocks */
    static cilist io___58 = { 0, 0, 0, 0, 0 };
    static cilist io___59 = { 0, 0, 0, 0, 0 };
    static cilist io___61 = { 0, 0, 0, 0, 0 };
    static cilist io___62 = { 0, 0, 0, 0, 0 };
    static cilist io___63 = { 0, 0, 0, 0, 0 };
    static cilist io___64 = { 0, 0, 0, 0, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___68 = { 0, 0, 0, 0, 0 };
    static cilist io___69 = { 0, 0, 0, 0, 0 };
    static cilist io___70 = { 0, 0, 0, 0, 0 };
    static cilist io___71 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___73 = { 0, 0, 0, 0, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___75 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___77 = { 0, 0, 0, 0, 0 };
    static cilist io___78 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___80 = { 0, 0, 0, 0, 0 };
    static cilist io___81 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___82 = { 0, 0, 0, 0, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___84 = { 0, 0, 0, 0, 0 };
    static cilist io___85 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___86 = { 0, 0, 0, 0, 0 };
    static cilist io___87 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___89 = { 0, 0, 0, 0, 0 };
    static cilist io___90 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___91 = { 0, 0, 0, 0, 0 };
    static cilist io___92 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___93 = { 0, 0, 0, 0, 0 };
    static cilist io___94 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___95 = { 0, 0, 0, 0, 0 };
    static cilist io___96 = { 0, 0, 0, fmt_640, 0 };


/* ...JetNet subroutine DUMP weights */
/* ...Dumps weights, threshold and other characteristics of the */
/* ...net to a file for use in other programs */
    if (*nf < 0) {
/* ...unformatted dump */
	jf = -(*nf);
	io___58.ciunit = jf;
	s_wsue(&io___58);
	do_uio(&c__1, (char *)&c__30, (ftnlen)sizeof(integer));
	e_wsue();
	io___59.ciunit = jf;
	s_wsue(&io___59);
	do_uio(&c__40, (char *)&jndat1_1.mstjn[0], (ftnlen)sizeof(integer));
	do_uio(&c__40, (char *)&jndat1_1.parjn[0], (ftnlen)sizeof(real));
	do_uio(&c__10, (char *)&jndat2_1.tinv[0], (ftnlen)sizeof(real));
	do_uio(&c__10, (char *)&jndat2_1.igfn[0], (ftnlen)sizeof(integer));
	do_uio(&c__10, (char *)&jndat2_1.etal[0], (ftnlen)sizeof(real));
	do_uio(&c__10, (char *)&jndat2_1.widl[0], (ftnlen)sizeof(real));
	do_uio(&c__10, (char *)&jndat2_1.satm[0], (ftnlen)sizeof(real));
	e_wsue();
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___61.ciunit = jf;
	    s_wsue(&io___61);
	    do_uio(&c__1, (char *)&jnint1_1.w[i__ - 1], (ftnlen)sizeof(real));
	    e_wsue();
/* L100: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___62.ciunit = jf;
	    s_wsue(&io___62);
	    do_uio(&c__1, (char *)&jnint1_1.t[i__ - 1], (ftnlen)sizeof(real));
	    e_wsue();
/* L110: */
	}
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___63.ciunit = jf;
	    s_wsue(&io___63);
	    do_uio(&c__1, (char *)&jnint1_1.nself[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_wsue();
/* L120: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___64.ciunit = jf;
	    s_wsue(&io___64);
	    do_uio(&c__1, (char *)&jnint1_1.ntself[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_wsue();
/* L130: */
	}
    } else {
/* ...Formatted dump */
	*(unsigned char *)&c__[1] = ' ';
	*(unsigned char *)&c__[0] = '*';
	nfsave = jndat1_1.mstjn[5];
	jndat1_1.mstjn[5] = *nf;
	io___67.ciunit = *nf;
	s_wsfe(&io___67);
	e_wsfe();
	jnhead_();
	jnstat_(&c__1);
	jnstat_(&c__2);
	jndat1_1.mstjn[5] = nfsave;
	io___68.ciunit = *nf;
	s_wsle(&io___68);
	e_wsle();
	io___69.ciunit = *nf;
	s_wsle(&io___69);
	e_wsle();
	io___70.ciunit = *nf;
	s_wsle(&io___70);
	e_wsle();
	if (jnint3_1.nxin == 0) {
	    io___71.ciunit = *nf;
	    s_wsfe(&io___71);
	    do_fio(&c__1, (char *)&c__0, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&c__1, (ftnlen)sizeof(integer));
	    e_wsfe();
	    i__1 = jnint2_1.m[0];
	    for (j = 1; j <= i__1; ++j) {
		io___73.ciunit = *nf;
		s_wsle(&io___73);
		e_wsle();
		io___74.ciunit = *nf;
		s_wsfe(&io___74);
		i__2 = jnint2_1.m[1];
		for (i__ = 1; i__ <= i__2; ++i__) {
		    do_fio(&c__1, (char *)&jnint1_1.w[jnindx_(&c__1, &i__, &j)
			     - 1], (ftnlen)sizeof(real));
		    do_fio(&c__1, c__ + jnint1_1.nself[jnindx_(&c__1, &i__, &
			    j) - 1], (ftnlen)1);
		}
		e_wsfe();
/* L200: */
	    }
	} else {
	    io___75.ciunit = *nf;
	    s_wsfe(&io___75);
	    e_wsfe();
	    i__1 = jnint3_1.nhprf;
	    for (ihprf = 1; ihprf <= i__1; ++ihprf) {
		io___77.ciunit = *nf;
		s_wsle(&io___77);
		e_wsle();
		io___78.ciunit = *nf;
		s_wsfe(&io___78);
		i__2 = jnint3_1.nrfw * ihprf;
		for (iw = jnint3_1.nrfw * (ihprf - 1) + 1; iw <= i__2; ++iw) {
		    do_fio(&c__1, (char *)&jnint1_1.w[iw - 1], (ftnlen)sizeof(
			    real));
		    do_fio(&c__1, c__ + jnint1_1.nself[iw - 1], (ftnlen)1);
		}
		e_wsfe();
/* L210: */
	    }
	    if (jnint3_1.nhrf * jnint3_1.nhprf < jnint2_1.m[1]) {
		io___80.ciunit = *nf;
		s_wsle(&io___80);
		e_wsle();
		io___81.ciunit = *nf;
		s_wsfe(&io___81);
		e_wsfe();
		i__1 = jnint2_1.m[0];
		for (j = 1; j <= i__1; ++j) {
		    io___82.ciunit = *nf;
		    s_wsle(&io___82);
		    e_wsle();
		    io___83.ciunit = *nf;
		    s_wsfe(&io___83);
		    i__2 = jnint2_1.m[1];
		    for (i__ = jnint3_1.nhrf * jnint3_1.nhprf + 1; i__ <= 
			    i__2; ++i__) {
			do_fio(&c__1, (char *)&jnint1_1.w[jnindx_(&c__1, &i__,
				 &j) - 1], (ftnlen)sizeof(real));
			do_fio(&c__1, c__ + jnint1_1.nself[jnindx_(&c__1, &
				i__, &j) - 1], (ftnlen)1);
		    }
		    e_wsfe();
/* L220: */
		}
	    }
	}
	io___84.ciunit = *nf;
	s_wsle(&io___84);
	e_wsle();
	io___85.ciunit = *nf;
	s_wsfe(&io___85);
	do_fio(&c__1, (char *)&c__1, (ftnlen)sizeof(integer));
	e_wsfe();
	io___86.ciunit = *nf;
	s_wsle(&io___86);
	e_wsle();
	io___87.ciunit = *nf;
	s_wsfe(&io___87);
	i__1 = jnint2_1.m[1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&jnint1_1.t[jnindx_(&c__1, &i__, &c__0) - 1]
		    , (ftnlen)sizeof(real));
	    do_fio(&c__1, c__ + jnint1_1.ntself[jnindx_(&c__1, &i__, &c__0) - 
		    1], (ftnlen)1);
	}
	e_wsfe();
	i__1 = jnint2_1.nl;
	for (il = 2; il <= i__1; ++il) {
	    io___89.ciunit = *nf;
	    s_wsle(&io___89);
	    e_wsle();
	    io___90.ciunit = *nf;
	    s_wsfe(&io___90);
	    i__2 = il - 1;
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&il, (ftnlen)sizeof(integer));
	    e_wsfe();
	    i__2 = jnint2_1.m[il - 1];
	    for (j = 1; j <= i__2; ++j) {
		io___91.ciunit = *nf;
		s_wsle(&io___91);
		e_wsle();
		io___92.ciunit = *nf;
		s_wsfe(&io___92);
		i__3 = jnint2_1.m[il];
		for (i__ = 1; i__ <= i__3; ++i__) {
		    do_fio(&c__1, (char *)&jnint1_1.w[jnindx_(&il, &i__, &j) 
			    - 1], (ftnlen)sizeof(real));
		    do_fio(&c__1, c__ + jnint1_1.nself[jnindx_(&il, &i__, &j) 
			    - 1], (ftnlen)1);
		}
		e_wsfe();
/* L310: */
	    }
	    io___93.ciunit = *nf;
	    s_wsle(&io___93);
	    e_wsle();
	    io___94.ciunit = *nf;
	    s_wsfe(&io___94);
	    do_fio(&c__1, (char *)&il, (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___95.ciunit = *nf;
	    s_wsle(&io___95);
	    e_wsle();
	    io___96.ciunit = *nf;
	    s_wsfe(&io___96);
	    i__2 = jnint2_1.m[il];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		do_fio(&c__1, (char *)&jnint1_1.t[jnindx_(&il, &i__, &c__0) - 
			1], (ftnlen)sizeof(real));
		do_fio(&c__1, c__ + jnint1_1.ntself[jnindx_(&il, &i__, &c__0) 
			- 1], (ftnlen)1);
	    }
	    e_wsfe();
/* L300: */
	}
    }
    return 0;
/* **** END OF JNDUMP **************************************************** */
} /* jndump_ */

/* *********************************************************************** */
/* Subroutine */ int jnerr_(integer *ierr)
{
    /* Format strings */
    static char fmt_600[] = "(\002 *** JETNET ERROR:\002,i2,\002 ***\002)";
    static char fmt_610[] = "(\002 Illegal number of layers (\002,i3,\002"
	    ")\002)";
    static char fmt_620[] = "(\002 Total number of nodes (\002,i6,\002) exce"
	    "eds limit (\002,i6,\002).\002)";
    static char fmt_630[] = "(\002 Total number of weights (\002,i6,\002) ex"
	    "ceeds limit (\002,i6,\002).\002)";
    static char fmt_640[] = "(\002 Number of nodes in output layer is incomp"
	    "atible \002,/,\002 with the dimension of the Potts-nodes.\002)";
    static char fmt_650[] = "(\002 JETNET must be initialized (with JNINIT o"
	    "r JNREAD) \002,\002before \002,a6,\002 can be called.\002)";
    static char fmt_660[] = "(\002 Total number of input nodes (\002,i6,\002"
	    ") exceeds limit (\002,i6,\002).\002)";
    static char fmt_670[] = "(\002 Total number of output nodes (\002,i6,"
	    "\002) exceeds limit (\002,i6,\002).\002)";
    static char fmt_680[] = "(\002 Undefined updating algorithm (\002,i2,"
	    "\002) chosen.\002)";
    static char fmt_690[] = "(\002 Inconsistent geometry for receptive field"
	    "s:\002,/,\002 (MSTJN(23) = \002,i4,\002, MSTJN(24) = \002,i4,"
	    "\002, MSTJN(25) = \002,i4,\002, MSTJN(26) = \002,i4,\002)\002)";
    static char fmt_700[] = "(\002 Too few input nodes (=\002,i4,\002) for r"
	    "eceptive fields\002,/,\002 (MSTJN(23) = \002,i4,\002 and MSTJN(2"
	    "4) = \002,i4,\002).\002)";
    static char fmt_710[] = "(\002 In JNSEFI: attempt to connect/disconnect "
	    "unconnectable\002,\002 nodes.\002)";
    static char fmt_720[] = "(\002 Cannot read file - wrong format. Try JNRO"
	    "LD instead.\002)";
    static char fmt_730[] = "(\002 Chopping not allowed on non-sigmoid funct"
	    "ions.\002)";
    static char fmt_740[] = "(\002 Undefined transfer function (\002,i2,\002"
	    ") in GJN.\002)";
    static char fmt_750[] = "(\002 Call to JNINIT after calling JMINIT\002)";
    static char fmt_760[] = "(\002 JNREAD cannot read data-file produced by "
	    "JMDUMP\002)";
    static char fmt_770[] = "(\002 JNROLD cannot read data-file produced by "
	    "JMDUMP\002)";
    static char fmt_780[] = "(\002 You cannot start learning by terminating "
	    "Conj. Grad.\002)";
    static char fmt_790[] = "(\002 Too many warnings issued by JETNET.\002)";
    static char fmt_800[] = "(\002 Nr. of restarts (\002,i4,\002) in Quickpr"
	    "op, line search, or \002)";
    static char fmt_805[] = "(\002 Scaled Conj. Grad. exceeds maximum MSTJN("
	    "36) = \002,i4)";
    static char fmt_810[] = "(\002 MSTJN(9) (\002,i3,\002) must be > 0\002)";
    static char fmt_820[] = "(\002 Layer \002,i2,\002 has no nodes\002)";
    static char fmt_830[] = "(\002 Nr. of weights (\002,i6,\002) exceeds lim"
	    "it in JNHESS.\002)";
    static char fmt_840[] = "(\002 Nr. of calls to JNHESS (\002,i5,\002) mus"
	    "t be an integer \002)";
    static char fmt_850[] = "(\002 multiple of MSTJN(9)*MSTJN(2) (\002,i5"
	    ",\002) if JNHEIG\002,\002 is invoked\002)";
    static char fmt_860[] = "(\002 Too many iterations in subroutine JNTQLI"
	    ".\002)";
    static char fmt_870[] = "(\002 Error function, MSTJN(4) = \002,i2,\002, "
	    "incompatible with\002)";
    static char fmt_880[] = "(\002 using output transfer function = \002,i3)";
    static char fmt_890[] = "(\002 Updating turned off, MSTJN(5) = 9, when c"
	    "alling JNINIT\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___97 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___98 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___99 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___100 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___101 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___102 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___103 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___104 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___105 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___106 = { 0, 0, 0, fmt_680, 0 };
    static cilist io___107 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___109 = { 0, 0, 0, fmt_700, 0 };
    static cilist io___110 = { 0, 0, 0, fmt_710, 0 };
    static cilist io___111 = { 0, 0, 0, fmt_720, 0 };
    static cilist io___112 = { 0, 0, 0, fmt_730, 0 };
    static cilist io___113 = { 0, 0, 0, fmt_740, 0 };
    static cilist io___114 = { 0, 0, 0, fmt_750, 0 };
    static cilist io___115 = { 0, 0, 0, fmt_760, 0 };
    static cilist io___116 = { 0, 0, 0, fmt_770, 0 };
    static cilist io___117 = { 0, 0, 0, fmt_780, 0 };
    static cilist io___118 = { 0, 0, 0, fmt_790, 0 };
    static cilist io___119 = { 0, 0, 0, fmt_800, 0 };
    static cilist io___120 = { 0, 0, 0, fmt_805, 0 };
    static cilist io___121 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___122 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___123 = { 0, 0, 0, fmt_810, 0 };
    static cilist io___124 = { 0, 0, 0, fmt_820, 0 };
    static cilist io___125 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___126 = { 0, 0, 0, fmt_830, 0 };
    static cilist io___127 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___128 = { 0, 0, 0, fmt_840, 0 };
    static cilist io___129 = { 0, 0, 0, fmt_850, 0 };
    static cilist io___130 = { 0, 0, 0, fmt_860, 0 };
    static cilist io___131 = { 0, 0, 0, fmt_870, 0 };
    static cilist io___132 = { 0, 0, 0, fmt_880, 0 };
    static cilist io___133 = { 0, 0, 0, fmt_890, 0 };


/* ...JetNet subroutine ERROR */
/* ...Writes out an error message and stops the execution */
    if (jndat1_1.mstjm[7] == 1) {
	jndat1_1.mstjn[5] = jndat1_1.mstjm[5];
    }
    io___97.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___97);
    do_fio(&c__1, (char *)&(*ierr), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ierr == 1) {
	io___98.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___98);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[0], (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 2) {
	io___99.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___99);
	do_fio(&c__1, (char *)&jnint2_1.mv0[jnint2_1.nl], (ftnlen)sizeof(
		integer));
	do_fio(&c__1, (char *)&c__2000, (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 3) {
	io___100.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___100);
	do_fio(&c__1, (char *)&jnint2_1.mm0[jnint2_1.nl], (ftnlen)sizeof(
		integer));
	do_fio(&c__1, (char *)&c_b168, (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 4) {
	io___101.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___101);
	e_wsfe();
    } else if (*ierr == 5) {
	io___102.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___102);
	do_fio(&c__1, "JNINDX", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 6) {
	io___103.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___103);
	do_fio(&c__1, "JNSEFI", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 7) {
	io___104.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___104);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[9], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&c__1000, (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 8) {
	io___105.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___105);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[jnint2_1.nl - 1], (ftnlen)
		sizeof(integer));
	do_fio(&c__1, (char *)&c__1000, (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 9) {
	io___106.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___106);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[4], (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 10) {
	io___107.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___107);
	for (i__ = 23; i__ <= 26; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
    } else if (*ierr == 11) {
	io___109.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___109);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[9], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&jndat1_1.mstjn[22], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&jndat1_1.mstjn[23], (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 12) {
	io___110.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___110);
	e_wsfe();
    } else if (*ierr == 13) {
	io___111.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___111);
	e_wsfe();
    } else if (*ierr == 14) {
	io___112.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___112);
	e_wsfe();
    } else if (*ierr == 15) {
	io___113.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___113);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[2], (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 16) {
	io___114.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___114);
	e_wsfe();
    } else if (*ierr == 17) {
	io___115.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___115);
	e_wsfe();
    } else if (*ierr == 18) {
	io___116.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___116);
	e_wsfe();
    } else if (*ierr == 19) {
	io___117.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___117);
	e_wsfe();
    } else if (*ierr == 20) {
	io___118.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___118);
	e_wsfe();
    } else if (*ierr == 21) {
	io___119.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___119);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[37], (ftnlen)sizeof(integer));
	e_wsfe();
	io___120.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___120);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[35], (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 22) {
	io___121.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___121);
	do_fio(&c__1, "JNTRAL", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 23) {
	io___122.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___122);
	do_fio(&c__1, "JNTEST", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 24) {
	io___123.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___123);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[8], (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 25) {
	io___124.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___124);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[6], (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 26) {
	io___125.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___125);
	do_fio(&c__1, "JNHESS", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 27) {
	io___126.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___126);
	i__1 = jnint2_1.mm0[jnint2_1.nl] + jnint2_1.mv0[jnint2_1.nl];
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 28) {
	io___127.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___127);
	do_fio(&c__1, "JNHDIA", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 29) {
	io___128.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___128);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[38], (ftnlen)sizeof(integer));
	e_wsfe();
	io___129.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___129);
	i__1 = jndat1_1.mstjn[8] * jndat1_1.mstjn[1];
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 30) {
	io___130.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___130);
	e_wsfe();
    } else if (*ierr == 31) {
	io___131.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___131);
	do_fio(&c__1, (char *)&jndat1_1.mstjn[3], (ftnlen)sizeof(integer));
	e_wsfe();
	io___132.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___132);
	do_fio(&c__1, (char *)&jndat2_1.igfn[jnint2_1.nl - 1], (ftnlen)sizeof(
		integer));
	e_wsfe();
    } else if (*ierr == 32) {
	io___133.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___133);
	e_wsfe();
    }
    if (*ierr > 0) {
	s_stop("0", (ftnlen)1);
    }
    return 0;
/* **** END OF JNERR ***************************************************** */
} /* jnerr_ */

/* *********************************************************************** */
/* Subroutine */ int jnfeed_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1;

    /* Local variables */
    static real beta;
    static integer i__, j, ihprf;
    static real dd;
    static integer ih, il, in, mi, jo, ix, iy, jy, iw, jx;
    extern doublereal gjn_(integer *, real *, integer *);
    static integer mij, inx, iny;

/* ...JetNet subroutine FEED signal through net */
/* ...Takes the the values of OIN and calculates the values of */
/* ...the output nodes without writing to OUT */
/* ...set beta in first layer */
    if (jndat2_1.tinv[0] == 0.f) {
	beta = jndat1_1.parjn[2];
    } else {
	beta = dabs(jndat2_1.tinv[0]);
    }
/* ...calculate nodes in first layer */
    if (jnint3_1.nxin == 0) {
/* ...normal first layer */
	i__1 = jnint2_1.m[1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.a[i__ - 1] = jnint1_1.t[i__ - 1];
	    mij = i__ - jnint2_1.m[1];
	    i__2 = jnint2_1.m[0];
	    for (j = 1; j <= i__2; ++j) {
		mij += jnint2_1.m[1];
		jnint1_1.a[i__ - 1] += jnint1_1.w[mij - 1] * jndat1_1.oin[j - 
			1];
/* L110: */
	    }
	    r__1 = beta * jnint1_1.a[i__ - 1];
	    jnint1_1.o[i__ - 1] = gjn_(&i__, &r__1, jnint2_1.ng);
/* L100: */
	}
    } else {
/* ...receptive fields in first layer */
	i__1 = jnint3_1.nhprf;
	for (ihprf = 1; ihprf <= i__1; ++ihprf) {
	    i__2 = jnint3_1.nyhrf;
	    for (iy = 1; iy <= i__2; ++iy) {
		ih = iy - jnint3_1.nyhrf + (ihprf - 1) * jnint3_1.nhrf;
		i__3 = jnint3_1.nxhrf;
		for (ix = 1; ix <= i__3; ++ix) {
		    ih += jnint3_1.nyhrf;
		    jnint1_1.a[ih - 1] = jnint1_1.t[ih - 1];
		    i__4 = jnint3_1.nyrf;
		    for (jy = 1; jy <= i__4; ++jy) {
			iw = jy - jnint3_1.nyrf + (ihprf - 1) * jnint3_1.nrfw;
			i__5 = jnint3_1.nxrf;
			for (jx = 1; jx <= i__5; ++jx) {
			    iw += jnint3_1.nyrf;
			    inx = ix + jx - 1;
			    if (inx > abs(jnint3_1.nxin)) {
				inx -= abs(jnint3_1.nxin);
			    }
			    iny = iy + jy - 1;
			    if (iny > abs(jnint3_1.nyin)) {
				iny -= abs(jnint3_1.nyin);
			    }
			    in = (inx - 1) * abs(jnint3_1.nyin) + iny;
			    jnint1_1.a[ih - 1] += jnint1_1.w[iw - 1] * 
				    jndat1_1.oin[in - 1];
/* L160: */
			}
/* L150: */
		    }
		    i__5 = jnint2_1.m[0];
		    for (in = (i__4 = jnint3_1.nxin * jnint3_1.nyin, abs(i__4)
			    ) + 1; in <= i__5; ++in) {
			iw = jnint3_1.nxrf * jnint3_1.nyrf + in - (i__4 = 
				jnint3_1.nxin * jnint3_1.nyin, abs(i__4)) + (
				ihprf - 1) * jnint3_1.nrfw;
			jnint1_1.a[ih - 1] += jnint1_1.w[iw - 1] * 
				jndat1_1.oin[in - 1];
/* L170: */
		    }
		    r__1 = beta * jnint1_1.a[ih - 1];
		    jnint1_1.o[ih - 1] = gjn_(&ih, &r__1, jnint2_1.ng);
/* L140: */
		}
/* L130: */
	    }
/* L120: */
	}
	i__1 = jnint2_1.m[1];
	for (ih = jnint3_1.nhrf * jnint3_1.nhprf + 1; ih <= i__1; ++ih) {
	    jnint1_1.a[ih - 1] = jnint1_1.t[ih - 1];
	    iw = jnint3_1.nhrf * jnint3_1.nhprf + ih - jnint2_1.m[1];
	    i__2 = jnint2_1.m[0];
	    for (in = 1; in <= i__2; ++in) {
		iw += jnint2_1.m[1];
		jnint1_1.a[ih - 1] += jnint1_1.w[iw - 1] * jndat1_1.oin[in - 
			1];
/* L190: */
	    }
	    r__1 = beta * jnint1_1.a[ih - 1];
	    jnint1_1.o[ih - 1] = gjn_(&ih, &r__1, jnint2_1.ng);
/* L180: */
	}
    }
/* ...calculate nodes in following layers */
    i__1 = jnint2_1.nl;
    for (il = 2; il <= i__1; ++il) {
/* ...set beta in layer IL */
	if (jndat2_1.tinv[il - 1] == 0.f) {
	    beta = jndat1_1.parjn[2];
	} else {
	    beta = (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	}
/* ...calculate nodes in layer IL */
	i__2 = jnint2_1.m[il];
	for (i__ = 1; i__ <= i__2; ++i__) {
	    mi = jnint2_1.mv0[il - 1] + i__;
	    jnint1_1.a[mi - 1] = jnint1_1.t[mi - 1];
	    mij = jnint2_1.mm0[il - 1] - jnint2_1.m[il] + i__;
	    i__3 = jnint2_1.mv0[il - 2] + jnint2_1.m[il - 1];
	    for (j = jnint2_1.mv0[il - 2] + 1; j <= i__3; ++j) {
		mij += jnint2_1.m[il];
		jnint1_1.a[mi - 1] += jnint1_1.w[mij - 1] * jnint1_1.o[j - 1];
/* L220: */
	    }
	    r__1 = beta * jnint1_1.a[mi - 1];
	    jnint1_1.o[mi - 1] = gjn_(&mi, &r__1, &jnint2_1.ng[il - 1]);
/* L210: */
	}
/* L200: */
    }
    if (jnint2_1.ipott < 2) {
	return 0;
    }
/* ...Special treatment of output layer if Potts-nodes */
    i__1 = jnint2_1.m[jnint2_1.nl] / jnint2_1.ipott;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dd = 0.f;
	jo = jnint2_1.mv0[jnint2_1.nl - 1] + (i__ - 1) * jnint2_1.ipott;
	i__2 = jnint2_1.ipott;
	for (j = 1; j <= i__2; ++j) {
	    dd += jnint1_1.o[jo + j - 1];
/* L310: */
	}
	i__2 = jnint2_1.ipott;
	for (j = 1; j <= i__2; ++j) {
	    jnint1_1.o[jo + j - 1] /= dd;
/* L320: */
	}
/* L300: */
    }
    return 0;
/* **** END OF JNFEED **************************************************** */
} /* jnfeed_ */

/* *********************************************************************** */
/* Subroutine */ int jnhead_(void)
{
    /* Format strings */
    static char fmt_600[] = "(14x,\002The Lund Neural Network Program - JETN"
	    "ET version 3.4\002)";
    static char fmt_610[] = "(14x,\002******   Latest date of change: July 6"
	    ", 1995  ******\002)";

    /* Builtin functions */
    integer s_wsle(cilist *), e_wsle(void), s_wsfe(cilist *), e_wsfe(void);

    /* Fortran I/O blocks */
    static cilist io___152 = { 0, 0, 0, 0, 0 };
    static cilist io___153 = { 0, 0, 0, 0, 0 };
    static cilist io___154 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___155 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___156 = { 0, 0, 0, 0, 0 };


/* ...JetNet subroutine write HEADer */
/* ...Writes a header on file number NF */
    if (jndat1_1.mstjm[7] == 1) {
	jndat1_1.mstjn[5] = jndat1_1.mstjm[5];
    }
    io___152.ciunit = jndat1_1.mstjn[5];
    s_wsle(&io___152);
    e_wsle();
    io___153.ciunit = jndat1_1.mstjn[5];
    s_wsle(&io___153);
    e_wsle();
    io___154.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___154);
    e_wsfe();
    io___155.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___155);
    e_wsfe();
    io___156.ciunit = jndat1_1.mstjn[5];
    s_wsle(&io___156);
    e_wsle();
    return 0;
/* **** END OF JNHEAD **************************************************** */
} /* jnhead_ */

/* *********************************************************************** */
/* Subroutine */ int jnheig_(integer *igrad)
{
    extern /* Subroutine */ int jnerr_(integer *);
    static integer nwgts;
    extern /* Subroutine */ int jntred_(integer *, integer *), jntqli_(
	    integer *, integer *);

/* ...JetNet subroutine Hessian EIGenvalues. */
/* ...Diagonalizes the Hessian matrix stored in D2E. The eigenvalues */
/* ...are placed in the vector OUT. If IGRAD=-1 then the */
/* ...eigenvectors of the Hessian are calculated. */
    if (jndat1_1.mstjn[7] == 0) {
	jnerr_(&c__28);
    }
    if (jndat1_1.mstjn[38] % (jndat1_1.mstjn[1] * jndat1_1.mstjn[8]) != 0 || 
	    jndat1_1.mstjn[38] <= 0) {
	jnerr_(&c__29);
    }
    nwgts = jnint2_1.mm0[jnint2_1.nl] + jnint2_1.mv0[jnint2_1.nl];
    jntred_(&nwgts, igrad);
    jntqli_(&nwgts, igrad);
    return 0;
/* **** END OF JNHEIG **************************************************** */
} /* jnheig_ */

/* *********************************************************************** */
/* Subroutine */ int jnhess_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    real r__1;

    /* Local variables */
    static real beta;
    static integer indx;
    static real term, fact2;
    static integer i__, j, k, l;
    static real q[90000]	/* was [300][300] */;
    extern /* Subroutine */ int jnerr_(integer *);
    static integer iofst, j2, j1, k1, k2, m1, nwgts, iofst1, iofst2;
    static real dd[300];
    static integer ii, il, mi, mj;
    extern /* Subroutine */ int jnfeed_(void);
    static real gamcap;
    static integer mk, ml, mm, iw, iv;
    static real factor;
    extern integer jnindx_(integer *, integer *, integer *);
    static integer il1, il2, mj2, mj1, mk1, mk2, il3;
    static real gam;
    static integer mij, mjk, mkl, mlm;
    static real sum;
    static integer mij1, mij2, mjk1, mjk2;

/* ...JetNet subroutine calculate HESSian */
/* ...Calculates the Hessian for the network. It assumes a summed square */
/* ...error (MSTJN(4)=0). */
    if (jndat1_1.mstjn[7] == 0) {
	jnerr_(&c__26);
    }
    if (jndat1_1.mstjn[8] <= 0) {
	jnerr_(&c__24);
    }
    nwgts = jnint2_1.mm0[jnint2_1.nl] + jnint2_1.mv0[jnint2_1.nl];
    if (nwgts > 300) {
	jnerr_(&c__27);
    }
    if (jndat1_1.mstjn[38] % (jndat1_1.mstjn[1] * jndat1_1.mstjn[8]) == 0) {
/* ...zero Hessian: */
	i__1 = nwgts;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = nwgts;
	    for (j = 1; j <= i__2; ++j) {
		jnint5_1.d2e[i__ + j * 300 - 301] = 0.f;
/* L110: */
	    }
/* L100: */
	}
    }
    ++jndat1_1.mstjn[38];
    jnfeed_();
/* ...rescale GPJN and GPPJN: */
    i__1 = jnint2_1.nl;
    for (il = 1; il <= i__1; ++il) {
	if (jndat2_1.tinv[il - 1] == 0.f) {
	    beta = jndat1_1.parjn[2];
	} else {
	    beta = (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	}
	i__2 = jnint2_1.m[il];
	for (i__ = 1; i__ <= i__2; ++i__) {
	    mi = jnindx_(&il, &i__, &c__0);
	    jnsigm_1.gpjn[mi - 1] *= beta;
/* Computing 2nd power */
	    r__1 = beta;
	    jnsigm_1.gppjn[mi - 1] *= r__1 * r__1;
/* L210: */
	}
/* L200: */
    }
/* ...compute Q: */
    for (il1 = jnint2_1.nl; il1 >= 2; --il1) {
	for (il2 = il1 - 1; il2 >= 1; --il2) {
	    i__1 = jnint2_1.m[il1];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		mi = jnindx_(&il1, &i__, &c__0);
		i__2 = jnint2_1.m[il2];
		for (j = 1; j <= i__2; ++j) {
		    mj = jnindx_(&il2, &j, &c__0);
		    if (il2 == il1 - 1) {
			mij = jnindx_(&il1, &i__, &j);
			q[mi + mj * 300 - 301] = jnint1_1.w[mij - 1];
		    } else {
			sum = 0.f;
			i__3 = jnint2_1.m[il2 + 1];
			for (j2 = 1; j2 <= i__3; ++j2) {
			    i__4 = il2 + 1;
			    mj2 = jnindx_(&i__4, &j2, &c__0);
			    i__4 = il2 + 1;
			    mjk = jnindx_(&i__4, &j2, &j);
			    sum += q[mi + mj2 * 300 - 301] * jnsigm_1.gpjn[
				    mj2 - 1] * jnint1_1.w[mjk - 1];
/* L260: */
			}
			q[mi + mj * 300 - 301] = sum;
		    }
/* L250: */
		}
/* L240: */
	    }
/* L230: */
	}
/* L220: */
    }
/* ...Loop over output units: */
    i__1 = jnint2_1.m[jnint2_1.nl];
    for (i__ = 1; i__ <= i__1; ++i__) {
	mi = jnindx_(&jnint2_1.nl, &i__, &c__0);
	gam = (jnint1_1.o[mi - 1] - jndat1_1.out[i__ - 1]) * jnsigm_1.gpjn[mi 
		- 1];
	gamcap = (jnint1_1.o[mi - 1] - jndat1_1.out[i__ - 1]) * 
		jnsigm_1.gppjn[mi - 1];
/* ...Diagonal weights - output layer: */
	jnint5_1.d2e[i__ + i__ * 300 - 301] += gamcap;
	if (jnint2_1.nl != 1) {
	    i__2 = jnint2_1.m[jnint2_1.nl - 1];
	    for (j1 = 1; j1 <= i__2; ++j1) {
		i__3 = jnint2_1.nl - 1;
		mj1 = jnindx_(&i__3, &j1, &c__0);
		mij1 = jnint2_1.m[jnint2_1.nl] + (i__ - 1) * jnint2_1.m[
			jnint2_1.nl - 1] + j1;
		term = gamcap * jnint1_1.o[mj1 - 1];
		jnint5_1.d2e[i__ + mij1 * 300 - 301] += term;
		i__3 = j1;
		for (j2 = 1; j2 <= i__3; ++j2) {
		    i__4 = jnint2_1.nl - 1;
		    mj2 = jnindx_(&i__4, &j2, &c__0);
		    mij2 = jnint2_1.m[jnint2_1.nl] + (i__ - 1) * jnint2_1.m[
			    jnint2_1.nl - 1] + j2;
		    jnint5_1.d2e[mij2 + mij1 * 300 - 301] += term * 
			    jnint1_1.o[mj2 - 1];
/* L320: */
		}
/* L310: */
	    }
	} else {
	    i__2 = jnint2_1.m[0];
	    for (j1 = 1; j1 <= i__2; ++j1) {
		mij1 = jnint2_1.m[jnint2_1.nl] + (i__ - 1) * jnint2_1.m[0] + 
			j1;
		term = gamcap * jndat1_1.oin[j1 - 1];
		jnint5_1.d2e[i__ + mij1 * 300 - 301] += term;
		i__3 = j1;
		for (j2 = 1; j2 <= i__3; ++j2) {
		    mij2 = jnint2_1.m[jnint2_1.nl] + (i__ - 1) * jnint2_1.m[0]
			     + j2;
		    jnint5_1.d2e[mij2 + mij1 * 300 - 301] += term * 
			    jndat1_1.oin[j2 - 1];
/* L340: */
		}
/* L330: */
	    }
	}
/* ...Diagonal weights - other layers: */
	iofst = 0;
	for (il = jnint2_1.nl - 1; il >= 1; --il) {
	    iofst += jnint2_1.m[il + 1] * (jnint2_1.m[il] + 1);
	    i__2 = jnint2_1.m[il];
	    for (j1 = 1; j1 <= i__2; ++j1) {
		mj1 = jnindx_(&il, &j1, &c__0);
		factor = gam * q[mi + mj1 * 300 - 301] * jnsigm_1.gppjn[mj1 - 
			1];
		jnint5_1.d2e[iofst + j1 + (iofst + j1) * 300 - 301] += factor;
		if (il > 1) {
		    i__3 = jnint2_1.m[il - 1];
		    for (k1 = 1; k1 <= i__3; ++k1) {
			i__4 = il - 1;
			mk1 = jnindx_(&i__4, &k1, &c__0);
			mjk1 = jnint2_1.m[il] + (j1 - 1) * jnint2_1.m[il - 1] 
				+ k1;
			term = factor * jnint1_1.o[mk1 - 1];
			jnint5_1.d2e[iofst + j1 + (iofst + mjk1) * 300 - 301] 
				+= term;
			i__4 = k1;
			for (k2 = 1; k2 <= i__4; ++k2) {
			    i__5 = il - 1;
			    mk2 = jnindx_(&i__5, &k2, &c__0);
			    mjk2 = jnint2_1.m[il] + (j1 - 1) * jnint2_1.m[il 
				    - 1] + k2;
			    jnint5_1.d2e[iofst + mjk2 + (iofst + mjk1) * 300 
				    - 301] += term * jnint1_1.o[mk2 - 1];
/* L380: */
			}
/* L370: */
		    }
		} else {
		    i__3 = jnint2_1.m[0];
		    for (k1 = 1; k1 <= i__3; ++k1) {
			mjk1 = jnint2_1.m[1] + (j1 - 1) * jnint2_1.m[0] + k1;
			term = factor * jndat1_1.oin[k1 - 1];
			jnint5_1.d2e[iofst + j1 + (iofst + mjk1) * 300 - 301] 
				+= term;
			i__4 = k1;
			for (k2 = 1; k2 <= i__4; ++k2) {
			    mjk2 = jnint2_1.m[1] + (j1 - 1) * jnint2_1.m[0] + 
				    k2;
			    jnint5_1.d2e[iofst + mjk2 + (iofst + mjk1) * 300 
				    - 301] += term * jndat1_1.oin[k2 - 1];
/* L400: */
			}
/* L390: */
		    }
		}
		i__3 = jnint2_1.m[il];
		for (j2 = 1; j2 <= i__3; ++j2) {
		    mj2 = jnindx_(&il, &j2, &c__0);
		    factor = gamcap * q[mi + mj1 * 300 - 301] * jnsigm_1.gpjn[
			    mj1 - 1] * q[mi + mj2 * 300 - 301] * 
			    jnsigm_1.gpjn[mj2 - 1];
		    if (il <= jnint2_1.nl - 2) {
			sum = 0.f;
			i__4 = il + 1;
			for (il2 = jnint2_1.nl - 1; il2 >= i__4; --il2) {
			    i__5 = jnint2_1.m[il2];
			    for (j = 1; j <= i__5; ++j) {
				mj = jnindx_(&il2, &j, &c__0);
				sum += q[mi + mj * 300 - 301] * 
					jnsigm_1.gppjn[mj - 1] * q[mj + mj1 * 
					300 - 301] * q[mj + mj2 * 300 - 301];
/* L430: */
			    }
/* L420: */
			}
			factor += sum * gam * jnsigm_1.gpjn[mj1 - 1] * 
				jnsigm_1.gpjn[mj2 - 1];
		    }
		    if (j2 >= j1) {
			jnint5_1.d2e[iofst + j1 + (iofst + j2) * 300 - 301] +=
				 factor;
		    }
		    if (il > 1) {
			i__4 = jnint2_1.m[il - 1];
			for (k1 = 1; k1 <= i__4; ++k1) {
			    i__5 = il - 1;
			    mk1 = jnindx_(&i__5, &k1, &c__0);
			    mjk1 = jnint2_1.m[il] + (j2 - 1) * jnint2_1.m[il 
				    - 1] + k1;
			    term = factor * jnint1_1.o[mk1 - 1];
			    jnint5_1.d2e[iofst + j1 + (iofst + mjk1) * 300 - 
				    301] += term;
			    if (j2 == j1) {
				i__5 = k1;
				for (k2 = 1; k2 <= i__5; ++k2) {
				    i__6 = il - 1;
				    mk2 = jnindx_(&i__6, &k2, &c__0);
				    mjk2 = jnint2_1.m[il] + (j1 - 1) * 
					    jnint2_1.m[il - 1] + k2;
				    jnint5_1.d2e[iofst + mjk2 + (iofst + mjk1)
					     * 300 - 301] += term * 
					    jnint1_1.o[mk2 - 1];
/* L440: */
				}
			    } else if (j2 > j1) {
				i__5 = jnint2_1.m[il - 1];
				for (k2 = 1; k2 <= i__5; ++k2) {
				    i__6 = il - 1;
				    mk2 = jnindx_(&i__6, &k2, &c__0);
				    mjk2 = jnint2_1.m[il] + (j1 - 1) * 
					    jnint2_1.m[il - 1] + k2;
				    jnint5_1.d2e[iofst + mjk2 + (iofst + mjk1)
					     * 300 - 301] += term * 
					    jnint1_1.o[mk2 - 1];
/* L450: */
				}
			    }
/* L431: */
			}
		    } else {
			i__4 = jnint2_1.m[0];
			for (k1 = 1; k1 <= i__4; ++k1) {
			    mjk1 = jnint2_1.m[1] + (j2 - 1) * jnint2_1.m[0] + 
				    k1;
			    term = factor * jndat1_1.oin[k1 - 1];
			    jnint5_1.d2e[iofst + j1 + (iofst + mjk1) * 300 - 
				    301] += term;
			    if (j2 == j1) {
				i__5 = k1;
				for (k2 = 1; k2 <= i__5; ++k2) {
				    mjk2 = jnint2_1.m[1] + (j1 - 1) * 
					    jnint2_1.m[0] + k2;
				    jnint5_1.d2e[iofst + mjk2 + (iofst + mjk1)
					     * 300 - 301] += term * 
					    jndat1_1.oin[k2 - 1];
/* L470: */
				}
			    } else if (j2 > j1) {
				i__5 = jnint2_1.m[0];
				for (k2 = 1; k2 <= i__5; ++k2) {
				    mjk2 = jnint2_1.m[1] + (j1 - 1) * 
					    jnint2_1.m[0] + k2;
				    jnint5_1.d2e[iofst + mjk2 + (iofst + mjk1)
					     * 300 - 301] += term * 
					    jndat1_1.oin[k2 - 1];
/* L480: */
				}
			    }
/* L460: */
			}
		    }
/* L410: */
		}
/* L360: */
	    }
/* L350: */
	}
/* ...End of diagonal weights. */
/* ...1st off-diagonal - output layer: */
	iofst2 = jnint2_1.m[jnint2_1.nl] * (jnint2_1.m[jnint2_1.nl - 1] + 1);
	i__2 = jnint2_1.m[jnint2_1.nl - 1];
	for (j = 1; j <= i__2; ++j) {
	    i__3 = jnint2_1.nl - 1;
	    mj = jnindx_(&i__3, &j, &c__0);
	    factor = gamcap * q[mi + mj * 300 - 301] * jnsigm_1.gpjn[mj - 1];
	    jnint5_1.d2e[i__ + (iofst2 + j) * 300 - 301] += factor;
	    mij = jnint2_1.m[jnint2_1.nl] + (i__ - 1) * jnint2_1.m[
		    jnint2_1.nl - 1] + j;
	    fact2 = gam * jnsigm_1.gpjn[mj - 1];
	    jnint5_1.d2e[mij + (iofst2 + j) * 300 - 301] += fact2;
	    if (jnint2_1.nl > 2) {
		i__3 = jnint2_1.m[jnint2_1.nl - 2];
		for (k = 1; k <= i__3; ++k) {
		    i__4 = jnint2_1.nl - 2;
		    mk = jnindx_(&i__4, &k, &c__0);
		    mjk = jnint2_1.m[jnint2_1.nl - 1] + (j - 1) * jnint2_1.m[
			    jnint2_1.nl - 2] + k;
		    jnint5_1.d2e[i__ + (iofst2 + mjk) * 300 - 301] += factor *
			     jnint1_1.o[mk - 1];
		    jnint5_1.d2e[mij + (iofst2 + mjk) * 300 - 301] += fact2 * 
			    jnint1_1.o[mk - 1];
/* L510: */
		}
	    } else {
		i__3 = jnint2_1.m[0];
		for (k = 1; k <= i__3; ++k) {
		    mjk = jnint2_1.m[1] + (j - 1) * jnint2_1.m[0] + k;
		    jnint5_1.d2e[i__ + (iofst2 + mjk) * 300 - 301] += factor *
			     jndat1_1.oin[k - 1];
		    jnint5_1.d2e[mij + (iofst2 + mjk) * 300 - 301] += fact2 * 
			    jndat1_1.oin[k - 1];
/* L520: */
		}
	    }
	    i__3 = jnint2_1.m[jnint2_1.nl - 1];
	    for (j2 = 1; j2 <= i__3; ++j2) {
		i__4 = jnint2_1.nl - 1;
		mj2 = jnindx_(&i__4, &j2, &c__0);
		mij2 = jnint2_1.m[jnint2_1.nl] + (i__ - 1) * jnint2_1.m[
			jnint2_1.nl - 1] + j2;
		term = factor * jnint1_1.o[mj2 - 1];
		jnint5_1.d2e[mij2 + (iofst2 + j) * 300 - 301] += term;
		if (jnint2_1.nl > 2) {
		    i__4 = jnint2_1.m[jnint2_1.nl - 2];
		    for (k = 1; k <= i__4; ++k) {
			i__5 = jnint2_1.nl - 2;
			mk = jnindx_(&i__5, &k, &c__0);
			mjk = jnint2_1.m[jnint2_1.nl - 1] + (j - 1) * 
				jnint2_1.m[jnint2_1.nl - 2] + k;
			jnint5_1.d2e[mij2 + (iofst2 + mjk) * 300 - 301] += 
				term * jnint1_1.o[mk - 1];
/* L540: */
		    }
		} else {
		    i__4 = jnint2_1.m[0];
		    for (k = 1; k <= i__4; ++k) {
			mjk = jnint2_1.m[1] + (j - 1) * jnint2_1.m[0] + k;
			jnint5_1.d2e[mij2 + (iofst2 + mjk) * 300 - 301] += 
				term * jndat1_1.oin[k - 1];
/* L550: */
		    }
		}
/* L530: */
	    }
/* L500: */
	}
/* ...1st off-diagonal - other layers: */
	iofst1 = 0;
	for (il = jnint2_1.nl - 1; il >= 2; --il) {
	    iofst1 += jnint2_1.m[il + 1] * (jnint2_1.m[il] + 1);
	    iofst2 += jnint2_1.m[il] * (jnint2_1.m[il - 1] + 1);
	    i__2 = jnint2_1.m[il];
	    for (j = 1; j <= i__2; ++j) {
		mj = jnindx_(&il, &j, &c__0);
		i__3 = jnint2_1.m[il - 1];
		for (k = 1; k <= i__3; ++k) {
		    i__4 = il - 1;
		    mk = jnindx_(&i__4, &k, &c__0);
		    factor = gamcap * q[mi + mj * 300 - 301] * jnsigm_1.gpjn[
			    mj - 1] * q[mi + mk * 300 - 301] * jnsigm_1.gpjn[
			    mk - 1] + gam * q[mi + mj * 300 - 301] * 
			    jnsigm_1.gppjn[mj - 1] * q[mj + mk * 300 - 301] * 
			    jnsigm_1.gpjn[mk - 1];
		    if (il <= jnint2_1.nl - 2) {
			sum = 0.f;
			i__4 = il + 1;
			for (il2 = jnint2_1.nl - 1; il2 >= i__4; --il2) {
			    i__5 = jnint2_1.m[il2];
			    for (j2 = 1; j2 <= i__5; ++j2) {
				mj2 = jnindx_(&il2, &j2, &c__0);
				sum += q[mi + mj2 * 300 - 301] * 
					jnsigm_1.gppjn[mj2 - 1] * q[mj2 + ml *
					 300 - 301] * q[mj2 + mj * 300 - 301];
/* L600: */
			    }
/* L590: */
			}
			factor += sum * gam * jnsigm_1.gpjn[ml - 1] * 
				jnsigm_1.gpjn[mj - 1];
		    }
		    jnint5_1.d2e[iofst1 + i__ + (iofst2 + j) * 300 - 301] += 
			    factor;
		    mjk = jnint2_1.m[il] + (j - 1) * jnint2_1.m[il - 1] + k;
		    fact2 = gam * q[mi + mj * 300 - 301] * jnsigm_1.gpjn[mj - 
			    1] * jnsigm_1.gpjn[mk - 1];
		    jnint5_1.d2e[iofst1 + mjk + (iofst2 + j) * 300 - 301] += 
			    fact2;
		    if (il - 1 > 1) {
			i__4 = jnint2_1.m[il - 2];
			for (l = 1; l <= i__4; ++l) {
			    i__5 = il - 2;
			    ml = jnindx_(&i__5, &l, &c__0);
			    mkl = jnint2_1.m[il - 1] + (k - 1) * jnint2_1.m[
				    il - 2] + l;
			    jnint5_1.d2e[iofst1 + j + (iofst2 + mkl) * 300 - 
				    301] += factor * jnint1_1.o[ml - 1];
			    jnint5_1.d2e[iofst1 + mjk + (iofst2 + mkl) * 300 
				    - 301] += fact2 * jnint1_1.o[ml - 1];
/* L610: */
			}
		    } else {
			i__4 = jnint2_1.m[0];
			for (l = 1; l <= i__4; ++l) {
			    mkl = jnint2_1.m[1] + (k - 1) * jnint2_1.m[0] + l;
			    jnint5_1.d2e[iofst1 + j + (iofst2 + mkl) * 300 - 
				    301] += factor * jndat1_1.oin[l - 1];
			    jnint5_1.d2e[iofst1 + mjk + (iofst2 + mkl) * 300 
				    - 301] += fact2 * jndat1_1.oin[l - 1];
/* L620: */
			}
		    }
		    i__4 = jnint2_1.m[il - 1];
		    for (k2 = 1; k2 <= i__4; ++k2) {
			i__5 = il - 1;
			mk2 = jnindx_(&i__5, &k2, &c__0);
			mjk2 = jnint2_1.m[il] + (k - 1) * jnint2_1.m[il - 1] 
				+ k2;
			term = factor * jnint1_1.o[mk2 - 1];
			jnint5_1.d2e[iofst1 + mjk2 + (iofst2 + k) * 300 - 301]
				 += term;
			if (il - 1 > 1) {
			    i__5 = jnint2_1.m[il - 2];
			    for (l = 1; l <= i__5; ++l) {
				i__6 = il - 2;
				ml = jnindx_(&i__6, &l, &c__0);
				mkl = jnint2_1.m[il - 1] + (k - 1) * 
					jnint2_1.m[il - 2] + l;
				jnint5_1.d2e[iofst1 + mjk2 + (iofst2 + mkl) * 
					300 - 301] += term * jnint1_1.o[ml - 
					1];
/* L640: */
			    }
			} else {
			    i__5 = jnint2_1.m[0];
			    for (l = 1; l <= i__5; ++l) {
				mkl = jnint2_1.m[1] + (k - 1) * jnint2_1.m[0] 
					+ l;
				jnint5_1.d2e[iofst1 + mjk2 + (iofst2 + mkl) * 
					300 - 301] += term * jndat1_1.oin[l - 
					1];
/* L650: */
			    }
			}
/* L630: */
		    }
/* L580: */
		}
/* L570: */
	    }
/* L560: */
	}
/* ...End of 1st off-diagonal. */
/* ...Higher off-diagonals - output layer: */
	iofst2 = jnint2_1.m[jnint2_1.nl] * (jnint2_1.m[jnint2_1.nl - 1] + 1);
	for (il = jnint2_1.nl - 2; il >= 1; --il) {
	    iofst2 += jnint2_1.m[il + 1] * (jnint2_1.m[il] + 1);
	    i__2 = jnint2_1.m[il];
	    for (k = 1; k <= i__2; ++k) {
		mk = jnindx_(&il, &k, &c__0);
		factor = gamcap * q[mi + mk * 300 - 301] * jnsigm_1.gpjn[mk - 
			1];
		jnint5_1.d2e[i__ + (iofst2 + k) * 300 - 301] += factor;
		i__3 = jnint2_1.m[jnint2_1.nl - 1];
		for (j = 1; j <= i__3; ++j) {
		    i__4 = jnint2_1.nl - 1;
		    mj = jnindx_(&i__4, &j, &c__0);
		    mjk = jnint2_1.m[jnint2_1.nl] + (i__ - 1) * jnint2_1.m[
			    jnint2_1.nl - 1] + j;
		    fact2 = gam * jnsigm_1.gpjn[mj - 1] * q[mj + mk * 300 - 
			    301] * jnsigm_1.gpjn[mk - 1];
		    jnint5_1.d2e[mjk + (iofst2 + k) * 300 - 301] += fact2;
		    if (il > 1) {
			i__4 = jnint2_1.m[il - 1];
			for (l = 1; l <= i__4; ++l) {
			    i__5 = il - 1;
			    ml = jnindx_(&i__5, &l, &c__0);
			    mkl = jnint2_1.m[il] + (k - 1) * jnint2_1.m[il - 
				    1] + l;
			    jnint5_1.d2e[i__ + (iofst2 + mkl) * 300 - 301] += 
				    factor * jnint1_1.o[ml - 1];
			    jnint5_1.d2e[mjk + (iofst2 + mkl) * 300 - 301] += 
				    fact2 * jnint1_1.o[ml - 1];
/* L700: */
			}
		    } else {
			i__4 = jnint2_1.m[0];
			for (l = 1; l <= i__4; ++l) {
			    mkl = jnint2_1.m[1] + (k - 1) * jnint2_1.m[0] + l;
			    jnint5_1.d2e[i__ + (iofst2 + mkl) * 300 - 301] += 
				    factor * jndat1_1.oin[l - 1];
			    jnint5_1.d2e[mjk + (iofst2 + mkl) * 300 - 301] += 
				    fact2 * jndat1_1.oin[l - 1];
/* L710: */
			}
		    }
/* L690: */
		}
/* L670: */
	    }
/* L660: */
	}
/* ...Higher off-diagonals - other layers: */
	iofst1 = 0;
	for (il1 = jnint2_1.nl - 1; il1 >= 2; --il1) {
	    iofst1 += jnint2_1.m[il1 + 1] * (jnint2_1.m[il1] + 1);
	    iofst2 = jnint2_1.m[jnint2_1.nl] * (jnint2_1.m[jnint2_1.nl - 1] + 
		    1);
	    for (il2 = il1 - 2; il2 >= 1; --il2) {
		iofst2 += jnint2_1.m[il2 + 1] * (jnint2_1.m[il2] + 1);
		i__2 = jnint2_1.m[il1];
		for (j = 1; j <= i__2; ++j) {
		    mj = jnindx_(&il1, &j, &c__0);
		    i__3 = jnint2_1.m[il2];
		    for (l = 1; l <= i__3; ++l) {
			ml = jnindx_(&il2, &l, &c__0);
			factor = gamcap * q[mi + ml * 300 - 301] * 
				jnsigm_1.gpjn[ml - 1] * q[mi + mj * 300 - 301]
				 * jnsigm_1.gpjn[mj - 1] + gam * q[mi + mj * 
				300 - 301] * jnsigm_1.gppjn[mj - 1] * q[mj + 
				ml * 300 - 301] * jnsigm_1.gpjn[ml - 1];
			if (il1 <= jnint2_1.nl - 2) {
			    sum = 0.f;
			    i__4 = il1 + 1;
			    for (il3 = jnint2_1.nl - 1; il3 >= i__4; --il3) {
				i__5 = jnint2_1.m[il3];
				for (j2 = 1; j2 <= i__5; ++j2) {
				    mj2 = jnindx_(&il3, &j2, &c__0);
				    sum += q[mi + mj2 * 300 - 301] * 
					    jnsigm_1.gppjn[mj2 - 1] * q[mj2 + 
					    mj * 300 - 301] * q[mj2 + ml * 
					    300 - 301];
/* L770: */
				}
/* L760: */
			    }
			    factor += sum * gam * jnsigm_1.gpjn[ml - 1] * 
				    jnsigm_1.gpjn[mj - 1];
			}
			jnint5_1.d2e[iofst1 + j + (iofst2 + l) * 300 - 301] +=
				 factor;
			i__4 = jnint2_1.m[il1 - 1];
			for (k = 1; k <= i__4; ++k) {
			    i__5 = il1 - 1;
			    mk = jnindx_(&i__5, &k, &c__0);
			    mkl = jnint2_1.m[il1] + (j - 1) * jnint2_1.m[il1 
				    - 1] + k;
			    fact2 = gam * q[mi + mj * 300 - 301] * 
				    jnsigm_1.gpjn[mj - 1] * jnsigm_1.gpjn[mk 
				    - 1] * q[mk + ml * 300 - 301] * 
				    jnsigm_1.gpjn[ml - 1];
			    jnint5_1.d2e[iofst1 + mkl + (iofst2 + l) * 300 - 
				    301] += fact2;
			    if (il2 > 1) {
				i__5 = jnint2_1.m[il2 - 1];
				for (m1 = 1; m1 <= i__5; ++m1) {
				    i__6 = il2 - 1;
				    mm = jnindx_(&i__6, &m1, &c__0);
				    mlm = jnint2_1.m[il2] + (l - 1) * 
					    jnint2_1.m[il2 - 1] + m1;
				    jnint5_1.d2e[iofst1 + j + (iofst2 + mlm) *
					     300 - 301] += factor * 
					    jnint1_1.o[mm - 1];
				    jnint5_1.d2e[iofst1 + mkl + (iofst2 + mlm)
					     * 300 - 301] += fact2 * 
					    jnint1_1.o[mm - 1];
/* L790: */
				}
			    } else {
				i__5 = jnint2_1.m[0];
				for (m1 = 1; m1 <= i__5; ++m1) {
				    mlm = jnint2_1.m[1] + (l - 1) * 
					    jnint2_1.m[0] + m1;
				    jnint5_1.d2e[iofst1 + j + (iofst2 + mlm) *
					     300 - 301] += factor * 
					    jndat1_1.oin[m1 - 1];
				    jnint5_1.d2e[iofst1 + mkl + (iofst2 + mlm)
					     * 300 - 301] += fact2 * 
					    jndat1_1.oin[m1 - 1];
/* L800: */
				}
			    }
/* L780: */
			}
/* L750: */
		    }
/* L740: */
		}
/* L730: */
	    }
/* L720: */
	}
/* L300: */
    }
/* ...End of loop over outputs. */
/* ...Add Jacobian part: */
    i__1 = jnint2_1.m[jnint2_1.nl];
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = jnint2_1.m[jnint2_1.nl];
	for (j = 1; j <= i__2; ++j) {
	    dd[j - 1] = 0.f;
	    i__3 = jnint2_1.m[jnint2_1.nl - 1];
	    for (k = 1; k <= i__3; ++k) {
		dd[jnint2_1.m[jnint2_1.nl] + (j - 1) * jnint2_1.m[jnint2_1.nl 
			- 1] + k - 1] = 0.f;
/* L902: */
	    }
/* L901: */
	}
	mi = jnindx_(&jnint2_1.nl, &i__, &c__0);
	jnint1_1.d__[mi - 1] = jnsigm_1.gpjn[mi - 1];
	for (il = jnint2_1.nl - 1; il >= 1; --il) {
	    i__2 = jnint2_1.m[il];
	    for (j = 1; j <= i__2; ++j) {
		mj = jnint2_1.mv0[il - 1] + j;
		sum = 0.f;
		if (il < jnint2_1.nl - 1) {
		    mij = jnint2_1.mm0[il] + (j - 1) * jnint2_1.m[il + 1];
		    i__3 = jnint2_1.mv0[il] + jnint2_1.m[il + 1];
		    for (ii = jnint2_1.mv0[il] + 1; ii <= i__3; ++ii) {
			++mij;
			sum += jnint1_1.d__[ii - 1] * jnint1_1.w[mij - 1];
/* L930: */
		    }
		    jnint1_1.d__[mj - 1] = sum * jnsigm_1.gpjn[mj - 1];
		} else {
		    mij = jnindx_(&jnint2_1.nl, &i__, &j);
		    jnint1_1.d__[mj - 1] = jnint1_1.d__[mi - 1] * jnint1_1.w[
			    mij - 1] * jnsigm_1.gpjn[mj - 1];
		}
/* L920: */
	    }
/* L910: */
	}
	dd[i__ - 1] = jnint1_1.d__[mi - 1];
	if (jnint2_1.nl == 1) {
	    i__2 = jnint2_1.m[0];
	    for (j = 1; j <= i__2; ++j) {
		dd[jnint2_1.m[1] + (i__ - 1) * jnint2_1.m[0] + j - 1] = 
			jnint1_1.d__[mi - 1] * jndat1_1.oin[j - 1];
/* L940: */
	    }
	} else {
	    i__2 = jnint2_1.m[jnint2_1.nl - 1];
	    for (j = 1; j <= i__2; ++j) {
		mj = jnint2_1.mv0[jnint2_1.nl - 2] + j;
		dd[jnint2_1.m[jnint2_1.nl] + (i__ - 1) * jnint2_1.m[
			jnint2_1.nl - 1] + j - 1] = jnint1_1.d__[mi - 1] * 
			jnint1_1.o[mj - 1];
		dd[jnint2_1.m[jnint2_1.nl] + jnint2_1.m[jnint2_1.nl] * 
			jnint2_1.m[jnint2_1.nl - 1] + j - 1] = jnint1_1.d__[
			mj - 1];
/* L950: */
	    }
	    iofst = jnint2_1.m[jnint2_1.nl] + jnint2_1.m[jnint2_1.nl] * 
		    jnint2_1.m[jnint2_1.nl - 1] + jnint2_1.m[jnint2_1.nl - 1];
	    for (il = jnint2_1.nl - 2; il >= 1; --il) {
		i__2 = jnint2_1.m[il];
		for (k = 1; k <= i__2; ++k) {
		    mk = jnint2_1.mv0[il - 1] + k;
		    indx = iofst + jnint2_1.m[il] * jnint2_1.m[il + 1] + k;
		    dd[indx - 1] = jnint1_1.d__[mk - 1];
		    i__3 = jnint2_1.m[il + 1];
		    for (j = 1; j <= i__3; ++j) {
			mj = jnint2_1.mv0[il] + j;
			indx = iofst + (j - 1) * jnint2_1.m[il] + k;
			dd[indx - 1] = jnint1_1.d__[mj - 1] * jnint1_1.o[mk - 
				1];
/* L980: */
		    }
/* L970: */
		}
		iofst = iofst + jnint2_1.m[il] * jnint2_1.m[il + 1] + 
			jnint2_1.m[il];
/* L960: */
	    }
	    i__2 = jnint2_1.m[0];
	    for (k = 1; k <= i__2; ++k) {
		i__3 = jnint2_1.m[1];
		for (j = 1; j <= i__3; ++j) {
		    mj = jnint2_1.mv0[0] + j;
		    indx = iofst + (j - 1) * jnint2_1.m[0] + k;
		    dd[indx - 1] = jnint1_1.d__[mj - 1] * jndat1_1.oin[k - 1];
/* L1000: */
		}
/* L990: */
	    }
	}
	i__2 = nwgts;
	for (iw = 1; iw <= i__2; ++iw) {
	    i__3 = nwgts;
	    for (iv = iw; iv <= i__3; ++iv) {
		jnint5_1.d2e[iw + iv * 300 - 301] += dd[iw - 1] * dd[iv - 1];
/* L1020: */
	    }
/* L1010: */
	}
/* L900: */
    }
/* L10: */
    if (jndat1_1.mstjn[38] % (jndat1_1.mstjn[1] * jndat1_1.mstjn[8]) == 0) {
/* ...Symmetrize and normalize the Hessian. */
	nwgts = jnint2_1.mm0[jnint2_1.nl] + jnint2_1.mv0[jnint2_1.nl];
	i__1 = nwgts;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint5_1.d2e[i__ + i__ * 300 - 301] /= (real) (jndat1_1.mstjn[1] *
		     jndat1_1.mstjn[8]);
	    i__2 = nwgts;
	    for (j = i__ + 1; j <= i__2; ++j) {
		jnint5_1.d2e[i__ + j * 300 - 301] /= (real) (jndat1_1.mstjn[1]
			 * jndat1_1.mstjn[8]);
		jnint5_1.d2e[j + i__ * 300 - 301] = jnint5_1.d2e[i__ + j * 
			300 - 301];
/* L30: */
	    }
/* L20: */
	}
    }
    return 0;
/* **** END OF JNHESS **************************************************** */
} /* jnhess_ */

/* *********************************************************************** */
integer jnindx_(integer *il, integer *i__, integer *j)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    extern /* Subroutine */ int jnerr_(integer *);
    static integer ix, iy, jx, jy, inx, iny;

/* ...JetNet function INDeX */
/* ...Gives the node vector index of node I in layer IL for J=0 */
/* ...else gives the weight vector index of weight between node */
/* ...I of layer IL and node J of layer IL-1 */
    if (jndat1_1.mstjn[7] == 0) {
	jnerr_(&c__5);
    }
    if (*j == 0) {
	ret_val = jnint2_1.mv0[*il - 1] + *i__;
    } else {
	if (jnint3_1.nxin == 0 || *il > 1) {
	    ret_val = jnint2_1.mm0[*il - 1] + (*j - 1) * jnint2_1.m[*il] + *
		    i__;
	} else {
	    if (*i__ <= jnint3_1.nhrf * jnint3_1.nhprf) {
		if (*j <= (i__1 = jnint3_1.nxin * jnint3_1.nyin, abs(i__1))) {
		    ix = (*i__ - 1) / jnint3_1.nyhrf + 1;
		    iy = (*i__ - 1) % jnint3_1.nyhrf + 1;
		    inx = (*j - 1) / abs(jnint3_1.nyin) + 1;
		    iny = (*j - 1) % abs(jnint3_1.nyin) + 1;
		    jx = inx - ix + 1;
		    if (jx <= 0) {
			jx += jnint3_1.nxrf;
		    }
		    if (jx <= 0) {
			jnerr_(&c__12);
		    }
		    if (jx > jnint3_1.nxrf) {
			jnerr_(&c__12);
		    }
		    jy = iny - iy + 1;
		    if (jy <= 0) {
			jy += jnint3_1.nyrf;
		    }
		    if (jy <= 0) {
			jnerr_(&c__12);
		    }
		    if (jy > jnint3_1.nyrf) {
			jnerr_(&c__12);
		    }
		    ret_val = (jx - 1) * jnint3_1.nyrf + jy + (*i__ - 1) / 
			    jnint3_1.nhrf * jnint3_1.nrfw;
		} else {
		    ret_val = jnint3_1.nxrf * jnint3_1.nyrf + *j - (i__1 = 
			    jnint3_1.nxin * jnint3_1.nyin, abs(i__1)) + (*i__ 
			    - 1) / jnint3_1.nhrf * jnint3_1.nrfw;
		}
	    } else {
		ret_val = jnint3_1.nhprf * jnint3_1.nrfw + (*j - 1) * 
			jnint2_1.m[1] + *i__ - jnint3_1.nxhrf * 
			jnint3_1.nyhrf;
	    }
	}
    }
    return ret_val;
/* **** END OF JNINDX **************************************************** */
} /* jnindx_ */

/* *********************************************************************** */
/* Subroutine */ int jninit_(void)
{
    /* Format strings */
    static char fmt_600[] = "(22x,\002Weights and thresholds set randomly"
	    "\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer idum, i__, j, ihprf;
    static real width;
    extern /* Subroutine */ int jnerr_(integer *);
    static integer il;
    extern /* Subroutine */ int jnhead_(void), jnsepa_(void), jnstat_(integer 
	    *);
    static real sumrfw;
    static integer mij;
    extern doublereal rjn_(integer *);

    /* Fortran I/O blocks */
    static cilist io___219 = { 0, 0, 0, fmt_600, 0 };


/* ...JetNet subroutine INITialize net */
/* ...Initializes a net according to switches and parameters in */
/* .../JNDAT1/ and /JNDAT2/ */
/* ...Check if JMINIT has been called */
    if (jndat1_1.mstjm[7] == 1) {
	jnerr_(&c__16);
    }
/* ...Set parameters in /JNINT2/ */
    jnsepa_();
/* ...Set initial values of weights and thresholds */
    i__1 = jnint2_1.nl;
    for (il = 1; il <= i__1; ++il) {
/* ...Set width in this layer */
	if (jndat2_1.widl[il - 1] <= 0.f) {
	    width = jndat1_1.parjn[3];
	} else {
	    width = jndat2_1.widl[il - 1];
	}
/* ...Initialize weights */
	i__2 = jnint2_1.mm0[il];
	for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__2; ++i__) {
	    idum = i__;
	    if (width >= 0.f) {
		jnint1_1.w[i__ - 1] = (rjn_(&idum) * 2.f - 1.f) * width;
	    } else {
		jnint1_1.w[i__ - 1] = -rjn_(&idum) * width;
	    }
/* L110: */
	}
/* ...Initialize thresholds */
	i__2 = jnint2_1.mv0[il];
	for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__2; ++i__) {
	    idum = i__;
	    if (width >= 0.f) {
		jnint1_1.t[i__ - 1] = (rjn_(&idum) * 2.f - 1.f) * width;
	    } else {
		jnint1_1.t[i__ - 1] = -rjn_(&idum) * width;
	    }
/* L120: */
	}
/* L100: */
    }
    if (jnint3_1.nxin != 0) {
	i__1 = jnint3_1.nhprf;
	for (ihprf = 1; ihprf <= i__1; ++ihprf) {
	    i__2 = jnint3_1.nhrf;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		jnint1_1.t[(ihprf - 1) * jnint3_1.nhrf + i__ - 1] = 
			jnint1_1.t[(ihprf - 1) * jnint3_1.nhrf];
/* L140: */
	    }
/* L130: */
	}
	if (jndat1_1.mstjn[26] < 0) {
	    i__1 = jnint2_1.m[2];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = jnint3_1.nhprf;
		for (ihprf = 1; ihprf <= i__2; ++ihprf) {
		    mij = jnint2_1.mm0[1] + (ihprf - 1) * jnint3_1.nhrf * 
			    jnint2_1.m[2] + i__;
		    sumrfw = jnint1_1.w[mij - 1];
		    i__3 = jnint3_1.nhrf;
		    for (j = 2; j <= i__3; ++j) {
			mij += jnint2_1.m[2];
			jnint1_1.w[mij - 1] = sumrfw;
/* L170: */
		    }
/* L160: */
		}
/* L150: */
	    }
	}
    }
/* ...Write statistics on output file */
    if (jndat1_1.mstjn[5] < 0) {
	return 0;
    }
    jnhead_();
    jnstat_(&c__1);
    io___219.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___219);
    e_wsfe();
    return 0;
/* **** END OF JNINIT **************************************************** */
} /* jninit_ */

/* *********************************************************************** */
/* Subroutine */ int jnlins_(void)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3;

    /* Builtin functions */
    double r_sign(real *, real *);

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int jnerr_(integer *);
    static real ac, bc, factor, de2, errnow, eta, tol;

/* ...JetNet subroutine do LINe Search. */
/* ...Performs a line search in the direction of G. */
/* ...The algorithm is a mixture between 'golden section search' and */
/* ...quadratic interpolation. Termination of the search is controlled */
/* ...by either of two criteria: (1) If the error has decreased */
/* ...sufficiently much - set by PARJN(24); or (2) if the predicted */
/* ...location of the error is within the preset - PARJN(25) - tolerance */
/* ...distance from the current best point. */
/* ...The first step is always equal to PARJN(1), but PARJN(1) is set to */
/* ...about half the size of the last successful step every time */
/* ...the algorithm finds a minimum (provided that this step size is */
/* ...smaller than the maximum allowed step size). */
/* ...ERRLN(1) = error value in current point. */
/* ...ERRLN(2-3) = error values in previous points. */
/* ...ERRLN(0) = error value in the starting point. */
/* ...STEPLN(1) = step to be taken (the current point is always at x=0). */
/* ...STEPLN(2-3) = distance to previous points. */
/* ...STEPLN(0) = distance to starting point. */
/* ...STEPMN = distance to best minimum so far. */
/* ...PARJN(26)=minimum allowed relative change in error. */
/* ...PARJN(27)=maximum allowed step size. */
/* ...ZEPS=machine precision. */
/* ...TINY=small number to prevent division by zero. */
    if (jndat1_1.mstjn[4] == 8) {
/* ...Back up to best point so far and terminate updating. */
	if (jndat1_1.parjn[7] <= jnint4_1.errmn) {
	    jnint4_1.errmn = jndat1_1.parjn[7];
	    jnint4_1.stepmn = 0.f;
	}
	jnint4_1.stepln[1] = -jnint4_1.stepmn;
/* ...Freeze updating: */
	jndat1_1.mstjn[4] = 9;
	jnint4_1.ilinon = 0;
	jnint4_1.nc = 0;
	jnint4_1.nit = 0;
	jnint4_1.nsc = 0;
	goto L20;
    }
    ++jnint4_1.nit;
    if (jnint4_1.nit > jndat1_1.mstjn[34]) {
/* ...Too many iterations -> restart from best minimum so far */
/* ...and rescale PARJN(1) according to previous success. */
	if (jndat1_1.mstjn[37] >= jndat1_1.mstjn[35]) {
	    jnerr_(&c__21);
	}
	++jndat1_1.mstjn[37];
	if (jnint4_1.ilinon > 0) {
	    if (jnint4_1.errln[1] <= jnint4_1.errln[0]) {
		if (dabs(jnint4_1.stepln[0]) > 1e-20f) {
/* Computing MIN */
/* Computing MAX */
		    r__2 = dabs(jnint4_1.stepln[0]), r__3 = jndat1_1.parjn[0] 
			    * 1.618034f;
		    r__1 = dmax(r__2,r__3);
		    jndat1_1.parjn[0] = dmin(r__1,jndat1_1.parjn[26]);
		} else {
/* Computing MIN */
		    r__1 = jndat1_1.parjn[0] * 1.618034f;
		    jndat1_1.parjn[0] = dmin(r__1,jndat1_1.parjn[26]);
		}
	    } else {
		jndat1_1.parjn[0] *= .381966f;
	    }
	} else {
	    jndat1_1.parjn[0] *= .381966f;
	}
	jndat1_1.parjn[0] = dmax(jndat1_1.parjn[0],1e-8f);
	jnint4_1.stepln[1] = -jnint4_1.stepmn;
	jnint4_1.ilinon = 0;
	jnint4_1.nc = 0;
	jnint4_1.nit = 0;
	goto L20;
    }
    eta = 0.f;
L10:
    if (jnint4_1.nit < 3) {
/* ...At least 3 points are needed */
	jnint4_1.errln[2] = jnint4_1.errln[1];
/* ...Store last updated error */
	jnint4_1.errln[1] = jndat1_1.parjn[7];
	if (jnint4_1.nit == 1) {
	    jnint4_1.errln[0] = jnint4_1.errln[1];
	    jnint4_1.errmn = jnint4_1.errln[1];
	    jnint4_1.stepln[0] = 0.f;
	    jnint4_1.stepmn = 0.f;
	    for (i__ = 2; i__ <= 3; ++i__) {
		jnint4_1.stepln[i__] = 0.f;
		jnint4_1.errln[i__] = 0.f;
/* L100: */
	    }
	    if (dabs(jnint4_1.derrln) > 1e-20f) {
		jnint4_1.stepln[1] = -r_sign(jndat1_1.parjn, &jnint4_1.derrln)
			;
	    } else {
		jnint4_1.stepln[1] = jndat1_1.parjn[0];
	    }
	} else {
	    jnint4_1.stepln[2] = -jnint4_1.stepln[1];
/* Computing 2nd power */
	    r__1 = jnint4_1.stepln[2];
	    de2 = (jnint4_1.errln[1] - jnint4_1.errln[2] + jnint4_1.stepln[2] 
		    * jnint4_1.derrln) * 2.f / (r__1 * r__1);
	    if ((r__1 = jnint4_1.derrln / jnint4_1.stepln[2], dabs(r__1)) < 
		    de2 * 10.f) {
		jnint4_1.stepln[1] = -jnint4_1.derrln / de2;
	    } else {
		jnint4_1.stepln[1] *= 1.618034f;
	    }
	}
	if (jnint4_1.errln[1] < jnint4_1.errmn) {
	    jnint4_1.stepmn = 0.f;
	    jnint4_1.errmn = jnint4_1.errln[1];
	}
    } else if (jnint4_1.ilinon > 0) {
/* ...Bracket the minimum */
/* ...Update error and step values: */
	jnint4_1.errln[3] = jnint4_1.errln[2];
	jnint4_1.errln[2] = jnint4_1.errln[1];
	jnint4_1.stepln[3] = -jnint4_1.stepln[1] + jnint4_1.stepln[2];
	jnint4_1.stepln[2] = -jnint4_1.stepln[1];
	jnint4_1.stepln[1] = 0.f;
	jnint4_1.errln[1] = jndat1_1.parjn[7];
	if (jnint4_1.errln[1] < jnint4_1.errmn) {
	    jnint4_1.stepmn = 0.f;
	    jnint4_1.errmn = jnint4_1.errln[1];
	}
/* ...Check if the search is improving - else take default step. */
	if ((r__1 = 1.f - jnint4_1.errln[0] / jnint4_1.errln[1], dabs(r__1)) <
		 jndat1_1.parjn[25]) {
	    jnint4_1.stepln[1] = -jnint4_1.stepln[2] * 1.618034f;
	    goto L20;
	}
/* ...Quadratic fit */
	if ((jnint4_1.errln[1] - jnint4_1.errln[3]) * 1e-20f > 
		jnint4_1.stepln[3]) {
	    factor = -1.f;
	} else {
	    bc = ((jnint4_1.errln[1] - jnint4_1.errln[3]) / (jnint4_1.stepln[
		    3] + 1e-20f) - (jnint4_1.errln[1] - jnint4_1.errln[2]) / (
		    jnint4_1.stepln[2] + 1e-20f)) / (jnint4_1.stepln[2] - 
		    jnint4_1.stepln[3] + 1e-20f);
	    ac = -bc * jnint4_1.stepln[2] - (jnint4_1.errln[1] - 
		    jnint4_1.errln[2]) / (jnint4_1.stepln[2] + 1e-20f);
	    eta = -ac / (bc * 2.f + 1e-20f);
	    if (dabs(eta) > 1e-20f) {
		factor = (jnint4_1.errln[1] - jnint4_1.errln[2]) / (
			jnint4_1.stepln[2] * (eta * 2.f - jnint4_1.stepln[2]) 
			+ 1e-20f);
	    } else {
		factor = -1.f;
	    }
	}
	if (jnint4_1.errln[1] < jnint4_1.errln[2]) {
	    if (jnint4_1.errln[1] < jnint4_1.errln[3]) {
		if (jnint4_1.stepln[2] * jnint4_1.stepln[3] < 0.f) {
/* ...Minimum is bracketed -> find it */
		    jnint4_1.ilinon = -1;
		    goto L10;
		} else {
/* ...Keep searching: */
		    if (factor > 0.f) {
/* ...Quadratic fit OK. */
			if ((r__1 = eta / jnint4_1.stepln[2], dabs(r__1)) < 
				10.f || (r__2 = eta / jnint4_1.stepln[3], 
				dabs(r__2)) < 10.f) {
			    jnint4_1.stepln[1] = eta;
			} else {
			    jnint4_1.stepln[1] = -jnint4_1.stepln[2] * 
				    1.618034f;
			}
		    } else {
			jnint4_1.stepln[1] = -jnint4_1.stepln[2] * 1.618034f;
		    }
		}
	    } else {
		if (jnint4_1.stepln[2] / jnint4_1.stepln[3] > 1.f) {
/* ...Back up to point (3) */
		    jnint4_1.stepln[1] = jnint4_1.stepln[3];
		} else {
/* ...Back up beyond point (3) */
		    jnint4_1.stepln[1] = jnint4_1.stepln[3] * 1.618034f;
		    jnint4_1.stepln[2] = jnint4_1.stepln[3];
		    jnint4_1.errln[2] = jnint4_1.errln[3];
		}
	    }
	} else if (jnint4_1.errln[1] > jnint4_1.errln[2]) {
	    if (jnint4_1.errln[2] < jnint4_1.errln[3]) {
		if (jnint4_1.stepln[3] / jnint4_1.stepln[2] > 1.f) {
/* ...Minimum is bracketed -> back up towards point (2) */
		    if (eta / jnint4_1.stepln[2] > 0.f && eta / 
			    jnint4_1.stepln[2] < 1.f) {
			jnint4_1.stepln[1] = eta;
		    } else {
			jnint4_1.stepln[1] = jnint4_1.stepln[2];
			jnint4_1.stepln[2] = jnint4_1.stepln[3];
			jnint4_1.errln[2] = jnint4_1.errln[3];
		    }
		    jnint4_1.ilinon = -2;
		} else {
		    jnint4_1.stepln[1] = jnint4_1.stepln[2] * 1.618034f;
		}
	    } else {
/* ...Rearrange and move beyond point (3) */
		jnint4_1.stepln[1] = jnint4_1.stepln[3] * 1.618034f;
		jnint4_1.stepln[2] = jnint4_1.stepln[3];
		jnint4_1.errln[2] = jnint4_1.errln[3];
	    }
	} else {
/* ...Take default step */
	    jnint4_1.stepln[1] = -jnint4_1.stepln[2] * 1.618034f;
	}
    } else if (jnint4_1.ilinon < 0) {
/* ...Find minimum (knowing that minimum is bracketed) */
	errnow = jndat1_1.parjn[7];
	if (errnow < jnint4_1.errmn) {
	    jnint4_1.stepmn = 0.f;
	    jnint4_1.errmn = errnow;
	}
/* ...Check bracket condition: */
	if (errnow >= jnint4_1.errln[2] || errnow >= jnint4_1.errln[3]) {
	    jnint4_1.ilinon = 1;
	    goto L10;
	}
	if (jnint4_1.errln[0] - errnow <= jndat1_1.parjn[23] * 
		jnint4_1.stepln[0] * jnint4_1.derrln) {
/* ...Satisfactory -> terminate search */
	    jnint4_1.nit = 0;
	    jnint4_1.ilinon = 0;
	    if (dabs(jnint4_1.stepln[0]) > 1e-8f) {
/* Computing MAX */
/* Computing MIN */
		r__2 = dabs(jnint4_1.stepln[0]) * .381966f;
		r__1 = dmin(r__2,jndat1_1.parjn[26]);
		jndat1_1.parjn[0] = dmax(r__1,1e-8f);
	    }
	    if (errnow > jnint4_1.errmn) {
/* ...Back up to best minimum so far */
		jnint4_1.stepln[1] = -jnint4_1.stepmn;
		goto L20;
	    } else {
		jndat1_1.mstjn[36] = 0;
		return 0;
	    }
	} else {
	    if (jnint4_1.ilinon != -1) {
/* ...Rearrange points: */
		if (errnow <= jnint4_1.errln[1]) {
		    if (jnint4_1.stepln[1] / (jnint4_1.stepln[2] + 1e-20f) > 
			    0.f) {
			jnint4_1.stepln[3] = -jnint4_1.stepln[1];
			jnint4_1.errln[3] = jnint4_1.errln[1];
			jnint4_1.stepln[2] -= jnint4_1.stepln[1];
			jnint4_1.stepln[1] = 0.f;
			jnint4_1.errln[1] = errnow;
		    } else {
			jnint4_1.stepln[3] -= jnint4_1.stepln[1];
			jnint4_1.stepln[2] = -jnint4_1.stepln[1];
			jnint4_1.errln[2] = jnint4_1.errln[1];
			jnint4_1.stepln[1] = 0.f;
			jnint4_1.errln[1] = errnow;
		    }
		} else {
		    if (jnint4_1.stepln[1] / (jnint4_1.stepln[2] + 1e-20f) > 
			    0.f) {
			jnint4_1.stepln[2] = jnint4_1.stepln[1];
			jnint4_1.errln[2] = errnow;
			jnint4_1.stepln[1] = 0.f;
		    } else {
			jnint4_1.stepln[3] = jnint4_1.stepln[1];
			jnint4_1.errln[3] = errnow;
			jnint4_1.stepln[1] = 0.f;
		    }
		}
	    } else {
		jnint4_1.stepln[1] = 0.f;
		jnint4_1.ilinon = -2;
	    }
/* ...Quadratic fit */
	    if ((jnint4_1.errln[1] - jnint4_1.errln[3]) * 1e-20f > 
		    jnint4_1.stepln[3]) {
		factor = -1.f;
	    } else {
		bc = ((jnint4_1.errln[1] - jnint4_1.errln[3]) / (
			jnint4_1.stepln[3] + 1e-20f) - (jnint4_1.errln[1] - 
			jnint4_1.errln[2]) / (jnint4_1.stepln[2] + 1e-20f)) / 
			(jnint4_1.stepln[2] - jnint4_1.stepln[3] + 1e-20f);
		ac = -bc * jnint4_1.stepln[2] - (jnint4_1.errln[1] - 
			jnint4_1.errln[2]) / (jnint4_1.stepln[2] + 1e-20f);
		eta = -ac / (bc * 2.f + 1e-20f);
		if (dabs(eta) > 1e-20f) {
		    factor = (jnint4_1.errln[1] - jnint4_1.errln[2]) / (
			    jnint4_1.stepln[2] * (eta * 2.f - jnint4_1.stepln[
			    2]) + 1e-20f);
		} else {
		    factor = -1.f;
		}
	    }
/* ..Tolerance: */
	    tol = dmax(jndat1_1.parjn[24],1e-8f);
	    if (factor > 0.f) {
/* ...Quadratic fit OK */
		if (eta / (jnint4_1.stepln[2] + 1e-20f) > 0.f) {
		    if (eta / (jnint4_1.stepln[2] + 1e-20f) < 1.f && (r__1 = 
			    eta - jnint4_1.stepln[2], dabs(r__1)) > tol) {
			jnint4_1.stepln[1] = eta;
		    } else {
			jnint4_1.stepln[1] = jnint4_1.stepln[2] * .381966f;
		    }
		} else if (eta / (jnint4_1.stepln[3] + 1e-20f) > 0.f) {
		    if (eta / (jnint4_1.stepln[3] + 1e-20f) < 1.f && (r__1 = 
			    eta - jnint4_1.stepln[3], dabs(r__1)) > tol) {
			jnint4_1.stepln[1] = eta;
		    } else {
			jnint4_1.stepln[1] = jnint4_1.stepln[3] * .381966f;
		    }
		} else {
/* ...Step too large -> decrease */
		    r__1 = dmin(jnint4_1.stepln[2],jnint4_1.stepln[3]);
		    jnint4_1.stepln[1] = r_sign(&r__1, &eta) * .381966f;
		}
	    } else {
/* ...Take step towards the most distant of points (2) and (3) */
		if (jnint4_1.stepln[2] / (jnint4_1.stepln[3] + 1e-20f) < -1.f)
			 {
		    jnint4_1.stepln[1] = jnint4_1.stepln[2] * .381966f;
		} else {
		    jnint4_1.stepln[1] = jnint4_1.stepln[3] * .381966f;
		}
	    }
	}
	if (dabs(jnint4_1.stepln[1]) <= tol) {
/* ...Predicted step less than tolerance from current point */
/* ...                                  -> terminate search */
	    jnint4_1.nit = 0;
	    jnint4_1.ilinon = 0;
	    if (dabs(jnint4_1.stepln[0]) > 1e-8f) {
/* Computing MAX */
/* Computing MIN */
		r__2 = dabs(jnint4_1.stepln[0]) * .381966f;
		r__1 = dmin(r__2,jndat1_1.parjn[26]);
		jndat1_1.parjn[0] = dmax(r__1,1e-8f);
	    }
	    if (errnow > jnint4_1.errmn) {
/* ...Back up to best minimum so far */
		jnint4_1.stepln[1] = -jnint4_1.stepmn;
	    } else {
		jndat1_1.mstjn[36] = abs(jnint4_1.ilinon);
		return 0;
	    }
	}
    }
L20:
/* ...Update weight vector: */
    i__1 = jnint2_1.mm0[jnint2_1.nl];
    for (i__ = 1; i__ <= i__1; ++i__) {
	jnint1_1.w[i__ - 1] += jnint4_1.stepln[1] * jnint1_1.g[i__ - 1] * (
		real) jnint1_1.nself[i__ - 1];
/* L110: */
    }
    i__1 = jnint2_1.mv0[jnint2_1.nl];
    for (i__ = 1; i__ <= i__1; ++i__) {
	jnint1_1.t[i__ - 1] += jnint4_1.stepln[1] * jnint1_1.g[i__ + 
		jnint2_1.mm0[jnint2_1.nl] - 1] * (real) jnint1_1.ntself[i__ - 
		1];
/* L120: */
    }
/* ...Keep track of starting point and best point up to now: */
    jnint4_1.stepln[0] += jnint4_1.stepln[1];
    jnint4_1.stepmn += jnint4_1.stepln[1];
    jndat1_1.mstjn[36] = abs(jnint4_1.ilinon);
    return 0;
/* **** END OF JNLINS **************************************************** */
} /* jnlins_ */

/* *********************************************************************** */
/* Subroutine */ int jnread_(integer *nf)
{
    /* Format strings */
    static char fmt_690[] = "(a)";
    static char fmt_650[] = "(tr63,f3.1)";
    static char fmt_601[] = "(tr11,6i7,f7.3,3i7)";
    static char fmt_600[] = "(tr11,10i7)";
    static char fmt_610[] = "(tr11,10f7.4)";
    static char fmt_620[] = "(6(f12.4,a1))";
    static char fmt_640[] = "(29x,\002Weights read from file\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_rsue(cilist *), do_uio(integer *, char *, ftnlen), e_rsue(void),
	     s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void)
	    , s_cmp(char *, char *, ftnlen, ftnlen), s_rsle(cilist *), e_rsle(
	    void);
    double pow_dd(doublereal *, doublereal *);
    integer s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static char line[100];
    static integer i__, j, ihprf;
    extern /* Subroutine */ int jnerr_(integer *);
    static real fvers;
    static integer ivers, jf, il;
    extern /* Subroutine */ int jnhead_(void);
    static integer iw;
    extern /* Subroutine */ int jnsepa_(void);
    static integer nfsave;
    extern integer jnindx_(integer *, integer *, integer *);
    extern /* Subroutine */ int jnstat_(integer *);
    static real par[21];
    static integer mst[6];
    static real trn;

    /* Fortran I/O blocks */
    static cilist io___230 = { 0, 0, 0, 0, 0 };
    static cilist io___235 = { 0, 0, 0, 0, 0 };
    static cilist io___236 = { 0, 0, 0, 0, 0 };
    static cilist io___237 = { 0, 0, 0, 0, 0 };
    static cilist io___238 = { 0, 0, 0, 0, 0 };
    static cilist io___239 = { 0, 0, 0, 0, 0 };
    static cilist io___240 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___242 = { 0, 0, 0, 0, 0 };
    static cilist io___243 = { 0, 0, 0, 0, 0 };
    static cilist io___244 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___246 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___247 = { 0, 0, 0, fmt_601, 0 };
    static cilist io___249 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___250 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___251 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___252 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___253 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___254 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___255 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___256 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___257 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___258 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___259 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___260 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___261 = { 0, 0, 0, 0, 0 };
    static cilist io___262 = { 0, 0, 0, 0, 0 };
    static cilist io___263 = { 0, 0, 0, 0, 0 };
    static cilist io___264 = { 0, 0, 0, 0, 0 };
    static cilist io___265 = { 0, 0, 0, 0, 0 };
    static cilist io___267 = { 0, 0, 0, 0, 0 };
    static cilist io___268 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___269 = { 0, 0, 0, 0, 0 };
    static cilist io___271 = { 0, 0, 0, 0, 0 };
    static cilist io___272 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___274 = { 0, 0, 0, 0, 0 };
    static cilist io___275 = { 0, 0, 0, 0, 0 };
    static cilist io___276 = { 0, 0, 0, 0, 0 };
    static cilist io___277 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___278 = { 0, 0, 0, 0, 0 };
    static cilist io___279 = { 0, 0, 0, 0, 0 };
    static cilist io___280 = { 0, 0, 0, 0, 0 };
    static cilist io___281 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___283 = { 0, 0, 0, 0, 0 };
    static cilist io___284 = { 0, 0, 0, 0, 0 };
    static cilist io___285 = { 0, 0, 0, 0, 0 };
    static cilist io___286 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___287 = { 0, 0, 0, 0, 0 };
    static cilist io___288 = { 0, 0, 0, 0, 0 };
    static cilist io___289 = { 0, 0, 0, 0, 0 };
    static cilist io___290 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___291 = { 0, 0, 0, fmt_640, 0 };


/* ...JetNet subroutine READ weights and parameters. */
/* ...Reads weights, thresholds and other statistics from a file NF and */
/* ...initializes the net */
    if (*nf < 0) {
/* ...unformatted read */
	nfsave = jndat1_1.mstjn[5];
	jf = -(*nf);
	io___230.ciunit = jf;
	s_rsue(&io___230);
	do_uio(&c__1, (char *)&ivers, (ftnlen)sizeof(integer));
	e_rsue();
	if (ivers < 0) {
	    jnerr_(&c__17);
	}
	if (ivers < 20) {
	    jnerr_(&c__13);
	}
	if (ivers / 10 == 2) {
/* ...New meanings for MSTJN(35-) and PARJN(20-) */
	    for (i__ = 1; i__ <= 6; ++i__) {
		mst[i__ - 1] = jndat1_1.mstjn[i__ + 33];
/* L400: */
	    }
	    for (i__ = 1; i__ <= 21; ++i__) {
		par[i__ - 1] = jndat1_1.parjn[i__ + 18];
/* L410: */
	    }
	}
	io___235.ciunit = jf;
	s_rsue(&io___235);
	do_uio(&c__40, (char *)&jndat1_1.mstjn[0], (ftnlen)sizeof(integer));
	do_uio(&c__40, (char *)&jndat1_1.parjn[0], (ftnlen)sizeof(real));
	do_uio(&c__10, (char *)&jndat2_1.tinv[0], (ftnlen)sizeof(real));
	do_uio(&c__10, (char *)&jndat2_1.igfn[0], (ftnlen)sizeof(integer));
	do_uio(&c__10, (char *)&jndat2_1.etal[0], (ftnlen)sizeof(real));
	do_uio(&c__10, (char *)&jndat2_1.widl[0], (ftnlen)sizeof(real));
	do_uio(&c__10, (char *)&jndat2_1.satm[0], (ftnlen)sizeof(real));
	e_rsue();
	if (ivers / 10 == 2) {
	    for (i__ = 1; i__ <= 6; ++i__) {
		jndat1_1.mstjn[i__ + 33] = mst[i__ - 1];
/* L420: */
	    }
	    for (i__ = 1; i__ <= 21; ++i__) {
		jndat1_1.parjn[i__ + 18] = par[i__ - 1];
/* L430: */
	    }
	}
	jnsepa_();
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___236.ciunit = jf;
	    s_rsue(&io___236);
	    do_uio(&c__1, (char *)&jnint1_1.w[i__ - 1], (ftnlen)sizeof(real));
	    e_rsue();
/* L100: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___237.ciunit = jf;
	    s_rsue(&io___237);
	    do_uio(&c__1, (char *)&jnint1_1.t[i__ - 1], (ftnlen)sizeof(real));
	    e_rsue();
/* L110: */
	}
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___238.ciunit = jf;
	    s_rsue(&io___238);
	    do_uio(&c__1, (char *)&jnint1_1.nself[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_rsue();
/* L120: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___239.ciunit = jf;
	    s_rsue(&io___239);
	    do_uio(&c__1, (char *)&jnint1_1.ntself[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_rsue();
/* L130: */
	}
	jndat1_1.mstjn[5] = nfsave;
    } else {
/* ...Formatted dump */
	io___240.ciunit = *nf;
	s_rsfe(&io___240);
	do_fio(&c__1, line, (ftnlen)100);
	e_rsfe();
	if (s_cmp(line + 26, " D", (ftnlen)2, (ftnlen)2) == 0) {
	    jnerr_(&c__17);
	}
	io___242.ciunit = *nf;
	s_rsle(&io___242);
	e_rsle();
	io___243.ciunit = *nf;
	s_rsle(&io___243);
	e_rsle();
	io___244.ciunit = *nf;
	s_rsfe(&io___244);
	do_fio(&c__1, (char *)&fvers, (ftnlen)sizeof(real));
	e_rsfe();
	ivers = (integer) (fvers * 10.f + .001f);
	if (ivers < 20) {
	    jnerr_(&c__13);
	}
	if (ivers / 10 == 2) {
/* ...New meanings for MSTJN(35-) and PARJN(20-) */
	    for (i__ = 1; i__ <= 6; ++i__) {
		mst[i__ - 1] = jndat1_1.mstjn[i__ + 33];
/* L440: */
	    }
	    for (i__ = 1; i__ <= 21; ++i__) {
		par[i__ - 1] = jndat1_1.parjn[i__ + 18];
/* L450: */
	    }
	}
L900:
	io___246.ciunit = *nf;
	s_rsfe(&io___246);
	do_fio(&c__1, line, (ftnlen)100);
	e_rsfe();
	if (s_cmp(line, "        I       1      2", (ftnlen)24, (ftnlen)24) !=
		 0) {
	    goto L900;
	}
	nfsave = jndat1_1.mstjn[5];
	io___247.ciunit = *nf;
	s_rsfe(&io___247);
	for (i__ = 1; i__ <= 6; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	do_fio(&c__1, (char *)&trn, (ftnlen)sizeof(real));
	for (i__ = 8; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	e_rsfe();
	d__1 = (doublereal) trn;
	jndat1_1.mstjn[6] = (integer) (pow_dd(&c_b469, &d__1) + .5f);
	io___249.ciunit = *nf;
	s_rsfe(&io___249);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ + 9], (ftnlen)sizeof(
		    integer));
	}
	e_rsfe();
	io___250.ciunit = *nf;
	s_rsfe(&io___250);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ + 19], (ftnlen)sizeof(
		    integer));
	}
	e_rsfe();
	io___251.ciunit = *nf;
	s_rsfe(&io___251);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ + 29], (ftnlen)sizeof(
		    integer));
	}
	e_rsfe();
	io___252.ciunit = *nf;
	s_rsfe(&io___252);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjn[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	io___253.ciunit = *nf;
	s_rsfe(&io___253);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjn[i__ + 9], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	io___254.ciunit = *nf;
	s_rsfe(&io___254);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjn[i__ + 19], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	io___255.ciunit = *nf;
	s_rsfe(&io___255);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjn[i__ + 29], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	d__1 = (doublereal) jndat1_1.parjn[21];
	jndat1_1.parjn[21] = pow_dd(&c_b469, &d__1);
	io___256.ciunit = *nf;
	s_rsfe(&io___256);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.igfn[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	e_rsfe();
	io___257.ciunit = *nf;
	s_rsfe(&io___257);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.tinv[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	io___258.ciunit = *nf;
	s_rsfe(&io___258);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.etal[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	io___259.ciunit = *nf;
	s_rsfe(&io___259);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.widl[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	io___260.ciunit = *nf;
	s_rsfe(&io___260);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.satm[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	io___261.ciunit = *nf;
	s_rsle(&io___261);
	e_rsle();
	jndat1_1.mstjn[5] = nfsave;
	if (ivers / 10 == 2) {
	    for (i__ = 1; i__ <= 6; ++i__) {
		jndat1_1.mstjn[i__ + 33] = mst[i__ - 1];
/* L460: */
	    }
	    for (i__ = 1; i__ <= 21; ++i__) {
		jndat1_1.parjn[i__ + 18] = par[i__ - 1];
/* L470: */
	    }
	}
	jnsepa_();
	io___262.ciunit = *nf;
	s_rsle(&io___262);
	e_rsle();
	io___263.ciunit = *nf;
	s_rsle(&io___263);
	e_rsle();
	io___264.ciunit = *nf;
	s_rsle(&io___264);
	e_rsle();
	if (jnint3_1.nxin == 0) {
	    io___265.ciunit = *nf;
	    s_rsle(&io___265);
	    e_rsle();
	    i__1 = jnint2_1.m[0];
	    for (j = 1; j <= i__1; ++j) {
		io___267.ciunit = *nf;
		s_rsle(&io___267);
		e_rsle();
		io___268.ciunit = *nf;
		s_rsfe(&io___268);
		i__2 = jnint2_1.m[1];
		for (i__ = 1; i__ <= i__2; ++i__) {
		    do_fio(&c__1, (char *)&jnint1_1.w[jnindx_(&c__1, &i__, &j)
			     - 1], (ftnlen)sizeof(real));
		    do_fio(&c__1, line + (i__ - 1), (ftnlen)1);
		}
		e_rsfe();
		i__2 = jnint2_1.m[1];
		for (i__ = 1; i__ <= i__2; ++i__) {
		    if (*(unsigned char *)&line[i__ - 1] == '*') {
			jnint1_1.nself[jnindx_(&c__1, &i__, &j) - 1] = 0;
		    } else {
			jnint1_1.nself[jnindx_(&c__1, &i__, &j) - 1] = 1;
		    }
/* L210: */
		}
/* L200: */
	    }
	} else {
	    io___269.ciunit = *nf;
	    s_rsle(&io___269);
	    e_rsle();
	    i__1 = jnint3_1.nhprf;
	    for (ihprf = 1; ihprf <= i__1; ++ihprf) {
		io___271.ciunit = *nf;
		s_rsle(&io___271);
		e_rsle();
		io___272.ciunit = *nf;
		s_rsfe(&io___272);
		i__2 = jnint3_1.nrfw * ihprf;
		for (iw = jnint3_1.nrfw * (ihprf - 1) + 1; iw <= i__2; ++iw) {
		    do_fio(&c__1, (char *)&jnint1_1.w[iw - 1], (ftnlen)sizeof(
			    real));
		    do_fio(&c__1, line + (iw - 1), iw - (iw - 1));
		}
		e_rsfe();
		i__2 = jnint3_1.nrfw * ihprf;
		for (iw = jnint3_1.nrfw * (ihprf - 1) + 1; iw <= i__2; ++iw) {
		    if (*(unsigned char *)&line[iw - 1] == '*') {
			jnint1_1.nself[iw - 1] = 0;
		    } else {
			jnint1_1.nself[iw - 1] = 1;
		    }
/* L230: */
		}
/* L220: */
	    }
	    if (jnint3_1.nhrf * jnint3_1.nhprf < jnint2_1.m[1]) {
		io___274.ciunit = *nf;
		s_rsle(&io___274);
		e_rsle();
		io___275.ciunit = *nf;
		s_rsle(&io___275);
		e_rsle();
		i__1 = jnint2_1.m[0];
		for (j = 1; j <= i__1; ++j) {
		    io___276.ciunit = *nf;
		    s_rsle(&io___276);
		    e_rsle();
		    io___277.ciunit = *nf;
		    s_rsfe(&io___277);
		    i__2 = jnint2_1.m[1];
		    for (i__ = jnint3_1.nhrf * jnint3_1.nhprf + 1; i__ <= 
			    i__2; ++i__) {
			do_fio(&c__1, (char *)&jnint1_1.w[jnindx_(&c__1, &i__,
				 &j) - 1], (ftnlen)sizeof(real));
			do_fio(&c__1, line + (i__ - 1), (ftnlen)1);
		    }
		    e_rsfe();
		    i__2 = jnint2_1.m[1];
		    for (i__ = jnint3_1.nhrf * jnint3_1.nhprf + 1; i__ <= 
			    i__2; ++i__) {
			if (*(unsigned char *)&line[i__ - 1] == '*') {
			    jnint1_1.nself[jnindx_(&c__1, &i__, &j) - 1] = 0;
			} else {
			    jnint1_1.nself[jnindx_(&c__1, &i__, &j) - 1] = 1;
			}
/* L250: */
		    }
/* L240: */
		}
	    }
	}
	io___278.ciunit = *nf;
	s_rsle(&io___278);
	e_rsle();
	io___279.ciunit = *nf;
	s_rsle(&io___279);
	e_rsle();
	io___280.ciunit = *nf;
	s_rsle(&io___280);
	e_rsle();
	io___281.ciunit = *nf;
	s_rsfe(&io___281);
	i__1 = jnint2_1.m[1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&jnint1_1.t[jnindx_(&c__1, &i__, &c__0) - 1]
		    , (ftnlen)sizeof(real));
	    do_fio(&c__1, line + (i__ - 1), (ftnlen)1);
	}
	e_rsfe();
	i__1 = jnint2_1.m[1];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (*(unsigned char *)&line[i__ - 1] == '*') {
		jnint1_1.ntself[jnindx_(&c__1, &i__, &c__0) - 1] = 0;
	    } else {
		jnint1_1.ntself[jnindx_(&c__1, &i__, &c__0) - 1] = 1;
	    }
/* L260: */
	}
	i__1 = jnint2_1.nl;
	for (il = 2; il <= i__1; ++il) {
	    io___283.ciunit = *nf;
	    s_rsle(&io___283);
	    e_rsle();
	    io___284.ciunit = *nf;
	    s_rsle(&io___284);
	    e_rsle();
	    i__2 = jnint2_1.m[il - 1];
	    for (j = 1; j <= i__2; ++j) {
		io___285.ciunit = *nf;
		s_rsle(&io___285);
		e_rsle();
		io___286.ciunit = *nf;
		s_rsfe(&io___286);
		i__3 = jnint2_1.m[il];
		for (i__ = 1; i__ <= i__3; ++i__) {
		    do_fio(&c__1, (char *)&jnint1_1.w[jnindx_(&il, &i__, &j) 
			    - 1], (ftnlen)sizeof(real));
		    do_fio(&c__1, line + (i__ - 1), (ftnlen)1);
		}
		e_rsfe();
		i__3 = jnint2_1.m[il];
		for (i__ = 1; i__ <= i__3; ++i__) {
		    if (*(unsigned char *)&line[i__ - 1] == '*') {
			jnint1_1.nself[jnindx_(&il, &i__, &j) - 1] = 0;
		    } else {
			jnint1_1.nself[jnindx_(&il, &i__, &j) - 1] = 1;
		    }
/* L320: */
		}
/* L310: */
	    }
	    io___287.ciunit = *nf;
	    s_rsle(&io___287);
	    e_rsle();
	    io___288.ciunit = *nf;
	    s_rsle(&io___288);
	    e_rsle();
	    io___289.ciunit = *nf;
	    s_rsle(&io___289);
	    e_rsle();
	    io___290.ciunit = *nf;
	    s_rsfe(&io___290);
	    i__2 = jnint2_1.m[il];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		do_fio(&c__1, (char *)&jnint1_1.t[jnindx_(&il, &i__, &c__0) - 
			1], (ftnlen)sizeof(real));
		do_fio(&c__1, line + (i__ - 1), (ftnlen)1);
	    }
	    e_rsfe();
	    i__2 = jnint2_1.m[il];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (*(unsigned char *)&line[i__ - 1] == '*') {
		    jnint1_1.ntself[jnindx_(&il, &i__, &c__0) - 1] = 0;
		} else {
		    jnint1_1.ntself[jnindx_(&il, &i__, &c__0) - 1] = 1;
		}
/* L330: */
	    }
/* L300: */
	}
    }
/* ...Write statistics on output file */
    if (jndat1_1.mstjn[5] < 0) {
	return 0;
    }
    jnhead_();
    jnstat_(&c__1);
    io___291.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___291);
    e_wsfe();
    return 0;
/* **** END OF JNREAD **************************************************** */
} /* jnread_ */

/* *********************************************************************** */
/* Subroutine */ int jnrold_(integer *nf)
{
    /* Format strings */
    static char fmt_690[] = "(a)";
    static char fmt_600[] = "(tr11,10i7)";
    static char fmt_610[] = "(tr11,10f7.4)";
    static char fmt_620[] = "(10(f7.4,a1))";
    static char fmt_630[] = "(10f8.4)";
    static char fmt_640[] = "(17x,\002Weights read from file produced with v"
	    "ersion 1\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_rsue(cilist *), do_uio(integer *, char *, ftnlen), e_rsue(void),
	     s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void)
	    , s_cmp(char *, char *, ftnlen, ftnlen), s_rsle(cilist *), e_rsle(
	    void), s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static char line[80];
    static real parj5;
    static integer mstj9, i__, j;
    static real parj11, parj12, parj13, parj14, parj15, parj16, parj17, 
	    parj18, parj19, parj20;
    extern /* Subroutine */ int jnerr_(integer *);
    static integer jf, il;
    extern /* Subroutine */ int jnhead_(void), jnsepa_(void);
    static integer nfsave;
    static real parjno[20];
    extern integer jnindx_(integer *, integer *, integer *);
    extern /* Subroutine */ int jnstat_(integer *);
    static integer mstjno[20];

    /* Fortran I/O blocks */
    static cilist io___305 = { 0, 0, 0, 0, 0 };
    static cilist io___309 = { 0, 0, 0, 0, 0 };
    static cilist io___310 = { 0, 0, 0, 0, 0 };
    static cilist io___311 = { 0, 0, 0, 0, 0 };
    static cilist io___312 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___314 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___316 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___317 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___318 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___319 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___320 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___321 = { 0, 0, 0, 0, 0 };
    static cilist io___322 = { 0, 0, 0, 0, 0 };
    static cilist io___323 = { 0, 0, 0, 0, 0 };
    static cilist io___325 = { 0, 0, 0, 0, 0 };
    static cilist io___326 = { 0, 0, 0, 0, 0 };
    static cilist io___328 = { 0, 0, 0, 0, 0 };
    static cilist io___329 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___330 = { 0, 0, 0, 0, 0 };
    static cilist io___331 = { 0, 0, 0, 0, 0 };
    static cilist io___332 = { 0, 0, 0, 0, 0 };
    static cilist io___333 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___334 = { 0, 0, 0, fmt_640, 0 };


/* ...JetNet subroutine Read weights from OLD versions */
/* ...(JETNET 1.0 and 1.1) */
/* ...Reads weights, thresholds and other statistics from a file NF and */
/* ...initializes the net */
/* ...MSTJN(9) has a new meaning: */
    mstj9 = jndat1_1.mstjn[8];
/* ...PARJN(5) has new meaning: */
    parj5 = jndat1_1.parjn[4];
/* ...PARJN(11-20) have new meanings: */
    parj11 = jndat1_1.parjn[10];
    parj12 = jndat1_1.parjn[11];
    parj13 = jndat1_1.parjn[12];
    parj14 = jndat1_1.parjn[13];
    parj15 = jndat1_1.parjn[14];
    parj16 = jndat1_1.parjn[15];
    parj17 = jndat1_1.parjn[16];
    parj18 = jndat1_1.parjn[17];
    parj19 = jndat1_1.parjn[18];
    parj20 = jndat1_1.parjn[19];
    if (*nf < 0) {
/* ...unformatted read */
	jf = -(*nf);
	io___305.ciunit = jf;
	s_rsue(&io___305);
	do_uio(&c__20, (char *)&mstjno[0], (ftnlen)sizeof(integer));
	do_uio(&c__20, (char *)&parjno[0], (ftnlen)sizeof(real));
	do_uio(&c__10, (char *)&jndat2_1.tinv[0], (ftnlen)sizeof(real));
	do_uio(&c__10, (char *)&jndat2_1.igfn[0], (ftnlen)sizeof(integer));
	e_rsue();
	for (i__ = 1; i__ <= 20; ++i__) {
	    jndat1_1.mstjn[i__ - 1] = mstjno[i__ - 1];
	    jndat1_1.parjn[i__ - 1] = parjno[i__ - 1];
/* L100: */
	}
	jndat1_1.mstjn[8] = mstj9;
	jndat1_1.parjn[4] = parj5;
	jndat1_1.parjn[10] = parj11;
	jndat1_1.parjn[11] = parj12;
	jndat1_1.parjn[12] = parj13;
	jndat1_1.parjn[13] = parj14;
	jndat1_1.parjn[14] = parj15;
	jndat1_1.parjn[15] = parj16;
	jndat1_1.parjn[16] = parj17;
	jndat1_1.parjn[17] = parj18;
	jndat1_1.parjn[18] = parj19;
	jndat1_1.parjn[19] = parj20;
	jnsepa_();
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___309.ciunit = jf;
	    s_rsue(&io___309);
	    do_uio(&c__1, (char *)&jnint1_1.w[i__ - 1], (ftnlen)sizeof(real));
	    e_rsue();
/* L110: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___310.ciunit = jf;
	    s_rsue(&io___310);
	    do_uio(&c__1, (char *)&jnint1_1.t[i__ - 1], (ftnlen)sizeof(real));
	    e_rsue();
/* L120: */
	}
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___311.ciunit = jf;
	    s_rsue(&io___311);
	    do_uio(&c__1, (char *)&jnint1_1.nself[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_rsue();
/* L130: */
	}
    } else {
/* ...Formatted dump */
	io___312.ciunit = *nf;
	s_rsfe(&io___312);
	do_fio(&c__1, line, (ftnlen)80);
	e_rsfe();
	if (s_cmp(line + 26, " D", (ftnlen)2, (ftnlen)2) == 0) {
	    jnerr_(&c__18);
	}
L900:
	io___314.ciunit = *nf;
	s_rsfe(&io___314);
	do_fio(&c__1, line, (ftnlen)80);
	e_rsfe();
	if (s_cmp(line, "         I      1      2", (ftnlen)24, (ftnlen)24) !=
		 0 && s_cmp(line, "          I      1      2", (ftnlen)25, (
		ftnlen)25) != 0) {
	    goto L900;
	}
	nfsave = jndat1_1.mstjn[5];
	io___316.ciunit = *nf;
	s_rsfe(&io___316);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&mstjno[i__ - 1], (ftnlen)sizeof(integer));
	}
	e_rsfe();
	io___317.ciunit = *nf;
	s_rsfe(&io___317);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&mstjno[i__ + 9], (ftnlen)sizeof(integer));
	}
	e_rsfe();
	io___318.ciunit = *nf;
	s_rsfe(&io___318);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&parjno[i__ - 1], (ftnlen)sizeof(real));
	}
	e_rsfe();
	io___319.ciunit = *nf;
	s_rsfe(&io___319);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.igfn[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	e_rsfe();
	io___320.ciunit = *nf;
	s_rsfe(&io___320);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.tinv[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	io___321.ciunit = *nf;
	s_rsle(&io___321);
	e_rsle();
	for (i__ = 1; i__ <= 20; ++i__) {
	    jndat1_1.mstjn[i__ - 1] = mstjno[i__ - 1];
	    jndat1_1.parjn[i__ - 1] = parjno[i__ - 1];
/* L200: */
	}
	jndat1_1.mstjn[5] = nfsave;
	jndat1_1.mstjn[8] = mstj9;
	jndat1_1.parjn[4] = parj5;
	jndat1_1.parjn[10] = parj11;
	jndat1_1.parjn[11] = parj12;
	jndat1_1.parjn[12] = parj13;
	jndat1_1.parjn[13] = parj14;
	jndat1_1.parjn[14] = parj15;
	jndat1_1.parjn[15] = parj16;
	jndat1_1.parjn[16] = parj17;
	jndat1_1.parjn[17] = parj18;
	jndat1_1.parjn[18] = parj19;
	jndat1_1.parjn[19] = parj20;
	jnsepa_();
	io___322.ciunit = *nf;
	s_rsle(&io___322);
	e_rsle();
	io___323.ciunit = *nf;
	s_rsle(&io___323);
	e_rsle();
	i__1 = jnint2_1.nl;
	for (il = 1; il <= i__1; ++il) {
	    io___325.ciunit = *nf;
	    s_rsle(&io___325);
	    e_rsle();
	    io___326.ciunit = *nf;
	    s_rsle(&io___326);
	    e_rsle();
	    i__2 = jnint2_1.m[il - 1];
	    for (j = 1; j <= i__2; ++j) {
		io___328.ciunit = *nf;
		s_rsle(&io___328);
		e_rsle();
		io___329.ciunit = *nf;
		s_rsfe(&io___329);
		i__3 = jnint2_1.m[il];
		for (i__ = 1; i__ <= i__3; ++i__) {
		    do_fio(&c__1, (char *)&jnint1_1.w[jnindx_(&il, &i__, &j) 
			    - 1], (ftnlen)sizeof(real));
		    do_fio(&c__1, line + (i__ - 1), (ftnlen)1);
		}
		e_rsfe();
		i__3 = jnint2_1.m[il];
		for (i__ = 1; i__ <= i__3; ++i__) {
		    if (*(unsigned char *)&line[i__ - 1] == '*') {
			jnint1_1.nself[jnindx_(&il, &i__, &j) - 1] = 0;
		    } else {
			jnint1_1.nself[jnindx_(&il, &i__, &j) - 1] = 1;
		    }
/* L230: */
		}
/* L220: */
	    }
	    io___330.ciunit = *nf;
	    s_rsle(&io___330);
	    e_rsle();
	    io___331.ciunit = *nf;
	    s_rsle(&io___331);
	    e_rsle();
	    io___332.ciunit = *nf;
	    s_rsle(&io___332);
	    e_rsle();
	    io___333.ciunit = *nf;
	    s_rsfe(&io___333);
	    i__2 = jnint2_1.m[il];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		do_fio(&c__1, (char *)&jnint1_1.t[jnindx_(&il, &i__, &c__0) - 
			1], (ftnlen)sizeof(real));
	    }
	    e_rsfe();
/* L210: */
	}
    }
/* ...Write statistics on output file */
    if (jndat1_1.mstjn[5] < 0) {
	return 0;
    }
    jnhead_();
    jnstat_(&c__1);
    io___334.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___334);
    e_wsfe();
    return 0;
/* **** END OF JNROLD **************************************************** */
} /* jnrold_ */

/* *********************************************************************** */
/* Subroutine */ int jnsatm_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer i__, il, mi;
    static real sum;

/* ...JetNet subroutine SATuration Measure */
/* ...Calculates the saturation measure "S" for each layer. */
/* ...Note: The response function for the layer must be a sigmoid. */
    i__1 = jnint2_1.nl;
    for (il = 1; il <= i__1; ++il) {
	if ((i__2 = jnint2_1.ng[il - 1], abs(i__2)) == 1 || (i__3 = 
		jnint2_1.ng[il - 1], abs(i__3)) == 5) {
	    sum = 0.f;
	    i__2 = jnint2_1.m[il];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		mi = jnint2_1.mv0[il - 1] + i__;
/* Computing 2nd power */
		r__1 = 1.f - jnint1_1.o[mi - 1] * 2.f;
		sum += r__1 * r__1;
/* L110: */
	    }
	    jnint2_1.sm[il - 1] += sum / (real) jnint2_1.m[il];
	} else if ((i__2 = jnint2_1.ng[il - 1], abs(i__2)) == 2) {
	    sum = 0.f;
	    i__2 = jnint2_1.m[il];
	    for (i__ = 1; i__ <= i__2; ++i__) {
		mi = jnint2_1.mv0[il - 1] + i__;
/* Computing 2nd power */
		r__1 = jnint1_1.o[mi - 1];
		sum += r__1 * r__1;
/* L120: */
	    }
	    jnint2_1.sm[il - 1] += sum / (real) jnint2_1.m[il];
	} else {
	    jnint2_1.sm[il - 1] = 0.f;
	}
/* L100: */
    }
    return 0;
/* **** END OF JNSATM **************************************************** */
} /* jnsatm_ */

/* *********************************************************************** */
/* Subroutine */ int jnscgr_(void)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real beta;
    static integer i__;
#define alpha ((real *)&jnint4_1 + 10)
    static real betak;
#define delta ((real *)&jnint4_1 + 6)
#define sigma ((real *)&jnint4_1 + 5)
    extern /* Subroutine */ int jnerr_(integer *);
#define lambda ((real *)&jnint4_1 + 11)
    static integer il;
    extern /* Subroutine */ int jncgbe_(real *, integer *);
#define cdelta ((real *)&jnint4_1 + 7)
    static real factor;

/* ...JetNet subroutine Scaled Conjugate GRadient */
/* ...Performs the Scaled Conjugate Gradient updating. */
/* ...The algorithm is described in: */
/* ...M. F. Moller, "A Scaled Conjugate Gradient Algorithm for Fast */
/* ...Supervised Learning", Neural Networks, Vol. 6, pp 525-533 (1993) */
/* ...The following notation is used (cf. Moller's article): */
/* ...S-vector = -(DW,DT) */
/* ...R-vector = (ODW,ODT) */
/* ...P-vector = G */
/* ...K = NSC */
/* ...MU=-DERRLN */
/* ...ZEPS=Machine precision */
    if (jndat1_1.mstjn[4] == 14) {
/* ...Move to last minimum and terminate the search. */
	if (jndat1_1.parjn[7] <= jnint4_1.errln[0]) {
	    jnint4_1.stepln[0] = 0.f;
	}
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.w[i__ - 1] -= jnint4_1.stepln[0] * jnint1_1.g[i__ - 1];
/* L400: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.t[i__ - 1] -= jnint4_1.stepln[0] * jnint1_1.g[i__ + 
		    jnint2_1.mm0[jnint2_1.nl] - 1];
/* L410: */
	}
	jnint4_1.nsc = 0;
	jnint4_1.nit = 0;
	jnint4_1.nc = 0;
	jnint4_1.ilinon = 0;
	jndat1_1.mstjn[4] = 9;
	jndat1_1.mstjn[36] = abs(jnint4_1.ilinon);
	return 0;
    }
    if (jnint4_1.ieval == 1) {
	goto L10;
    }
    if (jnint4_1.isucc > 0) {
/* ...Calculate 2nd order information: */
	++jnint4_1.nsc;
	jnint4_1.nc = 0;
	if (jnint4_1.icurve == 0) {
/* ...1st sweep -> Create new search direction and */
/* ...             Get curvature information */
	    jnint4_1.errln[0] = jndat1_1.parjn[7];
	    i__1 = (jnint4_1.nsc - 1) % (jnint2_1.mm0[jnint2_1.nl] + 
		    jnint2_1.mv0[jnint2_1.nl]);
	    jncgbe_(&betak, &i__1);
	    jnint4_1.derrln = 0.f;
	    beta = 1.f;
	    jnint4_1.gvec2 = 0.f;
	    for (il = jnint2_1.nl; il >= 1; --il) {
/* ...set effective beta in layer IL: */
		if (jndat2_1.tinv[il - 1] == 0.f) {
		    beta *= jndat1_1.parjn[2];
		} else {
		    beta *= (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
		}
		i__1 = jnint2_1.mm0[il];
		for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__1; ++i__) {
		    jnint1_1.g[i__ - 1] = betak * jnint1_1.g[i__ - 1] + 
			    jnint1_1.odw[i__ - 1] * (real) jnint1_1.nself[i__ 
			    - 1] * beta;
		    jnint4_1.derrln -= jnint1_1.odw[i__ - 1] * (real) 
			    jnint1_1.nself[i__ - 1] * beta * jnint1_1.g[i__ - 
			    1];
/* Computing 2nd power */
		    r__1 = jnint1_1.g[i__ - 1];
		    jnint4_1.gvec2 += r__1 * r__1;
/* L110: */
		}
		i__1 = jnint2_1.mv0[il];
		for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__1; ++i__) {
		    jnint1_1.g[i__ + jnint2_1.mm0[jnint2_1.nl] - 1] = betak * 
			    jnint1_1.g[i__ + jnint2_1.mm0[jnint2_1.nl] - 1] + 
			    jnint1_1.odt[i__ - 1] * (real) jnint1_1.ntself[
			    i__ - 1] * beta;
		    jnint4_1.derrln -= jnint1_1.odt[i__ - 1] * (real) 
			    jnint1_1.ntself[i__ - 1] * beta * jnint1_1.g[i__ 
			    + jnint2_1.mm0[jnint2_1.nl] - 1];
/* Computing 2nd power */
		    r__1 = jnint1_1.g[i__ + jnint2_1.mm0[jnint2_1.nl] - 1];
		    jnint4_1.gvec2 += r__1 * r__1;
/* L120: */
		}
/* L100: */
	    }
/* ...Initial value for lambda */
	    if (jnint4_1.nsc == 1) {
		*lambda = jndat1_1.parjn[28];
	    }
	    jnint4_1.nit = 1;
	    *sigma = jndat1_1.parjn[27] / (sqrt(jnint4_1.gvec2) + 1e-20f);
	    factor = (real) jndat1_1.mstjn[1];
	    i__1 = jnint2_1.mm0[jnint2_1.nl];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jnint1_1.dw[i__ - 1] = -jnint1_1.dw[i__ - 1] * factor;
		jnint1_1.w[i__ - 1] += *sigma * jnint1_1.g[i__ - 1];
/* L200: */
	    }
	    i__1 = jnint2_1.mv0[jnint2_1.nl];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jnint1_1.dt[i__ - 1] = -jnint1_1.dt[i__ - 1] * factor;
		jnint1_1.t[i__ - 1] += *sigma * jnint1_1.g[i__ + jnint2_1.mm0[
			jnint2_1.nl] - 1];
/* L210: */
	    }
	    jnint4_1.icurve = 1;
	    jndat1_1.mstjn[36] = abs(jnint4_1.ilinon);
	    jnint4_1.stepln[0] = *sigma;
	    return 0;
	} else {
/* ...2nd sweep -> Curvature information exists */
	    *delta = 0.f;
	    factor = *sigma * (real) jndat1_1.mstjn[1] + 1e-20f;
	    i__1 = jnint2_1.mm0[jnint2_1.nl];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jnint1_1.dw[i__ - 1] = -jnint1_1.dw[i__ - 1] / factor;
		*delta += jnint1_1.g[i__ - 1] * jnint1_1.dw[i__ - 1];
		jnint1_1.w[i__ - 1] -= *sigma * jnint1_1.g[i__ - 1];
/* L220: */
	    }
	    i__1 = jnint2_1.mv0[jnint2_1.nl];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jnint1_1.dt[i__ - 1] = -jnint1_1.dt[i__ - 1] / factor;
		*delta += jnint1_1.g[i__ + jnint2_1.mm0[jnint2_1.nl] - 1] * 
			jnint1_1.dt[i__ - 1];
		jnint1_1.t[i__ - 1] -= *sigma * jnint1_1.g[i__ + jnint2_1.mm0[
			jnint2_1.nl] - 1];
/* L230: */
	    }
	    jnint4_1.ilinon = 1;
	    jnint4_1.icurve = 0;
	    jnint4_1.stepln[0] = 0.f;
	}
    }
    if (*delta + *lambda * jnint4_1.gvec2 <= 0.f) {
/* ...Make Hessian positive definite: */
	*lambda = (*lambda - *delta / (jnint4_1.gvec2 + 1e-20f)) * 2.f;
    }
    *delta += *lambda * jnint4_1.gvec2;
/* ...Update weights to calculate comparison parameter: */
    *alpha = -jnint4_1.derrln / (*delta + 1e-20f);
    if (dabs(*alpha) <= 1e-8f || jnint4_1.nit >= jndat1_1.mstjn[34] || *
	    lambda > 1.f) {
/* ...Search is stuck! -> Restart. */
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.dw[i__ - 1] = 0.f;
/* L280: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.dt[i__ - 1] = 0.f;
/* L290: */
	}
	jnint4_1.stepln[0] = 0.f;
	jnint4_1.isucc = 1;
	jnint4_1.ieval = 0;
	jnint4_1.ilinon = 0;
	jnint4_1.nsc = 0;
	++jndat1_1.mstjn[37];
	if (jndat1_1.mstjn[37] > jndat1_1.mstjn[35]) {
	    jnerr_(&c__21);
	}
	jndat1_1.mstjn[36] = 0;
	return 0;
    }
    i__1 = jnint2_1.mm0[jnint2_1.nl];
    for (i__ = 1; i__ <= i__1; ++i__) {
	jnint1_1.w[i__ - 1] += *alpha * jnint1_1.g[i__ - 1];
/* L300: */
    }
    i__1 = jnint2_1.mv0[jnint2_1.nl];
    for (i__ = 1; i__ <= i__1; ++i__) {
	jnint1_1.t[i__ - 1] += *alpha * jnint1_1.g[i__ + jnint2_1.mm0[
		jnint2_1.nl] - 1];
/* L310: */
    }
    jnint4_1.stepln[0] = *alpha;
    jnint4_1.ieval = 1;
    ++jnint4_1.nit;
    jndat1_1.mstjn[36] = abs(jnint4_1.ilinon);
    return 0;
/* ...Come here if the comparison parameter is to be calculated: */
L10:
/* Computing 2nd power */
    r__1 = jnint4_1.derrln;
    *cdelta = *delta * 2.f * (jnint4_1.errln[0] - jndat1_1.parjn[7]) / (r__1 *
	     r__1 + 1e-20f);
    jnint4_1.ieval = 0;
    if (*cdelta >= 0.f) {
/* ...Successful reduction in error. */
	jnint4_1.isucc = 1;
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.dw[i__ - 1] = 0.f;
/* L320: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.dt[i__ - 1] = 0.f;
/* L330: */
	}
	jnint4_1.stepln[0] = 0.f;
	jnint4_1.ilinon = 0;
	if (*cdelta >= .75f) {
	    *lambda /= 4.f;
	}
    } else {
/* ...Not a successful error reduction -> move back and make new attempt */
	jnint4_1.isucc = 0;
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.w[i__ - 1] -= *alpha * jnint1_1.g[i__ - 1];
/* L340: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.t[i__ - 1] -= *alpha * jnint1_1.g[i__ + jnint2_1.mm0[
		    jnint2_1.nl] - 1];
/* L350: */
	}
	jnint4_1.stepln[0] = 0.f;
    }
    if (*cdelta < .25f) {
	*lambda += *delta * (1.f - *cdelta) / jnint4_1.gvec2;
    }
    jndat1_1.mstjn[36] = abs(jnint4_1.ilinon);
    return 0;
/* **** END OF JNSCGR **************************************************** */
} /* jnscgr_ */

#undef cdelta
#undef lambda
#undef sigma
#undef delta
#undef alpha


/* *********************************************************************** */
/* Subroutine */ int jnsefi_(integer *ila, integer *i1, integer *i2, integer *
	j1, integer *j2, integer *no)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer idum;
    static real width;
    extern /* Subroutine */ int jnerr_(integer *);
    static integer if__, jf, ii, jj, il, jl;
    extern integer jnindx_(integer *, integer *, integer *);
    extern doublereal rjn_(integer *);

/* ...JetNet subroutine SElect FIelds */
/* ...Switches the updating of the weights between nodes I1 to I2 in layer */
/* ...ILA and nodes J1 to J2 in layer ILA-1 on or off according to NO. If */
/* ...NO<=0 updating is turned off else it is turned on. In addition if */
/* ...NO=0 the weight is set to zero and if NO=1 the weight is */
/* ...reinitialized. */
    if (jndat1_1.mstjn[7] == 0) {
	jnerr_(&c__6);
    }
    if__ = max(*i1,1);
/* Computing MIN */
    i__1 = *i2, i__2 = jnint2_1.m[*ila];
    il = min(i__1,i__2);
    if (jndat2_1.widl[*ila - 1] <= 0.f) {
	width = jndat1_1.parjn[3];
    } else {
	width = jndat2_1.widl[*ila - 1];
    }
    if (*j1 != 0 || *j2 != 0) {
	jf = max(*j1,1);
/* Computing MIN */
	i__1 = *j2, i__2 = jnint2_1.m[*ila - 1];
	jl = min(i__1,i__2);
	i__1 = il;
	for (ii = if__; ii <= i__1; ++ii) {
	    i__2 = jl;
	    for (jj = jf; jj <= i__2; ++jj) {
		if (*no > 0) {
		    jnint1_1.nself[jnindx_(ila, &ii, &jj) - 1] = 1;
		} else {
		    jnint1_1.nself[jnindx_(ila, &ii, &jj) - 1] = 0;
		}
		if (*no == 1) {
		    idum = jj;
		    jnint1_1.w[jnindx_(ila, &ii, &jj) - 1] = (rjn_(&idum) * 
			    2.f - 1.f) * width;
		} else if (*no == 0) {
		    jnint1_1.w[jnindx_(ila, &ii, &jj) - 1] = 0.f;
		}
/* L110: */
	    }
/* L100: */
	}
    } else {
	i__1 = il;
	for (ii = if__; ii <= i__1; ++ii) {
	    if (*no > 0) {
		jnint1_1.ntself[jnindx_(ila, &ii, &c__0) - 1] = 1;
	    } else {
		jnint1_1.ntself[jnindx_(ila, &ii, &c__0) - 1] = 0;
	    }
	    if (*no == 1) {
		idum = ii;
		jnint1_1.t[jnindx_(ila, &ii, &c__0) - 1] = (rjn_(&idum) * 2.f 
			- 1.f) * width;
	    } else if (*no == 0) {
		jnint1_1.t[jnindx_(ila, &ii, &c__0) - 1] = 0.f;
	    }
/* L200: */
	}
    }
    return 0;
/* **** END OF JNSEFI **************************************************** */
} /* jnsefi_ */

/* *********************************************************************** */
/* Subroutine */ int jnsepa_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int jnerr_(integer *);
    static integer nrflw, il, it, iw;
    extern integer jnindx_(integer *, integer *, integer *);
    static real eta;

/* ...JetNet subroutine SEt PArameters */
/* ...Sets parameters in /JNINT2/ */
/* ...last layer: */
    jnint2_1.nl = jndat1_1.mstjn[0] - 1;
    if (jnint2_1.nl > 10 || jnint2_1.nl < 1) {
	jnerr_(&c__1);
    }
/* ...number of nodes in each layer */
    if (jndat1_1.mstjn[9] > 1000) {
	jnerr_(&c__7);
    }
    if (jndat1_1.mstjn[jnint2_1.nl + 9] > 1000) {
	jnerr_(&c__8);
    }
    i__1 = jndat1_1.mstjn[0];
    for (il = 1; il <= i__1; ++il) {
	if (jndat1_1.mstjn[il + 8] == 0) {
	    jndat1_1.mstjn[6] = il;
	    jnerr_(&c__25);
	}
/* L400: */
    }
    if (jndat1_1.mstjn[22] != 0) {
/* ...receptive fields will be used, check consistency and set indexes */
	if (jndat1_1.mstjn[23] == 0 || jndat1_1.mstjn[24] <= 0 || 
		jndat1_1.mstjn[25] <= 0) {
	    jnerr_(&c__10);
	}
	jnint3_1.nxin = jndat1_1.mstjn[22];
	jnint3_1.nyin = jndat1_1.mstjn[23];
	jnint3_1.nxrf = jndat1_1.mstjn[24];
	jnint3_1.nyrf = jndat1_1.mstjn[25];
	jnint3_1.nhprf = abs(jndat1_1.mstjn[26]);
	if (jndat1_1.mstjn[9] < (i__1 = jnint3_1.nxin * jnint3_1.nyin, abs(
		i__1))) {
	    jnerr_(&c__11);
	}
	if (jnint3_1.nxrf > abs(jnint3_1.nxin) || jnint3_1.nyrf > abs(
		jnint3_1.nyin)) {
	    jnerr_(&c__11);
	}
	if (jnint3_1.nxin > 0) {
	    jnint3_1.nxhrf = jnint3_1.nxin - jnint3_1.nxrf + 1;
	} else {
	    jnint3_1.nxhrf = -jnint3_1.nxin;
	}
	if (jnint3_1.nyin > 0) {
	    jnint3_1.nyhrf = jnint3_1.nyin - jnint3_1.nyrf + 1;
	} else {
	    jnint3_1.nyhrf = -jnint3_1.nyin;
	}
	jnint3_1.nhrf = jnint3_1.nxhrf * jnint3_1.nyhrf;
/* Computing MAX */
	i__1 = jndat1_1.mstjn[10], i__2 = jnint3_1.nhrf * jnint3_1.nhprf;
	jndat1_1.mstjn[10] = max(i__1,i__2);
	jnint3_1.nrfw = jnint3_1.nxrf * jnint3_1.nyrf + jndat1_1.mstjn[9] - (
		i__1 = jnint3_1.nxin * jnint3_1.nyin, abs(i__1));
	nrflw = jnint3_1.nrfw * jnint3_1.nhprf + (jndat1_1.mstjn[10] - 
		jnint3_1.nhrf * jnint3_1.nhprf) * jndat1_1.mstjn[9];
    } else {
	jnint3_1.nxin = 0;
	jnint3_1.nyin = 0;
	jnint3_1.nxrf = 0;
	jnint3_1.nyrf = 0;
	jnint3_1.nxhrf = 0;
	jnint3_1.nyhrf = 0;
	jnint3_1.nhrf = 0;
	jnint3_1.nrfw = 0;
	nrflw = 0;
    }
    i__1 = jnint2_1.nl;
    for (il = 0; il <= i__1; ++il) {
	jnint2_1.m[il] = jndat1_1.mstjn[il + 9];
/* L100: */
    }
/* ...offset index in node vectors and weight vectors */
    jnint2_1.mv0[0] = 0;
    jnint2_1.mm0[0] = 0;
    jnint2_1.mv0[1] = jnint2_1.m[1];
    jnint2_1.mm0[1] = jnint2_1.m[0] * jnint2_1.m[1];
    if (jnint3_1.nxin != 0) {
	jnint2_1.mm0[1] = nrflw;
    }
    i__1 = jnint2_1.nl + 1;
    for (il = 3; il <= i__1; ++il) {
	jnint2_1.mv0[il - 1] = jnint2_1.mv0[il - 2] + jnint2_1.m[il - 1];
	jnint2_1.mm0[il - 1] = jnint2_1.mm0[il - 2] + jnint2_1.m[il - 1] * 
		jnint2_1.m[il - 2];
/* L110: */
    }
    if (jnint2_1.mv0[jnint2_1.nl] > 2000) {
	jnerr_(&c__2);
    }
    if (jnint2_1.mm0[jnint2_1.nl] > 150000) {
	jnerr_(&c__3);
    }
/* ...check Potts-nodes */
    jnint2_1.ipott = jndat1_1.mstjn[3];
    if (jnint2_1.ipott >= 2) {
	if (jnint2_1.m[jnint2_1.nl] % jnint2_1.ipott != 0) {
	    jnerr_(&c__4);
	}
	jndat2_1.igfn[jnint2_1.nl - 1] = 3;
    }
    if (jnint2_1.ipott == 1) {
	jndat2_1.igfn[jnint2_1.nl - 1] = 5;
    }
/* ...set transfer functions to use */
    i__1 = jnint2_1.nl;
    for (il = 1; il <= i__1; ++il) {
	if (jndat2_1.igfn[il - 1] == 0) {
	    jnint2_1.ng[il - 1] = jndat1_1.mstjn[2];
	} else {
	    jnint2_1.ng[il - 1] = jndat2_1.igfn[il - 1];
	}
/* L120: */
    }
/* ...Check consistency between error measure and transfer function. */
    if (jndat1_1.mstjn[3] == 1) {
	if (jndat2_1.igfn[jnint2_1.nl - 1] == 2 || jndat2_1.igfn[jnint2_1.nl 
		- 1] == 4) {
	    jnerr_(&c__31);
	}
    }
/* ...Zero weight and threshold vectors */
    i__1 = jnint2_1.mv0[jnint2_1.nl];
    for (i__ = 1; i__ <= i__1; ++i__) {
	jnint1_1.dt[i__ - 1] = 0.f;
	jnint1_1.ntself[i__ - 1] = 1;
/* L200: */
    }
    i__1 = jnint2_1.mm0[jnint2_1.nl];
    for (i__ = 1; i__ <= i__1; ++i__) {
	jnint1_1.dw[i__ - 1] = 0.f;
	jnint1_1.nself[i__ - 1] = 1;
/* L210: */
    }
/* ...set precision chopping */
    if (jndat1_1.mstjn[27] > 0 || jndat1_1.mstjn[28] > 0 || jndat1_1.mstjn[29]
	     > 0) {
	jnint2_1.icpon = 1;
    }
/* ...If updating is turned off, stop here. */
    if (jndat1_1.mstjn[4] == 9) {
	jnerr_(&c__32);
    }
/* ...Initialize Quickprop, Rprop and Conjugate Gradient searches */
    if (jndat1_1.mstjn[4] >= 3 && jndat1_1.mstjn[4] <= 14) {
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.odw[i__ - 1] = 0.f;
	    jnint1_1.g[i__ - 1] = 0.f;
/* L300: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jnint1_1.odt[i__ - 1] = 0.f;
	    jnint1_1.g[jnint2_1.mm0[jnint2_1.nl] + i__ - 1] = 0.f;
/* L310: */
	}
    }
    if (jndat1_1.mstjn[4] == 8) {
	jnerr_(&c__19);
    }
    jndat1_1.mstjn[7] = 1;
/* ...Initialize Rprop learning rate: */
    for (il = jnint2_1.nl; il >= 1; --il) {
	if (jndat2_1.etal[il - 1] == 0.f) {
	    eta = jndat1_1.parjn[0] / (real) jndat1_1.mstjn[1];
	} else {
	    eta = jndat2_1.etal[il - 1] / (real) jndat1_1.mstjn[1];
	}
	i__1 = jnint2_1.m[il];
	for (i__ = 1; i__ <= i__1; ++i__) {
	    it = jnindx_(&il, &i__, &c__0);
	    jnint1_1.etav[jnint2_1.mm0[jnint2_1.nl] + it - 1] = eta;
	    i__2 = jnint2_1.m[il - 1];
	    for (j = 1; j <= i__2; ++j) {
		iw = jnindx_(&il, &i__, &j);
		jnint1_1.etav[iw - 1] = eta;
/* L520: */
	    }
/* L510: */
	}
/* L500: */
    }
/* ...Reset restart counter */
    jndat1_1.mstjn[37] = 0;
    return 0;
/* **** END OF JNSEPA **************************************************** */
} /* jnsepa_ */

/* *********************************************************************** */
/* Subroutine */ int jnstat_(integer *is)
{
    /* Format strings */
    static char fmt_600[] = "(22x,\002Initialized for a\002,i2,\002 layered "
	    "net with\002)";
    static char fmt_615[] = "(22x,i3,\002 nodes in layer number\002,i2,\002 "
	    "(output layer)\002)";
    static char fmt_620[] = "(22x,i3,\002 nodes in layer number\002,i2)";
    static char fmt_610[] = "(22x,i3,\002 nodes in layer number 0 (input lay"
	    "er)\002)";
    static char fmt_630[] = "(5x,\002with\002,i3,\002-dimensional Potts node"
	    "s in output layer.\002)";
    static char fmt_631[] = "(5x,\002receptive fields in first layer assumin"
	    "g the \002,i4,\002 first nodes in the\002,/,5x,\002input layer a"
	    "re organised in a plane of \002,i3,\002*\002,i3,\002 nodes, wher"
	    "e the\002,/,5x,\002receptive field nodes scan \002,i3,\002*\002,"
	    "i3,\002 input nodes each with \002,i3,\002 hidden\002/,5x,\002no"
	    "des per field.\002)";
    static char fmt_632[] = "(5x,\002The input layer is assumed to be cycl"
	    "ic \002,\002in the x-direction.\002)";
    static char fmt_633[] = "(5x,\002The input layer is assumed to be cycl"
	    "ic \002,\002in the y-direction.\002)";
    static char fmt_634[] = "(5x,\002The weights from equivalent nodes with "
	    "receptive \002,\002fields are clamped.\002)";
    static char fmt_635[] = "(22x,\002Using Cross-Entropy error.\002)";
    static char fmt_639[] = "(22x,\002Standard Back-Propagation updating."
	    "\002)";
    static char fmt_640[] = "(22x,\002Manhattan updating.\002)";
    static char fmt_645[] = "(22x,\002Langevin updating.\002)";
    static char fmt_646[] = "(22x,\002Quickprop updating.\002)";
    static char fmt_647[] = "(22x,\002Conjugate Gradient updating (Polak-Rib"
	    "iere).\002)";
    static char fmt_648[] = "(22x,\002Conjugate Gradient updating (Hestenes-"
	    "Stiefel).\002)";
    static char fmt_649[] = "(22x,\002Conjugate Gradient updating (Fletcher-"
	    "Reeves).\002)";
    static char fmt_651[] = "(22x,\002Conjugate Gradient updating (Shanno)"
	    ".\002)";
    static char fmt_652[] = "(22x,\002Scaled Conj. Grad. updating (Polak-Rib"
	    "iere).\002)";
    static char fmt_653[] = "(22x,\002Scaled Conj. Grad. updating (Hestenes-"
	    "Stiefel).\002)";
    static char fmt_654[] = "(22x,\002Scaled Conj. Grad. updating (Fletcher-"
	    "Reeves).\002)";
    static char fmt_655[] = "(22x,\002Scaled Conj. Grad. updating (Shanno)"
	    ".\002)";
    static char fmt_656[] = "(22x,\002Rprop updating.\002)";
    static char fmt_650[] = "(18x,\002Values of parameters and switches in J"
	    "ETNET\002)";
    static char fmt_660[] = "(a10,10i7)";
    static char fmt_661[] = "(a10,6i7,f7.3,3i7)";
    static char fmt_670[] = "(a10,10f7.4)";
    static char fmt_680[] = "(5x,\002Time factor for this net:\002,i10)";
    static char fmt_690[] = "(5x,\002Effective number of weights:\002,i7)";
    static char fmt_700[] = "(5x,\002The Hessian Matrix: (\002,i3,\002 x "
	    "\002,i3,\002)\002)";
    static char fmt_720[] = "(5x,\002Column: \002,i3)";
    static char fmt_710[] = "(5x,7(e9.2,1x))";
    static char fmt_730[] = "(5x,\002Diagonal elements only\002)";
    static char fmt_740[] = "(5x,\002Trace(H) = \002,f10.5)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3;

    /* Builtin functions */
    integer s_wsle(cilist *), e_wsle(void), s_wsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_wsfe(void);
    double r_lg10(real *);

    /* Local variables */
    static real par22;
    static integer jhop, nhop, i__, j, nwfac;
    static real trace;
    static integer nwgts, nwsum, il, iw, nextra;

    /* Fortran I/O blocks */
    static cilist io___364 = { 0, 0, 0, 0, 0 };
    static cilist io___365 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___366 = { 0, 0, 0, fmt_615, 0 };
    static cilist io___368 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___369 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___370 = { 0, 0, 0, 0, 0 };
    static cilist io___371 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___372 = { 0, 0, 0, fmt_631, 0 };
    static cilist io___373 = { 0, 0, 0, fmt_632, 0 };
    static cilist io___374 = { 0, 0, 0, fmt_633, 0 };
    static cilist io___375 = { 0, 0, 0, fmt_634, 0 };
    static cilist io___376 = { 0, 0, 0, fmt_635, 0 };
    static cilist io___377 = { 0, 0, 0, fmt_639, 0 };
    static cilist io___378 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___379 = { 0, 0, 0, fmt_645, 0 };
    static cilist io___380 = { 0, 0, 0, fmt_646, 0 };
    static cilist io___381 = { 0, 0, 0, fmt_647, 0 };
    static cilist io___382 = { 0, 0, 0, fmt_648, 0 };
    static cilist io___383 = { 0, 0, 0, fmt_649, 0 };
    static cilist io___384 = { 0, 0, 0, fmt_651, 0 };
    static cilist io___385 = { 0, 0, 0, fmt_652, 0 };
    static cilist io___386 = { 0, 0, 0, fmt_653, 0 };
    static cilist io___387 = { 0, 0, 0, fmt_654, 0 };
    static cilist io___388 = { 0, 0, 0, fmt_655, 0 };
    static cilist io___389 = { 0, 0, 0, fmt_656, 0 };
    static cilist io___390 = { 0, 0, 0, 0, 0 };
    static cilist io___392 = { 0, 0, 0, 0, 0 };
    static cilist io___393 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___394 = { 0, 0, 0, 0, 0 };
    static cilist io___395 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___397 = { 0, 0, 0, fmt_661, 0 };
    static cilist io___398 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___399 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___400 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___401 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___402 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___403 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___404 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___405 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___406 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___407 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___408 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___409 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___410 = { 0, 0, 0, 0, 0 };
    static cilist io___413 = { 0, 0, 0, fmt_680, 0 };
    static cilist io___414 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___416 = { 0, 0, 0, fmt_700, 0 };
    static cilist io___417 = { 0, 0, 0, 0, 0 };
    static cilist io___420 = { 0, 0, 0, fmt_720, 0 };
    static cilist io___422 = { 0, 0, 0, fmt_710, 0 };
    static cilist io___424 = { 0, 0, 0, fmt_710, 0 };
    static cilist io___425 = { 0, 0, 0, 0, 0 };
    static cilist io___428 = { 0, 0, 0, fmt_700, 0 };
    static cilist io___429 = { 0, 0, 0, fmt_730, 0 };
    static cilist io___430 = { 0, 0, 0, 0, 0 };
    static cilist io___431 = { 0, 0, 0, fmt_710, 0 };
    static cilist io___432 = { 0, 0, 0, fmt_710, 0 };
    static cilist io___433 = { 0, 0, 0, 0, 0 };
    static cilist io___434 = { 0, 0, 0, fmt_740, 0 };
    static cilist io___435 = { 0, 0, 0, 0, 0 };


/* ...JetNet subroutine output STATistics */
/* ...Statistics chosen by IS is written on the default file */
    if (*is == 1) {
/* ...Write out number of layers, units and receptive field status */
	io___364.ciunit = jndat1_1.mstjn[5];
	s_wsle(&io___364);
	e_wsle();
	io___365.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___365);
	i__1 = jnint2_1.nl + 1;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfe();
	io___366.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___366);
	do_fio(&c__1, (char *)&jnint2_1.m[jnint2_1.nl], (ftnlen)sizeof(
		integer));
	do_fio(&c__1, (char *)&jnint2_1.nl, (ftnlen)sizeof(integer));
	e_wsfe();
	for (il = jnint2_1.nl - 1; il >= 1; --il) {
	    io___368.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___368);
	    do_fio(&c__1, (char *)&jnint2_1.m[il], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&il, (ftnlen)sizeof(integer));
	    e_wsfe();
/* L100: */
	}
	io___369.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___369);
	do_fio(&c__1, (char *)&jnint2_1.m[0], (ftnlen)sizeof(integer));
	e_wsfe();
	io___370.ciunit = jndat1_1.mstjn[5];
	s_wsle(&io___370);
	e_wsle();
	if (jnint2_1.ipott > 1) {
	    io___371.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___371);
	    do_fio(&c__1, (char *)&jnint2_1.ipott, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	if (jndat1_1.mstjn[22] != 0) {
	    io___372.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___372);
	    i__2 = (i__1 = jndat1_1.mstjn[22] * jndat1_1.mstjn[23], abs(i__1))
		    ;
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	    i__3 = abs(jndat1_1.mstjn[22]);
	    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
	    i__4 = abs(jndat1_1.mstjn[23]);
	    do_fio(&c__1, (char *)&i__4, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[24], (ftnlen)sizeof(integer)
		    );
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[25], (ftnlen)sizeof(integer)
		    );
	    i__5 = abs(jndat1_1.mstjn[26]);
	    do_fio(&c__1, (char *)&i__5, (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	if (jndat1_1.mstjn[22] < 0) {
	    io___373.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___373);
	    e_wsfe();
	}
	if (jndat1_1.mstjn[23] < 0) {
	    io___374.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___374);
	    e_wsfe();
	}
	if (jndat1_1.mstjn[26] < 0) {
	    io___375.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___375);
	    e_wsfe();
	}
	if (jnint2_1.ipott == 1) {
	    io___376.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___376);
	    e_wsfe();
	}
	if (jndat1_1.mstjn[4] == 0) {
	    io___377.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___377);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 1) {
	    io___378.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___378);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 2) {
	    io___379.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___379);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 3) {
	    io___380.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___380);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 4) {
	    io___381.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___381);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 5) {
	    io___382.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___382);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 6) {
	    io___383.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___383);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 7) {
	    io___384.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___384);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 10) {
	    io___385.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___385);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 11) {
	    io___386.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___386);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 12) {
	    io___387.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___387);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 13) {
	    io___388.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___388);
	    e_wsfe();
	} else if (jndat1_1.mstjn[4] == 15) {
	    io___389.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___389);
	    e_wsfe();
	}
	io___390.ciunit = jndat1_1.mstjn[5];
	s_wsle(&io___390);
	e_wsle();
    } else if (*is == 2) {
/* ...Write out values of parameters and switches */
	par22 = jndat1_1.parjn[21];
	jndat1_1.parjn[21] = r_lg10(&par22);
	io___392.ciunit = jndat1_1.mstjn[5];
	s_wsle(&io___392);
	e_wsle();
	io___393.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___393);
	e_wsfe();
	io___394.ciunit = jndat1_1.mstjn[5];
	s_wsle(&io___394);
	e_wsle();
	io___395.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___395);
	do_fio(&c__1, "I ", (ftnlen)2);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	}
	e_wsfe();
	io___397.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___397);
	do_fio(&c__1, "MSTJN (I)", (ftnlen)9);
	for (i__ = 1; i__ <= 6; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
/* Computing MAX */
	r__3 = (real) jndat1_1.mstjn[6];
	r__2 = dmax(r__3,1.f);
	r__1 = r_lg10(&r__2);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	for (i__ = 8; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
	io___398.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___398);
	do_fio(&c__1, "(10+I)", (ftnlen)6);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ + 9], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
	io___399.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___399);
	do_fio(&c__1, "(20+I)", (ftnlen)6);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ + 19], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
	io___400.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___400);
	do_fio(&c__1, "(30+I)", (ftnlen)6);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjn[i__ + 29], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
	io___401.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___401);
	do_fio(&c__1, "PARJN (I)", (ftnlen)9);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjn[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_wsfe();
	io___402.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___402);
	do_fio(&c__1, "(10+I)", (ftnlen)6);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjn[i__ + 9], (ftnlen)sizeof(
		    real));
	}
	e_wsfe();
	io___403.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___403);
	do_fio(&c__1, "(20+I)", (ftnlen)6);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjn[i__ + 19], (ftnlen)sizeof(
		    real));
	}
	e_wsfe();
	io___404.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___404);
	do_fio(&c__1, "(30+I)", (ftnlen)6);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjn[i__ + 29], (ftnlen)sizeof(
		    real));
	}
	e_wsfe();
	io___405.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___405);
	do_fio(&c__1, "IGFN (I)", (ftnlen)8);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.igfn[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
	io___406.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___406);
	do_fio(&c__1, "TINV (I)", (ftnlen)8);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.tinv[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_wsfe();
	io___407.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___407);
	do_fio(&c__1, "ETAL (I)", (ftnlen)8);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.etal[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_wsfe();
	io___408.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___408);
	do_fio(&c__1, "WIDL (I)", (ftnlen)8);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.widl[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_wsfe();
	io___409.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___409);
	do_fio(&c__1, "SATM (I)", (ftnlen)8);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat2_1.satm[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_wsfe();
	io___410.ciunit = jndat1_1.mstjn[5];
	s_wsle(&io___410);
	e_wsle();
	jndat1_1.parjn[21] = par22;
    } else if (*is == 3) {
/* ...Write out time factor for net */
	nwfac = 0;
	if (jnint3_1.nxin == 0) {
	    nwfac = jnint2_1.mm0[jnint2_1.nl] + jnint2_1.mv0[jnint2_1.nl];
	    nwsum = nwfac;
	} else {
	    nwfac = jnint2_1.mm0[jnint2_1.nl] + jnint2_1.mv0[jnint2_1.nl] - 
		    jnint2_1.mm0[1] + jnint3_1.nhrf * jnint3_1.nrfw * 
		    jnint3_1.nhprf + (jndat1_1.mstjn[10] - jnint3_1.nhrf * 
		    jnint3_1.nhprf) * jndat1_1.mstjn[9];
	    nwsum = jnint2_1.mm0[jnint2_1.nl] + jnint2_1.mv0[jnint2_1.nl] - 
		    jnint2_1.mm0[1] + jnint3_1.nrfw * jnint3_1.nhprf - (
		    jnint3_1.nhrf - 1) * jnint3_1.nhprf + (jndat1_1.mstjn[10] 
		    - jnint3_1.nhrf * jnint3_1.nhprf) * jndat1_1.mstjn[9];
	    if (jndat1_1.mstjn[26] < 0) {
		nwsum -= (jnint3_1.nhrf - 1) * jnint3_1.nhprf * 
			jndat1_1.mstjn[11];
	    }
	}
	io___413.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___413);
	do_fio(&c__1, (char *)&nwfac, (ftnlen)sizeof(integer));
	e_wsfe();
	io___414.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___414);
	do_fio(&c__1, (char *)&nwsum, (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*is == 4) {
/* ...Write out Hessian Matrix */
	nwgts = jnint2_1.mm0[jnint2_1.nl] + jnint2_1.mv0[jnint2_1.nl];
	io___416.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___416);
	do_fio(&c__1, (char *)&nwgts, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nwgts, (ftnlen)sizeof(integer));
	e_wsfe();
	io___417.ciunit = jndat1_1.mstjn[5];
	s_wsle(&io___417);
	e_wsle();
	nhop = nwgts / 7;
	nextra = nwgts - nhop * 7;
	i__1 = nwgts;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io___420.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___420);
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	    e_wsfe();
	    i__2 = nhop;
	    for (jhop = 1; jhop <= i__2; ++jhop) {
		io___422.ciunit = jndat1_1.mstjn[5];
		s_wsfe(&io___422);
		for (j = 1; j <= 7; ++j) {
		    do_fio(&c__1, (char *)&jnint5_1.d2e[(jhop - 1) * 7 + j + 
			    i__ * 300 - 301], (ftnlen)sizeof(real));
		}
		e_wsfe();
/* L210: */
	    }
	    if (nextra > 0) {
		io___424.ciunit = jndat1_1.mstjn[5];
		s_wsfe(&io___424);
		i__2 = nextra;
		for (j = 1; j <= i__2; ++j) {
		    do_fio(&c__1, (char *)&jnint5_1.d2e[nhop * 7 + j + i__ * 
			    300 - 301], (ftnlen)sizeof(real));
		}
		e_wsfe();
	    }
	    io___425.ciunit = jndat1_1.mstjn[5];
	    s_wsle(&io___425);
	    e_wsle();
/* L200: */
	}
    } else if (*is == 5) {
/* ...Write out the diagonal elements of the Hessian Matrix and its Trace */
	trace = 0.f;
	nwgts = jnint2_1.mm0[jnint2_1.nl] + jnint2_1.mv0[jnint2_1.nl];
	i__1 = nwgts;
	for (iw = 1; iw <= i__1; ++iw) {
	    trace += jnint5_1.d2e[iw + iw * 300 - 301];
/* L220: */
	}
	io___428.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___428);
	do_fio(&c__1, (char *)&nwgts, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nwgts, (ftnlen)sizeof(integer));
	e_wsfe();
	io___429.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___429);
	e_wsfe();
	io___430.ciunit = jndat1_1.mstjn[5];
	s_wsle(&io___430);
	e_wsle();
	nhop = nwgts / 7;
	nextra = nwgts - nhop * 7;
	i__1 = nhop;
	for (jhop = 1; jhop <= i__1; ++jhop) {
	    io___431.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___431);
	    for (j = 1; j <= 7; ++j) {
		do_fio(&c__1, (char *)&jnint5_1.d2e[(jhop - 1) * 7 + j + ((
			jhop - 1) * 7 + j) * 300 - 301], (ftnlen)sizeof(real))
			;
	    }
	    e_wsfe();
/* L230: */
	}
	if (nextra > 0) {
	    io___432.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___432);
	    i__1 = nextra;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&jnint5_1.d2e[nhop * 7 + j + (nhop * 7 
			+ j) * 300 - 301], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	}
	io___433.ciunit = jndat1_1.mstjn[5];
	s_wsle(&io___433);
	e_wsle();
	io___434.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___434);
	do_fio(&c__1, (char *)&trace, (ftnlen)sizeof(real));
	e_wsfe();
	io___435.ciunit = jndat1_1.mstjn[5];
	s_wsle(&io___435);
	e_wsle();
    }
    return 0;
/* **** END OF JNSTAT **************************************************** */
} /* jnstat_ */

/* *********************************************************************** */
/* Subroutine */ int jntest_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int jnerr_(integer *), jnfeed_(void);
    extern integer jnindx_(integer *, integer *, integer *);

/* ...JetNet subroutine TEST */
/* ...Sets the values of OUT according to given pattern in OIN and */
/* ...current weights */
    if (jndat1_1.mstjn[7] == 0) {
	jnerr_(&c__23);
    }
    jnfeed_();
    i__1 = jnint2_1.m[jnint2_1.nl];
    for (i__ = 1; i__ <= i__1; ++i__) {
	jndat1_1.out[i__ - 1] = jnint1_1.o[jnindx_(&jnint2_1.nl, &i__, &c__0) 
		- 1];
/* L100: */
    }
    return 0;
/* **** END OF JNTEST **************************************************** */
} /* jntest_ */

/* *********************************************************************** */
/* Subroutine */ int jntral_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3, r__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double r_sign(real *, real *), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static real beta, olde;
    static integer idum;
    static real wmax;
    static integer i__;
    static real scale, width;
    extern /* Subroutine */ int jnerr_(integer *);
    extern doublereal errjn_(integer *);
    static integer il;
    extern /* Subroutine */ int jnfeed_(void);
    static integer it, iw;
    static real factor;
    extern /* Subroutine */ int jndelt_(void), jnchop_(integer *), jncogr_(
	    void);
    extern doublereal gausjn_(integer *);
    extern /* Subroutine */ int jnscgr_(void), jnsatm_(void);
    static real switch__, eta, err;
    extern doublereal rjn_(integer *);

/* ...JetNet subroutine TRaining ALgorithm */
/* ...Trains the net. */
    if (jndat1_1.mstjn[7] == 0) {
	jnerr_(&c__22);
    }
    if (jndat1_1.mstjn[8] <= 0) {
	jnerr_(&c__24);
    }
    ++jndat1_1.mstjn[6];
    jnfeed_();
    if (jnint4_1.ilinon == 0) {
	jndelt_();
    }
    err = errjn_(&c__0);
    err /= (real) jnint2_1.m[jnint2_1.nl];
    jndat1_1.parjn[6] = err;
    jnint2_1.er1 += err;
    jnint2_1.er2 += err;
    if (jndat1_1.mstjn[21] != 0) {
	jnsatm_();
    }
    if (jndat1_1.mstjn[6] % jndat1_1.mstjn[1] != 0) {
	return 0;
    }
/* ...update only every MSTJN(2) calls */
    jndat1_1.parjn[7] = jnint2_1.er1 / (real) jndat1_1.mstjn[1];
    jnint2_1.er1 = 0.f;
    if (jndat1_1.mstjn[20] > 0) {
/* ...Include pruning factors */
	beta = 1.f;
	for (il = jnint2_1.nl; il >= 1; --il) {
/* ...set beta in layer IL: */
	    if (jndat2_1.tinv[il - 1] == 0.f) {
		beta *= jndat1_1.parjn[2];
	    } else {
		beta *= (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	    }
/* Computing 2nd power */
	    r__1 = jndat1_1.parjn[17];
	    factor = (real) jndat1_1.mstjn[1] * 2.f * jndat1_1.parjn[13] * (
		    r__1 * r__1) / beta;
	    i__1 = jnint2_1.mm0[il];
	    for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		r__2 = jndat1_1.parjn[17];
/* Computing 2nd power */
		r__3 = jnint1_1.w[i__ - 1];
/* Computing 2nd power */
		r__1 = r__2 * r__2 + r__3 * r__3;
		jnint1_1.dw[i__ - 1] -= factor * jnint1_1.w[i__ - 1] / (r__1 *
			 r__1);
/* L120: */
	    }
	    i__1 = jnint2_1.mv0[il];
	    for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		r__2 = jndat1_1.parjn[17];
/* Computing 2nd power */
		r__3 = jnint1_1.t[i__ - 1];
/* Computing 2nd power */
		r__1 = r__2 * r__2 + r__3 * r__3;
		jnint1_1.dt[i__ - 1] -= factor * jnint1_1.t[i__ - 1] / (r__1 *
			 r__1);
/* L130: */
	    }
/* L110: */
	}
    }
    if (jndat1_1.mstjn[4] == 0) {
/* ...Normal updating: */
	beta = 1.f;
	for (il = jnint2_1.nl; il >= 1; --il) {
/* ...set beta in layer IL: */
	    if (jndat2_1.tinv[il - 1] == 0.f) {
		beta *= jndat1_1.parjn[2];
	    } else {
		beta *= (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	    }
/* ...set eta in layer IL: */
	    if (jndat2_1.etal[il - 1] == 0.f) {
		eta = jndat1_1.parjn[0] / (real) jndat1_1.mstjn[1] * beta;
	    } else {
		eta = jndat2_1.etal[il - 1] / (real) jndat1_1.mstjn[1] * beta;
	    }
	    i__1 = jnint2_1.mm0[il];
	    for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.w[i__ - 1] = (1.f - jndat1_1.parjn[4] * (real) 
			jnint1_1.nself[i__ - 1]) * jnint1_1.w[i__ - 1] + 
			jnint1_1.dw[i__ - 1] * (real) jnint1_1.nself[i__ - 1] 
			* eta;
		jnint1_1.dw[i__ - 1] *= jndat1_1.parjn[1];
/* L210: */
	    }
	    i__1 = jnint2_1.mv0[il];
	    for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__1; ++i__) {
		jnint1_1.t[i__ - 1] = (1.f - jndat1_1.parjn[4] * (real) 
			jnint1_1.ntself[i__ - 1]) * jnint1_1.t[i__ - 1] + 
			jnint1_1.dt[i__ - 1] * (real) jnint1_1.ntself[i__ - 1]
			 * eta;
		jnint1_1.dt[i__ - 1] *= jndat1_1.parjn[1];
/* L220: */
	    }
/* L200: */
	}
	jnint4_1.ilinon = 0;
	jnint4_1.nc = 0;
	jnint4_1.nsc = 0;
    } else if (jndat1_1.mstjn[4] == 1) {
/* ...Manhattan updating: */
/* ...set eta in layer IL: */
	i__1 = jnint2_1.nl;
	for (il = 1; il <= i__1; ++il) {
	    if (jndat2_1.etal[il - 1] == 0.f) {
		eta = jndat1_1.parjn[0] / (real) jndat1_1.mstjn[1];
	    } else {
		eta = jndat2_1.etal[il - 1] / (real) jndat1_1.mstjn[1];
	    }
	    i__2 = jnint2_1.mm0[il];
	    for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__2; ++i__) {
		jnint1_1.w[i__ - 1] = (1.f - jndat1_1.parjn[4] * (real) 
			jnint1_1.nself[i__ - 1]) * jnint1_1.w[i__ - 1] + 
			r_sign(&eta, &jnint1_1.dw[i__ - 1]) * (real) 
			jnint1_1.nself[i__ - 1];
		jnint1_1.dw[i__ - 1] *= jndat1_1.parjn[1];
/* L310: */
	    }
	    i__2 = jnint2_1.mv0[il];
	    for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__2; ++i__) {
		jnint1_1.t[i__ - 1] = (1.f - jndat1_1.parjn[4] * (real) 
			jnint1_1.ntself[i__ - 1]) * jnint1_1.t[i__ - 1] + 
			r_sign(&eta, &jnint1_1.dt[i__ - 1]) * (real) 
			jnint1_1.ntself[i__ - 1];
		jnint1_1.dt[i__ - 1] *= jndat1_1.parjn[1];
/* L320: */
	    }
/* L300: */
	}
	jnint4_1.ilinon = 0;
	jnint4_1.nc = 0;
	jnint4_1.nsc = 0;
    } else if (jndat1_1.mstjn[4] == 2) {
/* ...Langevin updating: */
	beta = 1.f;
	for (il = jnint2_1.nl; il >= 1; --il) {
/* ...set effective beta in layer IL: */
	    if (jndat2_1.tinv[il - 1] == 0.f) {
		beta *= jndat1_1.parjn[2];
	    } else {
		beta *= (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	    }
/* ...set eta in layer IL: */
	    if (jndat2_1.etal[il - 1] == 0.f) {
		eta = jndat1_1.parjn[0] / (real) jndat1_1.mstjn[1] * beta;
	    } else {
		eta = jndat2_1.etal[il - 1] / (real) jndat1_1.mstjn[1] * beta;
	    }
	    i__1 = jnint2_1.mm0[il];
	    for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__1; ++i__) {
		idum = i__;
		jnint1_1.w[i__ - 1] = (1.f - jndat1_1.parjn[4] * (real) 
			jnint1_1.nself[i__ - 1]) * jnint1_1.w[i__ - 1] + 
			jnint1_1.dw[i__ - 1] * (real) jnint1_1.nself[i__ - 1] 
			* eta + gausjn_(&idum) * jndat1_1.parjn[5];
		jnint1_1.dw[i__ - 1] *= jndat1_1.parjn[1];
/* L410: */
	    }
	    i__1 = jnint2_1.mv0[il];
	    for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__1; ++i__) {
		idum = i__;
		jnint1_1.t[i__ - 1] = (1.f - jndat1_1.parjn[4] * (real) 
			jnint1_1.ntself[i__ - 1]) * jnint1_1.t[i__ - 1] + 
			jnint1_1.dt[i__ - 1] * (real) jnint1_1.ntself[i__ - 1]
			 * eta + gausjn_(&idum) * jndat1_1.parjn[5];
		jnint1_1.dt[i__ - 1] *= jndat1_1.parjn[1];
/* L420: */
	    }
/* L400: */
	}
	jnint4_1.ilinon = 0;
	jnint4_1.nc = 0;
	jnint4_1.nsc = 0;
    } else if (jndat1_1.mstjn[4] == 3) {
/* ...Fahlman's Quickprop: */
	wmax = 0.f;
	beta = 1.f;
	for (il = jnint2_1.nl; il >= 1; --il) {
/* ...set beta in layer IL: */
	    if (jndat2_1.tinv[il - 1] == 0.f) {
		beta *= jndat1_1.parjn[2];
	    } else {
		beta *= (r__1 = jndat2_1.tinv[il - 1], dabs(r__1));
	    }
/* ...set eta in layer IL: */
	    if (jndat2_1.etal[il - 1] == 0.f) {
		eta = jndat1_1.parjn[0] / (real) jndat1_1.mstjn[1] * beta;
	    } else {
		eta = jndat2_1.etal[il - 1] / (real) jndat1_1.mstjn[1] * beta;
	    }
	    i__1 = jnint2_1.mm0[il];
	    for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__1; ++i__) {
/* Computing MAX */
/* Computing MIN */
		r__3 = jndat1_1.parjn[20], r__4 = jnint1_1.dw[i__ - 1] / (
			jnint1_1.odw[i__ - 1] - jnint1_1.dw[i__ - 1] + 1e-20f)
			;
		r__1 = -jndat1_1.parjn[20], r__2 = dmin(r__3,r__4);
		scale = dmax(r__1,r__2);
		r__1 = jnint1_1.odw[i__ - 1] * jnint1_1.dw[i__ - 1];
		switch__ = (real) jnint1_1.nself[i__ - 1] * (r_sign(&c_b852, &
			r__1) + .5f);
		jnint1_1.g[i__ - 1] = jnint1_1.dw[i__ - 1] * switch__ * eta + 
			scale * jnint1_1.g[i__ - 1];
		jnint1_1.w[i__ - 1] = (1.f - jndat1_1.parjn[4]) * jnint1_1.w[
			i__ - 1] + jnint1_1.g[i__ - 1];
		jnint1_1.odw[i__ - 1] = jnint1_1.dw[i__ - 1];
		jnint1_1.dw[i__ - 1] = 0.f;
		if ((r__1 = jnint1_1.w[i__ - 1], dabs(r__1)) > wmax) {
		    wmax = (r__2 = jnint1_1.w[i__ - 1], dabs(r__2));
		}
/* L510: */
	    }
	    i__1 = jnint2_1.mv0[il];
	    for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__1; ++i__) {
/* Computing MAX */
/* Computing MIN */
		r__3 = jndat1_1.parjn[20], r__4 = jnint1_1.dt[i__ - 1] / (
			jnint1_1.odt[i__ - 1] - jnint1_1.dt[i__ - 1] + 1e-20f)
			;
		r__1 = -jndat1_1.parjn[20], r__2 = dmin(r__3,r__4);
		scale = dmax(r__1,r__2);
		r__1 = jnint1_1.odt[i__ - 1] * jnint1_1.dt[i__ - 1];
		switch__ = (real) jnint1_1.ntself[i__ - 1] * (r_sign(&c_b852, 
			&r__1) + .5f);
		jnint1_1.g[jnint2_1.mm0[jnint2_1.nl] + i__ - 1] = jnint1_1.dt[
			i__ - 1] * switch__ * eta + scale * jnint1_1.g[
			jnint2_1.mm0[jnint2_1.nl] + i__ - 1];
		jnint1_1.t[i__ - 1] = (1.f - jndat1_1.parjn[4]) * jnint1_1.t[
			i__ - 1] + jnint1_1.g[jnint2_1.mm0[jnint2_1.nl] + i__ 
			- 1];
		jnint1_1.odt[i__ - 1] = jnint1_1.dt[i__ - 1];
		jnint1_1.dt[i__ - 1] = 0.f;
		if ((r__1 = jnint1_1.t[i__ - 1], dabs(r__1)) > wmax) {
		    wmax = (r__2 = jnint1_1.t[i__ - 1], dabs(r__2));
		}
/* L520: */
	    }
/* L500: */
	}
	if (wmax > jndat1_1.parjn[21]) {
/* ...Quickprop is stuck -> reset weights and restart */
	    i__1 = jnint2_1.nl;
	    for (il = 1; il <= i__1; ++il) {
		if (jndat2_1.widl[il - 1] <= 0.f) {
		    width = jndat1_1.parjn[3];
		} else {
		    width = jndat2_1.widl[il - 1];
		}
		i__2 = jnint2_1.mm0[il];
		for (i__ = jnint2_1.mm0[il - 1] + 1; i__ <= i__2; ++i__) {
		    idum = i__;
		    if (width >= 0.f) {
			jnint1_1.w[i__ - 1] = (rjn_(&idum) * 2.f - 1.f) * 
				width;
		    } else {
			jnint1_1.w[i__ - 1] = -rjn_(&idum) * width;
		    }
/* L540: */
		}
		i__2 = jnint2_1.mv0[il];
		for (i__ = jnint2_1.mv0[il - 1] + 1; i__ <= i__2; ++i__) {
		    idum = i__;
		    if (width >= 0.f) {
			jnint1_1.t[i__ - 1] = (rjn_(&idum) * 2.f - 1.f) * 
				width;
		    } else {
			jnint1_1.t[i__ - 1] = -rjn_(&idum) * width;
		    }
/* L550: */
		}
/* L530: */
	    }
	    ++jndat1_1.mstjn[37];
	    if (jndat1_1.mstjn[37] > jndat1_1.mstjn[35]) {
		jnerr_(&c__21);
	    }
	}
	jnint4_1.ilinon = 0;
	jnint4_1.nc = 0;
	jnint4_1.nsc = 0;
    } else if (jndat1_1.mstjn[4] >= 4 && jndat1_1.mstjn[4] <= 8) {
/* ...Conjugate Gradient updating: */
	jncogr_();
    } else if (jndat1_1.mstjn[4] == 9) {
/* ...Minimization terminated - don't update: */
	return 0;
    } else if (jndat1_1.mstjn[4] >= 10 && jndat1_1.mstjn[4] <= 14) {
/* ...Scaled Conjugate Gradient: */
	jnscgr_();
    } else if (jndat1_1.mstjn[4] == 15) {
/* ...Riedmiller's & Braun's Rprop: */
	i__1 = jnint2_1.mm0[jnint2_1.nl];
	for (iw = 1; iw <= i__1; ++iw) {
	    if (jnint1_1.dw[iw - 1] * jnint1_1.odw[iw - 1] > 0.f) {
/* Computing MIN */
/* Computing MAX */
		r__3 = jndat1_1.parjn[32], r__4 = jnint1_1.etav[iw - 1] * 
			jndat1_1.parjn[29];
		r__1 = jndat1_1.parjn[31], r__2 = dmax(r__3,r__4);
		jnint1_1.etav[iw - 1] = dmin(r__1,r__2);
	    } else if (jnint1_1.dw[iw - 1] * jnint1_1.odw[iw - 1] < 0.f) {
/* Computing MIN */
/* Computing MAX */
		r__3 = jndat1_1.parjn[32], r__4 = jnint1_1.etav[iw - 1] * 
			jndat1_1.parjn[30];
		r__1 = jndat1_1.parjn[31], r__2 = dmax(r__3,r__4);
		jnint1_1.etav[iw - 1] = dmin(r__1,r__2);
	    }
	    r__1 = jnint1_1.dw[iw - 1] * (real) jnint1_1.nself[iw - 1];
	    jnint1_1.w[iw - 1] += r_sign(&jnint1_1.etav[iw - 1], &r__1);
	    jnint1_1.odw[iw - 1] = jnint1_1.dw[iw - 1];
	    jnint1_1.dw[iw - 1] = 0.f;
/* L700: */
	}
	i__1 = jnint2_1.mv0[jnint2_1.nl];
	for (it = 1; it <= i__1; ++it) {
	    if (jnint1_1.dt[it - 1] * jnint1_1.odt[it - 1] > 0.f) {
/* Computing MIN */
/* Computing MAX */
		r__3 = jndat1_1.parjn[32], r__4 = jnint1_1.etav[jnint2_1.mm0[
			jnint2_1.nl] + it - 1] * jndat1_1.parjn[29];
		r__1 = jndat1_1.parjn[31], r__2 = dmax(r__3,r__4);
		jnint1_1.etav[jnint2_1.mm0[jnint2_1.nl] + it - 1] = dmin(r__1,
			r__2);
	    } else if (jnint1_1.dt[it - 1] * jnint1_1.odt[it - 1] < 0.f) {
/* Computing MIN */
/* Computing MAX */
		r__3 = jndat1_1.parjn[32], r__4 = jnint1_1.etav[jnint2_1.mm0[
			jnint2_1.nl] + it - 1] * jndat1_1.parjn[30];
		r__1 = jndat1_1.parjn[31], r__2 = dmax(r__3,r__4);
		jnint1_1.etav[jnint2_1.mm0[jnint2_1.nl] + it - 1] = dmin(r__1,
			r__2);
	    }
	    r__1 = jnint1_1.dt[it - 1] * (real) jnint1_1.ntself[it - 1];
	    jnint1_1.t[it - 1] += r_sign(&jnint1_1.etav[jnint2_1.mm0[
		    jnint2_1.nl] + it - 1], &r__1);
	    jnint1_1.odt[it - 1] = jnint1_1.dt[it - 1];
	    jnint1_1.dt[it - 1] = 0.f;
/* L710: */
	}
    } else {
	jnerr_(&c__9);
    }
/* ...do fixed precision weights */
    if (jnint2_1.icpon == 1) {
	jnchop_(&c__0);
    }
/* ...Scale temperature */
    if (jndat1_1.mstjn[21] >= 0) {
	d__1 = (doublereal) jndat1_1.parjn[12];
	d__2 = (doublereal) (1.f / (real) jndat1_1.mstjn[8]);
	scale = pow_dd(&d__1, &d__2);
	jndat1_1.parjn[2] /= scale;
	i__1 = jnint2_1.nl;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jndat2_1.tinv[i__ - 1] /= scale;
/* L600: */
	}
    }
    if (jndat1_1.mstjn[6] % (jndat1_1.mstjn[1] * jndat1_1.mstjn[8]) != 0) {
	return 0;
    }
/* ...Update some parameters every epoch */
    olde = jndat1_1.parjn[8];
    jndat1_1.parjn[8] = jnint2_1.er2 / (real) (jndat1_1.mstjn[1] * 
	    jndat1_1.mstjn[8]);
    jnint2_1.er2 = 0.f;
    if (jndat1_1.mstjn[20] > 0) {
/* ...Update pruning parameters */
	jndat1_1.parjn[9] = jndat1_1.parjn[15] * jndat1_1.parjn[9] + (1.f - 
		jndat1_1.parjn[15]) * jndat1_1.parjn[8];
	if (jndat1_1.parjn[8] < olde || jndat1_1.parjn[8] < jndat1_1.parjn[18]
		) {
	    jndat1_1.parjn[13] += jndat1_1.parjn[14];
	} else if (jndat1_1.parjn[8] < jndat1_1.parjn[9]) {
	    jndat1_1.parjn[13] -= jndat1_1.parjn[14];
	} else {
	    jndat1_1.parjn[13] = jndat1_1.parjn[16] * jndat1_1.parjn[13];
	}
    }
    if (jndat1_1.mstjn[21] != 0) {
/* ...Calculate saturation measures */
	i__1 = jnint2_1.nl;
	for (il = 1; il <= i__1; ++il) {
	    jndat2_1.satm[il - 1] = jnint2_1.sm[il - 1] / (real) (
		    jndat1_1.mstjn[1] * jndat1_1.mstjn[8]);
	    jnint2_1.sm[il - 1] = 0.f;
/* L610: */
	}
    }
    if (jndat1_1.mstjn[21] < 0) {
	i__1 = jnint2_1.nl;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (jndat2_1.tinv[i__ - 1] >= 0.f) {
		goto L620;
	    }
	    if (jndat2_1.satm[i__ - 1] > .75f) {
/* Computing 2nd power */
		r__1 = jndat2_1.satm[i__ - 1] - .5f;
		jndat2_1.tinv[i__ - 1] /= (jndat1_1.parjn[12] - 1.f) * 16.f * 
			(r__1 * r__1) + 1.f;
		goto L630;
	    } else if (jndat2_1.satm[i__ - 1] < .25f) {
/* Computing 2nd power */
		r__1 = .5f - jndat2_1.satm[i__ - 1];
		jndat2_1.tinv[i__ - 1] *= (jndat1_1.parjn[12] - 1.f) * 16.f * 
			(r__1 * r__1) + 1.f;
		goto L630;
	    }
L620:
	    ;
	}
L630:
	;
    }
/* ...Scale parameters: */
    if (jndat1_1.mstjn[4] <= 2) {
	if (jndat1_1.parjn[10] > 0.f) {
/* ...Change eta using 'bold driver': */
	    if (jndat1_1.parjn[8] >= olde) {
		jndat1_1.parjn[0] *= jndat1_1.parjn[10];
		for (i__ = 1; i__ <= 10; ++i__) {
		    jndat2_1.etal[i__ - 1] *= jndat1_1.parjn[10];
/* L640: */
		}
	    } else {
		jndat1_1.parjn[0] *= (1.f - jndat1_1.parjn[10]) * .1f + 1.f;
		for (i__ = 1; i__ <= 10; ++i__) {
		    jndat2_1.etal[i__ - 1] *= (1.f - jndat1_1.parjn[10]) * 
			    .1f + 1.f;
/* L650: */
		}
	    }
	} else if (jndat1_1.parjn[10] < 0.f) {
/* ...Decrease eta geometrically: */
	    jndat1_1.parjn[0] *= dabs(jndat1_1.parjn[10]);
	    for (i__ = 1; i__ <= 10; ++i__) {
		jndat2_1.etal[i__ - 1] *= dabs(jndat1_1.parjn[10]);
/* L660: */
	    }
	}
    }
/* ...Scale alpha: */
    jndat1_1.parjn[1] *= jndat1_1.parjn[11];
/* ...Scale Langevin noise: */
    jndat1_1.parjn[5] *= jndat1_1.parjn[19];
    return 0;
/* **** END OF JNTRAL **************************************************** */
} /* jntral_ */

/* *********************************************************************** */
/* Subroutine */ int jntred_(integer *n, integer *igrad)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);

    /* Local variables */
#define a ((real *)&jnint5_1)
    static real d__[300], e[300], f, g, h__;
    static integer i__, j, k, l;
    static real scale, hh;

/* ...JetNet subroutine Tridiagonal REduction. */
/* ...Householder reduction of the NxN Hessian stored in D2E. */
/* ...This routine is taken from "Numerical Recipies" by W.H.Press */
/* ...et. al., where it is called TRED2. It has been slightly changed to */
/* ...fit into JETNET. */
    if (*n > 1) {
	for (i__ = *n; i__ >= 2; --i__) {
	    l = i__ - 1;
	    h__ = 0.f;
	    scale = 0.f;
	    if (l > 1) {
		i__1 = l;
		for (k = 1; k <= i__1; ++k) {
		    scale += (r__1 = a[i__ + k * 300 - 301], dabs(r__1));
/* L11: */
		}
		if (scale == 0.f) {
		    e[i__ - 1] = a[i__ + l * 300 - 301];
		} else {
		    i__1 = l;
		    for (k = 1; k <= i__1; ++k) {
			a[i__ + k * 300 - 301] /= scale;
/* Computing 2nd power */
			r__1 = a[i__ + k * 300 - 301];
			h__ += r__1 * r__1;
/* L12: */
		    }
		    f = a[i__ + l * 300 - 301];
		    r__1 = sqrt(h__);
		    g = -r_sign(&r__1, &f);
		    e[i__ - 1] = scale * g;
		    h__ -= f * g;
		    a[i__ + l * 300 - 301] = f - g;
		    f = 0.f;
		    i__1 = l;
		    for (j = 1; j <= i__1; ++j) {
			if (*igrad != 0) {
			    a[j + i__ * 300 - 301] = a[i__ + j * 300 - 301] / 
				    h__;
			}
			g = 0.f;
			i__2 = j;
			for (k = 1; k <= i__2; ++k) {
			    g += a[j + k * 300 - 301] * a[i__ + k * 300 - 301]
				    ;
/* L13: */
			}
			if (l > j) {
			    i__2 = l;
			    for (k = j + 1; k <= i__2; ++k) {
				g += a[k + j * 300 - 301] * a[i__ + k * 300 - 
					301];
/* L14: */
			    }
			}
			e[j - 1] = g / h__;
			f += e[j - 1] * a[i__ + j * 300 - 301];
/* L15: */
		    }
		    hh = f / (h__ + h__);
		    i__1 = l;
		    for (j = 1; j <= i__1; ++j) {
			f = a[i__ + j * 300 - 301];
			g = e[j - 1] - hh * f;
			e[j - 1] = g;
			i__2 = j;
			for (k = 1; k <= i__2; ++k) {
			    a[j + k * 300 - 301] = a[j + k * 300 - 301] - f * 
				    e[k - 1] - g * a[i__ + k * 300 - 301];
/* L16: */
			}
/* L17: */
		    }
		}
	    } else {
		e[i__ - 1] = a[i__ + l * 300 - 301];
	    }
	    d__[i__ - 1] = h__;
/* L18: */
	}
    }
    if (*igrad != 0) {
	d__[0] = 0.f;
    }
    e[0] = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*igrad != 0) {
	    l = i__ - 1;
	    if (d__[i__ - 1] != 0.f) {
		i__2 = l;
		for (j = 1; j <= i__2; ++j) {
		    g = 0.f;
		    i__3 = l;
		    for (k = 1; k <= i__3; ++k) {
			g += a[i__ + k * 300 - 301] * a[k + j * 300 - 301];
/* L19: */
		    }
		    i__3 = l;
		    for (k = 1; k <= i__3; ++k) {
			a[k + j * 300 - 301] -= g * a[k + i__ * 300 - 301];
/* L20: */
		    }
/* L21: */
		}
	    }
	}
	d__[i__ - 1] = a[i__ + i__ * 300 - 301];
	if (*igrad != 0) {
	    a[i__ + i__ * 300 - 301] = 1.f;
	    if (l >= 1) {
		i__2 = l;
		for (j = 1; j <= i__2; ++j) {
		    a[i__ + j * 300 - 301] = 0.f;
		    a[j + i__ * 300 - 301] = 0.f;
/* L22: */
		}
	    }
	}
/* L23: */
    }
    return 0;
/* **** END OF JNTRED **************************************************** */
} /* jntred_ */

#undef a


/* *********************************************************************** */
/* Subroutine */ int jntqli_(integer *n, integer *igrad)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);

    /* Local variables */
    static integer iter;
    static real b, c__, d__[300], e[300], f, g;
    static integer i__, k, l, m;
    static real p, r__;
#define z__ ((real *)&jnint5_1)
    extern /* Subroutine */ int jnerr_(integer *);
    static real dd, ss;

/* ...JetNet subroutine Tridiagonal QL algorithm with Implicit shifts. */
/* ...QL algorithm with implicit shifts to determine the eigenvalues and */
/* ...eigenvectors of the NxN Hessian. The subroutine JNTRED must be */
/* ...called before JNTQLI is invoked (this is not checked for). */
/* ...Eigenvalues and eigenvectors are computed if IGRAD is non-zero. */
/* ...At return the eigenvectors are stored as columns in D2E and the */
/* ...eigenvalues are placed in the vector OUT. */
/* ...This routine is taken from "Numerical Recipies" by W.H.Press */
/* ...et. al., where it is called TQLI. It has been slightly modified */
/* ...to fit into JETNET. */
    if (*n > 1) {
	i__1 = *n;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    e[i__ - 2] = e[i__ - 1];
/* L11: */
	}
	e[*n - 1] = 0.f;
	i__1 = *n;
	for (l = 1; l <= i__1; ++l) {
	    iter = 0;
L1:
	    i__2 = *n - 1;
	    for (m = l; m <= i__2; ++m) {
		dd = (r__1 = d__[m - 1], dabs(r__1)) + (r__2 = d__[m], dabs(
			r__2));
		if ((r__1 = e[m - 1], dabs(r__1)) + dd == dd) {
		    goto L2;
		}
/* L12: */
	    }
	    m = *n;
L2:
	    if (m != l) {
		if (iter == 100) {
		    jnerr_(&c__30);
		}
		++iter;
		g = (d__[l] - d__[l - 1]) / (e[l - 1] * 2.f);
/* Computing 2nd power */
		r__1 = g;
		r__ = sqrt(r__1 * r__1 + 1.f);
		g = d__[m - 1] - d__[l - 1] + e[l - 1] / (g + r_sign(&r__, &g)
			);
		ss = 1.f;
		c__ = 1.f;
		p = 0.f;
		i__2 = l;
		for (i__ = m - 1; i__ >= i__2; --i__) {
		    f = ss * e[i__ - 1];
		    b = c__ * e[i__ - 1];
		    if (dabs(f) >= dabs(g)) {
			c__ = g / f;
/* Computing 2nd power */
			r__1 = c__;
			r__ = sqrt(r__1 * r__1 + 1.f);
			e[i__] = f * r__;
			ss = 1.f / r__;
			c__ *= ss;
		    } else {
			ss = f / g;
/* Computing 2nd power */
			r__1 = ss;
			r__ = sqrt(r__1 * r__1 + 1.f);
			e[i__] = g * r__;
			c__ = 1.f / r__;
			ss *= c__;
		    }
		    g = d__[i__] - p;
		    r__ = (d__[i__ - 1] - g) * ss + c__ * 2.f * b;
		    p = ss * r__;
		    d__[i__] = g + p;
		    g = c__ * r__ - b;
		    if (*igrad != 0) {
			i__3 = *n;
			for (k = 1; k <= i__3; ++k) {
			    f = z__[k + (i__ + 1) * 300 - 301];
			    z__[k + (i__ + 1) * 300 - 301] = ss * z__[k + i__ 
				    * 300 - 301] + c__ * f;
			    z__[k + i__ * 300 - 301] = c__ * z__[k + i__ * 
				    300 - 301] - ss * f;
/* L13: */
			}
		    }
/* L14: */
		}
		d__[l - 1] -= p;
		e[l - 1] = g;
		e[m - 1] = 0.f;
		goto L1;
	    }
/* L15: */
	}
    }
/* ...Put eigenvalues in OIN: */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jndat1_1.out[i__ - 1] = d__[i__ - 1];
/* L100: */
    }
    return 0;
/* **** END OF JNTQLI **************************************************** */
} /* jntqli_ */

#undef z__


/* **********************************************************************C */
/* PART TWO: SELF-ORGANIZING MAP NETWORK                                C */
/* **********************************************************************C */
doublereal gjm_(real *x, integer *n)
{
    /* System generated locals */
    real ret_val, r__1, r__2, r__3;

    /* Builtin functions */
    double tanh(doublereal), exp(doublereal);

    /* Local variables */
    extern /* Subroutine */ int jmerr_(integer *);

/* ...JetMap function G. */
/* ...Gives response function N with argument X. */
    if (*n == 1) {
	ret_val = (tanh(*x) + 1.f) * .5f;
    } else if (*n == 2) {
/* Computing MAX */
/* Computing MIN */
	r__3 = -(*x);
	r__1 = -50.f, r__2 = dmin(r__3,50.f);
	ret_val = exp((dmax(r__1,r__2)));
    } else {
	jndat1_1.mstjm[2] = *n;
	jmerr_(&c__11);
    }
    return ret_val;
/* **** END OF GJM ******************************************************* */
} /* gjm_ */

/* *********************************************************************** */
/* Subroutine */ int jmdump_(integer *nf)
{
    /* Format strings */
    static char fmt_600[] = "(26x,\002 Dump of weights generated by\002)";
    static char fmt_610[] = "(\002Unit \002,i2)";
    static char fmt_620[] = "(\002Unit (\002,i2,\002,\002,i2,\002)\002)";
    static char fmt_630[] = "(10(f8.4,1x))";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsue(cilist *), do_uio(integer *, char *, ftnlen), e_wsue(void),
	     s_wsfe(cilist *), e_wsfe(void), s_wsle(cilist *), e_wsle(void), 
	    do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer inod;
#define indw ((integer *)&jnint1_1 + 460000)
    static integer i__, j, k;
    extern /* Subroutine */ int jmerr_(integer *);
    static integer jf;
    extern /* Subroutine */ int jnhead_(void);
    static integer iw, nfsave;
    extern /* Subroutine */ int jmindx_(integer *, integer *, integer *), 
	    jmstat_(integer *);

    /* Fortran I/O blocks */
    static cilist io___481 = { 0, 0, 0, 0, 0 };
    static cilist io___482 = { 0, 0, 0, 0, 0 };
    static cilist io___484 = { 0, 0, 0, 0, 0 };
    static cilist io___486 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___488 = { 0, 0, 0, 0, 0 };
    static cilist io___491 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___492 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___493 = { 0, 0, 0, fmt_630, 0 };


/* ...JetMap subroutine DUMP weights */
/* ...Dumps weights and other characteristics of the */
/* ...net to file NF for use in other programs. */
    if (jndat1_1.mstjm[7] != 1) {
	jmerr_(&c__8);
    }
    if (*nf < 0) {
/* ...Unformatted dump */
	jf = -(*nf);
	io___481.ciunit = jf;
	s_wsue(&io___481);
	do_uio(&c__1, (char *)&c_n30, (ftnlen)sizeof(integer));
	e_wsue();
	io___482.ciunit = jf;
	s_wsue(&io___482);
	do_uio(&c__20, (char *)&jndat1_1.mstjm[0], (ftnlen)sizeof(integer));
	do_uio(&c__20, (char *)&jndat1_1.parjm[0], (ftnlen)sizeof(real));
	e_wsue();
	i__1 = indw[jmint1_1.nodes[3]] - 1;
	for (iw = 1; iw <= i__1; ++iw) {
	    io___484.ciunit = jf;
	    s_wsue(&io___484);
	    do_uio(&c__1, (char *)&jnint1_1.w[iw - 1], (ftnlen)sizeof(real));
	    e_wsue();
/* L100: */
	}
    } else {
/* ...Formatted dump */
	nfsave = jndat1_1.mstjm[5];
	jndat1_1.mstjm[5] = *nf;
	io___486.ciunit = *nf;
	s_wsfe(&io___486);
	e_wsfe();
	jnhead_();
	jmstat_(&c__1);
	jmstat_(&c__2);
	jndat1_1.mstjm[5] = nfsave;
	i__1 = jmint1_1.nodes[3];
	for (inod = 1; inod <= i__1; ++inod) {
	    io___488.ciunit = *nf;
	    s_wsle(&io___488);
	    e_wsle();
	    jmindx_(&inod, &i__, &j);
	    if (jmint1_1.ndim == 1) {
		io___491.ciunit = *nf;
		s_wsfe(&io___491);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		e_wsfe();
	    } else {
		io___492.ciunit = *nf;
		s_wsfe(&io___492);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    iw = indw[inod - 1] - 1;
	    io___493.ciunit = *nf;
	    s_wsfe(&io___493);
	    i__2 = jmint1_1.nodes[0];
	    for (k = 1; k <= i__2; ++k) {
		do_fio(&c__1, (char *)&jnint1_1.w[iw + k - 1], (ftnlen)sizeof(
			real));
	    }
	    e_wsfe();
/* L200: */
	}
    }
    return 0;
/* **** END OF JMDUMP **************************************************** */
} /* jmdump_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jmerr_(integer *ierr)
{
    /* Format strings */
    static char fmt_600[] = "(\002 *** JETMAP ERROR:\002,i2,\002 ***\002)";
    static char fmt_610[] = "(\002 Illegal number of dimensions (\002,i2,"
	    "\002)\002)";
    static char fmt_620[] = "(\002 Total number of input nodes (\002,i6,\002"
	    ") exceeds limit (\002,i6,\002)\002)";
    static char fmt_630[] = "(\002 Total number of network nodes (\002,i6"
	    ",\002) exceeds limit (\002,i6,\002)\002)";
    static char fmt_640[] = "(\002 The number of weights (\002,i6,\002) exce"
	    "eds limit (\002,i5,\002)\002)";
    static char fmt_650[] = "(\002 Network must be initialized (with JMINIT "
	    "or JMREAD) \002,\002before \002,a6,\002 can be called\002)";
    static char fmt_660[] = "(\002 Call to JMINIT after calling JNINIT\002)";
    static char fmt_670[] = "(\002 Undefined response function (\002,i2,\002"
	    ") in GJM.\002)";
    static char fmt_680[] = "(\002 JMREAD cannot read data-file produced by "
	    "JNDUMP\002)";
    static char fmt_690[] = "(\002 Too many warnings issued by JETMAP\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
#define indw ((integer *)&jnint1_1 + 460000)

    /* Fortran I/O blocks */
    static cilist io___496 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___497 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___498 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___499 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___500 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___501 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___502 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___503 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___504 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___505 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___506 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___507 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___508 = { 0, 0, 0, fmt_680, 0 };
    static cilist io___509 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___510 = { 0, 0, 0, fmt_650, 0 };


/* ...JetMap subroutine ERRor. */
/* ...Writes out an error message and stops the execution. */
    if (jndat1_1.mstjn[7] == 1) {
	jndat1_1.mstjm[5] = jndat1_1.mstjn[5];
    }
    io___496.ciunit = jndat1_1.mstjm[5];
    s_wsfe(&io___496);
    do_fio(&c__1, (char *)&(*ierr), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*ierr == 1) {
	io___497.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___497);
	do_fio(&c__1, (char *)&jndat1_1.mstjm[0], (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 2) {
	io___498.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___498);
	do_fio(&c__1, (char *)&jndat1_1.mstjm[9], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&c__1000, (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 3) {
	io___499.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___499);
	i__2 = (i__1 = jndat1_1.mstjm[10] * jndat1_1.mstjm[11], abs(i__1));
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&c__1000, (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 4) {
	io___500.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___500);
	i__1 = jmint1_1.nodes[3] * (jmint1_1.nodes[0] + 1);
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&c_b168, (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 5) {
	io___501.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___501);
	do_fio(&c__1, "JMTEST", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 6) {
	io___502.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___502);
	do_fio(&c__1, "JMTRAL", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 7) {
	io___503.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___503);
	do_fio(&c__1, "JMINDX", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 8) {
	io___504.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___504);
	do_fio(&c__1, "JMDUMP", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 9) {
	io___505.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___505);
	do_fio(&c__1, "JMINWE", (ftnlen)6);
	e_wsfe();
    } else if (*ierr == 10) {
	io___506.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___506);
	e_wsfe();
    } else if (*ierr == 11) {
	io___507.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___507);
	do_fio(&c__1, (char *)&jndat1_1.mstjm[2], (ftnlen)sizeof(integer));
	e_wsfe();
    } else if (*ierr == 12) {
	io___508.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___508);
	e_wsfe();
    } else if (*ierr == 13) {
	io___509.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___509);
	e_wsfe();
    } else if (*ierr == 14) {
	io___510.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___510);
	do_fio(&c__1, "JMNBHD", (ftnlen)6);
	e_wsfe();
    }
    if (*ierr > 0) {
	s_stop("0", (ftnlen)1);
    }
    return 0;
/* **** END OF JMERR ***************************************************** */
} /* jmerr_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jmfeed_(void)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
#define indw ((integer *)&jnint1_1 + 460000)
    static real odum, rsum;
    static integer j, k, iw;
    extern /* Subroutine */ int jmwarn_(integer *);
    extern doublereal gjm_(real *, integer *);

/* ...JetMap subroutine FEED signal to net. */
/* ...Feeds the input signal into the network and calculates */
/* ...MXNDJM and DW. */
    jndat1_1.mxndjm = 0;
    if (jmint1_1.isw[0] == 1) {
/* *** Sigmoidal unit *** */
	odum = 0.f;
	iw = 0;
	i__1 = jmint1_1.nodes[3];
	for (j = 1; j <= i__1; ++j) {
	    rsum = 0.f;
	    i__2 = jmint1_1.nodes[0];
	    for (k = 1; k <= i__2; ++k) {
		++iw;
		jnint1_1.dw[iw - 1] = jndat1_1.oin[k - 1] - jnint1_1.w[iw - 1]
			;
		rsum += jndat1_1.oin[k - 1] * jnint1_1.w[iw - 1];
/* L110: */
	    }
	    r__1 = (rsum - .5f) * jndat1_1.parjm[2];
	    jnint1_1.o[j - 1] = gjm_(&r__1, jmint1_1.isw);
	    if (jnint1_1.o[j - 1] > odum) {
		odum = jnint1_1.o[j - 1];
		jndat1_1.mxndjm = j;
	    }
/* L100: */
	}
    } else if (jmint1_1.isw[0] == 2) {
/* *** Gaussian unit *** */
	odum = exp(-49.f);
	iw = 0;
	i__1 = jmint1_1.nodes[3];
	for (j = 1; j <= i__1; ++j) {
	    rsum = 0.f;
	    i__2 = jmint1_1.nodes[0];
	    for (k = 1; k <= i__2; ++k) {
		++iw;
		jnint1_1.dw[iw - 1] = jndat1_1.oin[k - 1] - jnint1_1.w[iw - 1]
			;
/* Computing 2nd power */
		r__1 = jnint1_1.dw[iw - 1];
		rsum += r__1 * r__1;
/* L210: */
	    }
	    r__1 = rsum * jndat1_1.parjm[2];
	    jnint1_1.o[j - 1] = gjm_(&r__1, jmint1_1.isw);
	    if (jnint1_1.o[j - 1] > odum) {
		odum = jnint1_1.o[j - 1];
		jndat1_1.mxndjm = j;
	    }
/* L200: */
	}
    }
    if (jndat1_1.mxndjm == 0) {
	jmwarn_(&c__1);
    }
    return 0;
/* **** END OF JMFEED **************************************************** */
} /* jmfeed_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jmindx_(integer *inod, integer *i__, integer *j)
{
    /* Local variables */
#define indw ((integer *)&jnint1_1 + 460000)
    extern /* Subroutine */ int jmerr_(integer *);

/* ...JetMap subroutine INDeX. */
/* ...If INOD>0 it returns the (I,J)-coordinates for INOD, if */
/* ...INOD=0 it returns the INOD-number (in O array) for unit (I,J) */
/* ...in net. */
    if (jndat1_1.mstjm[7] != 1) {
	jmerr_(&c__7);
    }
    if (*inod > 0) {
/* ...INOD -> (I,J) */
	if (jmint1_1.ndim == 1) {
	    *j = 1;
	    *i__ = *inod;
	} else if (jmint1_1.ndim == 2) {
	    *j = *inod % jmint1_1.nodes[2];
	    *i__ = *inod / jmint1_1.nodes[2] + 1;
	    if (*j == 0) {
		*j = jmint1_1.nodes[2];
		--(*i__);
	    }
	}
    } else if (*inod == 0) {
/* ...(I,J) -> INOD */
	if (jmint1_1.ndim == 1) {
	    *inod = *i__;
	} else if (jmint1_1.ndim == 2) {
	    *inod = (*i__ - 1) * jmint1_1.nodes[2] + *j;
	}
    }
    return 0;
/* **** END OF JMINDX **************************************************** */
} /* jmindx_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jminit_(void)
{
    /* Format strings */
    static char fmt_600[] = "(26x,\002Self-organized Clustering\002)";
    static char fmt_610[] = "(25x,\002Learning Vector Quantization\002)";
    static char fmt_660[] = "(24x,\002LVQ with neighborhood function\002)";
    static char fmt_620[] = "(16x,\002Weights set randomly between 0.0 and"
	    " + \002,f6.3)";
    static char fmt_630[] = "(20x,\002Weights set randomly between +/- \002,"
	    "f6.3)";
    static char fmt_640[] = "(26x,\002Weights are not normalized\002)";
    static char fmt_650[] = "(28x,\002Weights are normalized\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), s_wsle(cilist *), e_wsle(void), 
	    do_fio(integer *, char *, ftnlen);

    /* Local variables */
    static integer idum;
#define indw ((integer *)&jnint1_1 + 460000)
    extern /* Subroutine */ int jmerr_(integer *), jnhead_(void);
    static integer iw;
    extern /* Subroutine */ int jmsepa_(void), jmnorm_(void), jmstat_(integer 
	    *);
    extern doublereal rjn_(integer *);

    /* Fortran I/O blocks */
    static cilist io___521 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___522 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___523 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___524 = { 0, 0, 0, 0, 0 };
    static cilist io___525 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___526 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___527 = { 0, 0, 0, 0, 0 };
    static cilist io___528 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___529 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___530 = { 0, 0, 0, 0, 0 };


/* ...JetMap subroutine INITialize net. */
/* ...Initializes the net according to switches and parameters in */
/* .../JMDAT1/ and /JMDAT2/. */
/* ...Check if JNINIT has been called: */
    if (jndat1_1.mstjn[7] == 1) {
	jmerr_(&c__10);
    }
/* ...Set parameters: */
    jmsepa_();
/* ...Set initial values of weights: */
    i__1 = indw[jmint1_1.nodes[3]] - 1;
    for (iw = 1; iw <= i__1; ++iw) {
	idum = iw;
	if (jndat1_1.mstjm[1] == 0) {
	    jnint1_1.w[iw - 1] = rjn_(&idum) * jndat1_1.parjm[3];
	} else if (jndat1_1.mstjm[1] == 1) {
	    jnint1_1.w[iw - 1] = (rjn_(&idum) * 2.f - 1.f) * jndat1_1.parjm[3]
		    ;
	}
/* L100: */
    }
/* ...Normalize weights: */
    if (jndat1_1.mstjm[6] == 1) {
	jmnorm_();
    }
/* ...Write statistics on output file: */
    if (jndat1_1.mstjm[5] < 0) {
	return 0;
    }
    jnhead_();
    jmstat_(&c__1);
    if (jndat1_1.mstjm[4] == 0) {
	io___521.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___521);
	e_wsfe();
    } else if (jndat1_1.mstjm[4] == 1) {
	io___522.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___522);
	e_wsfe();
    } else if (jndat1_1.mstjm[4] == 2) {
	io___523.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___523);
	e_wsfe();
    }
    io___524.ciunit = jndat1_1.mstjm[5];
    s_wsle(&io___524);
    e_wsle();
    if (jndat1_1.mstjm[1] == 0) {
	io___525.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___525);
	do_fio(&c__1, (char *)&jndat1_1.parjm[3], (ftnlen)sizeof(real));
	e_wsfe();
    } else if (jndat1_1.mstjm[1] == 1) {
	io___526.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___526);
	do_fio(&c__1, (char *)&jndat1_1.parjm[3], (ftnlen)sizeof(real));
	e_wsfe();
    }
    io___527.ciunit = jndat1_1.mstjm[5];
    s_wsle(&io___527);
    e_wsle();
    if (jndat1_1.mstjm[6] == 0) {
	io___528.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___528);
	e_wsfe();
    } else if (jndat1_1.mstjm[6] == 1) {
	io___529.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___529);
	e_wsfe();
    }
    io___530.ciunit = jndat1_1.mstjm[5];
    s_wsle(&io___530);
    e_wsle();
    return 0;
/* **** END OF JMINIT **************************************************** */
} /* jminit_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jminwe_(integer *inod)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
#define indw ((integer *)&jnint1_1 + 460000)
    static integer k;
    extern /* Subroutine */ int jmerr_(integer *);
    static integer iw;

/* ...JetMap subroutine INitial WEight. */
/* ...Sets the weight vector for unit INOD equal to the */
/* ...input pattern stored in OIN. */
    if (jndat1_1.mstjm[7] != 1) {
	jmerr_(&c__9);
    }
    i__1 = jmint1_1.nodes[0];
    for (k = 1; k <= i__1; ++k) {
	iw = indw[*inod - 1] + k - 1;
	jnint1_1.w[iw - 1] = jndat1_1.oin[k - 1];
/* L100: */
    }
    return 0;
/* **** END OF JMINWE **************************************************** */
} /* jminwe_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jmnbhd_(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Builtin functions */
    integer i_sign(integer *, integer *);
    double sqrt(doublereal);

    /* Local variables */
#define nbhd ((integer *)&jnint1_1 + 310000)
    static integer iend, jend, inod;
#define indw ((integer *)&jnint1_1 + 460000)
    static integer i__, j, nbnod;
    extern /* Subroutine */ int jmerr_(integer *);
    static integer nbnum;
    static real rdist;
    static integer istrt, jstrt, ic, jc;
    extern /* Subroutine */ int jmindx_(integer *, integer *, integer *), 
	    jmwarn_(integer *);
    static integer ind, jnd;

/* ...JetMap subroutine NeighBourHooD. */
/* ...Specifies the neighbourhood to update according to MSTJM(9). */
    if (jndat1_1.mstjm[7] != 1) {
	jmerr_(&c__14);
    }
    if (jmint1_1.ndim == 1 && abs(jndat1_1.mstjm[8]) > 121) {
	jmwarn_(&c__3);
	jndat1_1.mstjm[8] = i_sign(&c__121, &jndat1_1.mstjm[8]);
    } else if (jmint1_1.ndim == 2 && abs(jndat1_1.mstjm[8]) > 5) {
	jmwarn_(&c__2);
	jndat1_1.mstjm[8] = i_sign(&c__5, &jndat1_1.mstjm[8]);
    }
    if (jndat1_1.mstjm[8] == 0) {
	i__1 = jmint1_1.nodes[3];
	for (inod = 1; inod <= i__1; ++inod) {
	    nbhd[inod + 999] = inod;
	    nbhd[inod - 1] = 1;
/* L200: */
	}
	return 0;
    }
    if (jmint1_1.ndim == 1) {
	i__1 = jmint1_1.nodes[3];
	for (inod = 1; inod <= i__1; ++inod) {
	    if (jndat1_1.mstjm[10] >= 0) {
/* ...Non-periodic boundary: */
/* Computing MAX */
		i__2 = 1, i__3 = inod - abs(jndat1_1.mstjm[8]);
		istrt = max(i__2,i__3);
/* Computing MIN */
		i__2 = jndat1_1.mstjm[10], i__3 = inod + abs(jndat1_1.mstjm[8]
			);
		iend = min(i__2,i__3);
	    } else {
/* ...Periodic boundary: */
		istrt = inod - abs(jndat1_1.mstjm[8]);
		iend = istrt + (abs(jndat1_1.mstjm[8]) << 1);
	    }
	    nbnum = 0;
	    i__2 = iend;
	    for (i__ = istrt; i__ <= i__2; ++i__) {
		ind = (i__ - 1) % jndat1_1.mstjm[10] + 1;
		++nbnum;
		nbhd[inod + nbnum * 1000 - 1] = ind;
/* L110: */
	    }
	    nbhd[inod - 1] = nbnum;
/* L100: */
	}
    } else if (jmint1_1.ndim == 2) {
	if (jndat1_1.mstjm[8] > 0) {
/* ...Square neighbourhood: */
	    i__1 = jmint1_1.nodes[3];
	    for (inod = 1; inod <= i__1; ++inod) {
		jmindx_(&inod, &ic, &jc);
		if (jndat1_1.mstjm[10] >= 0) {
/* ...Non-periodic in dim. 1: */
/* Computing MAX */
		    i__2 = 1, i__3 = ic - jndat1_1.mstjm[8];
		    istrt = max(i__2,i__3);
/* Computing MIN */
		    i__2 = jndat1_1.mstjm[10], i__3 = ic + jndat1_1.mstjm[8];
		    iend = min(i__2,i__3);
		} else {
/* ...Periodic in dim. 1: */
		    istrt = abs(jndat1_1.mstjm[10]) + ic - jndat1_1.mstjm[8];
		    istrt = (istrt - 1) % abs(jndat1_1.mstjm[10]) + 1;
		    iend = istrt + (jndat1_1.mstjm[8] << 1);
		}
		if (jndat1_1.mstjm[11] >= 0) {
/* ...Non-periodic in dim.2: */
/* Computing MAX */
		    i__2 = 1, i__3 = jc - jndat1_1.mstjm[8];
		    jstrt = max(i__2,i__3);
/* Computing MIN */
		    i__2 = jndat1_1.mstjm[11], i__3 = jc + jndat1_1.mstjm[8];
		    jend = min(i__2,i__3);
		} else {
/* ...Periodic in dim. 2: */
		    jstrt = abs(jndat1_1.mstjm[11]) + jc - jndat1_1.mstjm[8];
		    jstrt = (jstrt - 1) % abs(jndat1_1.mstjm[11]) + 1;
		    jend = jstrt + (jndat1_1.mstjm[8] << 1);
		}
		nbnum = 0;
		i__2 = iend;
		for (i__ = istrt; i__ <= i__2; ++i__) {
		    i__3 = jend;
		    for (j = jstrt; j <= i__3; ++j) {
			ind = (i__ - 1) % jndat1_1.mstjm[10] + 1;
			jnd = (j - 1) % jndat1_1.mstjm[11] + 1;
			++nbnum;
			nbnod = 0;
			jmindx_(&nbnod, &ind, &jnd);
			nbhd[inod + nbnum * 1000 - 1] = nbnod;
/* L140: */
		    }
/* L130: */
		}
		nbhd[inod - 1] = nbnum;
/* L120: */
	    }
	} else if (jndat1_1.mstjm[8] < 0) {
/* ...Circular neighbourhood: */
	    i__1 = jmint1_1.nodes[3];
	    for (inod = 1; inod <= i__1; ++inod) {
		jmindx_(&inod, &ic, &jc);
		if (jndat1_1.mstjm[10] >= 0) {
/* ...Non-periodic in dim. 1: */
/* Computing MAX */
		    i__2 = 1, i__3 = ic + jndat1_1.mstjm[8];
		    istrt = max(i__2,i__3);
/* Computing MIN */
		    i__2 = jndat1_1.mstjm[10], i__3 = ic - jndat1_1.mstjm[8];
		    iend = min(i__2,i__3);
		} else {
/* ...Periodic in dim. 1: */
		    istrt = abs(jndat1_1.mstjm[10]) + ic + jndat1_1.mstjm[8];
		    istrt = (istrt - 1) % abs(jndat1_1.mstjm[10]) + 1;
		    iend = istrt - (jndat1_1.mstjm[8] << 1);
		    ic = istrt - jndat1_1.mstjm[8];
		}
		if (jndat1_1.mstjm[11] >= 0) {
/* ...Non-periodic in dim.2: */
/* Computing MAX */
		    i__2 = 1, i__3 = jc + jndat1_1.mstjm[8];
		    jstrt = max(i__2,i__3);
/* Computing MIN */
		    i__2 = jndat1_1.mstjm[11], i__3 = jc - jndat1_1.mstjm[8];
		    jend = min(i__2,i__3);
		} else {
/* ...Periodic in dim. 2: */
		    jstrt = abs(jndat1_1.mstjm[11]) + jc + jndat1_1.mstjm[8];
		    jstrt = (jstrt - 1) % abs(jndat1_1.mstjm[11]) + 1;
		    jend = jstrt - (jndat1_1.mstjm[8] << 1);
		    jc = jstrt - jndat1_1.mstjm[8];
		}
		nbnum = 0;
		i__2 = iend;
		for (i__ = istrt; i__ <= i__2; ++i__) {
		    i__3 = jend;
		    for (j = jstrt; j <= i__3; ++j) {
/* Computing 2nd power */
			i__4 = ic - i__;
/* Computing 2nd power */
			i__5 = jc - j;
			rdist = sqrt((real) (i__4 * i__4 + i__5 * i__5));
			if (rdist <= (real) (-jndat1_1.mstjm[8])) {
			    ind = (i__ - 1) % abs(jndat1_1.mstjm[10]) + 1;
			    jnd = (j - 1) % abs(jndat1_1.mstjm[11]) + 1;
			    ++nbnum;
			    nbnod = 0;
			    jmindx_(&nbnod, &ind, &jnd);
			    nbhd[inod + nbnum * 1000 - 1] = nbnod;
			}
/* L170: */
		    }
/* L160: */
		}
		nbhd[inod - 1] = nbnum;
/* L150: */
	    }
	}
    }
    jmint1_1.nbo = jndat1_1.mstjm[8];
    return 0;
/* **** END OF JMNBHD **************************************************** */
} /* jmnbhd_ */

#undef indw
#undef nbhd


/* *********************************************************************** */
/* Subroutine */ int jmnorm_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer inod;
#define indw ((integer *)&jnint1_1 + 460000)
    static integer iwend;
    static real rnorm;
    static integer iw, iwstrt;

/* ...JetMap subroutine NORMalize weights. */
/* ...Normalizes the weights. */
    i__1 = jmint1_1.nodes[3];
    for (inod = 1; inod <= i__1; ++inod) {
	rnorm = 0.f;
	iwstrt = indw[inod - 1];
	iwend = indw[inod - 1] + jmint1_1.nodes[0] - 1;
	i__2 = iwend;
	for (iw = iwstrt; iw <= i__2; ++iw) {
	    rnorm += jnint1_1.w[iw - 1] * jnint1_1.w[iw - 1];
/* L110: */
	}
	rnorm = sqrt(rnorm);
	i__2 = iwend;
	for (iw = iwstrt; iw <= i__2; ++iw) {
	    jnint1_1.w[iw - 1] /= rnorm;
/* L120: */
	}
/* L100: */
    }
    return 0;
/* **** END OF JMNORM **************************************************** */
} /* jmnorm_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jmread_(integer *nf)
{
    /* Format strings */
    static char fmt_690[] = "(a)";
    static char fmt_710[] = "(tr63,f2.1)";
    static char fmt_610[] = "(a,10i7)";
    static char fmt_620[] = "(a,10f7.4)";
    static char fmt_630[] = "(a5,i2)";
    static char fmt_640[] = "(10f9.4)";
    static char fmt_650[] = "(a6,i2,a1,i2,a1)";
    static char fmt_600[] = "(21x,\002Weights read from file\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_rsue(cilist *), do_uio(integer *, char *, ftnlen), e_rsue(void),
	     s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void)
	    , s_cmp(char *, char *, ftnlen, ftnlen), s_rsle(cilist *), e_rsle(
	    void), s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static char line[100];
    static integer inod;
#define indw ((integer *)&jnint1_1 + 460000)
    static integer ndum, i__, j, k;
    extern /* Subroutine */ int jmerr_(integer *);
    static real fvers;
    static integer ivers, jf;
    extern /* Subroutine */ int jnhead_(void);
    static integer iw;
    extern /* Subroutine */ int jmsepa_(void);
    static integer nfsave;
    extern /* Subroutine */ int jmstat_(integer *);

    /* Fortran I/O blocks */
    static cilist io___559 = { 0, 0, 0, 0, 0 };
    static cilist io___561 = { 0, 0, 0, 0, 0 };
    static cilist io___563 = { 0, 0, 0, 0, 0 };
    static cilist io___565 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___567 = { 0, 0, 0, 0, 0 };
    static cilist io___568 = { 0, 0, 0, 0, 0 };
    static cilist io___569 = { 0, 0, 0, fmt_710, 0 };
    static cilist io___571 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___572 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___574 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___575 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___576 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___577 = { 0, 0, 0, 0, 0 };
    static cilist io___579 = { 0, 0, 0, 0, 0 };
    static cilist io___580 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___581 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___583 = { 0, 0, 0, 0, 0 };
    static cilist io___584 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___586 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___587 = { 0, 0, 0, fmt_600, 0 };


/* ...JetMap subroutine READ weights and parameters. */
/* ...Reads weights and parameters from file NF and initializes the net. */
    ndum = jndat1_1.mstjm[5];
    if (*nf < 0) {
/* ...Unformatted read */
	jf = -(*nf);
	io___559.ciunit = jf;
	s_rsue(&io___559);
	do_uio(&c__1, (char *)&ivers, (ftnlen)sizeof(integer));
	e_rsue();
	if (ivers >= 0) {
	    jmerr_(&c__12);
	}
	io___561.ciunit = jf;
	s_rsue(&io___561);
	do_uio(&c__20, (char *)&jndat1_1.mstjm[0], (ftnlen)sizeof(integer));
	do_uio(&c__20, (char *)&jndat1_1.parjm[0], (ftnlen)sizeof(real));
	e_rsue();
	jmsepa_();
	i__1 = indw[jmint1_1.nodes[3]] - 1;
	for (iw = 1; iw <= i__1; ++iw) {
	    io___563.ciunit = jf;
	    s_rsue(&io___563);
	    do_uio(&c__1, (char *)&jnint1_1.w[iw - 1], (ftnlen)sizeof(real));
	    e_rsue();
/* L100: */
	}
    } else {
/* ...Formatted read */
	nfsave = jndat1_1.mstjm[5];
	io___565.ciunit = *nf;
	s_rsfe(&io___565);
	do_fio(&c__1, line, (ftnlen)100);
	e_rsfe();
	if (s_cmp(line + 26, "Du", (ftnlen)2, (ftnlen)2) == 0) {
	    jmerr_(&c__12);
	}
	io___567.ciunit = *nf;
	s_rsle(&io___567);
	e_rsle();
	io___568.ciunit = *nf;
	s_rsle(&io___568);
	e_rsle();
	io___569.ciunit = *nf;
	s_rsfe(&io___569);
	do_fio(&c__1, (char *)&fvers, (ftnlen)sizeof(real));
	e_rsfe();
	ivers = (integer) (fvers * 10.f + .001f);
L900:
	io___571.ciunit = *nf;
	s_rsfe(&io___571);
	do_fio(&c__1, line, (ftnlen)100);
	e_rsfe();
	if (s_cmp(line, "         I      1", (ftnlen)17, (ftnlen)17) != 0) {
	    goto L900;
	}
	io___572.ciunit = *nf;
	s_rsfe(&io___572);
	do_fio(&c__1, line, (ftnlen)10);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjm[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	e_rsfe();
	io___574.ciunit = *nf;
	s_rsfe(&io___574);
	do_fio(&c__1, line, (ftnlen)10);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjm[i__ + 9], (ftnlen)sizeof(
		    integer));
	}
	e_rsfe();
	io___575.ciunit = *nf;
	s_rsfe(&io___575);
	do_fio(&c__1, line, (ftnlen)10);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjm[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	io___576.ciunit = *nf;
	s_rsfe(&io___576);
	do_fio(&c__1, line, (ftnlen)10);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjm[i__ + 9], (ftnlen)sizeof(
		    real));
	}
	e_rsfe();
	io___577.ciunit = *nf;
	s_rsle(&io___577);
	e_rsle();
	jndat1_1.mstjm[5] = nfsave;
	jmsepa_();
	if (jmint1_1.ndim == 1) {
	    i__1 = jmint1_1.nodes[3];
	    for (inod = 1; inod <= i__1; ++inod) {
		io___579.ciunit = *nf;
		s_rsle(&io___579);
		e_rsle();
		io___580.ciunit = *nf;
		s_rsfe(&io___580);
		do_fio(&c__1, line, (ftnlen)100);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		e_rsfe();
		iw = indw[inod - 1] - 1;
		io___581.ciunit = *nf;
		s_rsfe(&io___581);
		i__2 = jmint1_1.nodes[0];
		for (k = 1; k <= i__2; ++k) {
		    do_fio(&c__1, (char *)&jnint1_1.w[iw + k - 1], (ftnlen)
			    sizeof(real));
		}
		e_rsfe();
/* L110: */
	    }
	} else if (jmint1_1.ndim == 2) {
	    i__1 = jmint1_1.nodes[3];
	    for (inod = 1; inod <= i__1; ++inod) {
		io___583.ciunit = *nf;
		s_rsle(&io___583);
		e_rsle();
		io___584.ciunit = *nf;
		s_rsfe(&io___584);
		do_fio(&c__1, line, (ftnlen)100);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, line, (ftnlen)100);
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		do_fio(&c__1, line, (ftnlen)100);
		e_rsfe();
		iw = indw[inod - 1] - 1;
		io___586.ciunit = *nf;
		s_rsfe(&io___586);
		i__2 = jmint1_1.nodes[0];
		for (k = 1; k <= i__2; ++k) {
		    do_fio(&c__1, (char *)&jnint1_1.w[iw + k - 1], (ftnlen)
			    sizeof(real));
		}
		e_rsfe();
/* L120: */
	    }
	}
    }
    jndat1_1.mstjm[5] = ndum;
/* ...Write statistics on output file: */
    if (jndat1_1.mstjm[5] < 0) {
	return 0;
    }
    jnhead_();
    jmstat_(&c__1);
    io___587.ciunit = jndat1_1.mstjm[5];
    s_wsfe(&io___587);
    e_wsfe();
/* L700: */
    return 0;
/* **** END OF JMREAD **************************************************** */
} /* jmread_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jmsepa_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer idim, inod;
#define indw ((integer *)&jnint1_1 + 460000)
    extern /* Subroutine */ int jmerr_(integer *), jmnbhd_(void);

/* ...JetMap subroutine SEt PArameters. */
/* ...Sets parameters in /JMINT1/. */
/* ...Number of dimensions: */
    jmint1_1.ndim = jndat1_1.mstjm[0];
    if (jmint1_1.ndim > 2 || jmint1_1.ndim <= 0) {
	jmerr_(&c__1);
    }
/* ...Number of nodes: */
    if (jndat1_1.mstjm[9] > 1000) {
	jmerr_(&c__2);
    }
    if (jmint1_1.ndim == 1) {
	jmint1_1.nodes[3] = abs(jndat1_1.mstjm[10]);
    } else if (jmint1_1.ndim == 2) {
	jmint1_1.nodes[3] = (i__1 = jndat1_1.mstjm[10] * jndat1_1.mstjm[11], 
		abs(i__1));
    }
    if (jmint1_1.nodes[3] > 1000) {
	jmerr_(&c__3);
    }
    for (idim = 1; idim <= 3; ++idim) {
	jmint1_1.nodes[idim - 1] = (i__1 = jndat1_1.mstjm[idim + 8], abs(i__1)
		);
/* L100: */
    }
/* *** Set internal switches: *** */
/* ...Response function: */
    jmint1_1.isw[0] = jndat1_1.mstjm[2];
/* ...Error measure: */
    jmint1_1.isw[1] = jndat1_1.mstjm[3];
/* *** Calculate weight pointers: *** */
    if (jmint1_1.nodes[3] * jmint1_1.nodes[0] > 150000) {
	jmerr_(&c__4);
    }
    i__1 = jmint1_1.nodes[3] + 1;
    for (inod = 1; inod <= i__1; ++inod) {
	indw[inod - 1] = (inod - 1) * jmint1_1.nodes[0] + 1;
/* L110: */
    }
    jndat1_1.mstjm[7] = 1;
/* ...Initialize neighbourhood: */
    jmnbhd_();
    return 0;
/* **** END OF JMSEPA **************************************************** */
} /* jmsepa_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jmstat_(integer *is)
{
    /* Format strings */
    static char fmt_600[] = "(22x,\002Initialized for a \002,i1,\002-dimensi"
	    "onal map with\002)";
    static char fmt_610[] = "(27x,i3,\002 nodes in dimension \002,i1)";
    static char fmt_620[] = "(27x,i3,\002 input nodes\002)";
    static char fmt_630[] = "(18x,\002Values of parameters and switches in J"
	    "ETNET\002)";
    static char fmt_640[] = "(a10,10i7)";
    static char fmt_650[] = "(a10,10f7.4)";
    static char fmt_660[] = "(5x,\002Time factor for this map:\002,i10)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), e_wsle(void), s_wsfe(cilist *), do_fio(integer *
	    , char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer idim;
#define indw ((integer *)&jnint1_1 + 460000)
    static integer i__, nwfac;

    /* Fortran I/O blocks */
    static cilist io___592 = { 0, 0, 0, 0, 0 };
    static cilist io___593 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___595 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___596 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___597 = { 0, 0, 0, 0, 0 };
    static cilist io___598 = { 0, 0, 0, 0, 0 };
    static cilist io___599 = { 0, 0, 0, 0, 0 };
    static cilist io___600 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___601 = { 0, 0, 0, 0, 0 };
    static cilist io___602 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___604 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___605 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___606 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___607 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___608 = { 0, 0, 0, 0, 0 };
    static cilist io___610 = { 0, 0, 0, fmt_660, 0 };


/* ...JetMap subroutine output STATistics. */
/* ...Statistics chosen by IS is written on file MSTJM(6). */
    if (*is == 1) {
	io___592.ciunit = jndat1_1.mstjm[5];
	s_wsle(&io___592);
	e_wsle();
	io___593.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___593);
	do_fio(&c__1, (char *)&jmint1_1.ndim, (ftnlen)sizeof(integer));
	e_wsfe();
	i__1 = jmint1_1.ndim;
	for (idim = 1; idim <= i__1; ++idim) {
	    io___595.ciunit = jndat1_1.mstjm[5];
	    s_wsfe(&io___595);
	    do_fio(&c__1, (char *)&jmint1_1.nodes[idim], (ftnlen)sizeof(
		    integer));
	    do_fio(&c__1, (char *)&idim, (ftnlen)sizeof(integer));
	    e_wsfe();
/* L100: */
	}
	io___596.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___596);
	do_fio(&c__1, (char *)&jmint1_1.nodes[0], (ftnlen)sizeof(integer));
	e_wsfe();
	io___597.ciunit = jndat1_1.mstjm[5];
	s_wsle(&io___597);
	e_wsle();
	io___598.ciunit = jndat1_1.mstjm[5];
	s_wsle(&io___598);
	e_wsle();
    } else if (*is == 2) {
	io___599.ciunit = jndat1_1.mstjm[5];
	s_wsle(&io___599);
	e_wsle();
	io___600.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___600);
	e_wsfe();
	io___601.ciunit = jndat1_1.mstjm[5];
	s_wsle(&io___601);
	e_wsle();
	io___602.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___602);
	do_fio(&c__1, "I", (ftnlen)1);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	}
	e_wsfe();
	io___604.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___604);
	do_fio(&c__1, "MSTJM I", (ftnlen)7);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjm[i__ - 1], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
	io___605.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___605);
	do_fio(&c__1, "MSTJM 10+I", (ftnlen)10);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.mstjm[i__ + 9], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
	io___606.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___606);
	do_fio(&c__1, "PARJM I", (ftnlen)7);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjm[i__ - 1], (ftnlen)sizeof(
		    real));
	}
	e_wsfe();
	io___607.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___607);
	do_fio(&c__1, "PARJM 10+I", (ftnlen)10);
	for (i__ = 1; i__ <= 10; ++i__) {
	    do_fio(&c__1, (char *)&jndat1_1.parjm[i__ + 9], (ftnlen)sizeof(
		    real));
	}
	e_wsfe();
	io___608.ciunit = jndat1_1.mstjm[5];
	s_wsle(&io___608);
	e_wsle();
    } else {
	nwfac = indw[jmint1_1.nodes[3]] - 1;
	io___610.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___610);
	do_fio(&c__1, (char *)&nwfac, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    return 0;
/* **** END OF JMSTAT **************************************************** */
} /* jmstat_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jmtest_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer inod;
#define indw ((integer *)&jnint1_1 + 460000)
    extern /* Subroutine */ int jmerr_(integer *), jmfeed_(void);

/* ...JetMap subroutine TEST */
/* ...Sets the values of OUT according to given pattern in OIN */
/* ...and current values of weights. */
    if (jndat1_1.mstjm[7] != 1) {
	jmerr_(&c__7);
    }
    jmfeed_();
    i__1 = jmint1_1.nodes[3];
    for (inod = 1; inod <= i__1; ++inod) {
	jndat1_1.out[inod - 1] = jnint1_1.o[inod - 1];
/* L100: */
    }
    return 0;
/* **** END OF JMTEST **************************************************** */
} /* jmtest_ */

#undef indw


/* *********************************************************************** */
/* Subroutine */ int jmtral_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
#define nbhd ((integer *)&jnint1_1 + 310000)
    static integer inbh, inod;
#define indw ((integer *)&jnint1_1 + 460000)
    static integer iwend;
    extern /* Subroutine */ int jmerr_(integer *), jmfeed_(void);
    static integer iw;
    extern /* Subroutine */ int jmnbhd_(void), jmnorm_(void);
    static integer iwstrt;

/* ...JetMap subroutine TRaining ALgorithm. */
/* ...Trains the net according to /JMDAT1/. */
    if (jndat1_1.mstjm[7] != 1) {
	jmerr_(&c__6);
    }
/* ...Update neighbourhood: */
    if (jndat1_1.mstjm[8] != jmint1_1.nbo) {
	jmnbhd_();
    }
    jmfeed_();
    if (jndat1_1.mxndjm == 0) {
	return 0;
    }
    if (jndat1_1.mstjm[4] == 0) {
/* ...Normal updating */
	i__1 = nbhd[jndat1_1.mxndjm - 1];
	for (inbh = 1; inbh <= i__1; ++inbh) {
	    iwstrt = indw[nbhd[jndat1_1.mxndjm + inbh * 1000 - 1] - 1];
	    iwend = indw[nbhd[jndat1_1.mxndjm + inbh * 1000 - 1]] - 1;
	    i__2 = iwend;
	    for (iw = iwstrt; iw <= i__2; ++iw) {
		jnint1_1.w[iw - 1] += jndat1_1.parjm[0] * jnint1_1.dw[iw - 1];
/* L110: */
	    }
/* L100: */
	}
    } else if (jndat1_1.mstjm[4] == 1) {
/* ...Learning Vector Quantization */
	iwstrt = indw[jndat1_1.mxndjm - 1];
	iwend = indw[jndat1_1.mxndjm] - 1;
	if (jndat1_1.out[jndat1_1.mxndjm - 1] < .5f) {
/* ...Wrong answer: */
	    i__1 = iwend;
	    for (iw = iwstrt; iw <= i__1; ++iw) {
		jnint1_1.w[iw - 1] -= jndat1_1.parjm[0] * jnint1_1.dw[iw - 1];
/* L200: */
	    }
	} else {
/* ...Correct answer: */
	    i__1 = iwend;
	    for (iw = iwstrt; iw <= i__1; ++iw) {
		jnint1_1.w[iw - 1] += jndat1_1.parjm[0] * jnint1_1.dw[iw - 1];
/* L210: */
	    }
	}
    } else if (jndat1_1.mstjm[4] == 2) {
/* ...LVQ with neighborhood function: */
	iwstrt = indw[jndat1_1.mxndjm - 1];
	iwend = indw[jndat1_1.mxndjm] - 1;
	i__1 = iwend;
	for (iw = iwstrt; iw <= i__1; ++iw) {
	    jnint1_1.w[iw - 1] += jndat1_1.parjm[0] * jnint1_1.dw[iw - 1] * 
		    jndat1_1.out[jndat1_1.mxndjm - 1];
/* L220: */
	}
    }
    i__1 = jmint1_1.nodes[3];
    for (inod = 1; inod <= i__1; ++inod) {
	jndat1_1.out[inod - 1] = jnint1_1.o[inod - 1];
/* L300: */
    }
/* ...If normalize: */
    if (jndat1_1.mstjm[6] == 1) {
	jmnorm_();
    }
    return 0;
/* **** END OF JMTRAL **************************************************** */
} /* jmtral_ */

#undef indw
#undef nbhd


/* *********************************************************************** */
/* Subroutine */ int jmwarn_(integer *iwarn)
{
    /* Format strings */
    static char fmt_600[] = "(\002 *** JETMAP WARNING:\002,i2,\002 ***\002)";
    static char fmt_610[] = "(\002 No response in net for presented input"
	    "\002)";
    static char fmt_620[] = "(\002 Illegal value of neighbourhood size (\002"
	    ",i2,\002)\002)";
    static char fmt_630[] = "(\002 absolute value of MSTJM(9) must be within"
	    " [0,5]\002)";
    static char fmt_640[] = "(\002 MSTJM(9) set to limit value 5\002)";
    static char fmt_650[] = "(\002 absolute value of MSTJM(9) must be within"
	    " [0,121]\002)";
    static char fmt_660[] = "(\002 MSTJM(9) set to limit value 121\002)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    extern /* Subroutine */ int jmerr_(integer *);

    /* Fortran I/O blocks */
    static cilist io___620 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___621 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___622 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___623 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___624 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___625 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___626 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___627 = { 0, 0, 0, fmt_660, 0 };


/* ...JetMap subroutine WARNing. */
/* ...Writes out a warning on file MSTJM(6). */
    ++jndat1_1.mstjn[33];
    jndat1_1.mstjn[32] = *iwarn;
    if (jndat1_1.mstjn[33] < jndat1_1.mstjn[31]) {
	io___620.ciunit = jndat1_1.mstjm[5];
	s_wsfe(&io___620);
	do_fio(&c__1, (char *)&(*iwarn), (ftnlen)sizeof(integer));
	e_wsfe();
	if (*iwarn == 1) {
	    io___621.ciunit = jndat1_1.mstjm[5];
	    s_wsfe(&io___621);
	    e_wsfe();
	} else if (*iwarn == 2) {
	    io___622.ciunit = jndat1_1.mstjm[5];
	    s_wsfe(&io___622);
	    do_fio(&c__1, (char *)&jndat1_1.mstjm[8], (ftnlen)sizeof(integer))
		    ;
	    e_wsfe();
	    io___623.ciunit = jndat1_1.mstjm[5];
	    s_wsfe(&io___623);
	    e_wsfe();
	    io___624.ciunit = jndat1_1.mstjm[5];
	    s_wsfe(&io___624);
	    e_wsfe();
	} else if (*iwarn == 3) {
	    io___625.ciunit = jndat1_1.mstjm[5];
	    s_wsfe(&io___625);
	    do_fio(&c__1, (char *)&jndat1_1.mstjm[8], (ftnlen)sizeof(integer))
		    ;
	    e_wsfe();
	    io___626.ciunit = jndat1_1.mstjm[5];
	    s_wsfe(&io___626);
	    e_wsfe();
	    io___627.ciunit = jndat1_1.mstjm[5];
	    s_wsfe(&io___627);
	    e_wsfe();
	}
    } else if (jndat1_1.mstjn[30] >= 1) {
	jmerr_(&c__13);
    }
    return 0;
/* **** END OF JMWARN **************************************************** */
} /* jmwarn_ */

/* *********************************************************************** */
/* Subroutine */ int jndata_(void)
{
    return 0;
} /* jndata_ */

/* ...JetNet block DATA */
/* ...Initial values for parameters and switches for JETNET. */
/* ...Brief explanation of parameters and switches : */
/* ... */
/* .../JNDAT1/: */
/* ...Feed-forward net: */
/* ...MSTJN(1) (D=3)      number of layers in net */
/* ...MSTJN(2) (D=10)     number of patterns per update in JNTRAL */
/* ...MSTJN(3) (D=1)      overall transfer function used in net */
/* ...        1 -> g(x)=1/(1+exp(-2x)) */
/* ...        2 -> g(x)=tanh(x) */
/* ...        3 -> g(x)=exp(x) (only used internally for Potts-nodes) */
/* ...        4 -> g(x)=x */
/* ...        5 -> g(x)=1/(1+exp(-2x)) (only used internally for */
/* ...             entropy error) */
/* ...MSTJN(4) (D=0)      error measure */
/* ...       -1 -> log-squared error:     E = -log(1-(o-t)**2) */
/* ...        0 -> summed square error:   E = 0.5*(o-t)**2 */
/* ...        1 -> entropy error:         E = -t*log(o) + (1-t)*log(1-o) */
/* ...      >=2 -> Kullback measure, using Potts nodes of dimension */
/* ...             MSTJN(4):              E = t*log(t/o) */
/* ...MSTJN(5) (D=0)      updating procedure */
/* ...        0 -> standard Back-Propagation updating */
/* ...        1 -> Manhattan updating */
/* ...        2 -> Langevin updating */
/* ...        3 -> Quickprop */
/* ...        4 -> Conjugate Gradient - Polak-Ribiere */
/* ...        5 -> Conjugate Gradient - Hestenes-Stiefel */
/* ...        6 -> Conjugate Gradient - Fletcher-Reeves */
/* ...        7 -> Conjugate Gradient - Shanno */
/* ...        8 -> Terminate Conjugate Gradient search */
/* ...        9 -> No updating */
/* ...       10 -> Scaled Conjugate Gradient - Polak-Ribiere */
/* ...       11 -> Scaled Conjugate Gradient - Hestenes-Stiefel */
/* ...       12 -> Scaled Conjugate Gradient - Fletcher-Reeves */
/* ...       13 -> Scaled Conjugate Gradient - Shanno */
/* ...       14 -> Terminate Scaled Conjugate Gradient Search */
/* ...       15 -> Rprop */
/* ...MSTJN(6) (D=6)      file number for output statistics */
/* ...MSTJN(7) (I)        number of calls to JNTRAL */
/* ...MSTJN(8) (I)        initialization done -> 0 = no */
/* ...MSTJN(9) (D=100)    number of updates per epoch */
/* ...MSTJN(10+I)         number of nodes in layer I (I=0 => input layer) */
/* ...MSTJN(10) (D=16) */
/* ...MSTJN(11) (D=8) */
/* ...MSTJN(12) (D=1) */
/* ...MSTJN(13-20) (D=0) */
/* ...MSTJN(21) (D=0)     pruning (>0 -> on) */
/* ...MSTJN(22) (D=0)     saturation measure (<>0 -> on) */
/* ...                    <0 -> update temperature to give measure ~0.5 */
/* ...MSTJN(23,24) (D=0)  geometry of input nodes for receptive fields */
/* ...MSTJN(25,26) (D=0)  geometry of receptive fields */
/* ...MSTJN(27)    (D=1)  number of hidden nodes per receptive field */
/* ...MSTJN(28-30) (D=0)  precision in bits (0 -> machine precision) for */
/* ...                    sigmoid functions (28), thresholds (29) and */
/* ...                    weights (30) */
/* ...MSTJN(31)    (D=1)  Warning procedure */
/* ...        0 -> No action is taken after a warning */
/* ...        1 -> The execution is stopped after the program */
/* ...             has experienced MSTJN(32) warnings */
/* ...             in any case only MSTJN(32) warning messages are printed */
/* ...             out. */
/* ...MSTJN(32) (D=10)    Maximum number of warning messages to be */
/* ...                    printed. As described above. */
/* ...MSTJN(33) (I)       code for latest warning issued by the program. */
/* ...MSTJN(34) (I)       Number of warnings issued by the program so far. */
/* ...MSTJN(35) (D=10)    Max. number of iterations allowed in line search. */
/* ...MSTJN(36) (D=10)    Max. number of allowed restarts in line search. */
/* ...MSTJN(37) (I)       Status of line search */
/* ...        0 -> Minimum found */
/* ...        1 -> Searching for minimum */
/* ...MSTJN(38) (I)       Number of restarts in Quickprop/ConjGr/ScConjGr */
/* ...MSTJN(39) (I)       Number of calls to JNHESS. */
/* ...MSTJN(40)           not used */
/* ... */
/* ... */
/* ...PARJN(1) (D=0.001)  learning parameter eta */
/* ...PARJN(2) (D=0.5)    momentum term alfa */
/* ...PARJN(3) (D=1.0)    overall inverse temperature beta */
/* ...PARJN(4) (D=0.1)    width of initial weights */
/* ...        > 0 -> [-width,+width] */
/* ...        < 0 -> [0,+width] */
/* ...PARJN(5) (D=0.0)    forgetting parameter epsilon */
/* ...PARJN(6) (D=0.0)    noise width in Langevin equation */
/* ...PARJN(7) (R)        last error per node */
/* ...PARJN(8) (R)        mean error in last update */
/* ...PARJN(9) (R)        mean error last epoch (equal to MSTJN(9) updates) */
/* ...PARJN(10)(R)        weighted mean average used in pruning */
/* ...PARJN(11) (D=1.0)   change in eta (scale factor per epoch) */
/* ...        > 0 -> Geometric with "bold driver" dynamics */
/* ...        < 0 -> Geometric decrease of eta */
/* ...PARJN(12) (D=1.0)   change in momentum alpha (scale factor per epoch) */
/* ...PARJN(13) (D=1.0)   change in temperature (scale factor per epoch) */
/* ...PARJN(14) (D=0.0)   pruning parameter lambda */
/* ...PARJN(15) (D=1.E-6) change in lambda */
/* ...PARJN(16) (D=0.9)   parameter gamma used for calculation of PARJN(10) */
/* ...PARJN(17) (D=0.9)   pruning "cut-off" */
/* ...PARJN(18) (D=1.0)   scale parameter W(0), used in pruning */
/* ...PARJN(19) (D=0.0)   target error when pruning */
/* ...PARJN(20) (D=1.0)   decrease in Langevin noise (scale factor per epoch) */
/* ...PARJN(21) (D=1.75)  maximum scale for Quickprop updating */
/* ...PARJN(22) (D=1000.) maximum allowed size of weights in Quickprop */
/* ...PARJN(23) (D=0.0)   constant added to g'(x) to avoid 'flat spot' */
/* ...PARJN(24) (D=0.1)   line search convergence parameter (0 < ... < 1) */
/* ...PARJN(25) (D=0.05)  tolerance of minimum in line search */
/* ...PARJN(26) (D=0.001) minimum allowed change in error in line search */
/* ...PARJN(27) (D=2.0)   maximum allowed step size in line search */
/* ...PARJN(28) (D=1.E-4) constant sigma_0 used in SCG */
/* ...PARJN(29) (D=1.E-6) initial value for lambda in SCG */
/* ...PARJN(30) (D=1.2)   scale-up factor used in Rprop */
/* ...PARJN(31) (D=0.5)   scale-down factor used in Rprop */
/* ...PARJN(32) (D=50.)   maximum scale-up factor in Rprop */
/* ...PARJN(33) (D=1.E-6) minimum scale-down factor in Rprop */
/* ...PARJN(34-40)        not used */
/* ... */
/* ... */
/* ...Self-organizing net: */
/* ...MSTJM(1)   (D=1)    number of dimensions in net */
/* ...MSTJM(2)   (D=0)    symmetry of initial weights */
/* ...        0 -> [0,+width] */
/* ...        1 -> [-width,+width] */
/* ...MSTJM(3)   (D=2)    response function */
/* ...        1 -> g(x)=0.5*(1.0+tanh(x) : for normalized data */
/* ...        2 -> g(x)=exp(-x)          : for unnormalized data */
/* ...MSTJM(4)   (D=1)    error measure */
/* ...        1 -> summed square error */
/* ...MSTJM(5)   (D=0)    updating procedure */
/* ...        0 -> unsupervized clustering & topological ordering */
/* ...        1 -> Learning Vector Quantization (LVQ 1) */
/* ...        2 -> as 1, but with neighborhood function. */
/* ...MSTJM(6)   (D=6)     output file number */
/* ...MSTJM(7)   (D=0)     normalize weights or not */
/* ...        0 -> unnormalized */
/* ...        1 -> normalized */
/* ...MSTJM(8) (I)         initialization done */
/* ...MSTJM(9)   (D=0)     neighbourhood size */
/* ...        0< -> square neighbourhood */
/* ...        <0 -> circular neighbourhood */
/* ...MSTJM(10)  (D=8)     number of input nodes */
/* ...MSTJM(11)  (D=10)    number of nodes in dimension 1. */
/* ...        <0 -> periodic boundary */
/* ...MSTJM(12)  (D=1)     number of nodes in dimension 2. */
/* ...        <0 -> periodic boundary */
/* ...MSTJM(13-20)         not used */
/* ... */
/* ... */
/* ...PARJM(1)   (D=0.001) learning parameter eta */
/* ...PARJM(2)   (D=0.0)   not used */
/* ...PARJM(3)   (D=0.01)  overall inverse temperature beta */
/* ...PARJM(4)   (D=0.5)   initial width of weights */
/* ...PARJM(5-20)          not used */
/* ... */
/* ... */
/* .../JNDAT2/: */
/* ...TINV(I) (D=0.0)     inverse temperature of layer I (if 0 use PARJN(3)) */
/* ... */
/* ...IGFN(I) (D=0)       sigmoid function for layer I (if 0 use MSTJN(3)) */
/* ... */
/* ...ETAL(I) (D=0.0)     learning parameter in layer I (if 0 use PARJN(1)) */
/* ... */
/* ...WIDL(I) (D=0.0)     initial width in layer I (if 0 use PARJN(4)) */
/* ... */
/* ...SATM(I) (R)         saturation measure "S" for layer I. */
/* ...      MSTJN(3)=1 -> S = sum[(1.-2.*O(J))**2] */
/* ...      MSTJN(3)=2 -> S = sum[O(J)**2] */
/* ... */
/* ...End of description */
/* **** END OF JNDATA **************************************************** */

/* *********************************************************************** */
/* Subroutine */ int jntdec_(integer *method)
{
    /* Format strings */
    static char fmt_600[] = "(31x,\002JETNET Test-Deck\002)";
    static char fmt_610[] = "(15x,\002Two overlapping Gaussian distributions"
	    " in \002,i2,\002 dimensions.\002)";
    static char fmt_620[] = "(15x,\002Their standard deviations are \002,f3."
	    "1,\002 and \002,f3.1)";
    static char fmt_621[] = "(15x,\002Their mean values are separated by "
	    "\002,f4.2)";
    static char fmt_625[] = "(15x,\002Generating training and test patterns."
	    "..\002)";
    static char fmt_626[] = "(15x,\002...done generating data.\002)";
    static char fmt_660[] = "(\002 Undefined training algorithm in call to J"
	    "NTDEC\002)";
    static char fmt_630[] = "(\002   Epoch   /  Training  /  General. \002)";
    static char fmt_640[] = "(i8,2x,2(\002 /\002,f9.3,2x))";
    static char fmt_650[] = "(\002 The optimal generalization performance "
	    "is \002,f4.1,\002%\002)";
    static char fmt_670[] = "(\002 Backprop should reach (81.0 +- 2.2)% in 1"
	    "00 epochs\002)";
    static char fmt_680[] = "(\002 Manhattan should reach (84.3 +- 0.6)% in "
	    "100 epochs\002)";
    static char fmt_690[] = "(\002 Langevin should reach (82.9 +- 1.8)% in 1"
	    "00 epochs\002)";
    static char fmt_700[] = "(\002 Quickprop should reach (82.8 +- 8.8)% in "
	    "100 epochs\002)";
    static char fmt_710[] = "(\002 Polak-Ribiere CG should reach (79.0 +- 7."
	    "0)% in 100\002,\002 epochs\002)";
    static char fmt_720[] = "(\002 Hestenes-Stiefel CG should reach (79.8 +-"
	    " 5.6)% in 100\002,\002 epochs\002)";
    static char fmt_730[] = "(\002 Fletcher-Reeves CG should reach (79.6 +- "
	    "5.6)% in 100\002,\002 epochs\002)";
    static char fmt_740[] = "(\002 Shanno CG should reach (71.7 +- 11.6)% in"
	    " 100 epochs\002)";
    static char fmt_750[] = "(\002 Polak-Ribiere SCG should reach (84.0 +- 1"
	    ".6)% in 100\002,\002 epochs\002)";
    static char fmt_760[] = "(\002 Hestenes-Stiefel SCG should reach (84.1 +"
	    "- 2.6)% in 100\002,\002 epochs\002)";
    static char fmt_770[] = "(\002 Fletcher-Reeves SCG should reach (81.4 +-"
	    " 5.2)% in 100\002,\002 epochs\002)";
    static char fmt_780[] = "(\002 Shanno SCG should reach (70.7 +- 8.1)% in"
	    " 100 epochs\002)";
    static char fmt_790[] = "(\002 Rprop should reach (83.5 +- 2.2)% in 100 "
	    "epochs\002)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen),
	     s_wsle(cilist *), e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer ipat, idum;
    static real test, tout[15000];
    static integer i__;
    static real train, trnmx;
    static integer ip, iepoch;
    extern doublereal gausjn_(integer *);
    extern /* Subroutine */ int jntral_(void), jninit_(void);
    static integer nright;
    extern /* Subroutine */ int jntest_(void);
    static real testmx;
    extern doublereal rjn_(integer *);
    static real tin[75000]	/* was [15000][5] */;

    /* Fortran I/O blocks */
    static cilist io___628 = { 0, 0, 0, fmt_600, 0 };
    static cilist io___629 = { 0, 0, 0, fmt_610, 0 };
    static cilist io___630 = { 0, 0, 0, fmt_620, 0 };
    static cilist io___631 = { 0, 0, 0, fmt_621, 0 };
    static cilist io___632 = { 0, 0, 0, 0, 0 };
    static cilist io___633 = { 0, 0, 0, fmt_625, 0 };
    static cilist io___639 = { 0, 0, 0, fmt_626, 0 };
    static cilist io___640 = { 0, 0, 0, fmt_660, 0 };
    static cilist io___641 = { 0, 0, 0, 0, 0 };
    static cilist io___642 = { 0, 0, 0, fmt_630, 0 };
    static cilist io___650 = { 0, 0, 0, fmt_640, 0 };
    static cilist io___651 = { 0, 0, 0, 0, 0 };
    static cilist io___652 = { 0, 0, 0, fmt_650, 0 };
    static cilist io___653 = { 0, 0, 0, fmt_670, 0 };
    static cilist io___654 = { 0, 0, 0, fmt_680, 0 };
    static cilist io___655 = { 0, 0, 0, fmt_690, 0 };
    static cilist io___656 = { 0, 0, 0, fmt_700, 0 };
    static cilist io___657 = { 0, 0, 0, fmt_710, 0 };
    static cilist io___658 = { 0, 0, 0, fmt_720, 0 };
    static cilist io___659 = { 0, 0, 0, fmt_730, 0 };
    static cilist io___660 = { 0, 0, 0, fmt_740, 0 };
    static cilist io___661 = { 0, 0, 0, fmt_750, 0 };
    static cilist io___662 = { 0, 0, 0, fmt_760, 0 };
    static cilist io___663 = { 0, 0, 0, fmt_770, 0 };
    static cilist io___664 = { 0, 0, 0, fmt_780, 0 };
    static cilist io___665 = { 0, 0, 0, fmt_790, 0 };


/* ...JetNet subroutine Test-DECk */
/* ...Runs a test-program using data from two overlapping Gaussian */
/* ...distributions in the input space. The test-program uses the */
/* ...method specified by METHOD. */
    io___628.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___628);
    e_wsfe();
    io___629.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___629);
    do_fio(&c__1, (char *)&c__5, (ftnlen)sizeof(integer));
    e_wsfe();
    io___630.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___630);
    do_fio(&c__1, (char *)&c_b15, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&c_b1164, (ftnlen)sizeof(real));
    e_wsfe();
    io___631.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___631);
    do_fio(&c__1, (char *)&c_b1168, (ftnlen)sizeof(real));
    e_wsfe();
    io___632.ciunit = jndat1_1.mstjn[5];
    s_wsle(&io___632);
    e_wsle();
/* ...Generate data: */
    io___633.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___633);
    e_wsfe();
    for (ipat = 1; ipat <= 15000; ++ipat) {
	idum = ipat;
	if (rjn_(&idum) > .5f) {
	    for (i__ = 1; i__ <= 5; ++i__) {
		tin[ipat + i__ * 15000 - 15001] = gausjn_(&idum) * 1.f;
/* L110: */
	    }
	    tout[ipat - 1] = 1.f;
	} else {
	    tin[ipat - 1] = gausjn_(&idum) * 2.f + 0.f;
	    for (i__ = 2; i__ <= 5; ++i__) {
		tin[ipat + i__ * 15000 - 15001] = gausjn_(&idum) * 2.f;
/* L120: */
	    }
	    tout[ipat - 1] = 0.f;
	}
/* L100: */
    }
    io___639.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___639);
    e_wsfe();
/* ...Set network architecture: MSTJN(1)-layered network with */
/* ...MSTJN(11) hidden nodes, MSTJN(12) output nodes and */
/* ...MSTJN(10) inputs. */
    jndat1_1.mstjn[0] = 3;
    jndat1_1.mstjn[9] = 5;
    jndat1_1.mstjn[10] = 10.f;
    jndat1_1.mstjn[11] = 1;
/* ...Set sigmoid function: */
    jndat1_1.mstjn[2] = 1;
/* ...Initial width of weights: */
    jndat1_1.parjn[3] = .5f;
/* ...Choose updating method */
    jndat1_1.mstjn[4] = *method;
    if (jndat1_1.mstjn[4] == 8 || jndat1_1.mstjn[4] == 9 || jndat1_1.mstjn[4] 
	    == 14 || jndat1_1.mstjn[4] < 0 || jndat1_1.mstjn[4] > 15) {
	io___640.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___640);
	e_wsfe();
	s_stop("0", (ftnlen)1);
    }
/* ...Initialize network: */
    jninit_();
/* ...Set parameters suitable for the given method of updating */
    if (jndat1_1.mstjn[4] == 0) {
/* ...Normal Backprop */
	jndat1_1.parjn[0] = 2.f;
	jndat1_1.parjn[1] = .5f;
	jndat1_1.parjn[10] = .999f;
    } else if (jndat1_1.mstjn[4] == 1) {
/* ...Manhattan */
	jndat1_1.parjn[0] = .05f;
	jndat1_1.parjn[1] = .5f;
	jndat1_1.parjn[10] = -.99f;
    } else if (jndat1_1.mstjn[4] == 2) {
/* ...Langevin */
	jndat1_1.parjn[0] = 1.f;
	jndat1_1.parjn[1] = .5f;
	jndat1_1.parjn[5] = .01f;
	jndat1_1.parjn[10] = .999f;
	jndat1_1.parjn[19] = .99f;
    } else if (jndat1_1.mstjn[4] == 3) {
/* ...Quickprop */
	jndat1_1.parjn[0] = 2.f;
	jndat1_1.parjn[1] = 0.f;
	jndat1_1.parjn[5] = 0.f;
	jndat1_1.parjn[10] = 1.f;
	jndat1_1.parjn[19] = 1.f;
	jndat1_1.mstjn[1] = 5000;
    } else if (jndat1_1.mstjn[4] >= 4 && jndat1_1.mstjn[4] <= 7) {
/* ...Conjugate Gradient */
	jndat1_1.parjn[0] = 1.f;
	jndat1_1.mstjn[1] = 5000;
    } else if (jndat1_1.mstjn[4] >= 10 && jndat1_1.mstjn[4] <= 13) {
/* ...Scaled Conjugate Gradient */
	jndat1_1.mstjn[1] = 5000;
    } else if (jndat1_1.mstjn[4] == 15) {
/* ...Rprop */
	jndat1_1.parjn[0] = 1.f;
	jndat1_1.mstjn[1] = 5000;
    }
/* ...Define the size of one epoch. Note that for batch training, the */
/* ...number of patterns per update, MSTJN(2), must be set to the */
/* ...total number of training patterns, and hence MSTJN(9), the */
/* ...number of updates per epoch must be set to one. */
/* Computing MAX */
    i__1 = 1, i__2 = 5000 / jndat1_1.mstjn[1];
    jndat1_1.mstjn[8] = max(i__1,i__2);
/* ...Other parameters keep their default values. */
    io___641.ciunit = jndat1_1.mstjn[5];
    s_wsle(&io___641);
    e_wsle();
    io___642.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___642);
    e_wsfe();
    testmx = 0.f;
    trnmx = 0.f;
/* ...Main loop over epochs: */
    for (iepoch = 1; iepoch <= 100; ++iepoch) {
/* ...Training loop: */
	nright = 0;
	for (ip = 1; ip <= 5000; ++ip) {
	    if (jndat1_1.mstjn[4] <= 2) {
/* ...Note that for non-batch training it is often a good idea to pick */
/* ...training patterns at random */
		ipat = (integer) (rjn_(&ip) * 5e3f) + 1;
	    } else {
		ipat = ip;
	    }
/* ...Put pattern into OIN: */
	    i__1 = jndat1_1.mstjn[9];
	    for (i__ = 1; i__ <= i__1; ++i__) {
		jndat1_1.oin[i__ - 1] = tin[ipat + i__ * 15000 - 15001];
/* L320: */
	    }
/* ...Put target output value into OUT: */
	    jndat1_1.out[0] = tout[ipat - 1];
/* ...Invoke training algorithm: */
	    jntral_();
/* ...Calculate performance on training set: */
	    if ((r__1 = jndat1_1.out[0] - tout[ipat - 1], dabs(r__1)) < .5f) {
		++nright;
	    }
/* L310: */
	}
	train = (real) nright / 5e3f;
	if (iepoch % 10 == 0) {
/* ...Testing loop: */
	    nright = 0;
	    for (ipat = 5001; ipat <= 15000; ++ipat) {
/* ...Put pattern into OIN: */
		i__1 = jndat1_1.mstjn[9];
		for (i__ = 1; i__ <= i__1; ++i__) {
		    jndat1_1.oin[i__ - 1] = tin[ipat + i__ * 15000 - 15001];
/* L340: */
		}
/* ...Get network output: */
		jntest_();
/* ...Calculate performance on test set (=generalization): */
		if ((r__1 = jndat1_1.out[0] - tout[ipat - 1], dabs(r__1)) < 
			.5f) {
		    ++nright;
		}
/* L330: */
	    }
	    test = (real) nright / 1e4f;
	    if (jndat1_1.mstjn[4] > 3 && jndat1_1.mstjn[4] < 15) {
		if (train > trnmx) {
		    trnmx = train;
		    testmx = test;
		}
		test = testmx;
		train = trnmx;
	    }
/* ...Display performance: */
	    io___650.ciunit = jndat1_1.mstjn[5];
	    s_wsfe(&io___650);
	    do_fio(&c__1, (char *)&iepoch, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&train, (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&test, (ftnlen)sizeof(real));
	    e_wsfe();
	}
/* ...Terminate CG and SCG training: */
	if (iepoch == 99) {
	    if (jndat1_1.mstjn[4] > 3 && jndat1_1.mstjn[4] < 15) {
		if (jndat1_1.mstjn[4] < 9) {
		    jndat1_1.mstjn[4] = 8;
		} else {
		    jndat1_1.mstjn[4] = 14;
		}
		trnmx = 0.f;
		testmx = 0.f;
	    }
	}
/* L300: */
    }
    io___651.ciunit = jndat1_1.mstjn[5];
    s_wsle(&io___651);
    e_wsle();
    io___652.ciunit = jndat1_1.mstjn[5];
    s_wsfe(&io___652);
    do_fio(&c__1, (char *)&c_b1194, (ftnlen)sizeof(real));
    e_wsfe();
    if (*method == 0) {
	io___653.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___653);
	e_wsfe();
    } else if (*method == 1) {
	io___654.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___654);
	e_wsfe();
    } else if (*method == 2) {
	io___655.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___655);
	e_wsfe();
    } else if (*method == 3) {
	io___656.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___656);
	e_wsfe();
    } else if (*method == 4) {
	io___657.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___657);
	e_wsfe();
    } else if (*method == 5) {
	io___658.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___658);
	e_wsfe();
    } else if (*method == 6) {
	io___659.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___659);
	e_wsfe();
    } else if (*method == 7) {
	io___660.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___660);
	e_wsfe();
    } else if (*method == 10) {
	io___661.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___661);
	e_wsfe();
    } else if (*method == 11) {
	io___662.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___662);
	e_wsfe();
    } else if (*method == 12) {
	io___663.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___663);
	e_wsfe();
    } else if (*method == 13) {
	io___664.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___664);
	e_wsfe();
    } else if (*method == 15) {
	io___665.ciunit = jndat1_1.mstjn[5];
	s_wsfe(&io___665);
	e_wsfe();
    }
    return 0;
/* **** END OF JNTDEC **************************************************** */
} /* jntdec_ */

/* *********************************************************************** */
doublereal rjn_(integer *idum)
{
    /* System generated locals */
    real ret_val;

    /* Local variables */
    static real runi;
    static integer i__, j, k, l, m;
    static real s, t, twom24;
    static integer i24, ii, ij, jj, kl;

/* ...JetNet function Random number generator. */
/* ...Generates random numbers uniformly distributed in ]0,1[ */
/* ...This function is taken from the Lund program JETSET */
/* ...written by T. Sjostrand. */
/* ...The algorithm is due to Marsaglia, Zaman and Tsang: */
/* ...Stat. Prob. Lett., vol. 9, (1990) */
/* ...This function is very much based on a routine */
/* ...written by F.James: F.James, Comp. Phys. Comm., vol 60 (1990). */
/* ...Names have been changed w.r.t. the JETSET function RLU and */
/* ...the switch that determines the current position has been */
/* ...removed. */
/* ...Initialize generation from given seed. */
    if (jndatr_1.mrjn[1] == 0) {
	ij = jndatr_1.mrjn[0] / 30082 % 31329;
	kl = jndatr_1.mrjn[0] % 30082;
	i__ = ij / 177 % 177 + 2;
	j = ij % 177 + 2;
	k = kl / 169 % 178 + 1;
	l = kl % 169;
	for (ii = 1; ii <= 97; ++ii) {
	    s = 0.f;
	    t = .5f;
	    for (jj = 1; jj <= 24; ++jj) {
		m = i__ * j % 179 * k % 179;
		i__ = j;
		j = k;
		k = m;
		l = (l * 53 + 1) % 169;
		if (l * m % 64 >= 32) {
		    s += t;
		}
/* L100: */
		t *= .5f;
	    }
/* L110: */
	    jndatr_1.rrjn[ii - 1] = s;
	}
	twom24 = 1.f;
	for (i24 = 1; i24 <= 24; ++i24) {
/* L120: */
	    twom24 *= .5f;
	}
	jndatr_1.rrjn[97] = twom24 * 362436.f;
	jndatr_1.rrjn[98] = twom24 * 7654321.f;
	jndatr_1.rrjn[99] = twom24 * 16777213.f;
	jndatr_1.mrjn[1] = 1;
	jndatr_1.mrjn[2] = 0;
	jndatr_1.mrjn[3] = 97;
	jndatr_1.mrjn[4] = 33;
    }
/* ...Generate next random number. */
L130:
    runi = jndatr_1.rrjn[jndatr_1.mrjn[3] - 1] - jndatr_1.rrjn[jndatr_1.mrjn[
	    4] - 1];
    if (runi < 0.f) {
	runi += 1.f;
    }
    jndatr_1.rrjn[jndatr_1.mrjn[3] - 1] = runi;
    --jndatr_1.mrjn[3];
    if (jndatr_1.mrjn[3] == 0) {
	jndatr_1.mrjn[3] = 97;
    }
    --jndatr_1.mrjn[4];
    if (jndatr_1.mrjn[4] == 0) {
	jndatr_1.mrjn[4] = 97;
    }
    jndatr_1.rrjn[97] -= jndatr_1.rrjn[98];
    if (jndatr_1.rrjn[97] < 0.f) {
	jndatr_1.rrjn[97] += jndatr_1.rrjn[99];
    }
    runi -= jndatr_1.rrjn[97];
    if (runi < 0.f) {
	runi += 1.f;
    }
    if (runi <= 0.f || runi >= 1.f) {
	goto L130;
    }
/* ...Update counters. Random number to output. */
    ++jndatr_1.mrjn[2];
    if (jndatr_1.mrjn[2] == 1000000000) {
	++jndatr_1.mrjn[1];
	jndatr_1.mrjn[2] = 0;
    }
    ret_val = runi;
    return ret_val;
/* **** END OF RJN ******************************************************* */
} /* rjn_ */

