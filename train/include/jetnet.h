#ifndef JETNET_H_
#define JETNET_H_

#include "jetnet34.h"

/* Common Block Declarations */

struct jndat1_1_ {
    integer mstjn[40];
    real parjn[40];
    integer mstjm[20];
    real parjm[20], oin[1000], out[1000];
    integer mxndjm;
};

#define jndat1_1 (*(struct jndat1_1_ *) &jndat1_)

struct {
    real o[2000], a[2000], d__[2000], t[2000], dt[2000], w[150000], dw[150000]
	    ;
    integer nself[150000], ntself[2000];
    real g[152000], odw[150000], odt[2000], etav[152000];
} jnint1_;

#define jnint1_1 jnint1_

struct jnint2_1_ {
    integer m[11], mv0[11], mm0[11], ng[10], nl, ipott;
    real er1, er2, sm[10];
    integer icpon;
};

#define jnint2_1 (*(struct jnint2_1_ *) &jnint2_)

struct jngaus_1_ {
    integer iset;
    real gasdev;
};

#define jngaus_1 (*(struct jngaus_1_ *) &jngaus_)

struct {
    real gpjn[2000], gppjn[2000];
} jnsigm_;

#define jnsigm_1 jnsigm_

struct jndat2_1_ {
    real tinv[10];
    integer igfn[10];
    real etal[10], widl[10], satm[10];
};

#define jndat2_1 (*(struct jndat2_1_ *) &jndat2_)

struct jnint4_1_ {
    integer ilinon, nc;
    real g2;
    integer nit;
    real errln[4], derrln, stepln[4], stepmn, errmn;
    integer ieval, isucc, icurve, nsc;
    real gvec2;
};

#define jnint4_1 (*(struct jnint4_1_ *) &jnint4_)

struct {
    integer nxin, nyin, nxrf, nyrf, nxhrf, nyhrf, nhrf, nrfw, nhprf;
} jnint3_;

#define jnint3_1 jnint3_

struct {
    real d2e[90000]	/* was [300][300] */;
} jnint5_;

#define jnint5_1 jnint5_

struct jmint1_1_ {
    integer ndim, isw[10], nodes[4], nbo;
};

#define jmint1_1 (*(struct jmint1_1_ *) &jmint1_)

struct jndatr_1_ {
    integer mrjn[5];
    real rrjn[100];
};

#define jndatr_1 (*(struct jndatr_1_ *) &jndatr_)

#endif /* JETNET_H_ */
