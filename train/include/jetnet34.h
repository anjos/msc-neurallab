/* Hello emacs, this is -*-c-*- */
/* $Id: jetnet34.h,v 1.1 2001/07/13 15:48:02 andre Exp $ */
/* André Rabello <Andre.Rabello@ufrj.br> */

/* The meat here is produced by f2c with flag -P */

#ifndef JETNET34_H_
#define JETNET34_H_

extern E_f errjn_(integer *idum);
extern E_f gausjn_(integer *idum);
extern E_f gjn_(integer *ind, real *x, integer *n);
extern int jncgbe_(real *betak, integer *iop);
extern int jnchop_(integer *ichp);
extern int jncogr_(void);
extern int jndelt_(void);
extern int jndump_(integer *nf);
extern int jnerr_(integer *ierr);
extern int jnfeed_(void);
extern int jnhead_(void);
extern int jnheig_(integer *igrad);
extern int jnhess_(void);
extern integer jnindx_(integer *il, integer *i__, integer *j);
extern int jninit_(void);
extern int jnlins_(void);
extern int jnread_(integer *nf);
extern int jnrold_(integer *nf);
extern int jnsatm_(void);
extern int jnscgr_(void);
extern int jnsefi_(integer *ila, integer *i1, integer *i2, integer *j1, integer *j2, integer *no);
extern int jnsepa_(void);
extern int jnstat_(integer *is);
extern int jntest_(void);
extern int jntral_(void);
extern int jntred_(integer *n, integer *igrad);
extern int jntqli_(integer *n, integer *igrad);
extern E_f gjm_(real *x, integer *n);
extern int jmdump_(integer *nf);
extern int jmerr_(integer *ierr);
extern int jmfeed_(void);
extern int jmindx_(integer *inod, integer *i__, integer *j);
extern int jminit_(void);
extern int jminwe_(integer *inod);
extern int jmnbhd_(void);
extern int jmnorm_(void);
extern int jmread_(integer *nf);
extern int jmsepa_(void);
extern int jmstat_(integer *is);
extern int jmtest_(void);
extern int jmtral_(void);
extern int jmwarn_(integer *iwarn);
extern int jndata_(void);
extern int jntdec_(integer *method);
extern E_f rjn_(integer *idum);
/* comlen jndat1_ 8484 */
/* comlen jnint1_ 3672000 */
/* comlen jnint2_ 232 */
/* comlen jngaus_ 8 */
/* comlen jnsigm_ 16000 */
/* comlen jndat2_ 200 */
/* comlen jnint4_ 80 */
/* comlen jnint3_ 36 */
/* comlen jnint5_ 360000 */
/* comlen jmint1_ 64 */
/* comlen jndatr_ 420 */

#endif /* JETNET34_H_ */
