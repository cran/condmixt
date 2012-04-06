#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define PI M_PI
#define MINSIGMA 0.1   /* 0.000001*/
#define MINPI 0.000001
#define MINXI 0.000001

/* Constants for quantile computations */
#define FACTOR 1.6
#define NTRY 50
#define IBIS 500

double lambertw(double z);

void lambertwR(double *z, double *w0);

double erfc(double x);

double erf(double x);

double sign(double x);

void hpnll(double *theta, double *x, int *n, double *nll,double *nllgrad);

double normpdf(double mu, double sigma, double x);

double normcdf(double mu, double sigma, double x);

double normlogpdf(double mu, double sigma, double x);

double gpdpdf(double xi, double beta, double x);

double gpdlogpdf(double xi, double beta, double x);

double hpdf(double xi, double mu, double sigma, double x);

double hcdf(double xi, double mu, double sigma, double x);

double hlogpdf(double xi, double mu, double sigma, double x);

double softplus(double x);

double softplusinv(double x);

void ummhfwd(double *theta, int *m, double *params);

void ummhbwd(double *params, int *m, double *theta);

void ummhprint(double *params, int *m);

double ummhpdf(double *params, int m, double x);

void ummgpdfR(double *params, int *m, double *x, int *n, double *p);

void ummgcdfR(double *params, int *m, double *x, int *n, double *p);

void ummhpdfR(double *params, int *m, double *x, int *n, double *p);

double ummhcdf(double *params, int m, double x);

void ummhcdfR(double *params, int *m, double *x, int *n, double *p);

void hpdfgrad(double xi, double mu, double sigma, double x,
	      double *dhdtheta);

void ummhnll(double *theta, int *m, double *x, int *n, double *nll, 
	     double *nllgrad);


void ummhnll_bimodal_tailpen(double *theta, int *m, double * lambda, double *w,
			     double * beta, double *mu, double *sigma, double *x, int *n, 
			     double *nll,double *nllgrad);


int ummqbrack(double (*cdf)(double *,int, double),
	      double *params, int m, double q, double *a, double *b);

void ummquant(double (*cdf)(double *,int, double),
	       double (*pdf)(double *,int, double),
	       double *params, int m, double q, double *aa, double *bb,
	      double tol, int itmax, double *xq);

double ummgpdf(double *params, int m, double x);

double ummgcdf(double *params, int m, double x);

void ummgbwd(double *params, int *m, double *theta);


#define HMAX 50 /* maximum number of hidden units */

/* structure containing the info for a conditional mixture model:

   - Linear component has nout(d+1) parameters arranged like this,
   b being the biases and w the linear weights:
          1     2         d+1     
   psi=[b(1) w(1,1) ... w(1,d) ...      ==> dedicated to output 1
            ...        ...
     (nout-1)*(d+1)+1          nout(d+1)
      b(nout) w(nout,1) ... w(nout,d)]  ==> dedicated to output nout

   - Non-linear component is an array of pointers of size hmax. Each 
   element of the array contains the weights related to one hidden
   unit. For instance, for hidden unit h we have, where c and v are 
   the hidden layer bias and weights and w the outer layer weights:
                1     2      d+1   d+2       d+nout+1
   theta[h-1]=[c(1) v(1) ... v(d)  w(1) ...  w(nout)] 
 */

typedef struct {
  double *psi; /* weights related to the linear component */
  double *theta[HMAX];  /* neural network weights, empty if h = 0 */
  int h;  /* number of hidden units */
  int m;  /* number of components */
  int s;  /* number of a second type of components (for double hybrid cmm) */
  int d;  /* dimension of the input */
  int nout; /* dimension of neural network output */
} CMM;


void cmmhprint(double *params,int *d,int *h,int *m);

void cmmprint (CMM *net);

void nnlin(CMM *net, double *x, double *a, double *z);

void cmmhfwd(CMM *net, double *x, double *params, double *a, double *z);

void cmmhfwdR(double *params,int *d,int *h,int *m, double *x, 
	      int *n, double *params_mixt, double *a, double *z);

void cmmhnll(CMM *net, double *x, double *y, int n, 
	     double *nll, double *nllgrad);

void cmmhnllR(double *params,int *d,int *h,int *m, double *x, 
	      double *y, int *n, double *nll, double *nllgrad);


void cmmhnll_bimodal_tailpen(CMM *net, double *x, double *y, int n, 
			     double *lambda,double *w, double *beta,double *mu,  
			     double *sigma, double *nll, double *nllgrad);

void cmmhnll_bimodal_tailpenR(double *params,int *d,int *h,int *m, double *x, 
			      double *y, int *n, double *lambda, double *w, 
			      double *beta, double *mu, double *sigma, double *nll, 
			      double *nllgrad);


void cmmhquant(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq);

void cmmhquant_trunc(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
		     double *xq);

void cmmhfwd_dirac(CMM *net, double *x, double *params, double *a, double *z);

void cmmhfwd_diracR(double *params,int *d,int *h,int *m, double *x, 
		    int *n, double *params_mixt, double *a, double *z);

void cmmhnll_bimodal_tailpen_dirac(CMM *net, double *x, double *y, int n, 
			     double *lambda,double *w, double *beta, double *mu, 
				   double *sigma, double *nll, double *nllgrad);


void cmmhnll_bimodal_tailpen_diracR(double *params,int *d,int *h,int *m, double *x, 
			      double *y, int *n, double *lambda, double *w, 
			      double *beta, double *mu, double *sigma, double *nll, 
				    double *nllgrad);


void cmmhnll_dirac(CMM *net, double *x, double *y, int n, double *nll);

void cmmhnll_diracR(double *params,int *d,int *h,int *m, double *x, 
		    double *y, int *n, double *nll);


void cmmhquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
		     double *xq);


void cmmhquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
		     double *xq);

void cmmhcquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
		      double *xq);

void cmmgfwd(CMM *net, double *x, double *params, double *a, double *z);

void cmmgnll(CMM *net, double *x, double *y, int n, 
	     double *nll, double *nllgrad);

void cmmgnllR(double *params,int *d,int *h,int *m, double *x, 
	      double *y, int *n, double *nll, double *nllgrad);

void cmmgquant(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq);

void cmmgquant_trunc(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
		     double *xq);

void cmmgfwdR(double *params,int *d,int *h,int *m, double *x, 
	      int *n, double *params_mixt, double *a, double *z);

void cmmgfwd_dirac(CMM *net, double *x, double *params, double *a, double *z);

void cmmgfwd_diracR(double *params,int *d,int *h,int *m, double *x, 
		    int *n, double *params_mixt, double *a, double *z);

void cmmgnll_dirac(CMM *net, double *x, double *y, int n, 
		   double *nll, double *nllgrad);

void cmmgnll_diracR(double *params,int *d,int *h,int *m, double *x, 
		    double *y, int *n, double *nll,double *nllgrad);

void cmmgquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
		     double *xq);

void cmmgcquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
		      double *xq);

void cmmlnll_dirac(CMM *net, double *x, double *y, int n, 
		   double *nll, double *nllgrad);

void cmmlnll_diracR(double *params,int *d,int *h,int *m, double *x, 
		    double *y, int *n, double *nll,double *nllgrad);

void cmmlquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
		     double *xq);


/* Gamma-Bernouilli conditional distribution: parameters are predicted by a neural network */
void cmmbergam_fwd(CMM *net, double *x, double *params, double *a, double *z);

void cmmbergam_fwdR(double *params,int *d,int *h, double *x, 
		int *n, double *params_bergam, double *a, double *z);


void cmmbergam_nll(CMM *net, double *x, double *y, int n, double *nll, double *nllgrad);

void cmmbergam_nllR(double *params,int *d,int *h, double *x, 
		    double *y, int *n, double *nll,double *nllgrad);
