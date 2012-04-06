#include "condmixt.h"

/* Lambert W function. */
double lambertw(double z){
  double tolrel, tolabs, w0, errabs, errrel, e0, wj;
  tolrel = pow(10.0,-6.0);
  tolabs = tolrel;
  w0 = 1.0/2.0;
  errabs = fabs(z - w0 * exp(w0));
  errrel = errabs/z;

  while ((errabs > tolabs) & (errrel > tolrel)){
    e0 = exp(w0);
    wj = w0-(w0*e0-z)/(e0*(w0+1)-(w0+2)*(w0*e0-z)/(2*w0+2));
    w0 = wj;
    /*    mexPrintf("w0 = %f\n",w0);*/
    errabs = fabs(z-w0*exp(w0));
    errrel = errabs/z;
  }
  return(w0);
}

/* Lambert W function for use with R */
void lambertwR(double *z, double *w0){
  double tolrel, tolabs, errabs, errrel, e0, wj;
  tolrel = pow(10.0,-6.0);
  tolabs = tolrel;
  *w0 = 1.0/2.0;
  errabs = fabs(*z - *w0 * exp(*w0));
  errrel = errabs/ *z;

  while ((errabs > tolabs) & (errrel > tolrel)){
    e0 = exp(*w0);
    wj = *w0-(*w0*e0-*z)/(e0*(*w0+1)-(*w0+2)*(*w0*e0-*z)/(2*(*w0)+2));
    *w0 = wj;
    /*    mexPrintf("w0 = %f\n",*w0);*/
    errabs = fabs(*z-*w0*exp(*w0));
    errrel = errabs/ *z;
  }

}

/* Complementary error function */
double erfc(double x){
  double t,z,ans;
  z=fabs(x); 
  t=1.0/(1.0+0.5*z); 
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+ 
       t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
      	t*(-0.82215223+t*0.17087277))))))))); 
	return x>= 0.0 ? ans : 2.0-ans; 		  
}

/* Error function */
double erf(double x){
  double t,z,ans;
  z=fabs(x); 
  t=1.0/(1.0+0.5*z); 
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+ 
       t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
      	t*(-0.82215223+t*0.17087277))))))))); 
  return x>= 0.0 ? 1.0-ans : ans-1.0; 		  
}

/* Sign function */
double sign(double x){
  double y;
  if (x > 0.0)
    y = 1.0;
  else{
    if (x < 0.0)
      y = -1.0;
    else
      y = 0.0;		
  }
  return(y);
}

void hpnll(double *theta, double *x, int *n, double *nll, double *nllgrad){
  int i;
  double w, gamma, beta, alpha, xi, mu, sigma;
  double dgdx, dbdx, dbds, dadx, dads;

  xi=theta[0];
  mu=theta[1];
  sigma=exp(theta[2]);

  w = lambertw((1.0 + xi)*(1.0 + xi)/(2.0*PI));
  gamma = 3.0/2.0 + erf(sign(1.0+xi)*sqrt(w/2.0))/2.0;
  alpha = sign(1.0+xi) * sigma * sqrt(w) + mu;

  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;
  nllgrad[0]=0;
  nllgrad[1]=0;
  nllgrad[2]=0;
  
    
  /* Compute partial derivatives that do not depend on data points:
     derivatives of gamma, beta and alpha with respect to xi, mu, sigma */
  dgdx=(exp(-w/2.0) * sqrt(w/(2.0*PI))) / (fabs(1.0 + xi) * (1.0 + w));
  dbdx=sign(1.0 + xi) * sigma * sqrt(w) / (1.0 + w);
  dbds=fabs(1.0 + xi) / sqrt(w);
  dadx=sigma * sqrt(w) / (fabs(1.0 + xi) * (1.0 + w));
  dads=sign(1.0 + xi) * sqrt(w);
 

  for (i=0; i<*n; i++){

    /* Compute neg-log-like  and its gradient */
    if (x[i] <= alpha){
	*nll+=log(gamma*sigma)+(x[i]-mu)*(x[i]-mu)/(2.0*sigma*sigma)+log(2*PI)/2;
	nllgrad[0]+=dgdx/gamma;
	nllgrad[1]+=-(x[i]-mu)/(sigma*sigma);
	nllgrad[2]+=1.0-(x[i]-mu)*(x[i]-mu)/(sigma*sigma);
    } 
    else {
      beta = sigma*fabs(1.0+xi)/sqrt(w);
      if (xi == 0.0){
	  *nll+= log(gamma*beta)+(x[i]-alpha)/beta;
	  nllgrad[0]+=dgdx/gamma + (dbdx*(1.0-(x[i]-alpha)/beta)-dadx)/beta;
	  nllgrad[1]+=-1.0/beta;
	  nllgrad[2]+=sigma/beta *((1.0-(x[i]-alpha)/beta)*dbds-dads);
      }
      else if (xi > 0.0 || (xi < 0.0 && (x[i] < alpha - beta/xi))){
	  *nll+=log(gamma*beta)+(1.0/xi + 1.0)*log(1.0+xi*(x[i]-alpha)/beta);
	  nllgrad[0]+=dgdx/gamma + dbdx/beta * 
	    (1.0-(xi+1.0)*(x[i]-alpha)/(beta+xi*(x[i]-alpha))) - 
	    log(1.0+xi*(x[i]-alpha)/beta)/(xi*xi) + 
	    (1.0+1.0/xi)*(x[i]-alpha)/(beta+xi*(x[i]-alpha)) -
	    (1.0+xi)/(beta+xi*(x[i]-alpha))*dadx;
	  nllgrad[1]+=-(1.0+xi)/(beta+xi*(x[i]-alpha));
	  nllgrad[2]+=(dbds/beta * 
		       (1.0-(xi+1.0)*(x[i]-alpha)/(beta+xi*(x[i]-alpha))) - 
		       (xi+1.0)/(beta+xi*(x[i]-alpha)) * dads)*sigma;
      }
      else{
	  *nll+=INFINITY;
	  nllgrad[0]+=INFINITY;
	  nllgrad[1]+=INFINITY;
	  nllgrad[2]+=INFINITY;
      }
    }

  } /* for i */
}



/* Univariate Normal density N(mu,sigma^2) evaluated at x. */
double normpdf(double mu, double sigma, double x){
  double z = (x-mu)/sigma;
  return(exp(-z*z/2.0)/(sqrt(2.0*PI)*sigma));
}


/* Univariate Normal distribution fonction N(mu,sigma^2) evaluated at x. */
double normcdf(double mu, double sigma, double x){
  double erfc(double x);
  double z = (x-mu)/sigma;
  return(0.5 * erfc(-z/sqrt(2)));
}


/* Logarithm of the univariate Normal density N(mu,sigma^2) evaluated at x. */
double normlogpdf(double mu, double sigma, double x){
  double z = (x-mu)/sigma;
  return(-z*z/2.0 - log(2.0*PI)/2.0 -log(sigma));
}



/* Generalized Pareto density G(x;xi,beta) evaluated at x. */
double gpdpdf(double xi, double beta, double x){
  double z;
  /* Check domain conditions. */
  if ((x < 0.0) | ((xi < 0.0) & (x > -beta/xi)) )
    return(0.0);
  z = x/beta;
  if (xi == 0.0)
    return(exp(-z)/beta);
  else{
    return(pow((1.0 + xi * z),(-1.0/xi - 1.0))/beta);
  }
}


/* Logarithm of the Generalized Pareto density G(xi,beta) evaluated at x. */
double gpdlogpdf(double xi, double beta, double x){
  double z;
  /* Check domain conditions. */
  if ((x < 0.0) | ((xi < 0.0) & (x > -beta/xi)) )
    return(INFINITY);
  z = x/beta;
  if (xi == 0.0)
    return(-z-log(beta));
  else{
    /*if (xi > 0.0)*/
    return((-1.0/xi - 1.0)*log(1.0 + xi * z)-log(beta));
  }
}

/* Hybrid Pareto density h(x; xi, mu, sigma) evaluated at one point. */
double hpdf(double xi, double mu, double sigma, double x){
  double w, gamma, beta, alpha,y;
  
  w = lambertw((1.0 + xi)*(1.0 + xi)/(2.0*PI));
  alpha = sign(1.0+xi) * sigma * sqrt(w) + mu;
  gamma = 3.0/2.0 + erf(sign(1.0+xi)*sqrt(w/2.0))/2.0;

  if (x <= alpha){
    y=normpdf(mu, sigma, x)/gamma;
  }
  else{
    beta = sigma*fabs(1.0+xi)/sqrt(w);
    y=gpdpdf(xi, beta, x-alpha)/gamma;    
  }
  
  return(y);
}


/* Hybrid Pareto cdf h(x; xi, mu, sigma) evaluated at one point. */
double hcdf(double xi, double mu, double sigma, double x){
  double w, gamma, beta, alpha,y;
  
  w = lambertw((1.0 + xi)*(1.0 + xi)/(2.0*PI));
  gamma = 3.0/2.0 + erf(sign(1.0+xi)*sqrt(w/2.0))/2.0;
  alpha = sign(1.0+xi) * sigma * sqrt(w) + mu;
  
  if (x <= alpha){
    y=normcdf(mu, sigma, x)/gamma;
  }
  else{
    beta = sigma*fabs(1.0+xi)/sqrt(w);
    if (xi == 0.0)
      y=(normcdf(mu,sigma,alpha)+1.0-exp(-(x-alpha)/beta))/gamma;
    else 
      y=(normcdf(mu,sigma,alpha)+1.0-pow(1.0+xi*(x-alpha)/beta,-1/xi))/gamma;
  }
  
  return(y);
}

/* Logarithm of the Hybrid Pareto density h(x; xi, mu, sigma) 
   evaluated at one point. */ 
double hlogpdf(double xi, double mu, double sigma, double x){
  double w, gamma, beta, alpha,y;
  
  w = lambertw((1.0 + xi)*(1.0 + xi)/(2.0*PI));
  gamma = 3.0/2.0 + erf(sign(1.0+xi)*sqrt(w/2.0))/2.0;
  alpha = sign(1.0+xi) * sigma * sqrt(w) + mu;
  
  if (x <= alpha){
    y=normlogpdf(mu, sigma, x)-log(gamma);
  }
  else{
    beta = sigma*fabs(1.0+xi)/sqrt(w);
    y=gpdlogpdf(xi, beta, x-alpha)-log(gamma);    
  }
  return(y);
}



/* --------------------- UMMH : Unconditional Mixture Model -------------------
   -----------------------  With Hybrid Pareto components   -------------------*/

/* Softplus function ensures positivity: y = log(1+exp(x)) */
double softplus(double x){
  double y;
  if(x>0.0){
    y = log(1.0+exp(-x));
    if (isfinite(y) & !isnan(y))
      y+=x;
    else
      y=x;
  }
  else {
    y= log(1.0+exp(x));
    if (!finite(y) | isnan(y))
      y=0.0;
  }
  return(y);
}

/* Softplus inverse function: y = log(exp(x)-1) */
double softplusinv(double x){
  double y;
  if(x>0.0){
    y = exp(x);
    if (isfinite(y) & !isnan(y))
      y=log(y-1.0);
    else
      y=x;
  }
  else {
    y= log(exp(x)-1.0);
  }
  return(y);
}


/* "Forward" pass in UMMH model */
void ummhfwd(double *theta, int *m, double *params){
  int i;
  double fprior=1.0;

  /* Apply transfer functions to mixture parameters. */
  for(i=*m-1;i>-1;i--){
    if(i>0){
      /* priors in [0,1] */
      params[i]=(1/(1+exp(-theta[4*i-1]))*(1.0-2*MINPI)+MINPI)*fprior;       
      fprior=fprior-params[i];
    }
    else
      params[i]=fprior;
    /* Rprintf("pi[%d] = %e \n",i+1,params[i]);*/
    params[*m+i]=softplus(theta[4*i])+MINXI;   /* tail parameters */
    params[2**m+i]=theta[4*i+1]; /* centers */
    params[3**m+i]=softplus(theta[4*i+2])+MINSIGMA; /* spread parameters */
  }

}

/* "Backward" pass in UMMH model */
void ummhbwd(double *params, int *m, double *theta){
  int i;
  double tmp,fprior=1.0;

  /* Apply inverse transfer functions to mixture parameters. */
  for(i=*m-1;i>-1;i--){
    if(i>0){
      tmp=params[i]/fprior;
      theta[4*i-1] = log(tmp/(1-tmp));
      fprior=fprior-params[i];
    }
    theta[4*i]=softplusinv(params[*m+i]-MINXI);  /* tail parameters */
    theta[4*i+1]=params[2**m+i]; /* centers */
    theta[4*i+2]=softplusinv(params[3**m+i]-MINSIGMA); /* spread parameters */
  }

}


/* Print the hybrid Pareto mixture model parameters */
void ummhprint(double *params, int *m){
  int j;
  /* double *params;
  params = (double*) R_alloc(4**m, sizeof(double));
  
  ummhfwd(theta,m,params);*/

  for (j=0; j<*m; j++){
    Rprintf("prior[%d] = %f, xi[%d] = %f, mu[%d] = %f, sigma[%d] = %f\n",
	      j+1,params[j],j+1,params[*m+j],j+1,params[2**m+j],j+1,params[3**m+j]);
  }

}


/* Hybrid Pareto mixture probability density function */
double ummhpdf(double *params, int m, double x){
  int j;
  double act,prob=0.0;;

  for(j=0; j<m; j++){     
    act=hpdf(params[m+j],params[2*m+j],params[3*m+j],x);
    prob+=params[j]*act;
  }		

  return(prob);
}


/* Hybrid Pareto mixture pdf at n points for use in R*/
void ummhpdfR(double *params, int *m, double *x, int *n, double *p){
  int i;

  for(i=0;i<*n;i++)
    p[i]=ummhpdf(params,*m,x[i]);
}


/* Hybrid Pareto mixture cumulative distribution function at n points*/
double ummhcdf(double *params, int m, double x){
  int j;
  double act,prob=0.0;

  for(j=0; j<m; j++){     
    act=hcdf(params[m+j],params[2*m+j],params[3*m+j],x);
    prob+=params[j]*act;
  }		
  return(prob);
}


/* Hybrid Pareto mixture cdf at n points for use in R*/
void ummhcdfR(double *params, int *m, double *x, int *n, double *p){
  int i;

  for(i=0;i<*n;i++)
    p[i]=ummhcdf(params,*m,x[i]);
}

/* Derivative of the hybrid Pareto pdf with respect to xi, mu and sigma 
   divided by the pdf: dhdxi/h */
void hpdfgrad(double xi, double mu, double sigma, double x,
	      double *dhdtheta){
  double w, gamma, alpha, beta, dh2dgdxi, dh2dgdbeta, dh2dgdalpha;
  
  w = lambertw((1.0 + xi)*(1.0 + xi)/(2.0*PI));
  gamma = 3.0/2.0 + erf(sign(1.0+xi)*sqrt(w/2.0))/2.0;
  alpha = sign(1.0+xi) * sigma * sqrt(w) + mu;
  
  /* Derivative with respect to xi, mu and sigma through gamma: dhdgammadxi. */
  dhdtheta[0] = -sqrt(w)*sign(1.0+xi)/
    (gamma*sqrt(2.0*PI)*exp(w/2.0)*(1.0+w)*(1.0+xi));	     /* dhdxi */
    
  if (x <= alpha){
    /* Derivative only depends on the Normal pdf */
    dhdtheta[1] = (x-mu)/ (sigma*sigma);			    /* dhdmu */
    dhdtheta[2] = (x-mu)*(x-mu)/(sigma*sigma*sigma)-1.0/sigma; /* dhdsigma */
  }
  else {
    /* Derivative only depends on the GP pdf */
    beta = sigma*fabs(1.0+xi)/sqrt(w);
    /* Compute the derivative of the GP pdf w/r to xi, beta and alpha. */
    if (xi != 0.0){
      dh2dgdxi=(1.0+xi)*(alpha-x)/(xi*(beta+xi*(x-alpha))) +
	log(1.0+ xi*(x-alpha)/beta)/(xi*xi);
      dh2dgdbeta=(x-beta-alpha)/(beta*(beta+ xi*(x-alpha)));
      dh2dgdalpha=(1.0+ xi)/(beta + xi*(x-alpha));
    }
    else{
      dh2dgdxi=0.0;
      dh2dgdbeta=(x-beta-alpha)/(beta*beta);
      dh2dgdalpha=1.0/beta;
    }
    /* Derivative w/r to xi has dependencies through GP
       pdf (xi, beta, alpha) */
    dhdtheta[0]+=(dh2dgdxi + dh2dgdbeta*(sign(1.0+ xi)*sigma/sqrt(w)-
					 sigma*fabs(1.0+xi)/(sqrt(w)*(1.0+w)*(1.0+xi)))+
		  dh2dgdalpha*sigma*sqrt(w)/((1.0+w)*fabs(1.0+xi)));  /* dhdxi */
    dhdtheta[1]=dh2dgdalpha;			      		  /* dhdmu */
    dhdtheta[2]=(dh2dgdalpha*sign(1.0+ xi)*sqrt(w) +
		  dh2dgdbeta*fabs(1.0+ xi)/sqrt(w));	       /* dhdsigma */
  }
}



/* Hybrid Pareto mixture negative log-likelihood at n points*/
void ummhnll(double *theta, int *m, double *x, int *n, double *nll, 
	       double *nllgrad){
  int i,j,nparams=4**m-1;
  double *params, *act, *post, *dhdtheta;
  double lprob,ppost,pprior,lprior; 
  
  
  /* Allocate memory for pointers. */
  params = (double*) R_alloc(4**m, sizeof(double));
  act = (double*) R_alloc(*m, sizeof(double));
  post = (double*) R_alloc(*m, sizeof(double));
  dhdtheta = (double*) R_alloc(3, sizeof(double));
  
	
  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;
  for(j=0;j<nparams;j++)
    nllgrad[j]=0;


  /* Compute mixture parameters from vector theta */  
  ummhfwd(theta,m,params);

  for (i=0; i<*n; i++){
    /* Compute model log-probability for pattern i */
    /* mth component */
    act[*m-1] = hlogpdf(params[2**m-1], params[3**m-1], params[4**m-1],x[i]);

    if(theta[4**m-5]>0){
      pprior=-log(1+exp(-theta[4**m-5])); /* re-use computation */
      lprob=pprior+act[*m-1]; /* accumulates log(sum priors act) */
      post[*m-1]=lprob; /* log(prior) + log(act) */
      lprior=-theta[4**m-5]+pprior; /* accumulates sum log(1-g(theta_{4*j-1}) */
    }
    else{
      pprior=-log(1+exp(theta[4**m-5]));
      lprob=theta[4**m-5]+pprior+act[*m-1];
      post[*m-1]=lprob;
      lprior=pprior;
    }
    for(j=*m-2; j>-1; j--){
      act[j] = hlogpdf(params[*m+j], params[2**m+j], params[3**m+j],x[i]);

      if(j>0){
	if(theta[4*j-1]>0){
	  pprior=-log(1+exp(-theta[4*j-1]));
	  post[j]=pprior+lprior+act[j];
	  lprior+=-theta[4*j-1]+pprior;
	}
	else{
	  pprior=-log(1+exp(theta[4*j-1]));
	  post[j]=theta[4*j-1]+pprior+lprior+act[j];
	  lprior+=pprior;
	}
      }
      else 
	post[j]=lprior+act[j];
      
      if (lprob > post[j])
	lprob=lprob + log(1+exp(post[j]-lprob));
      else
	lprob=post[j] + log(1+exp(lprob-post[j]));            
    } /* for j */
 
    *nll -= lprob;
    
    /* Compute gradient for pattern i */
    ppost=0.0; /* Maintains \sum_{i=0}^{j-1} post[i] */
    pprior=0.0; /* Maintains \sum_{i=0}^{j-1} params[i] */
    for (j=0; j<*m; j++){
      /* Compute deltas for priors */
      post[j]=exp(post[j]-lprob); /* posterior */
      
      if(j > 0){ /* there are m-1 priors */
	ppost+=post[j-1];	
	nllgrad[4*j-1] += (params[j]*ppost-pprior*post[j])/
	  (pprior+params[j])*(1.0-2*MINPI);
      }
      pprior+=params[j];

      /* Gradient w/r to the hybrid Pareto parameters */
      hpdfgrad(params[*m+j], params[2**m+j],params[3**m+j],
	       x[i],dhdtheta);

      /* Compute gradients for hybrid Pareto parameters */
      nllgrad[4*j] += - post[j] *dhdtheta[0]*
	(1-exp(MINXI-params[*m+j]));   /* xi_j */
      /* deltas[4*j] = - post[j] *dhdtheta[0]*params[net->m+j]*
	 (1-params[net->m+j]);	 */
      nllgrad[4*j+1] += - post[j] * dhdtheta[1];    /* mu_j */
      nllgrad[4*j+2] += - post[j] * dhdtheta[2]
	*(1-exp(MINSIGMA-params[3**m+j]));	   /* sigma_j */
    } /* for j loop */
  } /* for i loop */
  
}



/* Hybrid Pareto mixture negative log-likelihood at n points with a penalty on the tail indexes based on a mixture (an exponential and a gaussian) with two modes: at zero and at 0.5 */
void ummhnll_bimodal_tailpen(double *theta, int *m, double * lambda, double *w,
			     double * beta, double *mu, double *sigma, double *x, int *n, 
			     double *nll,double *nllgrad){
			      			     
  int i,j,nparams=4**m-1;
  double *params, *act, *post, *dhdtheta;
  double lprob,ppost,pprior,lprior,tailpen,loga,logb; 
  
  
  /* Allocate memory for pointers. */
  params = (double*) R_alloc(4**m, sizeof(double));
  act = (double*) R_alloc(*m, sizeof(double));
  post = (double*) R_alloc(*m, sizeof(double));
  dhdtheta = (double*) R_alloc(3, sizeof(double));
	
  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;
  for(j=0;j<nparams;j++)
    nllgrad[j]=0;


  /* Compute mixture parameters from vector theta */  
  ummhfwd(theta,m,params);

  for (i=0; i<*n; i++){
    /* Compute model log-probability for pattern i */
    /* mth component */
    act[*m-1] = hlogpdf(params[2**m-1], params[3**m-1], params[4**m-1],x[i]);

    if(theta[4**m-5]>0){
      pprior=-log(1+exp(-theta[4**m-5])); /* re-use computation */
      lprob=pprior+act[*m-1]; /* accumulates log(sum priors act) */
      post[*m-1]=lprob; /* log(prior) + log(act) */
      lprior=-theta[4**m-5]+pprior; /* accumulates sum log(1-g(theta_{4*j-1}) */
    }
    else{
      pprior=-log(1+exp(theta[4**m-5]));
      lprob=theta[4**m-5]+pprior+act[*m-1];
      post[*m-1]=lprob;
      lprior=pprior;
    }
    for(j=*m-2; j>-1; j--){
      act[j] = hlogpdf(params[*m+j], params[2**m+j], params[3**m+j],x[i]);

      if(j>0){
	if(theta[4*j-1]>0){
	  pprior=-log(1+exp(-theta[4*j-1]));
	  post[j]=pprior+lprior+act[j];
	  lprior+=-theta[4*j-1]+pprior;
	}
	else{
	  pprior=-log(1+exp(theta[4*j-1]));
	  post[j]=theta[4*j-1]+pprior+lprior+act[j];
	  lprior+=pprior;
	}
      }
      else 
	post[j]=lprior+act[j];
      
      if (lprob > post[j])
	lprob=lprob + log(1+exp(post[j]-lprob));
      else
	lprob=post[j] + log(1+exp(lprob-post[j]));            
    } /* for j */
 
    *nll -= lprob;
    
    /* Compute gradient for pattern i */
    ppost=0.0; /* Maintains \sum_{i=0}^{j-1} post[i] */
    pprior=0.0; /* Maintains \sum_{i=0}^{j-1} params[i] */
    for (j=0; j<*m; j++){
      /* Compute deltas for priors */
      post[j]=exp(post[j]-lprob); /* posterior */

      /* Rprintf("priors[%d] = %e \n",j+1,params[j]);*/
      
      if(j > 0){ /* there are m-1 priors */
	ppost+=post[j-1];	
	nllgrad[4*j-1] += (params[j]*ppost-pprior*post[j])/
	  (pprior+params[j])*(1.0-2*MINPI);
      }
      pprior+=params[j];

      /* Gradient w/r to the hybrid Pareto parameters */
      hpdfgrad(params[*m+j], params[2**m+j],params[3**m+j],
	       x[i],dhdtheta);

      /* Compute gradients for hybrid Pareto parameters */     
      nllgrad[4*j] += - post[j] *dhdtheta[0]*
	(1-exp(MINXI-params[*m+j]));   /* xi_j */
      /* deltas[4*j] = - post[j] *dhdtheta[0]*params[net->m+j]*
	 (1-params[net->m+j]);	 */
      nllgrad[4*j+1] += - post[j] * dhdtheta[1];    /* mu_j */
      nllgrad[4*j+2] += - post[j] * dhdtheta[2]
	*(1-exp(MINSIGMA-params[3**m+j]));	   /* sigma_j */
    } /* for j loop */
  } /* for i loop */

  /* Compute the tail penalty */
  for (j=0;j<*m;j++){
    loga=log(*w) +log(*beta) - *beta*params[*m+j];
    logb=log(1-*w)-(params[*m+j]-*mu)*(params[*m+j]-*mu)/(2**sigma**sigma) - log(2*PI)/2.0 - log(*sigma);
    *nll -= *lambda*(loga + log(1 + exp(logb-loga)));
    tailpen= *beta + ((params[*m+j]-*mu)/(*sigma**sigma)-*beta)/(exp(loga-logb)+1);
    nllgrad[4*j] += *lambda*tailpen*(1-exp(MINXI-params[*m+j])); 
  }


}




/* --------------------- UMM : Quantile computations  ----------------------*/

/* Bracket the quantile of level q in interval [a,b] given an initial 
   interval and the cdf of UMM. */
int ummqbrack(double (*cdf)(double *,int, double),
	      double *params, int m, double q, double *a, double *b){
  int j;
  double fa,fb;
  
  fa= (*cdf)(params,m,*a)-q;
  fb=(*cdf)(params,m,*b)-q;

  /*Rprintf("fa = %f\n",fa);
    Rprintf("fb = %f\n",fb);*/

  for(j=0;j<NTRY;j++){
    if(fa*fb < 0.0) return 1;
    if(fabs(fa) < fabs(fb)){
      *a+=FACTOR*(*a-*b);
      fa= (*cdf)(params,m,*a)-q;
    }
    else{
      *b+=FACTOR*(*b-*a);
      fb=(*cdf)(params,m,*b)-q;
    }
  }
  return 0;
}


/* Quantile computation given a cdf and pdf of UMM:
   numerical approximation by bisection and Newton's methods,
   first bracketing interval given by [a,b] */
void ummquant(double (*cdf)(double *,int, double),
	       double (*pdf)(double *,int, double),
	       double *params, int m, double q, double *aa, double *bb,
	       double tol, int itmax, double *xq){ 
  
  int j, brack;
  double xqi, fa,fq,a,b;


  /* Bracket the quantile. */
  brack=ummqbrack(cdf,params, m, q, aa, bb);
  a=*aa;
  b=*bb;
  /*mexPrintf("a = %f\n",a); */
  /* mexPrintf("b = %f\n",b); */

   if (brack == 1){ /* bracketing succeeded */
    /* start a few iterations of the bisection method */
    j=0;
    fa= (*cdf)(params,m,a)-q;
    while(j<IBIS){ /* IBIS = 100*/
      *xq=a+(b-a)/2;
      fq=(*cdf)(params,m,*xq)-q;
      if((fq==0.0) | ((b-a)/2 < tol))
	break;
      j++;
      if(fa*fq>0){
	a=*xq;
	fa=fq;
      }
      else
	b=*xq;
      /* mexPrintf("xq = %f\n",*xq);*/
    } /* bisection loop */

    j=0;
    while(j<itmax){ /* Newton's method with initial value *xq */
      xqi=*xq-((*cdf)(params,m,*xq)-q)/(*pdf)(params,m,*xq);
      if(fabs(xqi-*xq)/fabs(xqi)<tol){
	*xq=xqi;
	return;
      }
      j++;
      *xq=xqi;
    } /* Newton loop */
   } /* successful bracket if */
   else{
     Rprintf("Bracketing failed\n");
     *xq=NAN;
   }
}

/* Gaussian mixture probability density function at n points*/
double ummgpdf(double *params, int m, double x){
  int j;
  double act,prob=0.0;

  for(j=0; j<m; j++){     
    act=normpdf(params[m+j],params[2*m+j],x);
    prob+=params[j]*act;
  }		

  return(prob);
}

/* Gaussian mixture pdf at n points for use in R*/
void ummgpdfR(double *params, int *m, double *x, int *n, double *p){
  int i;

  for(i=0;i<*n;i++)
    p[i]=ummgpdf(params,*m,x[i]);
}


/* Gaussian mixture cumulative distribution function */
double ummgcdf(double *params, int m, double x){
  int j;
  double act,prob=0.0;

  for(j=0; j<m; j++){     
    act=normcdf(params[m+j],params[2*m+j],x);
    prob+=params[j]*act;
  }		
  return(prob);
}

/* Gaussian mixture cdf at n points for use in R*/
void ummgcdfR(double *params, int *m, double *x, int *n, double *p){
  int i;

  for(i=0;i<*n;i++)
    p[i]=ummgcdf(params,*m,x[i]);
}



/* "Backward" pass in UMMG model */
void ummgbwd(double *params, int *m, double *theta){
  int i;
  double tmp,fprior=1.0;

  /* Apply inverse transfer functions to mixture parameters. */
  for(i=*m-1;i>-1;i--){
    if(i>0){
      tmp=params[i]/fprior;
      theta[3*i-1] = log(tmp/(1-tmp));
      fprior=fprior-params[i];
    }
    theta[3*i]=params[*m+i]; /* centers */
    theta[3*i+1]=softplusinv(params[2**m+i]-MINSIGMA); /* spread parameters */
  }

}


/* Create a CMM with hybrid Pareto components from R inputs */
void cmmhprint(double *params,int *d,int *h,int *m){
  CMM net;
  double **ptr,**endptr;
  int l;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4-1;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+4*net.m;
    }        
  }
  cmmprint(&net);

  /* Change the value of the CMM weights 
  int nparams=net.nout*(net.d+1)+net.h*(1+net.d+net.nout);
  for (l=0;l<nparams;l++)
    net.psi[l]=l;

    cmmprint(&net);*/

}

/* Print content of CMM structure */
  void cmmprint (CMM *net){
    int j;
    double **ptr,**endptr;
    double *ptrlin, *endptrlin;
    
    Rprintf("Dimension of input = %d\n", net->d);
    Rprintf("Number of hidden units = %d\n",net->h);
    Rprintf("Number of components = %d\n",net->m);
    if (net->s > 0)
      Rprintf("Number of reverse components = %d\n",net->s);
    Rprintf("Number of nn output = %d\n",net->nout);
    Rprintf("Linear weights :\n");
    
    endptrlin=net->psi+net->nout*(net->d+1);
    j=1;
    for(ptrlin=net->psi;ptrlin<endptrlin;ptrlin++){
      Rprintf("%d : %f\n",j,*ptrlin);
      j++;
    }
    
  if(net->h > 0){
    Rprintf("Non-linear weights :\n");
    endptr=net->theta+net->h;
    for(j=0;j<net->d+net->nout+1;j++){
      for(ptr=net->theta;ptr<endptr;ptr++)
	Rprintf(" %f ",(*ptr)[j]);
      Rprintf("\n");
    }
  }
}


/* Compute the output of a neural network with linear output activation
   function and a linear component. */
void nnlin(CMM *net, double *x, double *a, double *z){
  int j,k;
  double *ptrlin, *endptrlin;

 /* linear component: compute a, z is not assigned */
  endptrlin=net->psi+net->nout*(net->d+1);
  j=0;
  k=0; /* output */
  for(ptrlin=net->psi;ptrlin<endptrlin;ptrlin++){
    if (j==0)
      a[k]=*ptrlin;
    else
      a[k]+=*ptrlin*x[j-1];
    j++;
    if(j == net->d+1){
      j=0;
      k++;
    }
  }

  if (net->h > 0) { /* there is a non-linear component */
    double **ptr, **endptr;
    endptr=net->theta+net->h;
    k=0; /* hidden unit */
    for(ptr=net->theta;ptr<endptr;ptr++){
      endptrlin=*ptr+net->nout+1+net->d;
      j=0;
      for(ptrlin=*ptr;ptrlin<endptrlin;ptrlin++){
	if (j < net->d+1){ /* compute hidden unit activation */
	  if(j==0)
	    z[k]=(*ptr)[j];  /* add input bias to hidden unit k */
	  else
	    z[k]+=(*ptr)[j]*x[j-1];
	  if(j==net->d){
	    z[k]=tanh(z[k]);
	  }
	}
	else{ /* compute kth hidden unit contribution to output j */
	  a[j-net->d-1]+=(*ptr)[j]*z[k];
	}
	j++;
      } /* loop over weights related to hidden unit k */
      k++;
    } /* loop over hidden units */
  }
}


/* Forward propagation through a conditional mixture of hybrid Paretos:
   prior of the first component determined from the remaining priors. 
   The priors are constrained to lie in the interval [MINPI, 1-MINPI]
   to ease computations. */
void cmmhfwd(CMM *net, double *x, double *params, double *a, double *z){
  int i;
  double fprior=1.0;

  /* Propagate forward  through nn */
  nnlin(net, x, a, z);  

  /* Apply transfer functions to nn outputs. */
  for(i=net->m-1;i>-1;i--){
    if(i>0){
      /* priors in [0,1] */
      params[i]=(1.0/(1.0+exp(-a[4*i-1]))*(1.0-2.0*MINPI)+MINPI)*fprior;  
      fprior=fprior-params[i];
    }
    else
      params[i]=fprior;
    params[net->m+i]=softplus(a[4*i])+MINXI;   /* tail parameters */
    params[2*net->m+i]=a[4*i+1]; /* centers */
    params[3*net->m+i]=softplus(a[4*i+2])+MINSIGMA; /* spread parameters */
  }

}




/* Hybrid Pareto conditional mixture negative log-likelihood at n points*/
void cmmhnll(CMM *net, double *x, double *y, int n, 
	       double *nll, double *nllgrad){
  double *params, *a, *z, *act, *post, *deltas, *dhdtheta, *dldz;
  int i,j,k,nparams=net->nout*(net->d+1)+net->h*(1+net->d+net->nout);
  double lprob,ppost,pprior,lprior;
  double *ptrlin,*endptrlin;

  /* Allocate memory for pointers. */
  params = (double*) R_alloc(net->nout+1, sizeof(double));
  a = (double*) R_alloc(net->nout, sizeof(double));
  z = (double*) R_alloc(net->h, sizeof(double));
  act = (double*) R_alloc(net->m, sizeof(double));
  post = (double*) R_alloc(net->m, sizeof(double));
  deltas = (double*) R_alloc(net->nout, sizeof(double));
  dhdtheta = (double*) R_alloc(3, sizeof(double));
  dldz = (double*) R_alloc(net->h, sizeof(double));


  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;
  for(j=0;j<nparams;j++)
    nllgrad[j]=0;


  for (i=0; i<n; i++){
    /* Get nn outputs, hidden unit activations and mixture parameters for
     pattern i */
    cmmhfwd(net,x+i*net->d,params,a,z);

    /* Compute model probability for pattern i */
    /* mth component */
    act[net->m-1] = hlogpdf(params[2*net->m-1], params[3*net->m-1],
		     params[4*net->m-1],y[i]);
    if (net->m > 1){
      if(a[4*net->m-5]>0){
	pprior=-log(1+exp(-a[4*net->m-5])); /* re-use computation */
	lprob=pprior+act[net->m-1]; /* accumulates log(sum priors act) */
	post[net->m-1]=lprob; /* log(prior) + log(act) */
	lprior=-a[4*net->m-5]+pprior; /* accumulates sum log(1-g(a_{4*j-1}) */
      }
      else{
	pprior=-log(1+exp(a[4*net->m-5]));
	lprob=a[4*net->m-5]+pprior+act[net->m-1];
	post[net->m-1]=lprob;
	lprior=pprior;
      }
      for(j=net->m-2; j>-1; j--){
	act[j] = hlogpdf(params[net->m+j], params[2*net->m+j],
			 params[3*net->m+j],y[i]);
	if(j>0){
	  if(a[4*j-1]>0){
	    pprior=-log(1+exp(-a[4*j-1]));
	    post[j]=pprior+lprior+act[j];
	    lprior+=-a[4*j-1]+pprior;
	  }
	  else{
	    pprior=-log(1+exp(a[4*j-1]));
	    post[j]=a[4*j-1]+pprior+lprior+act[j];
	    lprior+=pprior;
	  }
	}
	else 
	  post[j]=lprior+act[j];
	
	if (lprob > post[j])
	  lprob=lprob + log(1+exp(post[j]-lprob));
	else
	  lprob=post[j] + log(1+exp(lprob-post[j]));            
      } /* for j */

    }
    else{
      post[0] = act[0];
      lprob=act[0];
    }
    
    
    *nll -= lprob;


    /* Compute gradient for pattern i */
    ppost=0.0; /* Maintains \sum_{i=0}^{j-1} post[i] */
    pprior=0.0; /* Maintains \sum_{i=0}^{j-1} params[i] */
    for (j=0; j<net->m; j++){
      /* Compute deltas for priors */
      post[j]=exp(post[j]-lprob); /* posterior */
      
      if(j > 0){ /* there are m-1 priors */
	ppost+=post[j-1];	
	deltas[4*j-1] = (params[j]*ppost-pprior*post[j])/
	  (pprior+params[j])*(1.0-2.0*MINPI);
      }
      pprior+=params[j];
      
      /* Gradient w/r to the hybrid Pareto parameters */
      hpdfgrad(params[net->m+j], params[2*net->m+j],
	       params[3*net->m+j], y[i],dhdtheta);
      
      /* Compute deltas for hybrid Pareto parameters */
      deltas[4*j] = - post[j] *dhdtheta[0]*
	(1-exp(MINXI-params[net->m+j]));   /* xi_j */
      deltas[4*j+1] = - post[j] * dhdtheta[1];    /* mu_j */
      deltas[4*j+2] = - post[j] * dhdtheta[2]
	*(1-exp(MINSIGMA-params[3*net->m+j]));	   /* sigma_j */
    } /* for j loop */
    
    
      /* Compute gradient with respect to nn parameters */
    endptrlin=nllgrad+net->nout*(net->d+1);
    j=0;
    k=0; /* output */
    for(ptrlin=nllgrad;ptrlin<endptrlin;ptrlin++){
      if (j==0) /* bias */
	*ptrlin+=deltas[k];
      else
	*ptrlin+=deltas[k]*x[i*net->d+j-1];
      j++;
      if (j == net->d+1){
	j=0;
	k++;
      }
    }
      
      
    if (net->h > 0) {
      double **ptr, **endptr=net->theta+net->h;
      k=0;  /* hidden unit */
      int l = net->nout*(net->d+1); /* # of linear parameters */
      /* Derivative w/r to output layer weights */
      for(ptr=net->theta;ptr<endptr;ptr++){
	endptrlin=*ptr+net->d+4*net->m;
	j=0; /* output */
	dldz[k]=0;
	for(ptrlin=*ptr+net->d+1;ptrlin<endptrlin;ptrlin++){
	  nllgrad[l+k*(4*net->m+net->d)+j+net->d+1]+=
	    deltas[j] * z[k];
	  dldz[k]+=deltas[j]*(*ptrlin);
	  j++;
	}
	k++;
      }
	
      /* Derivative w/r to hidden layer weights */
      for(k=0;k<net->h;k++){
	nllgrad[l+k*(net->d+4*net->m)]+=(1-z[k]*z[k])*dldz[k];
	for(j=0;j<net->d;j++)
	  nllgrad[l+k*(net->d+4*net->m)+1+j]+=(1-z[k]*z[k])*dldz[k]*
	    x[i*net->d+j];
      }
    }
  }  /* for i loop */
  
}

/* Create a CMM with hybrid Pareto components from R inputs 
   and call cmmhnll() */
void cmmhnllR(double *params,int *d,int *h,int *m, double *x, 
	      double *y, int *n, double *nll, double *nllgrad){
  CMM net;
  double **ptr,**endptr;
  int l;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4-1;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+4*net.m;
    }        
  }
  cmmhnll(&net,x,y,*n,nll,nllgrad);
}




/* Create a CMM with hybrid Pareto components from R inputs 
   and call cmmhfwd() */
void cmmhfwdR(double *params,int *d,int *h,int *m, double *x, 
	      int *n, double *params_mixt, double *a, double *z){
  CMM net;
  double **ptr,**endptr;
  int l,i;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4-1;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+4*net.m;
    }        
  }
  for(i=0;i<*n;i++)
    cmmhfwd(&net,x+i*net.d,params_mixt+i*4*net.m,a+i*net.nout,z+i*net.h);
}



/* Hybrid Pareto conditional mixture negative log-likelihood at n points*/
/* Bimodal penalty mixture (an exponential and a gaussian) with two modes: 
   at zero and at mu (must be set)  */
void cmmhnll_bimodal_tailpen(CMM *net, double *x, double *y, int n, 
			     double *lambda,double *w, double *beta, double *mu, 
			     double *sigma, double *nll, double *nllgrad){
  
  double *params, *a, *z, *act, *post, *deltas, *dhdtheta, *dldz;
  int i,j,k,nparams=net->nout*(net->d+1)+net->h*(1+net->d+net->nout);
  double lprob,ppost,pprior,lprior,tailpen,loga,logb;
  double *ptrlin,*endptrlin;

  /* Allocate memory for pointers. */
  params = (double*) R_alloc(net->nout+1, sizeof(double));
  a = (double*) R_alloc(net->nout, sizeof(double));
  z = (double*) R_alloc(net->h, sizeof(double));
  act = (double*) R_alloc(net->m, sizeof(double));
  post = (double*) R_alloc(net->m, sizeof(double));
  deltas = (double*) R_alloc(net->nout, sizeof(double));
  dhdtheta = (double*) R_alloc(3, sizeof(double));
  dldz = (double*) R_alloc(net->h, sizeof(double));

  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;
  for(j=0;j<nparams;j++)
    nllgrad[j]=0;

  *lambda=*lambda/(n*net->m); /* Divide the tail penalty by the number of observations  */

  for (i=0; i<n; i++){
    /* Get nn outputs, hidden unit activations and mixture parameters for
     pattern i */
    cmmhfwd(net,x+i*net->d,params,a,z);

    /* Compute model probability for pattern i */
    /* mth component */
    act[net->m-1] = hlogpdf(params[2*net->m-1], params[3*net->m-1],
		     params[4*net->m-1],y[i]);

    if (net->m > 1){
      if(a[4*net->m-5]>0){
	pprior=-log(1+exp(-a[4*net->m-5])); /* re-use computation */
	lprob=pprior+act[net->m-1]; /* accumulates log(sum priors act) */
	post[net->m-1]=lprob; /* log(prior) + log(act) */
	lprior=-a[4*net->m-5]+pprior; /* accumulates sum log(1-g(a_{4*j-1}) */
      }
      else{
	pprior=-log(1+exp(a[4*net->m-5]));
	lprob=a[4*net->m-5]+pprior+act[net->m-1];
	post[net->m-1]=lprob;
	lprior=pprior;
      }
      for(j=net->m-2; j>-1; j--){
	act[j] = hlogpdf(params[net->m+j], params[2*net->m+j],
			 params[3*net->m+j],y[i]);

	if(j>0){
	  if(a[4*j-1]>0){
	    pprior=-log(1+exp(-a[4*j-1]));
	    post[j]=pprior+lprior+act[j];
	    lprior+=-a[4*j-1]+pprior;
	  }
	  else{
	    pprior=-log(1+exp(a[4*j-1]));
	    post[j]=a[4*j-1]+pprior+lprior+act[j];
	    lprior+=pprior;
	  }
	}
	else 
	  post[j]=lprior+act[j];
	
	if (lprob > post[j])
	  lprob=lprob + log(1+exp(post[j]-lprob));
	else
	  lprob=post[j] + log(1+exp(lprob-post[j]));            
      } /* for j */

    }
    else{
      post[0] = act[0];
      lprob=act[0];
    }
     
    *nll -= lprob;

    /* Rprintf("nll[%d] = %f\n",i+1,*nll); */

    /* Compute gradient for pattern i */
    ppost=0.0; /* Maintains \sum_{i=0}^{j-1} post[i] */
    pprior=0.0; /* Maintains \sum_{i=0}^{j-1} params[i] */

    for (j=0; j<net->m; j++){
      /* Compute deltas for priors */
      post[j]=exp(post[j]-lprob); /* posterior */
      
      if(j > 0){ /* there are m-1 priors */
	ppost+=post[j-1];	
	deltas[4*j-1] = (params[j]*ppost-pprior*post[j])/
	  (pprior+params[j])*(1.0-2*MINPI);
      }
      pprior+=params[j];

      /* Compute tail penalty together with gradient */
      loga=log(*w) +log(*beta) - *beta*params[net->m+j];
      logb=log(1.0-*w)-(params[net->m+j]-*mu)*(params[net->m+j]-*mu)/(2.0**sigma**sigma) - M_LN_SQRT_2PI - log(*sigma);
      *nll -= *lambda*(loga + log(1.0 + exp(logb-loga)));

      tailpen= *beta + ((params[net->m+j]-*mu)/(*sigma**sigma)-*beta)/(exp(loga-logb)+1.0);

      /* Gradient w/r to the hybrid Pareto parameters */
      hpdfgrad(params[net->m+j], params[2*net->m+j],
	       params[3*net->m+j], y[i],dhdtheta);
      
      /* Compute deltas for hybrid Pareto parameters */
      deltas[4*j] = (*lambda*tailpen - post[j] *dhdtheta[0])*
	(1-exp(MINXI-params[net->m+j]));   /* xi_j */
      deltas[4*j+1] = - post[j] * dhdtheta[1];    /* mu_j */
      deltas[4*j+2] = - post[j] * dhdtheta[2]
	*(1-exp(MINSIGMA-params[3*net->m+j]));	   /* sigma_j */
    } /* for j loop */
    
    
      /* Compute gradient with respect to nn parameters */
    endptrlin=nllgrad+net->nout*(net->d+1);
    j=0;
    k=0; /* output */
    for(ptrlin=nllgrad;ptrlin<endptrlin;ptrlin++){
      if (j==0) /* bias */
	*ptrlin+=deltas[k];
      else
	*ptrlin+=deltas[k]*x[i*net->d+j-1];
      j++;
      if (j == net->d+1){
	j=0;
	k++;
      }
    }
      
      
    if (net->h > 0) {
      double **ptr, **endptr=net->theta+net->h;
      k=0;  /* hidden unit */
      int l = net->nout*(net->d+1); /* # of linear parameters */
      /* Derivative w/r to output layer weights */
      for(ptr=net->theta;ptr<endptr;ptr++){
	endptrlin=*ptr+net->d+4*net->m;
	j=0; /* output */
	dldz[k]=0;
	for(ptrlin=*ptr+net->d+1;ptrlin<endptrlin;ptrlin++){
	  nllgrad[l+k*(4*net->m+net->d)+j+net->d+1]+=
	    deltas[j] * z[k];
	  dldz[k]+=deltas[j]*(*ptrlin);
	  j++;
	}
	k++;
      }
	
      /* Derivative w/r to hidden layer weights */
      for(k=0;k<net->h;k++){
	nllgrad[l+k*(net->d+4*net->m)]+=(1-z[k]*z[k])*dldz[k];
	for(j=0;j<net->d;j++)
	  nllgrad[l+k*(net->d+4*net->m)+1+j]+=(1-z[k]*z[k])*dldz[k]*
	    x[i*net->d+j];
      }
    }
  }  /* for i loop */
  
}

/* Create a CMM with hybrid Pareto components from R inputs 
   and call cmmhnll_bimodal_tailpen() */
void cmmhnll_bimodal_tailpenR(double *params,int *d,int *h,int *m, double *x, 
			      double *y, int *n, double *lambda, double *w, 
			      double *beta, double *mu, double *sigma, double *nll, 
			      double *nllgrad){
  CMM net;
  double **ptr,**endptr;
  int l;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4-1;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+4*net.m;
    }        
  }
  cmmhnll_bimodal_tailpen(&net,x,y,*n,lambda,w,beta,mu,sigma,nll,nllgrad);
}


/* Forward propagation through a conditional mixture of hybrid Paretos:
   prior of the first component determined from the remaining priors. 
   The priors are constrained to lie in the interval [MINPI, 1-MINPI]
   to ease computations. The first input is the probability of the Dirac */
void cmmhfwd_dirac(CMM *net, double *x, double *params, double *a, double *z){
  int i;
  double fprior=1.0;

  /* Propagate forward  through nn */
  nnlin(net, x, a, z);  

  /* First output is the Dirac probability (nout = 4m) */
  params[0] = 1.0/(1.0+exp(-a[0]));

  /* Apply transfer functions to nn outputs. */
  for(i=net->m;i>0;i--){
    if(i>1){
      /* priors in [0,1] */
      params[i]=(1.0/(1.0+exp(-a[4*(i-1)]))*(1.0-2.0*MINPI)+MINPI)*fprior;  
      fprior=fprior-params[i];
    }
    else
      params[i]=fprior;
    params[net->m+i]=softplus(a[4*(i-1)+1])+MINXI;   /* tail parameters */
    params[2*net->m+i]=a[4*(i-1)+2]; /* centers */
    params[3*net->m+i]=softplus(a[4*(i-1)+3])+MINSIGMA; /* spread parameters */
  }

}

/* Create a CMM with hybrid Pareto components from R inputs 
   and call cmmhfwd_dirac() */
void cmmhfwd_diracR(double *params,int *d,int *h,int *m, double *x, 
	      int *n, double *params_mixt, double *a, double *z){
  CMM net;
  double **ptr,**endptr;
  int l,i;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }
  for(i=0;i<*n;i++)
    cmmhfwd_dirac(&net,x+i*net.d,params_mixt+i*(4*net.m+1),a+i*net.nout,z+i*net.h);
}




/* Hybrid Pareto conditional mixture negative log-likelihood at n points*/
/* Bimodal penalty mixture (an exponential and a gaussian) with two modes: 
   at zero and at 0.2 and delta of Dirac probability of null observation  */
void cmmhnll_bimodal_tailpen_dirac(CMM *net, double *x, double *y, int n, 
			     double *lambda,double *w, double *beta, double *mu, 
			     double *sigma, double *nll, double *nllgrad){


  double *params, *a, *z, *act, *post, *deltas, *dhdtheta, *dldz;
  int i,j,k,nparams=net->nout*(net->d+1)+net->h*(1+net->d+net->nout);
  double lprob,ppost,pprior,lprior,tailpen,loga,logb;
  double *ptrlin,*endptrlin;

  /* Allocate memory for pointers. */
  params = (double*) R_alloc(net->nout+1, sizeof(double));
  a = (double*) R_alloc(net->nout, sizeof(double));
  z = (double*) R_alloc(net->h, sizeof(double));
  act = (double*) R_alloc(net->m, sizeof(double));
  post = (double*) R_alloc(net->m, sizeof(double));
  deltas = (double*) R_alloc(net->nout, sizeof(double));
  dhdtheta = (double*) R_alloc(3, sizeof(double));
  dldz = (double*) R_alloc(net->h, sizeof(double));

  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;
  for(j=0;j<nparams;j++)
    nllgrad[j]=0;

  *lambda=*lambda/(n*net->m); /* Divide the tail penalty by the number of observations  */

  for (i=0; i<n; i++){
    /* Get nn outputs, hidden unit activations and mixture parameters for
     pattern i */
    cmmhfwd_dirac(net,x+i*net->d,params,a,z);

    if (y[i] > 0){ /* positive obs. (it rained that day) */
      /* Compute model probability for pattern i */
      /* mth component */
      act[net->m-1] = hlogpdf(params[2*net->m], params[3*net->m],
			      params[4*net->m],y[i]);
      /* Rprintf("xi = %f, mu = %f, sigma = %f \n",params[2*net->m], params[3*net->m],params[4*net->m]); */

      if (net->m > 1){ /* more than one component */
	if(a[4*net->m-4]>0){
	  pprior=-log(1+exp(-a[4*net->m-4])); /* re-use computation */
	  lprob=pprior+act[net->m-1]; /* accumulates log(sum priors act) */
	  post[net->m-1]=lprob; /* log(prior) + log(act) */
	  lprior=-a[4*net->m-4]+pprior; /* accumulates sum log(1-g(a_{4*j-1}) */
	}
	else{
	  pprior=-log(1+exp(a[4*net->m-4]));
	  lprob=a[4*net->m-4]+pprior+act[net->m-1];
	  post[net->m-1]=lprob;
	  lprior=pprior;
	}
	for(j=net->m-2; j>-1; j--){
	  act[j] = hlogpdf(params[net->m+j+1], params[2*net->m+j+1],
			   params[3*net->m+j+1],y[i]);
	  
	  if(j>0){
	    if(a[4*j]>0){
	      pprior=-log(1+exp(-a[4*j]));
	      post[j]=pprior+lprior+act[j];
	      lprior+=-a[4*j]+pprior;
	    }
	    else{
	      pprior=-log(1+exp(a[4*j]));
	      post[j]=a[4*j]+pprior+lprior+act[j];
	      lprior+=pprior;
	    }
	  }
	  else 
	    post[j]=lprior+act[j];
	  
	  if (lprob > post[j])
	    lprob=lprob + log(1+exp(post[j]-lprob));
	  else
	    lprob=post[j] + log(1+exp(lprob-post[j]));            
	} /* for j */
	
      }
      else{ /* just one component */
	post[0] = act[0];
	lprob=act[0];
      }
      
      *nll -= (lprob + log(params[0]));
      
      /* Compute gradient for pattern i */
      ppost=0.0; /* Maintains \sum_{i=0}^{j-1} post[i] */
      pprior=0.0; /* Maintains \sum_{i=0}^{j-1} params[i] */
      
      for (j=0; j<net->m; j++){
	/* Compute deltas for priors */
	post[j]=exp(post[j]-lprob); /* posterior */
	
	if(j > 0){ /* there are m-1 priors */
	  ppost+=post[j-1];	
	  deltas[4*j] = (params[j+1]*ppost-pprior*post[j])/
	    (pprior+params[j+1])*(1.0-2*MINPI);
	}
	pprior+=params[j+1];
	
	/* Compute tail penalty together with gradient */
	loga=log(*w) +log(*beta) - *beta*params[net->m+j+1];
	logb=log(1-*w)-(params[net->m+j+1]-*mu)*(params[net->m+j+1]-*mu)/(2**sigma**sigma) - M_LN_SQRT_2PI - log(*sigma);
	*nll -= *lambda*(loga + log(1 + exp(logb-loga)));

	tailpen= *beta + ((params[net->m+j+1]-*mu)/(*sigma**sigma)-*beta)/(exp(loga-logb)+1);
	
	/* Gradient w/r to the hybrid Pareto parameters */
	hpdfgrad(params[net->m+j+1], params[2*net->m+j+1],
		 params[3*net->m+j+1], y[i],dhdtheta);
	
	/* Compute deltas for hybrid Pareto parameters */
	deltas[4*j+1] = (*lambda*tailpen - post[j] *dhdtheta[0])*
	  (1.0-exp(MINXI-params[net->m+j+1]));   /* xi_j */
	deltas[4*j+2] = - post[j] * dhdtheta[1];    /* mu_j */
	deltas[4*j+3] = - post[j] * dhdtheta[2]
	  *(1.0-exp(MINSIGMA-params[3*net->m+j+1]));	   /* sigma_j */
      } /* for j loop */

      /* Delta for the Dirac parameter */
      deltas[0] = params[0]-1.0;
    } /* positive observation condition */
    else { /* null observation (Dirac probability) */
      *nll -= log(1.0-params[0]); 
      deltas[0] = params[0];
      for (j=1; j<net->nout; j++)
	deltas[j] = 0.0;
    }
    
    /* Rprintf("nll[%d] = %f\n",i+1,*nll); */
       
      /* Compute gradient with respect to nn parameters */
    endptrlin=nllgrad+net->nout*(net->d+1);
    j=0;
    k=0; /* output */
    for(ptrlin=nllgrad;ptrlin<endptrlin;ptrlin++){
      if (j==0) /* bias */
	*ptrlin+=deltas[k];
      else
	*ptrlin+=deltas[k]*x[i*net->d+j-1];
      j++;
      if (j == net->d+1){
	j=0;
	k++;
      }
    }
      
      
    if (net->h > 0) {
      double **ptr, **endptr=net->theta+net->h;
      k=0;  /* hidden unit */
      int l = net->nout*(net->d+1); /* # of linear parameters */
      /* Derivative w/r to output layer weights */
      for(ptr=net->theta;ptr<endptr;ptr++){
	endptrlin=*ptr+net->d+net->nout+1;
	j=0; /* output */
	dldz[k]=0;
	for(ptrlin=*ptr+net->d+1;ptrlin<endptrlin;ptrlin++){
	  nllgrad[l+k*(net->nout+1+net->d)+j+net->d+1]+=
	    deltas[j] * z[k];
	  dldz[k]+=deltas[j]*(*ptrlin);
	  j++;
	}
	k++;
      }
	
      /* Derivative w/r to hidden layer weights */
      for(k=0;k<net->h;k++){
	nllgrad[l+k*(net->d+net->nout+1)]+=(1-z[k]*z[k])*dldz[k];
	for(j=0;j<net->d;j++)
	  nllgrad[l+k*(net->d+net->nout+1)+1+j]+=(1-z[k]*z[k])*dldz[k]*
	    x[i*net->d+j];
      }
    }
  }  /* for i loop */
  
}

/* Create a CMM with hybrid Pareto components from R inputs 
   and call cmmhnll_bimodal_tailpen_dirac() */
void cmmhnll_bimodal_tailpen_diracR(double *params,int *d,int *h,int *m, double *x, 
			      double *y, int *n, double *lambda, double *w, 
			      double *beta, double *mu, double *sigma, double *nll, 
			      double *nllgrad){
  CMM net;
  double **ptr,**endptr;
  int l;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }
  cmmhnll_bimodal_tailpen_dirac(&net,x,y,*n,lambda,w,beta,mu,sigma,nll,nllgrad);
}



/* Hybrid Pareto conditional mixture negative log-likelihood at n points*/
/* Delta of Dirac probability of null observation (dont compute gradient) */
void cmmhnll_dirac(CMM *net, double *x, double *y, int n, double *nll){


  double *params, *a, *z, *act, *post;
  int i,j;
  double lprob,pprior,lprior;

  /* Allocate memory for pointers. */
  params = (double*) R_alloc(net->nout+1, sizeof(double));
  a = (double*) R_alloc(net->nout, sizeof(double));
  z = (double*) R_alloc(net->h, sizeof(double));
  act = (double*) R_alloc(net->m, sizeof(double));
  post = (double*) R_alloc(net->m, sizeof(double)); 

  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;

  for (i=0; i<n; i++){
    /* Get nn outputs, hidden unit activations and mixture parameters for
     pattern i */
    cmmhfwd_dirac(net,x+i*net->d,params,a,z);

    if (y[i] > 0){ /* positive obs. (it rained that day) */
      /* Compute model probability for pattern i */
      /* mth component */
      act[net->m-1] = hlogpdf(params[2*net->m], params[3*net->m],
			      params[4*net->m],y[i]);
      /* Rprintf("xi = %f, mu = %f, sigma = %f \n",params[2*net->m], params[3*net->m],params[4*net->m]); */

      if (net->m > 1){ /* more than one component */
	if(a[4*net->m-4]>0){
	  pprior=-log(1+exp(-a[4*net->m-4])); /* re-use computation */
	  lprob=pprior+act[net->m-1]; /* accumulates log(sum priors act) */
	  post[net->m-1]=lprob; /* log(prior) + log(act) */
	  lprior=-a[4*net->m-4]+pprior; /* accumulates sum log(1-g(a_{4*j-1}) */
	}
	else{
	  pprior=-log(1+exp(a[4*net->m-4]));
	  lprob=a[4*net->m-4]+pprior+act[net->m-1];
	  post[net->m-1]=lprob;
	  lprior=pprior;
	}
	for(j=net->m-2; j>-1; j--){
	  act[j] = hlogpdf(params[net->m+j+1], params[2*net->m+j+1],
			   params[3*net->m+j+1],y[i]);
	  
	  if(j>0){
	    if(a[4*j]>0){
	      pprior=-log(1+exp(-a[4*j]));
	      post[j]=pprior+lprior+act[j];
	      lprior+=-a[4*j]+pprior;
	    }
	    else{
	      pprior=-log(1+exp(a[4*j]));
	      post[j]=a[4*j]+pprior+lprior+act[j];
	      lprior+=pprior;
	    }
	  }
	  else 
	    post[j]=lprior+act[j];
	  
	  if (lprob > post[j])
	    lprob=lprob + log(1+exp(post[j]-lprob));
	  else
	    lprob=post[j] + log(1+exp(lprob-post[j]));            
	} /* for j */
	
      }
      else{ /* just one component */
	post[0] = act[0];
	lprob=act[0];
      }
      
      *nll -= (lprob + log(params[0]));
                 
    } /* positive observation condition */
    else { /* null observation (Dirac probability) */
      *nll -= log(1.0-params[0]); 
    }
    
  }  /* for i loop */
  
}


/* Create a CMM with hybrid Pareto components from R inputs 
   and call cmmhnll_dirac() */
void cmmhnll_diracR(double *params,int *d,int *h,int *m, double *x, 
			      double *y, int *n, double *nll){
  CMM net;
  double **ptr,**endptr;
  int l;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }
  cmmhnll_dirac(&net,x,y,*n,nll);
}



/* Create a CMM with  hybrid Pareto components from R inputs 
   and perform quantile computations */
void cmmhquant(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq){
  CMM net;
  double **ptr,**endptr;
  int l,i,j;
  double tol=10^(-16); /* Precision on estimate */
  int itmax=0; /* No Newton iterations */

  net.psi=theta;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4-1;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+4*net.m;
    }        
  }

  double *params = (double *) R_alloc(4*net.m,sizeof(double));
  double *act = (double *) R_alloc(net.nout,sizeof(double));
  double *z = (double *) R_alloc(net.h,sizeof(double));


  for (i=0;i<*n;i++){
    cmmhfwd(&net,x+i*net.d,params,act,z);
    for(j=0;j<*nq;j++)
      /*mexPrintf("q[%d] = %f\n",j+1,q[j]); */
      ummquant(ummhcdf,ummhpdf,params,net.m,q[j],a,b,tol,itmax,xq+j+i**nq); 
  }

}


/* Create a CMM with hybrid Pareto components truncated below zero
   from R inputs and perform quantile computations */
void cmmhquant_trunc(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq){
  CMM net;
  double **ptr,**endptr;
  int l,i,j;
  double tol=10^(-16); /* Precision on estimate */
  int itmax=0; /* No Newton iterations */
  double F0; /* Probability of being smaller than 0 with CMMH */

  net.psi=theta;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4-1;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+4*net.m;
    }        
  }

  double *params = (double *) R_alloc(4*net.m,sizeof(double));
  double *act = (double *) R_alloc(net.nout,sizeof(double));
  double *z = (double *) R_alloc(net.h,sizeof(double));


  for (i=0;i<*n;i++){
    cmmhfwd(&net,x+i*net.d,params,act,z);
    F0=ummhcdf(params,net.m,0.0);
    for(j=0;j<*nq;j++)
      /*mexPrintf("q[%d] = %f\n",j+1,q[j]); */
      ummquant(ummhcdf,ummhpdf,params,net.m,q[j]*(1.0-F0),a,b,tol,itmax,xq+j+i**nq); 
  }

}


/* Create a CMM with hybrid Pareto components truncated below zero
   and one Dirac at zero from R inputs and perform quantile computations */
void cmmhquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq){
  CMM net;
  double **ptr,**endptr;
  int l,i,j;
  double tol=10^(-16); /* Precision on estimate */
  int itmax=0; /* No Newton iterations */
  double F0; /* Probability of being smaller than 0 with CMMH */

  net.psi=theta;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }

  double *params = (double *) R_alloc(net.nout+1,sizeof(double));
  double *params_cmmh = (double *) R_alloc(net.nout,sizeof(double));
  double *out = (double *) R_alloc(net.nout,sizeof(double));
  double *z = (double *) R_alloc(net.h,sizeof(double));


  for (i=0;i<*n;i++){
    cmmhfwd_dirac(&net,x+i*net.d,params,out,z);
    params_cmmh = params+1;
    for(j=0;j<*nq;j++){
      if (q[j] > 1.0-params[0]){
	F0=ummhcdf(params_cmmh,net.m,0.0);
	/*Rprintf("F0 = %f\n",F0);*/
	ummquant(ummhcdf,ummhpdf,params_cmmh,net.m,((q[j]-1.0)/params[0]+1.0)*(1.0-F0)+F0,a,b,tol,itmax,xq+j+i**nq); 
      } else
	xq[j+i**nq]=0.0;
    }
  }

}


/* Create a CMM with hybrid Pareto components truncated below zero and 
   one Dirac at zero from R inputs and perform conditional quantile 
   computations */
void cmmhcquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq){
  CMM net;
  double **ptr,**endptr;
  int l,i,j;
  double tol=10^(-16); /* Precision on estimate */
  int itmax=0; /* No Newton iterations */
  double F0; /* Probability of being smaller than 0 with CMMH */

  net.psi=theta;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*4;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }

  double *params = (double *) R_alloc(net.nout+1,sizeof(double));
  double *params_cmmh = (double *) R_alloc(net.nout,sizeof(double));
  double *out = (double *) R_alloc(net.nout,sizeof(double));
  double *z = (double *) R_alloc(net.h,sizeof(double));


  for (i=0;i<*n;i++){
    cmmhfwd_dirac(&net,x+i*net.d,params,out,z);
    params_cmmh = params+1;
    F0=ummhcdf(params_cmmh,net.m,0.0);
    for(j=0;j<*nq;j++){
      ummquant(ummhcdf,ummhpdf,params_cmmh,net.m,q[j]*(1.0-F0)+F0,a,b,tol,itmax,xq+j+i**nq); 
    }

  }
}



/*-----------------------------------------------------------------------*/
/*----------- Conditional mixture model with Gaussian components---------*/
/*-----------------------------------------------------------------------*/


/* Forward propagation through a conditional mixture Gaussian components*/
void cmmgfwd(CMM *net, double *x, double *params, double *a, double *z){
  int i;
  double fprior=1.0;

  /* Propagate forward  through nn */
  nnlin(net, x, a, z);  

  /* Apply transfer functions to nn outputs. */
  for(i=net->m-1;i>-1;i--){
    if(i>0){
      /* priors in [0,1] */
      params[i]=(1.0/(1.0+exp(-a[3*i-1]))*(1.0-2.0*MINPI)+MINPI)*fprior;  
      fprior=fprior-params[i];
    }
    else
      params[i]=fprior;
    params[net->m+i]=a[3*i];  /* means */
    params[2*net->m+i]=softplus(a[3*i+1])+MINSIGMA; /* standard deviations */
  }

}



/* Gaussian conditional mixture negative log-likelihood at n points*/
void cmmgnll(CMM *net, double *x, double *y, int n, 
	       double *nll, double *nllgrad){
  double *params, *a, *z, *act, *post, *deltas, *dldz;
  int i,j,k,nparams=net->nout*(net->d+1)+net->h*(1+net->d+net->nout);
  double lprob,ppost,pprior,lprior;
  double *ptrlin,*endptrlin;

  /* Allocate memory for pointers. */
  params = (double*) R_alloc(net->nout+1, sizeof(double));
  a = (double*) R_alloc(net->nout, sizeof(double));
  z = (double*) R_alloc(net->h, sizeof(double));
  act = (double*) R_alloc(net->m, sizeof(double));
  post = (double*) R_alloc(net->m, sizeof(double));
  deltas = (double*) R_alloc(net->nout, sizeof(double));
  dldz = (double*) R_alloc(net->h, sizeof(double));
  

  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;
  for(j=0;j<nparams;j++)
    nllgrad[j]=0; 

  for (i=0; i<n; i++){
    /* Get nn outputs, hidden unit activations and mixture parameters for
     pattern i */
    cmmgfwd(net,x+i*net->d,params,a,z);

    /* Compute model probability for pattern i */
    /* mth component */
    act[net->m-1] = normlogpdf(params[2*net->m-1], params[3*net->m-1], y[i]);
    if(a[3*net->m-4]>0){
      pprior=-log(1+exp(-a[3*net->m-4])); /* re-use computation */
      post[net->m-1]=pprior+act[net->m-1];  /* log(prior) + log(act) */
      lprob=post[net->m-1]; /* accumulates log(sum priors act) */
      lprior=-a[3*net->m-4]+pprior; /* accumulates sum log(1-g(a_{4*j-1}) */
    }
    else{
      pprior=-log(1+exp(a[3*net->m-4]));
      post[net->m-1]=a[3*net->m-4]+pprior+act[net->m-1];
      lprob=post[net->m-1];
      lprior=pprior;
    }
    for(j=net->m-2; j>-1; j--){
      act[j] = normlogpdf(params[net->m+j], params[2*net->m+j],y[i]);
      if(j>0){
	if(a[3*j-1]>0){
	  pprior=-log(1+exp(-a[3*j-1]));
	  post[j]=pprior+lprior+act[j];
	  lprior+=-a[3*j-1]+pprior;
	}
	else{
	  pprior=-log(1+exp(a[3*j-1]));
	  post[j]=a[3*j-1]+pprior+lprior+act[j];
	  lprior+=pprior;
	}
      }
      else 
	post[j]=lprior+act[j];
      
      if (lprob > post[j])
	lprob=lprob + log(1+exp(post[j]-lprob));
      else
	lprob=post[j] + log(1+exp(lprob-post[j]));            
    } /* for j */
    
    *nll -= lprob;
    /* Compute gradient for pattern i */
    ppost=0.0; /* Maintains \sum_{i=0}^{j-1} post[i] */
    pprior=0.0; /* Maintains \sum_{i=0}^{j-1} params[i] */
    for (j=0; j<net->m; j++){
      /* Compute deltas for priors */
      post[j]=exp(post[j]-lprob); /* posterior */
      
      if(j > 0){ /* there are m-1 priors */
	ppost+=post[j-1];	
	deltas[3*j-1] = (params[j]*ppost-pprior*post[j])/
	  (pprior+params[j])*(1.0-2.0*MINPI);
      }
      pprior+=params[j];
            
      /* Compute deltas for the Gaussian parameters */
      deltas[3*j] = - post[j] * /* mu_j */
	(y[i]-params[net->m+j])/(params[2*net->m+j]*params[2*net->m+j]); 
      deltas[3*j+1] = - post[j]/params[2*net->m+j] *
	((y[i]-params[net->m+j])*(y[i]-params[net->m+j])/
	 (params[2*net->m+j]*params[2*net->m+j])-1)
	*(1-exp(MINSIGMA-params[2*net->m+j]));	   /* sigma_j */
    } /* for j loop */
    
    
      /* Compute gradient with respect to nn parameters */
    endptrlin=nllgrad+net->nout*(net->d+1);
    j=0;
    k=0; /* output */
    for(ptrlin=nllgrad;ptrlin<endptrlin;ptrlin++){
      if (j==0) /* bias */
	*ptrlin+=deltas[k];
      else
	*ptrlin+=deltas[k]*x[i*net->d+j-1];
      j++;
      if (j == net->d+1){
	j=0;
	k++;
      }
    }
      
      
    if (net->h > 0) {
      double **ptr, **endptr=net->theta+net->h;
      k=0;  /* hidden unit */
      int l = net->nout*(net->d+1); /* # of linear parameters */
      /* Derivative w/r to output layer weights */
      for(ptr=net->theta;ptr<endptr;ptr++){
	endptrlin=*ptr+net->d+net->nout+1;
	j=0; /* output */
	dldz[k]=0;
	for(ptrlin=*ptr+net->d+1;ptrlin<endptrlin;ptrlin++){
	  nllgrad[l+k*(net->nout+1+net->d)+j+net->d+1]+=
	    deltas[j] * z[k];
	  dldz[k]+=deltas[j]*(*ptrlin);
	  j++;
	}
	k++;
      }
      
      /* Derivative w/r to hidden layer weights */
      for(k=0;k<net->h;k++){
	nllgrad[l+k*(net->d+net->nout+1)]+=(1-z[k]*z[k])*dldz[k];
	for(j=0;j<net->d;j++)
	  nllgrad[l+k*(net->d+net->nout+1)+1+j]+=(1-z[k]*z[k])*dldz[k]*
	    x[i*net->d+j];
      }
    }
  }  /* for i loop */  
}


/* Create a CMM with Gaussian components from R inputs 
   and call cmmgnll() */
void cmmgnllR(double *params,int *d,int *h,int *m, double *x, 
	      double *y, int *n, double *nll, double *nllgrad){
  CMM net;
  double **ptr,**endptr;
  int l;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*3-1;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }
  cmmgnll(&net,x,y,*n,nll,nllgrad);
}



/* Create a CMM with Gaussian components from R inputs 
   and perform quantile computations */
void cmmgquant(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq){
  CMM net;
  double **ptr,**endptr;
  int l,i,j;
  double tol=10^(-16); /* Precision on estimate */
  int itmax=0; /* No Newton iterations */

  net.psi=theta;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*3-1;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }

  double *params = (double *) R_alloc(3*net.m,sizeof(double));
  double *act = (double *) R_alloc(net.nout,sizeof(double));
  double *z = (double *) R_alloc(net.h,sizeof(double));


  for (i=0;i<*n;i++){
    cmmgfwd(&net,x+i*net.d,params,act,z);
    for(j=0;j<*nq;j++)
      ummquant(ummgcdf,ummgpdf,params,net.m,q[j],a,b,tol,itmax,xq+j+i**nq); 
  }

}

/* Create a CMM with Gaussian components truncated below zero
   from R inputs and perform quantile computations */
void cmmgquant_trunc(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq){
  CMM net;
  double **ptr,**endptr;
  int l,i,j;
  double tol=10^(-16); /* Precision on estimate */
  int itmax=0; /* No Newton iterations */
  double F0; /* Probability of being smaller than 0 with CMMG */

  net.psi=theta;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*3-1;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }

  double *params = (double *) R_alloc(3*net.m,sizeof(double));
  double *act = (double *) R_alloc(net.nout,sizeof(double));
  double *z = (double *) R_alloc(net.h,sizeof(double));


  for (i=0;i<*n;i++){
    cmmgfwd(&net,x+i*net.d,params,act,z);
    F0=ummgcdf(params,net.m,0.0);
    for(j=0;j<*nq;j++)
      ummquant(ummgcdf,ummgpdf,params,net.m,q[j]*(1.0-F0),a,b,tol,itmax,xq+j+i**nq); 
  }

}


/* Create a CMM with Gaussian components from R inputs 
   and call cmmgfwd() */
void cmmgfwdR(double *params,int *d,int *h,int *m, double *x, 
	      int *n, double *params_mixt, double *a, double *z){
  CMM net;
  double **ptr,**endptr;
  int l,i;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*3-1;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }
  for(i=0;i<*n;i++)
    cmmgfwd(&net,x+i*net.d,params_mixt+i*3*net.m,a+i*net.nout,z+i*net.h);
}


/* Forward propagation through a conditional mixture of Gaussians:
   prior of the first component determined from the remaining priors. 
   The priors are constrained to lie in the interval [MINPI, 1-MINPI]
   to ease computations. The first input is the probability of the Dirac */
void cmmgfwd_dirac(CMM *net, double *x, double *params, double *a, double *z){
  int i;
  double fprior=1.0;

  /* Propagate forward  through nn */
  nnlin(net, x, a, z);  

  /* First output is the Dirac probability (nout = 3m) */
  params[0] = 1.0/(1.0+exp(-a[0]));

  /* Apply transfer functions to nn outputs. */
  for(i=net->m;i>0;i--){
    if(i>1){
      /* priors in [0,1] */
      params[i]=(1.0/(1.0+exp(-a[3*(i-1)]))*(1.0-2.0*MINPI)+MINPI)*fprior;  
      fprior=fprior-params[i];
    }
    else
      params[i]=fprior;
    params[net->m+i]=a[3*(i-1)+1];   /* means */
    params[2*net->m+i]=softplus(a[3*(i-1)+2])+MINSIGMA; /*standard deviations*/
  }

}


/* Create a CMM with Gaussian components from R inputs 
   and call cmmgfwd_dirac() */
void cmmgfwd_diracR(double *params,int *d,int *h,int *m, double *x, 
	      int *n, double *params_mixt, double *a, double *z){
  CMM net;
  double **ptr,**endptr;
  int l,i;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*3;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }
  for(i=0;i<*n;i++)
    cmmgfwd_dirac(&net,x+i*net.d,params_mixt+i*(3*net.m+1),a+i*net.nout,z+i*net.h);
}



/* Gaussian conditional mixture negative log-likelihood at n points*/
/* First component is a delta Dirac for the probability of null observation*/
void cmmgnll_dirac(CMM *net, double *x, double *y, int n, 
	       double *nll, double *nllgrad){
  double *params, *a, *z, *act, *post, *deltas, *dldz;
  int i,j,k,nparams=net->nout*(net->d+1)+net->h*(1+net->d+net->nout);
  double lprob,ppost,pprior,lprior;
  double *ptrlin,*endptrlin;

  /* Allocate memory for pointers. */
  params = (double*) R_alloc(net->nout+1, sizeof(double));
  a = (double*) R_alloc(net->nout, sizeof(double));
  z = (double*) R_alloc(net->h, sizeof(double));
  act = (double*) R_alloc(net->m, sizeof(double));
  post = (double*) R_alloc(net->m, sizeof(double));
  deltas = (double*) R_alloc(net->nout, sizeof(double));
  dldz = (double*) R_alloc(net->h, sizeof(double));
  

  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;
  for(j=0;j<nparams;j++)
    nllgrad[j]=0;

  for (i=0; i<n; i++){
    /* Get nn outputs, hidden unit activations and mixture parameters for
     pattern i */
    cmmgfwd_dirac(net,x+i*net->d,params,a,z);

    if (y[i] > 0){ /* positive obs. (it rained that day) */
      /* Compute model probability for pattern i */
    /* mth component */
      act[net->m-1] = normlogpdf(params[2*net->m], params[3*net->m], y[i]);
      if(a[3*net->m-3]>0){
	pprior=-log(1+exp(-a[3*net->m-3])); /* re-use computation */
	post[net->m-1]=pprior+act[net->m-1];  /* log(prior) + log(act) */
	lprob=post[net->m-1]; /* accumulates log(sum priors act) */
	lprior=-a[3*net->m-3]+pprior; /* accumulates sum log(1-g(a_{4*j-1}) */
      }
      else{
	pprior=-log(1+exp(a[3*net->m-3]));
	post[net->m-1]=a[3*net->m-3]+pprior+act[net->m-1];
	lprob=post[net->m-1];
	lprior=pprior;
      }
      for(j=net->m-2; j>-1; j--){
	act[j] = normlogpdf(params[net->m+j+1], params[2*net->m+j+1],y[i]);
	if(j>0){
	  if(a[3*j]>0){
	    pprior=-log(1+exp(-a[3*j]));
	    post[j]=pprior+lprior+act[j];
	    lprior+=-a[3*j]+pprior;
	  }
	  else{
	    pprior=-log(1+exp(a[3*j]));
	    post[j]=a[3*j]+pprior+lprior+act[j];
	    lprior+=pprior;
	  }
	}
	else 
	  post[j]=lprior+act[j];
      
      if (lprob > post[j])
	lprob=lprob + log(1+exp(post[j]-lprob));
      else
	lprob=post[j] + log(1+exp(lprob-post[j]));            
    } /* for j */

      
      *nll -= (lprob+ log(params[0]));

      /* Compute gradient for pattern i */
      ppost=0.0; /* Maintains \sum_{i=0}^{j-1} post[i] */
      pprior=0.0; /* Maintains \sum_{i=0}^{j-1} params[i] */
      for (j=0; j<net->m; j++){
	/* Compute deltas for priors */
	post[j]=exp(post[j]-lprob); /* posterior */

	if(j > 0){ /* there are m-1 priors */
	  ppost+=post[j-1];	
	  deltas[3*j] = (params[j+1]*ppost-pprior*post[j])/
	    (pprior+params[j+1])*(1.0-2.0*MINPI);
	}
	pprior+=params[j+1];


	/* Compute deltas for the Gaussian parameters */
	deltas[3*j+1] = - post[j] * /* mu_j */
	  (y[i]-params[net->m+j+1])/(params[2*net->m+j+1]*params[2*net->m+j+1]); 
	deltas[3*j+2] = - post[j]/params[2*net->m+j+1] *
	  ((y[i]-params[net->m+j+1])*(y[i]-params[net->m+j+1])/
	   (params[2*net->m+j+1]*params[2*net->m+j+1])-1)
	  *(1-exp(MINSIGMA-params[2*net->m+j+1]));	   /* sigma_j */
      } /* for j loop */
      
      /* Delta for the Dirac parameter */
      deltas[0] = params[0]-1.0;
    }/* positive observation condition */
    else { /* null observation (Dirac probability) */
      *nll -= log(1.0-params[0]); 
      deltas[0] = params[0];
      for (j=1; j<net->nout; j++)
	deltas[j] = 0.0;
    }
    
    /* Compute gradient with respect to nn parameters */
    endptrlin=nllgrad+net->nout*(net->d+1);
    j=0;
    k=0; /* output */
    for(ptrlin=nllgrad;ptrlin<endptrlin;ptrlin++){
      if (j==0) /* bias */
	*ptrlin+=deltas[k];
      else
	*ptrlin+=deltas[k]*x[i*net->d+j-1];
      j++;
      if (j == net->d+1){
	j=0;
	k++;
      }
    }
        
    if (net->h > 0) {
      double **ptr, **endptr=net->theta+net->h;
      k=0;  /* hidden unit */
      int l = net->nout*(net->d+1); /* # of linear parameters */
      /* Derivative w/r to output layer weights */
      for(ptr=net->theta;ptr<endptr;ptr++){
	endptrlin=*ptr+net->d+net->nout+1;
	j=0; /* output */
	dldz[k]=0;
	for(ptrlin=*ptr+net->d+1;ptrlin<endptrlin;ptrlin++){
	  nllgrad[l+k*(net->nout+1+net->d)+j+net->d+1]+=
	    deltas[j] * z[k];
	  dldz[k]+=deltas[j]*(*ptrlin);
	  j++;
	}
	k++;
      }
      
      /* Derivative w/r to hidden layer weights */
      for(k=0;k<net->h;k++){
	nllgrad[l+k*(net->d+net->nout+1)]+=(1-z[k]*z[k])*dldz[k];
	for(j=0;j<net->d;j++)
	  nllgrad[l+k*(net->d+net->nout+1)+1+j]+=(1-z[k]*z[k])*dldz[k]*
	    x[i*net->d+j];
      }
    }
  } /* for i loop */ 
}  




/* Create a CMM with Gaussian components from R inputs 
   and call cmmgnll_dirac() */
void cmmgnll_diracR(double *params,int *d,int *h,int *m, double *x, 
		    double *y, int *n, double *nll,double *nllgrad){
  CMM net;
  double **ptr,**endptr;
  int l;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*3;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }
  cmmgnll_dirac(&net,x,y,*n,nll,nllgrad);
}


/* Create a CMM with Gaussian components truncated below zero 
   and one Dirac at zero from R inputs and perform quantile computations */
void cmmgquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq){
  CMM net;
  double **ptr,**endptr;
  int l,i,j;
  double tol=10^(-16); /* Precision on estimate */
  int itmax=0; /* No Newton iterations */
  double F0; /* Probability of being smaller than 0 with CMMG */

  net.psi=theta;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*3;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }

  double *params = (double *) R_alloc(net.nout+1,sizeof(double));
  double *params_cmmg = (double *) R_alloc(net.nout,sizeof(double));
  double *out = (double *) R_alloc(net.nout,sizeof(double));
  double *z = (double *) R_alloc(net.h,sizeof(double));


  for (i=0;i<*n;i++){
    cmmgfwd_dirac(&net,x+i*net.d,params,out,z);
    params_cmmg = params+1;
    for(j=0;j<*nq;j++){
      if (q[j] > 1.0-params[0]){
	ummquant(ummgcdf,ummgpdf,params_cmmg,net.m,((q[j]-1.0)/params[0]+1.0)*(1.0-F0)+F0,a,b,tol,itmax,xq+j+i**nq); 
      } else
	xq[j+i**nq]=0.0;
    }
  }

}


/* Create a CMM with Gaussian components truncated below zero
   and one Dirac at zero from R inputs and perform conditional 
   quantile computations */
void cmmgcquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq){
  CMM net;
  double **ptr,**endptr;
  int l,i,j;
  double tol=10^(-16); /* Precision on estimate */
  int itmax=0; /* No Newton iterations */
  double F0; /* Probability of being smaller than 0 with CMMG */

  net.psi=theta;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*3;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }

  double *params = (double *) R_alloc(net.nout+1,sizeof(double));
  double *params_cmmg = (double *) R_alloc(net.nout,sizeof(double));
  double *out = (double *) R_alloc(net.nout,sizeof(double));
  double *z = (double *) R_alloc(net.h,sizeof(double));


  for (i=0;i<*n;i++){
    cmmgfwd_dirac(&net,x+i*net.d,params,out,z);
    params_cmmg = params+1;
    for(j=0;j<*nq;j++){
      ummquant(ummgcdf,ummgpdf,params_cmmg,net.m,q[j]*(1.0-F0)+F0,a,b,tol,itmax,xq+j+i**nq); 
    }    
  }
}


/* Log-Normal conditional mixture negative log-likelihood at n points*/
/* First component is a delta Dirac for the probability of null observation*/
void cmmlnll_dirac(CMM *net, double *x, double *y, int n, 
		   double *nll, double *nllgrad){
  double *params, *a, *z, *act, *post, *deltas, *dldz;
  int i,j,k,nparams=net->nout*(net->d+1)+net->h*(1+net->d+net->nout);
  double lprob,ppost,pprior,lprior;
  double *ptrlin,*endptrlin;

  /* Allocate memory for pointers. */
  params = (double*) R_alloc(net->nout+1, sizeof(double));
  a = (double*) R_alloc(net->nout, sizeof(double));
  z = (double*) R_alloc(net->h, sizeof(double));
  act = (double*) R_alloc(net->m, sizeof(double));
  post = (double*) R_alloc(net->m, sizeof(double));
  deltas = (double*) R_alloc(net->nout, sizeof(double));
  dldz = (double*) R_alloc(net->h, sizeof(double));
  

  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;
  for(j=0;j<nparams;j++)
    nllgrad[j]=0;

  for (i=0; i<n; i++){
    /* Get nn outputs, hidden unit activations and mixture parameters for
     pattern i */
    cmmgfwd_dirac(net,x+i*net->d,params,a,z);

    if (y[i] > 0){ /* positive obs. (it rained that day) */
      /* Compute model probability for pattern i */
    /* mth component */
      act[net->m-1] = normlogpdf(params[2*net->m], params[3*net->m], log(y[i]));
      if(a[3*net->m-3]>0){
	pprior=-log(1+exp(-a[3*net->m-3])); /* re-use computation */
	post[net->m-1]=pprior+act[net->m-1];  /* log(prior) + log(act) */
	lprob=post[net->m-1]; /* accumulates log(sum priors act) */
	lprior=-a[3*net->m-3]+pprior; /* accumulates sum log(1-g(a_{4*j-1}) */
      }
      else{
	pprior=-log(1+exp(a[3*net->m-3]));
	post[net->m-1]=a[3*net->m-3]+pprior+act[net->m-1];
	lprob=post[net->m-1];
	lprior=pprior;
      }
      for(j=net->m-2; j>-1; j--){
	act[j] = normlogpdf(params[net->m+j+1], params[2*net->m+j+1],log(y[i]));
	if(j>0){
	  if(a[3*j]>0){
	    pprior=-log(1+exp(-a[3*j]));
	    post[j]=pprior+lprior+act[j];
	    lprior+=-a[3*j]+pprior;
	  }
	  else{
	    pprior=-log(1+exp(a[3*j]));
	    post[j]=a[3*j]+pprior+lprior+act[j];
	    lprior+=pprior;
	  }
	}
	else 
	  post[j]=lprior+act[j];
      
      if (lprob > post[j])
	lprob=lprob + log(1+exp(post[j]-lprob));
      else
	lprob=post[j] + log(1+exp(lprob-post[j]));            
    } /* for j */

      
      *nll -= (lprob+ log(params[0]));

      /* Compute gradient for pattern i */
      ppost=0.0; /* Maintains \sum_{i=0}^{j-1} post[i] */
      pprior=0.0; /* Maintains \sum_{i=0}^{j-1} params[i] */
      for (j=0; j<net->m; j++){
	/* Compute deltas for priors */
	post[j]=exp(post[j]-lprob); /* posterior */

	if(j > 0){ /* there are m-1 priors */
	  ppost+=post[j-1];	
	  deltas[3*j] = (params[j+1]*ppost-pprior*post[j])/
	    (pprior+params[j+1])*(1.0-2.0*MINPI);
	}
	pprior+=params[j+1];


	/* Compute deltas for the Gaussian parameters */
	deltas[3*j+1] = - post[j] * /* mu_j */
	  (log(y[i])-params[net->m+j+1])/(params[2*net->m+j+1]*params[2*net->m+j+1]); 
	deltas[3*j+2] = - post[j]/params[2*net->m+j+1] *
	  ((log(y[i])-params[net->m+j+1])*(log(y[i])-params[net->m+j+1])/
	   (params[2*net->m+j+1]*params[2*net->m+j+1])-1)
	  *(1-exp(MINSIGMA-params[2*net->m+j+1]));	   /* sigma_j */
      } /* for j loop */
      
      /* Delta for the Dirac parameter */
      deltas[0] = params[0]-1.0;
    }/* positive observation condition */
    else { /* null observation (Dirac probability) */
      *nll -= log(1.0-params[0]); 
      deltas[0] = params[0];
      for (j=1; j<net->nout; j++)
	deltas[j] = 0.0;
    }
    
    /* Compute gradient with respect to nn parameters */
    endptrlin=nllgrad+net->nout*(net->d+1);
    j=0;
    k=0; /* output */
    for(ptrlin=nllgrad;ptrlin<endptrlin;ptrlin++){
      if (j==0) /* bias */
	*ptrlin+=deltas[k];
      else
	*ptrlin+=deltas[k]*x[i*net->d+j-1];
      j++;
      if (j == net->d+1){
	j=0;
	k++;
      }
    }
        
    if (net->h > 0) {
      double **ptr, **endptr=net->theta+net->h;
      k=0;  /* hidden unit */
      int l = net->nout*(net->d+1); /* # of linear parameters */
      /* Derivative w/r to output layer weights */
      for(ptr=net->theta;ptr<endptr;ptr++){
	endptrlin=*ptr+net->d+net->nout+1;
	j=0; /* output */
	dldz[k]=0;
	for(ptrlin=*ptr+net->d+1;ptrlin<endptrlin;ptrlin++){
	  nllgrad[l+k*(net->nout+1+net->d)+j+net->d+1]+=
	    deltas[j] * z[k];
	  dldz[k]+=deltas[j]*(*ptrlin);
	  j++;
	}
	k++;
      }
      
      /* Derivative w/r to hidden layer weights */
      for(k=0;k<net->h;k++){
	nllgrad[l+k*(net->d+net->nout+1)]+=(1-z[k]*z[k])*dldz[k];
	for(j=0;j<net->d;j++)
	  nllgrad[l+k*(net->d+net->nout+1)+1+j]+=(1-z[k]*z[k])*dldz[k]*
	    x[i*net->d+j];
      }
    }
  } /* for i loop */ 
}  




/* Create a CMM with Log-Normal components from R inputs 
   and call cmmlnll_dirac() */
void cmmlnll_diracR(double *params,int *d,int *h,int *m, double *x, 
		    double *y, int *n, double *nll,double *nllgrad){
  CMM net;
  double **ptr,**endptr;
  int l;

  net.psi=params;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*3;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }
  cmmlnll_dirac(&net,x,y,*n,nll,nllgrad);
}


/* Create a CMM with LogNormal components and one Dirac at zero from 
   R inputs and perform quantile computations */
void cmmlquant_dirac(double *theta,int *d,int *h,int *m, double *x, 
	       int *n, double *q, int *nq, double *a, double *b,
	       double *xq){
  CMM net;
  double **ptr,**endptr;
  int l,i,j;
  double tol=10^(-16); /* Precision on estimate */
  int itmax=0; /* No Newton iterations */

  net.psi=theta;
  net.h=*h;
  net.m=*m;
  net.nout=net.m*3;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }

  double *params = (double *) R_alloc(net.nout+1,sizeof(double));
  double *params_cmmg = (double *) R_alloc(net.nout,sizeof(double));
  double *out = (double *) R_alloc(net.nout,sizeof(double));
  double *z = (double *) R_alloc(net.h,sizeof(double));


  for (i=0;i<*n;i++){
    cmmgfwd_dirac(&net,x+i*net.d,params,out,z);
    params_cmmg = params+1;
    for(j=0;j<*nq;j++){
      if (q[j] > 1.0-params[0]){
	ummquant(ummgcdf,ummgpdf,params_cmmg,net.m,(q[j]-1.0)/params[0]+1.0,a,b,tol,itmax,xq+j+i**nq);
	xq[j+i**nq]=exp(xq[j+i**nq]);
      } else
	xq[j+i**nq]=0.0;
    }
  }

}



/* ----------------------------------------------------------------------------------------*/
/* Gamma-Bernouilli conditional distribution: parameters are predicted by a neural network */
/* ----------------------------------------------------------------------------------------*/

/* Forward propagation through a conditional Bernouilli-Gamma mixture.
   The first input is the probability of success of the Bernouilli, then scale (theta)
   and shape (gamma) parameters of the Gamma distribution. */
void cmmbergam_fwd(CMM *net, double *x, double *params, double *a, double *z){
  int i;

  /* Propagate forward  through nn */
  nnlin(net, x, a, z);  

  /* First output is the probability of success of the Bernouilli */
  params[0] = 1.0/(1.0+exp(-a[0]));

  /* Second and third outputs are the Gamma parameters. */
  params[1] = softplus(a[1])+MINSIGMA; /* shape paramater k */
  params[2] = softplus(a[2])+MINSIGMA; /* scale paramater theta */
}


/* Create a CMM for a Bernouilli-Gamma mixture from R inputs 
   and call cmmbergam_fwd() */
void cmmbergam_fwdR(double *params,int *d,int *h, double *x, 
	      int *n, double *params_bergam, double *a, double *z){
  CMM net;
  double **ptr,**endptr;
  int l,i;

  net.psi=params;
  net.h=*h;
  net.m=0;
  net.nout=3;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }
  for(i=0;i<*n;i++)
    cmmbergam_fwd(&net,x+i*net.d,params_bergam+i*net.nout,a+i*net.nout,z+i*net.h);
}


/* Bernouilli-Gamma conditional mixture negative log-likelihood at n points*/
void cmmbergam_nll(CMM *net, double *x, double *y, int n, double *nll, 
		   double *nllgrad){


  double *params, *a, *z, *deltas, *dldz;
  int i,j,k,nparams=net->nout*(net->d+1)+net->h*(1+net->d+net->nout);
  double *ptrlin,*endptrlin;

  /* Allocate memory for pointers. */
  params = (double*) R_alloc(net->nout, sizeof(double));
  a = (double*) R_alloc(net->nout, sizeof(double));
  z = (double*) R_alloc(net->h, sizeof(double));
  deltas = (double*) R_alloc(net->nout, sizeof(double));
  dldz = (double*) R_alloc(net->h, sizeof(double));

  /* Make sure initialization of neg-log-like and its gradient
     is zero */
  *nll=0.0;
  for(j=0;j<nparams;j++)
    nllgrad[j]=0;

  for (i=0; i<n; i++){
    /* Get nn outputs, hidden unit activations and mixture parameters for
     pattern i */
    cmmbergam_fwd(net,x+i*net->d,params,a,z);

    if (y[i] > 0){ /* positive obs. (it rained that day) */      
      *nll -= ((params[1]-1.0)*log(y[i]) - params[1]*log(params[2]) - lgammafn(params[1])
	       - y[i]/params[2] + log(params[0]));
      
      deltas[0] = params[0]-1.0; /* Delta for the Dirac parameter */
      deltas[1] = (digamma(params[1])-log(y[i]/params[2]))*(1-exp(MINSIGMA-params[1]));
      deltas[2] = (params[1] - y[i]/params[2])*(1-exp(MINSIGMA-params[2]))/params[2];

    } /* positive observation condition */
    else { /* null observation  */
      *nll -= log(1.0-params[0]); 
      deltas[0] = params[0];
      deltas[1] = 0.0;
      deltas[2] = 0.0;
    }
    
    /* Rprintf("nll[%d] = %f\n",i+1,*nll); */
       
      /* Compute gradient with respect to nn parameters */
    endptrlin=nllgrad+net->nout*(net->d+1);
    j=0;
    k=0; /* output */
    for(ptrlin=nllgrad;ptrlin<endptrlin;ptrlin++){
      if (j==0) /* bias */
	*ptrlin+=deltas[k];
      else
	*ptrlin+=deltas[k]*x[i*net->d+j-1];
      j++;
      if (j == net->d+1){
	j=0;
	k++;
      }
    }
            
    if (net->h > 0) {
      double **ptr, **endptr=net->theta+net->h;
      k=0;  /* hidden unit */
      int l = net->nout*(net->d+1); /* # of linear parameters */
      /* Derivative w/r to output layer weights */
      for(ptr=net->theta;ptr<endptr;ptr++){
	endptrlin=*ptr+net->d+net->nout+1;
	j=0; /* output */
	dldz[k]=0;
	for(ptrlin=*ptr+net->d+1;ptrlin<endptrlin;ptrlin++){
	  nllgrad[l+k*(net->nout+1+net->d)+j+net->d+1]+=
	    deltas[j] * z[k];
	  dldz[k]+=deltas[j]*(*ptrlin);
	  j++;
	}
	k++;
      }
	
      /* Derivative w/r to hidden layer weights */
      for(k=0;k<net->h;k++){
	nllgrad[l+k*(net->d+net->nout+1)]+=(1-z[k]*z[k])*dldz[k];
	for(j=0;j<net->d;j++)
	  nllgrad[l+k*(net->d+net->nout+1)+1+j]+=(1-z[k]*z[k])*dldz[k]*
	    x[i*net->d+j];
      }
    }
  }  /* for i loop */
  
}


/* Create a CMM for a Bernouilli-Gamma mixture from R inputs 
   and call cmmbergam_nll() */
void cmmbergam_nllR(double *params,int *d,int *h, double *x, 
		    double *y, int *n, double *nll,double *nllgrad){
  CMM net;
  double **ptr,**endptr;
  int l;

  net.psi=params;
  net.h=*h;
  net.m=0;
  net.nout=3;
  net.d=*d;
  net.s=0;

  if(net.h > 0){
    /* if so, we need to split the parameter vector into 
       the linear and non-linear part */
    l = net.nout*(net.d+1); /* start of non-linear part */
    endptr=net.theta+net.h;
    for(ptr=net.theta;ptr<endptr;ptr++){
      *ptr=net.psi+l;
      l+=net.d+net.nout+1;
    }        
  }
  cmmbergam_nll(&net,x,y,*n,nll,nllgrad);
}
