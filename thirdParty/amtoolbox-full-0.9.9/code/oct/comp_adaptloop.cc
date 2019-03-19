/* Adapted to AMT by Piotr Majdak, 2016 */
/* development of the original files: 1999 - 2004 Stephan Ewert */

#include <octave/oct.h>
#include <stdlib.h>
#include "math.h"

typedef struct adaptloopstatevar
{
  int    loops;
  int    nsigs;
  double *state;
  double corr;
  double mult;
  double limit;
  double minlvl;
  double *a1;
  double *factor, *expfac, *offset;
} adaptloopstate;

void adaptloop_init(adaptloopstate *s, const int nsigs, const int loops)
{
   s->loops  = loops;
   s->nsigs  = nsigs;
   s->state  = (double*)malloc(loops*nsigs*sizeof(double));
   s->a1     = (double*)malloc(loops*sizeof(double));
   s->factor = (double*)malloc(loops*sizeof(double));
   s->expfac = (double*)malloc(loops*sizeof(double));
   s->offset = (double*)malloc(loops*sizeof(double));   
}

void adaptloop_free(adaptloopstate *s)
{
   free(s->offset);
   free(s->expfac);
   free(s->factor);
   free(s->a1);
   free(s->state);
}

void adaptloop_set(adaptloopstate *s, const int fs, const double limit,
		   const double minlvl, const double *tau)
{

   double *pstate;
   double maxvalue;
   int jj, w, loops;
         
   pstate = s->state;
   loops = s->loops;

   /* get the b0 and a1 of the RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
      and the steady state */
   for ( jj=0; jj<loops; jj++)
   {
      s->a1[jj]=exp(-1/(tau[jj]*((double)fs)));
      
      /* This is a clever way of avoiding exponents by taking the sqrt
       * of the previous entry */
      if ( jj == 0 )
      {
	 pstate[jj] = sqrt(minlvl);
      } else {
	 pstate[jj] = sqrt(pstate[jj-1]);
      }

      /* Copy the initial state to all the other subband signals */
      for (w=1; w<s->nsigs; w++)
      {
	 pstate[jj+w*loops]=pstate[jj];
      }
   }

   /* Calculate constants for overshoot limitation */
   for ( jj=0; jj<loops; jj++)
   {
      /* Max. possible output value */
      maxvalue = (1.0-pstate[jj]*pstate[jj])*limit-1.0;
      
      /* Factor in formula to speed it up  */
      s->factor[jj] = maxvalue*2; 			
      
      /*Exponential factor in output limiting function*/
      s->expfac[jj] = -2./maxvalue;
      s->offset[jj] = maxvalue - 1.0;
   }
   
   s->corr  = pstate[loops-1];
   s->mult  = 100.0/(1-s->corr);
   s->limit = limit;
   s->minlvl = minlvl;   
}

void adaptloop_run(adaptloopstate *s, const double *insig, const int siglen,
		   double *outsig)
{
   double tmp1;
   double *pstate, minlvl, *b0;
   int ii, jj, w;      

   /* Unpack for easier reference */
   pstate = s->state;
   minlvl = s->minlvl;

   b0 = (double*)malloc(s->loops*sizeof(double));

   for ( jj=0; jj<s->loops; jj++)
   {
      b0[jj]=1-s->a1[jj];
   }

   if (s->limit<=1.0)
   {
      /* Main loop, no overshoot limitation */

      for ( w=0; w<s->nsigs; w++)
      {
	 for ( ii=0; ii<siglen; ii++)
	 {	 
	    tmp1 = insig[ii+w*siglen];
	    
	    if ( tmp1 < minlvl )
	    {
	       tmp1 = minlvl;
	    }
	    
	    /* for each aloop */
	    for ( jj=0; jj<s->loops; jj++)
	    {	       
	       tmp1 /= pstate[jj];
	       pstate[jj] = s->a1[jj]*pstate[jj] + b0[jj]*tmp1;	 
	    }
	    
	    /* Scale to model units */		
	    outsig[ii+w*siglen] = (tmp1 - s->corr) * s->mult;	
	 }	
	 pstate+=s->loops;
      }
   }
   else
   {

      /* Main loop, overshoot limitation */
            
      for ( w=0; w<s->nsigs; w++)
      {
	 for ( ii=0; ii<siglen; ii++)
	 {	 
	    tmp1 = insig[ii+w*siglen];
	    
	    if ( tmp1 < minlvl )
	    {
	       tmp1 = minlvl;
	    }
	    
	    /* compute the adaptation loops */
	    for ( jj=0; jj<s->loops; jj++)
	    {
	       tmp1 /= pstate[jj];
	       
	       if (tmp1 > 1.0)
	       {
	       tmp1 = s->factor[jj]/(1+exp(s->expfac[jj]*(tmp1-1)))-s->offset[jj];
	       }
	       pstate[jj] = s->a1[jj]*pstate[jj] + b0[jj]*tmp1;	 
	    }
	    
	    /* Scale to model units */		
	    outsig[ii+w*siglen] = (tmp1 - s->corr) * s->mult;	
	 }
	 pstate+=s->loops;
      }
   }

   free(b0);

}

DEFUN_DLD (comp_adaptloop, args, ,
  "This function calls the C-library\n\
  c=comp_adaptloop(insig,fs,limit,minlvl);\n\
  Yeah.")
{
   
   const Matrix insig = args(0).matrix_value();
   
   const int siglen = insig.rows();
   const int nsigs  = insig.columns();

   const int fs     = args(1).int_value();
   const double limit  = args(2).double_value();
   const double minlvl = args(3).double_value();
   const Matrix tau    = args(4).matrix_value();

   const int nloops = tau.rows()*tau.columns();
   
   adaptloopstate s;

   Matrix outsig(siglen,nsigs);  
   
   adaptloop_init(&s, nsigs, nloops);
   adaptloop_set(&s, fs, limit, minlvl, (const double*)tau.data());
   adaptloop_run(&s, (double*)insig.data(), siglen, (double*)outsig.data());
   adaptloop_free(&s);

   //adaptloop((double*)insig.data(),fs,siglen,nsigs,limit,minlvl,(double*)outsig.data());
   return octave_value (outsig);
}

