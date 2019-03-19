/* Adapted to AMT by Piotr Majdak, 2016 */
/* Development of original files: 
  mhc.c --  Meddis haircell model for monaural and
        binaural models (peripheral preprocessing).
        This file can be compiled by the matlab
        mex compiler. (c) 2001 Jeroen Breebaart
*/

#include "mex.h"
#include "math.h"      
 
      
void meddishaircell(double insig[], int fs, int siglen, int nsigs, double theResult[])
{

   int ii, jj;

/* Parameters from Meddis' April 1990 JASA paper: */
   const double A=5;
   const double B=300;
   double g=2000;
   double y=5.05;
   double l=2500;
   double x=66.31;
   double r=6580;
   const double M=1;

   const double h=50000;


/* internal variables */

   double kt;
   double spont;
   double q;
   double w;
   double temp;
   double replenish;
   double eject;
   double reuptake;
   double reprocess;
   double loss;
   double srate;

   srate=fs;

   /* Rescale by the sampling rate */
   r=r/srate;
   x=x/srate;
   y=y/srate;
   l=l/srate;
   g=g/srate;

   for (jj=0; jj < nsigs; jj++)
   {      
      kt = g*A/(A+B);
      spont = M*y*kt/(l*kt+y*(l+r));
      q=spont*(l+r)/kt;
      w=spont*r/x;
            
      /* do the MHC thing! */
      for (ii=0; ii < siglen; ii++)
      {
	 temp=(insig[ii+jj*siglen]+A+abs(A+insig[ii+jj*siglen]))/2;
	 kt=g*temp/(temp+B);
	 replenish=(y*(M-q)+abs(y*(M-q)))/2;
	 eject=kt*q;
	 loss=l*spont;
	 reuptake=r*spont;
	 reprocess=x*w;
	 
	 q=q+replenish-eject+reprocess;
	 spont=spont+eject-loss-reuptake;
	 w=w+reuptake-reprocess;
	 theResult[ii+jj*siglen]=spont*h;	 
      }
   }
}
       
 
void mexFunction(int nlhs, mxArray  *plhs[],
		 int nrhs, const mxArray  *prhs[])
{

   int siglen,nsigs;            /* temporary array size holders */
   int fs;

   siglen = mxGetM(prhs[0]);
   nsigs  = mxGetN(prhs[0]);
   fs     = mxGetScalar(prhs[1]);

/*   Step 1: Error Checking Step 1a: is nlhs 1?  If not,
     generate an error message and exit mexample (mexErrMsgTxt
     does this for us!) */
   if (nlhs!=1)
      mexErrMsgTxt("mhc requires one output argument.");
      
   /*   Step 1b: is nrhs 2? */
   if (nrhs!=2)
      mexErrMsgTxt("mhc requires two input arguments, \
    prhs[0] and SAMPLERATE");   
   
   plhs[0] = mxCreateDoubleMatrix(siglen, nsigs, mxREAL);
      
   meddishaircell(mxGetPr(prhs[0]), fs, siglen, nsigs,
		  mxGetPr(plhs[0]));
 
}
