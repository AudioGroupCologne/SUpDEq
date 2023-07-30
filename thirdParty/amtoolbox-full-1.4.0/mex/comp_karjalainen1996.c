#include <stdio.h>
#include <math.h>
#include "mex.h"

double  *NFArray(int size)
{
	double	*p;
	p = (double *)mxCalloc(sizeof(*p),size);
	return p;
}

/*******************************************************/
/* Karja_network96                                     */
/* A multi-channel implementation of                   */
/* a nonlinear adaptation network (Karjalainen, 1996)  */
/*       By Aki Härmä 15. 12 1998  		       */
/*******************************************************/
void trans(double *sig, double fs, long int len, long int dim,
	   double *SLOW, double *FAST)
{
 long int q, w, ind;
 double kup1,kdown1,kup2,kdown2,thres,peak,wy1,wy2,wy,wx,k1,k2,y;
 double *y1,*y2;

y1=NFArray(dim); y2=NFArray(dim);
kup1  =   1/fs;
kdown1=   100/fs;
kup2  =   10/fs;
kdown2=   7.5/fs;
thres=0.15;

 for (q=0;q<len;q++)
   {
   for (w=0;w<dim;w++)
     {
     ind=q+w*len;
     peak=log(sig[ind])-3.5; /* Feedforward path */
     wy1=(y1[w]>peak)?y1[w]:peak;
     wy2=(y2[w]>peak-1)?y2[w]:peak-1;
     wy=(wy1+wy2)/2;
     wx=sig[ind]/exp(wy);    /* Input scaled by gain */
     wx=(wx<300)?wx:300;
     FAST[ind]=wx;           /* save FAST output */

     k1=(wx>wy1)?kup1:kdown1;
     y1[w]=k1*wx+(1-k1)*wy1;
     k2=(wx>wy2)?kup2:kdown2;
     y2[w]=k2*wx+(1-k2)*wy2;
     y=(11.5*0.5)*(y1[w]+y2[w]);
     SLOW[ind]=0.14*pow(10,y/40)-thres; /* save sone level */
     }
   }
return;
}

/************************  MEXFUNCTION *****************************/
void mexFunction(
	int		nlhs,
	mxArray	*plhs[],
	int		nrhs,
	const mxArray    *prhs[]
	)
{
  long int len, dim;
  double *sig, *SLOW, *FAST;
  double fs;
 
        sig = mxGetPr(prhs[0]);
      	len = (long int)mxGetM(prhs[0]);
	dim = (long int)mxGetN(prhs[0]);
	/*mexPrintf("Dim: %d   Len: %d",dim,len);*/
	fs = mxGetScalar(prhs[1]);

       plhs[0] = mxCreateDoubleMatrix(len, dim, mxREAL);
	SLOW = mxGetPr(plhs[0]);
       plhs[1] = mxCreateDoubleMatrix(len, dim, mxREAL);
	FAST = mxGetPr(plhs[1]);
	if(len<dim) mexPrintf("More channels than samples?");
	else trans(sig,fs,len,dim,SLOW,FAST);	
	return;
}


