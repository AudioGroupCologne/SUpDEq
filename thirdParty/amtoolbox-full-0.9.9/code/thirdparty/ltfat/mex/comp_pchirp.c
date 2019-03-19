
#include "mex.h"

#include "ltfat.h"

#define PI 3.1415926535897932384626433832795

static inline long long positiverem_long(long long a,long long b)
{
  const long long c = a%b;
  return(c<0 ? c+b : c);
}

/* Calling convention:
 *  pchirp(L,n);
 */
void mexFunction( int UNUSED(nlhs), mxArray *plhs[], 
		            int UNUSED(nrhs), const mxArray *prhs[] )
{ 
   const long long L=(long long) mxGetScalar(prhs[0]);
   const long long n=(long long) mxGetScalar(prhs[1]);

   plhs[0] = mxCreateDoubleMatrix(L, 1, mxCOMPLEX);
   double *gr = mxGetPr(plhs[0]);
   double *gi = mxGetPi(plhs[0]);


   const long long LL=2*L;
   const long long Lponen=positiverem_long((L+1)*n,LL);
   
   for (long long m=0;m<L;m++)
   {
      const long long idx = positiverem_long(
   	 positiverem_long(Lponen*m,LL)*m,LL);

      gr[m] = cos(PI*idx/L);
      gi[m] = sin(PI*idx/L);
   }

   return;
  
}


