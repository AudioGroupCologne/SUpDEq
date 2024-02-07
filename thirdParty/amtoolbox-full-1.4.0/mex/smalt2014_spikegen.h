#ifndef _SPIKEGEN_H
#define _SPIKEGEN_H
#include <math.h>
#include <stdlib.h>
/* #include <iostream.h> */
#include <time.h>

/* #define MAXSPIKES 100000 */

/* ------------------------------------------------------------------------------------ */
/* Pass the output of Synapse model through the Spike Generator */

/* The spike generator now uses a method coded up by B. Scott Jackson (bsj22@cornell.edu) 
   Scott's original code is available from Laurel Carney's web site at:
   http://www.urmc.rochester.edu/smd/Nanat/faculty-research/lab-pages/LaurelCarney/auditory-models.cfm
*/

int SpikeGenerator(double *synout, double tdres, int totalstim, double nrep, double *sptime) 
{  

	double  c0,s0,c1,s1,dead;
    int     nspikes,k,NoutMax,Nout,deadtimeIndex,randBufIndex;
    double  j, DT;
    double	deadtimeRnd, endOfLastDeadtime, refracMult0, refracMult1, refracValue0, refracValue1;
    double	Xsum, unitRateIntrvl, countTime;

    mxArray	*randInputArray[1], *randOutputArray[1];
    double *randNums, *randDims;    
    
    c0      = 0.5;
	s0      = 0.001;
	c1      = 0.5;
	s1      = 0.0125;
    dead    = 0.00075;
    
    DT = totalstim * tdres;  /* Total duration of the rate function */
    Nout = 0;
    NoutMax = (long) ceil(totalstim*tdres*nrep/dead);
    
    randInputArray[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
    randDims = mxGetPr(randInputArray[0]);
    randDims[0] = 1;
    randDims[1] = NoutMax+1;
    mexCallMATLAB(1, randOutputArray, 1, randInputArray, "rand");
    randNums = mxGetPr(randOutputArray[0]);
    randBufIndex = 0;
    
	/* Calculate useful constants */
	deadtimeIndex = (long) floor(dead/tdres);  /* Integer number of discrete time bins within deadtime */
	deadtimeRnd = deadtimeIndex*tdres;		   /* Deadtime rounded down to length of an integer number of discrete time bins */

	refracMult0 = 1 - tdres/s0;  /* If y0(t) = c0*exp(-t/s0), then y0(t+tdres) = y0(t)*refracMult0 */
	refracMult1 = 1 - tdres/s1;  /* If y1(t) = c1*exp(-t/s1), then y1(t+tdres) = y1(t)*refracMult1 */

	/* Calculate effects of a random spike before t=0 on refractoriness and the time-warping sum at t=0 */
	endOfLastDeadtime = log(randNums[randBufIndex++]) / synout[0] + dead;  /* End of last deadtime before t=0 */
	refracValue0 = c0*exp(endOfLastDeadtime/s0);                     /* Value of first exponential in refractory function */
	refracValue1 = c1*exp(endOfLastDeadtime/s1);                     /* Value of second exponential in refractory function */
	Xsum = synout[0] * (-endOfLastDeadtime + c0*s0*(exp(endOfLastDeadtime/s0)-1) + c1*s1*(exp(endOfLastDeadtime/s1)-1));  
        /* Value of time-warping sum */
		/*  ^^^^ This is the "integral" of the refractory function ^^^^ (normalized by 'tdres') */

	/* Calculate first interspike interval in a homogeneous, unit-rate Poisson process (normalized by 'tdres') */
	unitRateIntrvl = -log(randNums[randBufIndex++]) /tdres;
	    /* NOTE: Both 'unitRateInterval' and 'Xsum' are divided (or normalized) by 'tdres' in order to reduce calculation time.  
		This way we only need to divide by 'tdres' once per spike (when calculating 'unitRateInterval'), instead of 
		multiplying by 'tdres' once per time bin (when calculating the new value of 'Xsum').                         */

	countTime = tdres;
	k = 0;
	for (j=0; j<nrep; ++j)  /* Loop through "stimulus" repetitions */
		{
		for (; (k<totalstim) && (countTime<DT); ++k, countTime+=tdres, 
									refracValue0*=refracMult0, refracValue1*=refracMult1)  /* Loop through rate vector */
			{
			if (synout[k]>0.0)  /* Nothing to do for non-positive rates, i.e. Xsum += 0 for non-positive rates. */
				{
				Xsum += synout[k]*(1 - refracValue0 - refracValue1);  /* Add synout*(refractory value) to time-warping sum */
			
				if ( Xsum >= unitRateIntrvl )  /* Spike occurs when time-warping sum exceeds interspike "time" in unit-rate process */
					{
					sptime[Nout] = countTime; Nout = Nout+1;
					if (Nout >= NoutMax)  /* If the output buffer for the spike times is too small . . . */
						{
						NoutMax += 1000;
						sptime = (double *)mxRealloc(sptime, NoutMax*sizeof(double));  /* . . . allocate additional memory . . . */
						if (sptime == NULL)  mexErrMsgTxt("Out of Memory");
						}
				
					unitRateIntrvl = -log(randNums[randBufIndex++]) /tdres;  /* Next interspike "time" in unit-rate process */
					Xsum = 0;
				
					/* Increase index and time to the last time bin in the deadtime, and reset (relative) refractory function */
					k += deadtimeIndex;
					countTime += deadtimeRnd;
					refracValue0 = c0;
					refracValue1 = c1;
					}
				}
			} /* End of rate vector loop */
		
			/* Reset index and time to begin a new repetion.  Don't just set to zero, since new repetition may start within the 
				deadtime following the last spike in the previous repetition.                                               */
			countTime -= DT;
			k -= totalstim;
		} /* End of "stimulus" repetition loop */

	/* Delete spike(s) that occur after the last repetition of the rate function ends */
	for (; (Nout>0)&&(sptime[Nout-1]>DT); --Nout)  sptime[Nout-1] = mxGetNaN();
	nspikes = Nout;  /* Number of spikes that occurred. */

    mxDestroyArray(randInputArray[0]); mxDestroyArray(randOutputArray[0]);	
    
	return(nspikes);
}

#endif
