/* Adapted to AMT by Robert Baumgarnter (2015) and Piotr Majdak (2016) */


/* This is Version 5 of the code for auditory periphery model of:

    Zilany, M.S.A., Bruce, I.C., Nelson, P.C., and Carney, L.H. (2009). "A Phenomenological
        model of the synapse between the inner hair cell and auditory nerve : Long-term adaptation 
        with power-law dynamics," Journal of the Acoustical Society of America 126(5): 2390-2412.        

   with the modifications and simulation options described in:

    Zilany, M.S.A., Bruce, I.C., Ibrahim, R.A., and Carney, L.H. "Improved parameters
        and expanded simulation options for a model of the auditory periphery," 
        submitted to Journal of the Acoustical Society of America.

   Humanization in this version includes:
   - Human middle-ear filter, based on the linear middle-ear circuit model of Pascal et al. (JASA 1998)
   - Human BM tuning, based on Shera et al. (PNAS 2002)
   - Human frequency-offset of control-path filter (i.e., cochlear amplifier mechanism), based on Greenwood (JASA 1990)
   - Human latency vs CF function, based on Harte et al. (JASA 2009)

   The modifications to the BM tuning are described in:

        Ibrahim, R. A., and Bruce, I. C. (2010). "Effects of peripheral tuning on the auditory nerve's representation
            of speech envelope and temporal fine structure cues," in The Neurophysiological Bases of Auditory Perception,
            eds. E. A. Lopez-Poveda and A. R. Palmer and R. Meddis, Springer, NY, pp. 429438.

   Please cite these papers if you publish any research
   results obtained with this code or any modified versions of this code.

   See the file readme.txt for details of compiling and running the model.  
   
   %%% © M. S. Arefeen Zilany (msazilany@gmail.com), Ian C. Bruce (ibruce@ieee.org),
         Rasha A. Ibrahim, Paul C. Nelson, and Laurel H. Carney - December 2012 %%%
   
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>      /* Added for MS Visual C++ compatability, by Ian Bruce, 1999 */
#include <mex.h>
#include <time.h>
/* #include <iostream.h> */

/* #include "zilany2013_complex.hpp" */
#include "zilany2014_complex.c"

#define MAXSPIKES 1000000
#ifndef TWOPI
#define TWOPI 6.28318530717959
#endif

#ifndef __max
#define __max(a,b) (((a) > (b))? (a): (b))
#endif

#ifndef __min
#define __min(a,b) (((a) < (b))? (a): (b))
#endif

/* This function is the MEX "wrapper", to pass the input and output variables between the .dll or .mexglx file and Matlab */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	double *px, cf, tdres, fibertype, noiseType, implnt;
	int    nrep, pxbins, lp, outsize[2], totalstim;

	double *pxtmp, *cftmp, *nreptmp, *tdrestmp, *fibertypetmp, *noiseTypetmp, *implnttmp;
        
    double *meanrate, *varrate, *psth;
   
	void   SingleAN(double *, double, int, double, int, double, double, double, double *, double *, double *);
	
	/* Check for proper number of arguments */
	
	if (nrhs != 7) 
	{
		mexErrMsgTxt("model_Synapse requires 7 input arguments.");
	}; 

	if (nlhs != 3)  
	{
		mexErrMsgTxt("model_Synapse requires 3 output argument.");
	};
	
	/* Assign pointers to the inputs */

	pxtmp		= mxGetPr(prhs[0]);
	cftmp		= mxGetPr(prhs[1]);
	nreptmp		= mxGetPr(prhs[2]);
	tdrestmp	= mxGetPr(prhs[3]);
    fibertypetmp= mxGetPr(prhs[4]);
    noiseTypetmp= mxGetPr(prhs[5]);
    implnttmp	= mxGetPr(prhs[6]);
	
	/* Check with individual input arguments */

	pxbins = mxGetN(prhs[0]);
	if (pxbins==1)
		mexErrMsgTxt("px must be a row vector\n");
	
	cf = cftmp[0];
	
	nrep = (int)nreptmp[0];
	if (nreptmp[0]!=nrep)
		mexErrMsgTxt("nrep must an integer.\n");
	if (nrep<1)
		mexErrMsgTxt("nrep must be greater that 0.\n");

    tdres = tdrestmp[0];	
   
	fibertype  = fibertypetmp[0];  /* spontaneous rate of the fiber */
    
    noiseType  = noiseTypetmp[0];  /* fixed or variable fGn */
    
    implnt = implnttmp[0];  /* actual/approximate implementation of the power-law functions */

	/* Calculate number of samples for total repetition time */

	totalstim = (int)floor(pxbins/nrep);    

    px = (double*)mxCalloc(totalstim*nrep,sizeof(double)); 

	/* Put stimulus waveform into pressure waveform */
    
   	for (lp=0; lp<pxbins; lp++)
			px[lp] = pxtmp[lp];
	
	/* Create an array for the return argument */
    
    outsize[0] = 1;
	outsize[1] = totalstim;
	
	plhs[0] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
    
	/* Assign pointers to the outputs */
	
	meanrate	  = mxGetPr(plhs[0]);
    varrate = mxGetPr(plhs[1]);
    psth	  = mxGetPr(plhs[2]);
			
	/* run the model */

	/*mexPrintf("ANmodel: Zilany, Bruce, Ibrahim, and Carney : Auditory Nerve Model\n");*/

	SingleAN(px,cf,nrep,tdres,totalstim,fibertype,noiseType,implnt,meanrate,varrate,psth);

 mxFree(px);

}

void SingleAN(double *px, double cf, int nrep, double tdres, int totalstim, double fibertype, double noiseType, double implnt, double *meanrate, double *varrate, double *psth)
{	
        	
	/*variables for the signal-path, control-path and onward */
	double *synouttmp,*sptime;

	int    i,nspikes,ipst;
	double I,spont;
    double sampFreq = 10e3; /* Sampling frequency used in the synapse */
        
    /* Declarations of the functions used in the program */
	double Synapse(double *, double, double, int, int, double, double, double, double, double *);
	int    SpikeGenerator(double *, double, int, int, double *);
    
    /* Allocate dynamic memory for the temporary variables */
    synouttmp  = (double*)mxCalloc(totalstim*nrep,sizeof(double));
    sptime  = (double*)mxCalloc((long) ceil(totalstim*tdres*nrep/0.00075),sizeof(double));  	
	   
    /* Spontaneous Rate of the fiber corresponding to Fibertype */    
    if (fibertype==1) spont = 0.1;
    if (fibertype==2) spont = 4.0;
    if (fibertype==3) spont = 100.0;
    
    /*====== Run the synapse model ======*/    
    I = Synapse(px, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, synouttmp);
            
    /* Wrapping up the unfolded (due to no. of repetitions) Synapse Output */
    for(i = 0; i<I ; i++)
	{       
		ipst = (int) (fmod(i,totalstim));
        meanrate[ipst] = meanrate[ipst] + synouttmp[i]/nrep;        
	};
    /* Synapse Output taking into account the Refractory Effects (Vannucci and Teich, 1978) */
    for(i = 0; i<totalstim ; i++)
	{       
		varrate[i] = meanrate[i]/pow((1+0.75e-3*meanrate[i]),3); /* estimated instananeous variance in the discharge rate */
        meanrate[i]    = meanrate[i]/(1+0.75e-3*meanrate[i]);  /* estimated instantaneous mean rate */     
	};
    /*======  Spike Generations ======*/
    
	nspikes = SpikeGenerator(synouttmp, tdres, totalstim, nrep, sptime);
	for(i = 0; i < nspikes; i++)
	{        
		ipst = (int) (fmod(sptime[i],tdres*totalstim) / tdres);
        psth[ipst] = psth[ipst] + 1;       
	};

    /* Freeing dynamic memory allocated earlier */

    mxFree(sptime); mxFree(synouttmp); 

} /* End of the SingleAN function */
/* -------------------------------------------------------------------------------------------- */
/*  Synapse model: if the time resolution is not small enough, the concentration of
   the immediate pool could be as low as negative, at this time there is an alert message
   print out and the concentration is set at saturated level  */
/* --------------------------------------------------------------------------------------------*/
double Synapse(double *ihcout, double tdres, double cf, int totalstim, int nrep, double spont, double noiseType, double implnt, double sampFreq, double *synouttmp)
{    
    /* Initalize Variables */     
    int z, b;
    int resamp = (int) ceil(1/(tdres*sampFreq));
    double incr = 0.0; int delaypoint = (int) floor(7500/(cf/1e3));   
    
    double alpha1, beta1, I1, alpha2, beta2, I2, binwidth;
    int    k,j,indx,i;    
    double synstrength,synslope,CI,CL,PG,CG,VL,PL,VI;
	double cf_factor,PImax,kslope,Ass,Asp,TauR,TauST,Ar_Ast,PTS,Aon,AR,AST,Prest,gamma1,gamma2,k1,k2;
	double VI0,VI1,alpha,beta,theta1,theta2,theta3,vsat,tmpst,tmp,PPI,CIlast,temp;
            
    double *sout1, *sout2, *synSampOut, *powerLawIn, *exponOut, *TmpSyn;            
    double *m1, *m2, *m3, *m4, *m5;
	double *n1, *n2, *n3;
    
    mxArray	*randInputArray[5], *randOutputArray[1];
    double *randNums;
    
    mxArray	*IhcInputArray[3], *IhcOutputArray[1];
    double *sampIHC, *ihcDims;	  
        
    exponOut = (double*)mxCalloc((long) ceil(totalstim*nrep),sizeof(double));
    powerLawIn = (double*)mxCalloc((long) ceil(totalstim*nrep+3*delaypoint),sizeof(double));
    sout1 = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    sout2 = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    synSampOut  = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    TmpSyn  = (double*)mxCalloc((long) ceil(totalstim*nrep+2*delaypoint),sizeof(double));
      
    m1 = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    m2 = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    m3  = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double)); 
    m4 = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    m5  = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double)); 
    
    n1 = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    n2 = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));
    n3 = (double*)mxCalloc((long) ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq),sizeof(double));    
	
    /*----------------------------------------------------------*/    
    /*------- Parameters of the Power-law function -------------*/
    /*----------------------------------------------------------*/ 
    binwidth = 1/sampFreq;
    /*alpha1 = 5e-6*100e3; beta1 = 5e-4; I1 = 0;*/ /* older version, 2012 and before */
    alpha1 = 2.5e-6*100e3; beta1 = 5e-4; I1 = 0;
    alpha2 = 1e-2*100e3; beta2 = 1e-1; I2 = 0;       
    /*----------------------------------------------------------*/    
    /*------- Generating a random sequence ---------------------*/
    /*----------------------------------------------------------*/ 
    randInputArray[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(randInputArray[0])= ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq);
    randInputArray[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(randInputArray[1])= 1/sampFreq;
    randInputArray[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(randInputArray[2])= 0.9; /* Hurst index */
    randInputArray[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(randInputArray[3])= noiseType; /* fixed or variable fGn */
    randInputArray[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(randInputArray[4])= spont; /* high, medium, or low */
            
    mexCallMATLAB(1, randOutputArray, 5, randInputArray, "zilany2014_ffGn");
    randNums = mxGetPr(randOutputArray[0]);      
    /*----------------------------------------------------------*/
    /*----- Double Exponential Adaptation ----------------------*/
    /*----------------------------------------------------------*/    
       if (spont==100) cf_factor = __min(800,pow(10,0.29*cf/1e3 + 0.7));
       if (spont==4)   cf_factor = __min(50,2.5e-4*cf*4+0.2);
       if (spont==0.1) cf_factor = __min(1.0,2.5e-4*cf*0.1+0.15);              
	         
	   PImax  = 0.6;                /* PI2 : Maximum of the PI(PI at steady state) */
       kslope = (1+50.0)/(5+50.0)*cf_factor*20.0*PImax;            
       /* Ass    = 300*TWOPI/2*(1+cf/100e3); */  /* Older value: Steady State Firing Rate eq.10 */
       Ass    = 800*(1+cf/100e3);    /* Steady State Firing Rate eq.10 */

       if (implnt==1) Asp = spont*3.0;   /* Spontaneous Firing Rate if actual implementation */
       if (implnt==0) Asp = spont*2.75; /* Spontaneous Firing Rate if approximate implementation */
       TauR   = 2e-3;               /* Rapid Time Constant eq.10 */
       TauST  = 60e-3;              /* Short Time Constant eq.10 */
       Ar_Ast = 6;                  /* Ratio of Ar/Ast */
       PTS    = 3;                  /* Peak to Steady State Ratio, characteristic of PSTH */
   
       /* now get the other parameters */
       Aon    = PTS*Ass;                          /* Onset rate = Ass+Ar+Ast eq.10 */
       AR     = (Aon-Ass)*Ar_Ast/(1+Ar_Ast);      /* Rapid component magnitude: eq.10 */
       AST    = Aon-Ass-AR;                       /* Short time component: eq.10 */
       Prest  = PImax/Aon*Asp;                    /* eq.A15 */
       CG  = (Asp*(Aon-Asp))/(Aon*Prest*(1-Asp/Ass));    /* eq.A16 */
       gamma1 = CG/Asp;                           /* eq.A19 */
       gamma2 = CG/Ass;                           /* eq.A20 */
       k1     = -1/TauR;                          /* eq.8 & eq.10 */
       k2     = -1/TauST;                         /* eq.8 & eq.10 */
               /* eq.A21 & eq.A22 */
       VI0    = (1-PImax/Prest)/(gamma1*(AR*(k1-k2)/CG/PImax+k2/Prest/gamma1-k2/PImax/gamma2));
       VI1    = (1-PImax/Prest)/(gamma1*(AST*(k2-k1)/CG/PImax+k1/Prest/gamma1-k1/PImax/gamma2));
       VI  = (VI0+VI1)/2;
       alpha  = gamma2/k1/k2;       /* eq.A23,eq.A24 or eq.7 */
       beta   = -(k1+k2)*alpha;     /* eq.A23 or eq.7 */
       theta1 = alpha*PImax/VI; 
       theta2 = VI/PImax;
       theta3 = gamma2-1/PImax;
  
       PL  = ((beta-theta2*theta3)/theta1-1)*PImax;  /* eq.4' */
       PG  = 1/(theta3-1/PL);                        /* eq.5' */
       VL  = theta1*PL*PG;                           /* eq.3' */
       CI  = Asp/Prest;                              /* CI at rest, from eq.A3,eq.A12 */
       CL  = CI*(Prest+PL)/PL;                       /* CL at rest, from eq.1 */
   	
       if(kslope>=0)  vsat = kslope+Prest;                
       tmpst  = log(2)*vsat/Prest;
       if(tmpst<400) synstrength = log(exp(tmpst)-1);
       else synstrength = tmpst;
       synslope = Prest/log(2)*synstrength;
       
       k = 0;     
       for (indx=0; indx<totalstim*nrep; ++indx)
       {
            tmp = synstrength*(ihcout[indx]);
            if(tmp<400) tmp = log(1+exp(tmp));
            PPI = synslope/synstrength*tmp;           
         
            CIlast = CI; 
            CI = CI + (tdres/VI)*(-PPI*CI + PL*(CL-CI));
            CL = CL + (tdres/VL)*(-PL*(CL - CIlast) + PG*(CG - CL));
            if(CI<0)
            {
                temp = 1/PG+1/PL+1/PPI;
                CI = CG/(PPI*temp);
                CL = CI*(PPI+PL)/PL;
            };
            exponOut[k] = CI*PPI;
            k=k+1;
        }                 
        for (k=0; k<delaypoint; k++)
			powerLawIn[k] = exponOut[0];    
        for (k=delaypoint; k<totalstim*nrep+delaypoint; k++)
			powerLawIn[k] = exponOut[k-delaypoint];
        for (k=totalstim*nrep+delaypoint; k<totalstim*nrep+3*delaypoint; k++)
			powerLawIn[k] = powerLawIn[k-1];         
   /*----------------------------------------------------------*/ 
   /*------ Downsampling to sampFreq (Low) sampling rate ------*/   
   /*----------------------------------------------------------*/    
    IhcInputArray[0] = mxCreateDoubleMatrix(1, k, mxREAL);
    ihcDims = mxGetPr(IhcInputArray[0]);
    for (i=0;i<k;++i)
        ihcDims[i] = powerLawIn[i];    
    IhcInputArray[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(IhcInputArray[1])= 1;    
    IhcInputArray[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(IhcInputArray[2])= resamp;    
    mexCallMATLAB(1, IhcOutputArray, 3, IhcInputArray, "resample");
    sampIHC = mxGetPr(IhcOutputArray[0]);
    
    mxFree(powerLawIn); mxFree(exponOut);
   /*----------------------------------------------------------*/
   /*----- Running Power-law Adaptation -----------------------*/     
   /*----------------------------------------------------------*/
    k = 0;
    for (indx=0; indx<floor((totalstim*nrep+2*delaypoint)*tdres*sampFreq); indx++)
    {
          sout1[k]  = __max( 0, sampIHC[indx] + randNums[indx]- alpha1*I1); 
          /*sout1[k]  = __max( 0, sampIHC[indx] - alpha1*I1); */   /* No fGn condition */
          sout2[k]  = __max( 0, sampIHC[indx] - alpha2*I2); 
                                   
         if (implnt==1)    /* ACTUAL Implementation */
         {
              I1 = 0; I2 = 0; 
              for (j=0; j<k+1; ++j)
                  {
                      I1 += (sout1[j])*binwidth/((k-j)*binwidth + beta1);
                      I2 += (sout2[j])*binwidth/((k-j)*binwidth + beta2);              
                   }
         } /* end of actual */
              
         if (implnt==0)    /* APPROXIMATE Implementation */
         {              
                if (k==0)
                {
                    n1[k] = 1.0e-3*sout2[k];
                    n2[k] = n1[k]; n3[0]= n2[k];
                }
                else if (k==1)
                {
                    n1[k] = 1.992127932802320*n1[k-1]+ 1.0e-3*(sout2[k] - 0.994466986569624*sout2[k-1]);
                    n2[k] = 1.999195329360981*n2[k-1]+ n1[k] - 1.997855276593802*n1[k-1];
                    n3[k] = -0.798261718183851*n3[k-1]+ n2[k] + 0.798261718184977*n2[k-1];
                }
                else
                {			
                    n1[k] = 1.992127932802320*n1[k-1] - 0.992140616993846*n1[k-2]+ 1.0e-3*(sout2[k] - 0.994466986569624*sout2[k-1] + 0.000000000002347*sout2[k-2]);
                    n2[k] = 1.999195329360981*n2[k-1] - 0.999195402928777*n2[k-2]+n1[k] - 1.997855276593802*n1[k-1] + 0.997855827934345*n1[k-2];
                    n3[k] =-0.798261718183851*n3[k-1] - 0.199131619873480*n3[k-2]+n2[k] + 0.798261718184977*n2[k-1] + 0.199131619874064*n2[k-2];
                }   
                I2 = n3[k];       

                if (k==0)
                {
                    m1[k] = 0.2*sout1[k];
                    m2[k] = m1[k];	m3[k] = m2[k];			
                    m4[k] = m3[k];	m5[k] = m4[k];
                }
                else if (k==1)
                {
                    m1[k] = 0.491115852967412*m1[k-1] + 0.2*(sout1[k] - 0.173492003319319*sout1[k-1]);
                    m2[k] = 1.084520302502860*m2[k-1] + m1[k] - 0.803462163297112*m1[k-1];
                    m3[k] = 1.588427084535629*m3[k-1] + m2[k] - 1.416084732997016*m2[k-1];
                    m4[k] = 1.886287488516458*m4[k-1] + m3[k] - 1.830362725074550*m3[k-1];
                    m5[k] = 1.989549282714008*m5[k-1] + m4[k] - 1.983165053215032*m4[k-1];
                }        
                else
                {
                    m1[k] = 0.491115852967412*m1[k-1] - 0.055050209956838*m1[k-2]+ 0.2*(sout1[k]- 0.173492003319319*sout1[k-1]+ 0.000000172983796*sout1[k-2]);
                    m2[k] = 1.084520302502860*m2[k-1] - 0.288760329320566*m2[k-2] + m1[k] - 0.803462163297112*m1[k-1] + 0.154962026341513*m1[k-2];
                    m3[k] = 1.588427084535629*m3[k-1] - 0.628138993662508*m3[k-2] + m2[k] - 1.416084732997016*m2[k-1] + 0.496615555008723*m2[k-2];
                    m4[k] = 1.886287488516458*m4[k-1] - 0.888972875389923*m4[k-2] + m3[k] - 1.830362725074550*m3[k-1] + 0.836399964176882*m3[k-2];
                    m5[k] = 1.989549282714008*m5[k-1] - 0.989558985673023*m5[k-2] + m4[k] - 1.983165053215032*m4[k-1] + 0.983193027347456*m4[k-2];
                }   
                I1 = m5[k]; 
            } /* end of approximate implementation */
        
        synSampOut[k] = sout1[k] + sout2[k]; 
        k = k+1;                  
      }   /* end of all samples */
      mxFree(sout1); mxFree(sout2);  
      mxFree(m1); mxFree(m2); mxFree(m3); mxFree(m4); mxFree(m5); mxFree(n1); mxFree(n2); mxFree(n3); 
    /*----------------------------------------------------------*/    
    /*----- Upsampling to original (High 100 kHz) sampling rate --------*/  
    /*----------------------------------------------------------*/    
    for(z=0; z<k-1; ++z)
    {    
        incr = (synSampOut[z+1]-synSampOut[z])/resamp;
        for(b=0; b<resamp; ++b)
        {
            TmpSyn[z*resamp+b] = synSampOut[z]+ b*incr; 
        }        
    }      
    for (i=0;i<totalstim*nrep;++i)
        synouttmp[i] = TmpSyn[i+delaypoint];      
    
    mxFree(synSampOut); mxFree(TmpSyn);   
    mxDestroyArray(randInputArray[0]); mxDestroyArray(randOutputArray[0]);
    mxDestroyArray(IhcInputArray[0]); mxDestroyArray(IhcOutputArray[0]); mxDestroyArray(IhcInputArray[1]); mxDestroyArray(IhcInputArray[2]);
    mxDestroyArray(randInputArray[1]);mxDestroyArray(randInputArray[2]); mxDestroyArray(randInputArray[3]);
    mxDestroyArray(randInputArray[4]);
    return((long) ceil(totalstim*nrep));
}    
/* ------------------------------------------------------------------------------------ */
/* Pass the output of Synapse model through the Spike Generator */

/* The spike generator now uses a method coded up by B. Scott Jackson (bsj22@cornell.edu) 
   Scott's original code is available from Laurel Carney's web site at:
   http://www.urmc.rochester.edu/smd/Nanat/faculty-research/lab-pages/LaurelCarney/auditory-models.cfm
*/

int SpikeGenerator(double *synouttmp, double tdres, int totalstim, int nrep, double *sptime) 
{  
   	double  c0,s0,c1,s1,dead;
    int     nspikes,k,NoutMax,Nout,deadtimeIndex,randBufIndex;      
    double	deadtimeRnd, endOfLastDeadtime, refracMult0, refracMult1, refracValue0, refracValue1;
    double	Xsum, unitRateIntrvl, countTime, DT;    
    
    mxArray	*randInputArray[1], *randOutputArray[1];
    double *randNums, *randDims;    
    
    c0      = 0.5;
	s0      = 0.001;
	c1      = 0.5;
	s1      = 0.0125;
    dead    = 0.00075;
    
    DT = totalstim * tdres * nrep;  /* Total duration of the rate function */
    Nout = 0;
    NoutMax = (long) ceil(totalstim*nrep*tdres/dead);    
       
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
    endOfLastDeadtime = __max(0,log(randNums[randBufIndex++]) / synouttmp[0] + dead);  /* End of last deadtime before t=0 */
    refracValue0 = c0*exp(endOfLastDeadtime/s0);     /* Value of first exponential in refractory function */
	refracValue1 = c1*exp(endOfLastDeadtime/s1);     /* Value of second exponential in refractory function */
	Xsum = synouttmp[0] * (-endOfLastDeadtime + c0*s0*(exp(endOfLastDeadtime/s0)-1) + c1*s1*(exp(endOfLastDeadtime/s1)-1));  
        /* Value of time-warping sum */
		/*  ^^^^ This is the "integral" of the refractory function ^^^^ (normalized by 'tdres') */

	/* Calculate first interspike interval in a homogeneous, unit-rate Poisson process (normalized by 'tdres') */
    unitRateIntrvl = -log(randNums[randBufIndex++])/tdres;  
	    /* NOTE: Both 'unitRateInterval' and 'Xsum' are divided (or normalized) by 'tdres' in order to reduce calculation time.  
		This way we only need to divide by 'tdres' once per spike (when calculating 'unitRateInterval'), instead of 
		multiplying by 'tdres' once per time bin (when calculating the new value of 'Xsum').                         */

	countTime = tdres;
	for (k=0; (k<totalstim*nrep) && (countTime<DT); ++k, countTime+=tdres, refracValue0*=refracMult0, refracValue1*=refracMult1)  /* Loop through rate vector */
	{
		if (synouttmp[k]>0)  /* Nothing to do for non-positive rates, i.e. Xsum += 0 for non-positive rates. */
		{
		  Xsum += synouttmp[k]*(1 - refracValue0 - refracValue1);  /* Add synout*(refractory value) to time-warping sum */
			
			if ( Xsum >= unitRateIntrvl )  /* Spike occurs when time-warping sum exceeds interspike "time" in unit-rate process */
			{
				sptime[Nout] = countTime; Nout = Nout+1;								
				unitRateIntrvl = -log(randNums[randBufIndex++]) /tdres; 
                 Xsum = 0;
				
			    /* Increase index and time to the last time bin in the deadtime, and reset (relative) refractory function */
				k += deadtimeIndex;
				countTime += deadtimeRnd;
				refracValue0 = c0;
				refracValue1 = c1;
			}
		}
	} /* End of rate vector loop */			
            
    mxDestroyArray(randInputArray[0]); mxDestroyArray(randOutputArray[0]);	
	nspikes = Nout;  /* Number of spikes that occurred. */
	return(nspikes);
}
