/* This is the BEZ2018 version of the code for auditory periphery model from the Carney, Bruce and Zilany labs.
 * 
 * This release implements the version of the model described in:
 *
 *   Bruce, I.C., Erfani, Y., and Zilany, M.S.A. (2018). "A Phenomenological
 *   model of the synapse between the inner hair cell and auditory nerve: 
 *   Implications of limited neurotransmitter release sites," to appear in
 *   Hearing Research. (Special Issue on "Computational Models in Hearing".)
 *
 * Please cite this paper if you publish any research
 * results obtained with this code or any modified versions of this code.
 *
 * See the file readme.txt for details of compiling and running the model.
 *
 * %%% Ian C. Bruce (ibruce@ieee.org), Yousof Erfani (erfani.yousof@gmail.com),
 *     Muhammad S. A. Zilany (msazilany@gmail.com) - December 2017 %%%
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <time.h>

#include "bruce2018_complex.hpp"

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

/* This function is the MEX "wrapper", to pass the input and output variables between the .mex* file and Matlab or Octave */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double *px, cf, tdres, tabs, trel, noiseType, implnt, spont;
    int    nrep, pxbins, lp, totalstim;
    mwSize  outsize[2];
    
    double *pxtmp, *cftmp, *nreptmp, *tdrestmp, *noiseTypetmp, *implnttmp, *sponttmp, *tabstmp, *treltmp;
    
    double *meanrate, *varrate, *psth, *synout, *trd_vector, *trel_vector;
    
    void SingleAN(double *, double, int, double, int, double, double, double, double, double, double *, double *, double *, double *, double *, double *);
    
    /* Check for proper number of arguments */
    
    if (nrhs != 9)
    {
        mexErrMsgTxt("model_Synapse requires 9 input arguments.");
    };
    
    if (nlhs > 6)
    {
        mexErrMsgTxt("model_Synapse has a maximum of 6 output argument.");
    };
    
    /* Assign pointers to the inputs */
    
    pxtmp		 = mxGetPr(prhs[0]);
    cftmp		 = mxGetPr(prhs[1]);
    nreptmp		 = mxGetPr(prhs[2]);
    tdrestmp	 = mxGetPr(prhs[3]);
    noiseTypetmp = mxGetPr(prhs[4]);
    implnttmp	 = mxGetPr(prhs[5]);
    sponttmp	 = mxGetPr(prhs[6]);
    tabstmp      = mxGetPr(prhs[7]);
    treltmp      = mxGetPr(prhs[8]);

    /* Check individual input arguments */
    
    spont = sponttmp[0];
    if ((spont<1e-4)||(spont>180))
        mexErrMsgTxt("spont must in the range [1e-4,180]\n");   
    
    pxbins = (int)mxGetN(prhs[0]);
    if (pxbins==1)
        mexErrMsgTxt("px must be a row vector\n");
    
    cf = cftmp[0];
    
    nrep = (int)nreptmp[0];
    if (nreptmp[0]!=nrep)
        mexErrMsgTxt("nrep must an integer.\n");
    if (nrep<1)
        mexErrMsgTxt("nrep must be greater that 0.\n");
    
    tdres = tdrestmp[0];
    
    noiseType  = noiseTypetmp[0];  /* fixed or variable fGn */
    
    implnt = implnttmp[0];  /* actual/approximate implementation of the power-law functions */
    
    tabs = tabstmp[0];  /* absolute refractory period */
    if ((tabs<0)||(tabs>20e-3))
        mexErrMsgTxt("tabs must in the range [0,20e-3]\n");
    
    trel = treltmp[0];  /* baseline relative refractory period */
    if ((trel<0)||(trel>20e-3))
        mexErrMsgTxt("trel must in the range [0,20e-3]\n");
    
    /* Calculate number of samples for total repetition time */
    
    totalstim = (int)floor(pxbins/nrep);
    
    px = (double*)mxCalloc(totalstim*nrep,sizeof(double));
    
    /* Put stimulus waveform into pressure waveform */
    
    for (lp=0; lp<pxbins; lp++)
        px[lp] = pxtmp[lp];
    
    /* Create arrays for the return arguments */
    
    outsize[0] = 1;
    outsize[1] = totalstim;
    
    plhs[0] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
    
    outsize[1] = totalstim*nrep;
    plhs[3] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
    
    plhs[4] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
    plhs[5] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
    
    /* Assign pointers to the outputs */
    
    psth	  = mxGetPr(plhs[0]);
    meanrate	  = mxGetPr(plhs[1]);
    varrate = mxGetPr(plhs[2]);
    synout = mxGetPr(plhs[3]);
    trd_vector = mxGetPr(plhs[4]);
    trel_vector = mxGetPr(plhs[5]);
    
    /* run the model */
        
    SingleAN(px,cf,nrep,tdres,totalstim,noiseType,implnt,spont,tabs,trel,meanrate,varrate,psth,synout,trd_vector,trel_vector);
    
    mxFree(px);
    
}

void SingleAN(double *px, double cf, int nrep, double tdres, int totalstim, double noiseType, double implnt, double spont, double tabs, double trel, double *meanrate, double *varrate, double *psth, double *synout, double *trd_vector, double *trel_vector)
{
    
    /*variables for the signal-path, control-path and onward */
    double *sptime;
    double  MeanISI;
    double  SignalLength;
    long    MaxArraySizeSpikes;
    double tau,t_rd_rest,t_rd_init, t_rd_jump,trel_i;
    int    nSites;
    int    i,nspikes,ipst;
    double I;
    double sampFreq = 10e3; /* Sampling frequency used in the synapse */
    double  total_mean_rate;
    /* Declarations of the functions used in the program */
    double Synapse(double *, double, double, int, int, double, double,  double, double, double *);
    int SpikeGenerator(double *,double ,double, double, double,double , int , double , double, double , int , int , double,long,  double *, double *);
        
    /*====== Run the synapse model ======*/
    I = Synapse(px, tdres, cf, totalstim, nrep, spont, noiseType, implnt, sampFreq, synout);
    
    /* Calculate the overall mean synaptic rate */
    total_mean_rate = 0;
    for(i = 0; i<I ; i++)
    {   total_mean_rate= total_mean_rate+ synout[i]/I;
    };
    
    /*======  Synaptic Release/Spike Generation Parameters ======*/
    
    nSites = 4;      /* Number of synpatic release sites */
    
    t_rd_rest = 14.0e-3;   /* Resting value of the mean redocking time */
    t_rd_jump = 0.4e-3;  /* Size of jump in mean redocking time when a redocking event occurs */
    t_rd_init   = t_rd_rest+0.02e-3*spont-t_rd_jump;  /* Initial value of the mean redocking time */
    tau =  60.0e-3;  /* Time constant for short-term adaptation (in mean redocking time) */
    
    
    /* We register only the spikes at times after zero, the sufficient array size (more than 99.7 percent of cases) to register spike times  after zero is :/ */
      /*MaxN=signalLengthInSec/meanISI+ 3*sqrt(signalLengthInSec/MeanISI)= nSpike_average +3*sqrt(nSpike_average)*/
    MeanISI = (1/total_mean_rate)+ (t_rd_init)/nSites+tabs+trel;
    SignalLength = totalstim*nrep*tdres;
    MaxArraySizeSpikes= ceil ((long) (SignalLength/MeanISI + 3* sqrt(SignalLength/MeanISI)));
    
    sptime  = (double*)mxCalloc(MaxArraySizeSpikes,sizeof(double));
    nspikes=0;
    do {
        if  (nspikes<0) /* Deal with cases where the spike time array was not long enough */
        {   mxFree(sptime);
            MaxArraySizeSpikes = MaxArraySizeSpikes+100; /* Make the spike time array 100 elements larger */
            sptime  = (double*)mxCalloc(MaxArraySizeSpikes,sizeof(double));
        }
        
        nspikes =  SpikeGenerator(synout,  tdres, t_rd_rest, t_rd_init, tau, t_rd_jump, nSites, tabs, trel, spont, totalstim, nrep,  total_mean_rate,MaxArraySizeSpikes,  sptime, trd_vector) ;
        
    } while (nspikes<0);  /* Repeat if spike time array was not long enough */
        
    /* Calculate the analytical estimates of meanrate and varrate and wrapping them up based on no. of repetitions */
    for(i = 0; i<I ; i++)
    {
        
        ipst = (int) (fmod(i,totalstim));
        if (synout[i]>0)
        {
            trel_i = __min(trel*100/synout[i],trel);
            trel_vector[i] = trel_i;
            
            meanrate[ipst] = meanrate[ipst] + synout[i]/(synout[i]*(tabs + trd_vector[i]/nSites + trel_i) + 1)/nrep;  /* estimated instantaneous mean rate */
            varrate[ipst]  = varrate[ipst] + ((11*pow(synout[i],7)*pow(trd_vector[i],7))/2 + (3*pow(synout[i],8)*pow(trd_vector[i],8))/16 + 12288*pow(synout[i],2)*pow(trel_i,2) + trd_vector[i]*(22528*pow(synout[i],3)*pow(trel_i,2) + 22528*synout[i]) + pow(trd_vector[i],6)*(3*pow(synout[i],8)*pow(trel_i,2) + 82*pow(synout[i],6)) + pow(trd_vector[i],5)*(88*pow(synout[i],7)*pow(trel_i,2) + 664*pow(synout[i],5)) + pow(trd_vector[i],4)*(976*pow(synout[i],6)*pow(trel_i,2) + 3392*pow(synout[i],4)) + pow(trd_vector[i],3)*(5376*pow(synout[i],5)*pow(trel_i,2) + 10624*pow(synout[i],3)) + pow(trd_vector[i],2)*(15616*pow(synout[i],4)*pow(trel_i,2) + 20992*pow(synout[i],2)) + 12288)/(pow(synout[i],2)*pow(synout[i]*trd_vector[i] + 4,4)*(3*pow(synout[i],2)*pow(trd_vector[i],2) + 40*synout[i]*trd_vector[i] + 48)*pow(trd_vector[i]/4 + tabs + trel_i + 1/synout[i],3))/nrep; /* estimated instananeous variance in the discharge rate */
        }
               else
            trel_vector[i] = trel;
    };
    
    /* Generate PSTH */
    for(i = 0; i < nspikes; i++)
    {
        ipst = (int) (fmod(sptime[i],tdres*totalstim) / tdres);
        psth[ipst] = psth[ipst] + 1;
    };
        
} /* End of the SingleAN function */

/* --------------------------------------------------------------------------------------------*/
double Synapse(double *ihcout, double tdres, double cf, int totalstim, int nrep, double spont, double noiseType, double implnt, double sampFreq, double *synout)
{
    /* Initalize Variables */
    int z, b;
    int resamp = (int) ceil(1/(tdres*sampFreq));
    double incr = 0.0; int delaypoint = (int) floor(7500/(cf/1e3));
    
    double alpha1, beta1, I1, alpha2, beta2, I2, binwidth;
    int    k,j,indx,i;
    
    double cf_factor,cfslope,cfsat,cfconst,multFac;
    
    double *sout1, *sout2, *synSampOut, *powerLawIn, *mappingOut, *TmpSyn;
    double *m1, *m2, *m3, *m4, *m5;
    double *n1, *n2, *n3;
    
    mxArray	*randInputArray[5], *randOutputArray[1];
    double *randNums;
    
    mxArray	*IhcInputArray[3], *IhcOutputArray[1];
    double *sampIHC, *ihcDims;
    
    mappingOut = (double*)mxCalloc((long) ceil(totalstim*nrep),sizeof(double));
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
    alpha1 = 1.5e-6*100e3; beta1 = 5e-4; I1 = 0;
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
    
    mexCallMATLAB(1, randOutputArray, 5, randInputArray, "bruce2018_ffgn");
    randNums = mxGetPr(randOutputArray[0]);
    /*----------------------------------------------------------*/
    /*----- Mapping Function from IHCOUT to input to the PLA ---*/
    /*----------------------------------------------------------*/
    cfslope = pow(spont,0.19)*pow(10,-0.87);
    cfconst = 0.1*pow(log10(spont),2)+0.56*log10(spont)-0.84;
    cfsat = pow(10,(cfslope*8965.5/1e3 + cfconst));
    cf_factor = __min(cfsat,pow(10,cfslope*cf/1e3 + cfconst))*2.0;

    multFac = __max(2.95*__max(1.0,1.5-spont/100),4.3-0.2*cf/1e3);

    k = 0;
    for (indx=0; indx<totalstim*nrep; ++indx)
    {
        mappingOut[k] = pow(10,(0.9*log10(fabs(ihcout[indx])*cf_factor))+ multFac);
        if (ihcout[indx]<0) mappingOut[k] = - mappingOut[k];
        k=k+1;
    }
    for (k=0; k<delaypoint; k++)
        powerLawIn[k] = mappingOut[0]+3.0*spont;
    for (k=delaypoint; k<totalstim*nrep+delaypoint; k++)
        powerLawIn[k] = mappingOut[k-delaypoint]+3.0*spont;
    for (k=totalstim*nrep+delaypoint; k<totalstim*nrep+3*delaypoint; k++)
        powerLawIn[k] = powerLawIn[k-1]+3.0*spont;
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
    
    mxFree(powerLawIn); mxFree(mappingOut);
    /*----------------------------------------------------------*/
    /*----- Running Power-law Adaptation -----------------------*/
    /*----------------------------------------------------------*/
    k = 0;
    for (indx=0; indx<floor((totalstim*nrep+2*delaypoint)*tdres*sampFreq); indx++)
    {
        sout1[k]  = __max( 0, sampIHC[indx] + randNums[indx]- alpha1*I1);
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
    /*-- Upsampling to original (High 100 kHz) sampling rate ---*/
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
        synout[i] = TmpSyn[i+delaypoint];
    
    mxFree(synSampOut); mxFree(TmpSyn);
    mxDestroyArray(randInputArray[0]); mxDestroyArray(randOutputArray[0]);
    mxDestroyArray(IhcInputArray[0]); mxDestroyArray(IhcOutputArray[0]); mxDestroyArray(IhcInputArray[1]); mxDestroyArray(IhcInputArray[2]);
    mxDestroyArray(randInputArray[1]);mxDestroyArray(randInputArray[2]); mxDestroyArray(randInputArray[3]);
    mxDestroyArray(randInputArray[4]);
    return((long) ceil(totalstim*nrep));
}
/* ------------------------------------------------------------ */
/* Pass the output of Synapse model through the Spike Generator */

int SpikeGenerator(double *synout, double tdres, double t_rd_rest, double t_rd_init, double tau, double t_rd_jump, int nSites, double tabs, double trel, double spont, int totalstim, int nrep,double total_mean_rate,long MaxArraySizeSpikes, double *sptime, double *trd_vector)
{
    
    /* Initializing the variables: */
    
    double*  preRelease_initialGuessTimeBins;
    int*     unitRateInterval;
    double*  elapsed_time;
    double*  previous_release_times;
    double*  current_release_times;
    double*  oneSiteRedock;
    double*  Xsum;
    
    double  MeanInterEvents;
    long    MaxArraySizeEvents;
    
    /* Generating a vector of random numbers using mexCallMATLAB */
    mxArray *randInputArray[1], *randOutputArray[1];
    double *randDims, *randNums;
    long randBufIndex;
    long randBufLen;
    
    
    long     spCount; /* total numebr of spikes fired */
    
    long     k;  /*the loop starts from kInit */
    
    int i, siteNo, kInit;
    double Tref, current_refractory_period, trel_k;
    int t_rd_decay, rd_first;
    
    double previous_redocking_period,  current_redocking_period;
    int oneSiteRedock_rounded, elapsed_time_rounded ;
    
    mxArray *sortInputArray[1], *sortOutputArray[1];
    double *sortDims, *preReleaseTimeBinsSorted;
    
    preRelease_initialGuessTimeBins = (double*)mxCalloc(nSites, sizeof(double));
    unitRateInterval                = (int*)mxCalloc(nSites, sizeof(double));;
    elapsed_time                    = (double*)mxCalloc(nSites, sizeof(double));
    previous_release_times          = (double*)mxCalloc(nSites, sizeof(double));
    current_release_times           = (double*)mxCalloc(nSites, sizeof(double));
    oneSiteRedock                   = (double*)mxCalloc(nSites, sizeof(double));
    Xsum                            = (double*)mxCalloc(nSites, sizeof(double));
    
    /* Estimating Max number of spikes and events (including before zero elements)  */
    MeanInterEvents = (1/total_mean_rate)+ (t_rd_init)/nSites;
    /* The sufficient array size (more than 99.7% of cases) to register event times  after zero is :/
     * /*MaxN=signalLengthInSec/meanEvents+ 3*sqrt(signalLengthInSec/MeanEvents)*/
    
    MaxArraySizeEvents= ceil ((long) (totalstim*nrep*tdres/MeanInterEvents+ 3 * sqrt(totalstim*nrep*tdres/MeanInterEvents)))+nSites;
    
    /* Max random array Size:   nSites elements for oneSiteRedock initialization, nSites elements for preRelease_initialGuessTimeBins initialization
     * 1 element for Tref initialization, MaxArraySizeSpikes elements for Tref in the loop, MaxArraySizeEvents elements one  time for redocking, another time for rate intervals
     * Also, for before zero elements, Averageley add 2nSites  events (redock-and unitRate) and add nSites (Max) for Trefs:
     * in total : 3 nSpikes  */
    randBufLen = (long) ceil( 2*nSites+ 1 + MaxArraySizeSpikes + 2*MaxArraySizeEvents + MaxArraySizeSpikes+ 3*nSites);
    
    randInputArray[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
    randDims = mxGetPr(randInputArray[0]);
    randDims[0] = 1;
    randDims[1] = randBufLen;
    mexCallMATLAB(1, randOutputArray, 1, randInputArray, "rand");
    randNums = mxGetPr(randOutputArray[0]);
    randBufIndex = 0;
    
    
    /* Initial < redocking time associated to nSites release sites */
    for (i=0; i<nSites; i++)
    {
        oneSiteRedock[i]=-t_rd_init*log(randNums[randBufIndex++]);
    }
    
    /* Initial  preRelease_initialGuessTimeBins  associated to nsites release sites */
    
    for (i=0; i<nSites; i++)
    {
        preRelease_initialGuessTimeBins[i]= __max(-totalstim*nrep,ceil ((nSites/__max(synout[0],0.1) + t_rd_init)*log(randNums[randBufIndex++] ) / tdres));
        
    }
    
    
    /* Call Sort function using  */
    sortInputArray[0] = mxCreateDoubleMatrix(1, nSites, mxREAL);
    sortDims = mxGetPr(sortInputArray[0]);
    for (i=0;i<nSites; i++)
    {
        sortDims[i] = preRelease_initialGuessTimeBins[i];
        
    }
    
    mexCallMATLAB(1, sortOutputArray, 1, sortInputArray, "sort");
    
    /*Now Sort the four initial preRelease times and associate
     * the farthest to zero as the site which has also generated a spike */
    
    preReleaseTimeBinsSorted =  mxGetPr(sortOutputArray[0]);
    
    /* Consider the inital previous_release_times to be  the preReleaseTimeBinsSorted *tdres */
    for (i=0; i<nSites; i++)
    {
        previous_release_times[i] = ((double)preReleaseTimeBinsSorted[i])*tdres;
    }
    
    /* The position of first spike, also where the process is started- continued from the past */
    kInit = (int) preReleaseTimeBinsSorted[0];
    
    
    /* Current refractory time */
    Tref = tabs - trel*log( randNums[ randBufIndex++ ] );
    
    /*initlal refractory regions */
    current_refractory_period = (double) kInit*tdres;
    
    spCount = 0; /* total numebr of spikes fired */
    k = kInit;  /*the loop starts from kInit */
    
    /* set dynamic mean redocking time to initial mean redocking time  */
    previous_redocking_period = t_rd_init;
    current_redocking_period = previous_redocking_period;
    t_rd_decay = 1; /* Logical "true" as to whether to decay the value of current_redocking_period at the end of the time step */
    rd_first = 0; /* Logical "false" as to whether to a first redocking event has occurred */
    
    /* a loop to find the spike times for all the totalstim*nrep */
    while (k < totalstim*nrep){
        
        for (siteNo = 0; siteNo<nSites; siteNo++)
        {
            
            if ( k > preReleaseTimeBinsSorted [siteNo] )
            {
            
                /* redocking times do not necessarily occur exactly at time step value - calculate the
                 * number of integer steps for the elapsed time and redocking time */
                oneSiteRedock_rounded =  (int) floor(oneSiteRedock[siteNo]/tdres);
                elapsed_time_rounded =  (int) floor(elapsed_time[siteNo]/tdres);
                if ( oneSiteRedock_rounded == elapsed_time_rounded )
                {
                    /* Jump  trd by t_rd_jump if a redocking event has occurred   */
                    current_redocking_period  =   previous_redocking_period  + t_rd_jump;
                    previous_redocking_period =   current_redocking_period;
                    t_rd_decay = 0; /* Don't decay the value of current_redocking_period if a jump has occurred */
                    rd_first = 1; /* Flag for when a jump has first occurred */
                }
                
                /* to be sure that for each site , the code start from its
                 * associated  previus release time :*/
                elapsed_time[siteNo] = elapsed_time[siteNo] + tdres;
            };
            
            
            /*the elapsed time passes  the one time redock (the redocking is finished),
             * In this case the synaptic vesicle starts sensing the input
             * for each site integration starts after the redockinging is finished for the corresponding site)*/
            if ( elapsed_time[siteNo] >= oneSiteRedock [siteNo] )
            {
                Xsum[siteNo] = Xsum[siteNo] + synout[__max(0,k)] / nSites;
                
                /* There are  nSites integrals each vesicle senses 1/nosites of  the whole rate */
            }
            
            
            
            if  ( (Xsum[siteNo]  >=  unitRateInterval[siteNo]) &&  ( k >= preReleaseTimeBinsSorted [siteNo] ) )
            {  /* An event- a release  happened for the siteNo*/
                
                oneSiteRedock[siteNo]  = -current_redocking_period*log( randNums[randBufIndex++]);
                current_release_times[siteNo] = previous_release_times[siteNo]  + elapsed_time[siteNo];
                elapsed_time[siteNo] = 0;               
                
                if ( (current_release_times[siteNo] >= current_refractory_period) )
                {  /* A spike occured for the current event- release
                 * spike_times[(int)(current_release_times[siteNo]/tdres)-kInit+1 ] = 1;*/
                    
                    /*Register only non negative spike times */
                    if (current_release_times[siteNo] >= 0)
                    {
                        sptime[spCount] = current_release_times[siteNo]; spCount = spCount + 1;
                    }
                    
                    trel_k = __min(trel*100/synout[__max(0,k)],trel);

                    Tref = tabs-trel_k*log( randNums[randBufIndex++] );   /*Refractory periods */
                    
                    current_refractory_period = current_release_times[siteNo] + Tref;
                    
                }
                
                previous_release_times[siteNo] = current_release_times[siteNo];
                
                Xsum[siteNo] = 0;
                unitRateInterval[siteNo] = (int) (-log(randNums[randBufIndex++]) / tdres);
                
            };
            /* Error Catching */
            if ( (spCount+1)>MaxArraySizeSpikes  || (randBufIndex+1 )>randBufLen  )
            {     /* mexPrintf ("--------Array for spike times or random Buffer length not large enough, Rerunning the function.-----\n\n"); */
                spCount = -1;
                k = totalstim*nrep;
                siteNo = nSites;
            }
            
        };
        
        /* Decay the adapative mean redocking time towards the resting value if no redocking events occurred in this time step */
        if ( (t_rd_decay==1) && (rd_first==1) )
        {
            current_redocking_period =   previous_redocking_period  - (tdres/tau)*( previous_redocking_period-t_rd_rest );
            previous_redocking_period =  current_redocking_period;
        }
        else
        {
            t_rd_decay = 1;
        }
        
        /* Store the value of the adaptive mean redocking time if it is within the simulation output period */
        if ((k>=0)&&(k<totalstim*nrep))
            trd_vector [k] = current_redocking_period;
        
        k = k+1;
        
        
    };
    
    mxFree(preRelease_initialGuessTimeBins);
    mxFree(unitRateInterval);
    mxFree(elapsed_time);
    mxFree(previous_release_times);
    mxFree(current_release_times);
    mxFree(oneSiteRedock);
    mxFree(Xsum);
    mxDestroyArray(randInputArray[0]); mxDestroyArray(randOutputArray[0]);
    mxDestroyArray(sortInputArray[0]); mxDestroyArray(sortOutputArray[0]);
    return (spCount);
    
}


