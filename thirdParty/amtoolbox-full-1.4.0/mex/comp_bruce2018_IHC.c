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
	
	double *px, cf, tdres, reptime, cohc, cihc;
	int    nrep, pxbins, lp, totalstim, species;
    mwSize outsize[2];
    
	double *pxtmp, *cftmp, *nreptmp, *tdrestmp, *reptimetmp, *cohctmp, *cihctmp, *speciestmp;
    double *ihcout,  *c1out, *c2out;
   
	void   IHCAN(double *, double, int, double, int, double, double, int, double *, double *, double *);
	
	/* Check for proper number of arguments */
	
	if (nrhs != 8) 
	{
		mexErrMsgTxt("model_IHC requires 8 input arguments.");
	}; 

	if (nlhs !=3)  
	{
		mexErrMsgTxt("model_IHC requires 3 output argument.");
	};
	
	/* Assign pointers to the inputs */

	pxtmp		= mxGetPr(prhs[0]);
	cftmp		= mxGetPr(prhs[1]);
	nreptmp		= mxGetPr(prhs[2]);
	tdrestmp	= mxGetPr(prhs[3]);
	reptimetmp	= mxGetPr(prhs[4]);
    cohctmp		= mxGetPr(prhs[5]);
    cihctmp		= mxGetPr(prhs[6]);
    speciestmp	= mxGetPr(prhs[7]);
		
	/* Check individual input arguments */

	pxbins = (int) mxGetN(prhs[0]);
	if (pxbins==1)
		mexErrMsgTxt("px must be a row vector\n");

    species = (int) speciestmp[0];
	if (speciestmp[0]!=species)
		mexErrMsgTxt("species must an integer.\n");
	if (species<1 || species>3)
		mexErrMsgTxt("Species must be 1 for cat, or 2 or 3 for human.\n");
 
    
	cf = cftmp[0];
    if (species==1)
    {
	if ((cf<124.9)|(cf>40.1e3))
	{
		mexPrintf("cf (= %1.1f Hz) must be between 125 Hz and 40 kHz for cat model\n",cf);
		mexErrMsgTxt("\n");
    }
    }
    if (species>1)
    {
  	if ((cf<124.9)|(cf>20.1e3))
	{
		mexPrintf("cf (= %1.1f Hz) must be between 125 Hz and 20 kHz for human model\n",cf);
		mexErrMsgTxt("\n");
    }
    }
    
	nrep = (int) nreptmp[0];
	if (nreptmp[0]!=nrep)
		mexErrMsgTxt("nrep must an integer.\n");
	if (nrep<1)
		mexErrMsgTxt("nrep must be greater that 0.\n");

    tdres = tdrestmp[0];
	
	reptime = reptimetmp[0];
	if (reptime<pxbins*tdres)  /* duration of stimulus = pxbins*tdres */
		mexErrMsgTxt("reptime should be equal to or longer than the stimulus duration.\n");

    cohc = cohctmp[0]; /* impairment in the OHC  */
	if ((cohc<0)|(cohc>1))
	{
		mexPrintf("cohc (= %1.1f) must be between 0 and 1\n",cohc);
		mexErrMsgTxt("\n");
	}

	cihc = cihctmp[0]; /* impairment in the IHC  */
	if ((cihc<0)|(cihc>1))
	{
		mexPrintf("cihc (= %1.1f) must be between 0 and 1\n",cihc);
		mexErrMsgTxt("\n");
	}
	
   
	/* Calculate number of samples for total repetition time */

	/*totalstim = (int)floor((reptime*1e3)/(tdres*1e3)); */ /*older definition*/
    totalstim = (int)floor(reptime/tdres+0.5);

    px = (double*)mxCalloc(totalstim,sizeof(double)); 

	/* Put stimulus waveform into pressure waveform */

	for (lp=0; lp<pxbins; lp++)
			px[lp] = pxtmp[lp];
	
	/* Create an array for the return argument */
	
    outsize[0] = 1;
	outsize[1] = totalstim*nrep;
    
	plhs[0] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);
  plhs[2] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);		
	/* Assign pointers to the outputs */
	
	ihcout  = mxGetPr(plhs[0]);
  c1out = mxGetPr(plhs[1]);
  c2out = mxGetPr(plhs[2]);		
	/* run the model */

	IHCAN(px,cf,nrep,tdres,totalstim,cohc,cihc,species,ihcout, c1out, c2out);

 mxFree(px);

}

void IHCAN(double *px, double cf, int nrep, double tdres, int totalstim,
                double cohc, double cihc, int species, double *ihcout, double *c1out, double *c2out)
{	
    
    /*variables for middle-ear model */
	double megainmax;
    double *mey1, *mey2, *mey3, meout,c1filterouttmp,c2filterouttmp,c1vihctmp,c2vihctmp;
    double fp,C,m11,m12,m13,m14,m15,m16,m21,m22,m23,m24,m25,m26,m31,m32,m33,m34,m35,m36;

	/*variables for the signal-path, control-path and onward */
	double *ihcouttmp,*tmpgain , *c1outtmp, *c2outtmp;
	int    grd;

    double bmplace,centerfreq,gain,taubm,ratiowb,bmTaubm,fcohc,TauWBMax,TauWBMin,tauwb;
    double Taumin[1],Taumax[1],bmTaumin[1],bmTaumax[1],ratiobm[1],lasttmpgain,wbgain,ohcasym,ihcasym,delay;
	int    i,n,delaypoint,grdelay[1],bmorder,wborder;
	double wbout1,wbout,ohcnonlinout,ohcout,tmptauc1,tauc1,rsigma,wb_gain;
            
    /* Declarations of the functions used in the program */
	double C1ChirpFilt(double, double,double, int, double, double);
	double C2ChirpFilt(double, double,double, int, double, double);
    double WbGammaTone(double, double, double, int, double, double, int);

    double Get_tauwb(double, int, int, double *, double *);
	double Get_taubm(double, int, double, double *, double *, double *);
    double gain_groupdelay(double, double, double, double, int *);
    double delay_cat(double cf);
    double delay_human(double cf);

    double OhcLowPass(double, double, double, int, double, int);
    double IhcLowPass(double, double, double, int, double, int);
	double Boltzman(double, double, double, double, double);
    double NLafterohc(double, double, double, double);
	double ControlSignal(double, double, double, double, double);

    double NLogarithm(double, double, double, double);
    
    /* Allocate dynamic memory for the temporary variables */
	ihcouttmp  = (double*)mxCalloc(totalstim*nrep,sizeof(double));
	c1outtmp  = (double*)mxCalloc(totalstim*nrep,sizeof(double));
	c2outtmp  = (double*)mxCalloc(totalstim*nrep,sizeof(double));

	    
	mey1 = (double*)mxCalloc(totalstim,sizeof(double));
	mey2 = (double*)mxCalloc(totalstim,sizeof(double));
	mey3 = (double*)mxCalloc(totalstim,sizeof(double));

	tmpgain = (double*)mxCalloc(totalstim,sizeof(double));
    
	/** Calculate the center frequency for the control-path wideband filter
	    from the location on basilar membrane, based on Greenwood (JASA 1990) */

	if (species==1) /* for cat */
    {
        /* Cat frequency shift corresponding to 1.2 mm */
        bmplace = 11.9 * log10(0.80 + cf / 456.0); /* Calculate the location on basilar membrane from CF */
        centerfreq = 456.0*(pow(10,(bmplace+1.2)/11.9)-0.80); /* shift the center freq */
    }

	if (species>1) /* for human */
    {
        /* Human frequency shift corresponding to 1.2 mm */
        bmplace = (35/2.1) * log10(1.0 + cf / 165.4); /* Calculate the location on basilar membrane from CF */
        centerfreq = 165.4*(pow(10,(bmplace+1.2)/(35/2.1))-1.0); /* shift the center freq */
    }
    
	/*==================================================================*/
	/*====== Parameters for the gain ===========*/
    
	if(species==1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for cat */
    if(species>1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for human */
    /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/
    if(gain>60.0) gain = 60.0;  
    if(gain<15.0) gain = 15.0;
    
	/*====== Parameters for the control-path wideband filter =======*/
	bmorder = 3;
	Get_tauwb(cf,species,bmorder,Taumax,Taumin);
	taubm   = cohc*(Taumax[0]-Taumin[0])+Taumin[0];
	ratiowb = Taumin[0]/Taumax[0];
	/*====== Parameters for the signal-path C1 filter ======*/
	Get_taubm(cf,species,Taumax[0],bmTaumax,bmTaumin,ratiobm);
	bmTaubm  = cohc*(bmTaumax[0]-bmTaumin[0])+bmTaumin[0];
	fcohc    = bmTaumax[0]/bmTaubm;
    /*====== Parameters for the control-path wideband filter =======*/
	wborder  = 3;
    TauWBMax = Taumin[0]+0.2*(Taumax[0]-Taumin[0]);
	TauWBMin = TauWBMax/Taumax[0]*Taumin[0];
    tauwb    = TauWBMax+(bmTaubm-bmTaumax[0])*(TauWBMax-TauWBMin)/(bmTaumax[0]-bmTaumin[0]);
	
	wbgain = gain_groupdelay(tdres,centerfreq,cf,tauwb,grdelay);
	tmpgain[0]   = wbgain; 
	lasttmpgain  = wbgain;
  	/*===============================================================*/
    /* Nonlinear asymmetry of OHC function and IHC C1 transduction function*/
	ohcasym  = 7.0;    
	ihcasym  = 3.0;
  	/*===============================================================*/
    /*===============================================================*/
    /* Prewarping and related constants for the middle ear */
     fp = 1e3;  /* prewarping frequency 1 kHz */
     C  = TWOPI*fp/tan(TWOPI/2*fp*tdres);
     if (species==1) /* for cat */
     {
         /* Cat middle-ear filter - simplified version from Bruce et al. (JASA 2003) */
         m11 = C/(C + 693.48);                    m12 = (693.48 - C)/C;            m13 = 0.0;
         m14 = 1.0;                               m15 = -1.0;                      m16 = 0.0;
         m21 = 1/(pow(C,2) + 11053*C + 1.163e8);  m22 = -2*pow(C,2) + 2.326e8;     m23 = pow(C,2) - 11053*C + 1.163e8; 
         m24 = pow(C,2) + 1356.3*C + 7.4417e8;    m25 = -2*pow(C,2) + 14.8834e8;   m26 = pow(C,2) - 1356.3*C + 7.4417e8;
         m31 = 1/(pow(C,2) + 4620*C + 909059944); m32 = -2*pow(C,2) + 2*909059944; m33 = pow(C,2) - 4620*C + 909059944;
         m34 = 5.7585e5*C + 7.1665e7;             m35 = 14.333e7;                  m36 = 7.1665e7 - 5.7585e5*C;
         megainmax=41.1405;
     };
     if (species>1) /* for human */
     {
         /* Human middle-ear filter - based on Pascal et al. (JASA 1998)  */
         m11=1/(pow(C,2)+5.9761e+003*C+2.5255e+007);m12=(-2*pow(C,2)+2*2.5255e+007);m13=(pow(C,2)-5.9761e+003*C+2.5255e+007);m14=(pow(C,2)+5.6665e+003*C);             m15=-2*pow(C,2);					m16=(pow(C,2)-5.6665e+003*C);
         m21=1/(pow(C,2)+6.4255e+003*C+1.3975e+008);m22=(-2*pow(C,2)+2*1.3975e+008);m23=(pow(C,2)-6.4255e+003*C+1.3975e+008);m24=(pow(C,2)+5.8934e+003*C+1.7926e+008); m25=(-2*pow(C,2)+2*1.7926e+008);	m26=(pow(C,2)-5.8934e+003*C+1.7926e+008);
         m31=1/(pow(C,2)+2.4891e+004*C+1.2700e+009);m32=(-2*pow(C,2)+2*1.2700e+009);m33=(pow(C,2)-2.4891e+004*C+1.2700e+009);m34=(3.1137e+003*C+6.9768e+008);     m35=2*6.9768e+008;				m36=(-3.1137e+003*C+6.9768e+008);
         megainmax=2;
     };
  	for (n=0;n<totalstim;n++) /* Start of the loop */
    {    
        if (n==0)  /* Start of the middle-ear filtering section  */
		{
	    	mey1[0]  = m11*px[0];
            if (species>1) mey1[0] = m11*m14*px[0];
            mey2[0]  = mey1[0]*m24*m21;
            mey3[0]  = mey2[0]*m34*m31;
            meout = mey3[0]/megainmax ;
        }
            
        else if (n==1)
		{
            mey1[1]  = m11*(-m12*mey1[0] + px[1]       - px[0]);
            if (species>1) mey1[1] = m11*(-m12*mey1[0]+m14*px[1]+m15*px[0]);
			mey2[1]  = m21*(-m22*mey2[0] + m24*mey1[1] + m25*mey1[0]);
            mey3[1]  = m31*(-m32*mey3[0] + m34*mey2[1] + m35*mey2[0]);
            meout = mey3[1]/megainmax;
		}
	    else 
		{
            mey1[n]  = m11*(-m12*mey1[n-1]  + px[n]         - px[n-1]);
            if (species>1) mey1[n]= m11*(-m12*mey1[n-1]-m13*mey1[n-2]+m14*px[n]+m15*px[n-1]+m16*px[n-2]);
            mey2[n]  = m21*(-m22*mey2[n-1] - m23*mey2[n-2] + m24*mey1[n] + m25*mey1[n-1] + m26*mey1[n-2]);
            mey3[n]  = m31*(-m32*mey3[n-1] - m33*mey3[n-2] + m34*mey2[n] + m35*mey2[n-1] + m36*mey2[n-2]);
            meout = mey3[n]/megainmax;
		}; 	/* End of the middle-ear filtering section */   
     
		/* Control-path filter */

        wbout1 = WbGammaTone(meout,tdres,centerfreq,n,tauwb,wbgain,wborder);
        wbout  = pow((tauwb/TauWBMax),wborder)*wbout1*10e3*__max(1,cf/5e3);
  
        ohcnonlinout = Boltzman(wbout,ohcasym,12.0,5.0,5.0); /* pass the control signal through OHC Nonlinear Function */
		ohcout = OhcLowPass(ohcnonlinout,tdres,600,n,1.0,2);/* lowpass filtering after the OHC nonlinearity */
        
		tmptauc1 = NLafterohc(ohcout,bmTaumin[0],bmTaumax[0],ohcasym); /* nonlinear function after OHC low-pass filter */
		tauc1    = cohc*(tmptauc1-bmTaumin[0])+bmTaumin[0];  /* time -constant for the signal-path C1 filter */
		rsigma   = 1/tauc1-1/bmTaumax[0]; /* shift of the location of poles of the C1 filter from the initial positions */

		if (1/tauc1<0.0) mexErrMsgTxt("The poles are in the right-half plane; system is unstable.\n");

		tauwb = TauWBMax+(tauc1-bmTaumax[0])*(TauWBMax-TauWBMin)/(bmTaumax[0]-bmTaumin[0]);

	    wb_gain = gain_groupdelay(tdres,centerfreq,cf,tauwb,grdelay);
		
		grd = grdelay[0]; 

        if ((grd+n)<totalstim)
	         tmpgain[grd+n] = wb_gain;

        if (tmpgain[n] == 0)
			tmpgain[n] = lasttmpgain;	
		
		wbgain      = tmpgain[n];
		lasttmpgain = wbgain;
	 		        
        /*====== Signal-path C1 filter ======*/
         
		c1filterouttmp = C1ChirpFilt(meout, tdres, cf, n, bmTaumax[0], rsigma); /* C1 filter output */
		c1out[n] = c1filterouttmp;
	 
        /*====== Parallel-path C2 filter ======*/

		c2filterouttmp  = C2ChirpFilt(meout, tdres, cf, n, bmTaumax[0], 1/ratiobm[0]); /* parallel-filter output*/
		c2out[n] = c2filterouttmp;
	    /*=== Run the inner hair cell (IHC) section: NL function and then lowpass filtering ===*/

        c1vihctmp  = NLogarithm(cihc*c1filterouttmp,0.1,ihcasym,cf);
	     
		c2vihctmp = -NLogarithm(c2filterouttmp*fabs(c2filterouttmp)*cf/10*cf/2e3,0.2,1.0,cf); /* C2 transduction output */
            
        ihcouttmp[n] = IhcLowPass(c1vihctmp+c2vihctmp,tdres,3000,n,1.0,7);
   };  /* End of the loop */
   
    /* Stretched out the IHC output according to nrep (number of repetitions) */
   
    for(i=0;i<totalstim*nrep;i++)
	{
		ihcouttmp[i] = ihcouttmp[(int) (fmod(i,totalstim))];
  	};   
   	/* Adjust total path delay to IHC output signal */
    if (species==1)
        delay      = delay_cat(cf);
    if (species>1)
    {/*    delay      = delay_human(cf); */
        delay      = delay_cat(cf); /* signal delay changed back to cat function for version 5.2 */
    };
    delaypoint =__max(0,(int) ceil(delay/tdres));    
         
    for(i=delaypoint;i<totalstim*nrep;i++)
	{        
		ihcout[i] = ihcouttmp[i - delaypoint];
  	};   

    /* Freeing dynamic memory allocated earlier */

    mxFree(ihcouttmp);
    mxFree(mey1); mxFree(mey2); mxFree(mey3);	
    mxFree(tmpgain);

} /* End of the SingleAN function */
/* -------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------- */
/** Get TauMax, TauMin for the tuning filter. The TauMax is determined by the bandwidth/Q10
    of the tuning filter at low level. The TauMin is determined by the gain change between high
    and low level */

double Get_tauwb(double cf, int species, int order, double *taumax,double *taumin)
{
  double Q10,bw,gain,ratio;
    
  if(species==1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for cat */
  if(species>1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for human */
  /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/ /* older values */

  if(gain>60.0) gain = 60.0;  
  if(gain<15.0) gain = 15.0;
   
  ratio = pow(10,(-gain/(20.0*order)));       /* ratio of TauMin/TauMax according to the gain, order */
  if (species==1) /* cat Q10 values */
  {
    Q10 = pow(10,0.4708*log10(cf/1e3)+0.4664);
  }
  if (species==2) /* human Q10 values from Shera et al. (PNAS 2002) */
  {
    Q10 = pow((cf/1000),0.3)*12.7*0.505+0.2085;
  }
  if (species==3) /* human Q10 values from Glasberg & Moore (Hear. Res. 1990) */
  {
    Q10 = cf/24.7/(4.37*(cf/1000)+1)*0.505+0.2085;
  }
  bw     = cf/Q10;
  taumax[0] = 2.0/(TWOPI*bw);
   
  taumin[0]   = taumax[0]*ratio;
  
  return 0;
}
/* -------------------------------------------------------------------------------------------- */
double Get_taubm(double cf, int species, double taumax,double *bmTaumax,double *bmTaumin, double *ratio)
{
  double gain,factor,bwfactor;
    
  if(species==1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for cat */
  if(species>1) gain = 52.0/2.0*(tanh(2.2*log10(cf/0.6e3)+0.15)+1.0); /* for human */
  /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/ /* older values */

 
  if(gain>60.0) gain = 60.0;  
  if(gain<15.0) gain = 15.0;

  bwfactor = 0.7;
  factor   = 2.5;

  ratio[0]  = pow(10,(-gain/(20.0*factor))); 

  bmTaumax[0] = taumax/bwfactor;
  bmTaumin[0] = bmTaumax[0]*ratio[0];     
  return 0;
}
/* -------------------------------------------------------------------------------------------- */
/** Pass the signal through the signal-path C1 Tenth Order Nonlinear Chirp-Gammatone Filter */

double C1ChirpFilt(double x, double tdres,double cf, int n, double taumax, double rsigma)
{
    static double C1gain_norm, C1initphase; 
    static double C1input[12][4], C1output[12][4];

    double ipw, ipb, rpa, pzero, rzero;
	double sigma0,fs_bilinear,CF,norm_gain,phase,c1filterout;
	int i,r,order_of_pole,half_order_pole,order_of_zero;
	double temp, dy, preal, pimg;

	COMPLEX p[11]; 
	
	/* Defining initial locations of the poles and zeros */
	/*======== setup the locations of poles and zeros =======*/
	  sigma0 = 1/taumax;
	  ipw    = 1.01*cf*TWOPI-50;
	  ipb    = 0.2343*TWOPI*cf-1104;
	  rpa    = pow(10, log10(cf)*0.9 + 0.55)+ 2000;
	  pzero  = pow(10,log10(cf)*0.7+1.6)+500;

	/*===============================================================*/     
         
     order_of_pole    = 10;             
     half_order_pole  = order_of_pole/2;
     order_of_zero    = half_order_pole;

	 fs_bilinear = TWOPI*cf/tan(TWOPI*cf*tdres/2);
     rzero       = -pzero;
	 CF          = TWOPI*cf;
   
   if (n==0)
   {		  
	p[1].x = -sigma0;     

    p[1].y = ipw;

	p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

    p[2]   = compconj(p[1]);    p[4] = compconj(p[3]); p[6] = compconj(p[5]);

    p[7]   = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

	   C1initphase = 0.0;
       for (i=1;i<=half_order_pole;i++)          
	   {
           preal     = p[i*2-1].x;
		   pimg      = p[i*2-1].y;
	       C1initphase = C1initphase + atan(CF/(-rzero))-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	   };

	/*===================== Initialize C1input & C1output =====================*/

      for (i=1;i<=(half_order_pole+1);i++)          
      {
		   C1input[i][3] = 0; 
		   C1input[i][2] = 0; 
		   C1input[i][1] = 0;
		   C1output[i][3] = 0; 
		   C1output[i][2] = 0; 
		   C1output[i][1] = 0;
      }

	/*===================== normalize the gain =====================*/
    
      C1gain_norm = 1.0;
      for (r=1; r<=order_of_pole; r++)
		   C1gain_norm = C1gain_norm*(pow((CF - p[r].y),2) + p[r].x*p[r].x);
      
   };
     
    norm_gain= sqrt(C1gain_norm)/pow(sqrt(CF*CF+rzero*rzero),order_of_zero);
	
	p[1].x = -sigma0 - rsigma;

	if (p[1].x>0.0) mexErrMsgTxt("The system becomes unstable.\n");
	
	p[1].y = ipw;

	p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

    p[2] = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);

    p[7] = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

    phase = 0.0;
    for (i=1;i<=half_order_pole;i++)          
    {
           preal = p[i*2-1].x;
		   pimg  = p[i*2-1].y;
	       phase = phase-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	};

	rzero = -CF/tan((C1initphase-phase)/order_of_zero);

    if (rzero>0.0) mexErrMsgTxt("The zeros are in the right-half plane.\n");
	 
   /*%==================================================  */
	/*each loop below is for a pair of poles and one zero */
   /*%      time loop begins here                         */
   /*%==================================================  */
 
       C1input[1][3]=C1input[1][2]; 
	   C1input[1][2]=C1input[1][1]; 
	   C1input[1][1]= x;

       for (i=1;i<=half_order_pole;i++)          
       {
           preal = p[i*2-1].x;
		   pimg  = p[i*2-1].y;
		  	   
           temp  = pow((fs_bilinear-preal),2)+ pow(pimg,2);
		   

           /*dy = (input[i][1] + (1-(fs_bilinear+rzero)/(fs_bilinear-rzero))*input[i][2]
                                 - (fs_bilinear+rzero)/(fs_bilinear-rzero)*input[i][3] );
           dy = dy+2*output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg);

           dy = dy-output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);*/
		   
	       dy = C1input[i][1]*(fs_bilinear-rzero) - 2*rzero*C1input[i][2] - (fs_bilinear+rzero)*C1input[i][3]
                 +2*C1output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg)
			     -C1output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);

		   dy = dy/temp;

		   C1input[i+1][3] = C1output[i][2]; 
		   C1input[i+1][2] = C1output[i][1]; 
		   C1input[i+1][1] = dy;

		   C1output[i][2] = C1output[i][1]; 
		   C1output[i][1] = dy;
       }

	   dy = C1output[half_order_pole][1]*norm_gain;  /* don't forget the gain term */
	   c1filterout= dy/4.0;   /* signal path output is divided by 4 to give correct C1 filter gain */
	                   
     return (c1filterout);
}  

/* -------------------------------------------------------------------------------------------- */
/** Parallelpath C2 filter: same as the signal-path C1 filter with the OHC completely impaired */

double C2ChirpFilt(double xx, double tdres,double cf, int n, double taumax, double fcohc)
{
	static double C2gain_norm, C2initphase;
    static double C2input[12][4];  static double C2output[12][4];
   
	double ipw, ipb, rpa, pzero, rzero;

	double sigma0,fs_bilinear,CF,norm_gain,phase,c2filterout;
	int    i,r,order_of_pole,half_order_pole,order_of_zero;
	double temp, dy, preal, pimg;

	COMPLEX p[11]; 	
    
    /*================ setup the locations of poles and zeros =======*/

	  sigma0 = 1/taumax;
	  ipw    = 1.01*cf*TWOPI-50;
      ipb    = 0.2343*TWOPI*cf-1104;
	  rpa    = pow(10, log10(cf)*0.9 + 0.55)+ 2000;
	  pzero  = pow(10,log10(cf)*0.7+1.6)+500;
	/*===============================================================*/     
         
     order_of_pole    = 10;             
     half_order_pole  = order_of_pole/2;
     order_of_zero    = half_order_pole;

	 fs_bilinear = TWOPI*cf/tan(TWOPI*cf*tdres/2);
     rzero       = -pzero;
	 CF          = TWOPI*cf;
   	    
    if (n==0)
    {		  
	p[1].x = -sigma0;     

    p[1].y = ipw;

	p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

    p[2] = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);

    p[7] = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

	   C2initphase = 0.0;
       for (i=1;i<=half_order_pole;i++)         
	   {
           preal     = p[i*2-1].x;
		   pimg      = p[i*2-1].y;
	       C2initphase = C2initphase + atan(CF/(-rzero))-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	   };

	/*===================== Initialize C2input & C2output =====================*/

      for (i=1;i<=(half_order_pole+1);i++)          
      {
		   C2input[i][3] = 0; 
		   C2input[i][2] = 0; 
		   C2input[i][1] = 0;
		   C2output[i][3] = 0; 
		   C2output[i][2] = 0; 
		   C2output[i][1] = 0;
      }
    
    /*===================== normalize the gain =====================*/
    
     C2gain_norm = 1.0;
     for (r=1; r<=order_of_pole; r++)
		   C2gain_norm = C2gain_norm*(pow((CF - p[r].y),2) + p[r].x*p[r].x);
    };
     
    norm_gain= sqrt(C2gain_norm)/pow(sqrt(CF*CF+rzero*rzero),order_of_zero);
    
	p[1].x = -sigma0*fcohc;

	if (p[1].x>0.0) mexErrMsgTxt("The system becomes unstable.\n");
	
	p[1].y = ipw;

	p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

    p[2] = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);

    p[7] = p[1]; p[8] = p[2]; p[9] = p[5]; p[10]= p[6];

    phase = 0.0;
    for (i=1;i<=half_order_pole;i++)          
    {
           preal = p[i*2-1].x;
		   pimg  = p[i*2-1].y;
	       phase = phase-atan((CF-pimg)/(-preal))-atan((CF+pimg)/(-preal));
	};

	rzero = -CF/tan((C2initphase-phase)/order_of_zero);	
    if (rzero>0.0) mexErrMsgTxt("The zeros are in the right-hand plane.\n");
   /*%==================================================  */
   /*%      time loop begins here                         */
   /*%==================================================  */

       C2input[1][3]=C2input[1][2]; 
	   C2input[1][2]=C2input[1][1]; 
	   C2input[1][1]= xx;

      for (i=1;i<=half_order_pole;i++)          
      {
           preal = p[i*2-1].x;
		   pimg  = p[i*2-1].y;
		  	   
           temp  = pow((fs_bilinear-preal),2)+ pow(pimg,2);
		   
           /*dy = (input[i][1] + (1-(fs_bilinear+rzero)/(fs_bilinear-rzero))*input[i][2]
                                 - (fs_bilinear+rzero)/(fs_bilinear-rzero)*input[i][3] );
           dy = dy+2*output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg);

           dy = dy-output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);*/
		   
	      dy = C2input[i][1]*(fs_bilinear-rzero) - 2*rzero*C2input[i][2] - (fs_bilinear+rzero)*C2input[i][3]
                 +2*C2output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg)
			     -C2output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);

		   dy = dy/temp;

		   C2input[i+1][3] = C2output[i][2]; 
		   C2input[i+1][2] = C2output[i][1]; 
		   C2input[i+1][1] = dy;

		   C2output[i][2] = C2output[i][1]; 
		   C2output[i][1] = dy;

       };

	  dy = C2output[half_order_pole][1]*norm_gain;
	  c2filterout= dy/4.0;
	  
	  return (c2filterout); 
}   

/* -------------------------------------------------------------------------------------------- */
/** Pass the signal through the Control path Third Order Nonlinear Gammatone Filter */

double WbGammaTone(double x,double tdres,double centerfreq, int n, double tau,double gain,int order)
{
  static double wbphase;
  static COMPLEX wbgtf[4], wbgtfl[4];

  double delta_phase,dtmp,c1LP,c2LP,out;
  int i,j;
  
  if (n==0)
  {
      wbphase = 0;
      for(i=0; i<=order;i++)
      {
            wbgtfl[i] = compmult(0,compexp(0));
            wbgtf[i]  = compmult(0,compexp(0));
      }
  }
  
  delta_phase = -TWOPI*centerfreq*tdres;
  wbphase += delta_phase;
  
  dtmp = tau*2.0/tdres;
  c1LP = (dtmp-1)/(dtmp+1);
  c2LP = 1.0/(dtmp+1);
  wbgtf[0] = compmult(x,compexp(wbphase));                 /* FREQUENCY SHIFT */
  
  for(j = 1; j <= order; j++)                              /* IIR Bilinear transformation LPF */
  wbgtf[j] = comp2sum(compmult(c2LP*gain,comp2sum(wbgtf[j-1],wbgtfl[j-1])),
      compmult(c1LP,wbgtfl[j]));
  out = REAL(compprod(compexp(-wbphase), wbgtf[order])); /* FREQ SHIFT BACK UP */
  
  for(i=0; i<=order;i++) wbgtfl[i] = wbgtf[i];
  return(out);
}

/* -------------------------------------------------------------------------------------------- */
/** Calculate the gain and group delay for the Control path Filter */

double gain_groupdelay(double tdres,double centerfreq, double cf, double tau,int *grdelay)
{ 
  double tmpcos,dtmp2,c1LP,c2LP,tmp1,tmp2,wb_gain;

  tmpcos = cos(TWOPI*(centerfreq-cf)*tdres);
  dtmp2 = tau*2.0/tdres;
  c1LP = (dtmp2-1)/(dtmp2+1);
  c2LP = 1.0/(dtmp2+1);
  tmp1 = 1+c1LP*c1LP-2*c1LP*tmpcos;
  tmp2 = 2*c2LP*c2LP*(1+tmpcos);
  
  wb_gain = pow(tmp1/tmp2, 1.0/2.0);
  
  grdelay[0] = (int)floor((0.5-(c1LP*c1LP-c1LP*tmpcos)/(1+c1LP*c1LP-2*c1LP*tmpcos)));

  return(wb_gain);
}
/* -------------------------------------------------------------------------------------------- */
/** Calculate the delay (basilar membrane, synapse, etc. for cat) */
double delay_cat(double cf)
{  
  double A0,A1,x,delay;

  A0    = 3.0;  
  A1    = 12.5;
  x     = 11.9 * log10(0.80 + cf / 456.0);      /* cat mapping */
  delay = A0 * exp( -x/A1 ) * 1e-3;
  
  return(delay);
}

/* Calculate the delay (basilar membrane, synapse, etc.) for human, based
        on Harte et al. (JASA 2009) */
double delay_human(double cf) 
{  
  double A,B,delay;

  A    = -0.37;  
  B    = 11.09/2;
  delay = B * pow(cf * 1e-3,A)*1e-3;
  
  return(delay);
}

/* -------------------------------------------------------------------------------------------- */
/* Get the output of the OHC Nonlinear Function (Boltzman Function) */

double Boltzman(double x, double asym, double s0, double s1, double x1)
  {
	double shift,x0,out1,out;

    shift = 1.0/(1.0+asym);  /* asym is the ratio of positive Max to negative Max*/
    x0    = s0*log((1.0/shift-1)/(1+exp(x1/s1)));
	    
    out1 = 1.0/(1.0+exp(-(x-x0)/s0)*(1.0+exp(-(x-x1)/s1)))-shift;
	out = out1/(1-shift);

    return(out);
  }  /* output of the nonlinear function, the output is normalized with maximum value of 1 */
  
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the OHC Low Pass Filter in the Control path */

double OhcLowPass(double x,double tdres,double Fc, int n,double gain,int order)
{
  static double ohc[4],ohcl[4];

  double c,c1LP,c2LP;
  int i,j;

  if (n==0)
  {
      for(i=0; i<(order+1);i++)
      {
          ohc[i] = 0;
          ohcl[i] = 0;
      }
  }    
  
  c = 2.0/tdres;
  c1LP = ( c - TWOPI*Fc ) / ( c + TWOPI*Fc );
  c2LP = TWOPI*Fc / (TWOPI*Fc + c);
  
  ohc[0] = x*gain;
  for(i=0; i<order;i++)
    ohc[i+1] = c1LP*ohcl[i+1] + c2LP*(ohc[i]+ohcl[i]);
  for(j=0; j<=order;j++) ohcl[j] = ohc[j];
  return(ohc[order]);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the IHC Low Pass Filter  */

double IhcLowPass(double x,double tdres,double Fc, int n,double gain,int order)
{
  static double ihc[8],ihcl[8];
  
  double C,c1LP,c2LP;
  int i,j;

  if (n==0)
  {
      for(i=0; i<(order+1);i++)
      {
          ihc[i] = 0;
          ihcl[i] = 0;
      }
  }     
  
  C = 2.0/tdres;
  c1LP = ( C - TWOPI*Fc ) / ( C + TWOPI*Fc );
  c2LP = TWOPI*Fc / (TWOPI*Fc + C);
  
  ihc[0] = x*gain;
  for(i=0; i<order;i++)
    ihc[i+1] = c1LP*ihcl[i+1] + c2LP*(ihc[i]+ihcl[i]);
  for(j=0; j<=order;j++) ihcl[j] = ihc[j];
  return(ihc[order]);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the Control path using Nonlinear Function after OHC */

double NLafterohc(double x,double taumin, double taumax, double asym)
{    
	double R,dc,R1,s0,x1,out,minR;

	minR = 0.05;
    R  = taumin/taumax;
    
	if(R<minR) minR = 0.5*R;
    else       minR = minR;
    
    dc = (asym-1)/(asym+1.0)/2.0-minR;
    R1 = R-minR;

    /* This is for new nonlinearity */
    s0 = -dc/log(R1/(1-minR));
	
    x1  = fabs(x);
    out = taumax*(minR+(1.0-minR)*exp(-x1/s0));
	if (out<taumin) out = taumin; 
    if (out>taumax) out = taumax;
    return(out);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the IHC Nonlinear Function (Logarithmic Transduction Functions) */

double NLogarithm(double x, double slope, double asym, double cf)
{
	double corner,strength,xx,splx,asym_t;
	    
    corner    = 80; 
    strength  = 20.0e6/pow(10,corner/20);
            
    xx = log(1.0+strength*fabs(x))*slope;
    
    if(x<0)
	{
		splx   = 20*log10(-x/20e-6);
		asym_t = asym -(asym-1)/(1+exp(splx/5.0));
		xx = -1/asym_t*xx;
	};   
    return(xx);
}
/* -------------------------------------------------------------------------------------------- */

/* divide */
COMPLEX compdiv(COMPLEX ne,COMPLEX de)
{
  double d;
  COMPLEX z;
  d=de.x*de.x+de.y*de.y;
  z.x=(ne.x*de.x+ne.y*de.y)/d;
  z.y=(ne.y*de.x-ne.x*de.y)/d;
  return(z);
}
/* this returns a complex number equal to exp(i*theta) */
COMPLEX compexp(double theta)
{
  COMPLEX dummy;
  dummy.x = cos(theta);
  dummy.y = sin(theta);
  return dummy;
}
/* Multiply a complex number by a scalar */
COMPLEX compmult(double scalar, COMPLEX compnum)
{
 COMPLEX answer;
 answer.x = scalar * compnum.x;
 answer.y = scalar * compnum.y;
 return answer;
}
/* Find the product of 2 complex numbers */
COMPLEX compprod(COMPLEX compnum1, COMPLEX compnum2)
{
 COMPLEX answer;
 answer.x = (compnum1.x * compnum2.x) - (compnum1.y * compnum2.y);
 answer.y = (compnum1.x * compnum2.y) + (compnum1.y * compnum2.x);
 return answer;
}
/* add 2 complex numbers */
COMPLEX comp2sum(COMPLEX summand1, COMPLEX summand2)
{
 COMPLEX answer;
 answer.x = summand1.x + summand2.x;
 answer.y = summand1.y + summand2.y;
 return answer;
}
/* add three complex numbers */
COMPLEX comp3sum(COMPLEX summand1, COMPLEX summand2, COMPLEX summand3)
{
 COMPLEX answer;
 answer.x = summand1.x + summand2.x + summand3.x;
 answer.y = summand1.y + summand2.y + summand3.y;
 return answer;
}

/* subtraction: complexA - complexB */
COMPLEX compsubtract(COMPLEX complexA, COMPLEX complexB)
{
 COMPLEX answer;
 answer.x = complexA.x - complexB.x;
 answer.y = complexA.y - complexB.y;
 return answer;
}
/* Get the real part of the complex */
double REAL(COMPLEX compnum)
{ return(compnum.x); } 

/* Get the imaginary part of the complex */
double IMAG(COMPLEX compnum)
{ return(compnum.y); } 

/* Get the conjugate of the complex signal */
COMPLEX compconj(COMPLEX complexA)
{
  COMPLEX answer;
  answer.x = complexA.x;
  answer.y = -complexA.y;
  return (answer);
}