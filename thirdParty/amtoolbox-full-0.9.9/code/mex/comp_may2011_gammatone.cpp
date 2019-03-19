// Gammatone filtering ...
// Note that this implementaion is based on the source code from the ohio-state university
// http://www.cse.ohio-state.edu/pnl/shareware/roman-jasa03/


#include <stdio.h>
#include "math.h"
#include "mex.h"

/* Input Arguments */
#define   INPUT      			 prhs[0]  // input signal 
#define   FS                     prhs[1]  // Sampling rate of input signal in hertz
#define   LOWERFREQ   			 prhs[2]  // 
#define   UPPERFREQ              prhs[3]  // 
#define   NFILTER                prhs[4]  // 
#define   USEMIDDLEEAR      	 prhs[5]  // 
									      
#define   TIMEALIGN              prhs[6]  // Time-align gammaton channels
#define   INFOFLAG               prhs[7]  // Print gammatone parameter on screen

/* Output Arguments */
#define   BM     				 plhs[0]  // Basilar membrane displacement
#define   ENV     				 plhs[1]  // Envelope

/* Helper functions */
#define max(x, y)   ((x) > (y) ? (x) : (y))
#define	min(A, B)	((A) < (B) ? (A) : (B))
#define swap(A,B)   temp = (A); (A)=(B); (B) = temp;
#define getRound(x) ((x) >= 0?(long)((x)+0.5):(long)((x)-0.5))

#define PI					  (3.14159265358979323846)
// ERB bandwidth correction 4th order
#define bandwidthCorrection   1.019  

#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names))


/* Outer/middle ear */
#define MIDDLE_EAR_SIZE 29
#define DB 60.0

double f[MIDDLE_EAR_SIZE];
double af[MIDDLE_EAR_SIZE];
double bf[MIDDLE_EAR_SIZE];
double tf[MIDDLE_EAR_SIZE];


void usage()
{
	mexPrintf("\ngammatone");		
	mexPrintf("\n MEX-File implementation of an auditory filterbank.");
	mexPrintf("\n");
	mexPrintf("\n USAGE: ");
	mexPrintf("\n\t[BM,ENV] = gammatone(IN,FS,LOWFREQ,UPFREQ,NFILTER,BEAR,BALIGN,BINFO);");
	mexPrintf("\n");
	mexPrintf("\n INPUT PARAMETER:  ");
	mexPrintf("\n\t      IN - mono input signal [nSamples x 1]");
	mexPrintf("\n\t      FS - sampling frequency in Hz");
	mexPrintf("\n\t LOWFREQ - center frequency of lowest auditory filter ");
	mexPrintf("\n\t           (default, lowFreq = 80)");
	mexPrintf("\n\t  UPFREQ - center frequency of highest auditory filter");
	mexPrintf("\n\t           (default, upFreq = 5000)");
	mexPrintf("\n\t NFILTER - number of auditory filters which will be spaced linear in");
	mexPrintf("\n\t           the ERB domain. If NFILTER is not a scalar but a vector, the");
	mexPrintf("\n\t           first value is assumed to represent the number of auditory");
	mexPrintf("\n\t           channels and the following values represents the indices of");
	mexPrintf("\n\t           the filters which should be processed. This can be useful if");
	mexPrintf("\n\t           a large number of channels is required as MATLAB might run");
	mexPrintf("\n\t           out of memory for a longer exerpt if all filters should be");
	mexPrintf("\n\t           computed in one step.");
	mexPrintf("\n\t           (default, NFILTER = round(21.4*log10((FS/2)*0.00437 + 1.0)))");
	mexPrintf("\n\t    BEAR - adjust gain coefficients of the auditory channels to");
	mexPrintf("\n\t           incorporate middle ear effects. Note that this feature can");
	mexPrintf("\n\t           only be used if the center frequencies of the auditory");
	mexPrintf("\n\t           channels are above 20 hertz and below 12500 hertz.");
	mexPrintf("\n\t  BALIGN - phase-aligned gammatone output (produce non-causal output)");
	mexPrintf("\n\t           (default, BALIGN = false)");
	mexPrintf("\n\t   BINFO - info flag printing gammatone parameters on the screen");
	mexPrintf("\n\t           (default, BINFO = false)");
	mexPrintf("\n");
	mexPrintf("\n OUTPUT PARAMETER:  ");
	mexPrintf("\n\t      BM - basilar membrane displacement [nSamples x nFilter]");
	mexPrintf("\n\t     ENV - instantaneous envelope        [nSamples x nFilter]");
	mexPrintf("\n");
	mexPrintf("\n NOTE: ");
	mexPrintf("\n\tThis implementaion is based on the source code from the ohio-state");
	mexPrintf("\n\tuniversity: www.cse.ohio-state.edu/pnl/shareware/roman-jasa03/");
	mexPrintf("\n");	
	mexPrintf("\n\tImplementation by Tobias May, 2007-08");
	mexPrintf("\n\tPlease send bug reports to: tobias.may@philips.com\n");
	mexPrintf("\n"); 
}

struct gammaTone
{
	double cf, bw;
	double midEarCoeff;
	double delay;
	double phaseCorrection;
	double delta;
	int    align;
	double p[4];
	double q[4];
};

double erb(double x)
{
	// ERB   Equivalent rectangular bandwidth
	return 24.7*(0.00437 * x + 1.00);
}


double DBtoAmplitude(double dB)
{
  return pow(10.0,(dB/20.0));
}

double loudnessLevelInPhons(double dB, double freq)
/*
Uses linear interpolation of the look-up tables to compute the loudness level,
in phons, of a pure tone of frequency freq using the reference curve for sound
 pressure level dB.
The equation is taken from section 4 of BS3383.
*/
{
  int i=0;
  double afy, bfy, tfy;

  if ((freq<20.0) | (freq>12500.0)) {
    mexPrintf("Can't compute a outer/middle ear gain for that frequency\n");
    exit(0);
  }
  while (f[i] < freq) i++;
  afy=af[i-1]+(freq-f[i-1])*(af[i]-af[i-1])/(f[i]-f[i-1]);
  bfy=bf[i-1]+(freq-f[i-1])*(bf[i]-bf[i-1])/(f[i]-f[i-1]);
  tfy=tf[i-1]+(freq-f[i-1])*(tf[i]-tf[i-1])/(f[i]-f[i-1]);
  return 4.2+afy*(dB-tfy)/(1.0+bfy*(dB-tfy));
}

double HzToERBRate(double Hz)
{
	return( 21.4*log10(Hz*0.00437 + 1.0) );
}

double ERBRateToHz(double ERBRate)
{
	return( (pow(10, ERBRate/21.4) - 1) / 0.00437 );
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *input, *bm, *env, *filter;
	double lowerERB, upperERB;
	double dt, twoPiT, gain, z;
	double f1Down, f2Down, f1Up, f2Up, data;
	double x[4], y[4];
	int    nSamples, nChannels, lowerFreq, upperFreq, nFilter;
 	int    ii, nn, chan, fs, bExit = false, bUseMidEar = false;
	int    chanIdx, nFilter2Process, *processFilter, infoFlag;
	int    bTimeAlign = false, tc;

	// Define field names for gammatone structure
	//const char *field_names[] = {"fs","nFilter","lowerFreq","upperFreq","cf","bw","gain","delay"};
	

	// Check for proper number of arguments
  	if ((nrhs < 2) || (nlhs > 2)){
		usage();
		// Set exit flag to true because of incorrect function call
		bExit = true;   
	}
	// Check if input arguments are valid
	else{
		// Get dimension of input signal
		nSamples  = (int)mxGetM(INPUT);
	  	nChannels = (int)mxGetN(INPUT);

		// Sampling frequency
		fs = (int) mxGetScalar(FS);

		// Input signal must be mono
		if(nChannels != 1){
			mexPrintf("\n--------------------------------------------------------------------------------"  );
		    mexPrintf("\n ERROR: gammatone.dll is just supporting mono signals."										    );
		    mexPrintf("\n--------------------------------------------------------------------------------\n");
		    // Set exit flag to true
		    bExit = true;
		}

		// Set defaults ...
		if (nrhs < 3)
			lowerFreq = 80;
		else
			lowerFreq = (int)mxGetScalar(LOWERFREQ);

		if (nrhs < 4)
			upperFreq = 5000;
		else
			upperFreq = (int)mxGetScalar(UPPERFREQ);

		if (nrhs < 5){
			// Round the number of auditory filter to the next integer value
			nFilter = getRound(HzToERBRate(fs/2));
			// Total number of filters to process within this function call
			nFilter2Process = nFilter;
			// Create index vector for auditory filter processing
			processFilter = (int *) mxCalloc(nFilter2Process, sizeof(int));
			// Create vector which indicates the indicees of the processed auditory filters 
			for (chan=0; chan<nFilter; chan++){
				processFilter[chan] = chan;
			}
		}
		else{
			// Figure out how may auditory filters should be processed ...
			filter = mxGetPr(NFILTER);
			// Assume the first value to be the total number of auditory filters
			nFilter = (int) filter[0];
			// Determine number of auditory filters which should be processed
			nFilter2Process = max((int)mxGetM(NFILTER),(int)mxGetN(NFILTER))-1;

			// In case just the total number of auditory filters was specified, process all filters
			if (nFilter2Process == 0){
				nFilter2Process = nFilter;
				// Create index vector for auditory filter processing
				processFilter = (int *) mxCalloc(nFilter2Process, sizeof(int));

				// Create vector which indicates the indicees of the processed auditory filters 
				for (chan=0; chan<nFilter2Process; chan++){
					processFilter[chan] = chan;
				}
			}
			else{
				// Create index vector for auditory filter processing
				processFilter = (int *) mxCalloc(nFilter2Process, sizeof(int));
				// Create vector which indicates the indicees of the processed auditory filters 
				for (chan=0; chan<nFilter2Process; chan++){
					processFilter[chan] = (int) filter[chan + 1] - 1;
				}
			}
		}

		if (nrhs < 6)
			bUseMidEar = false;
		else{
			bUseMidEar = (int)mxGetScalar(USEMIDDLEEAR);
		}

		if (nrhs < 7)
			bTimeAlign = false;
		else{
			bTimeAlign = (int)mxGetScalar(TIMEALIGN);
		}


		// Check if lower and upper frequency limit fall within the valid range
		// of the outer/middle ear filter ...
		if (bUseMidEar && (lowerFreq <= 20  || upperFreq >= 12500 )){
			
			// Display error message
			mexPrintf("\n--------------------------------------------------------------------------------"  );
			mexPrintf("\n ERROR using gammatone.dll:                                                     "  );
			mexPrintf("\n The lower and upper frequency limit must be within 20 and 12500 Hz if the      "  );
			mexPrintf("\n outer/middle ear gain is used.                                                 "  );
			mexPrintf("\n--------------------------------------------------------------------------------\n");
		
			// Set exit flag to true because of incorrect function call
			bExit = true;  
		}
		else if(bUseMidEar){

	   		/* parameters of equal-loudness functions from BS3383,
			"Normal equal-loudness level contours for pure tones under 
			free-field listening conditions", table 1. f is the tone frequency
			af and bf are frequency-dependent coefficients tf is the threshold 
			sound pressure level of the tone, in dBs */
		
			f[0]=20.0;     af[0]=2.347;  bf[0]=0.00561;   tf[0]=74.3;
			f[1]=25.0;     af[1]=2.190;  bf[1]=0.00527;   tf[1]=65.0;
			f[2]=31.5;     af[2]=2.050;  bf[2]=0.00481;   tf[2]=56.3;
			f[3]=40.0;     af[3]=1.879;  bf[3]=0.00404;   tf[3]=48.4;
			f[4]=50.0;     af[4]=1.724;  bf[4]=0.00383;   tf[4]=41.7;
			f[5]=63.0;     af[5]=1.579;  bf[5]=0.00286;   tf[5]=35.5;
			f[6]=80.0;     af[6]=1.512;  bf[6]=0.00259;   tf[6]=29.8;
			f[7]=100.0;    af[7]=1.466;  bf[7]=0.00257;   tf[7]=25.1;
			f[8]=125.0;    af[8]=1.426;  bf[8]=0.00256;   tf[8]=20.7;
			f[9]=160.0;    af[9]=1.394;  bf[9]=0.00255;   tf[9]=16.8;
			f[10]=200.0;   af[10]=1.372; bf[10]=0.00254;  tf[10]=13.8;
			f[11]=250.0;   af[11]=1.344; bf[11]=0.00248;  tf[11]=11.2;
			f[12]=315.0;   af[12]=1.304; bf[12]=0.00229;  tf[12]=8.9;
			f[13]=400.0;   af[13]=1.256; bf[13]=0.00201;  tf[13]=7.2;
			f[14]=500.0;   af[14]=1.203; bf[14]=0.00162;  tf[14]=6.0;
			f[15]=630.0;   af[15]=1.135; bf[15]=0.00111;  tf[15]=5.0;
			f[16]=800.0;   af[16]=1.062; bf[16]=0.00052;  tf[16]=4.4;
			f[17]=1000.0;  af[17]=1.000; bf[17]=0.00000;  tf[17]=4.2;
			f[18]=1250.0;  af[18]=0.967; bf[18]=-0.00039; tf[18]=3.7;
			f[19]=1600.0;  af[19]=0.943; bf[19]=-0.00067; tf[19]=2.6;
			f[20]=2000.0;  af[20]=0.932; bf[20]=-0.00092; tf[20]=1.0;
			f[21]=2500.0;  af[21]=0.933; bf[21]=-0.00105; tf[21]=-1.2;
			f[22]=3150.0;  af[22]=0.937; bf[22]=-0.00104; tf[22]=-3.6;
			f[23]=4000.0;  af[23]=0.952; bf[23]=-0.00088; tf[23]=-3.9;
			f[24]=5000.0;  af[24]=0.974; bf[24]=-0.00055; tf[24]=-1.1;
			f[25]=6300.0;  af[25]=1.027; bf[25]=0.00000;  tf[25]=6.6;
			f[26]=8000.0;  af[26]=1.135; bf[26]=0.00089;  tf[26]=15.3;
			f[27]=10000.0; af[27]=1.266; bf[27]=0.00211;  tf[27]=16.4;
			f[28]=12500.0; af[28]=1.501; bf[28]=0.00488;  tf[28]=11.6;
		
		}

		if (nrhs < 8)
			infoFlag = false;
		else
			infoFlag = (int)mxGetScalar(INFOFLAG);
	}

	if(!bExit){
		
		// Asign pointers
		input = mxGetPr(INPUT);
	
		// Create a matrix for the return argument
		BM = mxCreateDoubleMatrix(nSamples, nFilter2Process, mxREAL);
		// Asign pointer
  		bm = mxGetPr(BM);

		// Envelope output
		if (nlhs>1){
			// Create a matrix for second return argument
			ENV = mxCreateDoubleMatrix(nSamples, nFilter2Process, mxREAL);
			// Asign pointer
			env = mxGetPr(ENV);
		}


		// Initialize gammatone filter structure 
		gammaTone* fChan;
		// Da der Kompiler meckert wenn nFilter keine Konstante ist, muss die
		// Gammatone-Struktur zur Laufzeit mit new initialisiert werden ...
		fChan = new gammaTone[nFilter];


		// Map frequency to erb domain
		lowerERB = HzToERBRate(lowerFreq);
		upperERB = HzToERBRate(upperFreq);
  
		// Calculate center frequencies
		fChan[0].cf = ERBRateToHz(lowerERB);

		for (chan = 1; chan<nFilter-1; chan++)
		{	
			// Linear space the gammatone center frequency in the ERB domain
			fChan[chan].cf = ERBRateToHz(lowerERB + chan * (upperERB-lowerERB)/((double)nFilter-1.0));
		}
		fChan[nFilter-1].cf = ERBRateToHz(upperERB);

		// Report parameter
		if (infoFlag){
			mexPrintf("\n \t |-------------------------------------| ");
			mexPrintf("\n \t |          GAMMATONE PARAMETER        | ");
			mexPrintf("\n \t |-------------------------------------| ");
			mexPrintf("\n \t | filter |  cf [Hz] |   bw   | delay  | ");
			mexPrintf("\n \t |-------------------------------------| ");
		}

		// Loop over number of auditory filters
		for (chan=0; chan<nFilter; chan++)
		{
			// Bandwidth
			fChan[chan].bw    = erb(fChan[chan].cf) * bandwidthCorrection;
			// Compute gammatone envelope delay in samples (for 4th order)
			fChan[chan].delay = (3.0 * fs) / (fChan[chan].bw * 2.0 * PI);
			// Compute phase correction to align peak in temporal fine structure
			fChan[chan].phaseCorrection = -fChan[chan].cf * 3.0/fChan[chan].bw;

			// Report parameters ...
			if (infoFlag){
				mexPrintf("\n \t |%5.0f   |%9.2f |%7.2f |%7.2f |",(double)chan+1.0,fChan[chan].cf,fChan[chan].bw,fChan[chan].delay);
			}
    	}

		if (infoFlag){
			mexPrintf("\n \t |-------------------------------------|\n");
		}


    	// Time resolution
		dt     = 1/double(fs);
		twoPiT = 2 * PI * dt;
			
		// Loop over number of auditory filters
		for (chan = 0; chan < nFilter2Process; chan++)
		{
			// Get index of current filter
			chanIdx = processFilter[chan];

			// Calculate gain depending on middle ear filter
			if (bUseMidEar){
				// Calculate middle ear coefficient
				fChan[chanIdx].midEarCoeff = DBtoAmplitude(loudnessLevelInPhons(DB,fChan[chanIdx].cf)-DB);
				// Calculate overall channel gain
				gain = fChan[chanIdx].midEarCoeff * pow(twoPiT * fChan[chanIdx].bw, 4) / 3.0;
			}
			else{
				gain = pow(twoPiT * fChan[chanIdx].bw, 4) / 3.0;
			}

			z = exp(-twoPiT * fChan[chanIdx].bw);


			if (bTimeAlign){
				// Integrate phase correction 
				f1Down = cos(fChan[chanIdx].cf * twoPiT + fChan[chanIdx].phaseCorrection) * z;
				f2Down = sin(fChan[chanIdx].cf * twoPiT + fChan[chanIdx].phaseCorrection) * z;

				// Delay
				tc = getRound(fChan[chan].delay);
			}
			else{
				f1Down = cos(fChan[chanIdx].cf * twoPiT) * z;
				f2Down = sin(fChan[chanIdx].cf * twoPiT) * z;

				// Delay
				tc = 0;
			}

			// Without phase correction
			f1Up = cos(fChan[chanIdx].cf * twoPiT) * z;
			f2Up = sin(fChan[chanIdx].cf * twoPiT) * z;

			// Initialize gammatone filter states
			for (ii = 0; ii < 4; ii++)
			{
				fChan[chanIdx].p[ii] = 0;
				fChan[chanIdx].q[ii] = 0;
			}

				// Loop over length of input signal
				for (nn = 0; nn < (nSamples + tc); nn++)
				{
					// *==========================================
					// * Basilar membrane displacement
					// *==========================================
					if (nn >= tc){
						// % Basilar membrane displacement
						bm[(nn-tc)+chan*nSamples] = fChan[chanIdx].p[3] * gain;
					}

					// *==========================================
					// * Envelope
					// *==========================================
					if (nlhs>1){
						if (nn >= tc){
							// % Basilar membrane displacement
							env[(nn-tc)+chan*nSamples] = sqrt(fChan[chanIdx].p[3] * fChan[chanIdx].p[3] + fChan[chanIdx].q[3] * fChan[chanIdx].q[3]) * gain;
						}
					}



					// Loop over the order of the gammatone filterbank
					for (ii = 0; ii < 4; ii++)
					{
						x[ii] = f1Up * fChan[chanIdx].p[ii] - f2Up * fChan[chanIdx].q[ii];
						y[ii] = f2Up * fChan[chanIdx].p[ii] + f1Up * fChan[chanIdx].q[ii];
					}

					// Zero-padding
					data = (nn < nSamples) ? input[nn] : 0.0;

					fChan[chanIdx].p[0] = data * f1Down + x[0];
					fChan[chanIdx].q[0] = data * f2Down + y[0];

					fChan[chanIdx].p[1] = fChan[chanIdx].p[0] + x[1];
					fChan[chanIdx].q[1] = fChan[chanIdx].q[0] + y[1];

					fChan[chanIdx].p[2] = fChan[chanIdx].p[1] + x[1] + x[2];
					fChan[chanIdx].q[2] = fChan[chanIdx].q[1] + y[1] + y[2];
		
					fChan[chanIdx].p[3] = fChan[chanIdx].p[2] + x[1] + 2*x[2] + x[3];
					fChan[chanIdx].q[3] = fChan[chanIdx].q[2] + y[1] + 2*y[2] + y[3];
				}
		}

		// Free the memory used for the gammatone structure
		delete fChan;

		} // end of if(!bExit)
		

}   /* end mexFunction() */
