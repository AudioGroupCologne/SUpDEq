/* Author: Tobias May (2009)

/*
% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 
*/

#include <math.h>
#include "mex.h"
#include "string.h"

/* Helper functions */
#define max(x, y)   ((x) > (y) ? (x) : (y))
#define	min(A, B)	((A) < (B) ? (A) : (B))
#define getRound(x) ((x) >= 0?(long)((x)+0.5):(long)((x)-0.5))

/* Input Arguments */
#define   INPUT      		 prhs[0]  // input signal [nSamples x nChannels]
#define   BLOCKSIZE		     prhs[1]  // blocksize
#define   HOPSIZE   		 prhs[2]  // hopsize
#define   WINDOW    		 prhs[3]  // window
#define   NFRAMES            prhs[4]  // number of frames

/* Output Arguments */
#define   OUTPUT			 plhs[0]  // framed signal matrix


void usage()
{
	mexPrintf("\n");		
	mexPrintf(" frameData   Frame input signal and apply analysis window. \n");
	mexPrintf("\n");		
	mexPrintf(" USAGE \n");
	mexPrintf(" frames = frameData(input,blockSize,hopSize,window,nFrames)\n");                                     
	mexPrintf("\n");                                                                      
	mexPrintf(" INPUT ARGUMENTS\n");                                                                
	mexPrintf("         input : mono input signal [nSamples x 1]\n");            
	mexPrintf("     blockSize : block size in samples\n");
    mexPrintf("       hopSize : step size in samples \n");
	mexPrintf("        window : analysis window of dimension [blockSize x 1]\n");
	mexPrintf("       nFrames : Number of frames \n");
	mexPrintf("\n"); 
	mexPrintf(" OUTPUT ARGUMENTS\n");                                                                
	mexPrintf("        frames : framed signal [blockSize x nFrames] \n");            
	mexPrintf("\n"); 
	mexPrintf(" Implementation by Tobias May © 2008, Philips Research Eindhoven \n");
	mexPrintf("\n"); 
}    


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *input, *output, *window;
	int blockSize, hopSize, overlap, nFrames;
	int hh, ii, nSamples, exit = false;
	int inIdx, outIdx;

  	/* Check for proper number of arguments */
  	if ((nrhs < 5) || (nlhs > 1) || (int) mxGetN(INPUT) > 1) {
    	usage();
    	exit = true;   
  	} 
  		
  	if(!exit){
  		/* Check the dimensions and type of INPUT */
	  	nSamples  = (int) mxGetM(INPUT);

	  	if (!mxIsNumeric(INPUT) || mxIsComplex(INPUT) || mxIsSparse(INPUT)  || !mxIsDouble(INPUT))
   			mexErrMsgTxt("frameData.dll requires a real input vector.");
		
 		// Assign pointers and values to the various parameters
		input     = mxGetPr(INPUT);
		blockSize = (int) mxGetScalar(BLOCKSIZE);
		hopSize   = (int) mxGetScalar(HOPSIZE);
		nFrames   = (int) mxGetScalar(NFRAMES);
		
		overlap   = blockSize - hopSize;
		window    = mxGetPr(WINDOW);

		// Determine dimensionality of output matrix
		// int dim[] = {blockSize, nFrames, nChannels};

		// OUTPUT = mxCreateNumericArray(3, dim, mxDOUBLE_CLASS, mxREAL);

    	// Create a matrix for the return arguments
		OUTPUT = mxCreateDoubleMatrix(blockSize, nFrames, mxREAL);
  		output = mxGetPr(OUTPUT);

        // Loop over number of frames
		for (hh=0; hh < nFrames; hh++){

			// Channel offset
			// inIdx  = 1 + (hh * hopSize) + nSamples;
			inIdx  = 1 + (hh * hopSize);
			// Frame offset
			// outIdx = (hh*blockSize)+(gg*blockSize*nFrames);
			outIdx = (hh*blockSize);
			
			// Loop over frame size
			for (ii=0; ii<blockSize; ii++){
				// Apply analysis window
				output[ii+outIdx] = input[ii+inIdx-1] * window[ii];
			}
		}
	}
}   /* end mexFunction() */
