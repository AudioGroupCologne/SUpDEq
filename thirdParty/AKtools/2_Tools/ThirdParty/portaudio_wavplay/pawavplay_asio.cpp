/*
pa_wavplay
  A matlab function for playing an audio buffer on multi-channel hardware.
  Matt Frear, MARCS Auditory Lab, Sydney, Australia
  m.frear@uws.edu.au

  with help with the recording code from:
  Paul Henderson, Rensselaer Polytechnic Institute, USA, Department of Architectural Acoustics
  hendep2@rpi.edu

  Modified July 2014 by Joseph Desloge, Sensimetrics Corp., Malden, MA USA
  desloge@sens.com
  
  mexFunction is the entry point from Matlab - parses arguments and converts input

  pa_wavplay opens the portaudio stream

  paWavCallback is the callback function for portaudio, plays the audio
 
 CHANGES:
 1.0 separated ASIO, DirectSound, and Win Audio into separate files.  Added WASAPI capability.
 0.21 User can specify range of channels to record from.
 0.20 Added recording and simultaneous recording and playback functionality
      Now uses floats internally instead of int16s
	  Much of the command parsing done in external .m files now.
 0.1 initial release
*/

#define API_NAME "ASIO"
//#define API_NAME "MME"
//#define API_NAME "Windows DirectSound"
//#define API_NAME "Windows WASAPI"

#define BLOCK 2048
//#define BLOCK 8192

#include <math.h>
#include <string.h>

#include "mex.h"
#include "portaudio.h"
#include "pa_win_wasapi.h"


/////////// my data types
typedef float SAMPLE; // format for portaudio - float = 32 bit

enum recordmode // play, or record, or both
{
	play,
	record,
	playrecord
};

// My struct for holding the audiobuffer and info about it. Received by the paWavCallback
typedef struct
{
	SAMPLE *buffer;      // PLAYBACK audio buffer //pdh
	SAMPLE *recbuffer;   // RECORDING audio buffer //pdh
	int bufpos;          // current play pos in the buffer
	int buflen;          // total length of the buffer
	int bufchannels;     // number of audio channels
	
	int recbuffirstchan; // first RECORDING audio channel
	int recbufflastchan; // last recording channel
	
	int recbufpos;       // current rec pos in the buffer
	int recbuflen;       // total length of the RECORDING buffer

	recordmode recmode;  // what we're doing - playing, recording, or both
}
paWavData;

/////////// function prototypes
static int paWavCallback( void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void *userData );
int pa_wavplayrec(SAMPLE *buffer, int buflen, int bufchannels, SAMPLE *recbuffer, int recbuflen, int firstrecchannel, int lastrecchannel, int deviceID, int recdevID, int samplerate);
void pa_printdevinfo();
void printdevice(int api_idx, int device);

// Convert functions - to convert an input buffer to our SAMPLE buffer
SAMPLE* convDouble(double *buffer, int buflen);
SAMPLE* convFloat(float *buffer, int buflen);
SAMPLE* convUchar(unsigned char *buffer, int buflen);

// Entry point from Matlab
// Parse arguments, convert buffer to my SAMPLE (float), then call pa_wavplay with parsed args
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	bool recording = true;
	bool playing = true;

	// Check for proper number of arguments 
/*	char *usage = "Usage: recordbuffer = pa_wavplay(playbuffer, playdevice, samplerate, recfirstchannel, reclastchannel, recnsamples, recdevice)\n pa_wavplayrecord with no arguments will list all the available devices\n";

	if (nrhs != 0 && nrhs != 6) 
		mexErrMsgTxt(usage);
*/
	if (nrhs == 0)
	{
//		mexPrintf(usage);
		pa_printdevinfo();
		return;
	}

	
	////// get deviceids. If deviceid < 0 don't play. if recdeviceid < 0 don't record ///////
	// get playback device
	int playdeviceid = 0;
	double* devptr = mxGetPr(prhs[1]);
	if (devptr != NULL)
		playdeviceid = (int)*devptr;

	if (playdeviceid < 0)
		playing = false;

	// get record device
	int recdeviceid = 0;
	devptr = mxGetPr(prhs[6]);
	if (devptr != NULL)
		recdeviceid = (int)*devptr;

	if (recdeviceid < 0)
		recording = false;

	if (!recording && !playing)
	{
		mexPrintf("playdevice < 0 and recdevice < 0. Nothing to be done");
//		mexErrMsgTxt(usage);
	}

	////// Next step: convert input buffer to SAMPLE buffer /////////
	if (mxIsComplex(prhs[0]))
		mexErrMsgTxt("audiobuffer must be noncomplex.");

	int bufrows = mxGetM(prhs[0]); // get number of rows in the buffer
	int bufcols = mxGetN(prhs[0]); // get number of columns
	
	SAMPLE *myBuf = NULL;

	if (playing)
	{
		// Float32 - sweet, same as SAMPLE, no conversion necessary
		if (mxIsSingle(prhs[0]))
		{
			SAMPLE *bufptr = (SAMPLE *)mxGetData(prhs[0]);

			if (bufptr == NULL)
				mexErrMsgTxt("audiobuffer is NULL");

			myBuf = bufptr;
		}

		// double
		else if (mxIsDouble(prhs[0]))
		{
			double *bufptr = mxGetPr(prhs[0]);

			if (bufptr == NULL)
				mexErrMsgTxt("audiobuffer is NULL");

			myBuf = convDouble(bufptr, bufrows*bufcols);
		}
		
		else
		{
			mexErrMsgTxt("audiobuffer is an invalid data type");
		}
	}

	// get samplerate
	double* sampleptr = mxGetPr(prhs[2]);

	if (sampleptr == NULL)
		mexErrMsgTxt("samplerate is NULL");

	double samplerate = sampleptr[0];

	int recordsamples = 0, firstrecordchannel = 0, lastrecordchannel = 0, Inputbufrows = 0, Inputbufcols = 0;

	SAMPLE *myInputBuf = NULL;

	if (recording)
	{
		// get number of record channels
		double* sptr = mxGetPr(prhs[3]);

		if (sptr == NULL)
			mexErrMsgTxt("firstrecordchannel is NULL");

		firstrecordchannel = (int)sptr[0];
		if (firstrecordchannel <=0)
			mexErrMsgTxt("invalid firstrecordchannel: <= 0");

		sptr = mxGetPr(prhs[4]);

		if (sptr == NULL)
			mexErrMsgTxt("lastrecchannel is NULL");

		lastrecordchannel = (int)sptr[0];
		if (lastrecordchannel <=0)
			mexErrMsgTxt("invalid lastrecordchannel <= 0");

		if (lastrecordchannel < firstrecordchannel)
			mexErrMsgTxt("invalid lastrecordchannel, < firstrecordchannel");

		Inputbufcols = lastrecordchannel - firstrecordchannel + 1;

		// get number of record samples
		sptr = mxGetPr(prhs[5]);

		if (sptr == NULL)
			mexErrMsgTxt("recsamples is NULL");

		recordsamples = (int)sptr[0];
	
		if (recordsamples <= 0)
			Inputbufrows = bufrows; // get number of rows in the buffer
		else
			Inputbufrows = recordsamples;

		///// allocate memory for the output matrix
		plhs[0] = mxCreateNumericMatrix(Inputbufrows,Inputbufcols, mxSINGLE_CLASS, mxREAL);
		SAMPLE* inbuff = (SAMPLE*)mxGetData(plhs[0]);
		myInputBuf = inbuff;
	}
	
	if (!playing)
	{
		playdeviceid = paNoDevice;
		bufrows = bufcols = 0;
	}
	
	if (!recording)
	{
		recdeviceid = paNoDevice;
		Inputbufcols = Inputbufrows = 0;
	}

	pa_wavplayrec(myBuf, bufrows*bufcols, bufcols, myInputBuf, firstrecordchannel, lastrecordchannel, Inputbufrows*Inputbufcols, playdeviceid, recdeviceid, (int)samplerate);

	if (!mxIsSingle(prhs[0]))
		delete myBuf; // delete myBuf if necessary

	return;
}

SAMPLE* convDouble(double *buffer, int buflen)
{
	// Convert to SAMPLEs for playback
	SAMPLE *myBuf = new SAMPLE[buflen]; // allocate buffer

	for (int i = 0; i < buflen; i++)
		myBuf[i] = (SAMPLE)(buffer[i]); 

	return myBuf;
}

int pa_wavplayrec(SAMPLE *buffer, int buflen, int bufchannels, SAMPLE *recbuffer, int recbuffirstchannel, int recbufflastchannel, int recbuflen, int deviceID, int recdevID, int samplerate)
{
	//mexPrintf("buflen=%i bufchannels=%i samplerate=%i\n", buflen, bufchannels, samplerate);

	// Make our wavData object, to pass to OpenStream. This gets given to us again in the
	// callback. I guess you could just have a global variable instead...
	paWavData wav;
	wav.buffer = buffer;
	wav.recbuffer = recbuffer;
	wav.buflen = buflen;
	wav.recbuflen = recbuflen;
	wav.bufchannels = bufchannels;
	wav.recbuffirstchan = recbuffirstchannel;
	wav.recbufflastchan = recbufflastchannel;
	wav.bufpos = 0;
	wav.recbufpos = 0;
	
	if (recbuffer == NULL)
		if (buffer == NULL)
			mexErrMsgTxt("recbuffer && buffer NULL in pa_wavplay");
		else
		{
			mexPrintf("Playing on device %i\n", deviceID);
			wav.recmode = play;
		}
	else
		if (buffer == NULL)
		{
			mexPrintf("Recording on device %i\n", recdevID);
			wav.recmode = record;
		}
		else
		{
			mexPrintf("Recording on device %i\n", recdevID);
			mexPrintf("Playing on device %i\n", deviceID);

			wav.recmode = playrecord;
		}

	PaStream *stream;
    PaError err;

    err = Pa_Initialize();
    if( err != paNoError ) 
		goto error;

	int dev_type_indx = -1;
	const PaHostApiInfo* Api_Info;
	int num_APIs = (int )Pa_GetHostApiCount();
	for (int i=0; i<num_APIs; i++) {
 		Api_Info = Pa_GetHostApiInfo(i);
		if (strcmp(Api_Info->name,API_NAME)==0)
			dev_type_indx = i;
	}
	if (dev_type_indx == -1) {
		goto error;
	}
	Api_Info = Pa_GetHostApiInfo(dev_type_indx);
	if (recdevID != -1)
		recdevID = Pa_HostApiDeviceIndexToDeviceIndex(dev_type_indx, recdevID);
	if (deviceID != -1)
		deviceID = Pa_HostApiDeviceIndexToDeviceIndex(dev_type_indx, deviceID);
	
	// Open an audio I/O stream. 
	PaWasapiStreamInfo api_info;
	api_info.size = sizeof(PaWasapiStreamInfo);
	api_info.version = 1;
	api_info.hostApiType = paWASAPI;
	api_info.flags = paWinWasapiExclusive;
	void *p_api_info;
	if ((strcmp(API_NAME, "Windows WASAPI")== 0)) {
		p_api_info = &api_info;
	} else {
		p_api_info = NULL;
	}
    PaStreamParameters inputParameters;		// Parameters governing portaudio input stream
	PaStreamParameters outputParameters;	// Parameters governing portaudio output stream
    inputParameters.channelCount = recbufflastchannel;
    inputParameters.device = recdevID;
    inputParameters.hostApiSpecificStreamInfo = p_api_info;
    inputParameters.sampleFormat = paFloat32;	// Input data returned as floating point
	if (recdevID != -1)
	    inputParameters.suggestedLatency = Pa_GetDeviceInfo(recdevID)->defaultLowInputLatency ;
	outputParameters.channelCount = bufchannels;
    outputParameters.device = deviceID;
    outputParameters.hostApiSpecificStreamInfo = p_api_info;
    outputParameters.sampleFormat = paFloat32;	// Output data returned as floating point
	if (deviceID != -1)
	    outputParameters.suggestedLatency = Pa_GetDeviceInfo(deviceID)->defaultLowInputLatency ;

	err = Pa_IsFormatSupported(NULL, &outputParameters, (double )samplerate);

	if (recbuffer == NULL)
		err = Pa_OpenStream(
			  &stream,
			  NULL,
			  &outputParameters,
			  samplerate,     // stream sample rate
			  BLOCK,            // frames per buffer 
			  paNoFlag,       // stream flag
			  (PaStreamCallback (__cdecl *))paWavCallback,  // callback function
			  &wav);          // user data
	else
		if (buffer == NULL)
			err = Pa_OpenStream(
				  &stream,
				  &inputParameters,
				  NULL,
				  samplerate,     // stream sample rate
				  BLOCK,            // frames per buffer 
				  paNoFlag,       // stream flag
				  (PaStreamCallback (__cdecl *))paWavCallback,  // callback function
				  &wav);          // user data
		else
			err = Pa_OpenStream(
				  &stream,
				  &inputParameters,
				  &outputParameters,
				  samplerate,     // stream sample rate
				  BLOCK,            // frames per buffer 
				  paNoFlag,       // stream flag
				  (PaStreamCallback (__cdecl *))paWavCallback,  // callback function
				  &wav);          // user data


    if( err != paNoError ) 
		goto error;

    err = Pa_StartStream( stream );
    
	if( err != paNoError ) 
		goto error;
    
	// sleep while stream's active
	while( Pa_IsStreamActive(stream) )
		Pa_Sleep( 50 );

	if (err = Pa_IsStreamStopped( stream ) == 0) {
	    err = Pa_StopStream( stream );
	    if( err != paNoError ) 
			goto error;
	} else if (err != paNoError )
		goto error;

	err = Pa_CloseStream( stream );
    
	if( err != paNoError ) 
		goto error;
    
	Pa_Terminate();

	return err;

	// wtf? A goto? Yeah, that's leftover from the portaudio example code
error:
    Pa_Terminate();
    mexPrintf( "An error occured while using the portaudio stream\n" );
    mexPrintf( "Error number: %d\n", err );
    mexPrintf( "Error message: %s\n", Pa_GetErrorText( err ) );
    return err;
	
}

// This routine will be called by the PortAudio engine when audio is needed.
// All we need to do is take the audio in my userData and put it into the output buffer.
// 
// The paWavData bufpos holds where we're up to in the audio buffer.
// Remember that we see the channels as stored one after the other, not interleaved.
// However, the portaudio outputbuffer needs each channel interleaved
// So we use the offset variable to calc where each channel starts in the buffer.
static int paWavCallback( void *inputBuffer, void *outputBuffer, 
	unsigned long framesPerBuffer, 
	const PaStreamCallbackTimeInfo* timeInfo, 
	PaStreamCallbackFlags statusFlags, void *userData )
{
    // Cast data passed through stream to my paWavData
    paWavData *data = (paWavData*)userData;

    SAMPLE *out = (SAMPLE *)outputBuffer;
    SAMPLE *in = (SAMPLE *)inputBuffer;

	//(void) outTime; // Prevent unused variable warnings.

	int offset = 0;
	if (data->bufchannels > 0)
		offset = data->buflen / data->bufchannels;
	
	int recchannelsize = data->recbuflen / (data->recbufflastchan - (data->recbuffirstchan - 1));
	

	// for each frame
    for(unsigned int i=0; i<framesPerBuffer; i++ )
    {
		if (data->recmode == play || data->recmode == playrecord)
		{
			// for each output channel
			for (int c = 0; c < data->bufchannels; c++)
			{
				// if we're not at the end of the buffer
				if (data->bufpos < (data->buflen / data->bufchannels))
				{
					SAMPLE *valptr = (data->buffer + data->bufpos);
					valptr = valptr + (c * offset); // add offset for this channel

					SAMPLE val = *valptr;

					*out = val;
					out++;
				}
				else
					*out++ = 0;
			}

			data->bufpos++;	// advance to next pos
		}
	
		// Hmm. Possible bug hereabouts.
		// had to add && in != null for when playing and recording simultaneously. Not sure why.
		if ((data->recmode == record || data->recmode == playrecord) && in != NULL)
		{
			// for each input channel
			int c = 0;
			for (int i = 0; i < data->recbufflastchan; i++)
			{
				// if we're not at the end of the buffer
				if (data->recbufpos < recchannelsize)
				{
					// set nothing for first unused channels, but increment in pointer
					if (i >= data->recbuffirstchan - 1) 
					{
						SAMPLE *valptr = (data->recbuffer + data->recbufpos);
						valptr = valptr + (c++ * recchannelsize); // add offset for this channel

						*valptr = *in;
					}

					in++;
				}
			}

			data->recbufpos++;	// advance to next pos
		}
    }
	
	// if we're past the end of the buffer
	if (data->recmode == play)
		if (data->bufpos >= (data->buflen / data->bufchannels))
			return 1; // return of non-zero signals stream end
	
	if (data->recmode == playrecord)
		// I used to change recmode from playrecord to record if there was nothing left to play
		// but this would emit an audible buzz as portaudio needs *out++ = 0;
		// if (data->bufpos >= (data->buflen / data->bufchannels))
		//	data->recmode = record; // signal no more to play
		if (data->recbufpos >= recchannelsize)
		{
			if (data->bufpos >= (data->buflen / data->bufchannels)) // if no more to play
				return 1;
			else
				data->recmode = play; // no more to record, play only
		}

	if (data->recmode == record)
		if (data->recbufpos >= recchannelsize)
			return 1; // return of non-zero signals stream end

	return 0;
}

// prints the device to stdout. Don't call this function, call pa_printdevinfo
// as Pa_Init needs to be done first.
void printdevice(int api_idx, int device)
{
	const PaDeviceInfo *pdi;
	pdi = Pa_GetDeviceInfo(Pa_HostApiDeviceIndexToDeviceIndex(api_idx, device));

		
	if (pdi != NULL)
	{
		char *format = NULL;
		
/*		switch (pdi->nativeSampleFormats)
		{
		case (0):
			format = strdup("float32");	break;
		case (1):
			format = strdup("int16"); break;
		case (2):
			format = strdup("int32"); break;
		case (3):
			format = strdup("int24"); break;
		case (4):
			format = strdup("packedint24");	break;
		case (5):
			format = strdup("int8"); break;
		case (6):
			format = strdup("uint8"); break;
		case (16):
			format = strdup("custom"); break;
		}

		mexPrintf("\nDevice %i\nDevice name: %s\nMax Input channels: %i  Max Output channels: %i\nNative Sample Format: %s\n"
			, device, pdi->name, pdi->maxInputChannels, pdi->maxOutputChannels, format);
*/
		mexPrintf("\nDevice %i\nDevice name: %s\nMax Input channels: %i  Max Output channels: %i\n"
			, device, pdi->name, pdi->maxInputChannels, pdi->maxOutputChannels);
		if (format != NULL)
			delete format;
	}
}

// print info about the audio device to stdout
void pa_printdevinfo()
{
	mexPrintf("Printing audio devices...\n");

	Pa_Initialize();		

	int dev_type_indx = -1;
	const PaHostApiInfo* Api_Info;
	int num_APIs = (int )Pa_GetHostApiCount();
	for (int i=0; i<num_APIs; i++) {
 		Api_Info = Pa_GetHostApiInfo(i);
		if (strcmp(Api_Info->name,API_NAME)==0)
			dev_type_indx = i;
	}
	if (dev_type_indx == -1) {
		return;
	}
	Api_Info = Pa_GetHostApiInfo(dev_type_indx);


	for( int d = 0; d < Api_Info->deviceCount; d++) {
		//PaDeviceIndex dev_indx = Pa_HostApiDeviceIndexToDeviceIndex(dev_type_indx, d);
		//printdevice((int )dev_indx);
		printdevice(dev_type_indx, d);
	}
	
	Pa_Terminate();
}