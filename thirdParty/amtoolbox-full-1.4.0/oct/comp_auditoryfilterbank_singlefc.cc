/*
 *=========================================================================
 * An efficient C implementation of the 4th order gammatone filter
 *-------------------------------------------------------------------------
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *-------------------------------------------------------------------------
%
%  [bm, env, instp, instf] = gammatone_c(x, fs, cf, hrect) 
%
%  x     - input signal
%  fs    - sampling frequency (Hz)
%  cf    - centre frequency of the filter (Hz)
%  hrect - half-wave rectifying if hrect = 1 (default 0)
%
%  bm    - basilar membrane displacement
%  env   - instantaneous envelope
%  instp - instantaneous phase (unwrapped radian)
%  instf - instantaneous frequency (Hz)
%
%
%  The gammatone filter is commonly used in models of the auditory system.
%  The algorithm is based on Martin Cooke's Ph.D work (Cooke, 1993) using 
%  the base-band impulse invariant transformation. This implementation is 
%  highly efficient in that a mathematical rearrangement is used to 
%  significantly reduce the cost of computing complex exponentials. For 
%  more detail on this implementation see
%  http://www.dcs.shef.ac.uk/~ning/resources/gammatone/
%
%  Once compiled in Matlab this C function can be used as a standard 
%  Matlab function:
%  >> mex gammatone_c.c
%  >> bm = gammatone_c(x, 16000, 200);
%
%  Ning Ma, University of Sheffield
%  n.ma@dcs.shef.ac.uk, 09 Mar 2006
% 
 * CHANGES:
 * 2012-05-30 Ning Ma <n.ma@dcs.shef.ac.uk>
 *   Fixed a typo in the implementation (a5). The typo does not make a lot
 *   of difference to the response. Thanks to Vijay Parsa for reporting
 *   the problem.
 *
 * 2010-02-01 Ning Ma <n.ma@dcs.shef.ac.uk>
 *   Clip very small filter coefficients to zero in order to prevent
 *   gradual underflow. Arithmetic operations may become very slow with
 *   subnormal numbers (those smaller than the minimum positive normal
 *   value, 2.225e-308 in double precision). This could happen if the 
 *   input signal cotains many zeros (e.g. impulse responses). Thanks to
 *   John Culling for reporting the problem.
 *=========================================================================
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"


/*=======================
 * Input arguments
 *=======================
 */
#define IN_x        prhs[0]
#define IN_fs       prhs[1]
#define IN_cf       prhs[2]
#define IN_hrect    prhs[3]

/*=======================
 * Output arguments
 *=======================
 */
#define OUT_bm      plhs[0] /* Basilar membrane response */
#define OUT_env     plhs[1] /* Instantaneous Hilbert envelope */
#define OUT_instp   plhs[2] /* Instantaneous phase */
#define OUT_instf   plhs[3] /* Instantaneous frequency */

/*=======================
 * Useful Const
 *=======================
 */
#define BW_CORRECTION      1.0190
#define VERY_SMALL_NUMBER  1e-200
#ifndef M_PI
#define M_PI               3.14159265358979323846
#endif

/*=======================
 * Utility functions
 *=======================
 */
#define myMax(x,y)     ( ( x ) > ( y ) ? ( x ) : ( y ) )
#define myMod(x,y)     ( ( x ) - ( y ) * floor ( ( x ) / ( y ) ) )
#define erb(x)         ( 24.7 * ( 4.37e-3 * ( x ) + 1.0 ) )

/*=======================
 * Main Function
 *=======================
 */
void mexFunction ( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
   double *x, *bm, *env, *instp, *instf;
   int i, j, t, fs, nsamples, hrect;
   double cf, a, tpt, tptbw, gain;
   double p0r, p1r, p2r, p3r, p4r, p0i, p1i, p2i, p3i, p4i;
   double a1, a2, a3, a4, a5, u0r, u0i; /*, u1r, u1i;*/
   double qcos, qsin, oldcs, coscf, sincf, oldphase, dp, dps;

  /*=========================================
   * Check argument numbers
   *=========================================
   */
   if ( nrhs < 3 ) {
      mexPrintf ( "??? Not enough input arguments.\n" );
      return;
   }
   if ( nrhs > 4 ) {
      mexPrintf ( "??? Too many input arguments.\n" );
      return;
   }
   if ( nlhs > 4 ) {
      mexPrintf ( "??? Too many output arguments.\n" );
      return;
   }
  /*=========================================
   * input arguments
   *=========================================
   */
   if ( nrhs < 4 ) {
      hrect = 0;
   }
   else {
      hrect = ( int ) mxGetScalar ( IN_hrect );
   }
   
   x = mxGetPr ( IN_x );   /* waveform */
   i = mxGetN ( IN_x ); 
   j = mxGetM ( IN_x );
   if ( i>1 && j>1 ) {
      mexPrintf ( "??? Input x must be a vector.\n" );
      return;
   }
   nsamples = myMax ( i, j );
   fs = ( int ) mxGetScalar ( IN_fs );   /* sampling rate */
   cf = mxGetScalar ( IN_cf );   /* centre frequency */

  /*=========================================
   * output arguments
   *=========================================
   */
   OUT_bm = mxCreateDoubleMatrix ( 1, nsamples, mxREAL );
   bm = mxGetPr ( OUT_bm );
   if ( nlhs > 1 ) {
      OUT_env = mxCreateDoubleMatrix ( 1, nsamples, mxREAL );
      env = mxGetPr ( OUT_env );
   }
   if ( nlhs > 2 ) {
      OUT_instp = mxCreateDoubleMatrix ( 1, nsamples, mxREAL );
      instp = mxGetPr ( OUT_instp );
   }
   if ( nlhs > 3 ) {
      OUT_instf = mxCreateDoubleMatrix ( 1, nsamples, mxREAL );
      instf = mxGetPr ( OUT_instf );
   }

  /*=========================================
   * Initialising variables
   *=========================================
   */
   oldphase = 0.0;
   tpt = ( M_PI + M_PI ) / fs;
   tptbw = tpt * erb ( cf ) * BW_CORRECTION;
   a = exp ( -tptbw );

  /* based on integral of impulse response */
   gain = ( tptbw*tptbw*tptbw*tptbw ) / 3;

  /* Update filter coefficients */
   a1 = 4.0*a; a2 = -6.0*a*a; a3 = 4.0*a*a*a; a4 = -a*a*a*a; a5 = a*a;
   p0r = 0.0; p1r = 0.0; p2r = 0.0; p3r = 0.0; p4r = 0.0;
   p0i = 0.0; p1i = 0.0; p2i = 0.0; p3i = 0.0; p4i = 0.0;
 
  /*===========================================================
   * exp(a+i*b) = exp(a)*(cos(b)+i*sin(b))
   * q = exp(-i*tpt*cf*t) = cos(tpt*cf*t) + i*(-sin(tpt*cf*t))
   * qcos = cos(tpt*cf*t)
   * qsin = -sin(tpt*cf*t)
   *===========================================================
   */
   coscf = cos ( tpt * cf );
   sincf = sin ( tpt * cf );
   qcos = 1; qsin = 0;   /* t=0 & q = exp(-i*tpt*t*cf)*/
   for ( t=0; t<nsamples; t++ )
   {
     /* Filter part 1 & shift down to d.c. */
      p0r = qcos*x[t] + a1*p1r + a2*p2r + a3*p3r + a4*p4r;
      p0i = qsin*x[t] + a1*p1i + a2*p2i + a3*p3i + a4*p4i;

     /* Clip coefficients to stop them from becoming too close to zero */
      if (fabs(p0r) < VERY_SMALL_NUMBER)
        p0r = 0.0F;
      if (fabs(p0i) < VERY_SMALL_NUMBER)
        p0i = 0.0F;
      
     /* Filter part 2 */
      u0r = p0r + a1*p1r + a5*p2r;
      u0i = p0i + a1*p1i + a5*p2i;

     /* Update filter results */
      p4r = p3r; p3r = p2r; p2r = p1r; p1r = p0r;
      p4i = p3i; p3i = p2i; p2i = p1i; p1i = p0i;
  
     /*==========================================
      * Basilar membrane response
      * 1/ shift up in frequency first: (u0r+i*u0i) * exp(i*tpt*cf*t) = (u0r+i*u0i) * (qcos + i*(-qsin))
      * 2/ take the real part only: bm = real(exp(j*wcf*kT).*u) * gain;
      *==========================================
      */
      bm[t] = ( u0r * qcos + u0i * qsin ) * gain;
      if ( 1 == hrect && bm[t] < 0 ) {
         bm[t] = 0;  /* half-wave rectifying */
      }
     /*==========================================
      * Instantaneous Hilbert envelope
      * env = abs(u) * gain;
      *==========================================
      */
      if ( nlhs > 1 ) {
         env[t] = sqrt ( u0r * u0r + u0i * u0i ) * gain;
      }
     /*==========================================
      * Instantaneous phase
      * instp = unwrap(angle(u));
      *==========================================
      */
      if (nlhs > 2) {
         instp[t] = atan2 ( u0i, u0r );
         /* unwrap it */
         dp = instp[t] - oldphase;
         if ( abs ( dp ) > M_PI ) {
            dps = myMod ( dp + M_PI, 2 * M_PI) - M_PI;
            if ( dps == -M_PI && dp > 0 ) {
               dps = M_PI;
            }
            instp[t] = instp[t] + dps - dp;
         }
         oldphase = instp[t];
      }
     /*==========================================
      * Instantaneous frequency
      * instf = cf + [diff(instp) 0]./tpt;
      *==========================================
      */
      if ( nlhs > 3 && t > 0 ) {
         instf[t-1] = cf + ( instp[t] - instp[t-1] ) / tpt;
      }

     /*====================================================
      * The basic idea of saving computational load:
      * cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
      * sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
      * qcos = cos(tpt*cf*t) = cos(tpt*cf + tpt*cf*(t-1))
      * qsin = -sin(tpt*cf*t) = -sin(tpt*cf + tpt*cf*(t-1))
      *====================================================
      */
      qcos = coscf * ( oldcs = qcos ) + sincf * qsin;
      qsin = coscf * qsin - sincf * oldcs;
   }
   if ( nlhs > 3 ) {
      instf[nsamples-1] = cf;
   }
   return;
}
  
/* end */
