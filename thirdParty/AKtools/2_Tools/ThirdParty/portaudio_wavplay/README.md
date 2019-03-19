licence PortAudio -> pa_wavplay
===============================
PortAudio Portable Real-Time Audio Library 
Copyright (c) 1999-2011 Ross Bencina and Phil Burk

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


pa-wavplay
==========

Play and record multi-channel audio in Matlab.

This project is a modification of one previously developed by Matt Frear and posted 
on both Source Forge and the Mathworks File Exchange.

The original code:

- played multi-channel audion from Matlab using ASIO, DirectSound, and Windows Audio.
- worked with 32-bit installations of Matlab only on Windows.

This update extends the original to:

- play multi-channel audio from Matlab using ASIO, WASAPI, DirectSound, and Windows Audio.
- works with both 32-bit and 64-bit installations of Matlab on Windows.  

The code is built upon the open source PortAudio API.

Requirements: Windows Vista or later.
