%% SUpDEq - Spatial Upsampling by Directional Equalization
%
% function win = supdeq_win(signalLength, windowLength)
%
% This function generates a single- or both-sided hann window which can be
% used for windowing in time domain
%
% Output:
% win           - Window in time domain
%
% Input:
% signalLength  - Length of signal to be windowed in samples 
% windowLength  - Values of head and tail window length in samples [a,b].
%                 If [a,0] is passed, a single sided (head-)window will be
%                 designed with head-length a. If [0,b] is passed, a single 
%                 sided (tail-)window will be designed with tail-length b.
%                 If only a single value is passed, a head window will be
%                 designed.
%
% Dependencies: -
%
% (C) 2018 by JMA, Johannes M. Arend
%             TH Köln - University of Applied Sciences
%             Institute of Communications Engineering
%             Department of Acoustics and Audio Signal Processing

function win = supdeq_win(signalLength, windowLength)

if length(windowLength) < 2
    winMode = 'headWin';
    headWinLength = windowLength(1);
end

if length(windowLength) >= 2
    if windowLength(1) == 0
        winMode = 'tailWin';
        tailWinLength = windowLength(2);
    elseif windowLength(2) == 0
        winMode = 'headWin';
        headWinLength = windowLength(1);
    else
        winMode = 'bothSidedWin';
        headWinLength = windowLength(1);
        tailWinLength = windowLength(2);
    end
end

%% Design window
if strcmp(winMode,'headWin')
   win_head = 0.5*(1 - cos(2*pi*(0:headWinLength*2-1)'/(headWinLength*2-1)));
   win = [win_head(1:headWinLength,:);ones(signalLength-headWinLength,1)];
end

if strcmp(winMode,'tailWin')
   win_tail = 0.5*(1 - cos(2*pi*(0:tailWinLength*2-1)'/(tailWinLength*2-1)));
   win = [ones(signalLength-tailWinLength,1);win_tail(tailWinLength+1:end,:)];
end

if strcmp(winMode,'bothSidedWin')
    win_head = 0.5*(1 - cos(2*pi*(0:headWinLength*2-1)'/(headWinLength*2-1)));
    win_tail = 0.5*(1 - cos(2*pi*(0:tailWinLength*2-1)'/(tailWinLength*2-1)));
    win = [win_head(1:headWinLength,:);ones(signalLength-headWinLength-tailWinLength,1);win_tail(tailWinLength+1:end,:)];
end

end

