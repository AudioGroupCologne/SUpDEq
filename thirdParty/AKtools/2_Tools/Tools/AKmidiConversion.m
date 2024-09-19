% [f0, midiNoteNumber, noteString] = AKmidiConversion(value, type, fA4)
% converts between MIDI note numbers, note strings, and fundamental
% frequencies f0.
%
% Examples:
% AKmidiConversion( 69,   'midiNoteNumber' )
% AKmidiConversion( 440,  'f0' )
% AKmidiConversion( 'A4', 'noteString' )
% all calls return f0=440 midiNoteNumber=69, and noteStr='A4'
%
%
% I N P U T
% value - numeric value or string
% type  - 'f0', 'midiNoteNumber', or 'noteString' to specify what value was
%         passed to AKmidiConversion
% fA4   - tuning frequency in Hz at A4 (default = 440)
%
% O U T P U T
% see examples
%
% 1/2018 - fabian.brinkmann@tu-berlin.de

% AKtools
% Copyright (C) 2016 Audio Communication Group, Technical University Berlin
% Licensed under the EUPL, Version 1.1 or as soon they will be approved by
% the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: 
% http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software 
% distributed under the License is distributed on an "AS IS" basis, 
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
% See the License for the specific language governing  permissions and
% limitations under the License. 
function [f0, midiNoteNumber, noteString] = AKmidiConversion(value, type, fA4)

% --- default parameter for tuning frequency
if nargin == 2
    fA4 = 440;
end

% --- get the fundamental frequency
if any(strcmpi(type, {'midiNoteNumber', 'm'})) || ( any(strcmpi(type, {'noteString', 'n'})) && ischar(value) )
    f0 = midi2freq(value, fA4);
else
    f0  = value;
end

% --- get the midiNote and noteString
[midiNoteNumber, noteString] = freq2midi(f0, fA4);

end

% ------------------------------------- sub-functions (self explainatory :)
function [midiNoteNumber, noteString] = freq2midi(f0, fA4)
    if (nargin == 1)
        fA4  = 440;
    end

    midiNoteNumber = round(69 + 12*log2(f0/fA4));
    noteString = midi2str( round(midiNoteNumber) );
end

function f0 = midi2freq(note, fA4)
    if (nargin == 1)
        fA4  = 440;
    end

    if ischar(note)
        midiNoteNumber = str2midi(note);
    else
        midiNoteNumber = note;
    end

    f0 = 2^( ( midiNoteNumber+12*log2(fA4)-69 ) / 12 );
end

function noteString = midi2str(midiNoteNumber)
    notes      = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
    noteString = [notes{mod(midiNoteNumber, 12)+1} num2str(floor(midiNoteNumber/12)-1)];
end

function midiNoteNumber = str2midi(noteString)

    notes = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};

    noteString = upper(noteString);

    idSharp =  strfind(noteString, '#');

    if idSharp
        noteOctave = str2double( noteString( idSharp+1:end ) );
        note       = noteString(1:idSharp);
    else
        noteOctave = str2double( noteString( 2:end ) );
        note       = noteString(1);
    end

    [~, noteID] = ismember(note, notes);

    midiNoteNumber = (noteOctave+1)*12 + noteID-1;
end






