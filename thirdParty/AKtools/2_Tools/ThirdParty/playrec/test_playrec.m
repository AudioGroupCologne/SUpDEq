function[] = test_playrec()
%TEST_PLAYREC runs through a set of tests to check if things seem to be
%working at the same time as providing an example how to use the utility.
%
%Robert Humphrey, January 2008

samplerates = [192000 96000 48000 44100 32000 16000 8000];

if exist('playrec', 'file') ~= 3
    error ('The Playrec MEX file does not yet exist in the current directory or on the search path');
end

buildDefines = regexpi(playrec('about'), 'Built with defines:\s*(.*?),?\s*(?:\n|$)', 'tokens');
deviceSummary = regexpi(playrec('about'), 'Available host API:\s*(.*?),?\s*(?:\n|$)', 'tokens');

if ~isempty(buildDefines)
    print_flush('Playrec was built with the following defines: %s\n', cell2mat(buildDefines{1}));    
end

if ~isempty(deviceSummary)
    print_flush('Playrec was built with the following host API: %s\n', cell2mat(deviceSummary{1}));    
end

if isempty(playrec('getDevices'))
    error ('There are no devices available using the selected host APIs.');
end

[playDevId, playDev] = select_play_device;

if playDevId == -1
    print_flush ('\nSkipping output tests as no device selected\n');
else
    print_flush ('\nStarting output tests using device %d (%s, %s):\n', ...
        playDevId, playDev.name, playDev.hostAPI);
    
    % Check if the device has already been initialised.  If it has then
    % this would mean it cannot be initialised with the required
    % configuration so it needs to be reset.  There are a set of get*
    % commands (eg getPlayDevice) to find out more about the current
    % configuration if necessary.
    if playrec('isInitialised')
        print_flush ('   Resetting playrec as previously initialised\n');
        playrec('reset');
    end
    
    %Initialise using a range of sample rates until one that works is
    %found.
    for rate = samplerates
        try
            playrec('init', rate, playDevId, -1);
            print_flush ('   Initialising device at %dHz succeeded\n', rate);
            break;
        catch
            print_flush ('   Initialising device at %dHz failed with error: %s\n', rate, lasterr);
        end
    end
    
    if ~playrec('isInitialised')
        warning ('Failed to initialise device at any sample rate');
    else
        %Is correctly initialised so lets do some testing!
        %Firstly play a 1s 1khz sine wave on each channel in turn
        test_wav = 0.8 * sin((1:rate)*2*pi*1000/rate)';
        
        print_flush('\n   Testing simple blocking output\n');
        print_flush('   A 1s long 1kHz tone should be heard on each output in turn.\n');
        
        %This goes through adding the samples for one channel only after
        %all previous pages (ie the previous channel) has finished
        for chan=1:playDev.outputChans
            print_flush('      Adding output on channel %d\n', chan);
            playrec('play', test_wav, chan);
            
            %This could use the block command instead
            while(playrec('isFinished') == 0)
            end
        end
        print_flush('      Output complete\n'); 
        
        pause(3);
        
        print_flush('\n   Testing simple non-blocking output\n');
        print_flush('   A 1s long 1kHz tone should be heard on each output in turn.\n');
        
        %This goes through adding the samples for one channel after
        %the other without waiting for completion.
        for chan=1:playDev.outputChans
            print_flush('      Adding output on channel %d\n', chan);
            playrec('play', test_wav, chan);
        end

        print_flush('      All samples added, waiting for output to complete\n');
        %This could use the block command instead
        while(playrec('isFinished') == 0)
        end
        print_flush('      Output complete\n'); 

        pause(3);
        
        print_flush('\n   Testing continuous non-blocking output on all channels, using pages 0.1s long\n');
        print_flush('   A different 5s long tone should be heard simultaneously on all outputs, with no glitches.\n');

        page_list = [];
        is_first = 1;
        % run for 5s
        for loopcount=0:49
            %Generate the next page of samples to output.  This is occuring
            %whilst the previous pages are still being output.
            sample_no = loopcount * rate/10 + (1:rate/10)';
            chan_freq = 500 + 100 * (1:playDev.outputChans);
            out_wav = 0.8 * sin(sample_no*chan_freq*2*pi/rate);
            print_flush('      Adding samples %d to %d\n', min(sample_no), max(sample_no));
            
            %Add the new page number to the end of the list
            page_list = [page_list, playrec('play', out_wav, 1:playDev.outputChans)];
            
            %If this is the first time through then reset the skipped
            %sample count
            if is_first
                is_first = 0;
                playrec('resetSkippedSampleCount');
            end
            
            %If there are more than 3 pages waiting, wait for the earliest
            %to end before continuing - any number of pages buffering
            %can be used, it is a trade off between the latency and the 
            %chance of audio glitches.  The number of samples in each page
            %also plays a part in deciding the best number to use.
            if(length(page_list) > 3)
                while(playrec('isFinished', page_list(1)) == 0)
                end 
                
                %Remove finished page from list
                page_list = page_list(2:end);
            end
        end

        print_flush('      All samples added, waiting for output to complete\n');
        
        %Although there are still samples to output, because they have been
        %supplied to Playrec glitches cannot occur between them so the 
        %number of samples skipped can be checked.       
        print_flush('      Glitches lasting %d samples occured\n', ...
            playrec('getSkippedSampleCount'));

        %This could use the block command instead
        while(playrec('isFinished') == 0)
        end
        print_flush('      Output complete\n'); 
                
    end
end

%
%recDev = select_rec_device;
%
%if recDev == -1
%    display ('\nSkipping record tests as no device selected');
%else
%    display ('\nStarting record tests...');
%    display ('\Record tests finished');
%end
%

function print_flush(varargin)

fprintf(varargin{:});
if is_octave
    fflush(stdout);
end

