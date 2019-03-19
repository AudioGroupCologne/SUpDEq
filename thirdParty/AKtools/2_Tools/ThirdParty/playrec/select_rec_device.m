function[deviceID, dev] = select_rec_device()
%SELECT_REC_DEVICE provides a prompt to select which device to record from
%using Playrec.  Returns the device ID followed by a struct containing
%all available information on the device.
%
%Robert Humphrey, January 2008

devs = playrec('getDevices');
validIDs = -1;

prompt = '\nAvailable input devices:\n -1) No Device\n';

for k=1:length(devs)
    if(devs(k).inputChans)
        prompt = [prompt, sprintf(' %2d) %s (%s) %d channels\n', ...
            devs(k).deviceID, devs(k).name, ...
            devs(k).hostAPI, devs(k).inputChans)];
        validIDs = [validIDs, devs(k).deviceID];
    end
end

fprintf([prompt, '\n']);

if is_octave
    fflush(stdout);
end

deviceID = input('Select which device to use [default -1]: ', 's');

while(~isempty(deviceID) && ~any(validIDs == str2double(deviceID)))
    deviceID = input('Invalid choice, please select which device to use [default -1]: ', 's');
end

if isempty(deviceID)
    deviceID = -1;
    dev = [];
else
    deviceID = str2double(deviceID);
    dev = devs([devs.deviceID] == deviceID);
end
