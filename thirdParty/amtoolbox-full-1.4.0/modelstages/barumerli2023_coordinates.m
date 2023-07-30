%BARUMERLI2023_COORDINATES Class useful to handl HRTF coordinates with different conventions
%   Usage: coordinates = barumerli2023_coordinates(SOFAobj);
%
%   Interface:
% 
%       barumerli2023_coordinates: constructor. The coordinates are converted into the 
%                    cartestian system is used.
%
%       return_positions: return the coordinates given a specific convention.
%
%       convert_positions: convert the stored coordinates given
%                    a specific convention.
%
%       find_position: return a position given the index of the coordinates
%                    matrix
%
%       normalize_distance: normalize the distance between receiver and
%                    source
% 
% 
%   Purpose:
%   Manage easily the coordinates system of the SOFA object.
% 
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%
%   Examples:
%   ---------
% 
%   coordinates = barumerli2023_coordinates(SOFAobj);
%   horpolar = coords.return_positions('horizontal-polar');
%   
%   Load coordinates from SOFA object and convert them into the
%   horizontal-polar system. The output is organized as in the SOFA object.
%
%   Url: http://amtoolbox.org/amt-1.4.0/doc/modelstages/barumerli2023_coordinates.php


%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: MATLAB SOFA M-STATISTICS M-Control M-Signal
%   #Author: Roberto Barumerli (2022)

% This file is licensed unter the GNU General Public License (GPL) either 
% version 3 of the license, or any later version as published by the Free Software 
% Foundation. Details of the GPLv3 can be found in the AMT directory "licences" and 
% at <https://www.gnu.org/licenses/gpl-3.0.html>. 
% You can redistribute this file and/or modify it under the terms of the GPLv3. 
% This file is distributed without any warranty; without even the implied warranty 
% of merchantability or fitness for a particular purpose. 

classdef barumerli2023_coordinates < handle

    properties (SetAccess = private)
        pos double% {mustBeReal, mustBeFinite}
        pos_type char% {mustBeMember(pos_type,{'horizontal-polar','spherical','cartesian', 'geodesic'})} = 'cartesian'
    end
    methods
        % constructor
        function obj = barumerli2023_coordinates(data, convention)
            if nargin == 1
                if strcmp(data.GLOBAL_Conventions, 'SOFA')
                    obj.pos = SOFAcalculateAPV(data);
                    obj.pos_type = data.SourcePosition_Type;
                else
                    error('barumerli2023_coordinates: SOFA object not valid!');
                end
            elseif nargin == 2
                if ismatrix(data)
                    assert(size(data, 2) == 3, 'Check matrix dimensions')
                    assert(ismember(convention,{'horizontal-polar','spherical','cartesian', 'geodesic'}), 'Check convention type')
                    obj.pos = data;
                    obj.pos_type = convention;
                else
                    error('barumerli2023_coordinates: position matrix not valid!');
                end
            else
                obj.pos = zeros(0,3);
                obj.pos_type = 'cartesian';
            end
            
            obj.convert_positions('cartesian');
        end
        
        % return coordinates with a specific 
        function r = return_positions(obj, pos_type)
            if strcmp(pos_type, obj.pos_type)
                r = obj.pos;
            else
                r = wrapperSOFAconvert(obj.pos, obj.pos_type, pos_type);
            end
        end
        
        function concatenate(obj, obj_new)
            obj.pos = vertcat(obj.pos, obj_new.return_positions(obj.pos_type));
        end
        
        function convert_positions(obj, pos_type)
            obj.pos = wrapperSOFAconvert(obj.pos,obj.pos_type, pos_type);
            obj.pos_type = pos_type;
        end
        
        function r = count_pos(obj)
            r = size(obj.pos, 1);
        end
        
        function r = find_position(obj, idx, pos_type)
            r = wrapperSOFAconvert(obj.pos(idx,:),obj.pos_type, pos_type);
        end
        
        function normalize_distance(obj)
            pos_type_temp = obj.pos_type;
            obj.convert_positions('spherical');
            obj.pos(:,3) = 1;
            obj.convert_positions(pos_type_temp);
        end
        
        function [idx, coords_new] = find_positions(obj, coords_search)
            [idx, coords_new] = local_SOFAfind(obj, coords_search);
        end    
        
        function plot(obj)
            r = return_positions(obj, 'cartesian');
            figure
            scatter3(r(:,1),r(:,2),r(:,3),20,0.5*ones(size(r, 1),1),'filled');
            view([1 0 0])
            axis equal;   
        end

    end
    
	methods (Access = private)
    end
end

function pos_new = wrapperSOFAconvert(pos, pos_type, pos_type_new)
    pos_new = local_SOFAconvertCoordinates(pos, pos_type, pos_type_new);
    
%     pos_inv = local_SOFAconvertCoordinates(pos_new, pos_type_new, pos_type);
%     assert(sum(abs(pos-pos_inv),'all')<1e-10)
   % assert(sum(abs(imag(pos_new)),'all')==0)
    assert(sum(abs(imag(pos_new(:)))) == 0)
    if strcmp(pos_type_new, 'horizontal-polar')
        c = sum(max(pos_new(abs(pos_new(:,1))<5,2)) > 270);
        if  c > 0
          warning('%i polar angles have been found to be greater than 270.',c)
        end
    end
end


function output = local_SOFAconvertCoordinates(input,input_type,output_type)
    %adapted from SOFAconvertCoordinates
   

    %% check input
    if strcmp(input_type,'cartesian')==0 && ...
            strcmp(input_type,'spherical')==0 && ...
            strcmp(input_type,'geodesic')==0 && ...
            strcmp(input_type,'horizontal-polar')==0
        error('Specified "input_type" is not supported');
    end
    if strcmp(output_type,'cartesian')==0 && ...
            strcmp(output_type,'spherical')==0 && ...
            strcmp(output_type,'geodesic')==0 && ...
            strcmp(output_type,'horizontal-polar')==0
        error('Specified "output_type" is not supported');
    end

    output=input;
    %% convert coordinates if necessary
    if strcmp(output_type,input_type)==0
        temp=input;
        switch input_type
            case 'cartesian'
                %do nothing
            case {'spherical','geodesic'}
                [temp(:,1),temp(:,2),temp(:,3)]=sph2cart(deg2rad(input(:,1)),deg2rad(input(:,2)),input(:,3));
            case 'horizontal-polar'
                [az, el] = hor2sph(input(:,1), input(:,2));
                [temp(:,1),temp(:,2),temp(:,3)]=sph2cart(deg2rad(az),deg2rad(el),input(:,3));
        end

        output=temp;
        switch output_type
            case 'cartesian'
                %do nothing
            case {'spherical','geodesic'}
                [output(:,1),output(:,2),output(:,3)]=cart2sph(temp(:,1),temp(:,2),temp(:,3));
                output(:,1:2)=rad2deg(output(:,1:2));
            case 'horizontal-polar'
                [output(:,1),output(:,2),output(:,3)]=cart2sph(temp(:,1),temp(:,2),temp(:,3));
                output(:,1:2)=rad2deg(output(:,1:2));
                [output(:, 1), output(:, 2)] = sph2hor(output(:,1), output(:,2));
        end
    end
end

function [idx, coords_new] = local_SOFAfind(coords, coords_search)
% local adaptation of SOFAfind
% [idx, azi, ele, r] = SOFAfind(Obj,azi,ele,r) finds the indecies to 
% the HRTFs from OBJ according to the trajectory given in AZI, ELE, R.
% Input: 
%		Obj: SOFA object containing HRTFs
%		azi, ele: direction (in degrees) for azimuth and elevation
%       r: optional radius. If not provided, radius will be ignored.
% 
% Output: 
%		idx: index of the filters (corresponds to AZI and ELE)
%		azi, ele: azimuth and elevation of the actual position (degrees)
%       r: actual radius

    
    %% create a 2D-grid with nearest positions
    pos = coords.return_positions('cartesian');
    pos_seach = coords_search.return_positions('cartesian');
    
    idx=zeros(size(pos_seach,1),1);
    
    for ii=1:size(pos_seach,1)
        dist = sum((pos-repmat(pos_seach(ii,:),size(pos,1),1)).^2,2);
        [~,idx(ii)]=min(dist);
    end

    %% Output
    % actually used angles
    coords_new = barumerli2023_coordinates(coords.pos(idx,:), coords.pos_type);
end



