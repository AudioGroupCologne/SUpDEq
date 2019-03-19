% calculate SH coefficients from the QSC-K8 loudspeaker directivite
% provided by the GRAS database
close all; clear; clc

%% ---------------------------------- load original data - takes some time!

H = readtable('QSCK8_1x1_64442_MPS_front_pole.csv');

% get absolute values from table
h = nan(size(H,1), 31);
for nn = 1:31
    re = table2array( H(:,(nn-1)*3+2) );
    im = table2array( H(:,(nn-1)*3+4) );
    im = strrep(im, 'i,', '');
    im = str2double(im);
    h(:,nn) = sqrt( re.^2 + im.^2);
end

h = h';

f = [20    25    31.5  40    50    63    80    100   125   160   200   250  ...
     315   400   500   630   800   1000  1250  1600  2000  2500  3150  4000 ...
     5000  6300  8000  10000 12500 16000 20000];

% get micrphone positions
H = H.Var1;

phi   = zeros(size(H));
theta = phi;

for nn = 1:numel(H)
    phi(nn)   = str2double(H{nn}(2:4));
    theta(nn) = str2double(H{nn}(6:8));
end

clear H nn re im

%% ---------------- get a subset of the grid to speed up further processing
sg      = AKgreatCircleGrid(90:-2:-90, 2, 90, 0);
sg(:,2) = 90 - sg(:,2);

h_sub = zeros(size(h, 1), size(sg,1));

for nn = 1:size(sg,1)
    id          = phi == sg(nn,1) & theta == sg(nn,2);
    h_sub(:,nn) = h(:, id);
end

% select some points for plotting intermediate, and final results
id = sg(:,1) == 0 & (sg(:,2) == 0 | sg(:,2) == 90 | sg(:,2) == 180);
id = find(id == 1);

clear nn

%% ------------------------------------------ spherical harmonics transform
clear Obj

% prepare struct for saving data first
Obj.GLOBAL_Title              = 'Loudspeaker directivity of an on axis equalized QSC-K8';
Obj.GLOBAL_Comment            = '';
Obj.GLOBAL_Origin             = 'part of AKtools: www.ak.tu-berlin.de/AKtools';
Obj.GLOBAL_DataType           = 'Spherical harmonics coefficients';
Obj.GLOBAL_SHconvention       = 'E. G. Williams: Fourier Acoustics. Academic Press (1999)';
Obj.GLOBAL_MetaData           = 'Most meta data entries are in accordance to AES69-2015: AES standard for file exchange - Spatial acoustic data file format. AES Standards Comittee, Audio Engineering Society, Inc. (2015)';
Obj.GLOBAL_History            = 'raw data from GRAS database (third octave pressure values), processed by QSC_K8_directivity.m';
Obj.GLOBAL_RoomType           = 'free field';
Obj.GLOBAL_AuthorContact      = 'fabian.brinkmann@tu-berlin.de';
Obj.GLOBAL_Organization       = 'Audio Communication Group, TU Berlin, Germany';
Obj.GLOBAL_License            = 'EUPL';
Obj.GLOBAL_DateCreated        = datestr(now, 'yyyy-mmm-dd HH:MM:SS');
Obj.GLOBAL_References         = 'GRAS database on www.depositOnce.tu-berlin.de';
Obj.GLOBAL_ApplicationName    = 'Matlab';
Obj.GLOBAL_ApplicationVersion = version;
Obj.ReceiverPosition          = [0 0 0];
Obj.ReceiverPosition_Type     = 'cartesian';
Obj.ReceiverPosition_Units    = 'metre';
Obj.N                         = f;
Obj.N_LongName                = 'frequency';
Obj.N_Units                   = 'Hz';

% SH parameters
Obj.SH.coeff            = 0;
Obj.SH.order            = 20;
[~, Obj.SH.n, Obj.SH.m] = AKsh(Obj.SH.order, [], 0, 0);
Obj.SH.f                = f;
Obj.SH.SHTmode          = 'complex';
Obj.SH.doFFT            = false;
Obj.SH.fs               = 44100;
Obj.SH.compact          = false;
Obj.SH.SHmode           = 'complex';
Obj.SH.isEven           = true;

% flip the azimuth of the front pole grid to work with a vertical polar
% grid in the next steps
sg2      = sg;
sg2(:,1) = mod(180 - sg2(:,1), 360);

% sh transform
f_nm = AKsht(h_sub, Obj.SH.doFFT, sg2, Obj.SH.order, Obj.SH.SHTmode, Obj.SH.fs, Obj.SH.compact, Obj.SH.SHmode);

% inverse transform to check if the order is sufficient
h_sh = real( AKisht(f_nm, Obj.SH.doFFT, sg2(id,:), Obj.SH.SHTmode) );

% plot
AKf(20,10)
semilogx(f, db( h_sub(:,id) ), 'k', 'linewidth', 2);
hold on
semilogx(f, db( h_sh ), '--r', 'linewidth', 2);

clear h_sh sg2

%% ------------------------------------------------- rotate to desired grid
%  the data were originally saved in a front pole grid, and is now
%  transformed to a vertical polar (often calles spherical) coordinate
%  system

g_nm = AKshRotate(f_nm, [0 90 0]);

% check if everything is in place
%      new az. new el. = old az. old el. 
sg2 = [     0      90         0       0    % frontal
           90      90        90      90    % left
          270      90       270      90    % right
          180      90         0     180    % back
            0       0         0      90    % up
            0     180       180      90    % down
          ];

AKf(40,20)
hold on
for nn = 1:size(sg2,1)
    % old sampling grid
    h_sh = AKisht(f_nm, Obj.SH.doFFT, [mod(180-sg2(nn,3), 360) sg2(nn,4)], Obj.SH.SHTmode, Obj.SH.isEven, Obj.SH.compact, Obj.SH.SHmode);
    semilogx(f, db( h_sh ), 'k', 'linewidth', 2);
    % new sampling grid
    h_sh = AKisht(g_nm, Obj.SH.doFFT, sg2(nn,1:2), Obj.SH.SHTmode, Obj.SH.isEven, Obj.SH.compact, Obj.SH.SHmode);
    semilogx(f, db( h_sh ), '--r', 'linewidth', 2);
end
set(gca, 'XScale', 'log')

clear nn h_sh sg2

Obj.SH.coeff = g_nm;

%% -------------------------------- normalize to 0 dB and frontal direction

hFront       = AKisht(Obj.SH.coeff, Obj.SH.doFFT, [0 90], Obj.SH.SHTmode);
gain         = mean(real(hFront(f>=100 & f<=1e3)));
Obj.SH.coeff = Obj.SH.coeff / gain;
hFront       = AKisht(Obj.SH.coeff, Obj.SH.doFFT, [0 90], Obj.SH.SHTmode);

AKf(20,10)
semilogx(f, db( hFront ), 'k', 'linewidth', 2);

clear gain hFront

%% --------------------------------------------------- save SH coefficients

save('QSC_K8_directivity.mat', '-struct', 'Obj')