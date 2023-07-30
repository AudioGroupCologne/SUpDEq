%demo_ExampleFiles - Prepare a set of SOFA example files with different conventions.
% 
% demo_ExampleFiles loads or creates SOFA files. Some (meta) data might be modified. The goal is to get a collection of SOFA files according to all conventions published in AES69-2022 (SOFA 2.0).

% #Author: Michael Mihocic (21.09.2022)
% 
% SOFA Toolbox
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

%% Prepare
SOFAstart;
out=[fullfile(mfilename('fullpath'), '..') '\examples\'];
mkdir(out);

%% Load and save SOFA files
% FreeFieldDirectivityTF_1.0
Obj=SOFAload('db://database/tu-berlin%20(directivity)/Trumpet_modern_a4_fortissimo.sofa');
SOFAsave([out Obj.GLOBAL_SOFAConventions '_' Obj.GLOBAL_SOFAConventionsVersion '.sofa'],Obj);
% FreeFieldDirectivityTF_1.1
Obj=SOFAload('file.sofa');
SOFAsave([out Obj.GLOBAL_SOFAConventions '_' Obj.GLOBAL_SOFAConventionsVersion '.sofa'],Obj);
% FreeFieldHRIR_1.0

% FreeFieldHRTF_1.0

% GeneralFIR-E_2.0

% GeneralFIRE_1.0

% GeneralFIR_1.0

% GeneralSOS_1.0

% GeneralString_0.2

% GeneralTF-E_1.0

% GeneralTF_1.0

% GeneralTF_2.0

% General_1.0

% SimpleFreeFieldHRIR_1.0
Obj=SOFAload('db://database/thk/HRIR_L2354.sofa');
SOFAsave([out Obj.GLOBAL_SOFAConventions '_' Obj.GLOBAL_SOFAConventionsVersion '.sofa'],Obj);

% SimpleFreeFieldHRSOS_1.0

% SimpleFreeFieldHRTF_1.0

% SimpleFreeFieldSOS_1.0

% SimpleHeadphoneIR_1.0

% SingleRoomMIMOSRIR_1.0

% SingleRoomSRIR_1.0
