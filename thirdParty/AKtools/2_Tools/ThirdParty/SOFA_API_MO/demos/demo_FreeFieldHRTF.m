% Demonstrates the usage of spherical harmonics (SH) for HRTF interpolation. 
% demo_FreeFieldHRTF loads an HRTF set, transforms to TF, then to SH, then 
% samples the horizontal and median plane in steps of 0.5 degrees. Finally, 
% the files are saved in SOFAdbPath as demo_FreeFieldHRTF_{TF, SH, TFrec}.sofa.

% SOFA API - demo script
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 



IR=SOFAload('db://database/thk/HRIR_L2354.sofa');
fs=IR.Data.SamplingRate;
IR.GLOBAL_APIVersion=SOFAgetVersion;
figure;
SOFAplotHRTF(IR,'magmedian'); title('SimpleFreeFieldHRIR (FIR) for reference');

SOFAsave(fullfile(SOFAdbPath,'sofa_api_mo_test','demo_FreeFieldHRTF_1_IR.sofa'),IR);

%% Convert to TF
TF=SOFAgetConventions('SimpleFreeFieldHRTF');
TF.ListenerPosition=IR.ListenerPosition;
TF.ListenerPosition_Type=IR.ListenerPosition_Type;
TF.ListenerPosition_Units=IR.ListenerPosition_Units;
TF.ListenerView=IR.ListenerView;
TF.ListenerView_Type=IR.ListenerView_Type;
TF.ListenerView_Units=IR.ListenerView_Units;
TF.ListenerUp=IR.ListenerUp;
TF.SourcePosition=IR.SourcePosition;
TF.SourcePosition_Type=IR.SourcePosition_Type;
TF.SourcePosition_Units=IR.SourcePosition_Units;
TF.EmitterPosition=IR.EmitterPosition;
TF.EmitterPosition_Type=IR.EmitterPosition_Type;
TF.EmitterPosition_Units=IR.EmitterPosition_Units;
TF.ReceiverPosition=IR.ReceiverPosition;
TF.ReceiverPosition_Type=IR.ReceiverPosition_Type;
TF.ReceiverPosition_Units=IR.ReceiverPosition_Units;

TF.Data.Real=zeros(IR.API.M,IR.API.R,IR.API.N+1);
TF.Data.Imag=zeros(IR.API.M,IR.API.R,IR.API.N+1);
for ii=1:IR.API.M
  for jj=1:IR.API.R
   sp=fft(squeeze(IR.Data.IR(ii,jj,:)),2*IR.API.N); % Delay not considered!
   TF.Data.Real(ii,jj,:)=real(sp(1:IR.API.N+1,:));
   TF.Data.Imag(ii,jj,:)=imag(sp(1:IR.API.N+1,:));
  end
end
TF.N=(0:fs/2/IR.API.N:fs/2)';

TF=SOFAupdateDimensions(TF);

SOFAsave(fullfile(SOFAdbPath,'sofa_api_mo_test','demo_FreeFieldHRTF_2_TF.sofa'),TF);

%% Plot median plane and horizontal planes for reference
figure;
SOFAplotHRTF(TF,'magmedian'); title('SimpleFreeFieldHRTF (TF) for reference');
figure;
SOFAplotHRTF(TF,'etchorizontal'); title ('SimpleFreeFieldHRTF (TF) for reference');

%% Convert to an emitter-based representation, TFE
TFE=TF; 
TFE.GLOBAL_SOFAConventions = 'GeneralTF-E';
TFE.GLOBAL_DataType = 'TF-E';
TFE.API.E=TF.API.M;
TFE.API.M=1;
TFE.Data=rmfield(TFE.Data,{'Real','Imag'});
TFE.Data.Real(1,:,:,:)=shiftdim(TF.Data.Real,1); % MRN --> 1RNM --> MRNE with M=1
TFE.API.Dimensions.Data.Real='MRNE';
TFE.Data.Imag(1,:,:,:)=shiftdim(TF.Data.Imag,1);
TFE.API.Dimensions.Data.Imag='MRNE';
TFE.EmitterPosition=TF.SourcePosition;
TFE.EmitterPosition_Type=TF.SourcePosition_Type;
TFE.EmitterPosition_Units=TF.SourcePosition_Units;
TFE.API.Dimensions.EmitterPosition='ECI';
TFE.SourcePosition=[0 0 0];
TFE.API.Dimensions.SourcePosition='IC';

TFE=SOFAupdateDimensions(TFE);

SOFAsave(fullfile(SOFAdbPath,'sofa_api_mo_test','demo_FreeFieldHRTF_3_TFE.sofa'),TFE);

%% Convert to SH
SH=TFE;
SH.GLOBAL_SOFAConventions = 'FreeFieldHRTF';

Lmax=floor(sqrt(size(SH.EmitterPosition,1))-1); % Max SH order
L=40; % actual SH order
[S, SH.API.E]=sph2SH(SH.EmitterPosition(:,1:2), L);

Sinv=pinv(S);
SH.Data.Real=zeros(1, SH.API.R, SH.API.N, SH.API.E);
SH.Data.Imag=zeros(1, SH.API.R, SH.API.N, SH.API.E);
for ii=1:TFE.API.R
  for jj=1:TFE.API.N
   SH.Data.Real(1,ii,jj,:)=Sinv*squeeze(TFE.Data.Real(1,ii,jj,:));
   SH.Data.Imag(1,ii,jj,:)=Sinv*squeeze(TFE.Data.Imag(1,ii,jj,:));
  end
end

SH.EmitterPosition=mean(SH.EmitterPosition);
SH.EmitterPosition_Type='Spherical Harmonics';

SH = SOFAupdateDimensions(SH);

SOFAsave(fullfile(SOFAdbPath,'sofa_api_mo_test','demo_FreeFieldHRTF_4_SH.sofa'),SH);

%% plot median and horitonal planes - spatially continuous
figure;
SOFAplotHRTF(SH,'magmedian'); title('FreeFieldHRTF (TFE) in Spherical Harmonics');
figure;
SOFAplotHRTF(SH,'etchorizontal'); title ('FreeFieldHRTF (TFE) in Spherical Harmonics');

%% plot spatially continuous geometry
SOFAplotGeometry(SH);
title ('FreeFieldHRTF (TFE) in Spherical Harmonics');

%% plot spatial spectra
figure;
for y=1:4
  subplot(2,2,y); hold on;
  fax=0:fs/(2*SH.API.N):(fs-(fs/(2*SH.API.N)))/2;
  r=1; plot(fax,20*log10(squeeze(abs(SH.Data.Real(1,r,:,y)+1i*SH.Data.Imag(1,r,:,y)))));
  r=2; plot(fax,20*log10(squeeze(abs(SH.Data.Real(1,r,:,y)+1i*SH.Data.Imag(1,r,:,y)))),'r');
  title(['SH spectra, coefficient index (ACN): ' num2str(y)]);
  xlabel('Frequency (Hz)');
  ylabel('Magnitude (dB)');
end
subplot(2,2,1);legend('First receiver','Second receiver');

%% plot all coefficients for a given frequency
figure; hold on;
f=10000; k=round(f/fs*2*SH.API.N);
plot(squeeze(20*log10(abs(SH.Data.Real(1,1,k,:)+1i*SH.Data.Imag(1,1,k,:)))),'ob');
hold on;
plot(squeeze(20*log10(abs(SH.Data.Real(1,2,k,:)+1i*SH.Data.Imag(1,2,k,:)))),'xr');
title(['SH representation, frequency: ' num2str(f) ' Hz']);
xlabel('Coefficients (ACN index)');
ylabel('Magnitude (dB)');
legend('First receiver','Second receiver');

%% interpolate for the horizontal and median planes to FreeFieldHRTF (TFE)
TFEint=SH;
elemin=-90;
ele=[elemin:1:90 89:-1:elemin zeros(1,length(1:355))]';
azi=[zeros(length(elemin:1:90),1); 180*ones(length(89:-1:elemin),1); (1:355)'];
radius=SH.EmitterPosition(:,3)*ones(size(ele));
TFEint.EmitterPosition=[azi ele radius];
TFEint.EmitterPosition_Type='spherical';
TFEint.EmitterPosition_Units=SH.EmitterPosition_Units;
Sint = sph2SH(TFEint.EmitterPosition(:,1:2), L);
TFEint.API.E=size(Sint,1);
TFEint.Data.Real=zeros(1,2,TFEint.API.N,TFEint.API.E);
TFEint.Data.Imag=zeros(1,2,TFEint.API.N,TFEint.API.E);
for ii=1:TFEint.API.R
  for jj=1:TFEint.API.N
    TFEint.Data.Real(1,ii,jj,:)=Sint*squeeze(SH.Data.Real(1,ii,jj,:));
    TFEint.Data.Imag(1,ii,jj,:)=Sint*squeeze(SH.Data.Imag(1,ii,jj,:));
  end
end

TFEint=SOFAupdateDimensions(TFEint);

SOFAsave(fullfile(SOFAdbPath,'sofa_api_mo_test','demo_FreeFieldHRTF_5_TFEint.sofa'),TFEint);

%% interpolate for the horizontal and median planes to SimpleFreeFieldHRTF (TF)
TFint=TF;
ele=[-90:0.5:90 89:-.5:-90 zeros(1,length(1:0.5:355))]';
azi=[zeros(length(-90:.5:90),1); 180*ones(length(89:-.5:-90),1); (1:0.5:355)'];
radius=1.2*ones(size(ele));
TFint.SourcePosition=[azi ele radius];
Sint = sph2SH(TFint.SourcePosition(:,1:2), sqrt(SH.API.E)-1);
TFint.API.M=size(Sint,1);
TFint.Data.Real=zeros(TFint.API.M,2,TFint.API.N);
TFint.Data.Imag=zeros(TFint.API.M,2,TFint.API.N);
for ii=1:TFint.API.R
  for jj=1:TFint.API.N
    TFint.Data.Real(:,ii,jj)=Sint*squeeze(SH.Data.Real(1,ii,jj,:));
    TFint.Data.Imag(:,ii,jj)=Sint*squeeze(SH.Data.Imag(1,ii,jj,:));
  end
end

TFint=SOFAupdateDimensions(TFint);

SOFAsave(fullfile(SOFAdbPath,'sofa_api_mo_test','demo_FreeFieldHRTF_6_TFint.sofa'),TFint);

%% compare
figure;
SOFAplotHRTF(TFint,'magmedian'); title('SimpleFreeFieldHRTF (TF): Interpolated');
figure;
SOFAplotHRTF(TFint,'etchorizontal'); title('SimpleFreeFieldHRTF (TF): Interpolated');



