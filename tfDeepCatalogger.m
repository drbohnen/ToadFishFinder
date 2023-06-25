function [predictedLabels, DetTime, Bcount, Ocount]=tfDeepCatalogger(v24,station,fileID,rootdir,classifier,predffreq,predffreq_uncert,deleteflag,fstartdatetime) 
%  
% [predictedLabels, DetTime, Bcount, Ocount]=tfDeepCatalogger(v24,station,...
%             fileID,rootdir,classifier,predffreq,predffreq_uncert,...
%             deleteflag,fstartdatetime) 
%
% INPUT 
% v24 = is a pressure corrected waveform sampled at 24kHz unfiltered  
% station = character string with name of station or site e.g., 'CB' 
% fileID = character string with ID for recording.  If feeding in 2 min 
% files you might use the file name from which the data were extracted 
% rootdir = output directrories for 'bwhistle' and 'other' spectrograms; 
%           e.g., 'detectionfolder' 
% classifier = the name of the classification object (pre-loaded into 
% workspace) 
%
% predffreq and predffreq_uncert are estimates of the fundamental frequency
% and the +/- range (Hz) over which you want to apply a weight of 1 to the
% kernel cross correlations.  Units of Hz. 
% 
% To apply w = 1 to all kernel frequencies, use the center frequency 
% of 220 and range of 120.  The lower cutoff predffreq-predffreq_uncert
% highcutoff = predffreq+predffreq_uncert; 
%
% deleteflag = 'bwhistle', 'other', 'all', 'none' - specify which
% spectrograms NOT to keep for review.  So, 'none', keeps everything.     
%
% OUTPUT 
% predictedLabels = labels returned for each detection considered 
% DetTime = event times in seconds since beginning of file 
% Bcount = total number of boat whistle calls 
% Ocount = total number of other signals 
% 
% USAGE 
% load tfClassifier_v4.mat;  % loads the model called classifier in workspace 
% [predictedLabels, DetTime, Bcount, Ocount]=tfDeepCatalogger(v24,'DB',...
%   'file1678012426.170819220002.wav','detectionfolder',classifier,...
%    predffreq,predffreq_uncert'none',[2017,10,20,2,00,0]);
% 
%  Note that the frequency range, sweep and standard deviation (width) of 
%  kernels, as well as detection threshold can be edit in the 'detector
%  paramater settings' block at the top of this function. 
% 
%     Defaults are: 
%     s=10;    % standard deviation of harmonic (Hz) width 
%     sweep=3; % sweep applied to F1 harmonic (hz)  
%     Frange=[100 339];   % min and max range of the first harmonic 
%     thres=0.25;         % detection threshold for peak detection 
%
%  Calls the fuctions tfMatchedFilterDet,tfMakespectro
%
% D. Bohnenstiehl 
% Toadfish Finder v.1.1
% June 2023  
 
%% detector parameter settings
s=10; %standard deviation of harmonic for detection kernel (Hz)
sweep=3; % sweep applied to F1 harmonic (hz) 
Frange=[100 339];   % min and max range of the first harmonic 
                      % The actual bin centers are selected based on 
                      % the frequency resolution of the spectrogram 
                      % The upper limit is Frange(2)*2 + 3*s Hz; 
                      % The lower limit is Frange(1) - 3*s Hz; 
                      % this is done in the tfMatchedFilterDet function
                      % call from this function. 
ploton=0;    % ploton=1 plots generated, ploton=0 plots not generated
thres=0.25;  % detection threshold for peak detection routine  
             % peak prominance for peak detection is set to thres/3 in 
             % matched filter detector: tfMatchedFilterDet.m 


%% make directories to hold the image files and outputs 
tic 
out_directory=fullfile(rootdir,[station,'_',strrep(fileID,'.','_')]); 
if ~exist(out_directory, 'dir'); mkdir(out_directory); 
else 
    rmdir(out_directory,'s'); 
    disp('     removing old directory')
    mkdir(out_directory);
end
out_directoryB=fullfile(out_directory,'bwhistle');  mkdir(out_directoryB) 
out_directoryO=fullfile(out_directory,'other');  mkdir(out_directoryO) 


%% run the match filter detector 
fs=24000;
[dettimesMF,ffreq,ccscore]=tfMatchedFilterDet(v24,Frange,s,sweep,thres,predffreq,predffreq_uncert, ploton);


if isempty(dettimesMF)==0
  DetInSamples=round(fs*dettimesMF);  % in number of samples 
  DetInSamples(DetInSamples+0.4*fs >= length(v24))=length(v24)-0.4*fs;  % deal with detection at end of data file 
  DetTime = DetInSamples/fs; % times of detection set to match samples 

  % if the fstartdatetme is a vector, convert to datetime type
  if isdatetime(fstartdatetime)==0 
  fstartdatetime=datetime(fstartdatetime); 
  end 

  % --  DetTime2=fstartdatetime+seconds(DetTime);   % should work but issue on 2021b linux  
  DetTime2=datenum(fstartdatetime)+DetTime/86400;   % bug fix for 2021b linux 
  DetTime2=datetime(DetTime2,'ConvertFrom','datenum');

  disp('     writing images'); 

%% for each detection 
disp(['     # detections: ' num2str(numel(DetInSamples))])  
parfor i=1:numel(DetInSamples)
% fprintf('time:  %06.2f  freq: %05.2f  ccscore: %0.3f\n',[DetTime(i) ffreq(i) ccscore(i)])

% get the data in a 850 ms winow and generate scalogram 
     winstart=DetInSamples(i)-0.4*fs;    % grab starting 400 ms before pick
     if winstart < 1; winstart = 1; end  % first index 1 in MATLAB 
     if winstart+0.851*fs >= length(v24); winstart=length(v24)-0.851*fs; end 
     ydet=v24(winstart:winstart+0.851*fs); % data segment 
    
    
% write and save image file 
im=tfMakespectro(ydet,0);            % function will trim to 850 ms 
basename=strcat(station,'_',fileID, sprintf('%011.7f',DetTime(i)));
basename=strrep(basename,'.wav','_');
imFileName=strcat('t_', basename,'.jpg'); 
imwrite(im,fullfile(out_directoryO,char(imFileName)));
filelist{i}=fullfile(char(out_directoryO),char(imFileName));  % write all to the other directory first 

end  % end of parfor 

%% get activations for each image in the out_directory  
  scalos2eval = imageDatastore(out_directoryO); 
  net = resnet50();     % Load pretrained network
  featureLayer = 'fc1000'; % define layer for feature extraction 
  ResNet50activations = activations(net,scalos2eval, featureLayer, ...
     'MiniBatchSize', 128, 'OutputAs', 'columns');

%% predict labels using classifier and activations 
   [predictedLabels, predictedScores] = predict(classifier, ResNet50activations, 'ObservationsIn', 'columns');


%% get indices for each class 
 bb=find(predictedLabels=='bwhistle');
 Bcount=length(bb); 
 oo=find(predictedLabels=='other');
 Ocount=length(oo); 

%% write tables with detections 
 btablename=['Btable_',station,'_',strrep(fileID,'.','_')]; 
 otablename=['Otable_',station,'_',strrep(fileID,'.','_')]; 


 %  save the bwhistle class table 
 if length(bb)> 1
 eval([btablename '=table(repmat(char(station),length(bb),1),predictedLabels(bb),ffreq(bb),ccscore(bb),double(predictedScores(bb,2)),DetTime(bb),DetTime2(bb),cellstr(repmat(char(fileID),length(bb),1)), ''VariableNames'',{''Site'',''CallID'',''ffreq'',''ccscore'',''pdscore'',''rel_time'',''abs_time'',''FileID''});']);   
 eval(['save ' fullfile(out_directory, 'Btable.mat') ' ' btablename]) 
 elseif length(bb)==1  % special case since you can't make a 1 row table as of 2022a 
 eval([btablename '={char(station),predictedLabels(bb),ffreq(bb),ccscore(bb),predictedScores(bb,2),DetTime(bb),DetTime2(bb),cellstr(char(fileID))};']) 
 eval([btablename '=cell2table(' btablename ');']); 
 eval([btablename '.Properties.VariableNames={''Site'',''CallID'',''ffreq'',''ccscore'',''pdscore'', ''rel_time'',''abs_time'',''FileID''};']) 
 eval([btablename '.Site=char(' btablename '.Site);']); 
 eval([btablename '.FileID=char(' btablename '.FileID);']);
 eval(['save ' fullfile(out_directory, 'Btable.mat') ' ' btablename]) 
 end


% save the other class table 
  if length(oo)> 1
eval([otablename '=table(repmat(char(station),length(oo),1),predictedLabels(oo),ffreq(oo),ccscore(oo),predictedScores(oo,2),DetTime(oo),DetTime2(oo),cellstr(repmat(char(fileID),length(oo),1)), ''VariableNames'',{''Site'',''CallID'',''ffreq'',''ccscore'',''pdscore'',''rel_time'',''abs_time'',''FileID''});']); 
eval(['save ' fullfile(out_directory, 'Otable.mat') ' ' otablename]) 
 elseif length(oo)==1
 eval([otablename '={char(station),predictedLabels(oo),ffreq(oo),ccscore(oo),predictedScores(oo,2), DetTime(oo),DetTime2(oo),cellstr(char(fileID))};']) 
 eval([otablename '=cell2table(' otablename ');']); 
 eval([otablename '.Properties.VariableNames={''Site'',''CallID'',''ffreq'',''ccscore'',''pdscore'', ''rel_time'',''abs_time'',''FileID''};']) 
 eval([otablename '.Site=char(' otablename '.Site);']); 
 eval([otablename '.FileID=char(' otablename '.FileID);']);
 eval(['save ' fullfile(out_directory, 'Otable.mat') ' ' otablename]) 
 end

%% move the bwhistle files to the bwhistle directory 
 parfor i=1:length(bb)
    movefile(char(filelist(bb(i))),out_directoryB)
 end

%% delete scalogram images if requested 
switch deleteflag 
case 'all' 
rmdir(out_directory,'s')
disp('deleteing all scalograms')
case 'bwhistle'
rmdir(out_directoryB,'s')
disp('deleteing perch scalograms')
case 'other' 
rmdir(out_directoryO,'s')
disp('deleteing other scalograms')

otherwise 
end

else  % goes with if isempty(dettimesMF)==0; 
predictedLabels='empty'; 
DetTime=NaN;  Bcount=0; Ocount=0;
end

toc

end

