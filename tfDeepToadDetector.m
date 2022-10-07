function [predictedLabels, etimes, Bcount, Ocount]=tfDeepToadDetector(v24,station,fileID,rootdir,classifier,deleteflag,fstartdatetime) 
% INPUT 
% v24 = is a pressure corrected waveform sampled at 24kHz unfiltered  
% station = character string with name of station or site e.g., 'CB' 
% fileID = character string with ID for recording.  If feeding in 2 min 
% files you might use the file name from which the data were extracted 
% rootdir = output directrories for 'boat' and 'other' spectrograms; 
% classifier = the name of the classification object (pre-loaded into 
% workspace) 
% deleteflag = 'boat', 'other', 'all', 'none' - specify which spectrogram
% not to keep for review  
%
% OUTPUT 
% predictedLabels = labels returned for each detection considered 
% etimes = event times in seconds since beginning of file 
% Bcount = total number of boatwhistle 
% Ocount = total number of other signals 
% 
% USAGE 
% load TFFclassifier_v1.mat;  % load the classifier 
% [predictedLabels, etimes, Bcount, Ocount]=tfDeepToadDetector(v24,'DB','1678012426.170819220002.wav','detectionfolder',classifier,'none',[2017,10,20,2,00,0]);

fs=24000; 
%% input detector parameter settings
s=10; %standard deviation of harmonic for detection kernel (Hz)
sweep=3; % sweep applied to F1 harmonic (hz) 
Frange=[100 337];   % min and max range of the first harmonic 
                      % The actual bin centers are selected based on 
                      % the frequency resolution of the spectrogram 
                      % note that the spectrogram upper limit is 690 Hz 
                      % so if you increase the max freq. here you need to adjust 
                      % that as well, so that the upper limit is
                      % Frange(2)*2 + ~20 Hz; 
ploton=0;    %ploton=1 plots generated, ploton=0 plots not generated
thres=0.25;  %detection threshold  
% note that the 'MinPeakProminence' is defined as thres/3 in the detection
% function, but this could be changed. 

%% make directories to hold the image files and outputs 
out_directory=fullfile(rootdir,[station,'_',strrep(fileID,'.','_')]); 
 if ~exist(out_directory, 'dir'); mkdir(out_directory); end
out_directoryB=fullfile(out_directory,'boat');  mkdir(out_directoryB) 
out_directoryO=fullfile(out_directory,'other');  mkdir(out_directoryO) 


% run the match filter detector 
  [tmp,ffreq,~,ccscore]=tfMatchedFilterDet24kHz_norm1_premoval(v24,Frange,s,sweep,thres, ploton);

  
if isempty(tmp)==0
  
  
  events=round(fs*tmp);  % in mumber of samples 
  events(events+0.4*fs >= length(v24))=length(v24)-0.4*fs; 
  etimes = events/fs; 

  if isdatetime(fstartdatetime)==0 
  fstartdatetime=datetime(fstartdatetime); 
  end 
  etimes2=fstartdatetime+seconds(etimes); 

  disp('writing images'); 



% for each detection 
size(events)
for i=1:size(events,1)
    fprintf('time:  %06.2f  freq: %05.2f  ccscore: %0.3f\n',[etimes(i) ffreq(i) ccscore(i)])

% get the data in a 850 ms winow and generate scalogram 
     winstart=events(i)-0.4*fs;  % grab starting 400 ms before pick
     if winstart < 1; winstart = 1; end 
     if winstart+0.851*fs >= length(v24); winstart=length(v24)-0.851*fs; end 
     ydet=v24(winstart:winstart+0.851*fs); 
     im=tfMakespectro(ydet,0);
    

basename=strcat(station,'_',fileID, sprintf('%011.7f',etimes(i)));
basename=strrep(basename,'.wav','_');
imFileName=strcat('t_', basename,'.jpg'); 
imwrite(im,fullfile(out_directoryO,char(imFileName)));
filelist{i}=fullfile(char(out_directoryO),char(imFileName));  % write all to the other directory first 

end

%% get activations for each image in the out_directory  
  scalos2eval = imageDatastore(out_directoryO); 
  net = resnet50();     % Load pretrained network
  featureLayer = 'fc1000'; % define layer for feature extraction 
  ResNet50activations = activations(net,scalos2eval, featureLayer, ...
     'MiniBatchSize', 128, 'OutputAs', 'columns');

%% predict labels using classifier and activations 
   [predictedLabels, predictedScores] = predict(classifier, ResNet50activations, 'ObservationsIn', 'columns');

   toc 

%% get indices for each class 
 bb=find(predictedLabels=='boat'); Bcount=length(bb); 
 oo=find(predictedLabels=='other'); Ocount=length(oo); 

%% write tables with detections 
 btablename=['Btable_',station,'_',strrep(fileID,'.','_')]; 
 otablename=['Otable_',station,'_',strrep(fileID,'.','_')]; 
 

 % 
 if length(bb)> 1
 eval([btablename '=table(repmat(char(station),length(bb),1),predictedLabels(bb),ffreq(bb),ccscore(bb),double(predictedScores(bb,2)),etimes(bb),etimes2(bb),repmat(char(fileID),length(bb),1), ''VariableNames'',{''Site'',''CallID'',''ffreq'',''ccscore'',''pdscore'',''rel_time'',''abs_time'',''FileID''});']);   
 eval(['save ' fullfile(out_directory, 'Btable.mat') ' ' btablename]) 
 elseif length(bb)==1  % special case since you can't make a 1 row table in 2022a 
 eval([btablename '={char(station),predictedLabels(bb),ffreq(bb),ccscore(bb),predictedScores(bb,2),etimes(bb),etimes2(bb),char(fileID)};']) 
 eval([btablename '=cell2table(' btablename ');']); 
 eval([btablename '.Properties.VariableNames={''Site'',''CallID'',''ffreq'',''ccscore'',''pdscore'', ''rel_time'',''abs_time'',''FileID''};']) 
 eval(['save ' fullfile(out_directory, 'Btable.mat') ' ' btablename]) 
 end

  if length(oo)> 1
eval([otablename '=table(repmat(char(station),length(oo),1),predictedLabels(oo),ffreq(oo),ccscore(oo),predictedScores(oo,2),etimes(oo),etimes2(oo),repmat(char(fileID),length(oo),1), ''VariableNames'',{''Site'',''CallID'',''ffreq'',''ccscore'',''pdscore'',''rel_time'',''abs_time'',''FileID''});']); 
eval(['save ' fullfile(out_directory, 'Otable.mat') ' ' otablename]) 
 elseif length(oo)==1
 eval([otablename '={char(station),predictedLabels(oo),ffreq(oo),ccscore(oo),predictedScores(oo,2), etimes(oo),etimes2(oo),char(fileID)};']) 
 eval([otablename '=cell2table(' otablename ');']); 
 eval([otablename '.Properties.VariableNames={''Site'',''CallID'',''ffreq'',''ccscore'',''pdscore'', ''rel_time'',''abs_time'',''FileID''};']) 
 eval(['save ' fullfile(out_directory, 'Otable.mat') ' ' otablename]) 
 end

%% move the perch files to the perch directory 
 parfor i=1:length(bb)
    movefile(char(filelist(bb(i))),out_directoryB)
 end



%% delete scalogram images if requested 
switch deleteflag 
case 'all' 
rmdir(out_directory,'s')
disp('deleteing all scalograms')
case 'boat'
rmdir(out_directoryB,'s')
disp('deleteing perch scalograms')
case 'other' 
rmdir(out_directoryO,'s')
disp('deleteing other scalograms')

otherwise 
end

else  % goes with if isempty(tmp)==0; 
predictedLabels='empty'; 
etimes=NaN;  Bcount=0; Ocount=0;
end


end

