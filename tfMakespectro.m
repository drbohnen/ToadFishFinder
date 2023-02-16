function im=tfMakespectro(v24,ploton)
%
% Function takes a vector of acoustic data 850 ms long sampled at 24 kHz 
% and returns a time frequency reassigned spectrogram image of the data 
% resample to 224 x 224 x 3; suitable for RESNET50   
%
% INPUT 
% v24 is a acoustic data segment sampled at 24 kHz. If longer than 
% 850 ms (20400 points) it will be trimmed. Input segment cannot to smaller. 
%
% ploton = 1 will display the center line column of the spectrogram 
%          This is used for scaling image amplitudes in the image 
%
% OUTPUT 
% im is a 224 x 224 x 3 image matrix 
%
%
% AUTHORS: 
% D. Bohnenstiehl (NCSU) 
% ToadFish Finder v.1 
% Oct 2022 


%% initial set up 
fs=24000; % sample rate 
data=v24(1:0.85*fs);  % trim the data to 30 ms 
[m,n] = size(v24);  % get length of data 

if m > n
data=data';  % transform to make a row vector 
end
L=length(data); 

%% spectrogram parameters 
frequencyLimits = [0 1200]; % Hz
leakage = 0; 
timeResolution = (2^12)/24000; 
overlapPercent = 85; 
reassignFlag = true;
timeValues = (0:length(v24)-1)./fs;

%% Compute spectral estimate
[P,F,T] = pspectrum(v24,timeValues, ...
    'spectrogram', ...
    'FrequencyLimits',frequencyLimits, ...
    'Leakage',leakage, ...
    'TimeResolution',timeResolution, ...
    'OverlapPercent',overlapPercent, ...
    'Reassign',reassignFlag);
P=10*log10(P);               % convert to dB scale 
P=flipud(P); F=flipud(F);    % flip upside down 
                             

% P is size 1024 x 30 
% very long wavelength background subtraction 
 PMod=P-(smooth(median(P,2),256)*ones(1,length(T))); 

%% image prep and return    
% the center line is used for scaling the image 
mid=floor(size(PMod,2)*0.55); 
scline=smooth(PMod(:,mid),17,'rlowess'); 

if ploton==1  % plotting   
figure(3); close(3); figure(3); 
plot(F,scline); hold on; 
end

[pamp,ploc]=findpeaks(scline,'MinPeakProminence',3);

mininput=-4; 
if ~isempty(pamp) 
if ploton==1; plot(F(ploc),pamp,'o'); end % plotting    
spamp=sort(pamp); 
maxinput=spamp(end-min([2,numel(spamp)-1])); 
else  
 maxinput=max(scline); 
end

if mininput >= maxinput
    mininput = min(scline);
end
    
    im = ind2rgb(im2uint8(rescale(PMod,0,1,'InputMin',mininput,'InputMax',maxinput)),parula);
    im=imresize(im,[224 224]);

end




