function tfLabeler(v24,station, filename, starttime, rootdir) 
% 
% Function lets you assign labels to detected candidate boatwhistle calls 
% and save labeled spectrograms in subfolders 'bwhistle', 'other', noisy' 
% 
% tfLabel(v24, station, filename, starttime, rootdirectory) 
%
% INPUTS
%   v24 = 24kHz sampled waveform; pressure corrected, unfiltered 
%           (e.g., 2 min record) 
%   station = e.g.,'C2'
%   fname = filename,e.g.,'1074286637.190619031500.wav' used to name output
%   starttime = start time of vector, e.g., [yyyy,mm,dd,hh,MM,ss] used to
%   name output 
%   rootdir = where to store the labeled subfolders with spectrograms 
%
% HOW TO USE 
%   y=audioread('1074286637.190619040000.wav');
%   y=resample(y,1,2);  % assume original 48kHz (or 1,4 if 96 kHz) 
%   gain = 169; % typical for ST instrument 
%   y=(y-mean(y))*(10^(gain/20));  % for soundtrap instrument
%   tfLabeler(y,'S1','1074286637.190619040000.wav',[2019,06,19,04,00,00],'labeledTF');
%   would create a folder 'labeldTF' in the local directory to store 'other'
%   and bwhistle subfolders holding label images for each class. 
%
%  Calls the fuctions bandbass_del, tfMatchedFilterDet,tfMakespectro
%
% AUTHOR: 
% D. Bohnenstiehl (NCSU) 
% ToadFish Finder v.1.1 
% June 2023 

close all
fs=24000; 

%% plot the data 
  np2=12; W=2^np2; 
  lL=50; uL=3000; % filter corner for display 

     v24f=bandpass_del(v24,lL,uL,fs,4); 
     [~,F,T,Pxx]=spectrogram(v24,W,floor(W*.95),2^np2,fs); 
     PxxMod=10*log10(Pxx); 
     PxxMod=PxxMod-(smooth(median(PxxMod,2),256)*ones(1,length(T))); 
     t=(0:1:length(v24)-1)*(1/fs);

     figure('Position',[100, 100, 500,1100]); AX(1)=subplot(2,1,1); 
     imagesc(T,F,PxxMod); colormap('jet'); %caxis([35,85])
     title([station ' start: ' datestr(starttime) ])
     ylabel('Hz');
     axis xy 
     ylim([0,1800])
     AX(2)=subplot(2,1,2); 
     plot(t,v24f)
     title(['Plot of BP : ' num2str(lL) ' - ' num2str(uL) ' Hz']); xlabel('seconds'); ylabel('uPa'); 
     linkaxes(AX,'x')
 

 %% detector parameter settings
s=10; %standard deviation of harmonic for detection kernel (Hz)
sweep=3; % sweep applied to F1 harmonic (hz) 
Frange=[100 380];     % min and max range of the first harmonic 
                      % The actual bin centers are selected based on 
                      % the frequency resolution of the spectrogram 
thres=0.25;  % detection threshold used for 'MinPeakHeight' of correlation coefficient  
             % note that the 'MinPeakProminence' is defined as thres/3 

  %% run the match filter detector 
  [tmp,ffreq,score]=tfMatchedFilterDet(v24,Frange,s,sweep,thres, 220,150, 0);
  events=round(fs*tmp');  % in number of samples 
  events(events+0.45*fs >= length(v24))=length(v24)-0.45*fs; 
  etimes = events/fs;     % time in seconds within file 

  fprintf('The number of picks found in this file is %1.0f\n', length(events)) 
  input('return to move forward') 

%% for each detection 
for i=1:numel(events)
    fprintf('time:  %06.2f  freq: %05.2f  score: %0.3f\n',[etimes(i) ffreq(i) score(i)])

   % reset the main window 
     figure(1); xlim([etimes(i)-0.5,etimes(i)+0.5 ])
     hold on; subplot(2,1,2);  xline(etimes(i),'-r','LineWidth',2);
     
   % get the data in a 850 ms winow and generate scalogram 
     winstart=events(i)-0.4*fs;  % grab starting 400 ms before pick
     if winstart < 1; winstart = 1; end 
     if winstart+0.851*fs >= length(v24); winstart=length(v24)-0.851*fs; end 
     ydet=v24(winstart:winstart+0.851*fs); 
     ydetf=v24f(winstart:winstart+0.851*fs); % trimmed to 0.85*fs by tfMakespectro 
     subplot(2,1,2); ylim([-max(abs(ydetf)),max(abs(ydetf)) ])
     im=tfMakespectro(ydet,1);
     f2=figure(2); f2.Position=[700,300,500,500];  imshow(im)

CLS=input('classify call: b = Boatwhistle; o = Other; n = Noisy Boatwhistle; return = Skip:\n','s');
basename=[station,'_',char(filename), sprintf('%011.7f',etimes(i)) ];
basename=strrep(basename,'.wav','_');

clear imgLoc 
switch CLS

    case 'b' 
    disp('boatwhistle')
    imgLoc=fullfile(rootdir,'bwhistle'); 
    imFileName=['bw_', basename,'.jpg'];

    if ~exist(imgLoc, 'dir'); mkdir(imgLoc); end
    imwrite(im,fullfile(imgLoc,imFileName));


     case 'o' 
     disp('other')
     imgLoc=fullfile(rootdir,'other'); 
     imFileName=['ot_', basename,'.jpg'];

    if ~exist(imgLoc, 'dir'); mkdir(imgLoc); end
    imwrite(im,fullfile(imgLoc,imFileName));


    case 'n' 
     disp('noisy boatwhistle')
     imgLoc=fullfile(rootdir,'noisy'); 
     imFileName=['nb_', basename,'.jpg'];

    if ~exist(imgLoc, 'dir'); mkdir(imgLoc); end
    imwrite(im,fullfile(imgLoc,imFileName));                      

    otherwise 
end  % end switch 

end
  
 end