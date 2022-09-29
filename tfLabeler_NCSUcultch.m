function tfLabeler_NCSUcultch(station, starttime, rootdir) 
% 
% Special version of ppLabeler for use with NCSU cultch reef data
% This version uses the 
% 
%
% Function will let you assign labels to detected calls 
% and save labeled scalograms in subfolders 'perch', 'other', noisy' 
% 
% ppLabel(v24, station, filename, starttime, rootdirectory) 
%
% INPUTS
%   station = e.g.,   'BB'
%   starttime = start time of vector, e.g., [yyyy,mm,dd,hh,MM,ss] 
%   rootdir = output directory where to store the labeled folders 
%
% HOW TO USE 
%
% AUTHORS: 
% D. Bohnenstiehl (NCSU) 
% perch picker v.1 
% 24 Aug 22 


[v,delay_rel_rq,start_time,file]=get_PSwaveform_utc(station, starttime, 0);
v24=resample(v,1,4); 

close all
fs=24000; 

%% plot the data 
  np2=12; W=2^np2; 
  lL=50; uL=3000; % filters 
   %v24f=bandpass(v24,[lL, uL],fs,'Steepness',0.95); 
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
     %plot(t,bandpass(v24,[2000,3500],fs,'Steepness',0.95)); 
    
     plot(t,v24f)
     title(['Plot of BP : ' num2str(lL) ' - ' num2str(uL) ' Hz']); xlabel('seconds'); ylabel('uPa'); 
     linkaxes(AX,'x')
     %xlim([0,10])

 %% run the detector 
  %[events] = tfKSpicker(v24f,2,3); 


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

 [tmp]=toadfish_detector96kHz_norm1_premoval(v24,Frange,s,sweep,thres, ploton);

events=round(fs*tmp'); 
 events(events+0.03*fs >= length(v24))=length(v24)-0.03*fs; 

  etimes = events/fs; 
  fprintf('The number of picks found in this file is %1.0f\n', length(events)) 
 
 input('return to move forward') 

%% for each detection 
size(events)
for i=1:size(events,1)
etimes(i)
   % reset the main window 
     figure(1); xlim([etimes(i,1)-0.5,etimes(i,1)+0.5 ])
     hold on; subplot(2,1,2);  xline(etimes(i),'-r','LineWidth',2);
     
   % get the data in a 300 ms winow and generate scalogram 
     winstart=events(i)-0.4*fs;  % grab starting 3 ms before pick
     if winstart < 1; winstart = 1; end 
     if winstart+0.851*fs >= length(v24); winstart=length(v24)-0.85*fs; end 
     ydet=v24(winstart:winstart+0.85*fs); 
     ydetf=v24f(winstart:winstart+0.85*fs); 
     subplot(2,1,2); ylim([-max(abs(ydetf)),max(abs(ydetf)) ])
     %im=tfMakescalo(ydet); 
     im=tfMakespectro(ydet);
     f2=figure(2); f2.Position=[700,300,500,500];  imshow(im)

CLS=input('classify the call: b = Boatwhistle; o = Not Boatwhistle; n = Noisy Boatwhistle; return = Skip this one:   ','s');


basename=[station,'_',char(file), sprintf('%011.7f',etimes(i)) ];
basename=strrep(basename,'.wav','_');


clear imgLoc 
switch CLS

    case 'b' 
    disp('boatwhistle')
    imgLoc=fullfile(rootdir,'boat'); 
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