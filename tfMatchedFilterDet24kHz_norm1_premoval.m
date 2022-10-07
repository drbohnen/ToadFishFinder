function  [det_time,det_freq1,det_freq2,det_score]=tfMatchedFilterDet24kHz_norm1_premoval(y,Frange,s,sweep,thres, ploton, usegpu)
% toadfish matched filter detector using 
% background subtraction, percussive filtering 
% and xcorr normalization 
% 
% INPUTS 
% y = demeaned & pressure corrected time series sampled at 24 kHz  
% Frange = is an a x 2 vector with [minF1  maxF1] (range of fundamental
% frequency kernels to run 
% sweep is the change in frequency (Hz) over the duration of the call 
% s is the width or std of call in Hz 
% thres = threshold for correlation to retain, e.g., 0.25 
% ploton = 1 will plot detections and spectrograms 
%
%
% OUTPUTS 
% det_time = time in seconds within the time series y 
% det_freq1 = frequency of first harmonic 
% det_freq2 = frequency of second harmonc = 2*det_freq1; 
% det_score = amplitude of the detection score 
% 
%  Example 
%  [det_time,det_freq1,det_freq2,det_score]=...
%  tfMatchedFilterDet24kHz_norm1_premoval(v24,[100,337],10,3,0.25,1);
% 
% Modified from Ricci et al. Plos One (2017)
% Oyster toadfish (Opsanus tau) boatwhistle call detection 
% and patterns within a large-scale oyster restoration site
%
% ToadFishFinder v1 
% Oct 2022 
% 

tic 

%% calculate PSD 
W=4096;  fs=24000; np2=13;  % zero pad to next  % zero pad to next 

%y=gpuArray(y);  % uncomment for GPU 

[~,F,T,Pxx]=spectrogram(y,W,floor(W*.8),2^np2,fs);   %  generate spectrogram 
a=find(F >=80 & F <= ceil(2*Frange(2)+2*s));  F=F(a); Pxx=Pxx(a,:);   % trim the frequency band requested  
[rowsF,~]=size(F); % number of frequencies in the spectrogram 


%% determine the harmonics for the kernels based on F1range 
F1a=F(F> Frange(1) & F < Frange(2));  F1b=F1a-sweep; F1=[F1a'; F1b']; 
F2=F1*2;  % define 2nd harmionic 
[~,nkerns]=size(F1);  % determine the number of kernels 

%% adjust the data spectrogram 
%PxxMod=10*log10(gather(Pxx)); % uncomment for GPU 
PxxMod=10*log10((Pxx)); % comment for GPU
% filter percussive signals. 
PxxModmed=median(PxxMod,1); 
PxxMod=PxxMod-repmat(PxxModmed,size(PxxMod,1),1);
% background subtraction 
PxxMod=PxxMod-(smooth(median(PxxMod,2),31)*ones(1,length(T))); 

% make a copy without percussive filtering for plotting 
if ploton==1 
%PxxMod2=10*log10(gather(Pxx));  % uncomment for GPU 
PxxMod2=10*log10((Pxx));
PxxMod2=PxxMod2-(smooth(median(PxxMod2,2),31)*ones(1,length(T)));  
end 

% figure out size of the kernel 
klengthsec=0.225;
t=0:diff(T(1:2)):klengthsec;  d=max(t); 
len_k=length(t);  % nunber of colums in kernel 

% preallocate 
cmatrix=nan(nkerns,length(T)+length(t)-1);   % length of T + len_k -1; 

%% now generate kernels 
for j=1:nkerns  % make each kernel
    % upper one 
fo=F2(1,j); f1=F2(2,j); 
k1=nan(length(F),length(t)); 
for i=1:length(t)  % for each time step 
x=F-(fo+(t(i)./d)*(f1-fo));
k1(:,i)=(1-(x.^2)/(s.^2)).*exp(-(x.^2)/(2*s^2)); 
end

% lower one 
fo=F1(1,j);  f1=F1(2,j); 
%t=0:diff(T(1:2)):0.225;  d=max(t); 
k2=nan(length(F),length(t)); 
for i=1:length(t)
x=F-(fo+(t(i)./d)*(f1-fo));
k2(:,i)=(1-(x.^2)/(s.^2)).*exp(-(x.^2)/(2*s^2)); 
end

k=k1+k2;  % add the result together to make the final c
knorm=sqrt(sum(k(:).^2));  

%% now cross correlate with data spectrogram 
[c]=xcorr2(PxxMod,k);    
c=c(rowsF,:);  % pulls the center value that represents only a time shift 

% norm to -1 to +1 
PxxModPad=padarray(PxxMod,[0, len_k-1],1,'both'); 

Pnorm=nan(1,length(cmatrix));
for jj=1:length(PxxModPad)-(len_k-1)
    Pnorm(jj)=sqrt(sum(reshape(PxxModPad(:,jj:jj+len_k-1),numel(k),1).^2));
end

cmatrix(j,:)=c./(Pnorm*knorm);   % stores it as a row 

end   % repeat for each of the nkerns 

[cout,iout]=max(cmatrix(:,len_k:end));  % returns the max value (cout) & which row gave the max (iout)

[det_score,LOCS]=findpeaks(cout,'MinPeakHeight',thres,'MinPeakDistance',4,'MinPeakProminence',thres/3); 
det_time=T(LOCS)'; 
det_freq1=F1(1,iout(LOCS))';
det_freq2=F2(1,iout(LOCS))';
det_score=det_score'; 


%% if plotting is on 
if ploton==1 
figure; ax(1)=subplot(2,1,1); imagesc(T,F,PxxMod2);  axis xy; hold on; 
colormap('jet'); caxis([-50,50]); 
plot(T-diff(T(1:2))/2,100+cout*500,'k','LineWidth',1); grid on; 
plot(T(LOCS)-diff(T(1:2))/2,det_freq1,'ok','MarkerSize',4,'LineWidth',1);
plot(T(LOCS)-diff(T(1:2))/2,det_freq2,'ok','MarkerSize',4,'LineWidth',1);
ax(2)=subplot(2,1,2);  plot(T-diff(T(1:2))/2,cout); 
hold on; plot(T(LOCS)-diff(T(1:2))/2,det_score,'or'); grid on; 
for j=1:length(LOCS)
text(T(LOCS(j)),det_score(j),char(sprintf('%3.1f\n',det_score(j))),'FontSize',14); 
end

linkaxes([ax(2),ax(1)],'x'); 

end

toc

