function  [det_time,det_freq1,det_score]=tfMatchedFilterDet(y,frange,s,sweep,thres, predffreq,predffreq_uncert, ploton)
% Toadfish detector using spectrogram cross correlation..
% 
% Spectrograms processed with background subtraction, 
% percussive filtering and xcorr normalization 
%
%  Example:
%  [det_time,det_freq1,det_score]=tfMatchedFilterDet(y,frange,s,...
%                    sweep,thres, predffreq,predffreq_uncert, ploton)
%  [det_time,det_freq1,det_score]=tfMatchedFilterDet(v24,[100,339],10,3,0.25,250,75,0);
% 
% INPUTS 
%     y = demeaned & pressure corrected time series sampled at 24 kHz  
%    frange = is an a 1 x 2 vector with [minF1  maxF1] range of 
%             fundamental frequencies to consider (e.g., [100, 400]).    
%    s = width of the harmonic in Hz; use to build kernel (e.g., 10) 
%    sweep = change in freq (Hz) of first harmonic over the 225 ms 
%            duration of kernel (e.g., 3)
%    thres = threshold for correlation to retain (e.g., 0.25) 
%    predffreq = predicted fo frequency Hz of first harmonic (e.g., 250)
%    predffreq_uncertainity = how certain about the fo prediction (e.g., 75)
%    ploton = 1 will plot detections and spectrograms (e.g., 0, for no plots).
% 
%  Correlations are downweighted beyond the predffreq +/- predffreq_uncert;
%  The predicted frequency of first harmonic can be determined using water temperature. 
%  
%  If you want to turn this off, set the prediction and its uncertainity 
%  to the median of your frequency range. For example: 250, 250. 
% 
% [det_time,det_freq1,det_score]=tfMatchedFilterDet(v24,[100,339],10,3,0.25,250,250,0);
%
% OUTPUTS 
%    det_time = time in seconds within the time series y 
%    det_freq1 = frequency of first harmonic 
%    det_score = normalized cross correlation detection score 
% 
% 
% Concept is after Ricci et al. Plos One (2017)
% Oyster toadfish (Opsanus tau) boatwhistle call detection 
% and patterns within a large-scale oyster restoration site
%
% ToadFishFinder v1 
% Del Bohnenstiehl  
% Jan 2023 



%% calculate PSD 
W=4096;  fs=24000; np2=13;  % zero pad to next  % zero pad to next 
y=gpuArray(y);              % uncomment for GPU 
[~,F,T,Pxx]=spectrogram(y,W,floor(W*.8),2^np2,fs);   %  generate spectrogram   
a=find(F >= ceil(frange(1)-3*s) & F <= ceil(2*frange(2)+3*s)); 
F=F(a); Pxx=Pxx(a,:);   % trim the frequency band requested  
[rowsF,~]=size(F);      % number of frequencies in the spectrogram 

%% determine the harmonics for the kernels based on F1range 
F1a=F(F> frange(1) & F < frange(2)); F1b=F1a-sweep; 
F1=[F1a'; F1b'];      % F1 2 x n matrix, with start and stop swept freqs of fo 
F2=F1*2;              % define 2nd harmionic 
[~,nkerns]=size(F1);  % determine the number of kernels 

% F1 is a 2 x n matrix, with start and stop 

%% adjust the data spectrogram 
PxxMod=10*log10(gather(Pxx)); % uncomment for GPU 
% PxxMod=10*log10((Pxx)); % comment for GPU
PxxModmed=median(PxxMod,1); % median of each column  
PxxMod=PxxMod-repmat(PxxModmed,size(PxxMod,1),1); % filter percussive signals. 
PxxMod=PxxMod-(smooth(median(PxxMod,2),31)*ones(1,length(T))); % additional background subtraction 

%% make a copy without percussive filtering for plotting 
if ploton==1 
PxxMod2=10*log10(gather(Pxx));  % uncomment for GPU 
% PxxMod2=10*log10((Pxx)); % comment for GPU
PxxMod2=PxxMod2-(smooth(median(PxxMod2,2),31)*ones(1,length(T)));  
end 

%% figure out size of the kernel 
klengthsec=0.225;
t=0:diff(T(1:2)):klengthsec;  d=max(t);  % t is time vector for kernel 
len_k=length(t);  % numnber of colums in kernel 

%% preallocate correlation matrix 
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
k2=nan(length(F),length(t)); 
for i=1:length(t) % for each time step
x=F-(fo+(t(i)./d)*(f1-fo));
k2(:,i)=(1-(x.^2)/(s.^2)).*exp(-(x.^2)/(2*s^2)); 
end

k=k1+k2;  % add the result together to make the final kernel 

knorm=sqrt(sum(k(:).^2)); 

%% now cross correlate with data spectrogram 
[c]=xcorr2(PxxMod,k);    
c=c(rowsF,:);  % pulls the center value that represents only a time shift 

%% norm to -1 to +1 
PxxModPad=padarray(PxxMod,[0, len_k-1],1,'both'); 
Pnorm=nan(1,length(cmatrix));
for jj=1:length(PxxModPad)-(len_k-1)
    Pnorm(jj)=sqrt(sum(reshape(PxxModPad(:,jj:jj+len_k-1),numel(k),1).^2));
end

cmatrix(j,:)=c./(Pnorm*knorm);   % stores it as a row for each kernel (or fund freq) 

end   % repeat for each of the nkerns 

cmatrix=cmatrix(:,len_k:end);    % trim to same size as T 

%% tapper the cmatrix based on prediction in frequency of boatwhistle 
if predffreq < frange(1); predffreq=frange(1); disp('resetting pred freq. to search range');  end 
if predffreq > frange(2); predffreq=frange(2); disp('resetting pred freq. to search range');  end 
if predffreq_uncert < 25; predffreq_uncert = 25; disp('uncertainity not realistic...resetting pred freq. uncert. to 25 Hz');  end 

lowcutoff = predffreq-predffreq_uncert;
highcutoff = predffreq+predffreq_uncert; 

xprd=1:20; yprd=exp(-xprd/15); xprd=xprd*2.94; % xprd/15 falls to 1/3 @ 50 Hz past cutoff corners. 
tapper4cmatrix=interp1([-12000, fliplr(lowcutoff-xprd),  lowcutoff, highcutoff, highcutoff+xprd  12000],[0.25,fliplr(yprd), 1, 1, yprd, 0.25],F1(1,:)); 
tapper4cmatrix=smooth(tapper4cmatrix,5);

figure; plot(F1(1,:),tapper4cmatrix)

cmatrix=cmatrix.*tapper4cmatrix;   
[cout,iout]=max(cmatrix);  % returns the max value (cout) & which row gave the max (iout)



[det_score,LOCS]=findpeaks(cout,'MinPeakHeight',thres,'MinPeakDistance',4,'MinPeakProminence',thres/3); 
det_time=gather(T(LOCS)'); 
det_time=gather(det_time); 
det_freq1=gather(F1(1,iout(LOCS)))';
det_freq2=gather(F2(1,iout(LOCS)))';
det_score=gather(det_score'); 


% %% if plotting is on 
if ploton==1 
figure; colormap('parula');
%
ax(1)=subplot(5,1,1); 
imagesc(T,F,PxxMod2); axis xy; hold on; caxis([-35,35]); 
ax(2)=subplot(5,1,2); 
imagesc(T,F,PxxMod); axis xy; hold on; caxis([-25,25]); 
plot(T(LOCS),det_freq1,'+k','MarkerSize',4,'LineWidth',1);
plot(T(LOCS),det_freq2,'+k','MarkerSize',4,'LineWidth',1); 
%
ax(3)=subplot(5,1,3);  plot(T,cout); 
hold on; plot(T(LOCS),det_score,'or'); grid on; ylim([0.1,0.6]); 
for j=1:length(LOCS)
text(gather(T(LOCS(j)))+0.15,gather(det_score(j)),sprintf('%3.1f',gather(det_score(j))),'FontSize',10); 
end
%
ax(4)=subplot(5,1,4); imagesc(T,F1(1,:),cmatrix);  axis xy; 
clim([-0.5,0.5]); 
%
ax(5)=subplot(5,1,5); hold on; 
for j=1:length(LOCS)
plot(gather(T(LOCS(j)))+cmatrix(:,LOCS(j)),F1(1,:));
end

linkaxes([ax(2),ax(1),ax(3),ax(4),ax(5)],'x'); 

end



