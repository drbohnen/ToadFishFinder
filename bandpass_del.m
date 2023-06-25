 function [d]=bandpass_del(c,flp,fhi,Fs,n) 
 % 
 % Function bandpass applies nth order butterworth filter 
 % [d]=bandpass(c,flp,fhi,Fs,n) 
 % 
 % INPUT 
 % c = input time series 
 % flp = lowpass corner frequency of filter 
 % fhi = hipass corner frequency 
 % Fs = sample rate 
 % n = filter order
 % 
 % OUTPUT 
 % d is the bandpassed waveform. 
 % 
 % Del Bohnenstiehl - NCSU 
 % drbohnen@ncsu.edu 
 
 if isempty(n) 
     n=2; 
 end
 fnq=0.5*Fs;              % Nyquist frequency 
 Wn=[flp/fnq fhi/fnq];    % non-dimensional frequency 
 [b,a]=butter(n,Wn);      % construct the filter 
 d=filtfilt(b,a,c);       % zero phase filter the data 
 return;
