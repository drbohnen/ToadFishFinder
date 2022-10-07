clear;
% close all
%
% This is a wrapper function for the NCSU cultch reef dataset that 
% the tfDeepToadDetector.m function 
% 
load STcalibration.mat          % calibration data 
                                % loads a variable called STcalib 
load PScultch1_dir2process.mat; % directory list for Pamlico Sound Cultch 
                                % loads a variable dir2process 
load tfclassifier_v1.mat     % trained classifier 
                                % loads a classification objected 
                                % named classifier 
% this is where you want the outou 
DirOut='/Volumes/G6/d_CultchTimeSeries/TF/toadfishdet/'; 

%% wrapper works by looping through the deployments in dir2process 

for H=1:3%:height(dir2process)   % choose which directories (i.e., site, deployments) 
DirIn= char(dir2process.DirIn(H)); % directory in
SiteOut=char(dir2process.Site(H)); 
DepOut=dir2process.Deployment(H); 
[filelist, fstart, fend] = mktableSTdir(char(DirIn)); % list of file in that directory 
a=find(fstart > datenum(dir2process.Sgate(H)) & fstart < datenum(dir2process.Egate(H))); % filter to time range specified 
filelist=filelist(a); fstart=fstart(a); fend=fend(a); % reset based on filtered time range 
site=char(dir2process.Site(H)); dep=dir2process.Deployment(H);  % define variable site and dep 

% The second deployment was 1 minutes recordings 
if dep==2 % for NC Cultch data, dep 2 was only recorded for 60s
    NSEC=60;
else
    NSEC=120;
end


%% Load File and Run perch detector on each file 
maxfiles=length(filelist); % 
Bcount=nan(maxfiles,1); Ocount=nan(maxfiles,1);  % reallocate 
for i=1:3 %maxfiles  
    fprintf('Processing %s\n', char(filelist(i).name));
    [y,fstart_UTC, fs, metadata]=readST(char(filelist(i).name),char(DirIn),NSEC,STcalib);   % DirIn cell, char changes to character for fileread function
    y=resample(y,1,4);  
% run the perch detector 
[~, ~,Bcount(i),Ocount(i)]=tfDeepToadDetector(y,char(dir2process.Site(H)),char(filelist(i).name), fullfile(DirOut,SiteOut,sprintf('%02.0f',DepOut)),classifier,'other',datetime(fstart_UTC,'ConvertFrom','datenum','Format','dd-MMM-uuuu HH:mm:ss.SSSSSSSS')); 
end


%% Make a big table wiht the time of all the boatwhistle detections 
F=dir(fullfile(DirOut,SiteOut,sprintf('%02.0f',DepOut),'*','Btable*.mat')); 
 clear('Btable*')

 if isempty(F)==0
for f=1:length(F)
    if f==1 
    load(fullfile(F(f).folder,F(f).name))
    V=whos('Btable*'); 
    eval([site '_' sprintf('%02.0f',dep) '_DetTable=' V.name ';']); 
    clear('Btable*')
    elseif f>1  
    load(fullfile(F(f).folder,F(f).name))
    V=whos('Btable*'); 
    eval([site '_' sprintf('%02.0f',dep) '_DetTable=[' site '_' sprintf('%02.0f',dep) '_DetTable; ' V.name '];']); 
    clear('Btable*')
    end
end

 end

%% File naming and saving 
filenamesout=char(filelist(1:length(Bcount)).name); 
filetimesout=datetime(fstart(1:length(Bcount)),'ConvertFrom','datenum'); 
eval([site '_' sprintf('%02.0f',dep) '_DetTab=table(Bcount, Ocount, filenamesout, filetimesout,''VariableNames'',{''perch'',''other'',''file'',''time''})'] ); 

 if isempty(F)==0
out1=fullfile(DirOut,[ site '_' sprintf('%02.0f',dep) '_DetTable.mat']); 
eval(['save ' out1 ' ' site '_' sprintf('%02.0f',dep) '_DetTable']); 
 end

out2=fullfile(DirOut,[ site '_' sprintf('%02.0f',dep) '_DetTab.mat']); 
eval(['save ' out2 ' ' site '_' sprintf('%02.0f',dep) '_DetTab']); 

%
clear Btable Bcount Ocount F 

end
 

