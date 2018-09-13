%	OVERVIEW:
%       This basic demo will allow you to load an ECG file in matlab 
%       compatible wfdb format, detect the locations of the R peaks,
%       perform signal quality (SQI) analysis and plot the results.
%
%   OUTPUT:
%       A figure with the loaded ECG signal and detected peaks will be
%       generated
%
%   DEPENDENCIES & LIBRARIES:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   REFERENCE: 
%       Vest et al. "An Open Source Benchmarked HRV Toolbox for Cardiovascular 
%       Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Giulia Da Poian   
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
close all
% Where are the data, in this demo they are located in a subfolder
InputFolder = [pwd filesep 'TestData' filesep 'mitdb-Arrhythmia']; % path to the folder where you data are located
SigName = '200m';

% load the ecg signa using load (it loads a variable called val)
load([InputFolder filesep SigName]);
% the signal has two channels, from now on we will use just one 
ecg = val(1,:);
% Get sampling frequency Fs from header file
sigInfo = readheader2([InputFolder filesep SigName '.hea']);
Fs = sigInfo.freq;
% time vector for visualization (in seconds)
tm = 0:1/Fs:(length(ecg)-1)/Fs;

% plot the signal
%xsfigure(1)
plot(tm,ecg);
xlabel('[s]');
ylabel('[mV]')


% Detection of the R-peaks using the jqrs.m function included in the
% toolbox, requires to set initialization parameters calling the
% InitializeHRVparams.m function

HRVparams = InitializeHRVparams('Demo');
% set the exact sampling frequency usign the one from the loaded signal
HRVparams.Fs = Fs;
% call the function that perform peak detection
r_peaks = jqrs(ecg,HRVparams);

% plot the detected r_peaks on the top of the ecg signal
figure(1)
hold on;
plot(r_peaks./Fs, ecg(r_peaks),'o');
legend('ecg signal', 'detected R peaks')
hold off
%%
% The following parts of the code are added by Nasim Katebi on 09/12/2018 
%
%   OVERVIEW:
%       First part which is 'HRV using RR intervals' is for the detection 
%       of heart rate variability.
%       For this purpose heart rate will be estimated using a median of RR 
%       intervals in each segment (window)of signal.
%       In second part,'Line up the beats', beats of ecg will be aligned
%       based on detected QRS locations.
%       
% 
%   OUTPUT:
%         hr: vector of heart rate (bpm) size of vector is equal to
%         numberofbeats=round(length(ecg)/(WINDOW*Fs)) 
%         alignedbeat: This is a matrix for aligning the beats based on qrs
%         location

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  HRV Using RR intervals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculating heart rate in each window of a signal.
WINDOW=10; %window size in second
for i=1:round(length(ecg)/(WINDOW*Fs))
    start=(i-1)*WINDOW*Fs+1;
    stop=i*WINDOW*Fs;
    r=r_peaks(r_peaks>start & r_peaks<stop); % finding r peaks in eachwindow
    hr{1,i}=60/(median(diff(sort(r/Fs)))); % heart rate (bpm) in each window using median of RR intervals
end
figure(2)
plot(cell2mat(hr),'Linewidth',2)
xlabel('window index')
ylabel('Heart Rate (bpm) ')
title('Heart rate variability using ecg signal')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% line up the beats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extracting the beats by dividing the signal in a location of half of RR intervals
beatlocation=round((r_peaks(1:end-1)+r_peaks(2:end))/2);
%qrslocation is the location of r peak in each beat 
qrslocation=nan(1,(numel(beatlocation)-1));

for i=1:numel(beatlocation)-1
beat{1,i}=(ecg(beatlocation(i):beatlocation(i+1))); %saving beats in a cell
qrslocation(1,i)=r_peaks(i+1)-beatlocation(i);
end

beatlength=beatlocation(2:end)-beatlocation(1:end-1);
lengthofarray=2*max(beatlength); 
alignedbeat=nan(length(beat),round(lengthofarray)); % alignedarray is a (number of beats)x(lengthofarray) dimension matrix which we created to align the beats based on qrs locations
qrs=lengthofarray/2; % beats will be aligned based on this location for qrs       
for i=1:length(beat)
  
    alignedbeat(i,((qrs-qrslocation(i)+1):qrs))=beat{1,i}(1:qrslocation(i));
    alignedbeat(i,((qrs:qrs+(length(beat{1,i})-qrslocation(i)))))=beat{1,i}(qrslocation(i):end);
    
end
figure(3)
hold on
for i=1:length(beat)
    signal=alignedbeat(i,:);
     plot(signal)  
end
title('lined up beats')
xlabel('sample')
ylabel('mV')






