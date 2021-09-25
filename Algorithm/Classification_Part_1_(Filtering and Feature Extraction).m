%%% This is the algorithm to calculate the features to create the dataset for the classification %%%


%% Reading collected data 

clc, clear all, close all
 
files = dir('*.mat')'; %read all the collected .mat files
L = length(files); 
for d = 1:L
    FZ(d) = load(files(d).name); %load files in a struct
end
F = struct2cell(FZ); %convert struct to cell 
FM = cell2mat(F); %convert cell to matrix 
sz = size(FM); %size of the matrix 

% The columns represents:  gyroscope signals, acelerometer signals, compass vectors in x, y, z axis


%% Time vector

fs = 20; %sampling rate [Hz] 
N = sz(1); % number of samples per signal 
S1 = sz(3); %number of signals
t = [0:N-1]./fs; %time axis in seconds (0-114.85 s)
t = t';

%% In this part of the algorithm the signals are divided in 6 matrices of interest (gyroscope signals and acelerometer signals in x, y, z axis) 

for k = 1:sz(3)*sz(2)
    AA(:,k) = FM(:,k);
end

GRX(:,1) = AA(:,1);
GRY(:,1) = AA(:,2);
GRZ(:,1) = AA(:,3);
ACX(:,1) = AA(:,4);
ACY(:,1) = AA(:,5);
ACZ(:,1) = AA(:,6);

for q = 1:sz(3)-1
        GRX(:,q+1) = AA(:,10+(q-1)*9);
        GRY(:,q+1) = AA(:,11+(q-1)*9);
        GRZ(:,q+1) = AA(:,12+(q-1)*9);
        ACX(:,q+1) = AA(:,13+(q-1)*9);
        ACY(:,q+1) = AA(:,14+(q-1)*9);
        ACZ(:,q+1) = AA(:,15+(q-1)*9);
end

%% Filtering 

% A 5th order butterworth passband filter (0.5-7.5 Hz) was used to reduce the noise and accomodate the disorders frequency bands:
% Parkinsonian tremor (4-6 Hz) 
% Essential Tremor (5-7 Hz)
% Diskynesia (1-3 Hz)

order = 5;
[b1,a1]=butter(order,[0.5,7.5]/(fs/2),'bandpass');  
for k = 1:S1
    GRX(:,k) = filter(b1,a1,GRX(:,k)');
    GRY(:,k) = filter(b1,a1,GRY(:,k)');
    GRZ(:,k) = filter(b1,a1,GRZ(:,k)');
    ACX(:,k) = filter(b1,a1,ACX(:,k)');
    ACY(:,k) = filter(b1,a1,ACY(:,k)');
    ACZ(:,k) = filter(b1,a1,ACZ(:,k)'); 
end


%%  Calculation of features

% Power of signals in time domain
for k = 1:S1
    aPowX1(k) = sum(abs(ACX(:,k)).^2)/(2*N+1);
    aPowY1(k) = sum(abs(ACY(:,k)).^2)/(2*N+1);
    aPowZ1(k) = sum(abs(ACZ(:,k)).^2)/(2*N+1);
    gPowX1(k) = sum(abs(GRX(:,k)).^2)/(2*N+1);
    gPowY1(k) = sum(abs(GRY(:,k)).^2)/(2*N+1);
    gPowZ1(k) = sum(abs(GRZ(:,k)).^2)/(2*N+1);
end

% frequency vector 
freq = fs*[0:N-1]./N;
freq = freq';
ind_freq = find(freq >= 4);

%fast fourier transform of the signals
aFT_X = fft(ACX);
aFT_Y = fft(ACY);
aFT_Z = fft(ACZ);
gFT_X = fft(GRX);
gFT_Y = fft(GRY);
gFT_Z = fft(GRZ);


%standard deviation of the signals
astd_X = std(ACX);
astd_Y = std(ACY);
astd_Z = std(ACZ);
gstd_X = std(GRX);
gstd_Y = std(GRY);
gstd_Z = std(GRZ);

%mean of the signals
amean_X = mean(ACX);
amean_Y = mean(ACY);
amean_Z = mean(ACZ);
gmean_X = mean(GRX);
gmean_Y = mean(GRY);
gmean_Z = mean(GRZ);


%Amplitude of the signals
for k = 1:S1 
    aAmp_X(k) = max(ACX(:,k))-min(ACX(:,k));
    aAmp_Y(k) = max(ACY(:,k))-min(ACY(:,k));
    aAmp_Z(k) = max(ACZ(:,k))-min(ACZ(:,k));
    gAmp_X(k) = max(GRX(:,k))-min(GRX(:,k));
    gAmp_Y(k) = max(GRY(:,k))-min(GRY(:,k));
    gAmp_Z(k) = max(GRZ(:,k))-min(GRZ(:,k));
end

Power of the signals in frequency domain
for k = 1:S1
    afPowX(k) = sum(abs(aFT_X(:,k)).^2)/(2*N+1);
    afPowY(k) = sum(abs(aFT_Y(:,k)).^2)/(2*N+1);
    afPowZ(k) = sum(abs(aFT_Z(:,k)).^2)/(2*N+1);
    gfPowX(k) = sum(abs(gFT_X(:,k)).^2)/(2*N+1);
    gfPowY(k) = sum(abs(gFT_Y(:,k)).^2)/(2*N+1);
    gfPowZ(k) = sum(abs(gFT_Z(:,k)).^2)/(2*N+1);
end

for k = 1:S1
    for m = 1:length(ind_freq)
        aFT_X_40(m,k) = aFT_X(m+ind_freq(1)-1,k);
        aFT_Y_40(m,k) = aFT_Y(m+ind_freq(1)-1,k);
        aFT_Z_40(m,k) = aFT_Z(m+ind_freq(1)-1,k);
        gFT_X_40(m,k) = gFT_X(m+ind_freq(1)-1,k);
        gFT_Y_40(m,k) = gFT_Y(m+ind_freq(1)-1,k);
        gFT_Z_40(m,k) = gFT_Z(m+ind_freq(1)-1,k);        
    end
end
for k = 1:S1
    afPowX_40(k) = sum(abs(aFT_X_40(:,k)).^2)/(2*N+1);
    afPowY_40(k) = sum(abs(aFT_Y_40(:,k)).^2)/(2*N+1);
    afPowZ_40(k) = sum(abs(aFT_Z_40(:,k)).^2)/(2*N+1);
    gfPowX_40(k) = sum(abs(gFT_X_40(:,k)).^2)/(2*N+1);
    gfPowY_40(k) = sum(abs(gFT_Y_40(:,k)).^2)/(2*N+1);
    gfPowZ_40(k) = sum(abs(gFT_Z_40(:,k)).^2)/(2*N+1);
end

%Percentage of power of the signals above 4 Hz

for k = 1:S1
    afPow_Perc4X(k) = afPowX_40(k)*100/afPowX(k);
    afPow_Perc4Y(k) = afPowY_40(k)*100/afPowY(k);
    afPow_Perc4Z(k) = afPowZ_40(k)*100/afPowZ(k);
    gfPow_Perc4X(k) = gfPowX_40(k)*100/gfPowX(k);
    gfPow_Perc4Y(k) = gfPowY_40(k)*100/gfPowY(k);
    gfPow_Perc4Z(k) = gfPowZ_40(k)*100/gfPowZ(k);
end


%% Correlation
for mm = 1:S1
    [aCRXX(:,mm), FL] = xcorr(ACX(:,mm),'biased');
    [aCRYY(:,mm), FL] = xcorr(ACY(:,mm),'biased');
    [aCRZZ(:,mm), FL] = xcorr(ACZ(:,mm),'biased');
    [aCRXY(:,mm), FL] = xcorr(ACX(:,mm),ACY(:,mm),'biased');
    [aCRXZ(:,mm), FL] = xcorr(ACX(:,mm),ACZ(:,mm),'biased');
    [aCRYZ(:,mm), FL] = xcorr(ACY(:,mm),ACZ(:,mm),'biased');
    [gCRXX(:,mm), FL] = xcorr(GRX(:,mm),'biased');
    [gCRYY(:,mm), FL] = xcorr(GRY(:,mm),'biased');
    [gCRZZ(:,mm), FL] = xcorr(GRZ(:,mm),'biased');
    [gCRXY(:,mm), FL] = xcorr(GRX(:,mm),GRY(:,mm),'biased');
    [gCRXZ(:,mm), FL] = xcorr(GRX(:,mm),GRZ(:,mm),'biased');
    [gCRYZ(:,mm), FL] = xcorr(GRY(:,mm),GRZ(:,mm),'biased');    
    %[CRA(:,mm), FL] = xcorr(A_rms(:,mm),'biased');
end   

t_corr = [-length(FL)/2:(length(FL)/2)-1]./fs;

%cross-correlation peak
for k = 1:S1 
    [aAmp_CRXY(k),aInd_CRXY(k)] = max(aCRXY(:,k));
    [aAmp_CRXZ(k),aInd_CRXZ(k)] = max(aCRXZ(:,k));
    [aAmp_CRYZ(k),aInd_CRYZ(k)] = max(aCRYZ(:,k));
    [gAmp_CRXY(k),gInd_CRXY(k)] = max(gCRXY(:,k));
    [gAmp_CRXZ(k),gInd_CRXZ(k)] = max(gCRXZ(:,k));
    [gAmp_CRYZ(k),gInd_CRYZ(k)] = max(gCRYZ(:,k));
end


for k = 1:S1
    a_tcorr_XY(k) = abs(t_corr(aInd_CRXY(k)));
    a_tcorr_XZ(k) = abs(t_corr(aInd_CRXZ(k)));
    a_tcorr_YZ(k) = abs(t_corr(aInd_CRYZ(k)));
    g_tcorr_XY(k) = abs(t_corr(gInd_CRXY(k)));
    g_tcorr_XZ(k) = abs(t_corr(gInd_CRXZ(k)));
    g_tcorr_YZ(k) = abs(t_corr(gInd_CRYZ(k)));
end


% autocorrelation peak
for k = 1:S1 
    [aAmp_CRXX(k),aInd_CRXX(k)] = max(aCRXX(:,k));
    [aAmp_CRYY(k),aInd_CRYY(k)] = max(aCRYY(:,k));
    [aAmp_CRZZ(k),aInd_CRZZ(k)] = max(aCRZZ(:,k));
    [gAmp_CRXX(k),gInd_CRXX(k)] = max(gCRXX(:,k));
    [gAmp_CRYY(k),gInd_CRYY(k)] = max(gCRYY(:,k));
    [gAmp_CRZZ(k),gInd_CRZZ(k)] = max(gCRZZ(:,k));
end

% autocorrelation peaks and lags
for k = 1:S1
    [aCRpkXX(1:length(findpeaks(aCRXX(:,k))),k)  aCRlagXX(1:length(findpeaks(aCRXX(:,k))),k)] = findpeaks(aCRXX(:,k));
    [aCRpkYY(1:length(findpeaks(aCRYY(:,k))),k)  aCRlagYY(1:length(findpeaks(aCRYY(:,k))),k)] = findpeaks(aCRYY(:,k));
    [aCRpkZZ(1:length(findpeaks(aCRZZ(:,k))),k)  aCRlagZZ(1:length(findpeaks(aCRZZ(:,k))),k)] = findpeaks(aCRZZ(:,k));

    [gCRpkXX(1:length(findpeaks(gCRXX(:,k))),k)  gCRlagXX(1:length(findpeaks(gCRXX(:,k))),k)] = findpeaks(gCRXX(:,k));
    [gCRpkYY(1:length(findpeaks(gCRYY(:,k))),k)  gCRlagYY(1:length(findpeaks(gCRYY(:,k))),k)] = findpeaks(gCRYY(:,k));
    [gCRpkZZ(1:length(findpeaks(gCRZZ(:,k))),k)  gCRlagZZ(1:length(findpeaks(gCRZZ(:,k))),k)] = findpeaks(gCRZZ(:,k));
end

% autcorrelation peak of the first lag
for k = 1:S1
    a1stlagXX(k) = aCRlagXX(1,k);
    a1stlagYY(k) = aCRlagYY(1,k);
    a1stlagZZ(k) = aCRlagZZ(1,k);
    g1stlagXX(k) = gCRlagXX(1,k);
    g1stlagYY(k) = gCRlagYY(1,k);
    g1stlagZZ(k) = gCRlagZZ(1,k);    
end

% amplitude of the first autcorrelation peak 
for k = 1:S1
    a1stpkXX(k) = aCRpkXX(1,k);
    a1stpkYY(k) = aCRpkYY(1,k);
    a1stpkZZ(k) = aCRpkZZ(1,k);
    g1stpkXX(k) = gCRpkXX(1,k);
    g1stpkYY(k) = gCRpkYY(1,k);
    g1stpkZZ(k) = gCRpkZZ(1,k);    
end

%sum of the peaks of autocorrelation
for k = 1:S1
    aSumpkXX(k) = sum(aCRpkXX(:,k));
    aSumpkYY(k) = sum(aCRpkYY(:,k));
    aSumpkZZ(k) = sum(aCRpkZZ(:,k));
    gSumpkXX(k) = sum(gCRpkXX(:,k));
    gSumpkYY(k) = sum(gCRpkYY(:,k));
    gSumpkZZ(k) = sum(gCRpkZZ(:,k));    
end

%number of autocorrelation peaks
for k = 1:S1
    aNpkXX(:,k) = sum((aCRpkXX(:,k)~= 0));
end

% frequency vector for the spectrum calculation 
T = FL/fs;
r = 1/(T(end)-T(1));
f = [0:length(T)-1]*r;
f=f';

% Spectrum
aSpXX = abs(fft(aCRXX));
aSpYY = abs(fft(aCRYY));
aSpZZ = abs(fft(aCRZZ));
gSpXX = abs(fft(gCRXX));
gSpYY = abs(fft(gCRYY));
gSpZZ = abs(fft(gCRZZ));

%% Cross-Spectrum
aSpXY = abs(fft(aCRXY));
aSpXZ = abs(fft(aCRXZ));
aSpYZ = abs(fft(aCRYZ));
gSpXX = abs(fft(gCRXX));
gSpXZ = abs(fft(gCRXZ));
gSpYZ = abs(fft(gCRYZ));

%% 
% dominant frequencies and PSD (Power Spectrum Density) of dominant frequency
for k = 1:S1 
    [aAmp_SpX(k),aInd_SpX(k)] = max(aSpXX(:,k));
    [aAmp_SpY(k),aInd_SpY(k)] = max(aSpYY(:,k));
    [aAmp_SpZ(k),aInd_SpZ(k)] = max(aSpZZ(:,k));
    [gAmp_SpX(k),gInd_SpX(k)] = max(gSpXX(:,k));
    [gAmp_SpY(k),gInd_SpY(k)] = max(gSpYY(:,k));
    [gAmp_SpZ(k),gInd_SpZ(k)] = max(gSpZZ(:,k));   

    af_dom1X(k) = f(aInd_SpX(k));
    af_dom1Y(k) = f(aInd_SpY(k)); 
    af_dom1Z(k) = f(aInd_SpZ(k));
    gf_dom1X(k) = f(gInd_SpX(k)); 
    gf_dom1Y(k) = f(gInd_SpY(k));
    gf_dom1Z(k) = f(gInd_SpZ(k));
end

% second dominant frequencies and PSD (Power Spectrum Density) of second dominant frequency 
for k = 1:S1
    aSp_ordX(:,k) = sort(aSpXX(:,k));
    aSp_ordY(:,k) = sort(aSpYY(:,k));
    aSp_ordZ(:,k) = sort(aSpZZ(:,k));
    gSp_ordX(:,k) = sort(gSpXX(:,k));
    gSp_ordY(:,k) = sort(gSpYY(:,k));
    gSp_ordZ(:,k) = sort(gSpZZ(:,k));
    
    aSp_2ndX(:,k) = aSp_ordX(length(aSpXX)-3,k);
    aSp_2ndY(:,k) = aSp_ordY(length(aSpYY)-3,k); 
    aSp_2ndZ(:,k) = aSp_ordZ(length(aSpZZ)-3,k); 
    gSp_2ndX(:,k) = gSp_ordX(length(gSpXX)-3,k); 
    gSp_2ndY(:,k) = gSp_ordY(length(gSpYY)-3,k); 
    gSp_2ndZ(:,k) = gSp_ordZ(length(gSpZZ)-3,k);     
    
    afindex_X(:,k) = find(aSpXX(:,k) == aSp_2ndX(:,k));
    afindex_Y(:,k) = find(aSpYY(:,k) == aSp_2ndY(:,k));
    afindex_Z(:,k) = find(aSpZZ(:,k) == aSp_2ndZ(:,k));
    gfindex_X(:,k) = find(gSpXX(:,k) == gSp_2ndX(:,k));
    gfindex_Y(:,k) = find(gSpYY(:,k) == gSp_2ndY(:,k));
    gfindex_Z(:,k) = find(gSpZZ(:,k) == gSp_2ndZ(:,k));    
    
    af_dom2X(:,k) =  f(afindex_X(1,k)); 
    af_dom2Y(:,k) =  f(afindex_Y(1,k)); 
    af_dom2Z(:,k) =  f(afindex_Z(1,k)); 
    gf_dom2X(:,k) =  f(gfindex_X(1,k));
    gf_dom2Y(:,k) =  f(gfindex_Y(1,k)); 
    gf_dom2Z(:,k) =  f(gfindex_Z(1,k)); 
    
end

%% plots

% Filtered signals in time domain 

% Parkinsonian tremor example
figure
subplot(3,2,1)
plot(t,GRX(:,17))
xlabel('tempo(s)')
ylabel('amplitude(º/s)')
ylim([-5, 5])
xlim([0 115])
title('GRX: TP')
subplot(3,2,3)
plot(t,GRY(:,17))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(º/s)')
title('GRY: TP')
subplot(3,2,5)
plot(t,GRZ(:,17))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(º/s)')
title('GRZ: TP')
subplot(3,2,2)
plot(t,ACX(:,17))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(g)')
title('ACX: TP')
subplot(3,2,4)
plot(t,ACY(:,17))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(g)')
title('ACY: TP')
subplot(3,2,6)
plot(t,ACZ(:,17))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(g)')
title('ACZ: TP')

%Essential Tremor example
figure
subplot(3,2,1)
plot(t,GRX(:,10))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(º/s)')
title('GRX: TE')
subplot(3,2,3)
plot(t,GRY(:,10))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(º/s)')
title('GRY: TE')
subplot(3,2,5)
plot(t,GRZ(:,10))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(º/s)')
title('GRZ: TE')
subplot(3,2,2)
plot(t,ACX(:,10))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(g)')
title('ACX: TE')
subplot(3,2,4)
plot(t,ACY(:,10))
ylim([-5, 5])
xlabel('tempo(s)')
ylabel('amplitude(g)')
xlim([0 115])
title('ACY: TE')
subplot(3,2,6)
plot(t,ACZ(:,10))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(g)')
title('ACZ: TE')

% Dyskinesia example
figure
subplot(3,2,1)
plot(t,GRX(:,55))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(º/s)')
title('GRX: Discinesia')
subplot(3,2,3)
plot(t,GRY(:,55))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(º/s)')
title('GRY: Discinesia')
subplot(3,2,5)
plot(t,GRZ(:,55))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(º/s)')
title('GRZ: Discinesia')
subplot(3,2,2)
plot(t,ACX(:,55))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(g)')
title('ACX: Discinesia')
subplot(3,2,4)
plot(t,ACY(:,55))
ylim([-5, 5])
xlabel('tempo(s)')
ylabel('amplitude(g)')
xlim([0 115])
title('ACY: Discinesia')
subplot(3,2,6)
plot(t,ACZ(:,55))
ylim([-5, 5])
xlim([0 115])
xlabel('tempo(s)')
ylabel('amplitude(g)')
title('ACZ: Discinesia')

% Power Spectrum Density 

figure
[gpsd_X f_gx] = pwelch(GRX,1000,300,1000,20);
[gpsd_Y f_gy] = pwelch(GRY,1000,300,1000,20);
[gpsd_Z f_gz] = pwelch(GRZ,1000,300,1000,20);
[apsd_X f_ax] = pwelch(ACX,1000,300,1000,20);
[apsd_Y f_ay] = pwelch(ACY,1000,300,1000,20);
[apsd_Z f_az] = pwelch(ACZ,1000,300,1000,20);


% Parkinsonian tremor example
subplot(3,2,1)
plot(f_gx,gpsd_X(:,17))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('GRX: TP')
subplot(3,2,3)
plot(f_gx,gpsd_Y(:,17))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('GRY: TP')
subplot(3,2,5)
plot(f_gx,gpsd_Z(:,17))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('GRZ: TP')
subplot(3,2,2)
plot(f_ax,apsd_X(:,17))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('ACX: TP')
subplot(3,2,4)
plot(f_ay,apsd_Y(:,17))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('ACY: TP')
subplot(3,2,6)
plot(f_az,apsd_Z(:,17))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('ACZ: TP')


%Essential Tremor example
figure
subplot(3,2,1)
plot(f_gx,gpsd_X(:,10))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('GRX: TE')
subplot(3,2,3)
plot(f_gx,gpsd_Y(:,10))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('GRY: TE')
subplot(3,2,5)
plot(f_gx,gpsd_Z(:,10))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('GRZ: TE')
subplot(3,2,2)
plot(f_ax,apsd_X(:,10))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('ACX: TE')
subplot(3,2,4)
plot(f_ay,apsd_Y(:,10))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('ACY: TE')
subplot(3,2,6)
plot(f_az,apsd_Z(:,10))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('ACZ: TE')


% Dyskinesia example

figure
subplot(3,2,1)
plot(f_gx,gpsd_X(:,55))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('GRX: Discinesia')
subplot(3,2,3)
plot(f_gx,gpsd_Y(:,55))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('GRY: Discinesia')
subplot(3,2,5)
plot(f_gx,gpsd_Z(:,55))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('GRZ: Discinesia')
subplot(3,2,2)
plot(f_ax,apsd_X(:,55))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('ACX: Discinesia')
subplot(3,2,4)
plot(f_ay,apsd_Y(:,55))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('ACY: Discinesia')
subplot(3,2,6)
plot(f_az,apsd_Z(:,55))
ylim([0, 0.5])
xlim([0 10])
xlabel('frequência (Hz)')
ylabel('PSD (W/Hz)')
title('ACZ: Discinesia')


