%%% testpitch.m
%%% BE3 SDI
%%% Stéphane Rossignol -- 2021

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Method 1 : Voicing+Pitch based on Analytic Signal (see BE1)

prename = "parenthese";  %%% trial on the file "parenthese.wav"
threshparam = [4 0.04 0.4];
[echtemps, voice1, f01] = voicingf0AS_BE(prename, threshparam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Method 2 : Voicing+Pitch based on Autocorrelation

fe = 8000; 
param = [0.98784626 0.4562832 0.0001777652 0.01163393 0.04 ...
         60 600 -0.7204987 6.8918675 249 0.7672432 ...
         0.98784626 0.4562832 0.0001777652 0.01163393];

pp1 = 2 * (param(6) - param(8)) / fe;  if (pp1 < 0), pp1 = 0; end
pp2 = 2 * (param(7) + param(9)) / fe;  if (pp2 > 1), pp2 = 1; end
aa = 1;
bb = fir2(param(10), ...
    [0 (param(6)-param(8))/fe*2 (param(7)+param(9))/fe*2 (param(7)+100)/fe*2 1], ...
    [0 1 param(11) 0 0]);
namefile = prename + '.wav';

[ttf0, voice, f0da] = pitcher_autoc(namefile, param, bb, aa);

% Interpolate results
[f02]   = interp1(ttf0, f0da, echtemps);
[voice2] = interp1(ttf0, voice, echtemps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Method 3 : Simple Pitch Tracking based on FFT

signal = audioread(namefile); 
[voice3, f03] = myPitchTrackingMethod(signal, fe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot pitch tracks

figure(1);  clf;
plot(echtemps, f01, 'r', 'DisplayName', 'Méthode analytique');
hold on;
plot(echtemps, f02, 'm', 'DisplayName', 'Autocorrélation');
plot(echtemps, f03, 'g', 'DisplayName', 'FFT (Méthode 3)');
plot(echtemps, signal / max(abs(signal)) * 300 * 0.9, 'b', ...
     'DisplayName', 'Signal (échelle normalisée)');
ylim([-300 300]);
title('Trajectoires de f_0 obtenues');
xlabel('Temps (s)');
ylabel('Fréquence fondamentale f_0 (Hz)');
legend('Location', 'best');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot voicing decisions

figure(2);  clf;
plot(echtemps, voice1, 'r', 'DisplayName', 'Voisement (analytique)');
hold on;
plot(echtemps, voice2 + 0.1, 'm', 'DisplayName', 'Voisement (autocorrélation)');
plot(echtemps, voice3 + 0.2, 'g', 'DisplayName', 'Voisement (FFT)');
plot(echtemps, signal / max(abs(signal)) * 0.9, 'b', ...
     'DisplayName', 'Signal (échelle normalisée)');
ylim([-1 1.2]);
title('Trajectoires de voisement obtenues');
xlabel('Temps (s)');
ylabel('Voisement (0: Non-voisé, 1: Voisé)');
legend('Location', 'best');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Combine pitch results

f0_final = combinePitchResults(f01, f02, f03);

figure(3);
plot(echtemps, f0_final, 'k');
title('Trajet final de f_0 combiné');
xlabel('Temps (s)');
ylabel('Fréquence fondamentale f_0 (Hz)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Optional: Plot combined pitch and normalized signal

figure(4); clf; hold on;
plot(echtemps, signal / max(abs(signal)) * 300 * 0.9, 'b', ...
    'DisplayName', 'Signal (échelle normalisée)');
plot(echtemps, f0_final, 'k', 'LineWidth', 2, ...
    'DisplayName', 'Trajet final f_0');
ylim([-300 300]);
title('Trajet final de f_0 superposé au signal');
xlabel('Temps (s)');
ylabel('Amplitude normalisée / f_0 (Hz)');
legend('Location','best');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions

function [voice3, f03] = myPitchTrackingMethod(audio, fs)
    frameSec = 0.04;  
    hopSec   = 0.01;  
    frameLen = round(frameSec * fs);
    hopLen   = round(hopSec   * fs);

    N = length(audio);
    numFrames = floor((N - frameLen) / hopLen) + 1;

    pitchPerFrame = zeros(numFrames,1);
    voicePerFrame = zeros(numFrames,1);
    timeFrames    = zeros(numFrames,1);

    w = hamming(frameLen);

    for k = 1:numFrames
        startIdx = (k-1)*hopLen + 1;
        endIdx   = startIdx + frameLen - 1;

        segment = audio(startIdx:endIdx).* w;
        X       = fft(segment);
        M       = abs(X(1:floor(end/2)+1));

        freqAxis = (0:length(M)-1)' * (fs/length(segment));

        [peakVal, peakIdx] = max(M);
        peakFreq = freqAxis(peakIdx);
        pitchPerFrame(k)    = peakFreq;
        voicePerFrame(k)    = (peakVal > 0.01);

        timeFrames(k) = (startIdx + endIdx)/(2*fs);
    end

    fullTime = (0:N-1)' / fs;
    f03      = interp1(timeFrames, pitchPerFrame, fullTime, 'nearest', 'extrap');
    voice3   = interp1(timeFrames, voicePerFrame, fullTime, 'nearest', 'extrap');
end

function f0_final = combinePitchResults(f01, f02, f03)
    % Make sure all are column vectors of the same length
    f01 = f01(:);
    f02 = f02(:);
    f03 = f03(:);
    
    % Simple average
    f0_final = (f01 + f02 + f03) / 3;

    % Remove outliers
    threshold = 2 * std(f0_final);
    meanVal   = mean(f0_final);
    mask      = abs(f0_final - meanVal) > threshold;
    f0_final(mask) = meanVal;
end
