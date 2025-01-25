%%% Silence detection
%%% 
%%% St\'ephane Rossignol -- 30/06/10

%%% you can use this code as a matlab/octave SCRIPT or as matlab/octave FUNCTION
%%% in order to change, comment/uncomment the lines just below

%%% STANDALONE
clear all;
close all;
[xxx,fe] = audioread('parenthese.wav');  % Fichier modifié pour 'parenthese.wav'

%%% FUNCTION
%function [signal] = sil_detect(xxx, fe)

ttt = [1:length(xxx)] / fe;

%%% parameters -- beginning
tsig = round(0.04 * fe);     %%% frame size
tsil = round(1.00 * fe);     %%% minimum length of a true silence
tvoice = round(0.04 * fe);   %%% minimum length of a voice activity segment
fen = blackman(tsig);        %%% ponderation window
thresh1 = 0.04;              %%% threshold value (re-estimated below, using 'lestim')
lestim = 2000;               %%% the 'lestim' first samples of the sound are supposed to be silence
                             %%% and are used to re-estimate the threshold
p1 = 1.2;
p2 = 0.99;                   %%% parameters used for the re-estimation of the threshold
%%% parameters -- end

nn = max(xxx);
figure(3);
clf;
plot(ttt, xxx);
hold on;
plot([1 length(xxx)] / fe, [thresh1 thresh1], 'r');
activity = 0; %%% activity not yet detected
sil = 1;      %%% the sound starts with a silence
fen = fen / sum(fen) * tsig;
actualstt = -1;

%%% threshold re-estimation -- my method
ii = tsig;
jj = 1;
while ii <= lestim
    yyy = xxx(ii - tsig + 1:ii) .* fen;
    ene(jj) = sqrt(sum(yyy .* yyy) / tsig);

    ii = ii + 1;
    jj = jj + 1;
end
enei = sort(ene);
thresh1 = p1 * enei(round(length(ene) * p2));

%%% threshold re-estimation -- your method
%thresh1 = ..........

%%% main loop
ii = tsig;
jj = 1;
ene = zeros(1, length(xxx) - tsig);
while ii <= length(xxx)
    yyy = xxx(ii - tsig + 1:ii) .* fen;
    ene(jj) = sqrt(sum(yyy .* yyy) / tsig);

    %%% voice activity detection
    if (ene(jj) > thresh1 && sil == 1)
        sil = 0;
        stt = ii;
    end
    if (ene(jj) < thresh1 && sil == 0)
        sil = 1;
        edd = ii;
        if (edd - stt < tvoice)
            sil = 0;   %%% too short voice activity segment; removed
        else
            if (actualstt < 0)
                actualstt = stt;
            end
            actualedd = edd;

            fprintf(1, 'voice activity length: %f (%d %d)\n', (edd - stt) / fe, stt, edd);

            activity = 1;  %%% the activity has been detected
            sttsil = edd;  %%% starting instant of the final silence

            plot([stt - round(tsig / 2) edd - round(tsig / 2)] / fe, 0.1 * [nn nn], 'k');
            drawnow;
        end
    end
    if (sil == 1 && activity == 1)
        if (ii - sttsil > tsil)
            fprintf(1, 'Processing stopped after a long silence\n');
            plot([ii ii] / fe, [min(xxx) max(xxx)]);
            ii = length(xxx);
        end
    end

    ii = ii + 1;
    jj = jj + 1;
end

actualstt = actualstt - round(tsig / 2);
actualedd = actualedd - round(tsig / 2);

actualstt = actualstt - round(0.1 * fe); %%% enlarging of the segment
actualedd = actualedd + round(0.1 * fe); %%% enlarging of the segment
plot([actualstt actualedd] / fe, 0.2 * [nn nn], 'm');  %%% line segment indicating which part of the signal is used

% Adjust energy padding to match ttt length
ene_padded = [zeros(1, round(tsig / 2)) ene zeros(1, length(ttt) - length(ene) - round(tsig / 2))];
plot(ttt, ene_padded, 'g', 'LineWidth', 2);

ylim([min(xxx) max(xxx)]);
hold off;
print -depsc2 detesil.eps

%%% the useful signal is saved in another variable
%%% => the ASR used afterwards should try to recognize only this useful part of the signal
signal = xxx(actualstt:actualedd);
