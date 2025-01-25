%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Voicing coefficient extraction -- Method: Analytic Signal
% f0 analysis added (see 'BE3')
%
% - prename        : file name without the .wav extension
% - threshparam(1) : thresholding method (Otsu=4)
% - threshparam(2) : length of the analysis window (in second)
% - threshparam(3) : thresholding parameter (in [0 1]; must be trained)
%
% => note: the sampling rate must be 16 kHz
%
% St\'ephane Rossignol                       -- 22/01/2014
%                                              (considering this version)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% you can use this code as a matlab/octave SCRIPT or as matlab/octave FUNCTION
%%% in order to change, comment/uncomment the lines just below

%%% standalone -- octave script
%clear all;close all;
%prename="parenthese";  %%% trial on the file "parenthese.wav"
%threshparam=[4 0.04 0.4];
%%% function -- octave function
function [echtemps,segmeAS,freqinst10] = voicingf0AS_BE(prename, threshparam)


%%% required by latest versions of "octave"
%pkg load signal



%%%%%% for BE1 -- voicing coefficient extraction



%%% parameters
threshmeth = threshparam(1); %%% thresholding method
KK         = 2.0;            %%% band-pass filter decay
hlength    = threshparam(2); %%% length of the analysis window (in second)
NLt        = 1.00*hlength;   %%% smoothing length (in seconds)


%%% for debugging purposes, some plots and comments are made
plt=0;


%%% some initializations
extsnrwin = [num2str(100) '-' num2str(round(hlength*1000))]; %%% where is saved the result

if ( exist(prename+extsnrwin+'-AS.txt','file')~=2 )
  fprintf(1,"1-");
end;


%%% the sound is read (note: fs=16000 is requested)
[signal, fs] = audioread(prename+'.wav');
lsig = length(signal);


%%% low-pass filtering + decimation
fcoupure=600; %%% frequency cutoff
lpfd = 2;     %%% if you time CPU time, not necessary to use 'lpfd = 1'
switch lpfd
  case 1     %%% both low-pass filtering and decimation are performed
    ndec = 5;                 %%% decimation factor
    N5   = round(hlength*fs); %%% filter length (in samples)
    h5   = fir1(N5, fcoupure/(fs/2));
    signal5 = filter(h5,1,signal);
    signal5 = signal5(1:ndec:length(signal5));

  case 2     %%% only low-pass filtering is performed
    N5 = round(hlength*fs); %%% filter length (in samples)
    h5 = fir1(N5, fcoupure/(fs/2));
    signal5 = filter(h5,1,signal);
    ndec = 1;

  otherwise  %%% nothing is performed
    N5   = 0;
    ndec = 1;      %%% there is no decimation
    signal5 = signal;
end;


%%% Hilbert filtering

%% initializations
N2t = hlength;        %%% filter length (in s)
N2  = round(N2t*fs);  %%% filter length (in samples)
N2  = N2/ndec;
win = blackman(N2+1)';


%% filter prototype
sf1 = 100.0;
sf2 = 310.0;

fint = (sf2-sf1)/6;
BU = [sf1  sf1+fint  sf1+2*fint  sf1+3*fint  sf1+4*fint  sf1+5*fint  sf2];   %%% useful bandwith (in Hz)

FF = [0.0   sf1-30   BU                   sf2+50   1000   fs/2]/(fs/2); %%% frequencies; 1=fs/2
MM = [0.0   1.0      (sf1^KK)./(BU.^KK)   0        0      0];           %%% magnitudes

%%
FF = FF/ndec;
reducf=2*(0:N2+1)/(N2+1);
Mint = interp1([FF(1:length(FF)-1) 2-fliplr(FF)], [MM(1:length(MM)-1) fliplr(MM)], reducf);

if (plt==1)
  figure(1);
  clf;
  plot(reducf/2*fs,Mint);
  title('magnitude frequency response (band pass filter)');
  xlabel('frequency in Hz');
  hold off;
end;

%% filters computation

h2r  = real(ifft(Mint,N2+1));
  % note: "Mint" is real+even, therefore the "ifft" should be real;
  % but due to numerical roundings, the imaginary part is not
  % absolutely null and must be thrown away
h2r  = [h2r(N2/2+2:N2+1) h2r(1:N2/2+1)];
h2r  = h2r.*win;  %%% in phase filter

Minti=[0 -i*Mint(2:N2/2+1) i*Mint(N2/2+2:N2+1)];
h2i  = real(ifft(Minti,N2+1));
  % note: "Minti" is imaginary+odd, therefore the "ifft" should be real;
  % but due to numerical roundings, the imaginary part is not
  % absolutely null and must be thrown away
h2i  = [h2i(N2/2+2:N2+1) h2i(1:N2/2+1)];
h2i  = h2i.*win;  %%% in quadrature filter


%%
if (plt==1)
  [hh2r, ww2r] = freqz(h2r,1,fs,'whole');
  [hh2i, ww2i] = freqz(h2i,1,fs,'whole');

  figure(2);
  clf;
  grid on;
  hold on;
  plot(ww2r/pi*fs/2,abs(hh2r),'b');
  plot(ww2i/pi*fs/2,abs(hh2i),'r');
  title('magnitude frequency response (b: in phase band pass filter; r: in quadrature band pass filter)');
  hold off;
end;

%% Hilbert filtering itself
dh        = 1;  %%% decimation during the Hilbert filtering (1: nothing; or 2)
signal10r = filter(h2r,1,signal5);
signal10i = filter(h2i,1,signal5);
signal10r = signal10r(1:dh:length(signal10r)); %%% decimation dh
signal10i = signal10i(1:dh:length(signal10i)); %%% decimation dh
signal10  = signal10r + i*signal10i;
lsig10    = length(signal10);
module10  = abs(signal10);

%%% time propagation until the Hilbert filtering
tpgh = 0;
if (lpfd>0)
  tpgh = (N5/2);
end;
tpgh = round(tpgh + ndec*(N2/2));


%%%
if (plt==1)
  figure(3);
  clf;
  plot([0:lsig-1]/fs,signal);
  grid on;
  title('original signal');
  hold off;

  figure(4);
  clf;
  plot([0:lsig-1]/fs,[real(signal10(tpgh+1:lsig));zeros(tpgh,1)],'b');
  hold on;
  plot([0:lsig-1]/fs,[imag(signal10(tpgh+1:lsig));zeros(tpgh,1)],'r');
  grid on;
  title('signal after Hilbert -- real part (b) and imaginary part (r)');
  hold off;
end;


%%% Smoothing

%% Smoothing filter
NL  = round(NLt*fs/(ndec*dh));  %%% filter length (in samples)
if (NL==0)
  NL = 1;
end;
smth     = blackman(NL)/(sum(blackman(NL)));   %%% shape of the smoothing filter

module10 = filter(smth, 1, module10);

%%

moduleAS = zeros(lsig,1);
moduleAS(1:ndec*dh:lsig) = module10;
moduleAS = filter(ones(ndec*dh,1),1,moduleAS);


%%% propagation time

tpgAS = 0;
if (lpfd>0)
  tpgAS = (N5/2);
end;
tpgAS = round(tpgAS + ndec*(N2/2) + ndec*dh*(NL/2));
moduleAS = [moduleAS(tpgAS+1:lsig); ones(tpgAS,1)*moduleAS(lsig)./[1:tpgAS]']; %%% the result is artificially extended


%%% saving files
save('-ascii', prename+extsnrwin+'-AS.txt',  'moduleAS');


%% thresholding
paramAS  = threshparam(3);
switch threshmeth
  case 1
    %%% 1 'empirical' parameter
    segmeAS  = (moduleAS>=paramAS*max(moduleAS));  %%% rough thresholding

  case 2
    %%% 2 'empirical' parameters
    param1 = 0.95;

    vvvAS     = sort(moduleAS);
    vchoiceAS = vvvAS(round(length(vvvAS)*param1));

    segmeAS   = (moduleAS>=paramAS*vchoiceAS);       %%% slightly more subtile

  case 3
    %%% 1 'empirical' parameter

    threshAS   = min([paramAS*mean(moduleAS) max(moduleAS)]);
      %%% the 'min' is used here in order to avoid to get a threshold
      %%% bigger than the max value

    segmeAS    = (moduleAS>=threshAS);  %%% thresholding

  case 4
    %%% Otsu method
    nhist = 50;

    [hhAS, ppAS] = hist(moduleAS, nhist);
    hhAS = hhAS/sum(hhAS);
    muT  = sum([0:nhist-1].*hhAS);
    maxsigm = 0;
    for jj=0:nhist-1-1  %%% second '-1': because 'w1' becomes = 0
      w0   = sum(hhAS(1:jj+1));
      w1   = 1-w0;
      mutt = sum([0:jj].*hhAS(1:jj+1));
      sigm = w0*w1*( (muT-mutt)/w1 - mutt/w0 )^2;
      if (sigm>maxsigm)
        maxsigm  = sigm;
        threshAS = ppAS(jj+1)*paramAS;
      end;
    end;
    segmeAS = (moduleAS>=threshAS);  %%% thresholding
end;

if (plt==1)
  figure(5);
  clf;
  plot([0:lsig-1]/fs,signal/max(abs(signal))*max(moduleAS)*0.9,'g');
  hold on;
  plot([0:lsig-1]/fs,moduleAS);
  plot([0 lsig-1]/fs,[threshAS threshAS]);
  plot([0:lsig-1]/fs,segmeAS*max(moduleAS),'r');
  title('signal; magnitude A; threshold position, segmentation result');
  xlim([0 lsig-1]/fs);
  hold off;
end;



%%%%%% for BE3 -- f0 analysis



signalAS=[signal10(tpgh+1:lsig);zeros(tpgh+1,1)];
for ii=1:length(signalAS)-1
  freqinst(ii) = (angle(signalAS(ii+1))-angle(signalAS(ii)))/2/pi*fs;

  %%% a slight correction for outliers
  if ( (freqinst(ii)<0. || freqinst(ii)>300.) && ii>1)
    freqinst(ii)=freqinst(ii-1);
  end;
end;

%%% some soft smothing
smthf0=blackman(300)/(sum(blackman(300)));   %%% shape of the smoothing filter
freqinst10 = filtfilt(smthf0, 1, freqinst);
freqinst10=freqinst10.*segmeAS';

echtemps=[0:lsig-1]/fs;

if (plt==1)
  figure(6);
  clf;
  plot([0:lsig-1]/fs,freqinst10,'r');
  hold on;
  plot([0:lsig-1]/fs,signal/max(abs(signal))*300*0.9,'b');
  ylim([-300 300]);
  hold off;
end;

