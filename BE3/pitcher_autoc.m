%%% Pitch-tracker based on the autocorrelation function
%%%
%%% usage:
%%%    [tt, voice, f0] = pitcher_autoc(namefile)
%%% where "namefile" is a ".wav" file
%%%
%%% S. Rossignol -- 19/01/2016 -- 2017 -- 2018



%%% you can use this code as a matlab/octave SCRIPT or as matlab/octave FUNCTION
%%% in order to change, comment/uncomment the lines just below

%%% standalone -- octave script
%%% standalone -- octave script
%clear all;
%close all;
%fe=8000;                   %%% sampling rate
%parametres=[0.98784626 0.4562832 0.0001777652 0.01163393 0.04 60 600 -0.7204987 6.8918675 249 0.7672432 0.98784626 0.4562832 0.0001777652 0.01163393];
%pp1=2*(parametres(6)-parametres(8))/fe;if (pp1<0.) pp1=0.; end;
%pp2=2*(parametres(7)+parametres(9))/fe;if (pp2>1.) pp2=1.; end;
%aa=1;
%bb = fir2 (parametres(10), [0 (parametres(6)-parametres(8))/fe*2 (parametres(7)+parametres(9))/fe*2 (parametres(7)+100)/fe*2 1], [0 1 parametres(11) 0 0]);
%namefile='parenthese.wav';
%%% function -- octave function
function [ttf0, voice, f0da] = pitcher_autoc(namefile, parametres, bb, aa)
fe=8000;                   %%% sampling rate



NN=ceil(fe*parametres(5)); %%% frame size
%stepsize=ceil(fe*0.01);    %%% stepsize between two successive frames
stepsize=1;                 %%% like for 'voicingf0AS_BE.m'

[sigs, fes] = audioread(namefile);

%%% the pitch-tracker works for signals sampled at 8000 Hz
if (fes~=fe)
  sigj = resample(sigs, fe, fes);
else
  sigj = sigs;
end;
lel = length(sigj);

%%% the pitch must belong to this interval
f0min=parametres(6);
f0max=parametres(7);

%%% band-pass filtering
sigi=filtfilt(bb,aa,sigj);

mintau=floor(fe/f0max)+1; %%% "tau=0" corresponds to the bin "mintau=1"
maxtau=ceil(fe/f0min)+1;  %%% "tau=0" corresponds to the bin "mintau=1"
if (maxtau>=NN)
  maxtau=NN-1;
end;


%%% 
thresha=parametres(12); %%% this threshold has been optimized on artificial signals
threshb=parametres(13); %%% threshold for voicing decision (not optimized) => distance to the segment
threshc=parametres(14); %%% threshold for voicing decision (not optimized) => energy
threshd=parametres(15); %%% threshold for post-processing 2

kk=0;
for jj=1:stepsize:lel-NN+1
  kk=kk+1;
  signal=sigi(jj:jj+NN-1);
  signal=signal-mean(signal);

  ttf0(kk)=(jj+NN/2)/fe;

  vari3 = xcorr(signal,'biased');
  ene(kk)=max(vari3);

  vari3a=vari3(NN:end);

  %%% short-term processing -- voicing
  voice(kk)=0;
  MM=vari3a(1);
  LL=length(vari3a)-1;

  dist = abs(vari3a(mintau:maxtau) - (-MM/LL*[mintau:maxtau]+MM)' );

  %%% voicing decision made here
  if (min(dist/vari3a(1)) < threshb && ene(kk)>threshc)
    voice(kk)=1;
  end;

  %%% short-term processing -- f0
  vari3a(1:mintau-1)=vari3(mintau);
  vari3a(maxtau+1:end)=vari3(maxtau);
  vari3m = (vari3a - min(vari3a))/(max(vari3a)-min(vari3a));

  [vv,ppv] = max(vari3m); %%% initialisation of f0da

  if (ppv>1)
    f0da(kk)=fe/(ppv-1);
  else
    f0da(kk)=fe/ppv;
  end;
  ii=mintau;
  while ii<=maxtau
    if ( ((vari3m(ii-1)<=vari3m(ii) && vari3m(ii+1)<vari3m(ii)) || ...
          (vari3m(ii-1)<vari3m(ii)  && vari3m(ii+1)<=vari3m(ii))) && ...
           vari3m(ii)>thresha )
      pp=polyfit([ii-1 ii ii+1],[vari3m(ii-1) vari3m(ii) vari3m(ii+1)],2);

      %%% it is slowing down things
      if (pp(1)~=0)
        iim=-pp(2)/2/pp(1);
      else
        iim=ii;
      end;

      f0da(kk) = fe/(iim-1); %%% "ii=1" corresponds to "tau=0"

      ii=maxtau+1;
    else
      ii=ii+1;
    end;
  end;
  f0dares(kk)=f0da(kk); %%% for post-processing purposes
  if (voice(kk)==0)
    f0da(kk)=0.;
  end;
end;

%%% post-processing 1
for jj=3:length(voice)-2
  if ( voice(jj)==0 && voice(jj-1)==1 && voice(jj-2)==1 && voice(jj+1)==1 && voice(jj+2)==1 )
    voice(jj)=1;
    f0da(jj)=f0dares(jj);
  end;
  if ( voice(jj)==1 && voice(jj-1)==0 && voice(jj-2)==0 && voice(jj+1)==0 && voice(jj+2)==0 )
    voice(jj)=0;
    f0da(jj)=0.;
  end;
end;

%%% post-processing 2
for jj=2:length(voice)-1
  if (f0da(jj-1)~=0 && f0da(jj)~=0)
    if ( abs(f0da(jj-1)-f0da(jj+1))/f0da(jj-1)<threshd && ...
         abs(f0da(jj)-f0da(jj+1))/f0da(jj)>threshd     && ...
         abs(f0da(jj)-f0da(jj-1))/f0da(jj)>threshd )
      f0da(jj) = ( f0da(jj+1) +  f0da(jj-1) )/2;
    end;
  end;
end;

