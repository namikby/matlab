clear 
close all
clc

%%%%% generating %%%%%% - Binary Random Generator

SGBT=randi([0,1],1,10);
fb=50000; % frequency of source in bps or Hz
Tb=1/fb; % period of one bit
t=0:Tb/20:(10*Tb)-Tb/20; 

figure(1)
stem(SGBT)
title('Izlaz slucajnog generatora binarnog toka')

%%%%% CODING %%%%% - LINE CODER

NRZ=nrz_coder(SGBT);
MAN=man_coder(SGBT); % Manchester coder makes double samples because of two states

% upsampling to 1 MHz
NRZ=up_sample(NRZ,20);
MAN=up_sample(MAN,10); % double the samples, half the upsampling

figure(2)
subplot(2,1,1)
plot(t,MAN)
title('Manchester')
subplot(2,1,2)
plot(t,NRZ)
title('NRZ')

%%%%% filter %%%%% - CHANNEL
A=1.95;
B=59;
delta=0.03;
l=2;
fc=1e6;
fs=2e6;
f=0:100:fs;

%filter design
gama=(1-(1i*delta/2))*(A*sqrt(f/fc)+1i*B*(f/fc));
H=exp(-gama*l);

% convolution with filter
chMAN=conv(fir2(10,f/fs,abs(H)),MAN);
chNRZ=conv(fir2(10,f/fs,abs(H)),NRZ);

lng=length(chMAN);

MAN=5*chMAN(5:lng-6); % amplifier 
NRZ=5*chNRZ(5:lng-6);

figure(3)
subplot(2,1,1)
plot(t,MAN)
title('Manchester on channel')
subplot(2,1,2)
plot(t,NRZ)
title('NRZ on channel')

%%%%% AWGN %%%%%
MAN=awgn(MAN,1); % SNR set to 1
NRZ=awgn(NRZ,1);

figure(4)
subplot(2,1,1)
plot(t,MAN)
title('Manchester on channel with AWGN')
subplot(2,1,2)
plot(t,NRZ)
title('NRZ on channel with AWGN')

%%%%% optimal filter %%%%%
x=[5*ones(10,1) -5*ones(10,1)];
H_MAN=conj(x(end:-1:1));
MAN_filt=filter(H_MAN,1,MAN);

MAN_filt=[MAN_filt(15:end) zeros(1,14)];
MAN=MAN_filt/100;

x=5*ones(20,1);
H_NRZ=conj(x(end:-1:1));
NRZ_filt=filter(H_NRZ,1,NRZ);

NRZ_filt=[NRZ_filt(10:end) zeros(1,9)];     % shifting left for moment of decision
NRZ=NRZ_filt/100;                           % scaling the output from (-500,500)

figure(5)
subplot(2,1,1)
plot(t,MAN)
title('Manchester after optimal filter')
subplot(2,1,2)
plot(t,NRZ)
title('NRZ after optimal filter')

%%%%% DETECTOR %%%%%%
MAN=detect_MAN(MAN);
NRZ=detect_NRZ(NRZ);

figure(6)
subplot(2,1,1)
stem(MAN)
title('Manchester after detection')
subplot(2,1,2)
stem(NRZ)
title('NRZ after detection')

%%%%% BER %%%%%

[NoNRZ,EbNRZ]=biterr(NRZ,SGBT);
[NoMAN,EbMAN]=biterr(MAN,SGBT);


%%%%%%%%%% functions %%%%%%%%%%%

function [output] = up_sample(in, n) % upsampling to get pseudoanalog signal
  output = upsample(in,n);
  for i = 2:length(output)
    if output(i) == 0
      output(i) = output(i-1);
    end
  end
end

function [out] = nrz_coder(in) % NRZ coder
  n = length(in);
  out = zeros(1,n); 
  for i=1:n
    if in(i) == 1
      out(i) = 5;
    else
      out(i) = -5;
    end
  end
end

function [out] = man_coder(in) % Manchester coder
  n = length(in);
  out = zeros(1,2*n);
  j=1;
  for i=1:n
      if in(i) == 1
          out(j) = 5;
          out(j+1) = -5;
      else
          out(j) = -5;
          out(j+1) = 5;
      end
      j=j+2;
  end
end

function [out] = detect_NRZ(in)  % detector for NRZ code
  n=length(in);
  out=zeros(1,n/20);
  j=1;
  for i=10:20:n
      if in(i)>0
          out(j)=1;
      else
          out(j)=0;
      end
      j=j+1;
  end
end

function [out] = detect_MAN(in)  % detector for Manchester code
  n=length(in);
  out=zeros(1,n/20);
  j=1;
  for i=5:20:n
      if in(i)>0
          out(j)=1;
      else
          out(j)=0;
      end
      j=j+1;
  end
end