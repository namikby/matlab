clear all
close all
clc

%%%%%%%%%%%%%%%
%%%%% TX %%%%%%
%%%%%%%%%%%%%%%

TX=randi([0,1],1,10);
Rb=1e6; % source 1Mbps
Tb=1/Rb; % period of one bit
fb=1/Tb; % source frequency
f=2*fb; % carrier frequency
t=0:Tb/40:(10*Tb)-Tb/40; % 40 samples per period
car=sqrt(2/fb)*(cos(2*pi*f*t)); % carrier signal

figure(1)
stem(TX)
title('Output of binary random source')

%%%%%%%%%%%%%%%%%%%%
%%%%% BPSK MOD %%%%%
%%%%%%%%%%%%%%%%%%%%

OUTPUT=nrz_coder(TX); % coding nrz signal with amplitude (1,-1)
OUTPUT=up_sample(OUTPUT,40); % repeating sample 40 times

figure(2)
plot(t,OUTPUT)
title('NRZ signal')
xlabel('Time(s)');
ylabel('Amplitude(V)');

OUTPUT=OUTPUT.*car; % modulation

figure(3)
plot(t,OUTPUT)
title('BPSK modulated signal')
xlabel('Time(s)');
ylabel('Amplitude(V)');

%%%%%%%%%%%%%%%%%
%%%%% CHANNEL %%%%%
%%%%%%%%%%%%%%%%%

h=[-0.015 0.058 -0.350 1.000 -0.350 0.058 -0.005]; % impulse responce of channel

OUTPUT=conv(OUTPUT,h); % convolution of signala with impulse response of channel
OUTPUT=OUTPUT(4:end-3); % irrelevant samples delete

figure(4)
plot(t,OUTPUT)
title('BPSK modulated signal in channel')
xlabel('Time(s)');
ylabel('Amplitude(V)');

%%%%%%%%%%%%%%%%
%%%%% AWGN %%%%%
%%%%%%%%%%%%%%%%

OUTPUT=awgn(OUTPUT,1,'measured'); % adding AWGN

figure(5)
plot(t,OUTPUT)
title('BPSK modulated signal in channel with AWGN')
xlabel('Time(s)');
ylabel('Amplitude(V)');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ZF EQUALIZATOR %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

X_ZF=[1.00 -0.350 0.058;
     -0.350 1.00 -0.350;
      0.058 -0.350 1.00];
Z_ZF=[0;1;0];
c_ZF=inv(X_ZF)*Z_ZF;
c_ZF=c_ZF';

ZF_eq=zeros(1,length(h));

for i=2:(length(h)-1)
    ZF_eq(i)=h(i+1)*c_ZF(1)+h(i-1)*c_ZF(3);
end

ZF_eq(1)=h(2)*c_ZF(1)+h(1)*c_ZF(2);
ZF_eq(7)=h(7)*c_ZF(2)+h(6)*c_ZF(3);

ZF_OUTPUT=conv(OUTPUT,ZF_eq);
ZF_OUTPUT=ZF_OUTPUT(4:end-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MMSE EQUALIZATOR %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_MMSE=[[-0.35 0 0];
       [1 -0.35 0];
       [-0.35 1 -0.35];
       [0 -0.35 1];
       [0 0 -0.35]];

Z_MMSE=[0;0;1;0;0];
   
Rxx=X_MMSE'*X_MMSE;
Rxz=X_MMSE'*Z_MMSE;

c_MMSE=inv(Rxx)*Rxz;
c_MMSE=c_MMSE';

MMSE_OUTPUT=conv(OUTPUT,c_MMSE);
MMSE_OUTPUT=MMSE_OUTPUT(1:end-2);

figure(6)
subplot(2,1,1)
plot(t,ZF_OUTPUT)
title('BPSK modulated signal in channel with AWGN and ZF eq')
xlabel('Vrijeme(s)');
ylabel('Amplituda(V)');
subplot(2,1,2)
plot(t,MMSE_OUTPUT)
title('BPSK modulated signal in channel with AWGN and MMSE eq')
xlabel('Time(s)');
ylabel('Amplitude(V)');

%%%%%%%%%%%%%%%%%%%%%%
%%%%% BPSK DEMOD %%%%%
%%%%%%%%%%%%%%%%%%%%%%

ZF_OUTPUT=car.*ZF_OUTPUT; % demodulation by multiplication with carrier
ZF_OUTPUT=1000*ZF_OUTPUT;
MMSE_OUTPUT=car.*MMSE_OUTPUT;
MMSE_OUTPUT=1000*MMSE_OUTPUT;

figure(7)
subplot(2,1,1)
plot(t,ZF_OUTPUT)
title('BPSK demodulated signal after ZF eq')
xlabel('Time(s)');
ylabel('Amplitude(V)');
subplot(2,1,2)
plot(t,MMSE_OUTPUT)
title('BPSK demodulated signal after MMSE eq')
xlabel('Time(s)');
ylabel('Amplitude(V)');

%%%%%%%%%%%%%%
%%%%% RX %%%%%
%%%%%%%%%%%%%%

ZF_RX=reshape(ZF_OUTPUT,40,[]); % one period in one column
ZF_RX=sum(ZF_RX); % summing each column to one period
ZF_RX=detector(ZF_RX); % is sume greater than 1

MMSE_RX=reshape(MMSE_OUTPUT,40,[]);
MMSE_RX=sum(MMSE_RX);
MMSE_RX=detector(MMSE_RX);

figure(8)
subplot(2,1,1)
stem(ZF_RX)
title('Recieved bits ZF')
subplot(2,1,2)
stem(MMSE_RX)
title('Recieved bits MMSE')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BER %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[No_ZF,Eb_ZF]=biterr(ZF_RX,TX);
[No_MMSE,Eb_MMSE]=biterr(MMSE_RX,TX);

function [out] = nrz_coder(in) % NRZ coder
  n = length(in);
  out = zeros(1,n); 
  for i=1:n
    if in(i) == 1
      out(i) = 1;
    else
      out(i) = -1;
    end
  end
end

function [output] = up_sample(in, n) % upsampling to get pseudoanalog signal
  output = upsample(in,n);
  for i = 2:length(output)
    if output(i) == 0
      output(i) = output(i-1);
    end
  end
end

function [out] = detector(in)  % detector
  n=length(in);
  out=zeros(1,n);
  for i=1:n
      if in(i)>=0
          out(i)=1;
      else
          out(i)=0;
      end   
  end
end