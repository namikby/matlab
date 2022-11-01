clear 
close all
clc

%%%%% generating 1e6 bits %%%%%
SGBT=[];
Nt=1e6;
Ns=1e4;
k_max=Nt/Ns;

for k=1:k_max
    SGBT=[SGBT randi([0,1],1,Ns)];
end

EbMAN1=zeros(1,5);
EbNRZ1=zeros(1,5);
NoNRZ1=zeros(1,5);
NoMAN1=zeros(1,5);

for i=1:5
clear NRZ* MAN* chMAN chNRZ x
SNR=[1,5,10,15,20];

NRZ=nrz_coder(SGBT);
MAN=man_coder(SGBT);

NRZ=up_sample(NRZ,20);
MAN=up_sample(MAN,10);

A=1.95;
B=59;
delta=0.03;
l=2;
fc=1e6;
fs=2e6;
f=0:100:fs;

gama=(1-(1i*delta/2))*(A*sqrt(f/fc)+1i*B*(f/fc));
H=exp(-gama*l);

chMAN=conv(fir2(10,f/fs,abs(H)),MAN);
chNRZ=conv(fir2(10,f/fs,abs(H)),NRZ);

lng=length(chMAN);

MAN=5*chMAN(5:lng-6);
NRZ=5*chNRZ(5:lng-6);

MAN=awgn(MAN,SNR(i));
NRZ=awgn(NRZ,SNR(i));

x=[5*ones(10,1) -5*ones(10,1)];
H_MAN=conj(x(end:-1:1));
MAN_filt=filter(H_MAN,1,MAN);
MAN_filt=[MAN_filt(15:end) zeros(1,14)];
MAN=MAN_filt/100;

x=5*ones(20,1);
H_NRZ=conj(x(end:-1:1));
NRZ_filt=filter(H_NRZ,1,NRZ);
NRZ_filt=[NRZ_filt(10:end) zeros(1,9)];
NRZ=NRZ_filt/100;

MAN=detect_MAN(MAN);
NRZ=detect_NRZ(NRZ);
[NoNRZ,EbNRZ]=biterr(NRZ,SGBT);
[NoMAN,EbMAN]=biterr(MAN,SGBT);

NoMAN1(i)=NoMAN;
EbMAN1(i)=EbMAN;
NoNRZ1(i)=NoNRZ;
EbNRZ1(i)=EbNRZ;

end

figure(1)
plot(EbMAN1,NoMAN1)
title('Manchester Pc=f(Eb/No)')

figure(2)
plot(EbNRZ1,NoNRZ1)
title('NRZ Pc=f(Eb/No)')

function [output] = up_sample(in, n)
  output = upsample(in,n);
  for i = 2:length(output)
    if output(i) == 0
      output(i) = output(i-1);
    end
  end
end

function [out] = nrz_coder(in)
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

function [out] = man_coder(in)
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

function [out] = detect_NRZ(in)
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

function [out] = detect_MAN(in)
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