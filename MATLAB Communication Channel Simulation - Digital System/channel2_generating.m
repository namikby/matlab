clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% GENERISANJE 1e5 BITS AND BER %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TX=randi([0,1],1,1e5);
Rb=1e6; 
Tb=1/Rb; 
fb=1/Tb; 
f=2*fb; 
t=0:Tb/40:(1e5*Tb)-Tb/40; 
car=sqrt(2/fb)*(cos(2*pi*f*t)); 


OUTPUT=nrz_coder(TX);
OUTPUT=up_sample(OUTPUT,40);
OUTPUT=OUTPUT.*car;

h=[-0.015 0.058 -0.350 1.000 -0.350 0.058 -0.005];

OUTPUT=conv(OUTPUT,h); 
OUTPUT1=OUTPUT(4:end-3); 

No_svi_ZF=[];
Eb_svi_ZF=[];
No_svi_MMSE=[];
Eb_svi_MMSE=[];

for i=1:6
    clear OUTPUT RX ZF_OUTPUT MMSE_OUTPUT ZF_RX MMSE_RX
    
    %%%%%%%%%%%%%%%
    %%%%% SNR %%%%%
    %%%%%%%%%%%%%%%
    
    SNR=[1,5,10,15,20,25];
    OUTPUT=awgn(OUTPUT1,SNR(i),'measured');
    
    %%%%%%%%%%%%%%%%%%
    %%%%% ZF EQ %%%%%%
    %%%%%%%%%%%%%%%%%%
    
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
    
    %%%%%%%%%%%%%%%%%%%
    %%%%% MMSE EQ %%%%%
    %%%%%%%%%%%%%%%%%%%
    
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
    
    ZF_OUTPUT=car.*ZF_OUTPUT;
    ZF_OUTPUT=1000*ZF_OUTPUT;
    MMSE_OUTPUT=car.*MMSE_OUTPUT;
    MMSE_OUTPUT=1000*MMSE_OUTPUT;
    
    ZF_RX=reshape(ZF_OUTPUT,40,[]);
    ZF_RX=sum(ZF_RX);
    ZF_RX=detector(ZF_RX);
    
    MMSE_RX=reshape(MMSE_OUTPUT,40,[]);
    MMSE_RX=sum(MMSE_RX);
    MMSE_RX=detector(MMSE_RX);
    
    [No_ZF,Eb_ZF]=biterr(TX,ZF_RX);
    [No_MMSE,Eb_MMSE]=biterr(TX,MMSE_RX);
    
    No_svi_ZF(i)=No_ZF;
    Eb_svi_ZF(i)=Eb_ZF;
    No_svi_MMSE(i)=No_MMSE;
    Eb_svi_MMSE(i)=Eb_MMSE;
end

figure(1)
subplot(2,1,1)
scatter(Eb_svi_ZF,No_svi_ZF)
title('ZF ekvalizator')
xlabel('Eb');
ylabel('N0');
subplot(2,1,2)
scatter(Eb_svi_MMSE,No_svi_MMSE)
title('MMSE ekvalizator')
xlabel('Eb');
ylabel('N0');

function [out] = nrz_coder(in) 
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

function [output] = up_sample(in, n) 
  output = upsample(in,n);
  for i = 2:length(output)
    if output(i) == 0
      output(i) = output(i-1);
    end
  end
end

function [out] = detector(in)  
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