close all;
clear all;
clc;

c=2.9979*10^8;  % Propagation velocity
fo=4.5e9;       % Carrier frequency (4.5GHz)
landa=c/fo;     % Wavelength (60cm for fo = 4.5e9)
theta_inc=pi/4; % Loook angle
Vplat=200;         % Velocity of platform
Rcenter=2000;  % Range distance to center of target area
l=0.5;           % Antenna Range length actual
%L=1.5;             %Antenna azimuth length
L=0.5;


Swath=(landa*Rcenter)/(l*cos(theta_inc));


Rnear=Rcenter-Swath/2;     % Range Near
Rfar=Rcenter+Swath/2;   % Rang Far
Tr=(2*Swath)/(c*9);     %Pulse duration=>window length is 10 Tr
Ts=(2*Rnear)/c-Tr/2;             % Start time of sampling
Tf=(2*Rfar)/c+Tr/2;           % End time of sampling

Bw=100e6;%therefore resolution is 1.5m
Kr=Bw/Tr;%Modulation Rate
Fs=2*Bw;
dt=1/Fs;
t=Ts:dt:Tf;
%%%%%target location%%%%%%%%%%%%%%%%
%R0=(Rcenter+Rnear)/2;
R0=Rcenter;
%R0=(Rcenter+Rfar)/2;
Timg=(landa*R0)/(L*Vplat);
eta0=0;
%eta0=-Timg/2;
%eta0=Timg/2;
%%%%%%%%%%%%%%%%%%

%max_PRF=(l*c)/(2*Rcenter*landa*tan(theta_inc));
max_PRF=(l*c)/(2*Rfar*landa*tan(theta_inc));
min_PRF=2*Vplat/L;
%PRF=max_PRF+min_PRF/2;
%PRF=900;
%PRF=600;
PRF=2*min_PRF;
Dta=1/PRF;
eta=-Timg/2:Dta:Timg/2;



R=R0+( ( (Vplat^2)*((eta-eta0).^2) )/(2*R0) );
Tdelay=2*R/c;
figure;
plot(eta,Tdelay);
hold on;

title('delay of received echo versus slow time');

counter=1;
Sig1=zeros(length(eta),length(t));
for Td=Tdelay
    Sig1(counter,:)=exp(-1j*2*pi*fo*Td)*exp( 1j*Kr*pi*((t-Td).^2) ) .* (-Tr/2+Td<=t & t<=Tr/2+Td);
    counter=counter+1;
end
figure;
imagesc(abs(Sig1));
 figure;
 mesh(abs(Sig1));

T=-Tr/2:dt:Tr/2;
St=exp(1j*Kr*pi*T.^2).*(-Tr/2<=T & T<=Tr/2);
len=length(t)+length(T)-1;
%Nfft=2^nextpow2(len);
Nfft=len;
St_f=fft(St,Nfft);
Ht=exp(-1j*Kr*pi*T.^2).*(-Tr/2<=T & T<=Tr/2);
H_f=fft(Ht,Nfft);

Sig2=zeros(length(eta),Nfft);

for counter=1:length(eta)
    Sig2(counter,:)=ifft( fft(Sig1(counter,:),Nfft) .* H_f , Nfft );
end

figure;
imagesc(abs(Sig2));
 figure;
 mesh(abs(Sig2));

Ka=(2*Vplat^2)/(landa*R0);
Sa=exp(-1j*Ka*pi*eta.^2);
lena=2*length(eta)-1;
%Nffta=2^nextpow2(lena);
Nffta=lena;

f=linspace(-PRF/2,PRF/2,Nffta);
RCM_range_dopler=((R0*landa^2)/(8*Vplat^2))*f.^2;
RCM_range_dopler=2*RCM_range_dopler/c;
RCM_range_dopler=RCM_range_dopler/dt;
RCM_range_dopler=round(RCM_range_dopler);

Sa_f=fft(Sa,Nffta);
Ha=exp(1j*Ka*pi*eta.^2);
Ha_f=fft(Ha,Nffta);
Ha_f=conj(Ha_f');
%Ha_f=(Sa_f)';


Sig3=zeros(Nffta,Nfft);

for counter=1:Nfft
    Sig3(:,counter)=fftshift(fft(Sig2(:,counter),Nffta));
end


for counter=1:Nffta
    Sig3(counter,:)=[ Sig3(counter,RCM_range_dopler(counter)+1:end-RCM_range_dopler(counter)),Sig3(counter,end-RCM_range_dopler(counter)+1:end),Sig3(counter,1:RCM_range_dopler(counter)) ];
end

figure;
imagesc(abs(Sig3));
figure;
mesh(abs(Sig3));

Sig4=zeros(Nffta,Nfft);

for counter=1:Nfft
    Sig4(:,counter)=ifft( Sig3(:,counter) .*fftshift( Ha_f ),Nffta);
end


 figure;
 imagesc(abs(Sig4));
 figure;
  mesh(abs(Sig4));


RCM=(0.886^2)*(R0*landa^2)/(8*L^2)



