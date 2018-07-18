%% Adaptive HW3
%
% Srilakshmi Alla
%
% 819663423
%%
clc;
clear all;
close all;
%% Question
% A communication channel is formed by a cascade of four processes as shown
% below. The matched filter operates with the shaping filter to simultaneously
% maximize signal to noise ratio (SNR) and minimize inter symbol interference (ISI).
% The channel induced distortion introduces ISI and the equalizer serves to minimize its effect.
%
% The input is 1000 samples of a complex QPSK random data sequence, the 
% shaping filter is a square root Nyquist filter with excess bandwidth $$ \alpha $$ = 0.5 
% operating at 4-samples per symbol, the channel impulse response is 
% {1  0  0  0  0  0.2  0  0  j*0.1}, and the noise is complex AWGN with $$ \alpha^2 $$=0.01. 

h=rcosine(1,4,'sqrt',0.5,6);        % Shaping Filter
h=h/max(h);                         % Normalize Shaping Filter
h2=reshape([h 0 0 0],4,13);         % Reshape Shaping Filter
 
N_dat=1000;
x0=(floor(2*rand(1,N_dat))-0.5)/0.5+j*(floor(2*rand(1,N_dat))-0.5)/0.5;
x1=zeros(1,4*N_dat);
reg=zeros(1,13);
 
% form Modulator Output
m=0;
for n=1:N_dat
    reg=[x0(n) reg(1:12)]
    for k=1:4
        x1(m+k)=reg*h2(k,:)'
    end
    m=m+4;
end
 
% form Demodulator Output
x4=filter(h,1,x1)/(h*h');   % no noise, No Channel
%% a
% On a single figure show the real part of 100 symbols of the time series
% formed at the output of the shaping filter and then at the output of the 
% matched filter without the channel and additive noise. Also show the
% constellation of the signal formed by the shaping filter and the matched filter.
figure;
subplot(3,2,[1 2]);
plot(0:1/4:100,real(x1(1:401)),'linewidth',2);
hold on;
stem(6:106,real(x0(1:101)),'linewidth',2);
hold off;
axis([0 100 -1.5 1.5]);
grid on;
title('Real part of timeseries at output of Shaping filter','fontsize',10);


subplot(3,2,[3 4]);
plot(0:1/4:100,real(x4(1:401)),'linewidth',2);
hold on;
stem(12:112,real(x0(1:101)),'linewidth',2);
hold off;
axis([0 100 -1.5 1.5]);
grid on;
title('Real part of timeseries at output of Matched filter without channel  and additive noise','fontsize',10);

% subplot (3,2,5)
% plot(x1(1:4:4*N_dat),'ro');
% grid on;
% axis([-1.5 1.5 -1.5 1.5]);
% title('Shaping filter constellation','fontsize',10)

subplot (3,2,5)
plot(x0(1:4:end),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Shaping filter constellation','fontsize',10)


subplot (3,2,6)
plot(x4(1:4:4*N_dat),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Matched filter constellation','fontsize',10)


%% b
% On a single figure show the real part of 100 symbols of the time series 
% formed at the output of the channel and then at the output of the matched
% filter without the additive noise. Also show the constellation of the signal 
% formed by the shaping filter, the channel output, and the matched filter.

x2=filter([1 0 0 0 0.2 0 0 j*0.1],1,x1);
x4=filter(h,1,x2)/(h*h'); % without the additive noise

figure;
subplot(3,3,[1 2 3]);
plot(0:1/4:100,real(x2(1:401)),'linewidth',2);
hold on;
stem(6:106,real(x0(1:101)),'linewidth',2);
hold off;
axis([0 100 -1.5 1.5]);
grid on;
title('Real part of timeseries at output of Channel filter','fontsize',10);


subplot(3,3,[4 5 6]);
plot(0:1/4:100,real(x4(1:401)),'linewidth',2);
hold on;
stem(12:112,real(x0(1:101)),'linewidth',2);
hold off;
axis([0 100 -1.5 1.5]);
grid on;
title('Real part of timeseries at output of Matched filter without additive noise','fontsize',10);

subplot (3,3,7)
plot(x1(1:4:4*N_dat),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Shaping filter constellation','fontsize',10)

subplot (3,3,8)
plot(x2(1:4:4*N_dat),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Channel constellation','fontsize',10)

subplot (3,3,9)
plot(x4(1:4:4*N_dat),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Matched filter constellation','fontsize',10)


x3=x2+0.00*(randn(1,4*N_dat)+j*randn(1,4*N_dat))/sqrt(2);
x4=filter(h,1,x3)/(h*h');   % with noise and Channel



figure;
subplot(3,3,[1 2 3]);
plot(0:1/4:100,real(x2(1:401)),'linewidth',2);
hold on;
stem(6:106,real(x0(1:101)),'linewidth',2);
hold off;
axis([0 100 -1.5 1.5]);
grid on;
title('Real part of timeseries at output of Channel filter','fontsize',10);


subplot(3,3,[4 5 6]);
plot(0:1/4:100,real(x4(1:401)),'linewidth',2);
hold on;
stem(12:112,real(x0(1:101)),'linewidth',2);
hold off;
axis([0 100 -1.5 1.5]);
grid on;
title('Real part of timeseries at output of Matched filter with noise and channel','fontsize',10);

subplot (3,3,7)
plot(x1(1:4:4*N_dat),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Shaping constellation','fontsize',10)

subplot (3,3,8)
plot(x2(1:4:4*N_dat),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Channel constellation','fontsize',10)

subplot (3,3,9)
plot(x4(1:4:4*N_dat),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Matched constellation','fontsize',10)

%% c
% Use a 40 tap decision directed LMS algorithm to equalize the channel. 
% The equalizer operates at 4-samples per symbol and is updated at symbol
% rate. On a single page, show the learning curve of the adaptation and show 
% 100 symbols of the time series and the constellation formed at the equalizer
% output. Use mu=0.002.
reg=zeros(1,40);
wts=zeros(1,40);
wts(4+0)=1;
mu=0.002;
 
m=1;
err_sv=zeros(1,N_dat);
for n=1:4*N_dat
    x5(n)=reg*wts';
    if n>40 && rem(n,4)==1 
        xd=sign(real(x5(n)))+j*sign(imag(x5(n)));
        xe=xd-x5(n);
        err_sv(m)=xe;
        m=m+1;
        wts=wts+mu*reg*conj(xe);
    end
    reg=[x4(n) reg(1:39)];
end


figure;
subplot(2,1,1);
plot(abs(err_sv),'linewidth',2);
grid on;
axis([0 N_dat 0 1.5])
title('Equalizer Learning curve,Linear Magnitude')

subplot(2,1,2);
plot(20*log10(abs(err_sv)),'linewidth',2);
grid on;
axis([0 N_dat -70 10])
title('Equalizer Learning curve,Log Magnitude')

figure;

subplot(2,3,[1 2 3]);
plot(0:1/4:100,real(x5(1:401)),'linewidth',2);
hold on;
stem(13:113,real(x0(1:101)),'linewidth',2);
hold off;
axis([0 100 -1.5 1.5]);
grid on;
title('Real part of timeseries at output of Equalizer filter','fontsize',10);


subplot(2,3,4);
plot(x2(1:4:4*N_dat),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Channel output constellation','fontsize',8)

subplot (2,3,5)
plot(x4(1:4:4*N_dat),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Matched filter output constellation','fontsize',8)

subplot (2,3,6)
plot(x5(1:4:4*N_dat),'ro');
grid on;
axis([-1.5 1.5 -1.5 1.5]);
title('Equalizer Filter output constellation','fontsize',8)


%% d
% On a single figure show the real part of the 1000 symbols for the channel 
% output, the matched filter output, and the Equalizer output
figure;
subplot(3,1,1);
plot(real(x2(1:4:4000)),'r.');
grid on;
title('Channel Output','fontsize',10);

subplot(3,1,2);
plot(real(x4(1:4:4000)),'r.');
grid on;
title('Matched filter Output','fontsize',10);

subplot(3,1,3);
plot(real(x5(1:4:4000)),'r.');
grid on;
title('Equalizer Output','fontsize',10);

%% Eye Diagrams
% Not part of homework.Just wanted to look at Eye Diagrams

figure;
subplot(3,1,1)
plot(0,0)
for nn=1:8:4000-8
    
    hold on;
    plot(-1:1/4:1,real(x2(nn:nn+8)))
    hold off;
    grid on;
end
title('Channel output Eye Diagram');

subplot(3,1,2)
plot(0,0)
for nn=1:8:4000-8
    
    hold on;
    plot(-1:1/4:1,real(x4(nn:nn+8)))
    hold off;
    grid on;
end
title('Matched Filter output Eye Diagram');

subplot(3,1,3)
plot(0,0)
for nn=1:8:4000-8
    
    hold on;
    plot(-1:1/4:1,real(x5(nn:nn+8)))
    hold off;
    grid on;
end
title('Equalizer Filter output Eye Diagram');
