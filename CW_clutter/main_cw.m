clear
close all
clc

clear
clc
close all
%   BL=[0:1:100]
%% Parameters
i=1;

BL=10; %baseline
baseline=floor(BL/0.6)*0.6; %adjust the baseline to the gid

freq=7.36e9; % frequency 
c = 3e8;
t_end=10;
sRate = 200; %sample rate
t=0:1/sRate:t_end; 
lambda=c/freq;
grid_dimensions=[baseline+100 baseline]; %set grid dimension

hTx=1;
blade_length=0.45; %blade length
hTg=1;           %fan height
hRx=hTx;

epsR=15; 
sigma_ground=0.005;


% 
Rx_position=[(grid_dimensions(1)-baseline)/2; 0; hTx];
Tx_position=[Rx_position(1)+baseline; 0; hRx];
[clutterPowerMatrix(:,:),y_gridPoints,x_gridPoints,kk,qq,maxCPM,Doppler_signature_clutter]=Clutter_pisa(hTx,freq,BL,t_end);
save clutterPowerMatrix.mat
save clutterPowerMatrixnorm.mat
% 
% f_l_p=2;
% [z1, p1, k1] = butter(9,f_l_p/(sRate),'low');                            % Hardware LPF of 5 order - filter coefficents
% [SOS1, gain1] = zp2sos(z1, p1, k1);                                        % Convert coefficients to SOS form
% h_LP = dfilt.df2tsos(SOS1, gain1);
% % fvtool(h_LP)
% 
% filtered_DS = filter(h_LP,(Doppler_signature));


%% plot

% Received power
% figure()
% plot(t,(real(Doppler_signature)),'k','LineWidth',2)
% grid on
% xlim([min(t) max(t)])
% xlabel('Time, s', 'FontSize', 12)
% ylabel('Received amplitude, Volt', 'FontSize', 12)
% set(gca,'FontSize',12);
% hold on
% plot(t,filtered_DS,'g','LineWidth',2)
% legend('Received signal', '2Hz LP filtered')

%Spectrogram
% Nhamming=256*4;
% figure()
                        

%% Power Spectral Density
[f, h]=pwelch(real(Doppler_signature_clutter),[],[],[],sRate,'onesided');
% [f_lpf h_lpf]=pwelch(real(filtered_DS),[],[],[],sRate,'onesided');%
figure()
semilogx(h,10*log10(f),'k','LineWidth',2);
xlabel('Frequency, Hz', 'FontSize', 12);
ylabel('PSD, dB/Hz', 'FontSize', 12);
grid on
xlim([min(h) max(h)])
set(gca,'FontSize',12);
% hold on
% semilogx(h_lpf,10*log10(f_lpf),'g','LineWidth',2);
% legend('Received signal', '2Hz LP filtered')
%  
% figure
% hist(sqrt(abs(Doppler_signature.^2)))
% histfit((sqrt(abs(Doppler_signature.^2))),[],'rayleigh')
% 
% 
% figure
% hist(real(Doppler_signature),sqrt(length(Doppler_signature)))
% 
% figure
% hist(imag(Doppler_signature),sqrt(length(Doppler_signature)))
% fft_DS=fftshift(fft(Doppler_signature));
% fft_DS_lpf=fftshift(fft(filtered_DS));
% fr=-sRate/2:sRate/length(Doppler_signature):sRate/2-sRate/length(Doppler_signature);

% figure
% plot(fr,10*log10(real(fft_DS)),'k','LineWidth',2)
% grid on
% hold on
% plot(fr,10*log10(real(fft_DS_lpf)),'g','LineWidth',2)
% xlabel('Frequency [Hz]')
% ylabel('Spectrum [dB]')

