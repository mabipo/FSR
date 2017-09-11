%Dopper target signature  Sr=-A(t) sin()

close all
clear all
clc
% f_c=150e6;
f_c=7.36e9;
c=3e8;
sRate=200;
hf=1;
hRx=1;
hTx=1;
BL=50;
grid_dimensions=[BL+40 BL/2];
tx=[(grid_dimensions(1)-BL)/2; 0; hTx];
rx=[tx(1)+BL; 0; hRx];
v=[0 2 0]; %[m/s]
x_tg=grid_dimensions(1)/2;
z_tg=hf;

t_end=grid_dimensions(2)/v(2);
tt=0:1/sRate:t_end;
y_tg=-grid_dimensions(2)/2;

tg0=[x_tg,y_tg,z_tg];
tg_inc=(v'*tt)';                                                          %receiver position increment
pos_tg=[tg0(1)+tg_inc(:,1) tg0(2)+tg_inc(:,2) tg0(3)+tg_inc(:,3)]; 

tx_tg=sqrt((tx(1)-pos_tg(:,1)).^2+(tx(2)-pos_tg(:,2)).^2+(tx(3)-pos_tg(:,3)).^2);
    rx_tg=sqrt((rx(1)-pos_tg(:,1)).^2+(rx(2)-pos_tg(:,2)).^2+(rx(3)-pos_tg(:,3)).^2);
    Doppler_signature_tg=-sin(2*pi*f_c/c*(tx_tg+rx_tg-BLnews ));   
Power_tg=10*log10(sum(abs(Doppler_signature_tg).^2));
% Power_tg_tot=10*log10(sum(10.^(nonzeros(Power_tg)./10)));

wind=window(@hamming,length(tt));
Doppler_signature=wind'.*Doppler_signature_tg;
 

figure()
plot(tt,Doppler_signature);
xlabel('Time, s')
ylabel('Amplitude, V')
xlim([min(tt) max(tt)])

%% Power Spectral Density
[f, h]=pwelch(real(Doppler_signature_tg),[],[],[],2000,'onesided');
figure()
semilogx(h,10*log10(f),'k','LineWidth',2);
xlabel('Frequency, Hz', 'FontSize', 12);
ylabel('PSD, dB/Hz', 'FontSize', 12);
grid on

xlim([min(h) max(h)])
set(gca,'FontSize',12);
