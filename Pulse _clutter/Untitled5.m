clear
clc
close all

f_c=7.36e9;
c=3e8;
sRate=200;
hRx=1;   %transmitter and receiver height
hTx=1;
BL=100;
lambda=c/f_c;

step=0.6;
rest=floor(BL/step);
BLnews=BL;
grid_dimensions=[2*BLnews BLnews]; %[horizontal vertical] grid dimentions to define the area to evaluate the clutter returns

tx=[(grid_dimensions(1)-BLnews)/2; 0; hTx];
rx=[tx(1)+BL; 0; hRx];
v=[0 2 0]; %[m/s]
x_tg=grid_dimensions(1)/2;
z_tg=1;
taop=[0.5e-9]%1e-9 2e-9 3e-9 4e-9 5e-9 6e-9 7e-9 8e-9 9e-9 10e-9];
t_oss= 5;
Tr=[0.5e-4]

%% tg
RCS_tg=0.00001

tt=0:1/sRate:t_oss;
y_tg=-b;
tg0=[x_tg,y_tg,z_tg];
tg_inc=(v'*tt)';                                                          %receiver position increment
pos_tg=[tg0(1)+tg_inc(:,1) tg0(2)+tg_inc(:,2) tg0(3)+tg_inc(:,3)];


tx_tg=sqrt((tx(1)-pos_tg(:,1)).^2+(tx(2)-pos_tg(:,2)).^2+(tx(3)-pos_tg(:,3)).^2);
rx_tg=sqrt((rx(1)-pos_tg(:,1)).^2+(rx(2)-pos_tg(:,2)).^2+(rx(3)-pos_tg(:,3)).^2);
Doppler_signature_tg=-sin(2*pi*f_c/c*(tx_tg+rx_tg-BL));

for ti=1:1:length(tt)
Power_tg_ist(ti)=(4*pi*RCS_tg*rx(3)^2*tx(3)^2*pos_tg(3)^4)/(lambda^2*tx_tg(ti)^4*rx_tg(ti)^4);

end 
i=0;
% while 
%     
% Power_tg= 10*log10(sum(Power_tg_ist((i*Tr):i*Tr+taop)));
% end

%     Power_tg(j)=10*log10(max(abs(Doppler_signature_tg).^2));
%     Power_tg_tot(j)=10*log10(sum(10.^(nonzeros(Power_tg)./10)));


wind=window(@hamming,length(tt));
Doppler_signature_tg_w=wind'.*Doppler_signature_tg;

%     figure()
%     plot(tt,Doppler_signature_tg_w);
%     xlabel('Time, s')
%     ylabel('Amplitude, V')
%     xlim([min(tt) max(tt)])

% Power Spectral Density
%     [f, h]=pwelch(real(Doppler_signature_tg),[],[],[],sRate,'onesided');
%     figure()
%     semilogx(h,10*log10(f),'k','LineWidth',2);
%     xlabel('Frequency, Hz', 'FontSize', 12);
%     ylabel('PSD, dB/Hz', 'FontSize', 12);
%     grid on
%     xlim([min(h) max(h)])
%     set(gca,'FontSize',12);

clear Doppler_signature Doppler_signature_tg


%% clutter

%Blade
blade_length=0.5;
% phase0=[0 2*pi/3 4*pi/3];                                                 %initial phase of the three blades
Tmean=1;                                                                   % blade rotation period mean
Tstd=0.25;                                                                   % blade rotation period std
Amean=blade_length;%/2
Astd=0.2;
hf_mean=1;
hf_std=0.5;                                                                       %fan origin height

%grid
% pointStep=[1.05 1.05];            %grid horizontal and vertical step
pointStep=[1.0 1.0]*0.6;            %grid horizontal and vertical step
x_gridPoints=0:pointStep(1):grid_dimensions(1);
if pointStep(2)<2*blade_length
    y_gridPos=blade_length+0.05:pointStep(2):grid_dimensions(2)/2+0.05;
    % y_gridPos=blade_length+0.7:pointStep(2):grid_dimensions(2)/2+0.05;
    y_gridNeg=-1*fliplr(y_gridPos);
    y_gridPoints=[y_gridNeg  y_gridPos];
    
    tx=[(grid_dimensions(1)-BLnews)/2; 0; hTx];
    rx=[tx(1)+BLnews; 0; hRx];
    
else
    if pointStep(2)>=2*blade_length
        y_gridPoints=-grid_dimensions(2)/2:pointStep(2):grid_dimensions(2)/2;
        tx=[(grid_dimensions(1)-BLnews)/2; pointStep(2)/2+0.05; hTx];
        rx=[tx(1)+BLnews; pointStep(2)/2+0.05; hRx];
    end
end

grid_points=zeros(length(x_gridPoints)*length(y_gridPoints),2);
grid_points(:,1)=kron(x_gridPoints.',ones(length(y_gridPoints),1));
grid_points(:,2)=repmat(y_gridPoints.',length(x_gridPoints),1);



% fan
rng('default')
rng(1)
hf_random=hf_mean+hf_std*rand(size(grid_points,1),1);
% hf_random=hf_mean*ones(size(grid_points,1),1);
fan_origin=[grid_points(:,1) grid_points(:,2) hf_random];   %place the fan origin alternatively in all the grid's points at heigth 1 ---> pivot on which the fan was placed
% theta_fan=0*ones(1,length(t));                                              %fan inclination with respect to the baseline: =0 when the fan is orthogonal to the baseline and defined clockwise
size(fan_origin,1)

%point where the scatter is supposed to be (center of the blade?)
%  Ground parameters for TRP model
epsR=7;    %Relative dielectric constant
sigma_ground=1.5e-2; %Conductivity


t=0:1/sRate:t_oss;
% Matrices initialization
Doppler_signature=zeros(1,length(t));
% maskDS=zeros(1,length(t));
clutterPowerMatrix=zeros(length(y_gridPoints),length(x_gridPoints));
% clutterPowerMatrixTRESH=zeros(length(y_gridPoints),length(x_gridPoints));
% sumDistancesMatrix=zeros(length(x_gridPoints),length(y_gridPoints));

rng('default')
rng(1)
tic
for fo=1:size(fan_origin,1)
    
    if mod(fo,length(y_gridPoints))==0
        rr=length(y_gridPoints);
        cc=fix(fo/length(y_gridPoints));
    else
        rr= mod(fo,length(y_gridPoints));
        cc=1+fix(fo/length(y_gridPoints));
    end
    %      grid(rr,cc)=fan_origin(fo,3);
    %     fan_pos=kron(fan_origin(fo,:).',ones(1,length(t)))+ blade1_rotation;
    
    %% Blades rotation in time+ fan inclination wrt the baseline
    theta_fan=pi.*rand(1)*ones(1,length(t));                                              %fan inclination with respect to the baseline: =0 when the fan is orthogonal to the baseline and defined clockwise
    period=Tmean+Tstd*rand(1);                                                                   % blade rotation period
    phase0=2*pi.*rand(3,1);
    A=Amean+Astd.*randn(1);
    %         A=Amean;
    %         period=Tmean;
    x_fan=[A*cos(2*pi*t/period+phase0(1)).*sin(theta_fan) ;  A*cos(2*pi*t/period+phase0(2)).*sin(theta_fan) ; A*cos(2*pi*t/period+phase0(3)).*sin(theta_fan)];
    y_fan=[A*cos(2*pi*t/period+phase0(1)).*cos(theta_fan) ;  A*cos(2*pi*t/period+phase0(2)).*cos(theta_fan) ; A*cos(2*pi*t/period+phase0(3)).*cos(theta_fan)];
    z_fan=[A*sin(2*pi*t/period+phase0(1)) ; A*sin(2*pi*t/period+phase0(2)) ; A*sin(2*pi*t/period+phase0(3))];
    
    blade1_rotation=[x_fan(1,:); y_fan(1,:) ;z_fan(1,:)];
    %     blade2_rotation=[x_fan(2,:); y_fan(2,:) ;z_fan(2,:)];
    %     blade3_rotation=[x_fan(3,:); y_fan(3,:) ;z_fan(3,:)];
    blades_rotation=[blade1_rotation];%; blade2_rotation; blade3_rotation];
    
    FSCS=0.1;
    %% ************* Doppler signature TRP calculation ************************
    fan_pos=kron(fan_origin(fo,:).',ones(1,length(t)))+ blades_rotation;
    
    Doppler_signature1=TRP_def(rx,tx,fan_pos,f_c,FSCS,epsR,sigma_ground,t_oss(j),sRate);%TRP
    Doppler_signature=Doppler_signature+Doppler_signature1;
    
    
    clutterPowerMatrix(rr,cc)=10*log10(max(abs(Doppler_signature1))^2);%/2
    
end
Power_clutter=10*log10((max(abs(Doppler_signature))^2));
%%
toc
% f_l_p=2/period;
% [z1, p1, k1] = butter(9,f_l_p*2/(sRate),'low');                            % Hardware LPF of 5 order - filter coefficents
% [SOS1, gain1] = zp2sos(z1, p1, k1);                                        % Convert coefficients to SOS form
% h_LP = dfilt.df2tsos(SOS1, gain1);
% fvtool(h_LP)
%
% filtered_DS = filter(h_LP,Doppler_signature);
% % filtered_maskDS = filter(h_LP,maskDS);
%
% fft_FilteredDS=fftshift(fft(filtered_DS));
% fft_FilteredMaskDS=fftshift(fft(filtered_maskDS));

%% PSD
[f,h]=pwelch(real(Doppler_signature),[],[],[],sRate,'onesided');%
% [fm hm]=pwelch(real(maskDS),[],[],[],sRate,'onesided');

%  [f_lpf, h_lpf]=pwelch(real(filtered_DS),[],[],[],sRate,'onesided');%
% [fm_lpf hm_lpf]=pwelch(real(filtered_maskDS),[],[],[],sRate,'onesided');%

figure()
semilogx(h,10*log10(f),'k','LineWidth',2);
% hold on
% semilogx(hm,10*log10(fm),'b--','LineWidth',2);
xlabel('Frequency, Hz', 'FontSize', 12);
ylabel('PSD, dB/Hz', 'FontSize', 12);
grid on

xlim([min(h) max(h)])
set(gca,'FontSize',12);
% legend('no threshold', 'threshold')
% Doppler_signature=Doppler_signature-mean(Doppler_signature);
%%Plots
% Doppler signature
figure()
plot(t,(real(Doppler_signature)),'k','LineWidth',2)
% hold on
% plot(t,(real(maskDS)),'b--','LineWidth',2)
grid on
xlim([min(t) max(t)])
xlabel('Time, s', 'FontSize', 12)
ylabel('Received amplitude, Volt', 'FontSize', 12)
set(gca,'FontSize',12);
% legend('no threshold', 'threshold')

% figure()
% plot(t,(real(filtered_DS)),'k','LineWidth',2)
% % hold on
% % plot(t,(real(filtered_maskDS)),'b--','LineWidth',2)
% grid on
% xlim([min(t) max(t)])
% xlabel('Time, s', 'FontSize', 12)
% ylabel('Received amplitude, Volt', 'FontSize', 12)
% set(gca,'FontSize',12);
% legend('filtered no threshold', 'filtered threshold')
%
% fft_DS=fftshift(fft(Doppler_signature));
% % fft_maskDS=fftshift(fft(maskDS));
% fr=-sRate/2:sRate/length(Doppler_signature):sRate/2-sRate/length(Doppler_signature);
%
% figure
% plot(fr,10*log10(abs(fft_DS)))
% hold on
% % plot(fr,10*log10(abs(fft_maskDS)))
% grid on
% % axis([-60 60 -80 0])
% % xlim([-200 200])
% xlabel('Frequency [Hz]')
% ylabel('Spectrum [dB]')


% %Spectrogram
% Nhamming=256*2;
% figure()
% [S,T,F,P] = twoSidedSpectrogram(real(Doppler_signature),sRate,Nhamming,Nhamming*10,Nhamming-20);
% surf(T+t(1),F,10*log10(P/max(max(P))),'edgecolor','none');axis tight;
% view(0,90);
% ylim([0 10])
% caxis([-50 0])
% xlabel('Time, s', 'FontSize', 12)
% ylabel('Doppler Frequency, Hz', 'FontSize', 12);
% colorbar
% set(gca,'FontSize',12);
% %


%
%
% a=abs(Doppler_signature);
% % b=abs(maskDS);

% figure
% histfit(a,[],'rayleigh')%rician
% xlabel('Power value')
% ylabel('N repetition')
% title('Doppler Signature Power')
%
% % figure
% % histfit(b,[],'rayleigh')%rician
% % title(['f=' num2str(freq/1e9) 'GHz, threshold=' num2str(treshold) ])
%
% figure
% histfit(real(Doppler_signature))
% xlabel('Amplitude value')
% ylabel('N repetition')
% title('I component')

% figure
% histfit(imag(Doppler_signature))
% xlabel('Amplitude value')
% ylabel('N repetition')
% title('Q component')
%
% figure
% histfit(real(filtered_DS))
% xlabel('Amplitude value')
% ylabel('N repetition')
% title('I component')
%
% figure
% histfit(imag(filtered_DS))
% xlabel('Amplitude value')
% ylabel('N repetition')
% title('Q component')


% if hTx==hRx
%     maxCP=max(max(clutterPowerMatrix));
%     save maxCP maxCP
% end

figure
surf(x_gridPoints,y_gridPoints,(clutterPowerMatrix),'Linestyle', 'none')%-maxCP
xlabel('x space [m]')
ylabel('y space [m]')
zlabel('Power [dB]')
ylim([min(y_gridPoints) max(y_gridPoints)])
xlim([min(x_gridPoints) max(x_gridPoints)])
title(['f=' num2str(f_c/1e9) 'GHz, [l_x , l_y]=[' num2str(pointStep(1)) ' , ' num2str(pointStep(2)) ']'])
% % title(['h_T_x=' num2str(hTx) 'm, h_R_x=' num2str(hRx) 'm, h_F_A_N= ' num2str(hf) 'm'])
caxis([max(max(clutterPowerMatrix))-50 max(max(clutterPowerMatrix))])%-maxCP
% % grid off
% % shading INTERP%FLAT


% figure
% surf(x_gridPoints,y_gridPoints,(clutterPowerMatrixTRESH),'Linestyle', 'none')%-maxCP
% xlabel('x space [m]')
% ylabel('y space [m]')
% zlabel('Power [dB]')
% ylim([min(y_gridPoints) max(y_gridPoints)])
% xlim([min(x_gridPoints) max(x_gridPoints)])
% title(['f=' num2str(freq/1e9) 'GHz, threshold=' num2str(treshold) ])
% title(['h_T_x=' num2str(hTx) 'm, h_R_x=' num2str(hRx) 'm, h_F_A_N= ' num2str(hf) 'm'])
% caxis([max(max(clutterPowerMatrix))-50 max(max(clutterPowerMatrix))])%-maxCP
%
%
% figure(300)
% hold on
% plot(x_gridPoints,clutterPowerMatrix(length(y_gridPoints)/2,:),'b','LineWidth',2)
% xlabel('x space [m]')
% hold on
% plot(x_gridPoints, treshold*ones(length(x_gridPoints)),'r--','LineWidth',2)
% ylabel('Power [dB]')
% grid on
% legend (['h_T_x=' num2str(hTx) ])
% legend (['h_F_A_N= ' num2str(hf) ])
% legend (['h_R_x=' num2str(hRx) ])

% figure
% plot(grid_points(:,1),grid_points(:,2),'c.','Markersize',10)
% hold on
% plot(tx(1),tx(2),'rs','Markersize',8)
% plot(rx(1),rx(2),'gs','Markersize',8)
% grid on
% xlabel('x')
% ylabel('y')
% title('Bush origin grid')
% ylim([min(y_gridPoints) max(y_gridPoints)])
% title(['[l_x , l_y]=[' num2str(pointStep(1)) ' , ' num2str(pointStep(2)) ']'])


% figure
% plot(x_gridPoints,y_gridPoints,(grid))
% hold on
% plot3(tx(1),tx(2),tx(3),'ks','Markersize',8)
% plot3(rx(1),rx(2),rx(3),'ko','Markersize',8)
% xlabel('x space [m]')
% ylabel('y space [m]')
% zlabel('grid ')

% figure
% C=contour(x_gridPoints,y_gridPoints,clutterPowerMatrix,[-70:2:-36 max(max(clutterPowerMatrix))]);
% clabel(C)%,'manual'

% figure
% % surf(x_gridPoints,y_gridPoints,(FSCS))
% % xlabel('x space [m]')
% % ylabel('y space [m]')
% % zlabel('FSCS [dB]')


for k=1: n
    [dTx_Int(k),int1(k,:),tx]=mainshell_uwb(BL, hTx,taop(k));
    theta=0:2*pi/1000:2*pi;
    x=real(int1(k,1)*cos(theta));
    y=real(int1(k,2)*sin(theta));
    nor=fix((x+(2*BLnews)/2)/0.6)*0.6; % N.B> to do check the area around the baseline
    nor=nor';
    for j=1:size(y,2)
        if y(j)>=0
            if y(j)<=0.55
                nor(j,2)=0.55;
            else
                nor(j,2)=floor(y(j)/0.6)*0.6+0.55;
            end
        else
            if y(j)>=-0.55
                nor(j,2)=-0.55;
            else
                nor(j,2)=floor(y(j)/0.6)*0.6-0.55;
            end
        end
        
    end
    nor=unique(nor,'rows');
    nor1=unique(nor(:,2));
    format long
    CP_ellips=-250*ones(size(clutterPowerMatrix,1),size(clutterPowerMatrix,2));
    
    for ii=2:size(nor1,1)/2+1
        [r, ~] = find(nor(:,2)==nor1(ii));
        na=nor(r,1);
        minna=double(min(na));
        maxna=double(max(na));
        cc1=find(abs(x_gridPoints-minna)<10^-6);
        cc2=find(abs(x_gridPoints-maxna)<10^-6);
        cc3=find(abs(y_gridPoints-nor1(ii))<10^-6);
        
        CP_ellips(cc3,cc1:cc2)=clutterPowerMatrix(cc3,cc1:cc2);
        CP_ellips(((size(clutterPowerMatrix,1))-cc3),cc1:cc2)=clutterPowerMatrix(((size(clutterPowerMatrix,1))-cc3),cc1:cc2);
        
        clear minna maxna
    end
    clear x y
    
    %
    PowerC_ellips(1,k)=10*log10(sum(sum(10.^(clutterPowerMatrix./10)))); %total power
    PowerC_ellips(2,k)=10*log10(sum(10.^(nonzeros(CP_ellips)./10)));%power of section
    
    
    % PowerC_ellips(3,k)=max(max(clutterPowerMatrix));
    % PowerC_ellips(4,k)=(sum(sum(10.^((clutterPowerMatrix-max(max(clutterPowerMatrix)))./10))));
    % PowerC_ellips(5,k)=sum(10.^((nonzeros(CP_ellips)-max(max(clutterPowerMatrix)))./10));
    [kk,qq]=find(clutterPowerMatrix==max(max(clutterPowerMatrix)));
    maxCPM= clutterPowerMatrix(kk,qq);
    fig3=figure;
    surf(x_gridPoints,y_gridPoints,(CP_ellips-max(max(clutterPowerMatrix))),'Linestyle', 'none')%
    caxis( [min(min(clutterPowerMatrix-maxCPM)) 0])
    
    set(gca,'FontSize',26)
    xlabel('x space [m]')
    ylabel('y space [m]')
    zlabel('[dB]')
    ylim([min(y_gridPoints) max(y_gridPoints)])
    xlim([min(x_gridPoints) max(x_gridPoints)])
    x=max(x_gridPoints)/2+int1(k,1)*cos(theta);
    y=int1(k,2)*sin(theta);
    hold on
    plot(x,y,'r','LineWidth',2);
    clear x y
    hold off
    title(['BL=' num2str(BL) 'm, h=' num2str(hTx) 'm,\tau =' num2str(taop(k)/1e-9) 'ns' ])
    colorbar
    colormap jet
    set(fig3,'units','pixel');
    set(fig3,'position',[0,0,960,760]);
    savefig(fig3,sprintf('Clutter%d%d',BL,k))
    
    
end

%     fig2=figure;
%     plot(taop/1e-9,dTx_Int,'bo-','LineWidth',2);
%     set(gca,'FontSize',26)
%     set(fig2,'units','pixel');
%     % set(fig2,'position',[0,0,2000,960]);
%     grid on
%     title (['Hidden space vs Pulse Width BL=' num2str(BL)])
%     xlabel('Pulse Width [ns]')
%     ylabel('Space [m]')
%     savefig(fig2,sprintf('HivsI%d',BL))

save PowerC_ellips
save dtX_Int
% % end

figure
dif=Power_tg-PowerC_ellips(2,:);
plot(taop,dif-min(dif))
grid on