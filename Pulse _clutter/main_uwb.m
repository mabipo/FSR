clear
close all
clc
%%
% BL=[500, 1000,1500];
BL=100;
freq=7.36e9;
c = 3e8;
hTx=1;
sRate=200;
taop=[0.5e-9 1e-9  2e-9 3e-9 4e-9 5e-9 6e-9 7e-9 8e-9 9e-9 10e-9];
PowerC_ellips=zeros(2,size(taop,2));
n=1;
a=(BL+n*c*taop)/2; %Semi-axes of ellipsoid
b=sqrt(a.^2- (BL/2).^2);
for i =1 :size(taop,2)
if b(i)^2<hTx^2
    error('The main shell is not defined');
end
end
% %
% for i=1:2%size(BL,2)
[clutterPowerMatrix(:,:),y_gridPoints,x_gridPoints,kk,qq,maxCPM,Doppler_signature_clutter]=Clutter_pisa(hTx,freq,BL);

save Doppler_signature_uwb.mat
save clutterPowerMatrix_uwb.mat

%%
extragrid=100;
[dTx_Int,int1]=mainshell_uwb(BL, hTx,taop);
theta=0:2*pi/1000:2*pi;

for k=1: size(taop,2)
x=real(int1(k,1)*cos(theta));
y=real(int1(k,2)*sin(theta));
nor=fix((x+(BL+extragrid)/2)/0.6)*0.6; % N.B> to do check the area around the baseline
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
%% Power Spectral Density
[f, h]=pwelch(real(Doppler_signature_clutter),[],[],[],sRate,'onesided');
figure()
semilogx(h,10*log10(f),'k','LineWidth',2);
xlabel('Frequency, Hz', 'FontSize', 12);
ylabel('PSD, dB/Hz', 'FontSize', 12);
grid on
xlim([min(h) max(h)])
set(gca,'FontSize',12);