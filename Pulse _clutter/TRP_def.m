function [u_tg]=TRP_def(rx_position,tx_position,target_position,freq,FSCS,epsR,sigma_ground,t_end,sRate,ref)

%% inizialization parameters
c=3e8;
lambda=c/freq;


t=0:1/sRate:t_end;
eps0=8.85e-12; %Dielectric constant

epsG=epsR-1j*(sigma_ground/(2*pi*freq*eps0)); %Complex relative dielectric permittivity of the ground 

baseline=sqrt((rx_position(1,:)-tx_position(1,:)).^2+(rx_position(2,:)-tx_position(2,:)).^2+(rx_position(3,:)-tx_position(3,:)).^2);
alpha2=atan((rx_position(3,:)+tx_position(3,:))./baseline); %Grazing angle 

gammaV_2=(epsG*sin(alpha2)-sqrt(epsG-cos(alpha2).^2))./(epsG*sin(alpha2)+sqrt(epsG-cos(alpha2).^2)); %Complex ground reflection coefficient 

%gammaH_2=(sin(alpha2)-sqrt(epsG-cos(alpha2).^2))./(sin(alpha2)+sqrt(epsG-cos(alpha2).^2));

%% Leakage signal
R2=sqrt(baseline.^2+(rx_position(3,:)+tx_position(3,:)).^2); %Reflected path from the ground Tx-RX
R1=sqrt(baseline^2+(rx_position(3,:)-tx_position(3,:))^2); % direct ray TR-RX
R1=repmat(R1,1,size(t,2));
R2=repmat(R2,1,size(t,2));
U1=lambda./(4*pi*R1); %Free-space loss direct ray
U2=lambda./(4*pi*R2); %Free-space loss refleceted ray
phi1=2*pi*R1./lambda; %phase direct ray
phi2=2*pi*R2./lambda;%phase refleceted ray

ul=U1.*exp(1j*phi1)+gammaV_2.*U2.*exp(1j*phi2); %Direct Leakage signal 
% 
% ul= ul*10^2.4; %guadagno calcolato
% %% Target signal
R3=sqrt((rx_position(1,:)-target_position(1,:)).^2 + (rx_position(2,:)-target_position(2,:)).^2 + (rx_position(3,:)-target_position(3,:)).^2 ); %Direct path Tx-target
R4=sqrt((rx_position(1,:)-target_position(1,:)).^2 + (rx_position(2,:)-target_position(2,:)).^2 + (rx_position(3,:)+target_position(3,:)).^2 ); %Reflected path from the ground Tx-target
R5=sqrt((tx_position(1,:)-target_position(1,:)).^2 + (tx_position(2,:)-target_position(2,:)).^2 + (tx_position(3,:)-target_position(3,:)).^2 ); %Direct path Rx-target
R6=sqrt((tx_position(1,:)-target_position(1,:)).^2 + (tx_position(2,:)-target_position(2,:)).^2 + (tx_position(3,:)+target_position(3,:)).^2 ); %Reflected path from the ground Rx-target


U3=lambda./(4*pi*R3); %Free-space loss
U4=lambda./(4*pi*R4);
U5=lambda./(4*pi*R5);
U6=lambda./(4*pi*R6);

phi3=2*pi*R3./lambda; %Path phase shift
phi4=2*pi*R4./lambda;
phi5=2*pi*R5./lambda;
phi6=2*pi*R6./lambda;

alpha4=atan((target_position(3,:)+rx_position(3,:))/R4); %Viewing angles 
alpha6=atan((target_position(3,:)+tx_position(3,:))/R6);

gammaV_4=(epsG*sin(alpha4)-sqrt(epsG-cos(alpha4).^2))./(epsG*sin(alpha4)+sqrt(epsG-cos(alpha4).^2)); %Complex ground reflection coefficient 
gammaV_6=(epsG*sin(alpha6)-sqrt(epsG-cos(alpha6).^2))./(epsG*sin(alpha6)+sqrt(epsG-cos(alpha6).^2));
% 
if ref==1
Ltx_tgt=U3.*exp(1j*phi3)+gammaV_4.*U4.*exp(1j*phi4); %propagation loss Tx-target
Ltgt_rx=U5.*exp(1j*phi5)+gammaV_6.*U6.*exp(1j*phi6); %propagation loss target-Rx
else
Ltx_tgt=U3.*exp(1j*phi3);
Ltgt_rx=U5.*exp(1j*phi5);
%  u_tg=FSCS.*Ltx_tgt.*Ltgt_rx;
end
u_tg=(sqrt(4*pi*FSCS)/lambda)*Ltx_tgt.*Ltgt_rx; %target signal
% u_tg=u_tg-mean(u_tg);
% u_tg= u_tg*10^2.4; %guadagno
u_tot=(ul-u_tg); %recived signal



