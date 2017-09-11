clear
clc
close all
%%
f_c=7.36e9;
c=3e8;
sRate=200;
hRx=1;
hTx=1;
BL=100;
lambda=c/f_c;

extragrid=40;
grid_dimensions=[BL+extragrid 100];
tx=[(grid_dimensions(1)-BL)/2; 0; hTx];
rx=[tx(1)+BL; 0; hRx];
v=[0 2 0]; %[m/s]
% x_tg=grid_dimensions(1)/2*0.6;
 x_tg=rx(1)-1;
z_tg=1;
% taop=[0.5e-9 1e-9 2e-9 3e-9 4e-9 5e-9 6e-9 7e-9 8e-9 9e-9 10e-9];
B=[3e9 2.5e9 2e9 1.5e9 1e9 0.5e9 0.4e9 0.3e9 0.2e9 0.1e9];
taop=1./B;
% taop=[1e-10 2 3e-10 6e-10 1e-9 3e-9 6e-9 10e-9 20e-9];
%% mainshell
n=1;
RCS_tg=10;
[int1]=mainshell_uwb(BL, hTx,taop);
t_end=(int1(:,2)*2)/v(2);
y_tg=-int1(:,2);
Power_tg=zeros(1,size(taop,2));
epsR=7;
sigma_ground=1.5e-2;
ref=0;
Tr=10^-6;
N=t_end./Tr;
% N=200;
ptz_db=zeros(size(taop,2),2);
ptz=ptz_db;
for i=1:length(taop)
    tg0=[x_tg,y_tg(i),z_tg];
    tt=0:1/sRate:t_end(i);
    tg_inc=(v'*tt)';                                                          %receiver position increment
    pos_tg=[tg0(1)+tg_inc(:,1) tg0(2)+tg_inc(:,2) tg0(3)+tg_inc(:,3)];
    %     tx_tg=sqrt((tx(1)-pos_tg(:,1)).^2+(tx(2)-pos_tg(:,2)).^2+(tx(3)-pos_tg(:,3)).^2);
    %     rx_tg=sqrt((rx(1)-pos_tg(:,1)).^2+(rx(2)-pos_tg(:,2)).^2+(rx(3)-pos_tg(:,3)).^2);
    %     Power_tg_ist=(4*pi*RCS_tg*rx(3)^2*tx(3)^2*pos_tg(3)^4)./(lambda^2*(tx_tg.^4).*(rx_tg.^4));
    %     Power_tg(i)= 10*log10(sum(Power_tg_ist));
    
    
    [u_tg]=TRP_def(rx,tx, pos_tg',f_c,RCS_tg,epsR,sigma_ground,t_end(i),sRate,ref);
%     figure
%     plot(tt, u_tg)
%     figure
%         plot(pos_tg(:,1) , pos_tg(:,2))
%         figure
%     plot(tt,(abs(u_tg)).^2)
    ptz(i,1)=(sum(abs(u_tg))^2)*taop(i)/Tr;
    %  ptz(i,1)=max(abs(u_tg))^2*taop(i)*N(i);
    ptz(i,2)=(sum(abs(u_tg))^2);
    ptz_db(i,1)=10*log10(ptz(i,1));
    ptz_db(i,2)=10*log10(ptz(i,2));
end
%% clutter
%
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

% nRow=grid_dimensions(2)/2/pointStep(2);
y_gridPos=blade_length+0.05:pointStep(2):grid_dimensions(2)/2+0.05;
y_gridNeg=-1*fliplr(y_gridPos);
y_gridPoints=[y_gridNeg  y_gridPos];
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
%% Ground parameters for TRP model
epsR=7;    %Relative dielectric constant
sigma_ground=1.5e-2; %Conductivity

theta=0:2*pi/1000:2*pi;
PowerC_ellips=zeros(2,size(taop,2));
for j=1:length(taop)
    t=0:1/sRate:t_end(j);
    %% Matrices initialization
    Doppler_signature=zeros(1,length(t));
    clutterPowerMatrix=zeros(length(y_gridPoints),length(x_gridPoints));
    
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
        ref=1;
        Doppler_signature1=TRP_def(rx,tx,fan_pos,f_c,FSCS,epsR,sigma_ground,t_end(j),sRate,ref);%TRP
        %         Doppler_signature=Doppler_signature+Doppler_signature1;
        %
        %         Power_clutter=10*log10((max(abs(Doppler_signature))^2));
   clutterPowerMatrix(rr,cc)=10*log10((max(abs(Doppler_signature1)))^2*(taop(i)*N(i))/Tr);%/2
%          clutterPowerMatrix(rr,cc)=10*log10(max(abs(Doppler_signature1))^2);%/2
        
    
    end
    [kk,qq]=find(clutterPowerMatrix==max(max(clutterPowerMatrix)));
    maxCPM= clutterPowerMatrix(kk,qq);
    
    n=1;
    a=(BL+n*c*taop(j))/2; %Semi-axes of ellipsoid
    b=sqrt(a.^2- (BL/2).^2);
    if b^2<hTx^2
        error('The main shell is not defined');
    end
    
    x=real(int1(j,1)*cos(theta));
    y=real(int1(j,2)*sin(theta));
    nor=fix((x+(BL+extragrid)/2)/0.6)*0.6; % N.B> to do check the area around the baseline
    nor=nor';
    for k=1:size(y,2)
        if y(k)>=0
            if y(k)<=0.55
                nor(k,2)=0.55;
            else
                nor(k,2)=floor(y(k)/0.6)*0.6+0.55;
            end
        else
            if y(k)>=-0.55
                nor(k,2)=-0.55;
            else
                nor(k,2)=floor(y(k)/0.6)*0.6-0.55;
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
    PowerC_ellips(1,j)=10*log10(sum(sum(10.^(clutterPowerMatrix./10)))); %total power
    PowerC_ellips(2,j)=10*log10(sum(10.^(nonzeros(CP_ellips)./10)));%power of section
    
    % PowerC_ellips(3,k)=max(max(clutterPowerMatrix));
    % PowerC_ellips(4,k)=(sum(sum(10.^((clutterPowerMatrix-max(max(clutterPowerMatrix)))./10))));
    % PowerC_ellips(5,k)=sum(10.^((nonzeros(CP_ellips)-max(max(clutterPowerMatrix)))./10));
    
    fig2=figure;
    surf(x_gridPoints,y_gridPoints,(CP_ellips-max(max(clutterPowerMatrix))),'Linestyle', 'none')%
    caxis( [min(min(clutterPowerMatrix-maxCPM)) 0])
    set(gca,'FontSize',26)
    xlabel('x space [m]')
    ylabel('y space [m]')
    zlabel('[dB]')
    ylim([min(y_gridPoints) max(y_gridPoints)])
    xlim([min(x_gridPoints) max(x_gridPoints)])
    x=max(x_gridPoints)/2+int1(j,1)*cos(theta);
    y=int1(j,2)*sin(theta);
    hold on
    plot(x,y,'r','LineWidth',2);
    clear x y
    hold off
    title(['BL=' num2str(BL) 'm, h=' num2str(hTx) 'm,\tau =' num2str(taop(j)/1e-9) 'ns' ])
    colorbar
    colormap jet
    set(fig2,'units','pixel');
    set(fig2,'position',[0,0,960,760]);
    savefig(fig2,sprintf('Clutter%d%d%d',BL,hTx,j))
    
end
fig1=figure;
surf(x_gridPoints,y_gridPoints,(clutterPowerMatrix-max(max(clutterPowerMatrix))),'Linestyle', 'none')%
caxis( [min(min(clutterPowerMatrix-maxCPM)) 0])
colorbar
colormap jet
set(gca,'FontSize',18)
xlabel('x space [m]')
ylabel('y space [m]')
zlabel('[Power_d_B]')
% set(fig1,'units','pixel');
% set(fig1,'position',[0,0,960,760]);
% savefig(fig1,sprintf('Clutter_tot%d%d',BL,hTx))




%%
dif=((10.^(PowerC_ellips(1,:)./10)-10.^(PowerC_ellips(2,:)./10))./10.^(PowerC_ellips(1,:)./10)).*100;


figure
plot(taop/1e-9,dif)
grid on
xlabel('Pulse Width [ns]')
ylabel('Clutter Reduction_%')


hspace=grid_dimensions(1)/2-int1(:,1)-tx(1);
figure
plot(taop/1e-9,(10.^(PowerC_ellips(2,:)./10))./10.^(PowerC_ellips(1,:)./10)*100)
hold on
grid on
% title('% power delate on total power ')% figure
xlabel('Pulse Width [ns]')
ylabel('Clutter_%')
yyaxis right

plot(taop/1e-9,hspace)
hold on
ylabel('Hidden spice [m]')
ylim ([0 max(hspace)])

ptz_tg=ptz_db(:,1)+10*log10(N);
SCR=ptz(:,1)./(10.^(PowerC_ellips(2,:)./10))';
figure
plot (taop/1e-9, 10*log10(SCR)+10*log10(N))
grid on
title('SCR')% figure
xlabel('Pulse Width [ns]')
ylabel('SCR_d_B')
hold on
figure
plot (t_end, ptz_db(:,2)-PowerC_ellips(1,:)')
grid on
title('SCR cw')% figure
xlabel('Toss [s]')
ylabel('SCR_d_B')
% xlim ([t_end(1) t_end(end)])

% figure
% plot(taop/1e-9,ptz_db(:,1))
% title('tg')% figure
% grid on
% figure
% plot(t_end,ptz_db(:,2))
% title('tg_cw')% figure
% grid on
%
% figure
% plot(taop,PowerC_ellips(2,:))
% title('cl')% figure
% grid on
% figure
% plot(t_end,PowerC_ellips(1,:))
% title('cl_cw')% figure
% grid on
%
%
%
%
% % fig2=figure;
% % plot(taop/1e-9,dTx_Int,'bo-','LineWidth',2);
% %   set(gca,'FontSize',26)
% % set(fig2,'units','pixel');
%
% %     title (['Hidden space vs Pulse Width BL=' num2str(BL)])
% %     xlabel('Pulse Width [ns]')
% %     ylabel('Space [m]')
% %     savefig(fig2,sprintf('HivsI%d',BL))


% fig=figure;
% 
% for i=1:size(int1,1)
%      x=max(x_gridPoints)/2+int1(i,1)*cos(theta);
%     y=int1(i,2)*sin(theta);
% 
% hold on
% h(i)=plot(x,y,'b','LineWidth',2);
% end
% plot(tx(1),tx(2),'rd','LineWidth',2,'MarkerSize',10)
% plot(rx(1),rx(2),'rd','LineWidth',2,'MarkerSize',10)
% set(gca,'FontSize',26)
% set(fig,'units','pixel');
% set(fig,'position',[0,0,2000,960]);
% title ('Main Shell')
% xlabel('x [m]')
% ylabel('y [m]')
% grid on
% axis equal
%legend(h(:),{'\tau = 0.20 ns','\tau = 0.36 ns','\tau = 0.52 ns','\tau = 0.68 ns','\tau = 0.84 ns','\tau = 1.00 ns'});  
% legend(h(:),{'\tau min','\tau2','\tau3','\tau4','\tau max '});  


