clear
clc
close all

f_c=1e9;
c=3e8;
sRate=200;
hRx=1;   %transmitter and receiver height
hTx=1;
BL=10;
lambda=c/f_c;

step=0.6;
rest=floor(BL/step); 
BLnews=BL;
grid_dimensions=[2*BLnews BLnews]; %[horizontal vertical] grid dimentions to define the area to evaluate the clutter returns
      
tx=[(grid_dimensions(1)-BLnews)/2; 0; hTx];
rx=[tx(1)+BLnews ; 0; hRx];
v=[0 2 0]; %[m/s]
x_tg=grid_dimensions(1)/2;
z_tg=1;
% taop=[0.5e-9 1e-9 2e-9 3e-9 4e-9 5e-9 6e-9 7e-9 8e-9 9e-9 10e-9];

%% tg
RCS_tg=-90;

    t_end=grid_dimensions(2)/v(2);
    tt=0:1/sRate:t_end;
    y_tg=-grid_dimensions(2)/2;
    tg0=[x_tg,y_tg,z_tg];
    tg_inc=(v'*tt)';                                                          %receiver position increment
    pos_tg=[tg0(1)+tg_inc(:,1) tg0(2)+tg_inc(:,2) tg0(3)+tg_inc(:,3)];
    tx_tg=sqrt((tx(1)-pos_tg(:,1)).^2+(tx(2)-pos_tg(:,2)).^2+(tx(3)-pos_tg(:,3)).^2);
    rx_tg=sqrt((rx(1)-pos_tg(:,1)).^2+(rx(2)-pos_tg(:,2)).^2+(rx(3)-pos_tg(:,3)).^2);
    Doppler_signature_tg=-sin(2*pi*f_c/c*(tx_tg+rx_tg-BLnews ));   
    for ti=1:1:length(tt)
       Power_tg_ist(ti)=(4*pi*RCS_tg*rx(3)^2*tx(3)^2*pos_tg(3)^4)/(lambda^2*tx_tg(ti)^4*rx_tg(ti)^4);
    end
%     Power_tg(j)=10*log10(max(abs(Doppler_signature_tg).^2));
%     Power_tg_tot(j)=10*log10(sum(10.^(nonzeros(Power_tg)./10)));
    Power_tg_ist_db=10*log10(Power_tg_ist);
    Power_tg= 10*log10(sum(Power_tg_ist));

    wind=window(@hamming,length(tt));
    Doppler_signature_tg_w=wind.*Doppler_signature_tg;
       
    figure()
    plot(tt,Doppler_signature_tg_w);
    xlabel('Time, s')
    ylabel('Amplitude, V')
    xlim([min(tt) max(tt)])
    
    % Power Spectral Density
%     [f, h]=pwelch(real(Doppler_signature_tg),[],[],[],sRate,'onesided');
%     figure()
%     semilogx(h,10*log10(f),'k','LineWidth',2);
%     xlabel('Frequency, Hz', 'FontSize', 12);
%     ylabel('PSD, dB/Hz', 'FontSize', 12);
%     grid on
%     xlim([min(h) max(h)])
%     set(gca,'FontSize',12);





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


    t=0:1/sRate:t_end;
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
        
        Doppler_signature1=TRP_def(rx,tx,fan_pos,f_c,FSCS,epsR,sigma_ground,t_end,sRate);%TRP
        Doppler_signature=Doppler_signature+Doppler_signature1;

       
        clutterPowerMatrix(rr,cc)=10*log10(max(abs(Doppler_signature1))^2);%/2
       
    end
     Power_clutter=10*log10((max(abs(Doppler_signature))^2));
      
    sgn_cl= Doppler_signature+Doppler_signature_tg';
     
    figure()
    plot(tt,Doppler_signature);
    xlabel('Time, s')
    ylabel('Amplitude, V')
    xlim([min(tt) max(tt)])
    grid on

     figure()
    plot(tt,sgn_cl);
    xlabel('Time, s')
    ylabel('Amplitude, V')
    xlim([min(tt) max(tt)])
    grid on
    
    figure
    plot(tt, Power_tg_ist_db-Power_clutter)
    grid on

