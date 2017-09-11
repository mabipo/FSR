function [clutterPowerMatrix,y_gridPoints,x_gridPoints,kk,qq,maxCPM,Doppler_signature]=Clutter_pisa(hTx,freq,BL,t_end)
%[clutterPowerMatrix,Doppler_signature1,y_gridPoints,x_gridPoints]=Clutter(hTx,freq,BL,taop,taopmin)


%% Parameters
sRate = 200;
t=0:1/sRate:t_end;

%Blade
blade_length=0.5;
period=1;
%theta_fan=0*ones(1,length(t));   %fan inclination
theta_fan=rand(1)*ones(1,length(t));
phase0=2*pi.*rand(3,1);         %initial phase
Amean=blade_length;%/2
Astd=0.2;

hf=1;  %fan origin height
hRx=hTx;

epsR=7; 
sigma_ground=1.5e-2;

%grid
extragrid=100;
grid_dimensions=[BL+extragrid BL]; %[horizontal vertical] grid dimentions to define the area to evaluate the clutter returns
pointStep=[1.0 1.0]*0.6;            %grid horizontal and vertical step
x_gridPoints=0:pointStep(1):grid_dimensions(1);

if pointStep(2)<2*blade_length
    y_gridPos=blade_length+0.05:pointStep(2):grid_dimensions(2)/2+0.05;
    y_gridNeg=-1*fliplr(y_gridPos);
    y_gridPoints=[y_gridNeg  y_gridPos];
    tx=[(grid_dimensions(1)-BL)/2; 0; hTx];
    rx=[tx(1)+BL; 0; hRx];
    
else if pointStep(2)>=2*blade_length
        y_gridPoints=-grid_dimensions(2)/2:pointStep(2):grid_dimensions(2)/2;
        tx=[(grid_dimensions(1)-BL)/2; pointStep(2)/2+0.05; hTx];
        rx=[tx(1)+BL; pointStep(2)/2+0.05; hRx];
    end
end


grid_points=zeros(length(x_gridPoints)*length(y_gridPoints),2);
grid_points(:,1)=kron(x_gridPoints.',ones(length(y_gridPoints),1));
grid_points(:,2)=repmat(y_gridPoints.',length(x_gridPoints),1);

% fan
fan_origin=[grid_points(:,1) grid_points(:,2) hf*ones(size(grid_points,1),1)];   %place the fan origin alternatively in all the grid's points at heigth 1 ---> pivot on which the fan was placed
% theta_fan=0*ones(1,length(t));                                              %fan inclination with respect to the baseline: =0 when the fan is orthogonal to the baseline and defined clockwise
size(fan_origin,1)
%point where the scatter is supposed to be (center of the blade?)


%% Matrices initialization
clutterPowerMatrix=zeros(length(y_gridPoints),length(x_gridPoints));
Doppler_signature=zeros(1,length(t));

% sumDistancesMatrix=zeros(length(x_gridPoints),length(y_gridPoints));
rng('default')
rng(1)

for fo=1:size(fan_origin,1)
    
    if mod(fo,length(y_gridPoints))==0
        rr=length(y_gridPoints);
        cc=fix(fo/length(y_gridPoints));
    else
        rr= mod(fo,length(y_gridPoints));
        cc=1+fix(fo/length(y_gridPoints));
    end
    
    %% Blades rotation in time+ fan inclination wrt the baseline

    
     A=Amean+Astd.*randn(1);
    
    x_fan=[A*cos(2*pi*t/period+phase0(1)).*sin(theta_fan) ;  A*cos(2*pi*t/period+phase0(2)).*sin(theta_fan) ; A*cos(2*pi*t/period+phase0(3)).*sin(theta_fan)];
    y_fan=[A*cos(2*pi*t/period+phase0(1)).*cos(theta_fan) ;  A*cos(2*pi*t/period+phase0(2)).*cos(theta_fan) ; A*cos(2*pi*t/period+phase0(3)).*cos(theta_fan)];
    z_fan=[A*sin(2*pi*t/period+phase0(1)) ; A*sin(2*pi*t/period+phase0(2)) ; A*sin(2*pi*t/period+phase0(3))];
    
    blade1_rotation=[x_fan(1,:); y_fan(1,:) ;z_fan(1,:)];
    %     blade2_rotation=[x_fan(2,:); y_fan(2,:) ;z_fan(2,:)];
    %     blade3_rotation=[x_fan(3,:); y_fan(3,:) ;z_fan(3,:)];
    blades_rotation=[blade1_rotation];%; blade2_rotation; blade3_rotation];
%     
    FSCS=0.25; %da dati
   
%     %

    %% ************* Doppler signature TRP calculation ************************
    fan_pos=kron(fan_origin(fo,:).',ones(1,length(t)))+ blades_rotation;
    Doppler_signature1=TRP_def(rx,tx,fan_pos,freq,FSCS,epsR,sigma_ground,t_end,sRate) ;%TRP
    
 Doppler_signature=Doppler_signature+Doppler_signature1; %somma delle
%     % doppler per ogni posizione della griglia
    clutterPowerMatrix(rr,cc)=10*log10(max(abs(Doppler_signature1))^2);
   
%     
 end
[kk,qq]=find(clutterPowerMatrix==max(max(clutterPowerMatrix)));
maxCPM= clutterPowerMatrix(kk,qq);

    fig=figure;
    surf(x_gridPoints,y_gridPoints,(clutterPowerMatrix-unique(maxCPM)),'Linestyle', 'none')%
%       surf(x_gridPoints,y_gridPoints,(clutterPowerMatrix),'Linestyle', 'none')%
    set(gca,'FontSize',26)
    xlabel('x space [m]')
    ylabel('y space [m]')
    zlabel('[dB]')
    ylim([min(y_gridPoints) max(y_gridPoints)])
    xlim([min(x_gridPoints) max(x_gridPoints)])
    caxis( [min(min(clutterPowerMatrix-maxCPM)) 0])
%        caxis( [min(min(clutterPowerMatrix)) 0])
    title(['BL=' num2str(BL) 'm, h=' num2str(hTx) 'm' ])
    colorbar
    colormap jet
    set(fig,'units','pixel');
    set(fig,'position',[0,0,960,760]);

  
%  end

end



