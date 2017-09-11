function [int1]=mainshell_uwb(BL,hTx,taop)
c = 3e8;
n=1; %main shell

tx=[BL, 0,hTx]; %position transmiter

% C=1; B=0; A=0;D=1;
% xx = [abs(Tx(1))+10 -abs(Tx(1))-10 -abs(Tx(1))-10 abs(Tx(1))+10]; % Generate data for x vertices
% yy= [abs(Tx(1))+10 abs(Tx(1))+10 -abs(Tx(1))-10 -abs(Tx(1))-10]; % Generate data for y vertices
% zz = -1/C*(A*xx + B*yy + D); % Solve
% taopmin=taop(1);
% % taopmax=taop(end);
% for i =1 :size(taop,2)
    a=(BL+n*c*taop)/2; %Semi-axes of ellipsoid
    b=sqrt(a.^2- (BL/2).^2);
%     [theta,phi]=meshgrid(0:2*pi/180:2*pi,0:2*pi/90:pi);
%     x=a*cos(theta).*sin(phi);
%     y=b.*sin(theta).*sin(phi);
%     z=1+b.*cos(phi);
    
    int1(:,1)=sqrt((1-tx(3)^2./b.^2).*a.^2);%x intersection
    int1(:,2)=sqrt((b.^2)-tx(3)^2);%y intersection
    
%     dTx_Int= BL/2-int1;
  
    
%     if taop==taopmin %||taop==taopmax
%         a=taopmin/10^-9;
%             if taop==taopmax/10^-9
%                 a=taopmax/10^-9;
%             else
%                 a=taopmin/10^-9;
%             end
%         
%         
%         fig=figure;
%         surf(x,y,z); axis equal;
%         hold on
%         patch(xx, yy, zz,'k','LineWidth',2);
%         xlabel('x space [m]')
%         ylabel('y space [m]')
%         taop=taop/1e-9;
%         title (['Itersection Main Shell-Ground BL=' num2str(BL) 'm, \tau =' num2str(taop(i)) 'ns'])
%         savefig(fig,sprintf('MainShell%d%',BL))
%                  saveas(fig,sprintf('Mainshell%d%d.jpeg',a,BL))
%         
%         
%     end
   
end



