clc
clear all
close all
profile on

addpath(genpath('../'));

[StateVar,VARIABLES,DOMAIN,BC,IBM,LS] = setUpVariables;

%% ==============================
% load 'initial at any time '


%% object location and IBM coefficients
% [IBM_coeffU,IBM_coeffV,IBM_coeffP]=IBMcoeffs(IBM,DOMAIN);

Flux = getDiffFlux(VARIABLES, DOMAIN,BC);




%% Level set at time t = dt
[LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN);

[IBM_coeffU,IBM_coeffV,IBM_coeffP] = LSIBMcoeffs(IBM,DOMAIN,LS);

% r= 0.5;
% ang=0:0.01:2*pi;
% xcir=r*cos(ang);
% ycir=r*sin(ang);
% xc = 3.5;
% yc = 3.5;
r= IBM.diamcyl/2;
ang=0:0.01:2*pi;

xcir=r*cos(ang);
ycir=r*sin(ang);
xc = IBM.xc;
yc = IBM.yc;

for i = 1:40
        
    if i == 101
        VARIABLES.dt = 2*VARIABLES.dt;
    end
    ControlVar.resi=1;
    ControlVar.ii=0;
    ControlVar.time = ControlVar.time+VARIABLES.dt;
%% =======================================================================
%                       SOLVE FOR FLOW
% ========================================================================

    [StateVar,ControlVar] = ... 
        SolveUVP (ControlVar,Flux,DOMAIN,VARIABLES,StateVar,IBM,IBM_coeffU,...
        IBM_coeffV,IBM_coeffP,BC);    
    
    StateVar.U_old = StateVar.U;
    StateVar.V_old = StateVar.V;
    StateVar.P_old = StateVar.P;

%% ========================================================================
%                       SCALAR TRANSPORT            
%  ========================================================================
    
    if ControlVar.flow_steady == 0
        [StateVar] = SolveTransportADRE(ControlVar,DOMAIN,VARIABLES,... 
            StateVar,IBM,IBM_coeffP,BC,Flux);
    end
    
   
%% ========================================================================
%                       SAVE THE DATA            
%  ========================================================================     
    
    CurrentStateVar =  struct('U',StateVar.U,'V',StateVar.V,...
            'P',StateVar.P,'phi',StateVar.phi,'psi',LS.psi,...
            'time',ControlVar.time);
    
  
    
%     if mod(i,ControlVar.savedat) == 0
%         flowfilename = strcat('dataRDE',num2str(i),'dt.mat');
%         save(flowfilename,'-v7.3','-struct','CurrentStateVar');
% 
%         figure(4)
%         contourf(DOMAIN.xp,DOMAIN.yp,LS.psi',200,'LineStyle','none');
%         colormap jet
%         colorbar
%         axis equal
%         for iCircle = 1:size(IBM.xc,1)
%             for jCircle = 1:size(IBM.yc,2)
% 
%                 hold on
%                 plot(IBM.xc(iCircle,jCircle)+xcir,...
%                     IBM.yc(iCircle,jCircle)+ycir,'w');
% 
%             end
%         end
% %         pause(0.1)
% %         
% 
% 
% 
%         
%     end    
%% ========================================================================
%                      LEVEL SET EQUATION            
%% ========================================================================
    [LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN);
%     [LS] = LSnormals(LS,DOMAIN);
% 
    [IBM_coeffU,IBM_coeffV,IBM_coeffP] = LSIBMcoeffs(IBM,DOMAIN,LS);
% 
    

    figure(13)
    contourf(DOMAIN.xp(1:end-10),DOMAIN.yp,(LS.u(1:end-10,:))',200,...
        'LineStyle','none');
    colormap jet
    axis equal
    hold on
    plot(xc+xcir,yc+ycir,'k');
    
    figure(1)
    contourf(DOMAIN.xp(1:end-10),DOMAIN.yp,(LS.nx(1:end-10,:))',200,...
        'LineStyle','none');
    colormap jet
    axis equal
    hold on
    plot(xc+xcir,yc+ycir,'k');
    
    figure(2)
    contourf(DOMAIN.xp(1:end-10),DOMAIN.yp,(LS.psi(1:end-10,:))',200,...
        'LineStyle','none');
    colormap jet
    axis equal
    hold on
    plot(xc+xcir,yc+ycir,'k');
    pause(0.05)


end

% [StateVar.PSI,StateVar.VOR]=PlotFlowField(CurrentStateVar,IBM,DOMAIN);
% 
% 
% profile viewer