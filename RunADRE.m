clc
clear all
close all
profile on



[StateVar,VARIABLES,DOMAIN,BC,IBM,LS]=SetUpVariables;

%% ==============================
% load 'initial at any time '

StateVar.U_old = StateVar.U;
StateVar.V_old = StateVar.V;
StateVar.P_old = StateVar.P;

% StateVar.U_star = StateVar.U;
% StateVar.V_star = StateVar.V;
% StateVar.U_star_old = StateVar.U_star;
% StateVar.V_star_old = StateVar.V_star;
        
%% ==============================
StateVar.phi = StateVar.phi_old;

%% object location and IBM coefficients
% [IBM_coeffU,IBM_coeffV,IBM_coeffP]=IBMcoeffs(IBM,DOMAIN);

%% U-momentum diffusive flux coeff
Flux.fl=1;
[Flux.Diffu_U]=DiffFlux(VARIABLES,Flux.fl,DOMAIN,BC);

%% V-momentum diffusive flux coeff
Flux.fl=-1;
[Flux.Diffu_V]=DiffFlux(VARIABLES,Flux.fl,DOMAIN);

%% Scalar Transport diffusive flux coeff
Flux.fl=0;
[Flux.Diffu_P]=DiffFlux(VARIABLES,Flux.fl,DOMAIN);



%% Control VAriables
ControlVar.PISO = 1;
ControlVar.time=0;
ControlVar.timedt = 50000;
ControlVar.savedat=50;
ControlVar.rat=floor(0.01/VARIABLES.dt);

ControlVar.tol=1e-3;
ControlVar.f=0;
ControlVar.tolbicg=10^(-5); %% max tol for bicon
ControlVar.maxit=2000;   %% max iter for bicon 

ControlVar.flow_steady= 1;
ControlVar.disc_scheme_vel = 2;
StateVar.P_cor_vec = zeros(1,(DOMAIN.imax-1)*(DOMAIN.jmax-1))';

ControlVar.flow_steady = 0;

ControlVar.tol_q=1e-5;
ControlVar.tolbicg_c=5e-5;
ControlVar.maxit_c=2000;

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