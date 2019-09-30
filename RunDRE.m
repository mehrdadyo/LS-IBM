clc
clear all
close all
% profile on



[StateVar,VARIABLES,DOMAIN,BC,IBM,LS]=SetUpVariables;


        
%% ==============================
StateVar.phi = StateVar.phi_old;


%% Scalar Transport diffusive flux coeff
Flux.fl=0;
[Flux.Diffu_P]=DiffFlux(VARIABLES,Flux.fl,DOMAIN);



%% Control VAriables

ControlVar.time=0;
ControlVar.timedt = 10000;
ControlVar.savedat=25;
ControlVar.rat=floor(0.001/VARIABLES.dt);

ControlVar.tol=1e-3;
ControlVar.f=0;
ControlVar.tolbicg=10^(-5); %% max tol for bicon
ControlVar.maxit=2000;   %% max iter for bicon 

ControlVar.flow_steady= 1;
ControlVar.disc_scheme_vel = 2;
StateVar.P_cor_vec = zeros(1,(DOMAIN.imax-1)*(DOMAIN.jmax-1))';


ControlVar.tol_q=1e-5;
ControlVar.tolbicg_c=5e-5;
ControlVar.maxit_c=2000;

%% Level set at time t = dt
[LS] = LSnormals(LS,DOMAIN);
[LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN);
% [LS] = LSnormals(LS,DOMAIN);

[~,~,IBM_coeffP] = LSIBMcoeffs(IBM,DOMAIN,LS);

r= 0.5;
ang=0:0.01:2*pi;
xcir=r*cos(ang);
ycir=r*sin(ang);
xc = 5;
yc = 5;


% profile viewer
for i = 1:ControlVar.timedt
        
%         if i == 21
%             VARIABLES.dt = 10*VARIABLES.dt;
%         end
    
    ControlVar.resi=1;
    ControlVar.ii=0;
    ControlVar.time = ControlVar.time+VARIABLES.dt;
%% =======================================================================
%                       SOLVE FOR FLOW
% ========================================================================
    % There is no flow in DRE.

%% ========================================================================
%                       SCALAR TRANSPORT            
%  ========================================================================
    
    [StateVar] = SolveTransportDRE(ControlVar,DOMAIN,VARIABLES,... 
        StateVar,IBM,IBM_coeffP,BC,Flux,LS);
    
   
%% ========================================================================
%                       SAVE THE DATA            
%  ========================================================================     
    
    CurrentStateVar =  struct('phi',StateVar.phi,'psi',LS.psi,...
            'time',ControlVar.time);
    
        if mod(i,100) == 0
            
        flowfilename = strcat('dataRDE',num2str(i),'dt.mat');
        save(flowfilename,'-v7.3','-struct','CurrentStateVar');
        
        end   
        
        
    if mod(i,25) == 0
        flowfilename = strcat('dataRDE',num2str(i),'dt.mat');
        save(flowfilename,'-v7.3','-struct','CurrentStateVar');
        figure(1)
        contourf(DOMAIN.xp(150:600),DOMAIN.yp(150:600),...
            LS.psi(150:600,150:600)',200,'LineStyle','none');
        colormap jet
        colorbar
        axis equal
        hold on
        plot(xc+xcir,yc+ycir,'w');
        pause(1)
        
        figure(2)
        contourf(DOMAIN.xp(150:600),DOMAIN.yp(150:600),...
            StateVar.phi(150:600,150:600)',200,'LineStyle','none');
        colormap jet
        colorbar
        axis equal
        hold on
        plot(xc+xcir,yc+ycir,'w');
        pause(1)


        
    end    
    
    
%% ========================================================================
%                      LEVEL SET EQUATION            
%% ========================================================================
%     [LS] = LSnormals(LS,DOMAIN);
    [LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN);
%      [LS] = LSnormals(LS,DOMAIN);

    [~,~,IBM_coeffP] = LSIBMcoeffs(IBM,DOMAIN,LS);

    



end

% [StateVar.PSI,StateVar.VOR]=PlotFlowField(CurrentStateVar,IBM,DOMAIN);


profile viewer