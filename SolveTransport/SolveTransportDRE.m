function [StateVar] = SolveTransportDRE(ControlVar,DOMAIN,VARIABLES,... 
    StateVar,IBM,IBM_coeffP,BC,Flux,LS)


        
        
frontCells = length(find(abs(LS.psi)<10*min(min(DOMAIN.dxp))));
noGhostCells = length(IBM_coeffP.I_g_p);
noSolidCells = length(IBM_coeffP.I_solid_p);
        
    [Soln.CM_phi] = COEFFPHIDRE(DOMAIN,Flux,IBM_coeffP);

         
% ===========================  The Hayase's Quick Loop ==================== 
        

        
    [Soln.RHS_phi] = RHSPHIDRE(StateVar,BC,Flux,DOMAIN,IBM_coeffP,IBM);
    ControlVar.setup.type = 'nofill';
    ControlVar.setup.milu = 'off';
    [Soln.L_phi,Soln.U_phi] = ilu(Soln.CM_phi,ControlVar.setup);
    [phi_vec,ff]=bicgstab(Soln.CM_phi,Soln.RHS_phi,...
        ControlVar.tolbicg_c,ControlVar.maxit_c,...
        Soln.L_phi,Soln.U_phi);
    [StateVar.phi_n] = FORMPHI(DOMAIN,phi_vec,BC,BC.phi_a,...
        IBM_coeffP.flag_p,0);
    StateVar.phi=StateVar.phi_n;
    StateVar.phi_old = StateVar.phi;
    
    if frontCells  == 0
            msg = ['   No of front cells =  ',num2str(frontCells), ...
                '  The LS shrank !!!!'  ];
            error(msg); 

            
    end
        
    
    disp([' Time = ',num2str(ControlVar.time), ...
        '  u =  ',num2str(max(max(abs(LS.u)))), ...
        '   No of front cells =  ',num2str(frontCells), ...
        '   No of ghoset cells =  ',num2str(noGhostCells), ...
        '   No of solid cells =  ',num2str(noSolidCells)]);  

%     if mod(i,savedat)==0
%         name_flowfile=strcat('flowexpongrid',num2str(i),'.mat');
%         save(name_flowfile,'-v7.3','U','V','P')
%     end
