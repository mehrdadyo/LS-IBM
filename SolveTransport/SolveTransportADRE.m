function [StateVar] = SolveTransportADRE(ControlVar,DOMAIN,VARIABLES,... 
    StateVar,IBM,IBM_coeffP,BC,Flux)


        
        ControlVar.flc = 0;
        [Flux.ConvF_phi]=ConvFlux...
            (StateVar,ControlVar,DOMAIN,VARIABLES);
        

        
       [Soln.CM_phi] =  COEFFPHIADRE(DOMAIN,Flux,IBM_coeffP,IBM);

         
% ===========================  The Hayase's Quick Loop ==================== 
        
        ControlVar.err_q=1;
        ControlVar.iter_qq=0;

        while ControlVar.err_q>ControlVar.tol_q && ControlVar.iter_qq<50
               
            ControlVar.iter_qq = ControlVar.iter_qq+1;
            [Soln.RHS_phi] = RHSPHIADRE(StateVar,BC,Flux,DOMAIN,IBM_coeffP,IBM);
            ControlVar.setup.type = 'nofill';
            ControlVar.setup.milu = 'off';
            [Soln.L_phi,Soln.U_phi] = ilu(Soln.CM_phi,ControlVar.setup);

            [phi_vec,ff]=bicgstab(Soln.CM_phi,Soln.RHS_phi,...
                ControlVar.tolbicg_c,ControlVar.maxit_c,...
                Soln.L_phi,Soln.U_phi);
            [StateVar.phi_n] = FORMPHI(DOMAIN,phi_vec,BC,BC.phi_a,...
                IBM_coeffP.flag_p,0);
            ControlVar.err_q=max(max(abs(StateVar.phi-StateVar.phi_n)));
            StateVar.phi=StateVar.phi_n;
        end
        disp([' Iteration = ',num2str(ControlVar.iter_qq),' Time = ',...
            num2str(ControlVar.time),'  Residual Quick  = ',...
            num2str(ControlVar.err_q)]);  
        
        StateVar.phi_old = StateVar.phi;
%     if mod(i,savedat)==0
%         name_flowfile=strcat('flowexpongrid',num2str(i),'.mat');
%         save(name_flowfile,'-v7.3','U','V','P')
%     end
