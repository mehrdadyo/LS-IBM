function [StateVar] = SolveTransportADRE(ControlVar,DOMAIN,VARIABLES,... 
    StateVar,IBM,IBM_coeffP,BC,Flux, LS)


        
        ControlVar.flc = 0;
        [Flux.ConvF_phi]=ConvFlux...
            (StateVar,ControlVar,DOMAIN,VARIABLES);
        

        
       [Soln.CM_phi, Soln.ap_p] =  COEFFPHIADRE(DOMAIN,Flux,IBM_coeffP,IBM, VARIABLES);

         
% ===========================  The Hayase's Quick Loop ==================== 
        
        ControlVar.err_q=1;
        ControlVar.iter_qq=0;
        err_q = 100;
        while ControlVar.err_q>ControlVar.tol_q && ControlVar.iter_qq<250
               
            ControlVar.iter_qq = ControlVar.iter_qq+1;
            [Soln.RHS_phi] = RHSPHIADRE(StateVar,BC,Flux,DOMAIN,...
                IBM_coeffP,IBM, VARIABLES, Soln);
            ControlVar.setup.type = 'nofill';
            ControlVar.setup.milu = 'off';
            [Soln.L_phi,Soln.U_phi] = ilu(Soln.CM_phi,ControlVar.setup);

            [phi_vec,ff]=bicgstab(Soln.CM_phi,Soln.RHS_phi,...
                ControlVar.tolbicg_c,ControlVar.maxit_c,...
                Soln.L_phi,Soln.U_phi);
            [StateVar.phi_n] = FORMPHI(DOMAIN,phi_vec,BC,BC.phi_a,...
                IBM_coeffP.flag_p,0, LS);
            ControlVar.err_q=max(max(abs(StateVar.phi-StateVar.phi_n)));
            StateVar.phi=StateVar.phi_n;
            
            if ControlVar.transport_steady
                if ~mod(ControlVar.iter_qq, 2)
                    disp(['The QUICK method is at iteration ',...
                        num2str(ControlVar.iter_qq), ' with residual = ',...
                        num2str(ControlVar.err_q)]);  
                end
            end
            
            
%             disp([' Iteration = ',num2str(ControlVar.iter_qq),...
%                 '  Residual Quick  = ',...
%                 num2str(ControlVar.err_q)]); 
            if err_q>ControlVar.err_q
                err_q = ControlVar.err_q;
                phi_conv = StateVar.phi_n;
            end
             
        end
        disp(['The QUICK method converged with ',num2str(ControlVar.iter_qq)...
            ' iterations ', 'with residual = ',...
            num2str(ControlVar.err_q)]);  
        
phi = phi_conv;%StateVar.phi;  
phi = (~IBM_coeffP.flag_p).* phi;

% % get indices of all cells which are in fluids and have negative phi and
% % larger than one phi
% negIndx = phi<0 & IBM_coeffP.flag_p==0;
% posIndx = phi>1 & IBM_coeffP.flag_p==0;
% % get number of cells that are inside fluid and have positive phi and less
% % than phi
% insideVol = sum(sum((phi>0 & IBM_coeffP.flag_p==0 & phi<1)));
% % get the sum of all the negative phis and larget than one phis
% totNegflux = sum(sum(phi(negIndx))) + sum(sum(phi(posIndx)));
% % % add the negative mass to the fluid cells
% phi(phi>0 & IBM_coeffP.flag_p==0 & phi<=1) = ...
%     phi(phi>0 & IBM_coeffP.flag_p==0 & phi<=1) + ...
%     totNegflux/insideVol;
% % reset the negative cells to zero and positives to one
% phi(negIndx) = 0;
% phi(posIndx) = 1;
phi(phi<0) = 0;
StateVar.phi = phi;  
% disp([' The difference in phi = ', ...
%     num2str(max(max(abs(StateVar.phi_old - StateVar.phi))))]);
StateVar.dphi_max = max(max(abs(StateVar.phi_old - StateVar.phi)));
StateVar.dphi_L2 = sum(sum((StateVar.phi_old - StateVar.phi).^2))/...
    nnz(StateVar.phi);

StateVar.phi_old = StateVar.phi;

%     if mod(i,savedat)==0
%         name_flowfile=strcat('flowexpongrid',num2str(i),'.mat');
%         save(name_flowfile,'-v7.3','U','V','P')
%     end
% negIndx = phi<0 & LS.psi>0;
