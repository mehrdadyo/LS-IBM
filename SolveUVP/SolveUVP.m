function [StateVar,ControlVar,Soln] = ... 
    SolveUVP (ControlVar,Flux,DOMAIN,VARIABLES,StateVar,IBM,IBM_coeffU,...
    IBM_coeffV,IBM_coeffP,BC)

    U_star = StateVar.U;
    V_star = StateVar.V;
    U_star_old = U_star;
    V_star_old = V_star;
    
    ControlVar.resi=1;
    ControlVar.ii=0;
    

    while (ControlVar.resi>ControlVar.tol || ControlVar.ii<2) && ...
            ((ControlVar.f==0)) 
%      for kk=1:10
        ControlVar.ii=ControlVar.ii+1;
        ControlVar.flc=1;
        [Flux.ConvF_U]=ConvFlux...
            (StateVar,ControlVar,DOMAIN,VARIABLES,BC);

        [Soln.CM_u,Soln.RHS_U,Soln.a_u_w,Soln.a_u_e,Soln.a_u_s,...
            Soln.a_u_n,Soln.a_u_p,Soln.d_u,Soln.A_g_sparse,Soln.M_u,...
            Soln.d_u_sec] = ...
            COEFFU(StateVar.U_old,StateVar.P,U_star,DOMAIN,...
            Flux.Diffu_U,Flux.ConvF_U,IBM_coeffU,IBM,VARIABLES,BC,...
            ControlVar.disc_scheme_vel);
        
        ControlVar.flc=-1;
        [Flux.ConvF_V]=ConvFlux...
            (StateVar,ControlVar,DOMAIN,VARIABLES);
    
        [Soln.CM_v,Soln.RHS_V,Soln.a_v_w,Soln.a_v_e,Soln.a_v_s,...
            Soln.a_v_n,Soln.a_v_p,Soln.d_v,Soln.d_v_sec] = ...
            COEFFV(StateVar.V_old,...
            StateVar.P,V_star,DOMAIN,Flux.Diffu_V,Flux.ConvF_V,...
            IBM_coeffV,IBM,VARIABLES,BC,...
           ControlVar.disc_scheme_vel);

        [Soln.CM_p,Soln.CM_p2, ap, aw, ae, as, an, Sp] = COEFFP(Soln.d_u,Soln.d_v,BC,DOMAIN);
        %% Solve U- momentum
        % preconditioning incomplete LU decomp
        ControlVar.setup.type = 'nofill';
        ControlVar.setup.milu = 'off';
        
        [Soln.L_u,Soln.U_u] = ilu(Soln.CM_u,ControlVar.setup);
        
        [U_star_vec,Soln.fu,Soln.relresu,Soln.iteru] = ...
            bicgstab(Soln.CM_u,Soln.RHS_U,ControlVar.tolbicg, ...
            ControlVar.maxit,Soln.L_u,Soln.U_u);                                           
        
        %% Solve V- momentum
%===================== old v momentum solver ==============================        
        % preconditioning incomplete LU decomp
        ControlVar.setup.type = 'nofill';
        ControlVar.setup.milu = 'off';
        
        [Soln.L_v,Soln.U_v] = ilu(Soln.CM_v,ControlVar.setup);
        
        [V_star_vec,Soln.fv,Soln.relresv,Soln.iterv] = ...
            bicgstab(Soln.CM_v,Soln.RHS_V,ControlVar.tolbicg, ...
            ControlVar.maxit,Soln.L_v,Soln.U_v);

        
        [U_star,V_star] = FORMUV(U_star_vec,...
            V_star_vec,BC,ControlVar.ii,DOMAIN,...
            IBM_coeffU,IBM_coeffV);

        [Soln.RHS_P,Soln.RHS_P2] = RHSP(U_star,V_star,...
            DOMAIN,IBM_coeffU,IBM_coeffV);

        %% Solve Pressure- correction (SIMPLE)      
%=========================================================================
        [StateVar.P_cor_vec,Soln.f,Soln.ress,Soln.it,Soln.resvec] = ...
            agmg(Soln.CM_p2,Soln.RHS_P2,[],1e-4,ControlVar.maxit,0,...
            StateVar.P_cor_vec);
%         StateVar.P_cor_vec = Soln.CM_p2\Soln.RHS_P2;
%         [Soln.L_p,Soln.U_p] = ilu(Soln.CM_p2,ControlVar.setup);      
%         [StateVar.P_cor_vec,Soln.f,Soln.ress,Soln.it] = ...
%             pcg(Soln.CM_p2,Soln.RHS_P2,ControlVar.tolbicg, ...
%             ControlVar.maxit,Soln.L_p,Soln.U_p);

        [Soln.PCOR] = FORMPCOR(StateVar.P_cor_vec,DOMAIN);

        [StateVar.U,StateVar.V,StateVar.P] = ...
            NEWUVP(U_star,V_star,StateVar.P, ...
            Soln,ControlVar,DOMAIN,BC, ...
            IBM_coeffU,IBM_coeffV,IBM_coeffP,VARIABLES,0);



        U_star_old = U_star;
        V_star_old = V_star;
        U_star = StateVar.U;
        V_star = StateVar.V;
%         P_star=P;

        if ControlVar.PISO == 1
            
             Soln.dU_star = StateVar.U - U_star_old;
             Soln.dV_star = StateVar.V - V_star_old;

            [Soln.RHS_P_pr,Soln.RHS_P2_pr] = RHSP_PISO(Soln,U_star,...
                V_star,DOMAIN,IBM_coeffU,IBM_coeffV);
            [StateVar.P_cor_vec,Soln.f,Soln.ress,Soln.it,Soln.resvec] = ...
                agmg(Soln.CM_p2,Soln.RHS_P2_pr,[],1e-4,ControlVar.maxit,0,...
                StateVar.P_cor_vec);
            [Soln.PCOR] = FORMPCOR(StateVar.P_cor_vec,DOMAIN);

            [StateVar.U,StateVar.V,StateVar.P] = ...
                NEWUVP(U_star,V_star,StateVar.P, ...
                Soln,ControlVar,DOMAIN,BC, ...
                IBM_coeffU,IBM_coeffV,IBM_coeffP,VARIABLES,ControlVar.PISO);


            [ControlVar.residual_vector,ControlVar.resi, ...
                ControlVar.messageFlow] = ConvergenceResiduals...
                (StateVar.U,StateVar.V,Soln.RHS_P2,ControlVar.ii, ...
                ControlVar.time,Soln.PCOR,Soln.CM_u,Soln.CM_v, ...
                Soln.RHS_U,Soln.RHS_V, ...
                DOMAIN,ControlVar.disc_scheme_vel,ControlVar.tol,BC,...
                ControlVar.PISO);

            break

            
        end

        
%% Error and time and residuals

        [ControlVar.residual_vector,ControlVar.resi, ...
            ControlVar.messageFlow] = ConvergenceResiduals...
            (StateVar.U,StateVar.V,Soln.RHS_P2,ControlVar.ii, ...
            ControlVar.time,Soln.PCOR,Soln.CM_u,Soln.CM_v, ...
            Soln.RHS_U,Soln.RHS_V, ...
            DOMAIN,ControlVar.disc_scheme_vel,ControlVar.tol,...
            BC, ControlVar.PISO);
        
        if (ControlVar.ii >100)
            disp([' The SIMPLE didnn reach the tol of = ', ...
                num2str(ControlVar.tol), ' in  ', ...
                num2str(ControlVar.ii), ' It reached ',  ...
                num2str(ControlVar.resi)]);
            break
        end
        
        
        
%         if (ControlVar.ii  >= 100) || (ControlVar.resi > 6)
%             
%             disp([' The SIMPLE did not reach the tol of = ', ...
%                 num2str(ControlVar.tol), ' in  ', ...
%                 num2str(ControlVar.ii), ' It reached ',  ...
%                 num2str(ControlVar.resi), ' reducing the under relaxation']);
% %             break
%             if ~mod(ControlVar.ii,100)
%                 VARIABLES.alpha_u = VARIABLES.alpha_u/2;
%                 VARIABLES.alpha_v = VARIABLES.alpha_v/2;
%             end
%             if ControlVar.residual_vector(1) > 6 || ControlVar.residual_vector(3) > 6
%                 VARIABLES.alpha_u = VARIABLES.alpha_u/2;
%             end
%             if ControlVar.residual_vector(2) > 6 || ControlVar.residual_vector(3) > 6
%                 VARIABLES.alpha_v = VARIABLES.alpha_v/2;
%             end
% %             VARIABLES.alpha_v = VARIABLES.alpha_v/2;
%         end

%         if (ControlVar.ii >150) || (ControlVar.resi > 20)
%             disp([' The SIMPLE did not reach the tol of = ', ...
%                 num2str(ControlVar.tol), ' in  ', ...
%                 num2str(ControlVar.ii), ' It reached ',  ...
%                 num2str(ControlVar.resi)]);
%             ControlVar.reduceTime = 1;
%             return;
%         end
        
        ControlVar.resiOld = ControlVar.resi;


    
     
    end
    
fU = double(IBM_coeffU.flag_u == 0);
StateVar.U = StateVar.U.*fU;

fV = double(IBM_coeffV.flag_v == 0);
StateVar.V = StateVar.V.*fV;

StateVar.dU_max = max(max(abs(StateVar.U_old - StateVar.U)));
StateVar.dU_L2 = sum(sum((StateVar.U_old - StateVar.U).^2))/...
    nnz(fU);

StateVar.dV_max = max(max(abs(StateVar.V_old - StateVar.V)));
StateVar.dV_L2 = sum(sum((StateVar.V_old - StateVar.V).^2))/...
    nnz(fV);



StateVar.U_old = StateVar.U;
StateVar.V_old = StateVar.V;
StateVar.P_old = StateVar.P;

