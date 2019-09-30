function [StateVar,ControlVar] = ... 
    SolveUVP1 (ControlVar,Flux,DOMAIN,VARIABLES,StateVar,IBM,IBM_coeffU,...
    IBM_coeffV,IBM_coeffP,BC)

    while (ControlVar.resi>ControlVar.tol || ControlVar.ii<2) && ...
            ((ControlVar.f==0)) 
%      for kk=1:10
        ControlVar.ii=ControlVar.ii+1;
        ControlVar.flc=1;
        [Flux.ConvF_U]=ConvFlux...
            (StateVar,ControlVar,DOMAIN,VARIABLES);

        [Soln.CM_u,Soln.RHS_U,Soln.a_u_w,Soln.a_u_e,Soln.a_u_s,...
            Soln.a_u_n,Soln.a_u_p,Soln.d_u,Soln.A_g_sparse,Soln.M_u] = ...
            COEFFU(StateVar.U_old,StateVar.P,StateVar.U_star,DOMAIN,...
            Flux.Diffu_U,Flux.ConvF_U,IBM_coeffU,IBM,VARIABLES,BC,...
            ControlVar.disc_scheme_vel);
        
        ControlVar.flc=-1;
        [Flux.ConvF_V]=ConvFlux...
            (StateVar,ControlVar,DOMAIN,VARIABLES);
    
        [Soln.CM_v,Soln.RHS_V,Soln.a_v_w,Soln.a_v_e,Soln.a_v_s,...
            Soln.a_v_n,Soln.a_v_p,Soln.d_v] = COEFFV(StateVar.V_old,...
            StateVar.P,StateVar.V_star,DOMAIN,Flux.Diffu_V,Flux.ConvF_V,...
            IBM_coeffV,IBM,VARIABLES,BC,...
           ControlVar.disc_scheme_vel);

        [Soln.CM_p,Soln.CM_p2] = COEFFP(Soln.d_u,Soln.d_v,BC,DOMAIN);
        %% Solve U- momentum
        % preconditioning incomplete LU decomp
        ControlVar.setup.type = 'nofill';
        ControlVar.setup.milu = 'off';
        
        [Soln.L_u,Soln.U_u] = ilu(Soln.CM_u,ControlVar.setup);
        
        [StateVar.U_star_vec,Soln.fu,Soln.relresu,Soln.iteru] = ...
            bicgstab(Soln.CM_u,Soln.RHS_U,ControlVar.tolbicg, ...
            ControlVar.maxit,Soln.L_u,Soln.U_u);                                           
        
        %% Solve V- momentum
%===================== old v momentum solver ==============================        
        % preconditioning incomplete LU decomp
        ControlVar.setup.type = 'nofill';
        ControlVar.setup.milu = 'off';
        
        [Soln.L_v,Soln.U_v] = ilu(Soln.CM_v,ControlVar.setup);
        
        [StateVar.V_star_vec,Soln.fv,Soln.relresv,Soln.iterv] = ...
            bicgstab(Soln.CM_v,Soln.RHS_V,ControlVar.tolbicg, ...
            ControlVar.maxit,Soln.L_v,Soln.U_v);

        
        [StateVar.U_star,StateVar.V_star] = FORMUV(StateVar.U_star_vec,...
            StateVar.V_star_vec,BC,ControlVar.ii,DOMAIN,...
            IBM_coeffU,IBM_coeffV);

        [Soln.RHS_P,Soln.RHS_P2] = RHSP(StateVar.U_star,StateVar.V_star,...
            DOMAIN,IBM_coeffU,IBM_coeffV);

        %% Solve Pressure- correction (SIMPLE)      
%=========================================================================
        [StateVar.P_cor_vec,Soln.f,Soln.ress,Soln.it,Soln.resvec] = ...
            agmg(Soln.CM_p2,Soln.RHS_P2,[],1e-4,ControlVar.maxit,0,...
            StateVar.P_cor_vec);

        [StateVar.PCOR] = FORMPCOR(StateVar.P_cor_vec,DOMAIN);

        [StateVar.U,StateVar.V,StateVar.P] = ...
            NEWUVP(StateVar.U_star,StateVar.V_star,StateVar.P, ...
            StateVar.PCOR,Soln.d_u,Soln.d_v,ControlVar.ii,DOMAIN,BC, ...
            IBM_coeffU,IBM_coeffV,IBM_coeffP,VARIABLES);



        StateVar.U_star_old = StateVar.U_star;
        StateVar.V_star_old = StateVar.V_star;
        StateVar.U_star = StateVar.U;
        StateVar.V_star = StateVar.V;
%         P_star=P;

        
%% Error and time and residuals

        [ControlVar.residual_vector,ControlVar.resi] = ConvergenceResiduals...
            (StateVar.U,StateVar.V,Soln.RHS_P2,ControlVar.ii, ...
            ControlVar.time,StateVar.PCOR,Soln.CM_u,Soln.CM_v, ...
            Soln.RHS_U,Soln.RHS_V, ...
            DOMAIN,ControlVar.disc_scheme_vel,ControlVar.tol);
        
        
        if (ControlVar.ii >100) || ...
                (ControlVar.ii>20 && ControlVar.resi>ControlVar.resiOld)
            disp([' The SIMPLE didnn reach the tol of = ', ...
                num2str(ControlVar.tol), ' in  ', ...
                num2str(ControlVar.ii), ' It reached ',  ...
                num2str(ControlVar.resi)]);
            break
        end
        
        ControlVar.resiOld = ControlVar.resi;


    
     
    end
    
fU = double(IBM_coeffU.flag_u == 0);
StateVar.U = StateVar.U.*fU;

fV = double(IBM_coeffV.flag_v == 0);
StateVar.V = StateVar.V.*fV;

