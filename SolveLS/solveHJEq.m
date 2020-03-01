function psi = solveHJEq(LS, VARIABLES, DOMAIN) 

    dt = VARIABLES.dtLS;
   
    u = LS.u;
    v = LS.v;
    psi_n = LS.psi;
 
    equation = "LevelSetEqn";
    scheme = VARIABLES.TimeSchemeLS;
    h = VARIABLES.LSband;

    %% convective fluxes
    [psinp] = LSFindDerivative(u,v,psi_n,DOMAIN,equation, h);

    %% Time integrations
    psi_n1 = psi_n - dt*( u.*psinp.psi_x + v.*psinp.psi_y );

    if scheme == "RK1"

        psi = psi_n1; 

    elseif scheme == "RK2" 

        [psinp] = LSFindDerivative(u,v,psi_n1,DOMAIN,equation, h);



        psi_n2 = psi_n1 - dt*( u.*psinp.psi_x + v.*psinp.psi_y );

        psi = 0.5 * (psi_n + psi_n2);


    elseif scheme == "RK3"

        [psinp] = LSFindDerivative(u,v,psi_n1,DOMAIN,equation, h);

        psi_n2 = psi_n1 - dt*( u.*psinp.psi_x + v.*psinp.psi_y );


        psi_n12 = 0.75 * psi_n + 0.25 * psi_n2;

        [psinp] = LSFindDerivative(u,v,psi_n1,DOMAIN,equation, h);

        psi_n32 = psi_n12 - dt * ( u.*psinp.psi_x + v.*psinp.psi_y );

        psi = (psi_n + 2* psi_n32)/3;

    end

end