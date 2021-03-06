function psi = LSreinitialization(VARIABLES, DOMAIN, LS, psi_o)
    
    equation = "ReinitializationEqn";
    epsSign = min(min(DOMAIN.dxp));

    scheme = VARIABLES.TimeSchemeRLS;
    n_iter = VARIABLES.n_iter_ReLS;
    dtau = VARIABLES.dtau;
    u = LS.u;
    v = LS.v;
    psi = LS.psi;
    psi_n = psi;
    h = VARIABLES.LSgamma;
    for i=1:n_iter

        [psinp] = LSFindDerivative(u,v,psi,DOMAIN,equation, h, psi_o);

        [G]= LSreinitilizationCoeff(psinp,psi_n);
        fl = double(G ~= 0);
        SignPsi = psi_n./sqrt(psi_n.^2+epsSign.^2);
    %     psi_m1 = psi - dtau*sign(psi_n).*(G-1);
        psi_m1 = psi - dtau*SignPsi.*(G-1).*fl;


        if scheme == "RK1"

            psi = psi_m1; 

        elseif scheme == "RK2" 

            [psinp] = LSFindDerivative(u,v,psi_m1,DOMAIN,equation, h, psi_o);
            [G]= LSreinitilizationCoeff(psinp,psi_n);
            fl = double(G ~= 0);
    %         psi_m2 = psi_m1 - dtau*sign(psi_n).*(G-1);
            psi_m2 = psi_m1 - dtau*SignPsi.*(G-1).*fl;

            psi = 0.5 * (psi + psi_m2);


        elseif scheme == "RK3"

            [psinp] = LSFindDerivative(u,v,psi_m1,DOMAIN,equation, h, psi_o);
            [G]= LSreinitilizationCoeff(psinp,psi_n);
            fl = double(G ~= 0);

    %         psi_m2 = psi_m1 - dtau*sign(psi_n).*(G-1);
            psi_m2 = psi_m1 - dtau*SignPsi.*(G-1).*fl;

            psi_m12 = 0.75 * psi + 0.25 * psi_m2;

            [psinp] = LSFindDerivative(u,v,psi_m12,DOMAIN,equation, h, psi_o);
            [G]= LSreinitilizationCoeff(psinp,psi_n);
            fl = double(G ~= 0);

    %         psi_m32 = psi_m12 - dtau*sign(psi_n).*(G-1);
            psi_m32 = psi_m12 - dtau*SignPsi.*(G-1).*fl;    

            psi = (psi + 2* psi_m32)/3;

        end


    end
    psi(psi > h) = h;
    psi(psi < -1*h) = -h;
end

