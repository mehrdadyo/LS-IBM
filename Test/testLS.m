%% test LS solver define a velocity profile and evolve the level set equation
r= IBM.diamcyl/2;
ang=0:0.01:2*pi;

xcir=r*cos(ang);
ycir=r*sin(ang);
xc = IBM.xc;
yc = IBM.yc;
for iTime = 1:250
        
    if iTime == 101
        VARIABLES.dt = 2*VARIABLES.dt;
    end
    ControlVar.resi=1;
    ControlVar.ii=0;
    ControlVar.time = ControlVar.time+VARIABLES.dt;
    time = ControlVar.time;
    time
    
  
    
%% ========================================================================
%                      LEVEL SET EQUATION            
%% ========================================================================
    [LS] = LSeqSolver(LS,StateVar,VARIABLES,DOMAIN);
% 
% 
    

%     figure(13)
%     contourf(DOMAIN.xp(1:end-10),DOMAIN.yp,(LS.u(1:end-10,:))',200,...
%         'LineStyle','none');
%     colormap jet
%     axis equal
%     hold on
%     plot(xc+xcir,yc+ycir,'k');
    if mod(iTime,10)
        figure(11)
        contourf(DOMAIN.xp(1:end-10),DOMAIN.yp,(LS.nx(1:end-10,:))',200,...
            'LineStyle','none');
        colormap jet
        axis equal
        hold on
        plot(xc+xcir,yc+ycir,'k');

        figure(12)
        contourf(DOMAIN.xp(1:end-10),DOMAIN.yp,(LS.psi(1:end-10,:)<0)',200,...
            'LineStyle','none');
        colormap jet
        axis equal
        hold on
        plot(xc+xcir,yc+ycir,'k');
    end

end


function LS = LSeqSolver(LS,StateVar,VARIABLES,DOMAIN)

    phi = StateVar.phi;
    psi_n = LS.psi;
    Da = VARIABLES.Da;
    Pe = VARIABLES.Pe;
    dt = VARIABLES.dt;
    n_iter = VARIABLES.n_iter_ReLS;
    dtau = VARIABLES.dtau;

    [u,v] = LSVelocity(phi,Da,Pe,DOMAIN,LS);

    LS.u = u;
    LS.v = v;
    equation = "LevelSetEqn";
    scheme = VARIABLES.TimeSchemeLS;
    % convective fluxes
    [psinp] = LSFindDerivative(u,v,psi_n,DOMAIN,equation);

    % psie = psinp.psie;
    % psiw = psinp.psiw;
    % psin = psinp.psin;
    % psis = psinp.psis;
    epsSign = min(min(DOMAIN.dxp))*2;

    psi_n1 = psi_n - dt*( u.*psinp.psi_x + v.*psinp.psi_y );

    if scheme == "RK1"

        psi = psi_n1; 

    elseif scheme == "RK2" 

        [psinp] = LSFindDerivative(u,v,psi_n1,DOMAIN,equation);



        psi_n2 = psi_n1 - dt*( u.*psinp.psi_x + v.*psinp.psi_y );

        psi = 0.5 * (psi_n + psi_n2);


    elseif scheme == "RK3"

        [psinp] = LSFindDerivative(u,v,psi_n1,DOMAIN,equation);

        psi_n2 = psi_n1 - dt*( u.*psinp.psi_x + v.*psinp.psi_y );


        psi_n12 = 0.75 * psi_n + 0.25 * psi_n2;

        [psinp] = LSFindDerivative(u,v,psi_n1,DOMAIN,equation);

        psi_n32 = psi_n12 - dt * ( u.*psinp.psi_x + v.*psinp.psi_y );

        psi = (psi_n + 2* psi_n32)/3;

    end

    psi_n = psi;

    %% Reinitialize the distance function
    equation = "ReinitializationEqn";

    scheme = VARIABLES.TimeSchemeRLS;

    for i=1:n_iter

        [psinp] = LSFindDerivative(u,v,psi,DOMAIN,equation);

        [G]= LSreinitilizationCoeff(psinp,psi_n);
        fl = double(G ~= 0);
        SignPsi = psi_n./sqrt(psi_n.^2+epsSign.^2);
    %     psi_m1 = psi - dtau*sign(psi_n).*(G-1);
        psi_m1 = psi - dtau*SignPsi.*(G-1).*fl;


        if scheme == "RK1"

            psi = psi_m1; 

        elseif scheme == "RK2" 

            [psinp] = LSFindDerivative(u,v,psi_m1,DOMAIN,equation);
            [G]= LSreInitilization(psinp,psi_n);
            fl = double(G ~= 0);
    %         psi_m2 = psi_m1 - dtau*sign(psi_n).*(G-1);
            psi_m2 = psi_m1 - dtau*SignPsi.*(G-1).*fl;

            psi = 0.5 * (psi + psi_m2);


        elseif scheme == "RK3"

            [psinp] = LSFindDerivative(u,v,psi_m1,DOMAIN,equation);
            [G]= LSreInitilization(psinp,psi_n);
            fl = double(G ~= 0);

    %         psi_m2 = psi_m1 - dtau*sign(psi_n).*(G-1);
            psi_m2 = psi_m1 - dtau*SignPsi.*(G-1).*fl;

            psi_m12 = 0.75 * psi + 0.25 * psi_m2;

            [psinp] = LSFindDerivative(u,v,psi_m12,DOMAIN,equation);
            [G]= LSreInitilization(psinp,psi_n);
            fl = double(G ~= 0);

    %         psi_m32 = psi_m12 - dtau*sign(psi_n).*(G-1);
            psi_m32 = psi_m12 - dtau*SignPsi.*(G-1).*fl;    

            psi = (psi + 2* psi_m32)/3;

        end


    end
    LS.psi = psi;
    [LS] = LSnormals(LS,DOMAIN);




end


function [u,v] = LSVelocity(phi,Da,Pe,DOMAIN,LS)


    %% ===================
    % This function extrapolates the velocity from the interface to a banded
    % neighborhood of the interface. This is necessary for evolving the
    % interface based on the level set equation.

    %% finding the velocity
    psi = LS.psi;

    nx = LS.nx;
    ny = LS.ny;
    
    tanTheta = ny./nx;
    theta = atan(tanTheta);
    theta = abs(theta);

    x = DOMAIN.xp;
    y = DOMAIN.yp;

    dx = min(min(DOMAIN.dxp));
    dy = min(min(DOMAIN.dyp));

    [Nx,Ny] = size(psi);
    [uI,vI] = deal(zeros(Nx,Ny));

    epsi = 10*dx;
    % i_vel = 0;
    for i=2:Nx-1
        for j=2:Ny-1
            if abs(psi(i,j))<epsi   % grids within the narrow band

%                 if psi(i,j) >= 0
% 
%                     %% outside the solid
%                     x_I = x(i) - nx(i,j)* psi(i,j) ;
%                     y_I = y(j) - ny(i,j)* psi(i,j) ;
% 
%                     x_pr = x_I + sqrt(2)*dx*nx(i,j);  % find x of x_pr
%                     y_pr = y_I + sqrt(2)*dy*ny(i,j);  % find y of x_pr
% 
%                     x_dpr = x_pr + sqrt(2)*dx*nx(i,j);
%                     y_dpr = y_pr + sqrt(2)*dy*ny(i,j);
% 
%                     d = sqrt(2)*dy;   % distance of x_pr from interface`
% 
% 
% 
%     %                 
% 
% 
%                 elseif psi(i,j)<0
%                     %% inside the solid
% 
%                     x_I = x(i) + nx(i,j)*abs( psi(i,j) );
%                     y_I = y(j) + ny(i,j)*abs( psi(i,j) );
% 
%                     x_pr = x_I + sqrt(2)*dx*nx(i,j);
%                     y_pr = y_I + sqrt(2)*dy*ny(i,j);
% 
%                     x_dpr = x_pr + sqrt(2)*dx*nx(i,j);
%                     y_dpr = y_pr + sqrt(2)*dy*ny(i,j);
% 
%                     d = sqrt(2)*dy;
% 
% 
%                 end



%                 [phi_x_pr] = biLinearInterpolation(x,y,phi,x_pr,y_pr);
%                 [phi_x_dpr] = biLinearInterpolation(x,y,phi,x_dpr,y_dpr);

                u_I = 1;% -(Da/Pe)*(2*phi_x_pr-1/2*phi_x_dpr)/(3/2+Da*d);


                uI(i,j) = u_I*nx(i,j);
                vI(i,j) = 0.1*u_I*ny(i,j);




            end
        end
    end

    % set the velocity at x equal to V and calculate the face velocities

    u = uI;
    v = vI;




end
