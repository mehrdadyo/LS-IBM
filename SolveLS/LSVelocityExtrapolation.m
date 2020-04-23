function [LS] = LSVelocityExtrapolation(VARIABLES, StateVar, DOMAIN, LS)

%% ===================
% This function extrapolates the velocity from the interface to a banded
% neighborhood of the interface. This is necessary for evolving the
% interface based on the level set equation.

%% finding the velocity
dissolution = VARIABLES.dissolution;
Pe = VARIABLES.Pe;

if isfield(LS, 'LS1')
    psi = LS.LS_s.psi;
    nx = LS.LS_s.nx;
    ny = LS.LS_s.ny;

    
    psi1 = LS.LS1.psi;
    psi2 = LS.LS2.psi;
    
    Da_save = VARIABLES.Da;
    
    [Nx,Ny] = size(psi);
    [LS.LS_s.u, LS.LS_s.v] = deal(zeros(Nx,Ny));
    [LS.LS1.u, LS.LS1.v] = deal(zeros(Nx,Ny));
    [LS.LS2.u, LS.LS2.v] = deal(zeros(Nx,Ny));
else
    psi = LS.psi;
    nx = LS.nx;
    ny = LS.ny;
    
    direc = -1;
    if dissolution
        direc = 1;
    end
    
    [Nx,Ny] = size(psi);
    [LS.u, LS.v] = deal(zeros(Nx,Ny));



end

gamma = VARIABLES.LSgamma;
beta = VARIABLES.LSbeta;



phi = StateVar.phi;


x = DOMAIN.xp;
y = DOMAIN.yp;

dx = min(min(DOMAIN.dxp));
dy = min(min(DOMAIN.dyp));




% epsi = 10*dx;
% i_vel = 0;

for i=10:Nx-9
    for j=2:Ny-1
        if abs(psi(i,j))< gamma   % grids within the narrow band
            
            if psi(i,j) >= 0
                            
                %% outside the solid
                x_I = x(i) - nx(i,j)* psi(i,j) ;
                y_I = y(j) - ny(i,j)* psi(i,j) ;
                
                x_pr = x_I + sqrt(2)*dx*nx(i,j);  % find x of x_pr
                y_pr = y_I + sqrt(2)*dy*ny(i,j);  % find y of x_pr
                
                x_dpr = x_pr + sqrt(2)*dx*nx(i,j);
                y_dpr = y_pr + sqrt(2)*dy*ny(i,j);

                d = sqrt(2)*dy;   % distance of x_pr from interface`
                
                
                
%                 
                
                
            elseif psi(i,j)<0
                %% inside the solid
                
                x_I = x(i) + nx(i,j)*abs( psi(i,j) );
                y_I = y(j) + ny(i,j)*abs( psi(i,j) );
                
                x_pr = x_I + sqrt(2)*dx*nx(i,j);
                y_pr = y_I + sqrt(2)*dy*ny(i,j);
                
                x_dpr = x_pr + sqrt(2)*dx*nx(i,j);
                y_dpr = y_pr + sqrt(2)*dy*ny(i,j);
                
                d = sqrt(2)*dy;
                
            
            end
            
            if isfield(LS, 'LS1')

            
                I = find(x<x_I);
                J = find(y<y_I);

                I=I(end);
                J=J(end);

                I1 = I;
                I2 = I1+1;
                I3 = I1;
                I4 = I1+1;

                J1 = J;
                J2 = J1;   
                J3 = J1+1;
                J4 = J1+1;

                isfluid = psi(I1,J1) || psi(I2,J2) || psi(I3,J3) || psi(I4,J4);
                issolid1 = psi1(I1,J1)<0 || psi1(I2,J2)<0 ...
                    || psi1(I3,J3)<0 || psi1(I4,J4)<0;
                issolid2 = psi2(I1,J1)<0 || psi2(I2,J2)<0 ...
                    || psi2(I3,J3)<0 || psi2(I4,J4)<0;
                
                direc = -1;

                
                if isfluid
                    if issolid1 
                        Da = Da_save.Da1;
                        
                        if dissolution.dissolution1
                            direc = 1;
                        end

                        
                    elseif issolid2
                        Da = Da_save.Da2;
                        if dissolution.dissolution2
                            direc = 1;
                        end
                        
                    end
                end
                
            end
                        


            [phi_x_pr] = biLinearInterpolation(x,y,phi,x_pr,y_pr);
            if i<Nx-9 && i>10
                [phi_x_dpr] = biLinearInterpolation(x,y,phi,x_dpr,y_dpr);
                u_I = direc *(Da/Pe)*(2*phi_x_pr-1/2*phi_x_dpr)/(3/2+Da*d);
            else
                u_I = direc*(Da/Pe)*(phi_x_pr)/(1+Da*d);
            end
            coeff = ((abs(psi(i,j)) - gamma)^2 ...
            *(2*abs(psi(i,j)) + gamma - 3 * beta))/(gamma - beta)^3;    
            c = (abs(psi(i,j))<= beta) + ...
                (abs(psi(i,j))<= gamma) * (abs(psi(i,j))> beta) * coeff +...
                (1 - (abs(psi(i,j))<= gamma));
%             c = 1;    
            u_I = c * u_I;
            
            if isfield(LS, 'LS1')
                LS.LS_s.u(i,j) = u_I*nx(i,j);
                LS.LS_s.v(i,j) = u_I*ny(i,j);
                
                if issolid1 
                    LS.LS1.u(i,j) = u_I*nx(i,j);
                    LS.LS1.v(i,j) = u_I*ny(i,j);

                elseif issolid2
                    LS.LS2.u(i,j) = u_I*nx(i,j);
                    LS.LS2.v(i,j) = u_I*ny(i,j);

                end
                
            else
                LS.u(i,j) = u_I*nx(i,j);
                LS.v(i,j) = u_I*ny(i,j);
            end
            
            uI(i,j) = u_I*nx(i,j);
            vI(i,j) = u_I*ny(i,j);
            



        end
    end
end

% set the velocity at x equal to V and calculate the face velocities

% LS.u = uI;
% LS.v = vI;
maxTravelU = max(max(abs(uI)))*VARIABLES.dt;
maxTravelV = max(max(abs(vI)))*VARIABLES.dt;
disp(['maxTravelU = ', num2str(maxTravelU),'     ','maxTravelV = ', ...
    num2str(maxTravelV)])
%% Staggered velocity
% u = 0.5*( uI(1:end-1,:)+uI(2:end,:) );
% v = 0.5*( vI(:,1:end-1)+vI(:,2:end) );
