function [u,v] = LSVelocityExtrapolation(VARIABLES, StateVar, DOMAIN, LS)

%% ===================
% This function extrapolates the velocity from the interface to a banded
% neighborhood of the interface. This is necessary for evolving the
% interface based on the level set equation.

%% finding the velocity
psi = LS.psi;

nx = LS.nx;
ny = LS.ny;

phi = StateVar.phi;
Da = VARIABLES.Da;
Pe = VARIABLES.Pe;


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
            
            

            [phi_x_pr] = biLinearInterpolation(x,y,phi,x_pr,y_pr);
            [phi_x_dpr] = biLinearInterpolation(x,y,phi,x_dpr,y_dpr);
            
            u_I = -(Da/Pe)*(2*phi_x_pr-1/2*phi_x_dpr)/(3/2+Da*d);
                

            uI(i,j) = u_I*nx(i,j);
            vI(i,j) = u_I*ny(i,j);
            



        end
    end
end

% set the velocity at x equal to V and calculate the face velocities

u = uI;
v = vI;

%% Staggered velocity
% u = 0.5*( uI(1:end-1,:)+uI(2:end,:) );
% v = 0.5*( vI(:,1:end-1)+vI(:,2:end) );
