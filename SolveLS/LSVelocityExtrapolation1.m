function [u,v] = LSVelocityExtrapolation1(phi,Da,Pe,DOMAIN,LS)

%% ===================
% This function extrapolates the velocity from the interface to a banded
% neighborhood of the interface. This is necessary for evolving the
% interface based on the level set equation.

%% finding the velocity
psi = LS.psi;

nx = LS.nx;
ny = LS.ny;

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
                x_pr = x(i) + sqrt(2)*dx*nx(i,j);  % find x of x_pr
                y_pr = y(j) + sqrt(2)*dy*ny(i,j);  % find y of x_pr

                d = psi(i,j) + sqrt(2)*dy;   % distance of x_pr from interface`
                
%                 if psi(i,j)<=5*dx
%                     x_pr = x(i) + dx*nx(i,j);  % find x of x_pr
%                     y_pr = y(j) + dy*ny(i,j);  % find y of x_pr
% 
%                     d = dy+psi(i,j);   % distance of x_pr from interface
%                 else
%                     x_pr = x(i) - dx*nx(i,j);  % find x of x_pr
%                     y_pr = y(j) - dy*ny(i,j);  % find y of x_pr
% 
%                     d = -dy+psi(i,j);   % distance of x_pr from interface
%                 end
                    
                
                
            elseif psi(i,j)<0
                %% inside the solid
                
                x_I = x(i) + nx(i,j)*abs( psi(i,j) );
                y_I = y(j) + ny(i,j)*abs( psi(i,j) );
                
                x_pr = x_I + sqrt(2)*dx*nx(i,j);
                y_pr = y_I + sqrt(2)*dy*ny(i,j);
                

                d = sqrt(2)*dy;
                
            
            end
            
            

            [phi_x_pr] = biLinearInterpolation(x,y,phi,x_pr,y_pr);

                
            u_I = -(Da/Pe)*phi_x_pr/(1+Da*d);

            % transfer u_I to a vector v = n* u_I
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
