%% =======================
% This function finds the first order derivatives in 5 point stencil
% The inputs are:
% - grid sizes in x and y direction
% - the level set matrix
% - the equation type: level set equation or reinitialization equation


%%


function [dir] = LSdirDerivates(DOMAIN,psi,u,v,equation)


dx = min(min(DOMAIN.dxp));
dy = min(min(DOMAIN.dyp));

% psi = LS.psi;
% u = LS.u;
% v = LS.psi;

%% x-dir 
[nx,ny] = size(psi);
% [v1,v2,v3,v4,v5] = deal(zeros(nx,ny));
%% direction


[vx1,vx2,vx3,vx4,vx5,vy1,vy2,vy3,vy4,vy5] = deal(zeros(nx,ny));
[vxn1,vxn2,vxn3,vxn4,vxn5,vyn1,vyn2,vyn3,vyn4,vyn5] = deal(zeros(nx,ny));
[vxp1,vxp2,vxp3,vxp4,vxp5,vyp1,vyp2,vyp3,vyp4,vyp5] = deal(zeros(nx,ny));

% tic;
if equation == "LevelSetEqn"
    for i=4:nx-3
        for j=2:ny-1
            if abs(psi(i,j))<10*dx   %% 4 grid points near the interface
                if u(i,j)>=0
                % x-direction
                    vx1(i,j) = 1/dx*(psi(i-2,j)-psi(i-3,j)); % 
                    vx2(i,j) = 1/dx*(psi(i-1,j)-psi(i-2,j));
                    vx3(i,j) = 1/dx*(psi(i,j)-psi(i-1,j));
                    vx4(i,j) = 1/dx*(psi(i+1,j)-psi(i,j));
                    vx5(i,j) = 1/dx*(psi(i+2,j)-psi(i+1,j));
                elseif u(i,j)<0
                    vx1(i,j) = 1/dx*(psi(i+3,j)-psi(i+2,j)); % 
                    vx2(i,j) = 1/dx*(psi(i+2,j)-psi(i+1,j));
                    vx3(i,j) = 1/dx*(psi(i+1,j)-psi(i,j));
                    vx4(i,j) = 1/dx*(psi(i,j)-psi(i-1,j));
                    vx5(i,j) = 1/dx*(psi(i-1,j)-psi(i-2,j));
                end
               
            end
        end
    end
    for i=2:nx-1
        for j=4:ny-3
            if abs(psi(i,j))<10*dx   %% 4 grid points near the interface
                
                % y-direction
                if v(i,j)>=0
                    vy1(i,j) = 1/dy*(psi(i,j-2)-psi(i,j-3)); % 
                    vy2(i,j) = 1/dy*(psi(i,j-1)-psi(i,j-2));
                    vy3(i,j) = 1/dy*(psi(i,j)-psi(i,j-1));
                    vy4(i,j) = 1/dy*(psi(i,j+1)-psi(i,j));
                    vy5(i,j) = 1/dy*(psi(i,j+2)-psi(i,j+1));
                elseif v(i,j)<0
                    vy1(i,j) = 1/dy*(psi(i,j+3)-psi(i,j+2)); % 
                    vy2(i,j) = 1/dy*(psi(i,j+2)-psi(i,j+1));
                    vy3(i,j) = 1/dy*(psi(i,j+1)-psi(i,j));
                    vy4(i,j) = 1/dy*(psi(i,j)-psi(i,j-1));
                    vy5(i,j) = 1/dy*(psi(i,j-1)-psi(i,j-2));

                end
            end
        end
    end
    
    dir.vx1 = vx1;
    dir.vx2 = vx2;
    dir.vx3 = vx3;
    dir.vx4 = vx4;
    dir.vx5 = vx5;

    dir.vy1 = vy1;
    dir.vy2 = vy2;
    dir.vy3 = vy3;
    dir.vy4 = vy4;
    dir.vy5 = vy5;
elseif equation == "ReinitializationEqn"
    for i=4:nx-3
        for j=2:ny-1
            if abs(psi(i,j))<10*dx   %% 4 grid points near the interface
                
                % x-direction

                
                vxn1(i,j) = 1/dx*(psi(i-2,j)-psi(i-3,j)); % 
                vxn2(i,j) = 1/dx*(psi(i-1,j)-psi(i-2,j));
                vxn3(i,j) = 1/dx*(psi(i,j)-psi(i-1,j));
                vxn4(i,j) = 1/dx*(psi(i+1,j)-psi(i,j));
                vxn5(i,j) = 1/dx*(psi(i+2,j)-psi(i+1,j));
                
                vxp1(i,j) = 1/dx*(psi(i+3,j)-psi(i+2,j)); % 
                vxp2(i,j) = 1/dx*(psi(i+2,j)-psi(i+1,j));
                vxp3(i,j) = 1/dx*(psi(i+1,j)-psi(i,j));
                vxp4(i,j) = 1/dx*(psi(i,j)-psi(i-1,j));
                vxp5(i,j) = 1/dx*(psi(i-1,j)-psi(i-2,j));
                


            end
        end
    end
    
    
    for i=2:nx-1
        for j=4:ny-3
            if abs(psi(i,j))<10*dx   %% 4 grid points near the interface
                

                
                % y-direction
                vyn1(i,j) = 1/dy*(psi(i,j-2)-psi(i,j-3)); % 
                vyn2(i,j) = 1/dy*(psi(i,j-1)-psi(i,j-2));
                vyn3(i,j) = 1/dy*(psi(i,j)-psi(i,j-1));
                vyn4(i,j) = 1/dy*(psi(i,j+1)-psi(i,j));
                vyn5(i,j) = 1/dy*(psi(i,j+2)-psi(i,j+1));
                
                vyp1(i,j) = 1/dy*(psi(i,j+3)-psi(i,j+2)); % 
                vyp2(i,j) = 1/dy*(psi(i,j+2)-psi(i,j+1));
                vyp3(i,j) = 1/dy*(psi(i,j+1)-psi(i,j));
                vyp4(i,j) = 1/dy*(psi(i,j)-psi(i,j-1));
                vyp5(i,j) = 1/dy*(psi(i,j-1)-psi(i,j-2));

            end
        end
    end
    
    dir.vxn1 = vxn1;
    dir.vxn2 = vxn2;
    dir.vxn3 = vxn3;
    dir.vxn4 = vxn4;
    dir.vxn5 = vxn5;
    
    dir.vxp1 = vxp1;
    dir.vxp2 = vxp2;
    dir.vxp3 = vxp3;
    dir.vxp4 = vxp4;
    dir.vxp5 = vxp5;

    dir.vyn1 = vyn1;
    dir.vyn2 = vyn2;
    dir.vyn3 = vyn3;
    dir.vyn4 = vyn4;
    dir.vyn5 = vyn5;
    
    dir.vyp1 = vyp1;
    dir.vyp2 = vyp2;
    dir.vyp3 = vyp3;
    dir.vyp4 = vyp4;
    dir.vyp5 = vyp5;
end




