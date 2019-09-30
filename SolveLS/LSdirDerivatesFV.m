%% =======================
% This function finds the first order derivatives in 5 point stencil
% The inputs are:
% - grid sizes in x and y direction
% - the level set matrix
% - the equation type: level set equation or reinitialization equation


%%


function [dir] = LSdirDerivatesFV(dx,psi,equation,StateVar)


%% 
[nx,ny] = size(psi);
u = StateVar.U;
v = StateVar.V;
%% direction

[psi1e,psi2e,psi3e,psi4e,psi5e,psi1w,psi2w,psi3w,psi4w,psi5w] = ...
    deal(zeros(nx,ny));

[psi1n,psi2n,psi3n,psi4n,psi5n,psi1s,psi2s,psi3s,psi4s,psi5s] = ...
    deal(zeros(nx,ny));

[psi1wN,psi2wN,psi3wN,psi4wN,psi5wN, ...
    psi1wP,psi2wP,psi3wP,psi4wP,psi5wP] = deal(zeros(nx,ny));

[psi1eN,psi2eN,psi3eN,psi4eN,psi5eN, ...
    psi1eP,psi2eP,psi3eP,psi4eP,psi5eP] = deal(zeros(nx,ny));

[psi1sN,psi2sN,psi3sN,psi4sN,psi5sN, ...
    psi1sP,psi2sP,psi3sP,psi4sP,psi5sP] = deal(zeros(nx,ny));

[psi1nN,psi2nN,psi3nN,psi4nN,psi5nN, ...
    psi1nP,psi2nP,psi3nP,psi4nP,psi5nP] = deal(zeros(nx,ny));

% tic;
if equation == "LevelSetEqn"
    for i=2:nx-1
        for j=2:ny-1
            if abs(psi(i,j))<4*dx   %% 4 grid points near the interface
                if u(i,j) >= 0 % East
                % x-direction
                    psi1e(i,j) = psi(i-2,j); % 
                    psi2e(i,j) = psi(i-1,j);
                    psi3e(i,j) = psi(i,j);
                    psi4e(i,j) = psi(i+1,j);
                    psi5e(i,j) = psi(i+2,j);
                elseif u(i,j)<0
                    psi1e(i,j) = psi(i+3,j); % 
                    psi2e(i,j) = psi(i+2,j);
                    psi3e(i,j) = psi(i+1,j);
                    psi4e(i,j) = psi(i,j);
                    psi5e(i,j) = psi(i+1,j);
                end
                
                if u(i-1,j) >= 0 % West
                % x-direction
                    psi1w(i,j) = psi(i-3,j); % 
                    psi2w(i,j) = psi(i-2,j);
                    psi3w(i,j) = psi(i-1,j);
                    psi4w(i,j) = psi(i,j);
                    psi5w(i,j) = psi(i+1,j);
                elseif u(i,j)<0
                    psi1w(i,j) = psi(i+2,j); % 
                    psi2w(i,j) = psi(i+1,j);
                    psi3w(i,j) = psi(i,j);
                    psi4w(i,j) = psi(i-1,j);
                    psi5w(i,j) = psi(i-2,j);
                end
                
                % y-direction
                if v(i,j) >= 0 % North
                % y-direction
                    psi1n(i,j) = psi(i,j-2); % 
                    psi2n(i,j) = psi(i,j-1);
                    psi3n(i,j) = psi(i,j);
                    psi4n(i,j) = psi(i,j+1);
                    psi5n(i,j) = psi(i,j+2);
                elseif v(i,j)<0
                    psi1n(i,j) = psi(i,j+3); % 
                    psi2n(i,j) = psi(i,j+2);
                    psi3n(i,j) = psi(i,j+1);
                    psi4n(i,j) = psi(i,j);
                    psi5n(i,j) = psi(i,j-1);
                end
                
                if v(i,j-1) >= 0 % South
                % y-direction
                    psi1s(i,j) = psi(i,j-3); % 
                    psi2s(i,j) = psi(i,j-2);
                    psi3s(i,j) = psi(i,j-1);
                    psi4s(i,j) = psi(i,j);
                    psi5s(i,j) = psi(i,j+1);
                elseif v(i,j-1)<0
                    psi1s(i,j) = psi(i,j+2); % 
                    psi2s(i,j) = psi(i,j+1);
                    psi3s(i,j) = psi(i,j);
                    psi4s(i,j) = psi(i,j-1);
                    psi5s(i,j) = psi(i,j-2);
                end

                
            end
        end
    end
    
    dir = struct(...
        'psi1w',psi1w,'psi2w',psi2w,'psi3w',psi3w,'psi4w',psi4w,...
        'psi5w',psi5w,...
        'psi1e',psi1e,'psi2e',psi2e,'psi3e',psi3e,'psi4e',psi4e,...
        'psi5e',psi5e,...
        'psi1s',psi1s,'psi2s',psi2s,'psi3s',psi3s,'psi4s',psi4s,...
        'psi5s',psi5s,...
        'psi1n',psi1n,'psi2n',psi2n,'psi3n',psi3n,'psi4n',psi4n,...
        'psi5n',psi5n);
            
    
elseif equation == "ReinitializationEqn"
    for i=2:nx-1
        for j=2:ny-1
            if abs(psi(i,j))<4*dx   %% 4 grid points near the interface
                
                % Esat face
                psi1eN(i,j) = psi(i-2,j); % 
                psi2eN(i,j) = psi(i-1,j);
                psi3eN(i,j) = psi(i,j);
                psi4eN(i,j) = psi(i+1,j);
                psi5eN(i,j) = psi(i+2,j);
                
                psi1eP(i,j) = psi(i+3,j);% 
                psi2eP(i,j) = psi(i+2,j);
                psi3eP(i,j) = psi(i+1,j);
                psi4eP(i,j) = psi(i,j);
                psi5eP(i,j) = psi(i+1,j);
                
                
                % West face
                psi1wN(i,j) = psi(i-3,j);  % 
                psi2wN(i,j) = psi(i-2,j);
                psi3wN(i,j) = psi(i-1,j);
                psi4wN(i,j) = psi(i,j);
                psi5wN(i,j) = psi(i+1,j);
                
                psi1wP(i,j) = psi(i+2,j);% 
                psi2wP(i,j) = psi(i+1,j);
                psi3wP(i,j) = psi(i,j);
                psi4wP(i,j) = psi(i-1,j);
                psi5wP(i,j) = psi(i-2,j);              
                            
                % North face
                psi1nN(i,j) = psi(i,j-2); % 
                psi2nN(i,j) = psi(i,j-1);
                psi3nN(i,j) = psi(i,j);
                psi4nN(i,j) = psi(i,j+1);
                psi5nN(i,j) = psi(i,j+2);
                
                psi1nP(i,j) = psi(i,j+3);% 
                psi2nP(i,j) = psi(i,j+2);
                psi3nP(i,j) = psi(i,j+1);
                psi4nP(i,j) = psi(i,j);
                psi5nP(i,j) = psi(i,j-1);
                
                
                % South face
                psi1sN(i,j) = psi(i,j-3); % 
                psi2sN(i,j) = psi(i,j-2);
                psi3sN(i,j) = psi(i,j-1);
                psi4sN(i,j) = psi(i,j);
                psi5sN(i,j) = psi(i,j+1);
                
                psi1sP(i,j) = psi(i,j+2);% 
                psi2sP(i,j) = psi(i,j+1);
                psi3sP(i,j) = psi(i,j);
                psi4sP(i,j) = psi(i,j-1);
                psi5sP(i,j) = psi(i,j-2);       

            end
        end
    end
    dir = struct(...
        'psi1eN',psi1eN,'psi2eN',psi2eN,'psi3eN',psi3eN,'psi4eN',psi4eN,...
        'psi5eN',psi5eN,...
        'psi1eP',psi1eP,'psi2eP',psi2eP,'psi3eP',psi3eP,'psi4eP',psi4eP,...
        'psi5eP',psi5eP,...
        'psi1wN',psi1wN,'psi2wN',psi2wN,'psi3wN',psi3wN,'psi4wN',psi4wN,...
        'psi5wN',psi5wN,...
        'psi1wP',psi1wP,'psi2wP',psi2wP,'psi3wP',psi3wP,'psi4wP',psi4wP,...
        'psi5wP',psi5wP,...
        'psi1nN',psi1nN,'psi2nN',psi2nN,'psi3nN',psi3nN,'psi4nN',psi4nN,...
        'psi5nN',psi5nN,...
        'psi1nP',psi1nP,'psi2nP',psi2nP,'psi3nP',psi3nP,'psi4nP',psi4nP,...
        'psi5nP',psi5nP,...
        'psi1sN',psi1sN,'psi2sN',psi2sN,'psi3sN',psi3sN,'psi4sN',psi4sN,...
        'psi5sN',psi5sN,...
        'psi1sP',psi1sP,'psi2sP',psi2sP,'psi3sP',psi3sP,'psi4sP',psi4sP,...
        'psi5sP',psi5sP);
   
end


