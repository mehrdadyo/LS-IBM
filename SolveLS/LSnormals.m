function [LS] = LSnormals(LS,DOMAIN)


imax = DOMAIN.imax;
jmax = DOMAIN.jmax;

psi = LS.psi;

[nx,ny] = deal(zeros(imax+1,jmax+1));

dx = min(min(DOMAIN.dxp));
dy = min(min(DOMAIN.dyp));
epsi = 20*dx ;
for i =2:imax
    for j = 2:jmax
        
        if abs(psi(i,j))< epsi

            nx(i,j) = (psi(i+1,j)-psi(i-1,j))/(2*dx);
            ny(i,j) = (psi(i,j+1)-psi(i,j-1))/(2*dy);
            grad = sqrt(nx(i,j)^2 + ny(i,j)^2);
            if (grad > 1e-15)
                nx(i,j) = nx(i,j)/grad;
                ny(i,j) = ny(i,j)/grad;
            end
            
%             if abs(nx(i,j))<1e-10
%                 nx(i,j) = 0;
%                 ny(i,j) = 1;
%             end
% 
%             if abs(ny(i,j))<1e-10
%                 ny(i,j) = 0;
%                 nx(i,j) = 1;
%             end
            
        end

    end
end
% grad = sqrt(nx.^2 + ny.^2);
% 
% nx = nx./grad;
% ny = ny./grad;

LS.nx = nx;
LS.ny = ny;
