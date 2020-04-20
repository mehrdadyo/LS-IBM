function psi = LSInitialize(DOMAIN, LSCase)

imax = DOMAIN.imax;
jmax = DOMAIN.jmax;

psi = zeros(imax+1,jmax+1);

if LSCase.case == 1 
    
    xc = LSCase.xc;
    yc = LSCase.yc;
    diamcyl = LSCase.diamcyl;
    
    for i=1:imax+1
        for j=1:jmax+1
            d = sqrt( (DOMAIN.xp(i)-xc).^2 + (DOMAIN.yp(j)-yc).^2 )-diamcyl/2;
            center = min(min(abs(d)));
            [I,J,~] = find(abs(d) == center);
            psi(i,j) = d(I(1),J(1));

        end
    end
    
elseif LSCase.case == 2  % fracture case
    ly = DOMAIN.ly ;     % width of fracture
    h = LSCase.h ;       % solid thickness
    d1 = DOMAIN.yp - h ; % distance from bottom wall of fracture
    d2 = (ly - h)- DOMAIN.yp; % distance from top wall of fracture
    
    d1 = ones(imax+1, jmax+1).*d1;  
    d2 = ones(imax+1, jmax+1).*d2;
    
    psi = double( abs(d1) <= abs(d2) ).* d1 + ...
        double( abs(d2) < abs(d1) ).* d2;
    
    
elseif LSCase.case == 3  % Simple rough fracture
    psi = simpleFractureInitialize(DOMAIN, LSCase.BoundaryCurve);
    
elseif LSCase.case == 4  % bi-mineral rough fracture
    psi = [];
    psi.psi1 = biMineralFractureInitialize(DOMAIN, LSCase.BoundaryCurve1);
    psi.psi2 = biMineralFractureInitialize(DOMAIN, LSCase.BoundaryCurve2);
end
    
end