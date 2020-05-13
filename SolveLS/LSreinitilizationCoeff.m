%% ============== 

%%

function [G]= LSreinitilizationCoeff(psinp,psi)



a = psinp.psi_xn;
b = psinp.psi_xp;
c = psinp.psi_yn;
d = psinp.psi_yp;


aP = max(a,0);
aM = min(a,0);

bP = max(b,0);
bM = min(b,0);

cP = max(c,0);
cM = min(c,0);

dP = max(d,0);
dM = min(d,0);

Gp = sqrt(max(aP.^2,bM.^2)+max(cP.^2,dM.^2));

Gn = sqrt(max(aM.^2,bP.^2)+max(cM.^2,dP.^2));

IdP = psi >= 0;
IdN = psi < 0;

G = Gp.*IdP + Gn.*IdN;





