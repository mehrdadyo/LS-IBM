function [StateVar] = CheckforNegative(ControlVar,DOMAIN,VARIABLES,...
    StateVar,IBM,BC,Flux, LS)

negIndx = StateVar.phi<0 & LS.psi>0;
min(min(StateVar.phi(LS.psi>0)));
if nnz(negIndx) == 0
    StateVar.phi_old = StateVar.phi;
    disp(strcat('number of negative cells: ', num2str(nnz(negIndx)),...
        ' minimum of concentration: ',...
        num2str(min(min(StateVar.phi(LS.psi>0)))) ))
else
    [~,~,IBM_coeffP] = LSIBMcoeffs(IBM,DOMAIN,LS, StateVar.phi);
    [StateVar] = SolveTransportADRE(ControlVar,DOMAIN,VARIABLES,...
        StateVar,IBM,IBM_coeffP,BC,Flux, LS);
    StateVar.phi_old = StateVar.phi;
    disp(strcat('number of negative cells: ', num2str(nnz(negIndx)),...
        ' minimum of concentration: ',...
        num2str(min(min(StateVar.phi(LS.psi>0)))) ))
end




        
  