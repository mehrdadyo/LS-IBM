


StateVar.U = U;
StateVar.V = V;
StateVar.P = P;
StateVar.phi = phi;
StateVar.phi_old = phi;

ControlVar.time = time;
StateVar.U_old = StateVar.U;
StateVar.V_old = StateVar.V;
StateVar.P_old = StateVar.P;

StateVar.U_star = StateVar.U;
StateVar.V_star = StateVar.V;
StateVar.U_star_old = StateVar.U_star;
StateVar.V_star_old = StateVar.V_star;
LS.psi =psi;

[LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN);

[IBM_coeffU,IBM_coeffV,IBM_coeffP] = LSIBMcoeffs(IBM,DOMAIN,LS);
