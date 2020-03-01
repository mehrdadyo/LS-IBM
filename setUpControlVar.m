function ControlVar = setUpControlVar(VARIABLES, DOMAIN)

%% Control VAriables
ControlVar.PISO = 1;
ControlVar.time=0;
ControlVar.timedt = 50000;
ControlVar.savedat=50;
ControlVar.rat=floor(0.01/VARIABLES.dt);

ControlVar.tol=1e-3;
ControlVar.f=0;
ControlVar.tolbicg=10^(-5); %% max tol for bicon
ControlVar.maxit=2000;   %% max iter for bicon 

ControlVar.flow_steady= 1;
ControlVar.disc_scheme_vel = 2;
StateVar.P_cor_vec = zeros(1,(DOMAIN.imax-1)*(DOMAIN.jmax-1))';

ControlVar.flow_steady = 0;

ControlVar.tol_q=1e-5;
ControlVar.tolbicg_c=5e-5;
ControlVar.maxit_c=2000;





end