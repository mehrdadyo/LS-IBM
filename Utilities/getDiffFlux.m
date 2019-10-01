function Flux = getDiffFlux(VARIABLES, DOMAIN,BC)


%% U-momentum diffusive flux coeff
Flux.fl=1;
[Flux.Diffu_U]=DiffFlux(VARIABLES,Flux.fl,DOMAIN,BC);

%% V-momentum diffusive flux coeff
Flux.fl=-1;
[Flux.Diffu_V]=DiffFlux(VARIABLES,Flux.fl,DOMAIN);

%% Scalar Transport diffusive flux coeff
Flux.fl=0;
[Flux.Diffu_P]=DiffFlux(VARIABLES,Flux.fl,DOMAIN);