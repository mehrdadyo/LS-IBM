function [ConvF] = ConvFlux(StateVar,ControlVar,DOMAIN,VARIABLES, BC)

U = StateVar.U;
V = StateVar.V;

[imax,jmax]=size(U);
jmax=jmax-1;


dt = VARIABLES.dt;
% flag=1 ==> U
% flag=-1 ==> V
% flag=0 ==> Scalar


flag = ControlVar.flc;

if flag==1          %CoEW=CoEWu  CoNS=CoNSu   dy=dyv   dx=dxv
    BC_e_p = BC.BC_e_p;
   
    CoEW = DOMAIN.CoEWu;
    CoNS = DOMAIN.CoNSu;
    dy = DOMAIN.dyv;
    dx = DOMAIN.dxv;
    dV = DOMAIN.dV_u;
    
    [Fe,Fw,Fn,Fs,alphae,alphaw,alphan,alphas]=deal(zeros(imax,jmax+1));

    Fe2=(CoEW(2:imax-1,2:jmax).*U(3:imax,2:jmax)+...
        (1-CoEW(2:imax-1,2:jmax)).*U(2:imax-1,2:jmax))...
            .*dy(2:imax-1,1:jmax-1);
    Fw2=(CoEW(1:imax-2,2:jmax).*U(2:imax-1,2:jmax)+...
        (1-CoEW(1:imax-2,2:jmax)).*U(1:imax-2,2:jmax))...
            .*dy(2:imax-1,1:jmax-1);
    Fn2=(CoNS(2:imax-1,2:jmax).*V(3:imax,2:jmax)+...
        (1-CoNS(2:imax-1,2:jmax)).*V(2:imax-1,2:jmax))...
            .*dx(2:imax-1,2:jmax);    
    Fs2=(CoNS(2:imax-1,2:jmax).*V(3:imax,1:jmax-1)+...
        (1-CoNS(2:imax-1,1:jmax-1)).*V(2:imax-1,1:jmax-1))...
            .*dx(2:imax-1,2:jmax);
        
    
        
        
%     alphae2=floor((sign(Fe2)+1)./2);
%     alphaw2=floor((sign(Fw2)+1)./2);
%     alphan2=floor((sign(Fn2)+1)./2);
%     alphas2=floor((sign(Fs2)+1)./2);
    
    Fe(2:end-1,2:end-1)=Fe2;
    Fw(2:end-1,2:end-1)=Fw2;
    Fn(2:end-1,2:end-1)=Fn2;
    Fs(2:end-1,2:end-1)=Fs2;
    
    if BC_e_p == 1
        Fe(end,2:end-1) = U(imax,2:jmax).*dy(imax,1:jmax-1);
        Fw(end,2:end-1) = (CoEW(imax-1,2:jmax).*U(imax,2:jmax)+...
            (1-CoEW(imax-1,2:jmax)).*U(imax-1,2:jmax))...
                .*dy(imax,1:jmax-1);
        Fn(end,2:end-1) = (CoNS(imax,2:jmax).*V(imax+1,2:jmax)+...
            (1-CoNS(imax,2:jmax)).*V(imax,2:jmax))...
            .*dx(imax,2:jmax);
        Fs(end,2:end-1) = (CoNS(imax,2:jmax).*V(imax+1,1:jmax-1)+...
            (1-CoNS(imax,1:jmax-1)).*V(imax,1:jmax-1))...
            .*dx(imax,2:jmax);
    end
    
    alphae = floor((sign(Fe)+1)./2);
    alphaw = floor((sign(Fw)+1)./2);
    alphan = floor((sign(Fn)+1)./2);
    alphas = floor((sign(Fs)+1)./2);
    
%     alphae(2:end-1,2:end-1)=alphae2; 
%     alphaw(2:end-1,2:end-1)=alphaw2;
%     alphan(2:end-1,2:end-1)=alphan2;
%     alphas(2:end-1,2:end-1)=alphas2;
    
    dF=Fe-Fw+Fn-Fs;
    A0_p=dV/dt;
    if ControlVar.flow_steady 
        A0_p = A0_p * 0;
    end
    
    ConvF = struct('Fe_u',Fe,'Fw_u',Fw,'Fn_u',Fn,'Fs_u',Fs,'dF_u',dF,...
        'A0_p_u',A0_p,'alphae_u',alphae,'alphaw_u',alphaw,...
        'alphan_u',alphan,'alphas_u',alphas);
elseif flag==-1     %CoEW=CoEWv  CoNS=CoNSv   dy=dyu   dx=dxu 
    CoEW = DOMAIN.CoEWv;
    CoNS = DOMAIN.CoNSv;
    dy = DOMAIN.dyu;
    dx = DOMAIN.dxu;
    dV = DOMAIN.dV_v;
    
    [Fe,Fw,Fn,Fs,alphae,alphaw,alphan,alphas]=deal(zeros(imax+1,jmax));
    
    Fe2=(CoEW(2:imax,2:jmax-1).*U(2:imax,3:jmax)+...
        (1-CoEW(2:imax,2:jmax-1)).*U(2:imax,2:jmax-1))...
            .*dy(2:imax,2:jmax-1);
    Fw2=(CoEW(2:imax,2:jmax-1).*U(1:imax-1,3:jmax)+...
        (1-CoEW(2:imax,2:jmax-1)).*U(1:imax-1,2:jmax-1))...
            .*dy(2:imax,2:jmax-1);
    Fn2=(CoNS(2:imax,2:jmax-1).*V(2:imax,3:jmax)+...
        (1-CoNS(2:imax,2:jmax-1)).*V(2:imax,2:jmax-1))...
            .*dx(1:imax-1,2:jmax-1);
    Fs2=(CoNS(2:imax,1:jmax-2).*V(2:imax,2:jmax-1)+...
        (1-CoNS(2:imax,1:jmax-2)).*V(2:imax,1:jmax-2))...
            .*dx(1:imax-1,2:jmax-1); 
        
    alphae2=floor((sign(Fe2)+1)./2);
    alphaw2=floor((sign(Fw2)+1)./2);
    alphan2=floor((sign(Fn2)+1)./2);
    alphas2=floor((sign(Fs2)+1)./2);
    
    Fe(2:end-1,2:end-1)=Fe2;
    Fw(2:end-1,2:end-1)=Fw2;
    Fn(2:end-1,2:end-1)=Fn2;
    Fs(2:end-1,2:end-1)=Fs2;
    
    alphae(2:end-1,2:end-1)=alphae2;
    alphaw(2:end-1,2:end-1)=alphaw2;    
    alphan(2:end-1,2:end-1)=alphan2;
    alphas(2:end-1,2:end-1)=alphas2;
    
    dF=Fe-Fw+Fn-Fs;
    A0_p=dV/dt;
    if ControlVar.flow_steady 
        A0_p = A0_p * 0;
    end
    
    
    ConvF = struct('Fe_v',Fe,'Fw_v',Fw,'Fn_v',Fn,'Fs_v',Fs,'dF_v',dF,...
        'A0_p_v',A0_p,'alphae_v',alphae,'alphaw_v',alphaw,...
        'alphan_v',alphan,'alphas_v',alphas);
    
elseif flag==0       %dy=dyv          dx=dxu
    
    dy = DOMAIN.dyv;
    dx = DOMAIN.dxu;
    dV = DOMAIN.dV_p;
    
    [Fe,Fw,Fn,Fs,alphae,alphaw,alphan,alphas]=deal(zeros(imax+1,jmax+1));
    Fe2=U(2:imax,2:jmax).*dy(2:imax,1:jmax-1);
    Fw2=U(1:imax-1,2:jmax).*dy(2:imax,1:jmax-1);
    Fn2=V(2:imax,2:jmax).*dx(1:imax-1,2:jmax);
    Fs2=V(2:imax,1:jmax-1).*dx(1:imax-1,2:jmax);
    
    alphae2=floor((sign(Fe2)+1)./2);
    alphaw2=floor((sign(Fw2)+1)./2);
    alphan2=floor((sign(Fn2)+1)./2);
    alphas2=floor((sign(Fs2)+1)./2);
    
    Fe(2:end-1,2:end-1)=Fe2;
    Fw(2:end-1,2:end-1)=Fw2;
    Fn(2:end-1,2:end-1)=Fn2;
    Fs(2:end-1,2:end-1)=Fs2;
   
    alphaw(2:end-1,2:end-1)=alphaw2;
    alphae(2:end-1,2:end-1)=alphae2; 
    alphan(2:end-1,2:end-1)=alphan2; 
    alphas(2:end-1,2:end-1)=alphas2; 
    
    dF=Fe-Fw+Fn-Fs;
    A0_p=dV/dt;
    
    if ControlVar.transport_steady
        A0_p = A0_p * 0;
    end
    
    ConvF = struct('Fe_p',Fe,'Fw_p',Fw,'Fn_p',Fn,'Fs_p',Fs,'dF_p',dF,...
        'A0_p_p',A0_p,'alphae_p',alphae,'alphaw_p',alphaw,...
        'alphan_p',alphan,'alphas_p',alphas);
   
end






