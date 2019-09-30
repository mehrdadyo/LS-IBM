function [Diffu]=DiffFlux(VARIABLES,flag,DOMAIN,BC)

dxu = DOMAIN.dxu;
dyu = DOMAIN.dyu;
dxv = DOMAIN.dxv;
dyv = DOMAIN.dyv;
dxp = DOMAIN.dxp;
dyp = DOMAIN.dyp;
imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
xu = DOMAIN.xu;
xp = DOMAIN.xp;
dV = DOMAIN.dV_p;
dt = VARIABLES.dt;


% flag=1 ==> U
% flag=-1 ==> V
% flag=0 ==> Scalar

if flag==1          %CoEW=CoEWu  CoNS=CoNSu
    BC_e_p = BC.BC_e_p;

    Gamma=VARIABLES.Re;
    
    De=dyv(2:imax-1,1:jmax-1)./(Gamma*dxu(2:imax-1,2:jmax));
    De=[zeros(1,size(De,2)+2);[zeros(size(De,1),1) De ...
        zeros(size(De,1),1)];zeros(1,size(De,2)+2)];  % add BC zeros around 
    

    
    Dw=dyv(2:imax-1,1:jmax-1)./(Gamma*dxu(1:imax-2,2:jmax));
    Dw=[zeros(1,size(Dw,2)+2);[zeros(size(Dw,1),1) Dw ...
        zeros(size(Dw,1),1)];zeros(1,size(Dw,2)+2)];  % add BC zeros around 
    
    Dn=dxv(2:imax-1,2:jmax)./(Gamma*dyu(2:imax-1,2:jmax));
    Dn=[zeros(1,size(Dn,2)+2);[zeros(size(Dn,1),1) Dn ...
        zeros(size(Dn,1),1)];zeros(1,size(Dn,2)+2)];  % add BC zeros around 
    
    Ds=dxv(2:imax-1,2:jmax)./(Gamma*dyu(2:imax-1,1:jmax-1));
    Ds=[zeros(1,size(Ds,2)+2);[zeros(size(Ds,1),1) Ds ...
        zeros(size(Ds,1),1)];zeros(1,size(Ds,2)+2)];  % add BC zeros around 
    
    if BC_e_p == 1
        %De(end,2:end-1) = dyv(imax,1:jmax-1)./(Gamma*dxu(imax,2:jmax));
        Dw(end,2:end-1) = dyv(imax,1:jmax-1)./(Gamma*dxu(imax-1,2:jmax));
        Dn(end,2:end-1) = dxv(imax,2:jmax)./(Gamma*dyu(imax,2:jmax));
        Ds(end,2:end-1) = dxv(imax,2:jmax)./(Gamma*dyu(imax,1:jmax-1));
    end

    
    Diffu = struct('De_u',De,'Dw_u',Dw,'Dn_u',Dn,'Ds_u',Ds);
    
elseif flag==-1         
    Gamma=VARIABLES.Re;
    
    De=dyu(2:imax,2:jmax-1)./(Gamma*dxv(2:imax,2:jmax-1));
    De=[zeros(1,size(De,2)+2);[zeros(size(De,1),1) De ...
        zeros(size(De,1),1)];zeros(1,size(De,2)+2)];  % add BC zeros around 
    
    Dw=dyu(2:imax,2:jmax-1)./(Gamma*dxv(1:imax-1,2:jmax-1));
    Dw=[zeros(1,size(Dw,2)+2);[zeros(size(Dw,1),1) Dw ...
        zeros(size(Dw,1),1)];zeros(1,size(Dw,2)+2)];  % add BC zeros around 
    
    Dn=dxu(1:imax-1,2:jmax-1)./(Gamma*dyv(2:imax,2:jmax-1));
    Dn=[zeros(1,size(Dn,2)+2);[zeros(size(Dn,1),1) Dn ...
        zeros(size(Dn,1),1)];zeros(1,size(Dn,2)+2)];  % add BC zeros around 
    
    Ds=dxu(1:imax-1,2:jmax-1)./(Gamma*dyv(2:imax,1:jmax-2));
    Ds=[zeros(1,size(Ds,2)+2);[zeros(size(Ds,1),1) Ds ...
        zeros(size(Ds,1),1)];zeros(1,size(Ds,2)+2)];  % add BC zeros around
  
    Diffu = struct('De_v',De,'Dw_v',Dw,'Dn_v',Dn,'Ds_v',Ds);
    
    
    
     
    
    
elseif flag==0
    Gamma=VARIABLES.D;
    
    De = Gamma*dyv(2:imax,1:jmax-1)./dxp(2:imax,2:jmax);
     De=[zeros(1,size(De,2)+2);[zeros(size(De,1),1) De ...
        zeros(size(De,1),1)];zeros(1,size(De,2)+2)];  % add BC zeros around 
     
     Dw = Gamma*dyv(2:imax,1:jmax-1)./dxp(1:imax-1,2:jmax);
     Dw=[zeros(1,size(Dw,2)+2);[zeros(size(Dw,1),1) Dw ...
        zeros(size(Dw,1),1)];zeros(1,size(Dw,2)+2)];  % add BC zeros around 
     
     Dn = Gamma*dxu(1:imax-1,2:jmax)./dyp(2:imax,2:jmax);
     Dn=[zeros(1,size(Dn,2)+2);[zeros(size(Dn,1),1) Dn ...
        zeros(size(Dn,1),1)];zeros(1,size(Dn,2)+2)];  % add BC zeros around 
     
     Ds = Gamma*dxu(1:imax-1,2:jmax)./dyp(2:imax,1:jmax-1);
     Ds=[zeros(1,size(Ds,2)+2);[zeros(size(Ds,1),1) Ds ...
        zeros(size(Ds,1),1)];zeros(1,size(Ds,2)+2)];  % add BC zeros around
    
    [D1_a,D2_a]=deal(zeros(1,jmax+1));

    
    D1=Gamma*(xu(1)-xp(2))/((xp(3)-xp(2))*(xp(3)-xu(1)));
    D2=Gamma*(xp(3)-2*xu(1)+xp(2))/((xp(2)-xu(1))*(xp(3)-xu(1)));  

    
    D1_a(1,2:jmax)= D1*dxv(1,1:jmax-1);
    D2_a(1,2:jmax)= D2*dxv(1,1:jmax-1);

    A0_p=dV/dt;
    
    Diffu = struct('De_p',De,'Dw_p',Dw,'Dn_p',Dn,'Ds_p',Ds,'D1_a',D1_a,...
            'D2_a',D2_a,'A0_p_p',A0_p);
    
           
end

