function psi = biMineralFractureInitialize(DOMAIN, BoundaryCurve)

xt = BoundaryCurve.xt;
xb = BoundaryCurve.xb;
yt = BoundaryCurve.yt;
yb = BoundaryCurve.yb;
fbTop = BoundaryCurve.boundfunctionTop;
fbBottom = BoundaryCurve.bFunctionBottom;

xt_s = BoundaryCurve.xt_s;
xt_e = BoundaryCurve.xt_e;
xb_s = BoundaryCurve.xb_s;
xb_e = BoundaryCurve.xb_e;

endBotIndx = BoundaryCurve.endBotIndx;

b = BoundaryCurve.b;
lambda = BoundaryCurve.lambda;
a = BoundaryCurve.a;
x_0 = BoundaryCurve.x_0;

X = DOMAIN.Xp;
Y = DOMAIN.Yp - DOMAIN.ly/2;
psi = zeros(size(X));
% profile on
% dispIndx = ceil(length(x)/20);
for i =1:length(DOMAIN.xp)
%     if ~mod(i,dispIndx)
%         display(strcat('...', num2str(i/length(x) *100), '.............'))
%     end
    indx = abs(xb-X(i,1)) < DOMAIN.lx/2;
    xb_temp = xb(indx);
    yb_temp = yb(indx);

    indx = abs(xt-X(i,1)) < DOMAIN.lx/2;
    xt_temp = xt(indx);
    yt_temp = yt(indx);
    %%%%% bottom
    dist = sqrt((xb_temp' - X(i,1:endBotIndx)).^2 + ...
        (yb_temp' - Y(i,1:endBotIndx)).^2 );
    y_curve = fbBottom(X(i,1), a, lambda, x_0, b);
    
    signNeg = (Y(i,1:endBotIndx) <= y_curve) * ...
        any((X(i,1)>= xb_s).*(X(i,1)<= xb_e));
    
    psi(i,1:endBotIndx) = min(dist) .* (1 - signNeg) - ...
         min(dist) .* signNeg ;

    %%%%% top
    dist = sqrt((xt_temp' - X(i,endBotIndx+1:end)).^2 + ...
        (yt_temp' - Y(i,endBotIndx+1:end)).^2 );
    y_curve = fbTop(X(i,1), a, lambda, x_0, b);
    
    signNeg = (Y(i,endBotIndx+1:end) >= y_curve) * ...
        any((X(i,1)>= xt_s).*(X(i,1)<= xt_e));
    
    psi(i,endBotIndx+1:end) = -min(dist) .* signNeg + ...
         min(dist) .* (1 - signNeg) ; 

end

