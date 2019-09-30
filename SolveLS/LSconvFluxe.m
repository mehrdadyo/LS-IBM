function [Fe,Fw,Fn,Fs] = LSconvFluxe(u,v,dx,dy,imax,jmax)


[Fe,Fw,Fn,Fs] = deal(zeros(imax+1,jmax+1));
Fe2 = u(2:imax,2:jmax).*dy(2:imax,1:jmax-1);
Fw2 = u(1:imax-1,2:jmax).*dy(2:imax,1:jmax-1);
Fn2 = v(2:imax,2:jmax).*dx(1:imax-1,2:jmax);
Fs2 = v(2:imax,1:jmax-1).*dx(1:imax-1,2:jmax);


Fe(2:end-1,2:end-1) = Fe2;
Fw(2:end-1,2:end-1) = Fw2;
Fn(2:end-1,2:end-1)= Fn2;
Fs(2:end-1,2:end-1) = Fs2;