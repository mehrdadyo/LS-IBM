function [g1,g2]=QuickInterp(x,xU,xUU,xD)
    
g1=((x-xU).*(x-xUU))./((xD-xU).*(xD-xUU));
g2=((x-xU).*(xD-x))./((xU-xUU).*(xD-xUU));

