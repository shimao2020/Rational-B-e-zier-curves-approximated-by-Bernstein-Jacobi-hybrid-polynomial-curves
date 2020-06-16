function [x,y]=DiscreCurves(px,py,n,Step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Discrete  curves
% Author Mao Shi  All Right Reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0;
for u=0:Step:1
    i=i+1;
    x(i)=deCasteljau(px,n,u);
    y(i)=deCasteljau(py,n,u);
end
end