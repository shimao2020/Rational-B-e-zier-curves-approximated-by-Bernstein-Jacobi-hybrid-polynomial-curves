function [x,y]=DiscreRationalCurves(px,py,w,n,Step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Discrete rational curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0;
for u=0:Step:1
    i=i+1;
    x(i)=RadeCasteljau(px,w,n,u);
    y(i)=RadeCasteljau(py,w,n,u);
end
end
