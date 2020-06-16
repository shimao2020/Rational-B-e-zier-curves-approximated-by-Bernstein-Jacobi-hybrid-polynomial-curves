%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% My Example 3 in Paper
%Stanislaw Lewanowicz 2012 Numer Algor Example3
%Qianqian Hu    Journal of Computational and Applied Mathematics 2013  Example3
% Author Mao Shi  All Right Reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
hold on
format long g

n=9;
px=[17,32,-23,33,-23,25,30,-5,-5,11];
py=[12,34,24,62,15,3,-2,-8,15,8];
w=[1,2,3,6,4,5,3,4,2,1];
m=18;
Px=px.*w;
Py=py.*w;
Step=0.001;

for i=0:m
   qx(i+1)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw Rational Curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[drax,dray]=DiscreRationalCurves(px,py,w,n,Step);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Condition of continuity between  rational Bezier curve and bezier curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=0;
v=0;

qx=Con_u_v(px,w,n,m,u,v,Px);
qy=Con_u_v(py,w,n,m,u,v,Py);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Approximation by Jacobi-Bernstein Hybrid curves without degree reduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[qx,qy]=Approximation(px,py,w,n,qx,qy,m,u,v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw  Curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[drx,dry]=DiscreCurves(qx,qy,m,Step);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Degree Reduction from m to mo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mo=10;
cqx=Cont_Curve_u_v(qx,m,mo,u,v);
cqy=Cont_Curve_u_v(qy,m,mo,u,v);

for i=0:m
    ww(i+1)=1;
end
[dqx,dqy]=Approximation(qx,qy,ww,m,cqx,cqy,mo,u,v);
disp('Control points of Bezier curve ')
sprintf('%0.7f,',dqx)
sprintf('%0.7f,',dqy)

[dqrx,dqry]=DiscreCurves(dqx,dqy,mo,Step);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Hausdoff distance for rational Bezier curve of degree n and Bezier curve
%  of degree mo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hd=HauDis(dqrx,dqry,drax,dray,Step);
 sprintf('Hausdoff distance Hd=%0.7f',Hd)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw graphic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 plot(drax,dray,'k',drx,dry,'G-.',dqrx,dqry,'R');
 legend('The given curve','Bezier curve with degree of m',' Bezier curve from degree reduction from m to mo','Location','SouthOutside');


 