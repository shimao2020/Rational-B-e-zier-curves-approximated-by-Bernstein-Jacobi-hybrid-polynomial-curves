function [qq]=Con_u_v(p,w,n,m,u,v,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Equations 9 and 10
%p--Given control point
%n--Degree of Given curve
%m--Degree of unknown curve
%u--Order of derivation
%v--Order of derivation
%P=p.*w--Multiplication of control points and weights
% u+v+1<n
%Theorem 1 in Paper
% Author Mao Shi  All Right Reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 q(1)=p(1);
 for r=1:u
    S1=0;
    for i=0:r-1
       S1=S1+(-1)^(r+i)*anchoosek(m+n,r)*anchoosek(r,i)*Ck(m,q,n,w,i);
      end
     S2=0;
     for i=max(0,r-n):r-1
          S2=S2+anchoosek(m,i)*anchoosek(n,r-i)*w(r-i+1)*q(i+1);
     end
      q(r+1)=(anchoosek(n,r)*FDeltaP(P,r,n)-S1-S2)/(w(1)*anchoosek(m,r));        
 end
 q(m+1)=p(n+1);
for s=1:v
    S1=0;
    for i=1:s
        S1=S1+(-1)^(s+i)*anchoosek(m+n,s)*anchoosek(s,i)*Ck(m,q,n,w,m+n-s+i);
    end
    S2=0;
    for i=1:min(s,n)
        S2=S2+(-1)^s*anchoosek(m,s-i)*anchoosek(n,i)*w(n-i+1)*q(m-s+i+1);
    end
    q(m-s+1)=(-1)^s*(anchoosek(n,s)*BDeltaP(P,s,n)-S1-S2)/(w(n+1)*anchoosek(m,s));
end
qq=q;
end