function [qx,qy]=Approximation(px,py,w,n,qx,qy,m,u,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Approximation by Jacobi-Bernstein Hybrid curves in paper Section 3
% Output Bezier points
% qx,qy are the point after satisfying the continuous condition
% % Author Mao Shi  All Right Reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=m-(u+v+2);

Px=px.*w;
Py=py.*w;

for k=0:M
    for ii=0:k
        A(k+1,ii+1)=(-1)^(k+ii)*nchoosek(k+u+1,ii)*nchoosek(k+v+1,k-ii)/nchoosek(k,ii);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first item on the left-hand side of equation (15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=0:M
    NNx(k+1)=0;
    NNy(k+1)=0;
     for i=0:n+k
       Nxx(k+1,i+1)=0;
       Nyy(k+1,i+1)=0;
       for j=max(0,i-n):min(k,i)
           Nxx(k+1,i+1)=Nxx(k+1,i+1)+nchoosek(k,j)*nchoosek(n,i-j)*A(k+1,j+1)*Px(i-j+1)/nchoosek(n+k,i);
           Nyy(k+1,i+1)=Nyy(k+1,i+1)+nchoosek(k,j)*nchoosek(n,i-j)*A(k+1,j+1)*Py(i-j+1)/nchoosek(n+k,i);
       end
       NNx(k+1)=NNx(k+1)+Nxx(k+1,i+1);
       NNy(k+1)=NNy(k+1)+Nyy(k+1,i+1);
     end
    NNx(k+1)=NNx(k+1)/(n+k+1);
    NNy(k+1)=NNy(k+1)/(n+k+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The second item on the left-hand side of equation (15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=0:m+n
    Cx(ii+1)=0;
    Cy(ii+1)=0;
    for jj=max(0,ii-n):min(m,ii)
        Cx(ii+1)=Cx(ii+1)+nchoosek(m,jj)*nchoosek(n,ii-jj)*qx(jj+1)*w(ii-jj+1)/nchoosek(m+n,ii);
        Cy(ii+1)=Cy(ii+1)+nchoosek(m,jj)*nchoosek(n,ii-jj)*qy(jj+1)*w(ii-jj+1)/nchoosek(m+n,ii);
    end
end

for k=0:M
    MMx(k+1)=0;
    MMy(k+1)=0;
     for i=0:n+m+k
       Mxx(k+1,i+1)=0;
       Myy(k+1,i+1)=0;
       for j=max(0,i-n-m):min(k,i)
           Mxx(k+1,i+1)=Mxx(k+1,i+1)+nchoosek(k,j)*nchoosek(n+m,i-j)*A(k+1,j+1)*Cx(i-j+1)/nchoosek(n+m+k,i);
           Myy(k+1,i+1)=Myy(k+1,i+1)+nchoosek(k,j)*nchoosek(n+m,i-j)*A(k+1,j+1)*Cy(i-j+1)/nchoosek(n+m+k,i);
       end
       MMx(k+1)=MMx(k+1)+Mxx(k+1,i+1);
       MMy(k+1)=MMy(k+1)+Myy(k+1,i+1);
     end
    MMx(k+1)=MMx(k+1)/(n+m+k+1);
    MMy(k+1)=MMy(k+1)/(n+m+k+1);
end

Lx=NNx-MMx;
Ly=NNy-MMy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The term on the right-hand side of equation (15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=0:M
    for i=0:n+k
       N(k+1,i+1)=0;
       for j=max(0,i-n):min(k,i)
           N(k+1,i+1)=N(k+1,i+1)+nchoosek(k,j)*nchoosek(n,i-j)*A(k+1,j+1)*w(i-j+1)/nchoosek(n+k,i);
       end
    end

    NNN=N(k+1,:);
    for l=0:k+n
        for j=0:M
            AA=0;
            for i=0:j
                AA=AA+(-1)^(i+j)*nchoosek(j+u+1,i)*nchoosek(j+v+1,j-i)*nchoosek(n+k,l)/(nchoosek(n+k+u+1+v+1+j,l+u+1+i));
            end
            BJ(l+1,j+1)=AA/(n+k+u+1+v+1+j+1);
        end
    end
    nmkk=NNN*BJ;
    nmk(k+1,1:length(nmkk))=nmkk;
end
BB=nmk;

CCX=Lx/BB;
CCY=Ly/BB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From Jacobi control points to Bezier control points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LL1=JacobiToBernstein(M,u+1,v+1);
 
 QQx=CCX*LL1;
 QQy=CCY*LL1;
 
 for i=0:M
     qx(u+1+i+1)=nchoosek(m-(u+v+2),i)*QQx(i+1)/nchoosek(m,u+1+i);
     qy(u+1+i+1)=nchoosek(m-(u+v+2),i)*QQy(i+1)/nchoosek(m,u+1+i);
 end
end

