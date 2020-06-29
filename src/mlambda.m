function [f,s,Z]=mlambda(a,Q,n)
% MLAMBDA : LAMBDA/MLAMBDA integer ambiguity resolution
%
% Resolve integer amibiguity by LAMBDA/MLAMBDA integer least-square estimation.
% Float solution is preprocessed by LAMBDA decorration/reduction and fixed
% (integer) solutions are searched in the transformed space by MLAMBDA method.
%
% Reference :
% [1] P.Teunissen, The least-square ambiguity decorrelation adjustment: a method
%     for fast GPS ambiguity estimation, J.Geodesy 70:65-82, 1995
% [2] X.-W.Chang et al., MLAMBDA: A modified LAMBDA method for integer least-
%     squares estimation, ION AM, 2005
%
%          Copyright (C) 2006 by T.TAKASU, All rights reserved.
%
% argin  : a   = float solution
%          Q   = float solution varience/covarience matrix
%          n   = number of fixed solutions
%
% argout : f   = fixed solutions [f1,f2,..,fn]
%          s   = residuals [s1,s2,...,sn](s1<s2<...<sn) si=((a-fi)'*Q^-1*(a-fi))
%          Z   = LAMBDA Z-transformation matrix
%
% version: $Revision: 2 $ $Date: 06/01/29 7:15 $
% history: 2006/01/26 1.1 new

if isempty(a)|isempty(Q), f=[]; s=[]; Z=[]; return, end
[L,D]=LD(Q);               % LDL' decomposition
[Z,L,D]=reduction(L,D);    % LAMBDA decorrelation/reduction
[z,s]=search(L,D,Z'*a,n);  % Modified LAMBDA search
f=Z'\z;

% LDL' decomposition : Q=L*diag(D)*L' ------------------------------------------
function [L,D]=LD(Q)
for i=size(Q,1):-1:1
   D(i,1)=Q(i,i);
   L(i,1:i)=Q(i,1:i)/sqrt(Q(i,i));
   for j=1:i-1, Q(j,1:j)=Q(j,1:j)-L(i,1:j)*L(i,j); end
   L(i,1:i)=L(i,1:i)/L(i,i);
end

% LAMBDA decorrelation/reduction : Z'*Q*Z=L'*diag(D)*L -------------------------
function [Z,L,D]=reduction(L,D)
n=size(L,1); Z=eye(n); j=n-1; k=n-1;
while j>0
   if j<=k
       for i=j+1:n
           mu=round(L(i,j));
           if mu~=0
               L(i:end,j)=L(i:end,j)-mu*L(i:end,i);
               Z(:,j)=Z(:,j)-mu*Z(:,i);
           end
       end
   end
   del=D(j)+L(j+1,j)^2*D(j+1);
   if del<D(j+1)
       eta=D(j)/del;
       lam=D(j+1)*L(j+1,j)/del;
       D(j:j+1)=[eta*D(j+1);del];
       L(j:j+1,1:j-1)=[-L(j+1,j),1;eta,lam]*L(j:j+1,1:j-1);
       L(j+1,j)=lam;
       L(j+2:end,[j,j+1])=L(j+2:end,[j+1,j]);
       Z(:,[j,j+1])=Z(:,[j+1,j]);
       k=j; j=n-1;
   else
       j=j-1;
   end
end

% Modified LAMBDA search -------------------------------------------------------
function [zn,fn]=search(L,D,zs,p)
zn=[]; fn=[];
maxdist=inf;
n=length(zs); S=zeros(n);
k=n; dist(k)=0;
zb(n)=zs(n);
z(n,1)=round(zb(n)); y=zb(n)-z(n); step(n)=sgn(y); imax=p;
while 1
    newdist=dist(k)+y^2/D(k);
    if newdist<maxdist
        if k~=1
            k=k-1;
            dist(k)=newdist;
            S(k,1:k)=S(k+1,1:k)+(z(k+1)-zb(k+1))*L(k+1,1:k);
            zb(k)=zs(k)+S(k,k);
            z(k)=round(zb(k)); y=zb(k)-z(k); step(k)=sgn(y);
        else
            if length(fn)<p-1
                zn=[zn,z];
                fn=[fn,newdist];
            else
                zn(:,imax)=z;
                fn(imax)=newdist;
                [f,imax]=max(fn);
                maxdist=f;
            end
            z(1)=z(1)+step(k);
            y=zb(1)-z(1);
            step(1)=-step(1)-sgn(step(1));
        end
    else
        if k==n, break
        else
            k=k+1;
            z(k)=z(k)+step(k);
            y=zb(k)-z(k);
            step(k)=-step(k)-sgn(step(k));
        end
    end
end
[fn,i]=sort(fn); zn=zn(:,i);

% sign operator ----------------------------------------------------------------
function s=sgn(x), if x<0, s=-1; else s=1; end
