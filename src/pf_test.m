function [x,pfx] = pf_test(x,rr,sref,sats,obsk,iobsk,mode)
C=299792458; 
lam=C/1.57542E9;
numSamples=100;
% rl=0.5;
% rtheta=0.02;
N=obs_dd(sats,obsk,iobsk,2)-obs_dd(sats,obsk,iobsk,1)/lam;
% N =InitfloatAmb(rr,iobsk(:,2))';
[z,R]=pf_obs_dds(sref,sats,obsk,iobsk,2);
% F=ones(length(N)+3,1); 
Q=zeros(length(N)+3,1);
for j = 1:numSamples
   if mode|isnan(x(1,j))
      Q(1:3)=1; 
%       F(1:3)=0; 
      x(1:3,j)=rr; 
   end    
   i=find((isnan(x(4:end,j)))&~isnan(N));
   Q(3+i)=100;
%    F(3+i)=0; 
   x(3+i,j)=N(i);
   for j = 4:length(x)
     k = find(i+3,j);
     if isempty(k)
         Q(k)=0.001;
     end   
   end
   x(:,j) = x(:,j)+sqrt(diag(Q).^2)*randn(length(x(:,j)),1);
   h = pf_measurement_update(x(:,j),sref,sats);
   pfdeta = z - h;
   pfq(j)=exp(-0.5*pfdeta'*inv(R)*pfdeta/10000000);   
end
   pfqsum = sum(pfq);
   pfq = pfq/pfqsum;
   %残差重采样方法
   m = size(x,1);
   pfx=zeros(m,1);
   for i=1:numSamples
      pfx=pfx+pfq(i)*x(:,i);
   end   
   [pfq,II] = sort(pfq);
    j = 1;
    i = numSamples; 
    while (i > 0)
         mmm = ceil(pfq(i)*numSamples);
         while (mmm > 0)
             pfx(:,j) = x(:,II(i));%+randn(1,ss)*noisebaseQ;
             j = j+1;
             mmm = mmm -1;
             if j>numSamples
                 break;
             end
          end
          if j>numSamples 
              break;
          end
          i = i -1;
    end
     x = pfx;    
end
function y=obs_dd(sats,obs,iobs,ch)
for n=1:length(sats)
    y(n,1)=obsdat(sats(n),obs,iobs,ch);
end
end
function y=obsdat(sat,obs,iobs,ch)
i=find(iobs(:,2)==sat);
if ~isempty(i)
    y=obs(i(1),ch); 
else
    y=nan; 
end
end
function h = pf_measurement_update(X,sref,sats)
C=299792458; lam=C/1.57542E9;
r11=geodist(X(1:3),sref);
i = find(sats,sref);
j = 1;
for n=1:length(sats)
    if n == i
        continue;
    end
    r12=geodist(X(1:3),sats(n));
    h(j,1)=(r11-r12)+lam*(X(sref+3)-X(sats(n)+3)); 
    j = j+1;
end
end
% geometric distance -----------------------------------------------------------
%Input: t0-date;t-epoch(obstime - recvr clkerr);rr-接收机坐标,nav
%Output:r-卫地距;e-观测方程系数阵，el-卫星高度角，dt-卫星钟差
function r=geodist(rr,sat)
% C=299792458; OMGE=7.292115167E-5;
sat_pos = [1.4118 2.8093 11.344
           1.9795 4.3758 11.3412
           1.1906 5.8387 11.3368
           -0.1482 6.2278 11.3593
           -1.4198 5.6642 11.3349
           -1.9875 4.1299 11.4931
           -1.239 2.6598 11.343
           0.1128 2.2419 11.3654];   
r=0;
% rk=1;
rs= sat_pos(sat,:);
rrs=rr-rs';
r=norm(rrs);
end
function [y,R]=pf_obs_dds(sref,sats,obs,iobs,ch)
    y1=pf_obsdat(sref,obs,iobs,ch);
    i = find(sats,sref);
    j=1;
for n=1:length(sats)
    C=299792458; lam=C/1.57542E9;
    y2=pf_obsdat(sats(n),obs,iobs,ch);
    if n == i
        continue;
    end    
    y(j,1)=lam*y1-lam*y2;
    j=j+1;
end
R=(ones(length(y(:,1)))+eye(length(y(:,1))))*0.01^2;
end
function y=pf_obsdat(sat,obs,iobs,ch)
i=find(iobs(:,2)==sat);
if ~isempty(i)
    y=obs(i(1),ch); 
else
    y=nan; 
end
end
function N =InitfloatAmb(rr,sats)
C=299792458; lam=C/1.57542E9;
for i = 1:length(sats)
    r12=geodist(rr,sats(i));
    N(i)=r12/lam;
end
end