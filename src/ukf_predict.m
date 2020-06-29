function [X,Pxx,zz,Pzz,Pxz]=ukf_predict(XX,wm,wc,nx,sref,iobsk,mode)%f--linearmean
sats = iobsk(:,2);
X=zeros(nx,1);
lengthx = length(sats)-1;
zz=zeros(lengthx,1);
Pzz=zeros(lengthx,lengthx);
Pxx=zeros(nx,nx);
Pxz=zeros(nx,lengthx);
for k=1:1:2*nx+1
    X=X+wm(k)*XX(:,k);
    z(:,k) = ukf_measurement_update(XX,sref,sats);
    i = find(~isnan(z(:,k)));
    zz=zz+wm(k)*z(i,k);
end
for k=1:1:2*nx+1    
    Pxx=Pxx+wc(k)*((XX(:,k)-X)*(XX(:,k)-X)');
    Pzz=Pzz+wc(k)*((z(i,k)-zz)*(z(i,k)-zz)');
    Pxz=Pxz+wc(k)*((XX(:,k)-X)*(z(i,k)-zz)');
end
R=(ones(length(z(i,1)))+eye(length(z(i,1))))*0.01^2;
% Pxx=Pxx+Q;%是否需要？？
Pzz=R+Pzz;
end
function h = ukf_measurement_update(X,sref,sats)
C=299792458; lam=C/1.57542E9;
h = zeros(8,1);
r11=geodist(X(1:3),sref);
i = find(sats==sref);
for n=1:length(sats)
    if n == i
        continue;
    end
    r12=geodist(X(1:3),sats(n));
    h(sats(n),1)=(r11-r12)+lam*(X(sref+3)-X(sats(n)+3)); 
end
h(find(h==0))=nan;
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
rs= sat_pos(sat,:);
rrs=rr-rs';
r=norm(rrs);
end