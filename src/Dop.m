figure;   %  Estimate Bias of Reference Point
clear;
clc;
% satpos = [1.4118 2.8093 11.344
%            1.9795 4.3758 11.3412
%            1.1906 5.8387 11.3368
%            -0.1482 6.2278 11.3593
%            -1.4198 5.6642 11.3349
%            -1.9875 4.1299 11.4931
%            -1.2390 2.6598 11.3430
%            0.1128 2.2419 11.3654];
sat_pos = [1.4118	2.8093	11.344
1.9795	4.3758	11.3412
1.1906	5.8387	11.3368
-0.1482	6.2278	11.3593
-1.4198	5.6642	11.3349
-1.9875	4.1299	11.4931
-1.239	2.6598	11.343
0.1128	2.2419	11.3654];
       [m,n] = size(sat_pos);
for j = 1:n
    center_xyz(j) = mean(sat_pos(:,j));
end
M_Dist_all=[];
for x =  -15:2:15
  for y =  -15:2:15  
      for z =  -15:2:15
       M_Dist_all=[M_Dist_all; x y z];
      end
  end
end
% M_Dist_all(:,1) = -20:2:20;
% M_Dist_all(:,2) = -20:2:20;
% M_Dist_all(:,3) = -60:2:-20;
DtX=M_Dist_all(:,1)-center_xyz(1);
DtY=M_Dist_all(:,2)-center_xyz(2);
dop =prmodel([M_Dist_all(:,1) M_Dist_all(:,2) M_Dist_all(:,3)]',1:8);
% dop = dop';
[X,Y]=meshgrid(min(DtX):0.05:max(DtX),min(DtY):0.05:max(DtY));%将x、y轴网格化，重构用于画图x、y轴数据
Z=griddata(DtX,DtY,dop,X,Y);%插值，重构用于画图的Z轴数据
mesh(X,Y,Z,'FaceLighting','gouraud','LineWidth',0.3);
% mesh(X,Y,Z,'FaceLighting','gouraud','LineWidth',0.3);
  grid on;
  xlabel('X/m');ylabel('Y/m');zlabel('PDOP');
  xlim([-15 15]);ylim([-15 15]);set(gca,'XTick',-15:2:15);
  set(gca,'YTick',-15:2:15);
%  set(gca,'ZTick',0:0.05:0.5);
  colorbar;
set(gca,'Fontname', 'Palatino Linotype','FontSize',12);
hold on 
% [m,n]=size(satpos);
scatter3(sat_pos(:,1),sat_pos(:,2),sat_pos(:,3));
% for j = 1:m
%     [x y z]=sphere();
% %     surf(x/4+satpos(j,1),y/4+satpos(j,2),z/4+satpos(j,3));%绘制半径为2的球
%     plot3(satpos(j,1),satpos(j,2),satpos(j,3));
% %     scatter3(satpos(j,1),satpos(j,2),satpos(j,3),1,'filled');
% end
function dop =prmodel(rr,sats)
[i,j]=size(rr);
for k=1:j
    for n=1:length(sats)
    [r,e]=geodist(rr(:,k),sats(n));
    H(n,:)=[e(1:2)',1];%加入卫星钟差和对流层改正
    end
    HTH = inv(H'*H);
    dop(k)=sqrt(trace(HTH(1:2,1:2)));
end
% for k=1:j
%     for n=1:length(sats)
%     [r,e]=geodist(rr(:,k),sats(n));
%     H(n,:)=[e',1];%加入卫星钟差和对流层改正
% end
%     HTH = inv(H'*H);
%     dop(1,k)=sqrt(trace(HTH));
%     dop(2,k)=sqrt(trace(HTH(1:3,1:3)));
%     dop(3,k)=sqrt(trace(HTH(1:2,1:2)));
% end
end
function [r,e]=geodist(rr,sat)
% sat_pos = [1.4118 2.8093 11.344
%            1.9795 4.3758 11.3412
%            1.1906 5.8387 11.3368
%            -0.1482 6.2278 11.3593
%            -1.4198 5.6642 11.3349
%            -1.9875 4.1299 11.4931
%            -1.239  2.6598 11.343
%            0.1128 2.2419 11.3654];   

sat_pos = [1.4118	2.8093	11.344
1.9795	4.3758	11.3412
1.1906	5.8387	11.3368
-0.1482	6.2278	11.3593
-1.4198	5.6642	11.3349
-1.9875	4.1299	11.4931
-1.239	2.6598	11.343
0.1128	2.2419	11.3654];
r=0;
rs= sat_pos(sat,:);
rrs=rr-rs';
r=norm(rrs);
e=rrs/r;
end