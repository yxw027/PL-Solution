function [xs,xres] =ppp_SD_LS(time,obs,iobs,mode,type)
% fid=fopen('output_debug.txt', 'wt');;
MAXNUM= 8;
NB = [];
NL = [];
C=299792458; 
lam=C/1.57542E9; 
if type == 1
    numSamples=100;
    for i= 1:numSamples
        x(1:3+MAXNUM,i)=nan; 
    end   
else
    x(1:3+MAXNUM,1)=nan; 
    P=zeros(length(x),length(x));  
end     
for k=1:length(time)
    if k>6000
        break;
    end
    if isempty(find(10*round(iobs(:,1))==10*round(time(k))))
        continue;
    end
     sref=1; sats=1:8; 
     i=find(round(10*iobs(:,1))==round(10*time(k)));     
     j = find(iobs(i,2)<=8&iobs(i,2)>0);
     iobskk=iobs(i,:);obskk=obs(i,:);
     obsk=obskk(j,:); iobsk=iobskk(j,:);  
     rr = zeros(3,1);
    if length(j)<4
        continue;
    end
%     [rr,sat,dop]=pointp_SD(obsk,iobsk);
%     kk = find(~isnan(rr));
%     if kk>0   
%         continue;
%     end
    slip=zeros(MAXNUM,1);
%     if k==1, sref=sat; 
%         sats=1:31; 
%         sats(sat)=[]; 
%     end
%     if k>1 
%          i_pre=find(round(iobs(:,1))==round(time(k-1)));
%           j_pre = find(iobs(i_pre,2)<=8);
%         iobs_pre = iobs(j_pre,:);obs_pre=obs(j_pre,:);
%         slip = detslp_dop(obs,iobs,obs_pre,iobs_pre,1,lam);
%     end
    % temporal update of states
%     rr=[2.6777 6.059 -1.7463]';true value
if k<60
    mode = 0;
else
    mode = 1;
end
       y=obs_dds(sref,sats,obsk,iobsk,2);%L1
       y_code=obs_dds(sref,sats,obsk,iobsk,1);%L1
%        rr=[0.025 3.317 -0.530]';
        rr=[-2.9701 5.7737 -0.3562]';%GPS-2 123
%          rr=[3.0223 0.7152 -0.3879]';%GPS-3-7         
%          rr = [2.6777 6.0590 -0.3719]'%GPS-6
      if type == 0
          [x,P] = ukf_test(x,P,rr,sref,sats,obsk,iobsk,mode);
      elseif type == 1
          [x,pfx] = pf_test(x,rr,sref,sats,obsk,iobsk,mode);
      else
           x = udstate(sats,x,rr,obsk,iobsk,mode,lam);              
             % observables/measurement model                       
            if ~mode
                while 1
                [h,H,R]=measmodel(x(1:3),x(4:end),sref,sats,lam,mode);
               i=find(~isnan(y)&~isnan(h)); 
               H=H(i,:);
%                if isempty(NB)
%                  NB = H'*R(i,i)*H;
%                  NL = H'*R(i,i)*(y(i)-h(i)); 
%                else
%                  NB = H'*R(i,i)*H+NB;
%                  NL = H'*R(i,i)*(y(i)-h(i))+NL; 
%                end
               NB = H'*R(i,i)*H;
               NL = H'*R(i,i)*(y(i)-h(i)); 
               x1 = NB\NL;
               P(4:end,4:end)=inv(NB);
               x(4:end) = x1+x(4:end);              
                if norm(x1)<0.002
                  break;
                 end
                end
           else
              r_res = Optimalsearch(x,sref,sats,lam,y);
               x(1:2) = r_res(1,1:2);
               r_res(1,1)
%              i=find(~isnan(y)&~isnan(h)); 
%              H=H(i,:);
%              HTPH = H'*R(i,i)*H;
%              HTPL = H'*R(i,i)*(y(i)-h(i)); 
%              x1 = HTPH\HTPL;
%              P(1:3,1:3) = inv(HTPH); 
%              x(1:3) =  x1+x(1:3);
%               x(3) = rr(3);
%              if norm(x1(1:3))<0.01
%                 break;
          end     
%           while 1
%              x = udstate(sats,x,rr,obsk,iobsk,mode,lam);    
%              % observables/measurement model          
%              [h,H,R]=measmodel(x(1:3),x(4:end),sref,sats,lam,mode);  
%              NB = H'*R*H+NB;
%              NL = H'*R*(y-h)+NL;
%              % measurment update of states
%              [x1,P1]=filt(x,y,h,H,R); 
%              t = x1-x;
%              if norm(t(1:3))<2
%                  x = x1;
%                  P = P1;
%                  break;
%              end
%           end              
      end      
    % [xa,Pa]=ddmat(x,P);
    % ambiguity resolution
    % xf(k,:)=fixamb(xa,Pa)';     
    if(type == 1)
       res = measmodel_res(pfx(1:3),pfx(4:end),sref,sats,lam,y);
       xs(k,:)=pfx(1:3)';   
       xres(k,:)=res;
    else
       res = measmodel_res(x(1:3),x(4:end),sref,sats,lam,y);
       xs(k,:)=x(1:3)';   
       xres(k,:)=res;
    end
       Na = Na_dds(sref,sats,x);%L1
       phase_res(k,:)=y';
       code_res(k,:)=y_code';
       N_res(k,:) = Na';       
end
    tt = 1:length(xs(:,1));
    figure(1);
    subplot(4,1,1)    
    plot(tt,xs(:,1),'LineWidth',3)
    title('xs');
    subplot(4,1,2)
    plot(tt,xs(:,2),'LineWidth',3)
    subplot(4,1,3)
    plot(tt,xs(:,3),'LineWidth',3)
    subplot(4,1,4)
    plot(xs(:,1),xs(:,2),'LineWidth',3)
    ylabel('PlaneTrack(m)','fontsize',10);
    figure(2)
    subplot(2,1,1) 
    [m,n]=size(phase_res);
     plot(1:m,phase_res(:,1),'LineWidth',3);
     hold on;
    for i = 2:n
        plot(1:m,phase_res(:,i),'LineWidth',3);
    end
     title('phase_res');
     ylabel('phase-res (m)','fontsize',20);
     subplot(2,1,2)
     [m,n]=size(code_res);
     plot(1:m,code_res(:,1),'LineWidth',3);
     hold on;
     for i = 2:n
        plot(1:m,code_res(:,i),'LineWidth',3);
     end
     ylabel('code-res (m)','fontsize',20);
     figure(3);
     [m,n]=size(N_res);
     plot(1:m, N_res,'LineWidth',3);
     ylabel('amb-res (m)','fontsize',20);
end
function Na = Na_dds(sref,sats,x)
i = find(sats,sref);k = 1;
for j = 1:length(sats)
    if i == j
        continue;
    end
    Na(k)=x(3+i)-x(3+j);
    k=k+1;
end
end
% single point positioning -----------------------------------------------------
%输入：t0-日期，obs-观测值，iobs-观测时刻，nav-电文，inav-PRN
%输出：rr-坐标（ref&rov），t-接收时刻-接收机钟差(ref&rov)，sat-参考星rov，dop-DOP
function [rr,sat,dop]=pointp_SD(obs,iobs)
C=299792458; f1=1.57542E9; f2=1.2276E9;
    i=find(iobs(:,3)==1); 
    sats=iobs(i,2);
    y=obs(i,1); % ion-free pseudorange
    x=zeros(3,1); xk=ones(3,1);
    while norm(x-xk)>0.1
        HH = [];H=[];h=[];yy=[];
        [h,H,el]=prmodel_SD(x(1:3),sats);
        i=find(~isnan(y)&~isnan(h)); 
        if length(i)<4
            x(:)=nan; 
            break;
        end
        HH=H(2:end,:)-H(1,:);
        hh = h(2:end,:)-h(1,:);
        yy = y(2:end,:)-y(1,:);    
        xk = x;
        x=x+(HH'*HH)\HH'*(yy-hh);
    end
    rr=x;
    [e,i]=max(el); sat=sats(i); 
    dop=sqrt(trace(inv(H'*H))); 
end
function [rr,sat,dop]=pointp_SD_cp(obs,iobs)
C=299792458; f1=1.57542E9;
lam = C/f1;
    i=find(iobs(:,3)==1);sats=iobs(i,2);
    yc=obs(i,1); % ion-free pseudorange
    yp=obs(i,2); 
    n = length(yp);
    x=zeros(3+n,1); xk=ones(3+n,1);
    while norm(x-xk)>0.1
        [h,H,el]=prmodel_SD(x(1:3),sats);
        i=find(~isnan(yc)&~isnan(h)&~isnan(yp)); 
         if length(i)<4
            x(:)=nan; 
            break;
         end
         HH=[];
        for j=2:n
            a = H(1,:)'-H(j,:)';
            b = [sats(1)+3 sats(j)+3];
            HH(j-1,[1:3,b])=[a',[lam,-lam]];
        end
        Hc = zeros(n-1,3+n);
        Hc(:,1:3)=H(2:end,:)-H(1,:);
        HH = [HH;Hc];
        yy1 = lam*yp(1)-lam*yp(2:end);
        yy2 = yc(1)-yc(2:end);
        hh = h(1)-h(2:end);
        hh = [hh;hh];
        yy = [yy1;yy2];
        x=x+(HH'*HH)\HH'*(yy-hh);
    end
    rr=x(1:3);
    [e,i]=max(el); sat=sats(i); 
    dop=sqrt(trace(inv(H'*H))); 
end
function [rr,t,sat,dop]=pointp(obs,iobs)
C=299792458; f1=1.57542E9; f2=1.2276E9;
    i=find(iobs(:,3)==1); 
    tr=iobs(i(1),1); 
    sats=iobs(i,2);
    y=obs(i,1); % ion-free pseudorange
    x=zeros(4,1); 
    xk=ones(4,1);
    while norm(x-xk)>0.1
        [h,H,el]=prmodel(x(1:3),sats);
        i=find(~isnan(y)&~isnan(h)); 
        H=H(i,:);
        if length(i)<4
            x(:)=nan; 
            break;
        end
        xk=x; x=x+(H'*H)\H'*(y(i)-h(i));
    end
    rr=x(1:3); 
    t=tr-x(4)/C; % t=tr-dtr
    [e,i]=max(el);
    sat=sats(i); 
    dop=sqrt(trace(inv(H'*H)));
end
% pseudorange model (ionosphere-free) ------------------------------------------
function [h,H,el]=prmodel_SD(rr,sats)
for n=1:length(sats)
    [r,e,el(n)]=geodist(rr,sats(n));
    h(n,1)=r; 
    H(n,:)=e';%加入卫星钟差和对流层改正
end
end

function [h,H,el]=prmodel(rr,sats)
for n=1:length(sats)
    [r,e,el(n)]=geodist(rr,sats(n));
    h(n,1)=r; 
%     h(n,1)=r+2.4/sin(el(n)); 
    H(n,:)=[e',1];%加入卫星钟差和对流层改正
end
end
function N =InitfloatAmb(rr,sats)
C=299792458; lam=C/1.57542E9;
for i = 1:length(sats)
    [r12,e12]=geodist(rr,sats(i));
    N(i)=r12/lam;
end
end
% temporal update of states ----------------------------------------------------
function x=udstate(sats,x,rr,obs,iobs,mode,lam)
if ~mode & isnan(x(1))
    x(1:3)=rr; 
end
% x(3)=-0.327;
%  N =InitfloatAmb(x(1:3),sats)';
N=obs_dd(sats,obs,iobs,2)-obs_dd(sats,obs,iobs,1)/lam;
% i=find((isnan(x(4:end))|~isnan(slip(:,1))|~isnan(slip(:,2)))&~isnan(N)); 
i=find((isnan(x(4:end)))&~isnan(N)); 
x(3+i)=N(i); 
end
% double difference of observables ---------------------------------------------
%计算双差值
function y=obs_dds(sref,sats,obs,iobs,ch)
C=299792458; lam=C/1.57542E9;
y = zeros(8,1);
    y1=obsdat(sref,obs,iobs,ch);
    i = find(sats,sref);
for n=1:length(sats)
    y2=obsdat(sats(n),obs,iobs,ch);
    if n == i
        continue;
    end   
    y(sats(n),1)=lam*(y1-y2);
end
y(find(y==0))=nan;
end
function y=obs_dd(sats,obs,iobs,ch)
% i = iobs(:,2);
% for n=1:i
%     y(n,1) = obs(n,ch);
% end
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
function res = measmodel_res(rr1,N,sref,sats,lam,obs)
[r11,e11]=geodist(rr1,sref);
i = find(sats,sref);
j = 1;
H = [];
for n=1:length(sats)
    if n == i
        continue;
    end
    [r12,e12]=geodist(rr1,sats(n));
    H(sats(n),:)=e11'-e12';
    h(sats(n))=(r11-r12)+lam*(N(sref)-N(sats(n))); 
    res(sats(n))=obs(sats(n)) - h(sats(n));
end    
    res(find(h == 0)) = nan;
end
% single difference of phase model ---------------------------------------------
function [h,H,R]=measmodel(rr1,N,sref,sats,lam,mode)
[r11,e11]=geodist(rr1,sref);
h = zeros(8,1);
i = find(sats,sref);
 for n=1:length(sats)
    if n == i
        continue;
    end
     Na = N(sref)-N(sats(n));
    if mode
        Na = round(N(sref)-N(sats(n)));
    end
    [r12,e12]=geodist(rr1,sats(n));
    if mode
      h(sats(n),1)=(r11-r12)+lam*Na; 
      H(sats(n),:)= e11'-e12';
    else
      h(sats(n),1)=(r11-r12)+lam*Na; 
      H(sats(n),[sref sats(n)])=[lam,-lam];
    end
 end
     R=(ones(length(h))+eye(length(h)))*0.01^2;
     h(find(h==0))=nan;
end

function [xa,Pa]=ddmat(x,P)
i = length(x);
for j = 4:i
    if ~isnan(x(j))
        break
    end
end
k = find(~isnan(x(j+1:end)));
xa = x(k)-x(j);
t = length(k);
dd = zeros(t,i);
dd(1,x(j)) = -1;
for h = 1:t
    dd(h,x(h)) = 1;
end
Pa = dd*P*dd';
end
% ambiguity resolution ---------------------------------------------------------
function x=fixamb(x,P);
i=1:3; j=4:length(x); 
j=j(~isnan(x(j))&diag(P(j,j))<10^2);
[N,s]=mlambda(x(j),P(j,j),2);
if isempty(N)|s(2)/s(1)<3, x(:)=nan; 
    return, end  % ratio-test
x([i,j])=[x(i)-P(i,j)/P(j,j)*(x(j)-N(:,1));N(:,1)];
end
% geometric distance -----------------------------------------------------------
%Input: t0-date;t-epoch(obstime - recvr clkerr);rr-接收机坐标,nav
%Output:r-卫地距;e-观测方程系数阵，el-卫星高度角，dt-卫星钟差
function [r,e,el]=geodist(rr,sat)
% C=299792458; OMGE=7.292115167E-5;
% sat_pos = [1.4118 2.8093 11.344
%            1.9795 4.3758 11.3412
%            1.1906 5.8387 11.3368
%            -0.1482 6.2278 11.3593
%            -1.4198 5.6642 11.3349
%            -1.9875 4.1299 11.4931
%            -1.239 2.6598 11.343
%            0.1128 2.2419 11.3654];   %1-8
% sat_pos = [-2.827 0.3581 11.3545
%            3.1962 0.761	11.3499
%            2.7086 8.1584 11.3408
%           -3.2637 7.7569 11.3347
%            6.459 3.501  6.077
%            -0.256 9.834  6.396
%            -6.443 2.428  6.015
%            0.083 -3.771  6.484];
sat_pos = [1.4118	2.8093	11.344
1.9795	4.3758	11.3412
1.1906	5.8387	11.3368
-0.1482	6.2278	11.3593
-1.4198	5.6642	11.3349
-1.9875	4.1299	11.4931
-1.239	2.6598	11.343
0.1128	2.2419	11.3654];
r=0;
% rk=1;
rs= sat_pos(sat,:);
rrs=rr-rs';
r=norm(rrs);
% while abs(r-rk)>1E-4
%     [rs,dt]=satpos(t0,t-r/C,nav);
%     rrs=rr-Rz(OMGE*r/C)*rs; rk=r; r=norm(rrs);
% end
e=rrs/r;
if norm(rr)>0
    el=-asin(rr'*e/norm(rr)); 
else
    el=pi/2; 
end
end
function slip = detslp_dop(obs,iobs,obs_pre,iobs_pre,rcv,lam)
    DTTOL= 0.005; 
    MAXACC =30.0;
    i=find(iobs(:,3)==rcv); tr=iobs(i(1),1); sats=iobs(i,2);
    i_pre=find(iobs_pre(:,3)==rcv); tr_pre=iobs_pre(i_pre(1),1); sats_pre=iobs_pre(i_pre,2);
    for j = 1:length(sats)
        slip(sats(j),1) = sats(j);
        if (~find(iobs(:,2)==sats(j))||~find(iobs_pre(:,2)==sats(j)))
            continue;
       jj= find(iobs(:,2)==sats(j));   
       jj_pre= find(iobs_pre(:,2)==sats(j));
       dph=obs(jj,3)-obs_pre(jj_pre,3);
       tt=tr-tr_pre;
%  /* cycle slip threshold (cycle) */
        thres=MAXACC*tt*tt/2.0/lam+fabs(tt)*4.0;
       if abs(tt)<DTTOL
           continue;
       end
       dpt=-obs(jj,3)*tt;
        if fabs(dph-dpt)<=thres
            continue;
        else           
            slip(sats(j),2) = 1;
        end     
    end
    slip(slip(:,2)==0)=nan;
    end
end
function res = Optimalsearch(X,sref,sats,lam,obs)
%   n = 1;
%   m = 1;
  rr1 = X(1:3);
  N= X(4:end);
  count = 1;
  for i = X(1)-2:0.1:X(1)+2
      for j = X(2)-2:0.1:X(2)+2
          rr1(1) = i;
          rr1(2) = j;
          ref = measmodel_res(rr1,N,sref,sats,lam,obs);
          k  = find(~isnan(ref));
          LS = norm(ref(k));
%           res(m,n)=LS;
          res_f(count,:)=[LS rr1']';
%           n = n+1;
          count = count+1;
      end
%       m = m+1;
  end
  res_f = sortrows(res_f,1);
  res = res_f(1:2,:);  
end