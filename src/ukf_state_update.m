function [x,P]=ukf_state_update(sats,x,P,rr,obs,iobs,mode)
C=299792458; lam=C/1.57542E9;
F=ones(length(x),1); Q=zeros(length(x),1);
if ~mode & isnan(x(1))
    x(1:3)=rr; 
    F(1:3)=0; 
    Q(1:3)=0.1; 
elseif mode
    F(1:3)=0; 
    Q(1:3)=2; 
end
N=obs_dd(sats,obs,iobs,2)-obs_dd(sats,obs,iobs,1)/lam;
% i=find((isnan(x(4:end))|~isnan(slip(:,1))|~isnan(slip(:,2)))&~isnan(N)); 
i=find((isnan(x(4:end)))&~isnan(N));
x(3+i)=N(i); 
F(3+i)=0; 
Q(3+i)=10;
% if mode
%       F(4:end)=0;
%     for j = 4:length(x) 
%       Q(j)=0.1;
%     end
% end
P=diag(F)*P*diag(F)'+diag(Q).^2;
P   %pretect
mode
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