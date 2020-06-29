function testobs(time,obs,iobs)
    for k=1:length(time)
        if isempty(find(10*round(iobs(:,1))==10*round(time(k))))
          continue;
        end
     sref=1; sats=1:8; 
     i=find(round(10*iobs(:,1))==round(10*time(k)));     
     j = find(iobs(i,2)<=8&iobs(i,2)>0);
     iobskk=iobs(i,:);obskk=obs(i,:);
     obsk=obskk(j,:); iobsk=iobskk(j,:);  
    y=obs_dds(sref,sats,obsk,iobsk,2);%L1
    y_code=obs_dds(sref,sats,obsk,iobsk,1);%L1
    phase_res(k,:)=y';
    code_res(k,:)=y_code';
end
     figure(2);
     subplot(2,1,1);
    [m,n]=size(phase_res);
     plot(1:m,phase_res(:,1),'LineWidth',3);
     hold on;
    for i = 2:n
        plot(1:m,phase_res(:,i),'LineWidth',3);
    end
    title('phase_res');
    ylabel('phase-res (m)','fontsize',20);
     subplot(2,1,2);
     [m,n]=size(code_res);
     plot(1:m,code_res(:,1),'LineWidth',3);
     hold on;
     for i = 2:n
        plot(1:m,code_res(:,i),'LineWidth',3);
     end
     ylabel('code-res (m)','fontsize',20);
     title('code_res');
end
       % double difference of observables ---------------------------------------------
%º∆À„À´≤Ó÷µ
function y=obs_dds(sref,sats,obs,iobs,ch)
y = zeros(8,1);
C=299792458; lam=C/1.57542E9;
    y1=obsdat(sref,obs,iobs,ch);
    i = find(sats==sref);
for n=1:length(sats)
    y2=obsdat(sats(n),obs,iobs,ch);
    if n == i
        continue;
    end
    y(sats(n),1)=(y1-y2);
    if ch == 2
        y(sats(n),1)=lam*y(sats(n),1);
    end   
end
    y(find(y == 0))=nan;
end
function y=obsdat(sat,obs,iobs,ch)
i=find(iobs(:,2)==sat);
if ~isempty(i)
    y=obs(i(1),ch); 
else
    y=nan; 
end
end