function y=obs_dds(sref,sats,obs,iobs,ch)
    y1=obsdat(sref,obs,iobs,ch);
    i = find(sats,sref);
for n=1:length(sats)
    C=299792458; lam=C/1.57542E9;
    y2=obsdat(sats(n),obs,iobs,ch);
    if n == i
        continue;
    end    
    y(sats(n),1)=lam*y1-lam*y2;
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