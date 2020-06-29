function [obs,iobs,time]=readrnx(files,rcv)
%     type=0;trnx='';tobs=[];ts=0;te=0;t0=0;
    f=fopen(files,'rt');
    if f<0
        error(['file open error :',files]); 
    end
    disp(['reading rinex file ... : ',files]);
    [trnx,tobs,t0,ts,te,type]=readrnxh(f);
    if trnx=='O'
        [obs,iobs,time]=readrnxo(t0,ts,te,f,tobs,type,rcv);
%   -------------------------------------     
    end
    fclose(f);

% read rinex header -----------------------------------------------------------
function [trnx,tobs,t0,ts,te,type]=readrnxh(f)
trnx=''; tobs={};te=[];
while 1
    s=fgetl(f); if ~isstr(s), break, end
    label=subs(s,61,20);
    if findstr(label,'RINEX VERSION / TYPE')
        trnx=subs(s,21,1);
        type = str2num(subs(s,6,4));
    elseif findstr(label,'# / TYPES OF OBSERV')
%     elseif findstr(label,'SYS / # / OBS TYPES  ')
%     elseif findstr(label,'# / TYPES OF OBSERV')
%     elseif findstr(label,'SYS / # / OBS TYPES ')
       
    elseif findstr(label,'# / TYPES OF OBSERV')
        if(type < 3)
          p=11;
          for i=1:s2n(s,1,6)
            if p>=59, s=fgetl(f); p=11; end
              tobs{i}=subs(s,p,2); p=p+6;
          end
        else
          p=8;
          if strcmp(s(1),'G')
           for i=1:s2n(s,2,5)
            if p>=59
                s=fgetl(f); 
                p=8; 
            end
              tobs{i}=subs(s,p,2); p=p+4;
          end
          end
        end
    elseif findstr(label,'SYS / # / OBS TYPES')
        if(type < 3)
          p=11;
          for i=1:s2n(s,1,6)
            if p>=59, s=fgetl(f); p=11; end
              tobs{i}=subs(s,p,2); p=p+6;
          end
        else
          p=8;
          if strcmp(s(1),'G')
           for i=1:s2n(s,2,5)
            if p>=59
                s=fgetl(f); 
                p=8; 
            end
              tobs{i}=subs(s,p,2); p=p+4;
          end
          end
        end
    elseif findstr(label,'TIME OF FIRST OBS')
        year =subs(s,3,4);
        mon  =subs(s,11,2);
        day  =subs(s,17,2);
        t0=datenum(str2num(year),str2num(mon),str2num(day));
        hour =subs(s,23,3);
        min  =subs(s,29,2);
        sec  =subs(s,34,10);
        ts   =str2num(hour)*3600+str2num(min)*60+str2num(sec);
    elseif findstr(label,'TIME OF LAST OBS')
        year =subs(s,3,4);
        mon  =subs(s,11,2);
        day  =subs(s,17,2);
        hour =subs(s,23,3);
        min  =subs(s,29,2);
        sec  =subs(s,34,10);
        t0=datenum(str2num(year),str2num(mon),str2num(day));
        te   =str2num(hour)*3600+str2num(min)*60+str2num(sec);
    elseif findstr(label,'END OF HEADER')
        break;
    end
end
% read rinex observation data --------------------------------------------------
function [obs,iobs,time]=readrnxo(t0,ts,te,f,tobs,type,rcv)
to={'C1','L1','D1','S1'}; ind=zeros(1,length(tobs));time = [];
for n=1:4
    i=find(strcmp(to{n},tobs)); 
    if ~isempty(i)
        ind(i)=n; 
    end
end
obs=repmat(nan,400000,4); iobs=zeros(400000,3); n=0;
while 1
    s=fgetl(f); if ~isstr(s), break, end
    if type<3
        e=s2e(s,1,26);
    if length(e)>=6
        if e(1)<80
            y=2000; 
        else
            y=1900; 
        end
        tt=(datenum(e(1)+y,e(2),e(3))-t0)*86400+e(4:6)*[3600;60;1];
        if round(tt)>te
            break;
%         else
%             flag=round(tt)>=ts;         
        end
        flag=1;
        ns=s2n(s,30,3); 
        iobs(n+(1:ns),1)=tt; 
        if rcv == 1
             time = [time tt];
        end      
        p=33;
        for i=1:ns
            if p>=69
                s=fgetl(f); 
                p=33; 
            end
            if flag&any(subs(s,p,1)==' G')
                iobs(n+i,2)=s2n(s,p+1,2); 
            end
            p=p+3;
        end
        for i=1:ns
            s=fgetl(f); p=1;
            for j=1:length(tobs)
                if p>=80, s=fgetl(f); p=1; end
                if flag&ind(j)>0, obs(n+i,ind(j))=s2n(s,p,14); end, p=p+16;
            end
        end
        if flag
            n=n+ns; 
        end
    end
    else
    e=s2e(s,2,28);
    if length(e)>=6
        tt=(datenum(e(1),e(2),e(3))-t0)*86400+e(4:6)*[3600;60;1];
        if round(tt)>te, break, else flag=round(tt)>=ts; end
        ns=s2n(s,33,3); 
        iobs(n+(1:ns),1)=tt; 
        if rcv == 1
             time = [time tt];
        end       
        for i=1:ns
            s=fgetl(f); p=1;
           if flag&any(subs(s,p,1)=='G')
               iobs(n+i,2)=s2n(s,p+1,2); 
               p = p+3;
           end
            if flag&any(s(1)=='G')
            for j=1:length(tobs)
                if p>=67, s=fgetl(f); p=1; end
                if flag&ind(j)>0
                    obs(n+i,ind(j))=s2n(s,p,14); 
                end
                p=p+16;
            end
            end
        end
        if flag, n=n+ns; end
    end
    end
end
obs=obs(1:n,:); obs(obs==0)=nan; 
% iobs(iobs(:,2)==0)=nan;
iobs=iobs(1:n,:); 
iobs(1:n,3)=rcv;
% j=iobs(iobs(:,2)==0);
% iobs(j,2) = nan;
% % j =  find(iobs(:,2)~=0);
% i = find(~isnan(iobs(:,2)));
% obs=obs(i,:);iobs=iobs(i,:);

% string to number/epoch/substring  --------------------------------------------
function a=s2n(s,p,n), a=sscanf(subs(s,p,n),'%f'); if isempty(a), a=nan; end
function e=s2e(s,p,n), e=sscanf(subs(s,p,n),'%d %d %d %d %d %f',6)';
function s=subs(s,p,n), s=s(p:min(p+n-1,length(s)));
