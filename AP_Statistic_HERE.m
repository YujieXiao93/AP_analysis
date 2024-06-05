function [spkwaveform,ttspan,...
    HW, AHP, Threshold, Amp,max_slope,min_slope,...
    HWidx,AHPidx, Thresholdidx,max_slopeidx,min_slopeidx] = AP_Statistic_HERE(data,peaktime,tspan,Goodidx)
%INPUT
% data :   a vector
% tspan:  a vector,have the same elememt with data
% NumSpk: successive number of  spike to be analized
%%
dt=tspan(end)-tspan(end-1);

spktime=peaktime(Goodidx);
%
pretime=0.04;
postime=0.08;

spkwaveform=[];

for i=1:length(spktime)
    idx=round(spktime(i)/dt-pretime/dt:spktime(i)/dt+postime/dt);
    ttspan=linspace(-pretime,postime,numel(idx));
    spkwaveform=[spkwaveform,data(idx)];
end

% figure,clf
% plot(spkwaveform)%确定提取出来的波形是正确的
%%

for i=1:length(spktime)
    data=spkwaveform(:,i);
    slope=diff(data)./dt./1000;
    slope=[0;slope];
    
    if i<length(peaktime)-1
        endtime=(peaktime(i+1)-peaktime(i));
    else
        endtime=postime*0.85;
    end
    
    max_slopeidx(i)=find(slope(ttspan<endtime)==max(slope(ttspan<endtime)),1,'first');
    max_slope(i)=slope(max_slopeidx(i));
    min_slopeidx(i)=find(slope(ttspan<endtime)==min(slope(ttspan<endtime)),1,'first');
    min_slope(i)=slope(min_slopeidx(i));
    
%     figure,plot(slope)
    idxA=find(slope(:)>20);
    idxB=find(data(:)>-60);
    idxC=intersect(idxA(:),idxB(:));
    Thresholdidx(i)=idxC(1);
    Threshold(i)=data(Thresholdidx(i));
    AHPidx(i)=find(data(ttspan>0&ttspan<endtime)==min(data(ttspan>0&ttspan<endtime)),1,'first')+sum(ttspan<0);
    AHP(i)=data(Thresholdidx(i))-data(AHPidx(i));
    
    Amp(i)=max(data)-data(Thresholdidx(i));
    
    %     HWidx(i,:) = eventedgefinder_HERE(data,data(Thresholdidx(i))+Amp(i)/2,1/dt,0.0001,0.0001,1,1);
    HWidx(i,:) = eventedgefinder_HERE(data,data(Thresholdidx(i))+Amp(i)/2,1,1);
    HW(i)=(HWidx(i,2)-HWidx(i,1))*dt*1000;
    
end

end


function d=deriv(a)
% First derivative of vector using 2-point central difference.
n=length(a);
d(1)=a(2)-a(1);
d(n)=a(n)-a(n-1);
for j = 2:n-1;
    d(j)=(a(j+1)-a(j-1)) ./ 2;
end
end

function d=deriv2(a)
n=length(a);
for j = 2:n-1;
    d(j)=a(j+1) - 2.*a(j) + a(j-1);
end
d(1)=d(2);
d(n)=d(n-1);
end

function d3=deriv3(a)
% Third derivative of vector a
%  T. C. O'Haver, 2008.
n=length(a);
for j = 3:n-2;
    d3(j)=a(j+2) - 2.*a(j+1) + 2.*a(j-1) - a(j-2);
end
d3(1:2)=d3(3);
d3(n-1:n)=d3(n-2);
end