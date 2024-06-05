clear all;clc
fidx=dir('*.abf');

StimAmp_DP=0:50:900;
TimeStim=[0.22,0.83];

%%
for k=1:length(fidx)
    filename=fidx(k).name;
    [data,si,~]=abfload(filename,'start',0,'stop','e','channels','a');
    %    [data,si,~]=abfload(filename);
    dt=si*1E-6;
    dataV=data(:,1,:);dataV=squeeze(dataV);%第三个是SWEEP数
    tspanV=linspace(0,size(dataV,1)*dt,size(dataV,1));
    
    %% count peak num
    waveformV=[];
    StimAmp=[0:50:900];
    spkFreq=[];
    Vrest=[];
    
%     StimAmp_H=[];
%     Rin=[];
%     tau=[];
%     sagRatio=[];
%     
%     Fun=@(tau,A,B,x) A.*exp(-x./tau)+B;
%     fo= fitoptions('Method','NonlinearLeastSquares',...
%         'Lower',[0,1,-100],...
%         'Upper',[10,100,10], ...
%         'Startpoint',[0.1,15,-70]);
%     
    for i=1:size(dataV,2)
%         for k=1:size(StimAmp)
            % ========================================= 去极化刺激的分析
            peakidx=peakfinder(dataV(:,i), -20, 10, 1, 0); %(x0, sel, thresh, extrema, include_endpoints)%i表示每一行
            spkTime{i}=peakidx.*dt;
            spkFreq=[spkFreq;numel(peakidx)/range(TimeStim)];
            waveformV=[waveformV,dataV(:,i)];
            Vrest=[Vrest;median(dataV(round(TimeStim(1)/dt-0.2/dt:TimeStim(1)/dt),i))];
            
            figure(1),clf
            plot(tspanV,dataV(:,i)),hold on
            plot(tspanV(peakidx),dataV(peakidx,i),'ro')
            ylim([-130,100])
            title(['I = ',num2str(StimAmp(i)),'   Freq = ',num2str(spkFreq(i))])
            drawnow
%         end               
    end
    
    %% 保存结果
    
    save(['DATA_',filename(1:end-4),'.mat'],'waveformV','tspanV',...
        'StimAmp','spkFreq','spkTime','TimeStim',...
        'Vrest');
    
end