clear all;clc

Fidx{1}=dir('DATA_C*.mat');
Fidx{2}=dir('DATA_T*.mat');
GroupName={'C','T'};

Num=7;
NPlusMinus=4;
thAP=1;

Binwidth=50;
BinL=0:Binwidth:900-Binwidth;
BinR=BinL+Binwidth;
BinsFIcurve{1}=[BinL(:),BinR(:)];%
BinsFIcurve{2}=[BinL(:),BinR(:)];%

colormap=[0.5,0.54,0.53;0.12,0.56,1];
%深灰色[0.5,0.54,0.53];红色[0.69,0.19,0.12];橙色[1,0.38,0];浅蓝色[0.12,0.56,1];
%%

for kk=1:length(Fidx)
    for k=1:length(Fidx{kk})
        filename=Fidx{kk}(k).name;
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        dt=tspanV(2)-tspanV(1);
        FilenameRecord{kk}{k,1}=filename;
        
        %% fitting ReLU
        Fun=@(a,b,x) max([zeros(length(x),1),a.*x(:)+b],[],2);
        fo= fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0.001,-1000],...
            'Upper',[10,1000], ...
            'Startpoint',[0.1,-10]);
%         
%         uniqueStim=unique(StimAmp);
%         StimAmp_act=zeros(size(uniqueStim))';
%         spkFreq_act=zeros(size(uniqueStim))';
        for i=1:length(StimAmp)
%             StimAmp_act(i)=StimAmp(i);
%             spkFreq_act(i)=spkFreq(i);
        end
        
        if numel(StimAmp)>5
            [tempx,i]=sort(StimAmp);
            tempy=spkFreq(i);
            maxidx=find(tempy==max(tempy),1,'first');
            if maxidx<numel(tempy)-1&tempy(maxidx)>20
                maxFreq=tempy(maxidx);
                maxStim=tempx(maxidx);
            else
                maxFreq=nan;
                maxStim=nan;
            end
            
            idx=find(tempy>=12,1,'first');%大于12mV视为AP
            x=[0;tempx(1:idx)'];
            y=[0;tempy(1:idx)];
            
            [Model,~]=fit(x(:),y(:),fittype(Fun),fo);
            xx=linspace(0,x(end),100);
            yy=feval(Model,xx);
            
                    figure(1),clf
                    plot(x,y,'ko'),hold on
                    plot(xx,yy,'r')
                    drawnow
%                     pause
            FIcurveParamRecord{kk}(k,:)=[-(Model.b/Model.a),Model.a,maxStim,maxFreq];
        else
            FIcurveParamRecord{kk}(k,:)=nan(1,4);
        end


        %% ------------- Bin FIcurve ---------------------------
        
        for i=1:size(BinsFIcurve{kk},1)
            idx=find(StimAmp>=BinsFIcurve{kk}(i,1)&StimAmp<BinsFIcurve{kk}(i,2));
            if ~isempty(idx)
                spkFreqBin{kk}(k,i)=mean(spkFreq(idx));
            else
                spkFreqBin{kk}(k,i)=nan;
            end
        end
        
                figure(100)
                plot(StimAmp,spkFreq,'o','markeredgecolor',colormap(kk,:),'markersize',6),hold on
                plot(xx,yy,'m')
                xlabel('Current(pA)');ylabel('Frequency(Hz)');
                drawnow
%                 pause
        
        %% ------------- Bin spkshape,1~5 AP first spike ---------------------------
        
        Wantidx=find((spkFreq>=(Num-NPlusMinus)*2&spkFreq<=(Num+NPlusMinus)*2)&(StimAmp(:)<1550));

        if ~isempty(Wantidx)
            spkwaveform_temp=[];
            HWtemp=[];
            AHPtemp=[];
            Thresholdtemp=[];
            Amptemp=[];
            max_slopetemp=[];
            min_slopetemp=[];
            
            for jj=1:length(Wantidx)
                data=waveformV(:,Wantidx(jj));
                [spkwaveform,ttspanspk,HW, AHP, Threshold, Amp,max_slope,min_slope,...
                 HWidx,AHPidx, Thresholdidx,...
                 max_slopeidx,min_slopeidx] = AP_Statistic_HERE(data,spkTime{Wantidx(jj)},tspanV,1:Num-NPlusMinus);%AP_Statistic_HERE里的pretime和postime会影响ttspanspk   
             
                spkwaveform_temp(jj,:,:)=spkwaveform;
                HWtemp(:,jj)=HW;
                AHPtemp(:,jj)=AHP;
                Thresholdtemp(:,jj)=Threshold;
                Amptemp(:,jj)=Amp;
                max_slopetemp(:,jj)=max_slope;
                min_slopetemp(:,jj)=min_slope;

                figure(1),clf
                H=[];
                for j=1%:size(spkwaveform,2)
%                     H(j)=subplot(1,size(spkwaveform,2),j);
                    plot(ttspanspk,spkwaveform(:,j)),hold on
                    plot(ttspanspk(HWidx(j,:)),spkwaveform(HWidx(j,:),j),'ro')
                    plot(ttspanspk(AHPidx(j)),spkwaveform(AHPidx(j),j),'ko')
                    plot(ttspanspk(Thresholdidx(j)),spkwaveform(Thresholdidx(j),j),'co')
                    axis tight
                    drawnow
                end
%                 linkaxes(H,'x')
%                 pause
            end          
    
            spkparamRecord{kk}{1}(:,k)=mean(HWtemp,2);
            spkparamRecord{kk}{2}(:,k)=mean(AHPtemp,2);
            spkparamRecord{kk}{3}(:,k)=mean(Thresholdtemp,2);
            spkparamRecord{kk}{4}(:,k)=mean(Amptemp,2);
            spkparamRecord{kk}{5}(:,k)=mean(max_slopetemp,2);
            spkparamRecord{kk}{6}(:,k)=mean(min_slopetemp,2);
            
            if size(spkwaveform_temp,2)==10000%%这个值需要根据ttspanspk修改
                spkwaveformRecord{kk}(k,:,:)=squeeze(mean(spkwaveform_temp,1));
                ttspanspkGood=ttspanspk;
            else
                spkwaveformRecord{kk}(k,:,:)=nan(10000,Num-NPlusMinus);
            end
        else
            spkparamRecord{kk}{1}(:,k)=nan(Num-NPlusMinus,1);%HW
            spkparamRecord{kk}{2}(:,k)=nan(Num-NPlusMinus,1);%AHP
            spkparamRecord{kk}{3}(:,k)=nan(Num-NPlusMinus,1);%Threshold
            spkparamRecord{kk}{4}(:,k)=nan(Num-NPlusMinus,1);%Amp
            spkparamRecord{kk}{5}(:,k)=nan(Num-NPlusMinus,1);%maxslope
            spkparamRecord{kk}{6}(:,k)=nan(Num-NPlusMinus,1);%minslope
            spkwaveformRecord{kk}(k,:,:)=nan(10000,Num-NPlusMinus);
        end
        
    end
end
% ttspanspk=ttspanspkGood;
FIcurveParamName={'rheobase','gain','maxStim','maxFreq'};
spkParamName={'HW', 'AHP', 'Threshold', 'Amp','max_slope','min_slope'};

save('Results_FIcurveStatistic.mat','Fidx','GroupName',...
    'spkwaveformRecord','spkparamRecord','ttspanspk',...
    'spkFreqBin','BinsFIcurve','FIcurveParamRecord',...
    'Num','thAP',...
    'FIcurveParamName','spkParamName','GroupName')

for kk=1:length(Fidx)
    
    xlswrite('Results_Statistic.xls',FilenameRecord{kk},['Param_',GroupName{kk}],'A2')
    xlswrite('Results_Statistic.xls',FIcurveParamName,['Param_',GroupName{kk}],'B1')
    xlswrite('Results_Statistic.xls',FIcurveParamRecord{kk},['Param_',GroupName{kk}],'B2')
    
    xlswrite('Results_Statistic.xls',spkParamName,['Param_',GroupName{kk}],'F1')
    xlswrite('Results_Statistic.xls',[spkparamRecord{kk}{1}(thAP,:)',...
                                      spkparamRecord{kk}{2}(thAP,:)',...
                                      spkparamRecord{kk}{3}(thAP,:)',...
                                      spkparamRecord{kk}{4}(thAP,:)',...
                                      spkparamRecord{kk}{5}(thAP,:)',...
                                      spkparamRecord{kk}{6}(thAP,:)'],['Param_',GroupName{kk}],'F2')
    
    xlswrite('Results_Statistic.xls',BinsFIcurve{kk}(:,1),['FIcurve_',GroupName{kk}],'A1')
    xlswrite('Results_Statistic.xls',spkFreqBin{kk}',['FIcurve_',GroupName{kk}],'B1')
end

%%

clear all;clc
load('Results_FIcurveStatistic.mat')
colormap=[0.5,0.54,0.53;0.12,0.56,1];
% ==================== FIcurve ====================
figure(1),clf
for kk=1:length(Fidx)
    
    x=BinsFIcurve{kk}(:,1);
    y=nanmean(spkFreqBin{kk},1);
    N=sum(~isnan(spkFreqBin{kk}),1);
    sem=nanstd(spkFreqBin{kk},[],1)./sqrt(N);
    
    h=errorbar(x,y,sem,'r-o','color',colormap(kk,:),'markerfacecolor',colormap(kk,:),'markersize',8);hold on
    h.CapSize=0;
    
%     plot(x,spkFreqBin{kk},'r-o','color',colormap(kk,:),'markerfacecolor','w','markersize',4),hold on
    
end
legend(GroupName)
ylabel('Frequency (Hz)')
xlabel('Current injection (pA)')

% for i=1:size(spkFreqBin{1},2)
%     A=spkFreqBin{1}(:,i);A=A(~isnan(A));
%     B=spkFreqBin{2}(:,i);B=B(~isnan(B));
%     if ~isempty(A)&~isempty(B)
%         %         [~,p] = ttest2(A,B);
%         p = ranksum(A,B);
%
%         figure(1)
%         if p<0.05&p>=0.01;
%             text(BinsFIcurve{1}(i,1),y(i),'*','fontsize',20);
%         elseif p<0.01&p>=0.001;
%             text(BinsFIcurve{1}(i,1),y(i),'**','fontsize',20);
%         elseif p<0.001
%             text(BinsFIcurve{1}(i,1),y(i),'***','fontsize',20);
%         end
%     end
% end


% ====================  statistic ====================
figure(4),clf
figure(400),clf
for k=1:length(FIcurveParamName)
    
    dataFortest=[];
    groupFortest=[];
    
    Str={};
    for kk=1:length(Fidx)
        data{kk}=FIcurveParamRecord{kk}(:,k);N=sum(~isnan(data{kk}));
        Mean=nanmean(data{kk});
        Sem=nanstd(data{kk})./sqrt(N);
        figure(4)
        subplot(1,length(FIcurveParamName),k)
        bar((2.5+2.5*(kk-1)),Mean,'facecolor',colormap(kk,:));hold on
        errorbar((2.5+2.5*(kk-1)),Mean,Sem,'color',colormap(kk,:))
        ylabel(FIcurveParamName{k})
        text(((2.5+2.5*(kk-1))-0.3),Mean/2,num2str(N),'color','k','fontsize',12)
        set(gca,'tickdir','out','box','off','ticklength',[0.02,0.02])
        Str{kk}=[GroupName{kk},' ',num2str(Mean),' \pm ',num2str(Sem)];
        plot((1.25+2.5*(kk-1))+0.1*rand(size(data{kk})),data{kk},'bo','markerfacecolor','none','MarkerEdgeColor',colormap(kk,:),'markersize',6,'lineWidth',1)

        
        figure(400)
        subplot(length(Fidx),length(FIcurveParamName),k+(kk-1)*length(FIcurveParamName))
        [n,x]=hist(data{kk},10);%bin number
        bar(x,n,'facecolor',colormap(kk,:))
        ylabel(FIcurveParamName{k})
        set(gca,'tickdir','out','box','off','ticklength',[0.02,0.02])
        
        dataFortest=[dataFortest;data{kk}];
        groupFortest=[groupFortest;kk.*ones(size(data{kk}));];        
        
    end
%     % ========== test
%     % 更加严格
%     %  [p,tbl,stats] = kruskalwallis(dataFortest,groupFortest,'off');
%     %  c = multcompare(stats,'display','off')
%     %  pRecord=[mean([c(1,6),c(4,6)]),mean([c(3,6),c(6,6)])];
%     
%     [p,tbl,stats] = anova1(dataFortest,groupFortest,'off');
%     c = multcompare(stats,'display','off')
%     pRecord=[mean([c(1,6)]),mean([c(3,6),c(5,6),c(6,6)])];
%     % ==========
%     Str{kk+1}=['p = ',num2str(p)];
%     figure(4)
%     title(Str)

end

%
figure(5),clf
figure(500),clf
for k=1:length(spkParamName)
    
    for kk=1:length(Fidx)
        data{kk}=spkparamRecord{kk}{k}(thAP,:);N=sum(~isnan(data{kk}));%统计7±3的AP串中第几个AP的波形参数
        Mean=nanmean(data{kk});
        Sem=nanstd(data{kk})./sqrt(N);
        figure(5)
        subplot(1,length(spkParamName),k)
        bar((2.5+2.5*(kk-1)),Mean,'facecolor',colormap(kk,:));hold on
        errorbar((2.5+2.5*(kk-1)),Mean,Sem,'color',colormap(kk,:))
        text(((2.5+2.5*(kk-1))-0.3),Mean/2,num2str(N),'color','k','fontsize',12)
        ylabel(spkParamName{k})
        plot((1.25+2.5*(kk-1))+0.1*rand(size(data{kk})),data{kk},'bo','markerfacecolor','none','MarkerEdgeColor',colormap(kk,:),'markersize',6,'lineWidth',1)   
       
        figure(500)
        subplot(length(Fidx),length(spkParamName),k+(kk-1)*length(spkParamName))
        [n,x]=hist(data{kk},10);%bin number
        bar(x,n,'facecolor',colormap(kk,:))
        ylabel(spkParamName{k})
    end
%     ========== test
% try%加这个try和end可以尝试test
%         A=data{1};B=data{end};
%         if swtest(A,0.05)==0&swtest(B,0.05)==0
%             [~,p]=ttest2(A,B);
%             tmethod='two samples T-test';
%         else
%             [p,~]=ranksum(A,B);
%             tmethod='Wilcoxon rank sum ';% or Mann-Whitney U-test
%         end
%         title({['No ',num2str(nanmean(A)),'\pm',num2str(nanstd(A)./sqrt(length(A)))]
%             ['YES ',num2str(nanmean(B)),'\pm',num2str(nanstd(B)./sqrt(length(B)))]
%             [tmethod]
%             ['  p = ',num2str(p)]})
% end
    % ==========
end


% ====================  waveform  ====================
figure(6),clf

xlim_max=0.008;
xlim_min=-0.005;

for kk=1:length(Fidx)
    
    for j=1:size(spkwaveformRecord{kk},3)
        subplot(2,size(spkwaveformRecord{kk},3),j);
        data=squeeze(nanmean(spkwaveformRecord{kk}(:,:,j),1));% from averaged
%         data=squeeze(nanmean(spkwaveformRecord{kk}(5,:,j),1));% from single cell
        plot(ttspanspk,data,'-','color',colormap(kk,:)),hold on%%ttspanspk与数据采样率相关
        axis([xlim_min,xlim_max,-60,60])
        xlabel('Time(s)');ylabel('Vm (mV)');
        
        subplot(2,size(spkwaveformRecord{kk},3),j+size(spkwaveformRecord{kk},3))
        tempidx=find(ttspanspk>-0.004&ttspanspk<0.0040);
        plot(data(tempidx(2:end)),diff(data(tempidx))*100,'color',colormap(kk,:)),hold on
        axis([-60,60,-150,500])
        xlabel('Vm (mV)');ylabel('dV/dt');
    end
    
end
legend(GroupName,'location','best')

%}