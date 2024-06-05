close all;
clear all;clc
% dbstop if error
% dbstop at 45 in AP_Statistic

Fidx{1}=dir('DATA_C*.mat');
Fidx{2}=dir('DATA_T*.mat');

GroupName={'C','T'};
% 
% WantState=2;%1/2：分析with/without holding状态下的数据
NumbWant=1;

BinWidth=50;
BinL=0:BinWidth:900-BinWidth;
BinR=BinL+BinWidth;
Bin=[BinL(:),BinR(:)];
%%
colormap=[0.5,0.54,0.53;0.12,0.56,1;
    1,0.38,0;0.69,0.19,0.12];
%深灰色[0.5,0.54,0.53];红色[0.69,0.19,0.12];橙色[1,0.38,0];浅蓝色[0.12,0.56,1];

figure(100),clf
for kk=1:length(Fidx)
    for k=1:length(Fidx{kk})
        
        filename=Fidx{kk}(k).name;
        load(filename);
        dt=tspanV(2)-tspanV(1);
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        FilenameRecord{kk}{k,1}=filename;
        FilenameRecord_row{kk}{1,k}=filename;
       
        %%
%         if WantState==1
%             idx=abs(Ihold)>=10;
%         else
%              idx=abs(Ihold)<10;

%         end
         idx=[]
        TimeStim=TimeStim(idx);
        %% ----------------------------------------
        
        spkwaveform_temp=[];
        spkparam_temp=[];
        
        for jj=1:length(TimeStim)
            data=waveformV(:,jj);
             ihold=0;
            [spkwaveform,ttspanspk,HW, AHP, Threshold, Amp,max_slope,min_slope,AIS_maxslope,Soma_maxslope...
                HWidx,AHPidx, Thresholdidx,max_slopeidx,min_slopeidx,max_slope2idx] = AP_Statistic(data,spkTime{jj},tspanV,1:NumbWant);
            spkwaveform_temp(:,:,jj)=spkwaveform;
            spkparam_temp(1,:,jj)=HW;
            spkparam_temp(2,:,jj)=AHP;
            spkparam_temp(3,:,jj)=Threshold;
            spkparam_temp(4,:,jj)=Amp;
            spkparam_temp(5,:,jj)=max_slope;
            spkparam_temp(6,:,jj)=min_slope;
            spkparam_temp(7,:,jj)=AIS_maxslope;
            spkparam_temp(8,:,jj)=Soma_maxslope;
            
            figure(100)
            for ii=1:NumbWant
                H(1)=subplot(2,NumbWant,ii);
                plot(ttspanspk,spkwaveform_temp(:,ii,jj)),hold on
                plot(ttspanspk(HWidx(ii,:)),spkwaveform_temp(HWidx(ii,:),ii,jj),'ro')
                plot(ttspanspk(AHPidx(ii)),spkwaveform_temp(AHPidx(ii),ii,jj),'ko')
                plot(ttspanspk(Thresholdidx(ii)),spkwaveform_temp(Thresholdidx(ii),ii,jj),'co')
                
                H(2)=subplot(2,NumbWant,ii+NumbWant);
                slope=[nan;diff(spkwaveform_temp(:,ii,jj))/(dt*1000)];
                slope2=[nan;diff(slope)/(dt*1000)];
                slope2=smooth(slope2);
                plot(ttspanspk,slope2','k'),hold on
                plot(ttspanspk,slope','r'),hold on
                plot(ttspanspk,slope2','b'),hold on
                plot(ttspanspk(max_slope2idx(ii)),slope2(max_slope2idx(ii)),'ro')
            end
            axis tight
            drawnow
        end
        %%
        spkparam=[];
        spkwaveform=[];
        Otherparam=[];
        for bb=1:size(Bin,1)
            idx=find(TimeStim>Bin(bb,1)&TimeStim<=Bin(bb,2));
            for pp=1:size(spkparam_temp,1)
                spkparam(bb,pp,:)=mean(spkparam_temp(pp,:,idx),3);
                spkwaveform(bb,:,:)=mean(spkwaveform_temp(:,:,idx),3);
            end
            Otherparam(bb,1)=mean(Vrest(idx));
            Otherparam(bb,2)=mean(spkFreq(idx));
        end
        
        %%
        spkparamRecord{kk}(:,:,:,k)=spkparam;
        spkwaveformRecord{kk}(:,:,:,k)=spkwaveform;%第一个spike waveform
        OtherparamRecord{kk}(:,:,k)=Otherparam;
    end
    spkParamName={'HW', 'AHP', 'Threshold',...
        'Amp','max slope','min slope',...
        'AIS_maxslope','Soma_maxslope'};
    
    OtherparamName={'Vrest','spkFreq'};
    
    save('Results_FIcurve.mat','Fidx','GroupName','colormap',...
        'tspanV','TimeStim',...
        'spkwaveformRecord',...
        'spkparamRecord','spkParamName',...
        'OtherparamRecord','OtherparamName',...
        'NumbWant','Bin',...
        'spkParamName','GroupName')
end
%%

figure(1),clf
for pp=1:length(spkParamName)
    for kk=1:length(GroupName)
        for nn=1%:NumbWant% nth spike
            
            subplot(NumbWant,length(spkParamName),pp+(nn-1)*length(spkParamName))
            
            data=squeeze(spkparamRecord{kk}(:,pp,nn,:));
            
            Mean=mean(data,2);
            N=sum(~isnan(data),2);
            Sem=std(data,[],2)./sqrt(N);
            errorbar(mean(Bin,2),Mean,Sem,'color',colormap(kk,:)),hold on
            ylabel(spkParamName{pp})
        end
        
    end
end

figure(2),clf
for pp=1:length(OtherparamName)
    for kk=1:length(GroupName)
           
            subplot(1,length(OtherparamName),pp)
            
            data=squeeze(OtherparamRecord{kk}(:,pp,:));
            
            Mean=mean(data,2);
            N=sum(~isnan(data),2);
            Sem=std(data,[],2)./sqrt(N);
            errorbar(mean(Bin,2),Mean,Sem,'color',colormap(kk,:)),hold on
            ylabel(OtherparamName{pp})
            xlabel('Time (s)')
        
    end
end


