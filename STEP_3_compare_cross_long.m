clear all
close all

%select IDPname before running
IDPname='GMV'; % GMV, CTh, SA or FA

addpath('Dependencies/');
inpath1=('variables_for_normative_modeling/');
inpath2=('output_normative_modeling/');

load([inpath1,IDPname,'_vars_for_R_cross'],'site','age1');

%read in normative modeling results: cross-sectional 
yfitF=nan(1000,5,length(unique(site)));
yfitM=nan(1000,5,length(unique(site)));
for i=1:length(unique(site))
    
    load([inpath2,...
        'GAMLSS_',IDPname,'_cross_predicted_sex0_site',num2str(i-1),'_WHOLE_SAMPLE.mat'],'predictions_quantiles');
    yfitF(:,:,i)=predictions_quantiles; %estimated females at site i (assumes females = 0)
    
    load([inpath2,...
        'GAMLSS_',IDPname,'_vars_for_R_cross_predicted_sex1_site',num2str(i-1),'_WHOLE_SAMPLE.mat'],'age','predictions_quantiles');
    yfitM(:,:,i)=predictions_quantiles; %estimated males at site i (assumes males = 1)
end


%read in normative modeling results: longitudinal (rate of change)
yfitF_rate=nan(1000,5,length(unique(site)));
yfitM_rate=nan(1000,5,length(unique(site)));
for i=1:length(unique(site))
    
    load([inpath2,...
        'GAMLSS_',IDPname,'_vars_for_R_long_predicted_sex0_site',num2str(i-1),'_WHOLE_SAMPLE.mat'],'age','predictions_quantiles');
    yfitF_rate(:,:,i)=predictions_quantiles; %estimated females at site 0
    
    load([inpath2,...
        'GAMLSS_',IDPname,'_vars_for_R_long_predicted_sex1_site',num2str(i-1),'_WHOLE_SAMPLE.mat'],'age','predictions_quantiles');
    yfitM_rate(:,:,i)=predictions_quantiles; %estimated males at site 0
    
end

 
%Estimate rate of change from cross-sectional measurements
J=length(age);
delta=age(2)-age(1);
est_rateF=nan(999,5,length(unique(site)));
est_rateM=nan(999,5,length(unique(site)));
est_rateF_percent=nan(999,5,length(unique(site)));
est_rateM_percent=nan(999,5,length(unique(site)));
agemid=(age(2:end)+age(1:end-1))/2; %midpoint
for i=1:length(unique(site))
    
    %Estimate rate of change from cross-sectional measurements
    est_rateF(:,:,i)=diff(yfitF(:,:,i))/delta;
    est_rateM(:,:,i)=diff(yfitM(:,:,i))/delta;
    est_rateF_percent(:,:,i)=est_rateF(:,:,i)./yfitF(1:end-1,:,i)*100;
    est_rateM_percent(:,:,i)=est_rateM(:,:,i)./yfitM(1:end-1,:,i)*100;
    
end

%plot results (average across males and females and across sites)
hf=figure; hf.Color='w'; hf.Position=[50,50,800,400];
subplot(1,2,1)
yfit=(mean(yfitM,3)+mean(yfitF,3))/2;
scatter(age1,y1); 
plot(age,yfit(:,3),'k','linewidth',2.5); %plot 50th centile
plot(age,yfit(:,[1,2,4,5]),'k:','linewidth',2); %plot other centiles
ylabel(IDPname); 
xlabel('Age');
title('Cross-sectional');

subplot(2,2,1)
yfit_rate=(mean(yfitM_rate,3)+mean(yfitF_rate,3))/2;
est_rate=(est_rateM_percent+est_rateF_percent)/2;
observed=yfit_rate./yfit*100;

S=500; %smoothing span, for visualization only 
plot(age,observed,'k','linewidth',2.5); hold on; 
plot(agemid,smooth(est_rate,S),'k--','linewidth',2); hold on; 
L(1) = plot(nan, nan, 'k','linewidth',1); L(2) = plot(nan, nan, 'k--','linewidth',2.5);
legend(L, {'Longitudinal', 'Cross-sectional prediction'},'Location','SouthWest','Box','off'); 
ylabel('Rate (%/year)'); 
title('Rate (%/year)');
xlabel('Age');


