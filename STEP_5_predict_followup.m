clear all
close all

%select IDPname before running
IDPname='GMV'; % GMV, CTh, SA or FA

addpath('Dependencies/');
inpath1=('variables_for_normative_modeling/');
inpath2=('output_normative_modeling/');

%load cross-validation structure (index of subjects in training and test)
instr=strcat(['ukb_',IDPname,'_vars_for_R_cross_TRAIN_1']);
load([inpath1,instr],'c','sex','site')
K_folds=2;

N=length(sex);

%initialise
ypredNaive=nan(N,1);
ypredX=nan(N,1);
ypredX_PERCENTILE=nan(N,1);

errorNaive=nan(N,1);
errorX=nan(N,1);
errorX_PERCENTILE=nan(N,1);

age1=nan(N,1);
age2=nan(N,1);
y1=nan(N,1);
y2=nan(N,1);
time_interval=nan(N,1);

for k=1:K_folds
    
    %read in normative modeling results: cross-sectional
    yfitF=nan(1000,5,length(unique(site)));
    yfitM=nan(1000,5,length(unique(site)));
    for i=1:length(unique(site))
        
        load([inpath2,...
            'GAMLSS_',IDPname,'_vars_for_R_cross_predicted_sex0_site',num2str(i-1),'_TRAIN_',num2str(k),'.mat'],'age','predictions_quantiles');
        yfitF(:,:,i)=predictions_quantiles; %estimated females at site 0
        
        
        load([inpath2,...
            'GAMLSS_',IDPname,'_vars_for_R_cross_predicted_sex1_site',num2str(i-1),'_TRAIN_',num2str(k),'.mat'],'age','predictions_quantiles');
        yfitM(:,:,i)=predictions_quantiles; %estimated males at site 0
    end
    
    
    
    
    %read in normative modeling results: longitudinal (rate of change)
    yfitF_rate=nan(1000,5,length(unique(site)));
    yfitM_rate=nan(1000,5,length(unique(site)));
    for i=1:length(unique(site))
        load([inpath2,...
            'GAMLSS_',IDPname,'_vars_for_R_long_predicted_sex0_site',num2str(i-1),'_TRAIN_',num2str(k),'.mat'],'age','predictions_quantiles');
        yfitF_rate(:,:,i)=predictions_quantiles; %estimated females at site 0
        
        load([inpath2,...
            'GAMLSS_',IDPname,'_vars_for_R_long_predicted_sex1_site',num2str(i-1),'_TRAIN_',num2str(k),'.mat'],'age','predictions_quantiles');
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
        
        est_rateF(:,:,i)=diff(yfitF(:,:,i))/delta;
        est_rateM(:,:,i)=diff(yfitM(:,:,i))/delta;
        est_rateF_percent(:,:,i)=est_rateF(:,:,i)./yfitF(1:end-1,:,i)*100;
        est_rateM_percent(:,:,i)=est_rateM(:,:,i)./yfitM(1:end-1,:,i)*100;
        
    end
    
    %% PREDICTION %%
    instr=strcat([IDPname,'_vars_for_R_cross_TRAIN_',num2str(k)]);
    load([inpath1,instr]);
    
    %naive prediction of no change
    CV_test_ypredNaive=CV_test_y1;
    CV_test_errorNaive=CV_test_y2-CV_test_y1;
    
    %%load in 100 percentile data for individual prediction of cross-sectionally
    %inferred rates of change
    yfitF_100=nan(1000,99,length(unique(site)));
    yfitM_100=nan(1000,99,length(unique(site)));
    for i=1:length(unique(site))
        load([inpath2,...
            'GAMLSS_100_',IDPname,'_vars_for_R_cross_predicted_sex0_site',num2str(i-1),'_TRAIN_',num2str(k),'.mat'],'predictions_quantiles');
        yfitF_100(:,:,i)=predictions_quantiles; %estimated females at site 0
        
        
        load([inpath2,...
            'GAMLSS_100_',IDPname,'_vars_for_R_cross_predicted_sex1_site',num2str(i-1),'_TRAIN_',num2str(k),'.mat'],'predictions_quantiles');
        yfitM_100(:,:,i)=predictions_quantiles; %estimated males at site 0
        
        
    end
    
    %Estimate rate of change from cross-sectional measurements using
    %individual percentile at baseline
    est_rateF_100=nan(999,99,length(unique(site)));
    est_rateM_100=nan(999,99,length(unique(site)));
    est_rateF_percent_100=nan(999,99,length(unique(site)));
    est_rateM_percent_100=nan(999,99,length(unique(site)));
    for i=1:length(unique(site))
        est_rateF_100(:,:,i)=diff(yfitF_100(:,:,i))/delta;
        est_rateM_100(:,:,i)=diff(yfitM_100(:,:,i))/delta;
        est_rateF_percent_100(:,:,i)=est_rateF_100(:,:,i)./yfitF_100(1:end-1,:,i)*100;
        est_rateM_percent_100(:,:,i)=est_rateM_100(:,:,i)./yfitM_100(1:end-1,:,i)*100;
        
    end
    
    %Predict value at follow-up using cross-sectionally predicted rate of change
    CV_test_ypredX=zeros(length(CV_test_age1),1);
    CV_test_ypredX_PERCENTILE=zeros(length(CV_test_age1),1);
    ind_percentile=zeros(length(CV_test_age1),1);
    diff_perc=zeros(length(CV_test_age1),99);
    for i=1:length(CV_test_age1)
        
        ind=site(i)+1;
        
        est_rateF_sub=est_rateF(:,:,ind);
        est_rateM_sub=est_rateM(:,:,ind);
        est_rateF_sub_100=est_rateF_100(:,:,ind);
        est_rateM_sub_100=est_rateM_100(:,:,ind);
        
        start_age=CV_test_age1(i);
        end_age=CV_test_age2(i);
        [~,ind0]=min(abs(agemid-start_age));
        [~,ind1]=min(abs(agemid-end_age));
        CV_test_sex=sex(c.test(k));
        
        if CV_test_sex(i)==0 %females
            
            %find difference between baseline value and each percentile
            diff_perc(i,:)=CV_test_y1(i)-yfitF_100(ind0,:,ind);
            
            %find which percentile the individual belongs to at baseline
            ind_percentile(i)=find(abs(diff_perc(i,:))==min(abs(diff_perc(i,:))));
            
            CV_test_ypredX_PERCENTILE(i)=trapz(agemid(ind0:ind1),est_rateF_sub_100(ind0:ind1,ind_percentile(i)))+CV_test_y1(i);
            CV_test_ypredX(i)=trapz(agemid(ind0:ind1),est_rateF_sub(ind0:ind1,3))+CV_test_y1(i);
        
        else %males
            
            %find diff_percerence between baseline value and each percentile
            diff_perc(i,:)=CV_test_y1(i)-yfitM_100(ind0,:,ind);
            
            %find which percentile the individual belongs to at baseline
            ind_percentile(i)=find(abs(diff_perc(i,:))==min(abs(diff_perc(i,:))));
            
            CV_test_ypredX_PERCENTILE(i)=trapz(agemid(ind0:ind1),est_rateM_sub_100(ind0:ind1,ind_percentile(i)))+CV_test_y1(i);
            CV_test_ypredX(i)=trapz(agemid(ind0:ind1),est_rateM_sub(ind0:ind1,3))+CV_test_y1(i);
        end
    end
    
    
    CV_test_errorX=CV_test_ypredX-CV_test_y2;
    CV_test_errorX_PERCENTILE=CV_test_ypredX_PERCENTILE-CV_test_y2;
   
    
    %error
    errorNaive(c.test(k))=CV_test_errorNaive; %Naive model of no change
    errorX(c.test(k))=CV_test_errorX; %cross-sectional prediction (50th percentile)
    errorX_PERCENTILE(c.test(k))=CV_test_errorX_PERCENTILE; %cross-sectional prediction (individual percentile)
    
    %predicted value
    ypredNaive(c.test(k))=CV_test_ypredNaive; %Naive model of no change
    ypredX(c.test(k))=CV_test_ypredX; %cross-sectional prediction (50th percentile)
    ypredX_PERCENTILE(c.test(k))=CV_test_ypredX_PERCENTILE; %cross-sectional prediction (individual percentile)
    
    %reorganise data
    age1(c.test(k))=CV_test_age1;
    age2(c.test(k))=CV_test_age2;
    sex(c.test(k))=CV_test_sex;
    y1(c.test(k))=CV_test_y1;
    y2(c.test(k))=CV_test_y2;
    CV_test_time_interval=CV_test_age2-CV_test_age1;
    time_interval(c.test(k))=CV_test_time_interval;
end




%plot errors
hf=figure; hf.Color='w'; hf.Position=[50,50,300,400];
hb=bar([mean(abs(errorNaive)),mean(abs(errorX)),mean(abs(errorX_PERCENTILE))]);
ylabel('Mean absolute error');
title('Females');
ax=gca;
ax.XTickLabel{1}='Naive';
ax.XTickLabel{2}='Cross-sectional';
ax.XTickLabel{3}='Cross-sectional indv percentiles';
ax.XTickLabelRotation=45;



%plot predicted vs observed for the 3 models 
hf=figure; hf.Color='w'; hf.Position=[100,100,1600,500];
subplot(1,3,1);
for i=1:length(age1)
    plot([0;age2(i)-age1(i)],[0;y2(i)-y1(i)],'r');
    if i==1
        hold on;
    end
end
for i=1:length(age1)
    plot([0;age2(i)-age1(i)],[0;ypredNaive(i)-y1(i)],'k');
end
xlabel('Time between baseline and followup (years)');
ylabel('Change (mm)');
title('red=measured, black=predicted (Naive)');
ax=gca;
ax.XTickLabelRotation=45;

subplot(1,3,2);
for i=1:length(age1)
    plot([0;age2(i)-age1(i)],[0;y2(i)-y1(i)],'r');
    if i==1
        hold on;
    end
end
for i=1:length(age1)
    plot([0;age2(i)-age1(i)],[0;ypredX(i)-y1(i)],'k');
end
xlabel('Time between baseline and followup (years)');
ylabel('Change (mm)');
title('red=measured, black=predicted (cross-sectional)');
ax=gca;
ax.XTickLabelRotation=45;


subplot(1,3,3);
for i=1:length(age1)
    plot([0;age2(i)-age1(i)],[0;y2(i)-y1(i)],'r');
    if i==1
        hold on;
    end
end
for i=1:length(age1)
    plot([0;age2(i)-age1(i)],[0;ypredX_PERCENTILE(i)-y1(i)],'k');
end
xlabel('Time between baseline and followup (years)');
ylabel('Change (mm)');
title('red=measured, black=predicted (cross percentile)');
ax=gca;
ax.XTickLabelRotation=45;




