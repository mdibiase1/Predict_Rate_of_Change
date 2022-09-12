function STEP_1_prepare_vars_for_norm_models(IDPname, y1, y2, age1, age2, sex, site, covars_visit_1, covars_visit_2)

addpath('Dependencies/');

% Inputs:
% IDPname = GMV, CTh, SA or FA
% y1 = IDP at baseline
% y2 = IDP at follow-up
% age1 = age at baseline
% age2 = age at follow-up
% sex = biological attribute of sex
% site = scanner/site
% neuroimaging covariates (matrix)

% Outputs:
% variables prepared for normative modeling 

%% STEP 1: main normative models (as shown in Fig. 1) %%

%normalise covariates
%visit 1
covars_visit_1_norm=zeros(size(covars_visit_1));
for i=1:size(covars_visit_1,2)
    covars_visit_1_norm(:,i)=inormal(covars_visit_1(:,i));
end

%visit 2
covars_visit_2_norm=zeros(size(covars_visit_2));
for i=1:size(covars_visit_2,2)
    covars_visit_2_norm(:,i)=inormal(covars_visit_2(:,i));
end


%derive neuroimaging covariate components for cross-sectional normative models
[score1] = pca(covars_visit_1_norm); %run PCA
covars_1_pca=score1(:,1:3); %extract first 3 components


%derive neuroimaging covariate components for longitudinal normative models
%(utilise covariates from both timepoints)
[score2] = pca([covars_visit_1_norm covars_visit_2_norm]);  %run PCA
covars_2_pca=score2(:,1:3); %extract first 3 components

%save variables for cross-sectional normative models
dir='variables_for_normative_modeling/';
mkdir(dir);
outstr=strcat([IDPname,'_vars_for_R_cross']);
save([dir,outstr], 'y1','age1','sex','site','covars_1_pca')

%Compute rate of change using longtudinal measurements
rate=(y2-y1)./(age2-age1);
ageM=(age1+age2)/2; %average age between baseline and follow-up

%save variables for longitudinal normative models
%need to include >1 site variable if site changes across timepoints
outstr=strcat([IDPname,'_vars_for_R_long']);
save([dir,outstr], 'rate','ageM','sex','site','covars_2_pca')



%% STEP 2: Prepare data for cross-validation for prediction analyses %%

N=size(y1,1); %number of subjects
K=2; %number of folds
c=cvpartition(N,'KFold',K); %generate a 2-fold partition


for i=1:K
    
    %CROSS-SECTIONAL
    %training set
    CV_y1=y1(c.training(i));
    CV_y2=y2(c.training(i));
    CV_age1=age1(c.training(i));
    CV_age2=age2(c.training(i));
    CV_sex=sex(c.training(i));
    CV_site=site(c.training(i));
    CV_covars_1_pca=covars_1_pca(c.training(i),:);
    
    %test set
    CV_test_y1=y1(c.test(i));
    CV_test_y2=y2(c.test(i));
    CV_test_age1=age1(c.test(i));
    CV_test_age2=age2(c.test(i));
    CV_test_sex=sex(c.test(i));
    CV_test_site=site(c.test(i));
    CV_test_covars_1_pca=covars_1_pca(c.test(i),:);
    
    %save variables for cross-sectional normative models with cross-validation
    outstr=strcat([IDPname,'_vars_for_R_cross_TRAIN_']);
    outstr=strcat(outstr,num2str(i));
    save([dir,outstr], 'c','CV_y1','CV_y2','CV_age1','CV_age2','CV_sex','CV_site','CV_covars_1_pca',...
        'CV_test_y1','CV_test_y2','CV_test_age1','CV_test_age2','CV_test_sex','CV_test_site','CV_test_covars_1_pca');
    
    %LONGITUDINAL
    %train
    CV_rate=rate(c.training(i));
    CV_ageM=ageM(c.training(i));
    CV_covars_2_pca=covars_2_pca(c.training(i),:);
    
    %test
    CV_test_rate=rate(c.test(i));
    CV_test_ageM=ageM(c.test(i));
    CV_test_covars_2_pca=covars_2_pca(c.test(i),:);
    
    %save variables for longitudinal normative models with cross-validation
    outstr=strcat(['ukb_',IDPname,'_vars_for_R_rate_long_TRAIN_']);
    outstr=strcat(outstr,num2str(i));
    save([dir,outstr], 'CV_rate','CV_ageM','CV_sex','CV_site','CV_covars_2_pca',...
        'CV_test_rate','CV_test_ageM','CV_test_sex','CV_test_site','CV_test_covars_2_pca')
    
end

