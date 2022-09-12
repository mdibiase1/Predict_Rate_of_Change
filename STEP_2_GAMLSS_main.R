rm(list = ls())
setwd("variables_for_normative_modeling/")
source('Dependencies/helpers.R') 
library(ggpointdensity)
library(cowplot)
library(patchwork)
library(qqplotr)
library(R.matlab)
library(Hmisc)
library(rstatix)
library(ggpubr)
library(scoring)
library(gamlss)
library(Hmisc)
library(pracma)

# read data from mat file
modality<-c("GMV","CTh","SA","FA")
type<-c("cross","long")

  # run analysis for each modality and each type (i.e., cross-sectional and longitudinal)
  for (idx_m in 1:length(modality)){ 
     for (idx_type in 1:length(type)){
      fname<-paste(modality[idx_m],"_vars_for_R_",type[idx_type],".mat",sep="")
      M<-readMat(fname)
      
      # read data
      if (idx_type==1){
        # baseline data
        print(paste("baseline raw data of ",modality[idx_m],sep=""))
        mydata <- as.data.frame(cbind(M$age1,M$sex,M$site,M$covars.1.pca,M$y1))
        names(mydata)<-c("age","sex","site","covar1","covar2","covar3","phenotype");
     
        } else if (idx_type==2){
          
        # rate of change
        print(paste("rate of change data of ",modality[idx_m],sep=""))
        mydata <- as.data.frame(cbind(M$ageM,M$sex,M$site,M$covars.2.pca,M$rate))
        names(mydata)<-c("age","sex","site","covar1","covar2","covar3","phenotype");
        
      } 
      
      
      # Model training
        mdl <- gamlss(phenotype ~ fp(age,npoly=1) + sex + random(as.factor(site)) +pb(covar1) + pb(covar2) + pb(covar3), # formula for mu
                      sigma.fo=~fp(age,npoly=1),
                      nu.fo=~1,
                      nu.tau=~1,
                      family=BCT(),
                      data=mydata,
                      control = gamlss.control(n.cyc = 100))
        
     
      if (mdl$converged) {
        print("Model Converged : YES")
        converged<-1
      } else {
        print("Model Cnverged: NO")
        converged<-0
      }
      
      
      # check model fitting in the entire sample
      gvd1<-deviance(mdl,what="G") # validation global deviance
      iterations<-mdl$iter # number of iterations taken to converge
      

        params<-predictAll(mdl,data=mydata, newdata=mydata, 
                           output='matrix',type="response",
                           y.value="median",what=c("mu", "sigma", "nu", "tau"))
      
      
      
      ## Compute z-score ##
      # https://rdrr.io/cran/gamlss.dist/man/BCt.html
      quantiles <- pBCT(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]); #for BCT  
      z_randomized<-qnorm(quantiles, mean = 0, sd =1);
      ks_out<-ks.test(z_randomized,"pnorm")
      path<-"output_normative_modeling";
      mkdir(path)
      outname2<-paste("output_normative_modeling/fitstats_GAMLSS_",modality[idx_m],"_",type[idx_type],"_WHOLE_SAMPLE.mat",sep="")
      print(paste("KS test, p=",ks_out$p.value,sep=""))
      writeMat(outname2,  AIC1=mdl$aic, gvd1=gvd1,  BIC=mdl$sbc, zscores=z_randomized, ks_p=ks_out$p.value,ks_d=ks_out$statistic, converged=converged, P_deviance=mdl$P.deviance)
      
      ### predict hypothetical data ###
      for (gender in 0:1) {
        site1_unique<-unique(M$site)
        for (idx_site in 1:length(site1_unique)){ 
          print(paste("gender",gender,"site",idx_site,sep=" "))
          min_age<-min(data$age)
          max_age<-max(data$age)
          age_test<-linspace(min_age,max_age,1000)
          site_test<-matrix(data=site,nrow=1000,ncol=1)
          sex_test<-matrix(data=gender,nrow=1000,ncol=1)
          covar1_test<-matrix(data=mean(mydata$covar1),nrow=1000,ncol=1)
          covar2_test<-matrix(data=mean(mydata$covar2),nrow=1000,ncol=1)
          covar3_test<-matrix(data=mean(mydata$covar3),nrow=1000,ncol=1)
          time_test<-matrix(data=0,nrow=1000,ncol=1)
          Subject_test<-matrix(data=500,nrow=1000,ncol=1)
          y_test<-matrix(data=mean(mydata$phenotype),nrow=1000,ncol=1) # for the sake of the data structure required for the function
          # but not contribute to the prediction
          
        
            data_test<-as.data.frame(cbind(age_test,sex_test,site_test,covar1_test,covar2_test,covar3_test,y_test))
            names(data_test)<-c("age","sex","site","covar1","covar2","covar3","phenotype")
            params<-predictAll(mdl,data=mydata,
                               newdata=data_test,output='matrix',type="response",
                               y.value="median",what=c("mu", "sigma", "nu","tau"))
            
          
          #main percentiles
          plot(phenotype~age, data=mydata,  col="lightgray", pch=10)
          predictions_quantiles<-matrix(data=0,ncol=5,nrow=1000)
          quantiles <- pnorm(c(-2:2))
          for (i in 1:length(quantiles)){
            Qua <- getQuantile(mdl, quantile=quantiles[i],term="age",fixed.at=list(sex=gender,site=site))
            out<-curve(Qua, min_age, max_age,  lwd=2, lty=1, add=T,col="red",n = 1000) 
            predictions_quantiles[,i]=as.vector(out$y)
          }
          
          
          # save to mat files
          outname<-paste("output_normative_modeling/GAMLSS_",modality[idx_m],"_",type[idx_type],"_predicted_sex",gender,"_site",idx_site,"_WHOLE_SAMPLE.mat",sep="")
          writeMat(outname, predictions_quantiles=as.matrix(predictions_quantiles), age=age_test)
          
        }
      }
    }
  }

