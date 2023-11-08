####################################################################################################################################################
####     Developped by Ben Aabdelkrim and co-workers for the software R. See :
####     "A lactation curve model with explicit representation of perturbations as a phenotyping tool for dairy livestock precision farming"
####################################################################################################################################################

#install.packages('nls.multstart')
library(nls.multstart)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
#library(nlstools) #this was commented out before
setwd("/Users/mihirjagtap/Documents/Dairy_Brain/PLM-perturbations")

#setwd("/Users/victor/Library/CloudStorage/GoogleDrive-vcabrera@wisc.edu/My Drive/_10_UW-MADISON/23Sabbatical/Technical/A.Gallo/PLM-Perturbation")
#Cow.data <- read.csv("Cow.data.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#head(Cow.data)


############################ The perturbation model function
# This function calculates the effect of perturbation for a given time ktps with the parameters k0, k1, k2 of the ith perturbation . 


perturb = function(tt,ktps,k0, k1, k2){
  if(k1==k2){return(rep(1, length(tt)));}
  is_enabled = as.numeric(tt>=ktps)
  delta_t = is_enabled*( tt-ktps); # Determination of the elapsed time since the beginning of the perturbation. 
  x = ((k0*k1)/(k1-k2));#First part of the PLM equation
  y = exp(-k2*delta_t)-exp(-k1*delta_t);#Second part of the PLM equation
  r =  (is_enabled*x*y);# The PLM equation
  r=1-r
  return(r);
}


############################ Fitting procedure
# This function calculates the rate loss milk based in the difference between the unperturbed data (New_wood) and the perturbed data (pred) results calculated using the PLM model.

Loss=function(New_wood,pred){
  ecart=sum(New_wood)-sum(pred)
  r=(ecart*100)/(sum(New_wood))
  return(r);
}

# This function represents the model of lactation proposed by Wood (1967) (a, b and c : wood Parameters, t: time)

wood = function(t, a,b,c){ 
  (a*(t^b))*exp((-c)*t)
}

# PLM fitting function for NMAX perturbations
compute_plm=function(Cow.data, NMAX) {
  TIME = Cow.data$time
  # Creating of dataframe to save the news parameters
  PARAM = as.data.frame(matrix(nrow=NMAX,ncol=11)); 
  colnames(PARAM)=c("ID","N", "ktps","k0","k1","k2","a","b","c","Loss","Parity")
  data_fit = as.data.frame(matrix(nrow=length(TIME),ncol=0));
  PARAM[1,"ID"]=as.character(unique(Cow.data$ID));
  PARAM[1,"Parity"]=unique(Cow.data$parite)
  PARAM[1,"ktps"]=0;PARAM[1,"k0"]=0;PARAM[1,"k1"]=0;PARAM[1,"k2"]=0;
  # Creating a list object to save the results of fitting
  MODELS = list();
  supp_errors="Y";
  # fitting procedure for N perturbations
  for(i in 1:NMAX){
    MODELS[[i]]=list();
    # This function allows to apply "PLM" on N perturbations affecting the lactation curve
    PROD_perturb=function(time, N){
      r = rep(1,length(time))
      if(N==1){return (r);}
      for(j in 1:(N-1)){
        r=r*(perturb(time, PARAM[j,"ktps"],  PARAM[j,"k0"], PARAM[j,"k1"], PARAM[j,"k2"]));
      }
      return(r);
    }
    
    # lower boundaries for the start parameters of Wood model
    
    a0_inf = 0;
    b0_inf = 0;
    c0_inf = 0;
    # upper boundaries for the start parameters of Wood model
    a0_sup = 100;
    b0_sup = 1;
    c0_sup = 1;
    # Wood parameter estimates a, b, and c to characterize the unperturbed lactation curve  
    opt_wood=function(a,b,c){
      w = wood(TIME, a,b,c);
      p =  PROD_perturb(TIME,i);
      return (w*p);
    }     
    mod_wood<- nls_multstart(
      MY_obs ~ opt_wood(a,b,c),
      data =data.frame(Cow.data),
      iter =100000,
      start_lower = c(a=a0_inf , b=b0_inf, c=c0_inf),
      start_upper = c(a=a0_sup , b=b0_sup, c=c0_sup),
      supp_errors = supp_errors,
      lower       = c(a=a0_inf , b=b0_inf, c=c0_inf),
      upper       = c(a=a0_sup , b=b0_sup, c=c0_sup),
      convergence_count =100
    );
    
    # Save the results of Wood fitting
    MODELS[[i]][["wood"]]=mod_wood
    mod_wood_summary = summary(mod_wood)
    PARAM[i,"ID"]=as.character(unique(Cow.data$ID))
    PARAM[i,"Parity"]=unique(Cow.data$parite)
    PARAM[i,"N"]=i
    PARAM[i,"a"] = mod_wood_summary$coefficients[1]
    PARAM[i,"b"] = mod_wood_summary$coefficients[2]
    PARAM[i,"c"] = mod_wood_summary$coefficients[3]
    # lower boundaries for the start parameters of PLM model 
    ktps_inf = min(TIME,na.rm=T)+3
    k0_inf = 0
    k1_inf = 0   
    k2_inf = 0
    # upper boundaries for the start parameters of PLM model
    ktps_sup = max(TIME,na.rm=T)  
    k0_sup = 1
    k1_sup = 10
    k2_sup = 10
    
    # PLM parameter estimates a, b, and c to characterize the perturbed lactation curve
    old =  PROD_perturb(TIME,i)*wood(TIME, PARAM[i,"a"], PARAM[i,"b"], PARAM[i,"c"]);
    opt_perturb=function(ktps,k0,k1,k2){
      z = perturb(TIME,ktps,k0,k1,k2);
      return(old*z);
    }
    mod_perturb<- nls_multstart(
      MY_obs ~ opt_perturb(ktps,k0,k1,k2),
      data = Cow.data,
      iter = 100000,
      start_lower = c(ktps = ktps_inf, k0 = k0_inf, k1=k1_inf, k2=k2_inf),
      start_upper = c(ktps = ktps_sup, k0 = k0_sup, k1=k1_sup, k2=k2_sup),
      supp_errors = supp_errors,
      lower       = c(ktps = ktps_inf, k0 = k0_inf, k1=k1_inf, k2=k2_inf),
      upper       = c(ktps = ktps_sup, k0 = k0_sup, k1=k1_sup, k2=k2_sup),
      convergence_count =100
    );
    
    # Save the results of PLM fitting
    MODELS[[i]][["perturb"]]=mod_perturb
    mod_perturb_summary = summary(mod_perturb)
    PARAM[i,"ktps"] = mod_perturb_summary$coefficients[1]
    PARAM[i,"k0"]   = mod_perturb_summary$coefficients[2]
    PARAM[i,"k1"]   = mod_perturb_summary$coefficients[3]
    PARAM[i,"k2"]   = mod_perturb_summary$coefficients[4]
    pred=predict(mod_perturb)
    MODELS[[i]][["predict"]]=pred
    MODELS[[i]][["time"]]=time
    New_wood=wood(TIME,PARAM[i,"a"], PARAM[i,"b"], PARAM[i,"c"])
    PARAM[i,"Loss"]=Loss(New_wood,pred)
    MODELS[[i]][["New_wood"]]=New_wood
  }
  
  # Predict Perturbed curve with higher sample rate
  TIME2 = seq(min(TIME),max(TIME),by=0.5);
  pred_final =  PROD_perturb(TIME2,i)*wood(TIME2, PARAM[i,"a"], PARAM[i,"b"], PARAM[i,"c"]);
  MODELS[[i]][["predict2"]]=pred_final
  MODELS[[i]][["time2"]]=TIME2
  data_fit["ID"]=as.character(Cow.data$ID)
  data_fit["Parity"]=Cow.data$parite
  data_fit["time"]=TIME
  wn<-wood(TIME,PARAM[i,"a"], PARAM[i,"b"], PARAM[i,"c"])
  data_fit["Wood_new"]=wn
  data_fit["fit"]=MODELS[[i]][["predict"]]
  data_fit["Prpert"]=(wn-MODELS[[i]][["predict"]])/(wn)
  # plot
  plot(TIME,Cow.data$MY_obs,main =paste(titre_graph,", Np:",i,", Pr:",round(PARAM[i,"Loss"]),"%"), xlab= "Time (days)",
       col="blue", ylim=c(0,max(Cow.data$MY_obs)),ylab = "Milk yield (kg)",pch=21,cex=0.5, axes = FALSE)
  axis(1)
  axis(2)
  lines(TIME,MODELS[[i]][["predict"]],col="red",lwd=2)
  lines(Cow.data$time,wood(TIME,PARAM[i,"a"], PARAM[i,"b"], PARAM[i,"c"]),col="black",lwd=2)
  legend("topright",box.lty=0, legend=c("data","Unperturbed","Perturbed"), cex=0.8,lty=c(NA,1,1),pch=c(21,NA,NA), lwd=c(1,1.5,1.5),col=c("blue","black","red"))
  return(list(param=PARAM, data_fit= data_fit, models=MODELS))
}


