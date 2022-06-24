# Develop the decision rule for BOCA-DA design 
# borrowing information from missing samples
install.packages(R2jags)
install.packages(bayestestR)
library(rjags)
library(bayestestR)
BOCA_rule=function(ncohort,n_pos,n_neg, n_mis,r0, r1,r2,pil,pis,pif,p_c){
  
  
  decision=NULL ## final decision
  D_n_p = NULL ## neg marker is promising
  D_p_p = NULL ## pos marker is promising
  p1_post = NULL ## posterior probability for pos marker
  p0_post = NULL ## posterior probability for neg marker
  enrich = NULL ##if the enter enriched stage
  earlystop=0
  t <- length(n_pos)
  
  BOCAMCMC<-"model{
        r1~dbinom(p1, n_pos)
        r0~dbinom(p0, n_neg)
       r2~dbinom(p,n_mis)
       p <-p1*(n_pos/(n_pos+n_neg))+p0*(1-n_pos/(n_pos+n_neg))
p0~dbeta(0.5,0.5)
p1~dbeta(0.5,0.5)
        }"
  
  for(i in 1: t){
    model <- jags.model(textConnection(BOCAMCMC), 
                        data = list(n_neg=sum(n_neg[1:i]),n_pos=sum(n_pos[1:i]),
                                    n_mis=sum(n_mis[1:i]),r0=sum(r0[1:i]),
                                    r1=sum(r1[1:i]),r2=sum(r2[1:i])),quiet=T)
    update(model, 20, progress.bar="none") ## Markov Chain for burn-in period of 2000 iterations
    re <- jags.samples(model, n.iter=50,variable.names=c("p0", "p1"), progress.bar="none")
    p1_array <- re$p1[1:50]
    p0_array <- re$p0[1:50]
    
    if( mean(p1_array> p_c)<pil){
      earlystop=1
      
      decision = "Stop the trial for futility: Marker(+):N; Marker(-):N" ## early stop for futility before enriched stage
      D_p_p = 0
      D_n_p = 0
      
     
      enrich='Not in enriched stage'
      break
    }
    else if (mean(p0_array> p_c)>pis){
      earlystop=1
      
      decision = "Stop the trial for superiority: Marker(+):Y; Marker(-):Y" ## early stop for superiority before enriched stage
      D_p_p = 1
      D_n_p = 1
      
     
      enrich='Not in enriched stage'
      break
      
    } else if (mean(p1_array> p_c)>pis & mean(p0_array> p_c)<pil ){
      earlystop=1
      decision = "Stop the trial:Marker(+):Y; Marker(-):N" ## early stop for superiority for pos and futility for neg
      D_p_p = 1
      D_n_p = 0
      
     
      enrich='Not in enriched stage'
      break
      
    }
    else if (mean(p1_array> p_c)>pis){
      earlystop=0
      enrich = '## enriched stage, only enroll (-) in trial' 
      ## enriched stage
      for (k in i:t){
        if (k==t & t<ncohort){decision = "only enroll Marker(-) in trial" 
        break}
        else if (k<t){
          model <- jags.model(textConnection(BOCAMCMC), 
                              data = list(n_neg=sum(n_neg[1:i]),n_pos=sum(n_pos[1:i]),
                                          n_mis=sum(n_mis[1:i]),r0=sum(r0[1:i]),
                                          r1=sum(r1[1:i]),r2=sum(r2[1:i])),quiet=T)
          update(model, 20, progress.bar="none") ## Markov Chain for burn-in period of 2000 iterations
          re <- jags.samples(model, n.iter=50,variable.names=c("p0", "p1"), progress.bar="none")
          p1_array <- re$p1[1:50]
          p0_array <- re$p0[1:50]
          
          if(mean(p0_array> p_c)<pil ){
            earlystop=1
            decision="Stop the trial for futility: Marker(+):N; Marker(-):N" ## early stop for (-) futility during enriched stage
            D_p_p = 1
            D_n_p = 0
            
            
            break
          }
          else if (mean(p0_array> p_c) > pis){
            earlystop=1
            decision="Stop the trial for superiority: Marker(+):Y; Marker(-):Y" ## early stop for (-) superiority during enriched stage
            D_p_p = 1
            D_n_p = 1
            
            
            break
          }
          
          else {
            earlystop=0
            # continue enriched stage until maximizing sample size
          }
        }
      }
      
      # reach maximum sample size in enriched stage
      if(t==ncohort & earlystop==0 & mean(p0_array> p_c)<pif){ # final analysis
        earlystop=11
        decision='Final analysis: Marker(+):Y; Marker(-):N' # pos promising but neg are not promising
        D_p_p = 1
        D_n_p = 0
        
        
        break
      }
      
      else if (t==ncohort & earlystop==0){
        earlystop=11
        decision='Final analysis: Marker(+):Y; Marker(-):Y' # both pos and neg are promising
        D_p_p = 1
        D_n_p = 1
        
        
        break
      }
    } ## end of enriched stage
    
    
    else if(   mean(p0_array> p_c)<pil){
      earlystop=0
      enrich = '## enriched stage, only enroll (+) in trial' 
      for (k in i:t){
        if (k==t & t<ncohort){decision = "only enroll Marker(+) in trial" 
        break}
        else if (k<t){
          model <- jags.model(textConnection(BOCAMCMC), 
                              data = list(n_neg=sum(n_neg[1:i]),n_pos=sum(n_pos[1:i]),
                                          n_mis=sum(n_mis[1:i]),r0=sum(r0[1:i]),
                                          r1=sum(r1[1:i]),r2=sum(r2[1:i])),quiet=T)
          update(model, 20, progress.bar="none") ## Markov Chain for burn-in period of 2000 iterations
          re <- jags.samples(model, n.iter=50,variable.names=c("p0", "p1"), progress.bar="none")
          p1_array <- re$p1[1:50]
          p0_array <- re$p0[1:50]
          
          if(  mean(p1_array> p_c) > pis ){
            earlystop=1
            decision="Stop the trial for superiority. Marker(+):Y; Marker(-):N" ## early stop for (+) superiority  during enriched stage
            D_p_p = 1
            D_n_p = 0
            
            
            break
          }
          else if( mean(p1_array> p_c)<pil  ){
            earlystop=1
            decision="Stop the trial for superiority. Marker(+):N; Marker(-):N" ## early stop for (+) futility during enriched stage
            D_p_p = 0
            D_n_p = 0
            
            
            break
          }
          
          else {
            earlystop=0
          }
          # continue enriched stage until maximizing sample size
        }
      }
      
      # reach maximum sample size in enriched stage
      if(t==ncohort & earlystop==0 &  mean(p1_array> p_c)<pif){ # final analysis
        earlystop=11
        decision='Final analysis: Marker(+):N; Marker(-):N' # neither pos nor neg is promising
        D_p_p = 0
        D_n_p = 0
        
        
        break
      }
      else if (t==ncohort & earlystop==0){
        earlystop=11
        decision='Final analysis: Marker(+):Y; Marker(-):N' # pos promising but neg are not promising
        D_p_p = 1
        D_n_p = 0
        
        
        break
      }
      
    } ## end of enriched stage
    
    else { ##  neither (+) nor (-) has evidence of futility or superiority at this interim
      earlystop=0
    }
    
  } # end of interim
  
  if (ncohort==t)   { 
    # reach maximum sample size at end of interims
    if( mean(p1_array> p_c)<pif){ # final analysis
      decision="final analysis: neither Marker(+) nor Marker(-) is promising" # neither pos nor neg is promising
      D_p_p = 0
      D_n_p = 0
      
      
      enrich='Not in enriched stage'
    }
    else if ( mean(p0_array> p_c)<pif) {
      decision="final analysis: Marker(+) is promising and Marker(-) is not" # pos promising but neg are not promising
      D_p_p = 1
      D_n_p = 0
      
      
      enrich='Not in enriched stage'
    }
    else {
      decision="final analysis: both Marker(+) and Marker(-) are promising" # both pos and neg are promising
      D_p_p = 1
      D_n_p = 1
      
      
      enrich='Not in enriched stage'
    }
  } 
  else if (is.null(decision)) {
    decision="enroll Marker(+) and Marker(-) in trial"
  }
  
  return(list("BOCA-DA decisions",decision,
              #round(quantile(p1_array, 0.25),digits=2),
              round(median(p1_array),digits=2),
              #round(quantile(p1_array, 0.75),digits=2),
              round(quantile(p1_array, probs = c(0.025, 0.975)),digits=2),
            #round(quantile(p0_array, 0.25),digits=2),
            round(median(p0_array),digits=2),
            #round(quantile(p0_array, 0.75),digits=2),
            round(quantile(p0_array, probs = c(0.025, 0.975)),digits=2),
            #ci(p0_array, method = "HDI"),
            round(mean(p1_array> p_c),digits=2),
            round(mean(p0_array> p_c),digits=2)
            )
         )
         
         }




# BOCA-DA design decision rule examples:

#ncohort:"# of cohorts",

#n_pos: "Marker(+) sample size",

#r1: "Marker(+) response",

#n_neg: "Marker(-) sample size",

#r0: "Marker(-) response", 

#n_mis: "Marker missing sample size",

#r2: "Marker(-) response", 

#pil:"P (Fut), cutoff of early stop for futility",

#pis: "P (Sup), cutoff of early stop for superiority",

#pif: "P (Final), cutoff of final analysis",

#p_c: "Target probability"

# example cohorts A:
BOCA_r1a <-BOCA_rule(ncohort=3,n_pos=c(20),n_neg=c(8),r1=c(10),r0=c(2),n_mis=c(2),r2=c(0),
                    pil=0.1,pis=0.99,pif=0.94,p_c=0.2)
BOCA_r1a
BOCA_r2a <-BOCA_rule(ncohort=3,n_pos=c(20,0),n_neg=c(8,8),r1=c(10,0),r0=c(2,1),n_mis=c(2,2),r2=c(0,1),
                     pil=0.1,pis=0.99,pif=0.94,p_c=0.2)
BOCA_r2a

BOCA_r3a <-BOCA_rule(ncohort=3,n_pos=c(20,0,0),n_neg=c(8,8,6),
                     r1=c(10,0,0),r0=c(2,1,1),n_mis=c(2,2,4),r2=c(0,1,2),
                     pil=0.1,pis=0.99,pif=0.94,p_c=0.2)
BOCA_r3a



BOCADA_rules <- as.data.frame(rbind(BOCA_r1a,
                                     BOCA_r2a,
                                     BOCA_r3a))


colnames(BOCADA_rules) <- c("Design","Decision",
                            "median of Marker(+) posterior response probability",
                            "95% credible intervals of Marker(+) posterior response probability",
                            "median of Marker(-) posterior response probability",
                            "95% credible intervals of Marker(-) posterior response probability",
                            "Power for Marker(+)",
                            "Power for Marker(-)")

BOCADA_rules

# example cohorts B:
BOCA_r1b <-BOCA_rule(ncohort=3,n_pos=c(20),n_neg=c(8),r1=c(6),r0=c(4),n_mis=c(2),r2=c(0),
                    pil=0.1,pis=0.99,pif=0.94,p_c=0.2)
BOCA_r1b
BOCA_r2b <-BOCA_rule(ncohort=3,n_pos=c(20,20),n_neg=c(8,8),r1=c(10,15),r0=c(4,4),n_mis=c(2,2),r2=c(0,1),
                     pil=0.1,pis=0.99,pif=0.94,p_c=0.2)
BOCA_r2b




BOCADA_rules <- as.data.frame(rbind(BOCA_r1b,
                                    BOCA_r2b))


colnames(BOCADA_rules) <- c("Design","Decision",
                            "median of Marker(+) posterior response probability",
                            "95% credible intervals of Marker(+) posterior response probability",
                            "median of Marker(-) posterior response probability",
                            "95% credible intervals of Marker(-) posterior response probability",
                            "Power for Marker(+)",
                            "Power for Marker(-)")

BOCADA_rules
