## BOCA-DA design simulation

## borrow info from marker missing samples

install.packages(R2jags)
library(rjags)

BOCA_m2=function(s,t,p_c,p0,p1,phi,pil,pis,pif,sm,mis){
  
  
  decision=NULL ## final decision
  D_n_p = NULL ## neg marker is promising
  D_p_p = NULL ## pos marker is promising
  ss_neg = NULL ## final sample size for neg marker in trial
  ss_pos = NULL ## final sample size for pos marker in trial
  ss_trial = NULL ## final sample size in trial of one simulation = ss_neg + ss_pos
  ss_neg_1 = NULL ## final sample size for neg marker 
  ss_pos_1 = NULL ## final sample size for pos marker 
  ss_all = NULL ## final total sample size = ss_neg_1+ss_pos_1
  p1_post = NULL ## posterior probability for pos marker
  p0_post = NULL ## posterior probability for neg marker
  p=p0*(1-phi)+p1*phi ## response rate for overall population
  
  BOCAMCMC<-"model{
        r1~dbinom(p1, n_pos)
        r0~dbinom(p0, n_neg)
       nm1~dbinom(phi,n_pos+n_neg+m0+m1)
       r2~dbinom(p,n_mis)
       p <-p1*phi+p0*(1-phi)
p0~dbeta(0.5,0.5)
p1~dbeta(0.5,0.5)
        }"
  
  for (l in 1: sm){
    n_mis=rbinom(1,s,mis) ## missing samples
    n_pos=rbinom(1,s-n_mis,phi) ## pos samples
    n_neg=s-n_pos-n_mis ## neg samples
    r0=rbinom(1,n_neg,p0) ## 1st stage response samples (-)
    r1=rbinom(1,n_pos,p1) ## 1st stage response samples (+)
    r2=rbinom(1,n_mis,p)
    m0=0 ## marker (-) samples not in trial
    m1=0 ## marker (+) samples not in trial
    m2=0 ## marker missing samples not in trial
    earlystop=0
    for(i in 1: (t-1)){
      model <- jags.model(textConnection(BOCAMCMC), 
                          data = list(n_neg=n_neg,n_pos=n_pos,n_mis=n_mis,r0=r0,r1=r1,r2=r2, nm1=n_pos+m1,m0=m0, m1=m1,phi=phi ),quiet=T)
      update(model, 2000, progress.bar="none") ## Markov Chain for burn-in period of 2000 iterations
      re <- jags.samples(model, n.iter=5000,variable.names=c("p0", "p1"), progress.bar="none")
      p1_array <- re$p1[1:5000]
      p0_array <- re$p0[1:5000]
      
      if( mean(p1_array> p_c)<pil){
        earlystop=1
        ss_pos[l] = n_pos
        ss_neg[l] = n_neg
        ss_trial[l] = n_pos + n_neg + n_mis
        ss_pos_1[l] = n_pos
        ss_neg_1[l] = n_neg
        ss_all[l] = n_pos + n_neg + n_mis
        decision[l] = 0 ## early stop for futility before enriched stage
        D_p_p[l] = 0
        D_n_p[l] = 0
        p1_post[l] = mean(p1_array)
        p0_post[l] = mean(p0_array)
        break
      }
      else if (mean(p0_array> p_c)>pis){
        earlystop=1
        ss_pos[l] = n_pos
        ss_neg[l] = n_neg
        ss_trial[l] = n_pos + n_neg + n_mis
        ss_pos_1[l] = n_pos
        ss_neg_1[l] = n_neg
        ss_all[l] = n_pos + n_neg + n_mis
        decision[l] = 00 ## early stop for superiority before enriched stage
        D_p_p[l] = 1
        D_n_p[l] = 1
        p1_post[l] = mean(p1_array)
        p0_post[l] = mean(p0_array)
        break
        
      } else if (mean(p1_array> p_c)>pis & mean(p0_array> p_c)<pil  ){
        earlystop=1
        ss_pos[l] = n_pos
        ss_neg[l] = n_neg
        ss_trial[l] = n_pos + n_neg + n_mis
        ss_pos_1[l] = n_pos
        ss_neg_1[l] = n_neg
        ss_all[l] = n_pos + n_neg + n_mis
        decision[l] = 100 ## early stop for superiority for pos and futility for neg
        D_p_p[l] = 1
        D_n_p[l] = 0
        p1_post[l] = mean(p1_array)
        p0_post[l] = mean(p0_array)
        break
        
      }
      else if (mean(p1_array> p_c)>pis){
        m1 = 0 ## marker (+) samples not in trial
        m2 = 0 ## marker (missing) samples not in trial
        ## enriched stage
        for (k in i:(t-1)){
          n_mis1=rbinom(1,s,mis) ## missing samples
          n_neg1=rbinom(1,s-n_mis1,1-phi)
          n_pos1=s-n_neg1-n_mis1 ## neg samples
          n_neg = n_neg+ n_neg1## marker (+) samples
          m1 = m1 + n_pos1 ## enriched stage marker (+) samples  not in trial
          m2 = m2 + n_mis1 ## enriched stage missing marker samples  not in trial
          r0=r0+rbinom(1,n_neg1,p0) ## enriched stage response samples (-)
          
          model <- jags.model(textConnection(BOCAMCMC), 
                              data = list(n_neg=n_neg,n_pos=n_pos,n_mis=n_mis,r0=r0,r1=r1,r2=r2,nm1=n_pos+m1,m0=m0,m1=m1,phi=phi  ),quiet=T)
          update(model, 2000, progress.bar="none") ## Markov Chain for burn-in period of 2000 iterations
          re <- jags.samples(model, n.iter=5000,variable.names=c("p0"), progress.bar="none")
          p0_array <- re$p0[1:5000]
          
          if(mean(p0_array> p_c)<pil ){
            earlystop=1
            decision[l]=11 ## early stop for (-) futility during enriched stage
            D_p_p[l] = 1
            D_n_p[l] = 0
            ss_pos[l] = n_pos
            ss_neg[l] = n_neg
            ss_trial[l] = n_pos + n_neg + n_mis
            ss_pos_1[l] = n_pos + m1
            ss_neg_1[l] = n_neg
            ss_all[l] = n_pos + n_neg + n_mis + m1 + m2
            p1_post[l] = mean(p1_array)
            p0_post[l] = mean(p0_array)
            break
          }
          else if (mean(p0_array> p_c) > pis ){
            earlystop=1
            decision[l]=22 ## early stop for (-) superiority during enriched stage
            D_p_p[l] = 1
            D_n_p[l] = 1
            ss_pos[l] = n_pos
            ss_neg[l] = n_neg
            ss_trial[l] = n_pos + n_neg + n_mis
            ss_pos_1[l] = n_pos + m1
            ss_neg_1[l] = n_neg
            ss_all[l] = n_pos + n_neg  + n_mis + m1 + m2
            p1_post[l] = mean(p1_array)
            p0_post[l] = mean(p0_array)
            break
          }
          
          else {
            earlystop=0
            
            # continue enriched stage until maximizing sample size
          }
        }
        
        # reach maximum sample size in enriched stage
        if(earlystop==0 & mean(p0_array> p_c)<pif){ # final analysis
          earlystop=1
          decision[l]=2 # pos promising but neg are not promising
          D_p_p[l] = 1
          D_n_p[l] = 0
          ss_pos[l] = n_pos
          ss_neg[l] = n_neg
          ss_pos_1[l] = n_pos + m1
          ss_neg_1[l] = n_neg
          ss_trial[l] = n_pos + n_neg +n_mis
          ss_all[l] = n_pos + n_neg  + n_mis +  m1  + m2 # should be equal to total sample size
          p1_post[l] = mean(p1_array)
          p0_post[l] = mean(p0_array)
          break
        }
        
        else if (earlystop==0){
          earlystop=1
          decision[l]=7 # both pos and neg are promising
          D_p_p[l] = 1
          D_n_p[l] = 1
          ss_pos[l] = n_pos
          ss_neg[l] = n_neg
          ss_trial[l] = n_pos + n_neg + n_mis
          ss_pos_1[l] = n_pos + m1
          ss_neg_1[l] = n_neg
          ss_all[l] = n_pos + n_neg   + n_mis +  m1  + m2# should be equal to total sample size
          p1_post[l] = mean(p1_array)
          p0_post[l] = mean(p0_array)
          break
        }
        
      } ## end of enriched stage
      
      
      else if(  mean(p0_array> p_c)<pil ){
        m0 = 0 ## marker (-) samples not in trial
        m2 = 0 ## marker missing samples not in trial
        ## enriched stage
        for (k in i:(t-1)){
          n_mis1=rbinom(1,s,mis) ## missing samples
          n_pos1 = rbinom(1,s-n_mis1,phi)
          n_pos = n_pos+ n_pos1## marker (+) samples
          m0 = m0 + s - n_mis1 - n_pos1 ## enriched stage marker (-) samples that not in trial
          m2 = m2 + n_mis1
          r1=r1+rbinom(1,n_pos1,p1) ## enriched stage response samples (+)
          
          model <- jags.model(textConnection(BOCAMCMC), 
                              data = list(n_neg=n_neg,n_pos=n_pos,n_mis=n_mis,r0=r0,r1=r1,r2=r2,nm1=n_pos+m1,m0=m0,m1=m1,phi=phi  ),quiet=T)
          update(model, 2000, progress.bar="none") ## Markov Chain for burn-in period of 2000 iterations
          re <- jags.samples(model, n.iter=5000,variable.names=c("p1"), progress.bar="none")
          p1_array <- re$p1[1:5000]
          
          if( mean(p1_array> p_c) > pis ){
            earlystop=1
            decision[l]=222 ## early stop for (+) superiority  during enriched stage
            D_p_p[l] = 1
            D_n_p[l] = 0
            ss_pos[l] = n_pos
            ss_neg[l] = n_neg
            ss_trial[l] = n_pos + n_neg + n_mis
            ss_pos_1[l] = n_pos
            ss_neg_1[l] = n_neg + m0
            ss_all[l] = n_pos + n_neg  + n_mis + m0 + m2
            p1_post[l] = mean(p1_array)
            p0_post[l] = mean(p0_array)
            break
          }
          else if( mean(p1_array> p_c)<pil ){
            earlystop=1
            decision[l]=111 ## early stop for (+) futility during enriched stage
            D_p_p[l] = 0
            D_n_p[l] = 0
            ss_pos[l] = n_pos
            ss_neg[l] = n_neg
            ss_trial[l] = n_pos + n_neg + n_mis
            ss_pos_1[l] = n_pos
            ss_neg_1[l] = n_neg + m0
            ss_all[l] = n_pos + n_neg  + n_mis + m0 + m2
            p1_post[l] = mean(p1_array)
            p0_post[l] = mean(p0_array)
            break
          }
          
          else {
            earlystop=0
            
            # continue enriched stage until maximizing sample size
          }
        }
        
        # reach maximum sample size in enriched stage
        if(earlystop==0 &  mean(p1_array> p_c)<pif){ # final analysis
          earlystop=1
          decision[l]=4 # neither pos nor neg is promising
          D_p_p[l] = 0
          D_n_p[l] = 0
          ss_pos[l] = n_pos
          ss_neg[l] = n_neg
          ss_pos_1[l] = n_pos
          ss_neg_1[l] = n_neg + m0
          ss_trial[l] = n_pos + n_neg  + n_mis
          ss_all[l] = n_pos + n_neg  + n_mis + m0  + m2# should be equal to total sample size
          p1_post[l] = mean(p1_array)
          p0_post[l] = mean(p0_array)
          break
        }
        else if (earlystop==0){
          earlystop=1
          decision[l]=2 # pos promising but neg are not promising
          D_p_p[l] = 1
          D_n_p[l] = 0
          ss_pos[l] = n_pos
          ss_neg[l] = n_neg
          ss_trial[l] = n_pos + n_neg + n_mis
          ss_pos_1[l] = n_pos
          ss_neg_1[l] = n_neg + m0
          ss_all[l] = n_pos + n_neg  + n_mis + m0  + m2# should be equal to total sample size
          p1_post[l] = mean(p1_array)
          p0_post[l] = mean(p0_array)
          break
        }
        
      } ## end of enriched stage
      
      else { ##  neither (+) nor (-) has evidence of futility or superiority at this interim
        earlystop=0
        n_mis1=rbinom(1,s,mis) ## missing samples
        n_mis=n_mis+n_mis1
        n_pos_1 = rbinom(1,s-n_mis1,phi) 
        n_pos=n_pos+n_pos_1 ## 1st stage marker (+) samples
        # n_neg=(i+1)*s-n_pos ## 1st stage marker (-) samples
        n_neg=n_neg+s-n_mis1-n_pos_1
        r1=r1+rbinom(1,n_pos_1,p1) ## 1st stage response samples (+)
        r0=r0+rbinom(1,s-n_mis1-n_pos_1,p0) ## 1st stage response samples (-)
        r2=r2+rbinom(1,n_mis1,p) ## 1st stage response samples (missing)
      } 
    } # end of interim
    
    if (earlystop==0)   { 
      # reach maximum sample size at end of interims
      if( mean(p1_array> p_c)<pif){ # final analysis
        decision[l]=5 # neither pos nor neg is promising
        D_p_p[l] = 0
        D_n_p[l] = 0
        ss_pos[l] = n_pos
        ss_neg[l] = n_neg
        ss_trial[l] = n_pos + n_neg + n_mis
        ss_pos_1[l] = n_pos
        ss_neg_1[l] = n_neg
        ss_all[l] = n_pos + n_neg + n_mis # should be equal to total sample size s*t
        p1_post[l] = mean(p1_array)
        p0_post[l] = mean(p0_array)
      }
      else if ( mean(p0_array> p_c)<pif) {
        decision[l]=6 # pos promising but neg are not promising
        D_p_p[l] = 1
        D_n_p[l] = 0
        ss_pos[l] = n_pos
        ss_neg[l] = n_neg
        ss_trial[l] = n_pos + n_neg + n_mis
        ss_pos_1[l] = n_pos
        ss_neg_1[l] = n_neg
        ss_all[l] = n_pos + n_neg+n_mis # should be equal to total sample size s*t
        p1_post[l] = mean(p1_array)
        p0_post[l] = mean(p0_array)
      }
      else {
        decision[l]=7 # both pos and neg are promising
        D_p_p[l] = 1
        D_n_p[l] = 1
        ss_pos[l] = n_pos
        ss_neg[l] = n_neg
        ss_trial[l] = n_pos + n_neg + n_mis
        ss_pos_1[l] = n_pos
        ss_neg_1[l] = n_neg
        ss_all[l] = n_pos + n_neg+n_mis # should be equal to total sample size s*t
        p1_post[l] = mean(p1_array)
        p0_post[l] = mean(p0_array)
      }
    } 
  } # end of simulation
  
  return(c("BOCA-DA",round(mean(D_p_p)*100,digits=2) ,round(mean(D_n_p)*100,digits=2) ,
           round(mean(ss_trial),digits=2),
           round(mean(ss_all),digits=2),
           round(mean(p1_post),digits=2),
           round(mean(p0_post),digits=2),
           p_c, p0,p1, phi,pil,pis,pif,mis
  ))
}


#BOCA-DA simulation examples
pp_c=0.2
pp0=0.05
pphi=0.25
mis=0.3

set.seed (12345)
BOCA_m2p1 <- BOCA_m2(s=10,t=8,p_c=pp_c,p0=pp0,p1=pp0,phi=pphi,pil=0.1,pis=0.99,pif=0.94,sm=5000,mis=mis)
set.seed (12345)
BOCA_m2p2 <- BOCA_m2(s=10,t=8,p_c=pp_c,p0=pp0,p1=pp0+0.1,phi=pphi,pil=0.1,pis=0.99,pif=0.94,sm=5000,mis=mis)
set.seed (12345)
BOCA_m2p3 <- BOCA_m2(s=10,t=8,p_c=pp_c,p0=pp0,p1=pp0+0.35,phi=pphi,pil=0.1,pis=0.99,pif=0.94,sm=5000,mis=mis)
set.seed (12345)
BOCA_m2p4 <- BOCA_m2(s=10,t=8,p_c=pp_c,p0=pp0+0.1,p1=pp0+0.1,phi=pphi,pil=0.1,pis=0.99,pif=0.94,sm=5000,mis=mis)
set.seed (12345)
BOCA_m2p5 <- BOCA_m2(s=10,t=8,p_c=pp_c,p0=pp0+0.1,p1=pp0+0.25,phi=pphi,pil=0.1,pis=0.99,pif=0.94,sm=5000,mis=mis)
set.seed (12345)
BOCA_m2p6 <- BOCA_m2(s=10,t=8,p_c=pp_c,p0=pp0+0.25,p1=pp0+0.35,phi=pphi,pil=0.1,pis=0.99,pif=0.94,sm=5000,mis=mis)
set.seed (12345)
BOCA_m2p7 <- BOCA_m2(s=10,t=8,p_c=pp_c,p0=pp0+0.35,p1=pp0+0.35,phi=pphi,pil=0.1,pis=0.99,pif=0.94,sm=5000,mis=mis)


BOCADA_SIM <- as.data.frame(rbind(BOCA_m2p1,
                                  BOCA_m2p2,
                                  BOCA_m2p3,
                                  BOCA_m2p4,
                                  BOCA_m2p5,
                                  BOCA_m2p6,
                                  BOCA_m2p7))

colnames(BOCADA_SIM) <- c("Design","Success in Marker (+) Group (%)","Success in Marker (-) Group (%)",
                          "Mean n","Mean Screened","posteria prob (+)","posteria prob (-)",
                          "Target","P (-)","P (+)","(+) Prevalence","P (Fut)","P (Sup)","P (Fin)","P (Missing)")

