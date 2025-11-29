# A very lazy way of tracing out the QRE -- just use the Stan file to do it!

library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

model<-"code_SFEMQRE/SFEMQRE_Trace.stan" |>
  stan_model()

strgList<-list(
  AD = list(start = 0,
            actions = c(0,0),
            transitions = rbind(c(0,0),c(0,0))
  ),
  AC = list(start = 1,
            actions = c(1,1),
            transitions = rbind(c(1,1),c(1,1))
  ),
  Grim = list(start = 1,
              actions = c(0,1),
              transitions = rbind(c(0,0),c(0,1))
  ),
  TFT = list(start = 1,
             actions = c(0,1),
             transitions = rbind(c(0,1),c(0,1))
  ),
  STFT = list(start = 0,
              actions = c(0,1),
              transitions = rbind(c(0,1),c(0,1))
  )#,
  #WSLS =list(start = 1,
  #           actions = c(0,1),
  #           transitions = rbind(c(1,0),c(0,1))
  #) 
)

SL<-names(strgList)

numFollow<-array(-1,c(length(strgList),length(strgList),4,4))
nextState<-array(-1,c(length(strgList),length(strgList),4,4))
firstState<-matrix(-1,length(strgList),length(strgList))
states<-rbind(c(0,0),c(0,1),c(1,0),c(1,1))
outcomes<-states

for (SS1 in 1:length(strgList)) {
  
  for (SS2 in 1:length(strgList)) {
    S1<-strgList[[SS1]]
    S2<-strgList[[SS2]]
    ns<-matrix(NA,4,4)
    for (ss in 1:4) {
      
      t1<-S1$transition[states[ss,1]+1,]
      t2<-S2$transitions[states[ss,2]+1,]
      
      for (oo in 1:dim(outcomes)[1]) {
        a1<-S1$actions[1+states[ss,1]]
        a2<-S2$actions[1+states[ss,2]]
        numFollow[SS1,SS2,ss,oo]<-1*(a1==outcomes[oo,1])+1*(a2==outcomes[oo,2])
        a1<-outcomes[oo,1]
        a2<-outcomes[oo,2]
        ns[ss,oo]<-2*t1[a2+1]+t2[a1+1]+1
      }
    }
    nextState[SS1,SS2,,]<-ns
    firstState[SS1,SS2]<-2*S1$start+S2$start+1
  }
  
}

TRACE<-tibble()

gld<-"code_SFEMQRE/01_estimate_treatments.csv" |>
  read.csv()

for (epsilon in c(0.01,0.05,0.1,0.15,0.2)) {
  
  print(epsilon)

    dStan<-list(
      nstrg = length(strgList),
      
      nTreatment = dim(gld)[1],
      g = gld$g,
      l = gld$l,
      delta = gld$delta,
      
      epsilon=epsilon,
      
      numFollow = numFollow,
      nextState = nextState,
      firstState = firstState,
      
      # not my real priors, but chosen to get a good spread of lambda
      prior_lambda = c(log(30),2),
      
      lGridSize = 50,
      nc = 100,
      
      ftol = 1e-6
      
      
    )
    
    trace<-model |>
      sampling(
        data = dStan,seed=42,
        iter=110, warmup=10 
        # warmup doesn't matter because we're not doing inference here
      )
    
    lambda<-extract(trace)$lambda
    mix<-extract(trace)$MIX
    
    for (ss in 1:dStan$nstrg) {
      for (tt in 1:dStan$nTreatment) {
        TRACE<-TRACE |>
          rbind(
            tibble(
              strg = SL[ss],
              g = gld$g[tt],
              l = gld$l[tt],
              delta = gld$delta[tt],
              r = gld$r[tt],
              lambda = lambda,
              epsilon = dStan$epsilon,
              mix = mix[,tt,ss]
              
            )
          )
      }
    }
}

TRACE<-TRACE |>
  arrange(
    r,delta,epsilon,lambda
  )

TRACE |>
  write.csv("code_SFEMQRE/02_traceQRE.csv")




(
  ggplot(TRACE,aes(x=lambda,y=mix,color=strg,linetype=as.factor(epsilon)))
  +geom_path()
  +facet_grid(r~delta)
  +scale_x_continuous(trans="log10")
  +theme_bw()
)

