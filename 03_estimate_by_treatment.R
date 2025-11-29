library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(bridgesampling)

RERUN<-TRUE

Dstrg<-"code_SFEMQRE/00_loaddata.csv" |>
  read.csv() |>
  group_by(treatmentID,r,g,l,delta,id) |>
  summarize(
    strgAD = sum(strgAD==coop),
    strgAC = sum(strgAC==coop),
    strgGrim = sum(strgGrim==coop),
    strgTFT = sum(strgTFT==coop),
    strgSTFT = sum(strgSTFT==coop),
    count = n()
  ) |>
  ungroup() |>
  mutate(
    TID = treatmentID |> as.factor() |> as.integer()
  )

gld<-Dstrg |>
  group_by(TID) |>
  summarize(
    r = mean(r),
    g = mean(g),
    l = mean(l),
    delta = mean(delta)
  )

gld |>
  write.csv("code_SFEMQRE/01_estimate_treatments.csv")

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

ESTIMATES<-tibble()

POSTPROBS<-tibble()

for (tt in 1:6) {
  
  print(gld[tt,])
  
  dstrg<-Dstrg |> filter(TID==tt)
  
  dStan<-list(
    N = dim(dstrg)[1],
    nstrg = length(strgList),
    follow = dstrg |> select(contains("strg")),
    treatmentID = rep(1,dim(dstrg)[1]),
    nTreatment = dstrg$TID |> unique() |> length(),
    
    g = gld$g[tt] |> array(),
    l = gld$l[tt] |> array(),
    delta = gld$delta[tt] |> array(),
    
    count = dstrg$count,
    
    numFollow = numFollow,
    nextState = nextState,
    firstState = firstState,
    
    prior_lambda = c(1.15,1.76),
    prior_MIX = rep(1,5),
    
    lGridSize = 50,
    nc = 100,
    
    ftol = 1e-6,
    UseData = rep(1,dim(dstrg)[1])
    
    
  )
  
  SFEMQRE<-"code_SFEMQRE/SFEMQRE_pooled.stan" |>
    stan_model()
  
  Fit<-SFEMQRE |>
    sampling(data=dStan,
             #chains = 1,iter = 10,
             seed=42)
  
  bs_SFEMQRE<-bridge_sampler(Fit)
  
  l<-loo(Fit)
  
  ESTIMATES<-ESTIMATES |>
    rbind(
      summary(Fit)$summary |>
        data.frame() |>
        rownames_to_column(var = "par") |>
        mutate(
          model = "QRE-SFEM",
          r = gld$r[tt],
          delta= gld$delta[tt],
          elpd_loo = l$estimates[1,1],
          elpd_loo_se = l$estimates[1,2]
        )
    )
  
  SFEM<-"code_SFEMQRE/SFEM_pooled.stan" |>
    stan_model()
  
  Fit<-SFEM |>
    sampling(data=dStan,
             #chains = 1,iter = 10,
             seed=42)
  
  bs_SFEM<-bridge_sampler(Fit)
  
  l<-loo(Fit)
  
  ESTIMATES<-ESTIMATES |>
    rbind(
      summary(Fit)$summary |>
        data.frame() |>
        rownames_to_column(var = "par") |>
        mutate(
          model = "SFEM",
          r = gld$r[tt],
          delta= gld$delta[tt],
          elpd_loo = l$estimates[1,1],
          elpd_loo_se = l$estimates[1,2]
        )
    )
  
  postprobs<-post_prob(bs_SFEMQRE,bs_SFEM)
  
  POSTPROBS<-rbind(
    POSTPROBS,
    tibble(
      r = gld$r[tt],
      delta= gld$delta[tt],
      postprob.SFEMQRE=postprobs[1]
    )
  )
  
}

ESTIMATES |>
  write.csv("code_SFEMQRE/03_estimate_by_treatment.csv")

POSTPROBS |>
  write.csv("code_SFEMQRE/03_estimate_by_treatment_postprobs.csv")