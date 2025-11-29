library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bridgesampling)

RERUN<-TRUE


# A list of all strategies
strg<-tibble()

for (ss in 0:(2^7-1)) {
  
  s<-intToBits(ss)[1:7]
  
  strg<-rbind(strg,
              tibble(
                x0  = s[1] |> as.numeric(),
                a0  = s[2] |> as.numeric(),
                x00 = s[3] |> as.numeric(),
                x01 = s[4] |> as.numeric(),
                a1  = s[5] |> as.numeric(),
                x10 = s[6] |> as.numeric(),
                x11 = s[7] |> as.numeric()
              )
  )
  
}

dim(strg)[1]



strg<-strg |>
  # remove everything with two D states or two C states
  filter(
    a0!=a1
  )  |>
  # remove everything with state0=1 and state1=0
  # At this point state 0 ==> defect and state 1 ==> cooperate
  filter(
    !(a0==1 & a1==0)
  ) |>
  # remove strategies that never get to state 1
  filter(
    !(x0 ==0 & x00==0 & x01==0 )
  ) |>
  #remove strategies that never get to state 0
  filter(
    !(x0 ==1 & x10==1 & x11==1 )
  ) |>
  
  #add back in AD and AC 
  rbind(
    tibble(
      x0  = 0,
      a0  = 0,
      x00 = 0,
      x01 = 0,
      a1  = 0,
      x10 = 0,
      x11 = 0 
    ),
    tibble(
      x0  = 1,
      a0  = 1,
      x00 = 1,
      x01 = 1,
      a1  = 1,
      x10 = 1,
      x11 = 1 
    )
  )


dim(strg[1])

# put these into list form

strgList<-list()

for (ss in 1:dim(strg)[1]) {
  
  strgList[[ss]]<-list(
    code=strg[ss,] |> matrix() |> as.vector() |> unlist(),
    start = strg$x0[ss],
    actions = c(strg$a0[ss],strg$a1[ss]),
    transitions = rbind(
                  c(strg$x00[ss],strg$x01[ss]),
                  c(strg$x10[ss],strg$x11[ss])
                  )
  )
  
}

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


apply.strg<-function(S,coop,coop_other) {
  
  state<-S$start
  
  a<-c()
  
  for (tt in 1:length(coop)) {
    
    a<-c(a,S$actions[state+1])
    
    state<-S$transitions[state+1,coop_other[tt]+1]
    
  }
  
  return(a)
  
  
}



D<- "code_SFEMQRE/00_loaddata.csv" |>
  read.csv() |>
  select(-X) |>
  select(-contains("strg")) |>
  group_by(treatmentID,id,match) 

D.counts<-D |>
  group_by(treatmentID,r,g,l,delta,id) |>
  summarize(
    count = n()
  ) |>
  ungroup() |>
  mutate(TID = treatmentID |> as.factor() |> as.numeric()) 
  

gld<-D.counts |>
  group_by(TID) |>
  summarize(
    r = mean(r),
    g = mean(g),
    l = mean(l),
    delta = mean(delta)
  )


for (ss in 1:length(strgList)) {
  
  print(paste(ss,"of",length(strgList)))
    
  d<-D |>
    group_by(treatmentID,id,match) |>
    mutate(
      x = apply.strg(strgList[[ss]],coop,coop_other)
    ) |>
    group_by(r,g,l,delta,id) |>
    summarize(
      follow = sum(coop==x)
    )
  
  name<-paste0("strg",ss)
  
  D.counts[,name]<-d$follow
  
}



#-------------------------------------------------------------------------------
# Estimation starts here
#-------------------------------------------------------------------------------

SFEMQRE<-"code_SFEMQRE/SFEMQRE_pooled.stan" |>
  stan_model()

SFEM<-"code_SFEMQRE/SFEM_pooled.stan" |>
  stan_model()

#-------------------------------------------------------------------------------
# Fully pooled QRE model
#-------------------------------------------------------------------------------

file<-"code_SFEMQRE/05_manystrategies_pooledQRE.csv"

if (!file.exists(file)| RERUN) {
  
  dStan<-list(
    N = dim(D.counts)[1],
    nstrg = length(strgList),
    follow = D.counts |> select(contains("strg")),
    treatmentID = D.counts$TID,
    nTreatment = D.counts$TID |> unique() |> length(),
    
    g = gld$g,
    l = gld$l,
    delta = gld$delta,
    
    count = D.counts$count,
    
    numFollow = numFollow,
    nextState = nextState,
    firstState = firstState,
    
    prior_lambda = c(1.15,1.76),
    prior_MIX = rep(1,length(strgList)),
    
    lGridSize = 50,
    nc = 100,
    
    ftol = 1e-6,
    UseData = rep(1,dim(D.counts)[1])
    
    
  )
  
  Fit.SFEMQRE<-SFEMQRE |>
    sampling(
      data=dStan,seed=42
    )
  
  bs<-bridge_sampler(Fit.SFEMQRE)
  
  summary(Fit.SFEMQRE)$summary |>
    data.frame() |>
    rownames_to_column(var = "par") |>
    mutate(
      logml = bs$logml
    ) |>
    write.csv(file)
  
}


