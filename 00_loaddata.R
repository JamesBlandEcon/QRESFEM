library(tidyverse)



D<-"code_SFEMQRE/DBF2011.txt" |>
  read.table(skip=11, header=TRUE,sep="\t") |>
  group_by(group,round) |>
  mutate(
    coop_other = sum(coop)-coop
  ) |>
  ungroup() |>
  # code up strategies
  arrange(r,delta,id,match,round)  |>
  group_by(
    r,delta,id
  ) |>
  mutate(
    interaction = 1:n()
  ) |>
  filter(interaction>110) |>
  group_by(r,delta,id,match) |>
  filter(
    min(round)==1
  ) |>
  ungroup() |>
  mutate(
    S = 12,D=50,P=25,
    g = (D-P)/(r-P)-1, l=-(S-P)/(r-P),
    treatmentID = paste(round(g,2),round(l,2),delta)
  )

Dstrg<-D |>
  group_by(treatmentID,r,g,l,delta,id,match) |>
  mutate(
    strgAD = 0,
    strgAC = 1,
    strgGrim = ifelse(round==1,1,ifelse(lag(cumsum(coop_other))==round,1,0)),
    strgTFT  = ifelse(round==1,1,lag(coop_other)),
    strgSTFT = ifelse(round==1,0,lag(coop_other))
  )

Dstrg |>
  write.csv("code_SFEMQRE/00_loaddata.csv")

