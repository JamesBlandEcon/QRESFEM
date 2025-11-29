library(tidyverse)
library(kableExtra)
library(latex2exp)
library(Cairo)
library(viridis)
library(bridgesampling)

library(rstan)

treatments<-"code_SFEMQRE/01_estimate_treatments.csv" |>
  read.csv()

strgList<-c("AD","AC","Grim","TFT","STFT")

estimates<-rbind(
  summary("code_SFEMQRE/01_estimate_SFEMQRE.rds" |> readRDS())$summary |>
    data.frame()|>
    rownames_to_column(var="par") |>
    mutate(
      model = "QRE-SFEM"
    ),
  summary("code_SFEMQRE/01_estimate_SFEM.rds" |> readRDS())$summary |>
    data.frame()|>
    rownames_to_column(var="par") |>
    mutate(
      model = "SFEM"
    ) 
  
)

#-------------------------------------------------------------------------------
# Posterior distribution of lambda and epsilon from the pooled model
#-------------------------------------------------------------------------------

Fit<-"code_SFEMQRE/01_estimate_SFEMQRE.rds" |> readRDS()

LambdaEpsilon<-tibble(
  lambda = extract(Fit)$lambda,
  epsilon = extract(Fit)$epsilon
) |>
  pivot_longer(
    cols = lambda:epsilon
  ) |>
  mutate(
    prior_density = ifelse(name=="lambda",dlnorm(value,mean=log(10),sd=1.17),2)
  ) |>
  mutate(
    param = name |> factor(levels = c("lambda","epsilon"))
  )

plt<-(
  ggplot(LambdaEpsilon,aes(x=value))
  +stat_density()
  +theme_bw()
  +geom_line(aes(y=prior_density),color="red",linewidth=1)
  +facet_wrap(~param,scales="free",labeller = label_parsed)
  +xlab("")+ylab("density")
)
cairo_pdf(filename=paste("code_SFEMQRE/pooledPosterior.pdf",sep=""),width=7,height=4,family="CMU Serif")
(plt) 
dev.off()
plt |> print()


#-------------------------------------------------------------------------------
# Pooled QRE-SFEM estimates
#-------------------------------------------------------------------------------

fmt<-"%.3f"

TAB<-estimates |>
  filter(grepl("MIX",par))|>
  filter(model=="QRE-SFEM") |>
  mutate(
    TID = par |> str_split_i(",",1) |> parse_number(),
    strgID = par |> str_split_i(",",2) |> parse_number(),
    strg = strgList[strgID]
  ) |>
  full_join(treatments,by="TID") |>
  arrange(r,delta,model) |>
  pivot_longer(
    cols = c(mean,sd),
    names_to = "name",
    values_to = "value"
  ) |>
  mutate(
    msd = paste0(ifelse(name=="sd","(",""),sprintf(fmt,value),ifelse(name=="sd",")",""))
  ) |>
  pivot_wider(
    id_cols = c(strg,name),
    names_from = c(r,delta,model),
    values_from = msd
  ) |>
  group_by(strg) |>
  mutate(
    strg = ifelse(is.na(lag(strg)),strg,ifelse(strg==lag(strg),"",strg))    
  ) |>
  select(-name)|>
  as.matrix()

colnames(TAB)<-c(c("","$R=32$","$R=40$","$R=48$","$R=32$","$R=40$","$R=48$"))

TAB |>
  kbl() |>
  kable_classic(full_width=FALSE) |>
  add_header_above(c("","$\\delta=0.50$"=3,"$\\delta=0.75$"=3))

TAB |>
  kbl(booktabs=TRUE,escape=FALSE,format = "latex", linesep = c("","\\addlinespace")) |>
  add_header_above(c("","$\\\\delta=0.50$"=3,"$\\\\delta=0.75$"=3),escape=FALSE) |>
  save_kable(file = "code_SFEMQRE/pooled_SFEMQRE.tex",escape=FALSE)
  

#-------------------------------------------------------------------------------
# QRE trace plots
#-------------------------------------------------------------------------------

trace<-"code_SFEMQRE/02_traceQRE.csv" |>
  read.csv() |>
  mutate(
    epstxt = paste0("\u03b5 = ",sprintf("%.2f",epsilon)),
    deltatxt = paste0("\u03b4 = ",sprintf("%.2f",delta)),
    Rtxt =  paste0("R = ",r)
  )

plotThis<-trace |>
  filter(epsilon == 0.01 | epsilon==0.1)

plt<-(
  ggplot(plotThis,aes(x=lambda,y=mix,color=strg,linetype=epstxt))
  +geom_path()
  +facet_grid(Rtxt~deltatxt)
  +scale_x_continuous(trans="log10")
  +theme_bw()
  +xlab(expression(lambda))
  +ylab("Strategy frequency")
  +labs(color=NULL,linetype=NULL)
)
cairo_pdf(filename=paste("code_SFEMQRE/QREtrace.pdf",sep=""),width=7,height=4,family="CMU Serif")
(plt) 
dev.off()
plt |> print()

#-------------------------------------------------------------------------------
# separate model per treatment
#-------------------------------------------------------------------------------

fmt<-"%.3f"

ESTIMATES <- "code_SFEMQRE/03_estimate_by_treatment.csv" |>
  read.csv() |>
  filter(
    !grepl("logmix",par) & !grepl("log_lik",par) & par !="lp__"
  ) |>
  mutate(
    parameter = strgList[par |> str_split_i(",",2) |> parse_number()],
    parameter = ifelse(grepl("epsilon",par),"$\\epsilon$",parameter),
    parameter = ifelse(grepl("lambda",par),"$\\lambda$",parameter)
  ) |>
  arrange(delta,r,model) |>
  pivot_longer(
    cols = c(mean,sd)
  ) |>
  mutate(
    msd = paste0(ifelse(name=="sd","(",""),sprintf(fmt,value),ifelse(name=="sd",")",""))
  ) |>
  pivot_wider(
    id_cols = c(parameter,name),
    names_from = c(model,r,delta),
    values_from = msd,
    values_fill = ""
  ) |>
  ungroup() |>
  mutate(
    parameter = ifelse(is.na(lag(parameter)),parameter,ifelse(parameter==lag(parameter),"",parameter)) 
  ) |>
  select(-name) |>
  as.matrix() 
  

colnames(ESTIMATES)<-ifelse(grepl("QRE",colnames(ESTIMATES)),"QRE","SFEM")
colnames(ESTIMATES)[1]<-""


ESTIMATES |>
  kbl()|>
  kable_classic(full_width=FALSE) |>
  add_header_above(
    c("","$R=32$"=2,"$R=40$"=2,"$R=48$"=2,"$R=32$"=2,"$R=40$"=2,"$R=48$"=2)
  ) |>
  add_header_above(
    c("","$\\delta=0.50$"=6,"$\\delta=0.75$"=6)
  )
ESTIMATES |>
  kbl(format="latex",booktabs=TRUE,escape=FALSE,linesep = c("","\\addlinespace"))|>
  add_header_above(
    c("","$R=32$"=2,"$R=40$"=2,"$R=48$"=2,"$R=32$"=2,"$R=40$"=2,"$R=48$"=2),escape=FALSE
  ) |>
  add_header_above(
    c("","$\\\\delta=0.50$"=6,"$\\\\delta=0.75$"=6),escape=FALSE
  ) |>
  save_kable(
    file = "code_SFEMQRE/separateModels.tex",escape=FALSE
  )

# posterior probabilities

POSTPROBS<-"code_SFEMQRE/03_estimate_by_treatment_postprobs.csv" |>
  read.csv() |>
  select(-X) |>
  arrange(delta,r) |>
  mutate(
    delta = paste0("$\\delta=",delta,"$")
  )|>
  rename(
    `$R$`=r,
    `posterior probability` = postprob.SFEMQRE
  )|>
  pivot_wider(
    id_cols = `$R$`,
    names_from = delta,
    values_from = `posterior probability`
  )

POSTPROBS |>
  kbl(digits = 4,booktabs=TRUE,escape=FALSE,format="latex") |>
  save_kable("code_SFEMQRE/separateModelsPostprobs.tex")

#-------------------------------------------------------------------------------
# Leave-out analysis
#-------------------------------------------------------------------------------

SFEM<- "code_SFEMQRE/03_estimate_by_treatment.csv" |>
  read.csv() |>
  filter(!grepl("epsilon",par)) |>
  filter(
    !grepl("logmix",par) & !grepl("log_lik",par) & par !="lp__" & model=="SFEM"
  )  |>
  full_join(
    gld,
    by = c("r","delta")
  ) |>
  mutate(
    strg = strgList[par |> str_split_i(",", 2) |> parse_number()]
  ) |>
  rename(
    SFEM= mean
  ) |>
  select(TID,strg,SFEM) 

d<-"code_SFEMQRE/04_leaveout.csv" |>
  read.csv() |>
  filter(
    grepl("MIX",par)
  ) |>
  mutate(
    leaveout = par |> str_split_i(",", 1) |> parse_number(),
    strg = strgList[par |> str_split_i(",", 2) |> parse_number()]
  ) |>
  filter(
    TID==leaveout
  ) |>
  full_join(
    SFEM,
    by =c("strg","TID")
  ) |>
  mutate(
    treatment = paste0("R =",r,", \u03b4 =", sprintf("%.2f",delta))
  )

plt<-(
  ggplot(d,aes(x=SFEM,y=mean,label=strg,color = treatment))
  +geom_text()
  #+geom_smooth()
  +xlab("SFEM") 
  +ylab("QRE-SFEM (out of sample)")
  +geom_abline(slope=1,intercept=0,linetype="dashed",color="black")
  +theme_bw()
)
cairo_pdf(filename=paste("code_SFEMQRE/leaveout.pdf",sep=""),width=7,height=4,family="CMU Serif")
(plt) 
dev.off()
plt |> print()


#-------------------------------------------------------------------------------
# Pull out some posterior probabilities
#-------------------------------------------------------------------------------

(postprobs<-"code_SFEMQRE/01_estimate_postprobs.rds" |> readRDS())
