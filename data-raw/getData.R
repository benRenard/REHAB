library(airGR)
library(dplyr)

D=read.table(file.path('data-raw','V5064010_SAFRAN.txt'),skip=29,header=TRUE,sep=';')
dates=as.POSIXlt(as.character(D$Date),format="%Y%m%d")
P=D$Ptot
ETP=D$ETP
Qobs=1000*((D$Q*0.001*60*60*24)/(2240*1000000))
theta=c(X1=125,X2=0.6,X3=160,X4=1.5)

inputs <- CreateInputsModel(FUN_MOD=RunModel_GR4J,
                            DatesR=dates,Precip=P,PotEvap=ETP)
options <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                            InputsModel=inputs,IndPeriod_Run=1:length(dates),
                            warning=FALSE,verbose=FALSE)
output <- RunModel(InputsModel = inputs,
                   RunOptions = options, Param = theta,
                   FUN_MOD = RunModel_GR4J)

ArdecheRiver=data.frame(date=dates,obs=Qobs,sim=output$Qsim) %>%
  filter(between(date,as.Date('1998-09-01'),as.Date('2008-08-31')))
plot(ArdecheRiver$date,ArdecheRiver$obs,log='y')
lines(ArdecheRiver$date,ArdecheRiver$sim,col='red')
plot(log(ArdecheRiver$sim),ArdecheRiver$obs-ArdecheRiver$sim)

save(ArdecheRiver,file=file.path('data','ArdecheRiver.RData'))
