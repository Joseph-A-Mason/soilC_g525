#This script was written to model OC in two horizons, each with a single 
#pool and input from roots or surface litter only. It also models 14C in
#soil organic matter and d13C change after a vegetation change from
#prairie to forest (could also be prairie to a C3 crop like wheat)

#load SoilR and ggplot2 packagesggplot()
library("SoilR")
library("ggplot2")
        
#This part controls input to the model of atmospheric 14C content
#Can use the reconstructed content as it varied over time from
#3000 years ago up to just before nuclear bomb testing, or 
#the average over that time.
#1.Use IntCal20. 2.Use Delta-14C = Average for time period
#Change IntCal20_3000 to other time series as need.
Intcal.switch <- 2
IntCal20_3000<-IntCal20$Delta.14C[6502:9501]
MeanD14C <- mean(IntCal20_3000)

#Here you enter the name of the profile and horizons.
#Substitute profile name for P98 if you use another one
ProfName<-"P98 Single-Pool"
H1<-"A1"
H2<-"Bw"

#Decomposition constant (k) for each of the two horizons. When you open this project, the
#values will be set for the A and Bw horizons of P-98. You will change them as a step in
#carrying out the project. If you want to start with values for a different profile, some
#are provided in the instructions.

Ok<-vector(length=2)
#First (= upper, fast-cycling) horizon
Ok[1]<--0.1
#Second (= lower, slow-cycling) horizon
Ok[2]<--0.00078

Nk<-vector(length=2)
Nk[1]<-Ok[1]
Nk[2]<-Ok[2]

#Input factor (default value = 1.0). In this model, the input of SOC from roots and leaf 
#litter is calculated to produce a steady state of SOC, based on the
#initial SOC inventory. However, if you want you can reduce or
#increase the input factor to represent change in vegetation or
#environmental conditions. For example, you could make the factor 0.8
#to represent the assumption that input from forest vegetation is 80%
#of input from prairie vegetation
In_factor<-1.0


#obsd13C is a list of the observed (modern) values of d13C (normalized
#ratio of 13C to 12C) in H1 and H2, in that order. You should change these
#if you use a profile other than P98. Values for other profiles are in the
#project instructions.

obsd13C<-(c(-25.67, -21.28))

#Now we need to enter the initial OC inventory of the two parts of 
#the single pool of OC in each horizon (Mg/ha, mass per unit area of
#the land surface).This is where understanding the model gets a little
#complicated. We are assuming there is one pool of OC in each horizon,
#to make modeling simple. But to model 14C content and d13C after vegetation
#change, we need to think of that single pool as having two parts, one
#that is OC from the older prairie vegetation, and the other that is new
#OC from forest. The old prairie-derived OC is initially all of the OC in
#the pool, but over time it's replaced by new forest-derived OC. 

#The initial OC inventories you see when you first open this script, 
#combined with the initial k values (decomposition constants) will 
#result in a steady state model in which input of OC = output of OC 
#and total OC in each pool does not change over time. You will then
#need to try out changes in the initial values of OC inventory and
#k (decomposition constants) and see how the results change.
#observed bulk d13C values, Horizon 1, Horizon 2

#OCInv is a vector, a list of values, for initial OC in each part of
#the single pool of OC in Horizons 1 and 2. The values are inside the
#parentheses ("c(...)"). In order, the values there are Horizon 1
#(prairie-derived) C, Horizon 1 New (forest-derived) C, 
#Horizon 2 Old C, and Horizon 2 New C. You will leave the two zeros,
#initial values of new forest-derived C, and change the first and third
#values representing initial prairie-derived C.

OCinv<-c(4920, 0.0, 834, 0.0)

#Initial Delta-14C (essebtially the 14C content relative to other 
#isotopes of C) of each horizon. To use horizons from another profile
#replace the four values listed for InitC14t with values for that
#profile given in the project instructions. Just as with OC inventory,
#the values here are list in the order, Horizon 1 Old C, 
#Horizon 1 New C, Horizon 2 Old C, Horizon 2 New C

InitC14t<-c(-8.87, -1.24, -140.9,-1.24)

Old<-"-Old"
New<-"-New"
Total<-"Total"
Bulk<-"Bulk"

#The function "d13C_w_decomp" models d13C change for a single pool of SOC from one of two
#C sources: Old, from prairie, now usually with 0 input, and new, from forest, with ongoing input
#Includes effect of decomposition, etc., based on observed trends of d13C with increasing Delta-14C
#May also include effect of transfer between pools or horizons, indirectly, if that transfer is
#included in basic models of total SOC and Delta-14C, e.g. C transfer from overlying horizons
#or faster-cycling pools will lead to less negative Delta-14C in pool being modeled, shifting
#d13C of the pool toward more negative values reflecting less decomposition
d13C_w_decomp <- function(nyr, Ind13C, Intercd13C, C14t, SlopeC14t){
  d13C1<- matrix(0.0, nyr, 1)
  d13C1[1]<- Ind13C
  for(m in 2:nyr) {
    d13C1[m]<-(SlopeC14t*C14t[m])+Intercd13C
  }
  return(d13C1)
}

#This section calculates the d13C value for old (prairie-derived) and new
#(forest-derived) OC, based on its 14C content, which reflects its residence
#time in the soil. You don't need to change anything here.

#d13C with Delta-14C intercept for C from old and new vegetation
#Horizon 1 Old C, Horizon 1 New C, Horizon 2 Old C, Horizon 2 New C
Intercd13C<-c(-18.0,-26.7, -18.0, -26.7) 
#slope of d13C with Delta-14C for C from old and new vegetation
#Horizon 1 Old C, Horizon 1 New C, Horizon 2 Old C, Horizon 2 New C
SlopeC14t<-c(-0.003, -0.008, -0.003, -0.008)

#From here on, leave the code unchanged for the project (you can try making changes
#if you are interested, of course)

#Root input from new and old vegetation (old is typically zero but doesn't have to be)
NRootIn<-c((In_factor*-1.0*Ok[1]*OCinv[1]),(In_factor*-1.0*Ok[2]*OCinv[3]))
ORootIn<-c(0,0)

#timescale for simulation
numyr<-3000
t_start<-1
t_end<-3000
timestep<-1
t=seq(t_start, t_end, timestep)

#Time to reach observed d13C, for the two horizons (will be recalculated)
H1_Time<-0.0
H2_Time<-0.0

#Create output matrix and add year numbers
d13Clist<- matrix(0.0,nrow=numyr, ncol=23)
d13Clist[1, 1] <- 1
for (i in 2:numyr) {
  d13Clist[i, 1] <- d13Clist[i-1, 1]+1
}

#Create A matrix for both horizons 
A=new(Class="BoundLinDecompOp",
      t_start,
      t_end,
      function(t0){
        matrix(nrow=4, ncol=4, byrow=TRUE,
          c(Ok[1], 0, 0, 0,
            0, Nk[1], 0, 0,
            0, 0, Ok[2], 0,
            0, 0, 0, Nk[2])
          )
      }
)

#Initial DC14 values, structured for use in SoilR General Model
F0<-ConstFc(c(InitC14t[1], InitC14t[2],InitC14t[3], InitC14t[4]), "Delta14C")

#Time, atmospheric 14C, and C input series
if(Intcal.switch==1){
  prebatmosdf<-data.frame(t[1:numyr], IntCal20_3000)
  inputFc=BoundFc(prebatmosdf, format="Delta14C")
  inputFluxes=new(
    "TimeMap",
    t_start,t_end,
    function(t0){matrix(nrow=4, ncol=1, c(ORootIn[1], NRootIn[1], ORootIn[2], NRootIn[2]))}
  )
} else if(Intcal.switch==2){
  prebatmosdf<-data.frame(t[1:numyr], MeanD14C)
  inputFc=BoundFc(prebatmosdf, format="Delta14C")
  inputFluxes=new(
    "TimeMap",
    t_start,t_end,
    function(t0){matrix(nrow=4, ncol=1, c(ORootIn[1], NRootIn[1], ORootIn[2], NRootIn[2]))}
  )
}

#Run one pool (two horizon) 13C models for C from old vegetation, filling results into output matrix
HorizModel<-GeneralModel_14(t=t,A=A, ivList=OCinv, initialValF = F0,inputFluxes = inputFluxes, inputFc = inputFc, di=-0.0001209681, solverfunc = deSolve.lsoda.wrapper,
                              pass = TRUE)
#these two lines replace very small Ct values with 0.0001--only needed for old input to allow calculation of d13C
d13Clist[,2:5]<-getC(HorizModel)
d13Clist[,2:5]<-replace(d13Clist[,2:5], d13Clist[,2:5]<0.0001, 0.0001)

#D14C time series
d13Clist[,6:9]<-getF14(HorizModel)

#Replaces NaN for initial D14C of new C
d13Clist[1,6:9]<- InitC14t[1:4]

#these loops replace erratic D14C values generated after pool becomes very small, with stable values--only needed for old input
for(e in 6:9){
  for(f in 2:numyr){
    if(d13Clist[f,e-4]<=0.0001){d13Clist[f,e]<-d13Clist[f-1,e]}
  }
}

#Fills in d13C values for new and old C in each horizon, adjusting for D14C
Ind13C<-c(0.0, 0.0, 0.0, 0.0)
for (j in 1:4) {
    #for (k in 1:4){
    Ind13C[j]<-(SlopeC14t[j]*d13Clist[1, j+5])+Intercd13C[j]
    #}
    d13Clist[,j+9] <- d13C_w_decomp(nyr=numyr, Ind13C=Ind13C[j], Intercd13C=Intercd13C[j], C14t=d13Clist[,j+5], SlopeC14t=SlopeC14t[j])
}

#Calculate bulk values, based on simple mixing model

#Total C for the two horizons
d13Clist[,14]<-d13Clist[,2]+d13Clist[,3]
d13Clist[,15]<-d13Clist[,4]+d13Clist[,5]
#Delta-14C for the two horizons
d13Clist[,16]<-(d13Clist[,6]*(d13Clist[,2]/d13Clist[,14]))+(d13Clist[,7]*(d13Clist[,3]/d13Clist[,14]))
d13Clist[,17]<-(d13Clist[,8]*(d13Clist[,4]/d13Clist[,15]))+(d13Clist[,9]*(d13Clist[,5]/d13Clist[,15]))
#d13C for the two horizons
d13Clist[,18]<-(d13Clist[,10]*(d13Clist[,2]/d13Clist[,14]))+(d13Clist[,11]*(d13Clist[,3]/d13Clist[,14]))
d13Clist[,19]<-(d13Clist[,12]*(d13Clist[,4]/d13Clist[,15]))+(d13Clist[,13]*(d13Clist[,5]/d13Clist[,15]))

#Difference from observed d13C for the two horizons
d13Clist[,20]<- abs(d13Clist[,18]-obsd13C[1])
d13Clist[,21]<- abs(d13Clist[,19]-obsd13C[2])

#Conventional radiocarbon age for the two horizons
d13Clist[,22]<--8033*log(AbsoluteFractionModern_from_Delta14C(d13Clist[,16]))
d13Clist[,23]<--8033*log(AbsoluteFractionModern_from_Delta14C(d13Clist[,17]))

#Write results and parameters to dataframe
d13Cdf <- as.data.frame(d13Clist)
names(d13Cdf)<- c("Year", "H1_OCt","H1_NCt","H2_OCt","H2_NCt","H1_OC14t","H1_NC14t",
                  "H2_OC14t","H2_NC14t", "H1_Od13C","H1_Nd13C","H2_Od13C","H2_Nd13C",
                  "H1_TCt","H2_TCt","H1_TC14t", "H2_TC14t","H1_Td13C","H2_Td13C",
                  "H1_Diffd13C", "H2_Diffd13C", "H1_14CAge", "H2_14CAge")

#Output of files listing parameters and results would go here

#create a set of plots of the results

ggplot(d13Cdf, aes(Year)) + 
  geom_line(aes(y=H1_TCt,linetype="H1"), colour="black",size=0.75)+ 
  geom_line(aes(y=H2_TCt, linetype="H2"), colour="black", size=0.75)+ 
  scale_linetype_manual("Horizon", values=c("H1"=2, "H2"=1), 
                        labels=c(H1, H2))+
  labs(x="Years", y="Total C")+
  theme(aspect.ratio = 0.3,
        panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major=element_line(colour="gray"),
        panel.grid.minor=element_line(colour="gray"),
        axis.title.y = element_text(size = rel(1), angle = 90),
        axis.title.x = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        legend.position="right",
        legend.key.width=unit(1.5,"cm"),
        legend.key =element_rect(fill="transparent"))

ggplot(d13Cdf, aes(Year)) + 
  geom_line(aes(y=H1_TC14t,linetype="H1"), colour="black",size=0.75)+ 
  geom_line(aes(y=H2_TC14t, linetype="H2"), colour="black", size=0.75)+ 
  scale_linetype_manual("Horizon", values=c("H1"=2, "H2"=1), 
                        labels=c(H1, H2))+
  labs(x="Years", y=expression(paste(Delta^14,"C")))+
  theme(aspect.ratio = 0.3,
        panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major=element_line(colour="gray"),
        panel.grid.minor=element_line(colour="gray"),
        axis.title.y = element_text(size = rel(1), angle = 90),
        axis.title.x = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        legend.position="right",
        legend.key.width=unit(1.5,"cm"),
        legend.key =element_rect(fill="transparent"))

ggplot(d13Cdf, aes(Year)) + 
  geom_line(aes(y=H1_14CAge,linetype="H1"), colour="black",size=0.75)+ 
  geom_line(aes(y=H2_14CAge, linetype="H2"), colour="black", size=0.75)+ 
  scale_linetype_manual("Horizon", values=c("H1"=2, "H2"=1), 
                        labels=c(H1, H2))+
  labs(x="Years", y=expression(paste(""^14*"C Age (yr)")))+
  theme(aspect.ratio = 0.3,
        panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major=element_line(colour="gray"),
        panel.grid.minor=element_line(colour="gray"),
        axis.title.y = element_text(size = rel(1), angle = 90),
        axis.title.x = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        legend.position="right",
        legend.key.width=unit(1.5,"cm"),
        legend.key =element_rect(fill="transparent"))
  
ggplot(d13Cdf, aes(Year)) + 
  geom_segment(aes(x=0, y=obsd13C[1], xend=numyr, yend=obsd13C[1]), colour="gray", linetype=4, size=1.2)+
  geom_segment(aes(x=0, y=obsd13C[2], xend=numyr, yend=obsd13C[2]), colour="gray", linetype=3, size=1.2)+
  geom_line(aes(y=H1_Td13C,linetype="H1"), colour="black",size=0.75)+ 
  geom_line(aes(y=H2_Td13C, linetype="H2"), colour="black", size=0.75)+ 
  scale_linetype_manual("Horizon", values=c("H1"=4, "H2"=3), 
                        labels=c(H1, H2))+
  labs(x="Years", y=expression(paste(delta^13,"C")))+
  theme(aspect.ratio = 0.3,
        panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major=element_line(colour="gray"),
        panel.grid.minor=element_line(colour="gray"),
        axis.title.y = element_text(size = rel(1), angle = 90),
        axis.title.x = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        legend.position="right",
        legend.key.width=unit(1.5,"cm"),
        legend.key =element_rect(fill="transparent"))

#Find and print time for d13C to reach observed value for each horizon
H1_Time <- d13Cdf$Year[which.min(d13Cdf$H1_Diffd13C)]
H2_Time <- d13Cdf$Year[which.min(d13Cdf$H2_Diffd13C)]
print(paste("Time for Horizon 1 to reach observed d13C:", H1_Time, "years"))
print(paste("Time for Horizon 2 to reach observed d13C:", H2_Time, "years"))

#Check on steady state. You can check this first by looking at the plot
#of total C to make sure it remains constant (horizontal lines) throughout
#the simulation time. As an added check for radiocarbon steady state, 
#run these lines. If there is minimal change in start and end values 
#for D14C in each horizon, then steady state for radiocarbon as well.

d13Cdf$H1_TC14t[1]
d13Cdf$H1_TC14t[3000]
d13Cdf$H2_TC14t[1]
d13Cdf$H2_TC14t[3000]

