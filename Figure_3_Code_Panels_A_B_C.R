######################################################
################# Fig 3, Panels D-G ##################
######################################################


############################
####### Load Packages ######
############################
library(tidyverse)
library(smplot2)
library(cowplot)
library(ggbreak) 
library(patchwork)
library(ggh4x)
library(ggplot2)

# This will reproduce Fig 3, Panels A,B,C


############################
######### Panel A #########
############################

# note that different functions will be used for each panel 

############################
##### Define functions #####
############################

# Hawk_Dove_Change: returns an single timepstep of the stochastic Hawk-Dove game with reflecting boundary
### Parameters###
#D number of Doves, H number of Hawks, b0 baseline births per timestep, d0 mortality probability per timestep
# B benefit of getting resource, C Cost of fighting, s strength of selection, a strength of density-dependence
# u mutation rate (not used in any simulations), k (a dummy parameter; never used, but was important for testing; kindly ignore it) 
Hawk_Dove_Change <- function(H,D,b0,d0,B,C,s,a,u,k){
  
  
  H2 <- rbinom(1,H,1-d0) # surviving Hawks 
  D2 <- rbinom(1,D,1-d0) # surviving Doves 
  
  if(H2==0) H2 <- 1 # reflecting boundary
  if(D2==0) D2 <- 1 # reflecting boundary
  
  H_prop <- (H2)/(H2+D2) # proportion Hawk
  D_Prop <- (D2)/(H2+D2) # proportion Dove 
  
  if(is.na(H_prop)) H_prop<-0 # make any potential divide by zero -> 0 (shouldn't ever happen with boundary)
  if(is.na(D_Prop)) D_Prop<-0 # make any potential divide by zero -> 0 (shouldn't ever happen with boundary)
  
  DD <- exp(-a*(H2+D2)) # Density dependence term
  
  Birth_Hawk   <-  rpois(1, max(0, H2*b0 *exp(s*(H_prop*(B-C)*0.5 + D_Prop*B) )*DD     ) ) # Hawk births 
  Birth_Dove   <-  rpois(1, max(0, D2*b0 *exp(s*(D_Prop*B*0.5))*DD )) # Dove Births
  
  H3 <- H2 + Birth_Hawk # New Hawk abundance
  D3 <- D2 + Birth_Dove # New Dove abundance
  
  return(c(H3,D3))
}

# Simulation: stochastic simulations of Hawk-Dove game for a set Timesteps 
Simulation <- function(Timesteps,H0,D0,b0,d0,B,C,s,a,u,k){
  
  Hawks    <- rep(NA,Timesteps) # make empty vector
  Hawks[1] <- H0 # Hawk initial
  
  Doves    <- rep(NA,Timesteps) # make empty vector
  Doves[1] <- D0 # Dove initial 
  
  for(Time in 2:Timesteps){ # for all timesteps 
    
    New_Abundances <- Hawk_Dove_Change(Hawks[Time-1],Doves[Time-1],b0,d0,B,C,s,a,u,k) # run stochastic model 
    Hawks[Time]    <- New_Abundances[1] # Save in vector 
    Doves[Time]    <- New_Abundances[2] # Save in vector 
  }
  
  return(list(Hawks,Doves)) # return as list 
}

# Hawk_Dove_Change_det: returns an single timepstep of the deterministic Hawk-Dove game 
Hawk_Dove_Change_det <- function(H,D,b0,d0,B,C,s,a,u,k){
  
  # Basically, the same as "Hawk_Dove_Change", but fully deterministic 
  
  H2 <- (1-d0)*H
  D2 <- (1-d0)*D
  
  H_prop <- (H2)/(H2+D2)
  D_Prop <- (D2)/(H2+D2)
  
  DD <- exp(-a*(H2+D2))
  
  Birth_Hawk   <-  H2*b0 *exp(s*(H_prop*(B-C)*0.5 + D_Prop*B) )*DD     
  Birth_Dove   <-  D2*b0 *exp(s*(D_Prop*B*0.5))*DD 
  
  Hawk_To_Dove  <- Birth_Hawk*u
  Dove_To_Hawk  <- Birth_Dove*u
  
  H3 <- H2 + Birth_Hawk + Dove_To_Hawk - Hawk_To_Dove 
  D3 <- D2 + Birth_Dove + Hawk_To_Dove - Dove_To_Hawk 
  
  
  return(c(H3,D3))
}


# Simulation: deterministic simulations of Hawk-Dove game for a set Timesteps 
Simulation_det <- function(Timesteps,H0,D0,b0,d0,B,C,s,a,u,k){
  # The same as "Simulation", but fully deterministic 
  
  Hawks    <- rep(NA,Timesteps)
  Hawks[1] <- H0
  
  Doves    <- rep(NA,Timesteps)
  Doves[1] <- D0
  
  for(Time in 2:Timesteps){
    
    New_Abundances <- Hawk_Dove_Change_det(Hawks[Time-1],Doves[Time-1],b0,d0,B,C,s,a,u,k)
    Hawks[Time]    <- New_Abundances[1]
    Doves[Time]    <- New_Abundances[2]
  }
  
  return(list(Hawks,Doves))
}


############################
##### Set parameters #####
############################

Timesteps <- 75000/2.4 # timesteps 
H0 <- 10 #initial conditions Hawk 
D0 <- 10 # initial conditions Dove
b0 <- .05 # baseline birth 
d0 <- .01 # death probability 
B <- 2 # benefit of gaining resource 
C <- 4.4 # cost of fighting for Hawks 
s <- .85 # strength of selection 
a <- .1 # density-dependence parameter 
u <- 0 # mutation rate (set to zero; ignore)
k=1 # dummy parameter (ignore)


set.seed(524) # set seed for simulation
Out <- Simulation(Timesteps,H0,D0,b0,d0,B,C,s,a,u,k) # stochastic simulation 
T_Sample <- seq(1,Timesteps,by=100) # vector; will sample Time series of outputs every 100 steps 
MX <- max(c(Out[[1]],Out[[2]])) # store max of each
Out_DET <- Simulation_det(Timesteps,H0,D0,b0,d0,B,C,s,a,u,k) # deterministic simulation 
HAWKEQ <- Out_DET[[1]][Timesteps] # Equilibrium abundance of deterministic sim for Hawks
DOVEEQ <- Out_DET[[2]][Timesteps] # Equilibrium abundance of deterministic sim for Doves


# some dataframe maniulation follows below 
TS_Data <- data.frame(cbind(c(1:Timesteps),Out[[1]],Out[[2]]),Out_DET[[1]],Out_DET[[2]])
colnames(TS_Data) <- c("Time","Hawk","Dove","Hawk_Det","Dove_Det")
Sampler <- c(seq(1,nrow(TS_Data), by=50))
TS_Data_1 <- data.frame(TS_Data$Time[Sampler],TS_Data$Hawk[Sampler],TS_Data$Dove[Sampler],TS_Data$Hawk_Det[Sampler],TS_Data$Dove_Det[Sampler])
colnames(TS_Data_1) <- c("Time","Hawk","Dove","Hawk_Det","Dove_Det")
TS_Data_2  <- TS_Data_1

# make Time series plot for Panel A
Hawk_Dove_TS_Plot <- TS_Data_2 %>% ggplot(aes(x = Time, y = Hawk)) +
  geom_line(size = 1,  color = sm_color('skyblue')) +
  geom_line(aes(y = Dove), size = 1, color = "darkolivegreen3") +  # Add the Hawk line
  sm_hgrid(legends = FALSE) +
  theme(  panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
             axis.title.x=element_blank(),
             axis.title.y=element_blank(), 
          plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=16,face="bold",color="black"))  +
  geom_line(aes(y = Hawk_Det), size = 1, color = "skyblue4",alpha=.7) +  # Add the Hawk line
  geom_line(aes(y = Dove_Det), size = 1, color = "darkolivegreen4",alpha=.7) +  # Add the Hawk line
  scale_y_continuous( expand = c(0.00, .0),limits=c(0,33))+ 
  scale_x_continuous( expand = c(0.01, .0))

# Save file 
ggsave(filename = "TS_Hawk_Dove_Plot_With_Det.svg", plot = Hawk_Dove_TS_Plot, device = "svg", width = 18/1.2, height = 3.3/1.2, units = "in")




############################
######### Panel B #########
############################

####################
#### Funcitons #####
###################

# Hawk_Dove_Change: returns an single timepstep of the stochastic Hawk-Dove game 
Hawk_Dove_Change <- function(H,D,b0,d0,B,C,s,a,u,k){
  
  H2 <- rbinom(1,H,1-d0) # surviving Hawks 
  D2 <- rbinom(1,D,1-d0) # surviving Doves 
  
  H_prop <- (H2)/(H2+D2) # proportion Hawk
  D_Prop <- (D2)/(H2+D2) # proportion Dove 
  
  if(is.na(H_prop)) H_prop<-0 # make any potential divide by zero -> 0
  if(is.na(D_Prop)) D_Prop<-0 # make any potential divide by zero -> 0
  
  DD <- exp(-a*(H2+D2)) # Density dependence term
  
  Birth_Hawk   <-  rpois(1, max(0, H2*b0 *exp(s*(H_prop*(B-C)*0.5 + D_Prop*B) )*DD     ) ) # Hawk births 
  Birth_Dove   <-  rpois(1, max(0, D2*b0 *exp(s*(D_Prop*B*0.5))*DD )) # Dove Births
  
  H3 <- H2 + Birth_Hawk # New Hawk abundance
  D3 <- D2 + Birth_Dove # New Dove abundance
  
  
  return(c(H3,D3))
}

# Simulation: stochastic simulations of Hawk-Dove game for a set Timesteps 
Simulation <- function(Timesteps,H0,D0,b0,d0,B,C,s,a,u,k){
  
  Hawks    <- rep(NA,Timesteps) # make empty vector
  Hawks[1] <- H0 # Hawk initial
  
  Doves    <- rep(NA,Timesteps) # make empty vector
  Doves[1] <- D0 # Dove initial 
  
  for(Time in 2:Timesteps){ # for all timesteps 
    
    New_Abundances <- Hawk_Dove_Change(Hawks[Time-1],Doves[Time-1],b0,d0,B,C,s,a,u,k) # run stochastic model 
    Hawks[Time]    <- New_Abundances[1] # Save in vector 
    Doves[Time]    <- New_Abundances[2] # Save in vector 
  }
  
  return(list(Hawks,Doves)) # return as list 
}

# Hawk_Dove_Change_det: returns an single timepstep of the deterministic Hawk-Dove game 
Hawk_Dove_Change_det <- function(H,D,b0,d0,B,C,s,a,u,k){
  
  # Basically, the same as "Hawk_Dove_Change", but fully deterministic 
  
  H2 <- (1-d0)*H
  D2 <- (1-d0)*D
  
  H_prop <- (H2)/(H2+D2)
  D_Prop <- (D2)/(H2+D2)
  
  DD <- exp(-a*(H2+D2))
  
  Birth_Hawk   <-  H2*b0 *exp(s*(H_prop*(B-C)*0.5 + D_Prop*B) )*DD     
  Birth_Dove   <-  D2*b0 *exp(s*(D_Prop*B*0.5))*DD 
  
  Hawk_To_Dove  <- Birth_Hawk*u
  Dove_To_Hawk  <- Birth_Dove*u
  
  H3 <- H2 + Birth_Hawk + Dove_To_Hawk - Hawk_To_Dove 
  D3 <- D2 + Birth_Dove + Hawk_To_Dove - Dove_To_Hawk 
  
  
  return(c(H3,D3))
}

# Simulation_det: deterministic simulations of Hawk-Dove game for a set Timesteps 
Simulation_det <- function(Timesteps,H0,D0,b0,d0,B,C,s,a,u,k){
  # The same as "Simulation", but fully deterministic 
  
  Hawks    <- rep(NA,Timesteps)
  Hawks[1] <- H0
  
  Doves    <- rep(NA,Timesteps)
  Doves[1] <- D0
  
  for(Time in 2:Timesteps){
    
    New_Abundances <- Hawk_Dove_Change_det(Hawks[Time-1],Doves[Time-1],b0,d0,B,C,s,a,u,k)
    Hawks[Time]    <- New_Abundances[1]
    Doves[Time]    <- New_Abundances[2]
  }
  
  return(list(Hawks,Doves))
}

# Simulation_Inv_Estab_PLOT: will simulate the Hawk-Dove game several times and return a time series of a specified number of failures (3) and one establishment  
Simulation_Inv_Estab_PLOT <- function(Timesteps,Trials,H0,D0,b0,d0,B,C,s,a,u,k,Invader){
  
  
  if(Invader=="Dove")   { # set initial conditions if Dove is invading
    INVADER  <- 2
    RESIDENT <- 1
    Out_D      <-  Simulation_det(10000,10,0,b0,d0,B,C,s,a,u,k) # run deterministic model with only resident (get equilibrium)
    Out_D2     <-  Simulation_det(10000,10,10,b0,d0,B,C,s,a,u,k) # run deterministic model with both (get equilibriria)
    H_initial <-   Out_D[[RESIDENT]][10000] # initial abundance of resident set at deterministic equilibrium
    D_initial <-  1 # initial abundance of invader set at 1 
    ESTABLISHMENT <- Out_D2[[INVADER]][10000] # equilibirum of invader in deterministic model in 2 type equilibrium 
    RES_EQUIL     <- Out_D2[[RESIDENT]][10000] # save equilibrium of resident in 2 type equilibrium 
    
  }
  
  if(Invader=="Hawk")   { # same as above code, but invader is  Hawk 
    INVADER  <- 1
    RESIDENT <- 2
    Out_D     <-  Simulation_det(10000,0,10,b0,d0,B,C,s,a,u,k)
    Out_D2     <-  Simulation_det(10000,10,10,b0,d0,B,C,s,a,u,k)
    D_initial <-  Out_D[[RESIDENT]][10000]
    H_initial <-  1
    ESTABLISHMENT <- Out_D2[[INVADER]][10000]
    RES_EQUIL     <- Out_D2[[RESIDENT]][10000]
    
  }
  
  Estab_Probs <- rep(NA,Trials) # vector; keep track of establishment / failure
  Estab_Abuns <- rep(NA,Trials) # keep track of abundances 
  
  fails <- 0 # initial number of failures to reach establishment 
  SUCCESS = "NO"
  
  for(TRIAL in 1:Trials){ # for loops goes for number of "trials"; trials is a dummy variable that should be set to a high value
    
    Hawks    <- rep(NA,Timesteps) # Hawk abundance vector
    Hawks[1] <- round(H_initial) # initial condition
    
    Doves    <- rep(NA,Timesteps) # Dove abundance vector 
    Doves[1] <- round(D_initial) # initial conditions 
    
    Abundances <- cbind(Hawks,Doves) # combine Abundances into 2 column matrix
    if(SUCCESS=="YES") break # if a sucessful establishment has been recorded, break
    
    for(Time in 2:Timesteps){ # for timesteps (Timesteps should be arbitrarily large; dummy variable)
      New_Abundances <- Hawk_Dove_Change(Abundances[Time-1,1],Abundances[Time-1,2],b0,d0,B,C,s,a,u,k) # new abundances 
      Abundances[Time,]    <- c(New_Abundances[1],New_Abundances[2]) # store in matrix 
      
      
      if(New_Abundances[RESIDENT]==0) { # if resident goes extinct 
        Abundances[Time,RESIDENT] <- round(RES_EQUIL/2) # add back into community 
      }
      
      if(New_Abundances[INVADER]==0 &&  max(Abundances[,INVADER],na.rm=T) <  ESTABLISHMENT ){ # establishment fails 
        
        Estab_Probs[TRIAL] <- 0 # note failure
        break
      }
      
      if(New_Abundances[INVADER]==0 &&  max(Abundances[,INVADER],na.rm=T) >=  ESTABLISHMENT ){  # if establishment succeeds
        
        Estab_Probs[TRIAL] <- 1 # note success 
        #Estab_Abuns[TRIAL] <- sm_auc(1:Time, Abundances[,INVADER]  ) # time integral (will calculate abundance)
        break
      }
      
      
    }

    
    if(Estab_Probs[TRIAL]==1 && fails==3){ # 1st success after 3rd failure 
      ABUNS0 <- cbind(rep("SUCC",Time+10),Abundances[1:(Time+10),INVADER])   # store abundances with added dummy timesteps 
      ABUNS <- rbind(ABUNS,ABUNS0 ) # turn to vector 
      SUCCESS <- "YES" # note success 
    }
    
    
    if(Estab_Probs[TRIAL]==0 && fails==2 && max(Abundances[,INVADER],na.rm=T)>2){ # record 3rd failure 
      fails <- fails + 1 
      ABUNS0 <- cbind(rep("fail_3",Time+150),Abundances[1:(Time+150),INVADER])  
      ABUNS <- rbind(ABUNS,ABUNS0 )
    }
    
    
    if(Estab_Probs[TRIAL]==0 && fails==1 && max(Abundances[,INVADER],na.rm=T)>2){# record 2nd failure 
      fails <- fails + 1 
      ABUNS0 <- cbind(rep("fail_2",Time+150),Abundances[1:(Time+150),INVADER])  
      
      ABUNS <- rbind(ABUNS,ABUNS0 )
    }
    
    if(Estab_Probs[TRIAL]==0 && fails==0 &&max(Abundances[,INVADER],na.rm=T)>2  ){ # record 1st failure 
      fails <- fails + 1 
      ABUNS <- cbind(rep("fail_1",Time+150),Abundances[1:(Time+150),INVADER])  
    }
    
    
    
  }    
  

  
  return(ABUNS) # return abundance time series data
  
}

##########################
###### Set parameters ####
##########################

Timesteps <- 100000 # timesteps (should be large enough to ensure extinction or establishment occurs)
Trials <- 1250 # dummy variable; given value is easily large enough
H0 <- 10 # initial Hawk abundance 
D0 <- 10 # initial Dove abundance 
b0 <- .05 # baseline births
d0 <- .01 # death probability 
B <- 2 # benefit of resource 
C <- 3.2 # cost of fighting 
s <- .5 # strength of selectino 
a <- .1 # density-dependence parameter
u <- 0 # mutation rate 
k=1 # dummy parameer (ignore )

#######################
#### Run simulation ###
#######################

set.seed(4) # set seed
TS_INV_Dove<- Simulation_Inv_Estab_PLOT(Timesteps,Trials,H0,D0,b0,d0,B,C,s,a,u,k,"Dove") # run sim 
TS_INV_Dove <- data.frame(TS_INV_Dove) # save output as data frame 
TS_INV_Dove$TIME <- seq(1:nrow(TS_INV_Dove)) # add timesteps 
colnames(TS_INV_Dove) <- c("Color","y","x") # column names for plotting 
TS_INV_Dove$y <- as.numeric(TS_INV_Dove$y) # make sure values are numeric 

# some data wrangling that puts the time series into the correct order 
TS_INV_Dove_EST <- subset(TS_INV_Dove,Color=="SUCC")
ZEROS <- 1:as.numeric(rownames(TS_INV_Dove_EST)[1])
LENZ  <- as.numeric(rownames(TS_INV_Dove_EST)[1])
DFBIND <- data.frame(cbind(rep("SUCC",LENZ), rep(-1,LENZ),ZEROS))
colnames(DFBIND) <- c("Color","y","x")
TS_INV_Dove_EST <- rbind(DFBIND,TS_INV_Dove_EST)
TS_INV_Dove_EST$y <- as.numeric(TS_INV_Dove_EST$y)
TS_INV_Dove_EST$x <- as.integer(TS_INV_Dove_EST$x)
TS_INV_Dove <- head(TS_INV_Dove,-11)
TS_INV_Dove <- subset(TS_INV_Dove,Color!="SUCC")
MEAN_ABUN     <- mean(TS_INV_Dove_EST$y[which(TS_INV_Dove_EST$y>-1)])
MIN_SOJ_time  <- nrow(TS_INV_Dove)+1
MAX_SOJ_time  <-nrow(TS_INV_Dove_EST)
MIDDLE <- MAX_SOJ_time - (MAX_SOJ_time  -MIN_SOJ_time)/2

# make plot
Dove <- TS_INV_Dove_EST %>% ggplot(aes(x = x, y = y)) +
  geom_area(fill =  "darkolivegreen3", alpha = 0.4) +
  
  geom_line(data = TS_INV_Dove, aes(x, y), color = "grey65",size=1) +# Second plot overlaid
  geom_line(size = 1.1,  color ="darkolivegreen3") +
  geom_segment(aes(x = MIN_SOJ_time+20, xend = MAX_SOJ_time, y = MEAN_ABUN, yend = MEAN_ABUN),
               linetype = "dashed", color = "black", size = 1.5)+

  geom_segment(aes(x = MIN_SOJ_time-10, xend = MIN_SOJ_time-10, y = 0, yend = MEAN_ABUN),
               linetype = "solid", color = "black", size = 1.5)+ 
  geom_segment(aes(x = MAX_SOJ_time-20, xend = MAX_SOJ_time-20, y = 0, yend = MEAN_ABUN),
               linetype = "solid", color = "black", size = 1.5)+  

  sm_hgrid(legends = FALSE) +

  theme(  panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
          plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=16,face="bold",color="black"))  +
  scale_y_continuous( expand = c(0.00, .0),limits=c(0,23))+ #limits = c(0,3),breaks=c(0,1 ,2,3)) + 
  scale_x_continuous( expand = c(0.01, .0),limits=c(0,4150),breaks=c(1250,2500,3750))#, breaks=c(1,1000,2000),labels = paste0(c(1,1000,2000))) 



ggsave(filename = "Dove_TS.svg", plot = Dove, device = "svg", width = 9, height =4, units = "in")

############################
######### Panel C #########
############################

# almost everything is the same as the above (Panel A) but for the Hawk rather than the Dove

Timesteps <- 100000
Trials <- 1250
H0 <- 10
D0 <- 10
b0 <- .05
d0 <- .01
B <- 2
C <- 3.2
s <- .5
a <- .1
u <- 0.0000
k=1
Simulation_Inv_Estab_PLOT <- function(Timesteps,Trials,H0,D0,b0,d0,B,C,s,a,u,k,Invader){
  
  if(Invader=="Dove")   {
    INVADER  <- 2
    RESIDENT <- 1
    Out_D      <-  Simulation_det(10000,10,0,b0,d0,B,C,s,a,u,k)
    Out_D2     <-  Simulation_det(10000,10,10,b0,d0,B,C,s,a,u,k)
    H_initial <-   Out_D[[RESIDENT]][10000]
    D_initial <-  1
    ESTABLISHMENT <- Out_D2[[INVADER]][10000]
    RES_EQUIL     <- Out_D2[[RESIDENT]][10000]
    
  }
  
  if(Invader=="Hawk")   {
    INVADER  <- 1
    RESIDENT <- 2
    Out_D     <-  Simulation_det(10000,0,10,b0,d0,B,C,s,a,u,k)
    Out_D2     <-  Simulation_det(10000,10,10,b0,d0,B,C,s,a,u,k)
    D_initial <-  Out_D[[RESIDENT]][10000]
    H_initial <-  1
    ESTABLISHMENT <- Out_D2[[INVADER]][10000]
    RES_EQUIL     <- Out_D2[[RESIDENT]][10000]
    
  }
  
  
  Estab_Probs <- rep(NA,Trials)
  Estab_Abuns <- rep(NA,Trials)
  
  fails <- 0 
  SUCCESS = "NO"
  
  for(TRIAL in 1:Trials){
    
    Hawks    <- rep(NA,Timesteps)
    Hawks[1] <- round(H_initial)
    
    Doves    <- rep(NA,Timesteps)
    Doves[1] <- round(D_initial)
    
    Abundances <- cbind(Hawks,Doves)
    if(SUCCESS=="YES") break
    
    for(Time in 2:Timesteps){
      New_Abundances <- Hawk_Dove_Change(Abundances[Time-1,1],Abundances[Time-1,2],b0,d0,B,C,s,a,u,k)
      Abundances[Time,]    <- c(New_Abundances[1],New_Abundances[2]) 
      
      
      if(New_Abundances[RESIDENT]==0) {
        Abundances[Time,RESIDENT] <- round(RES_EQUIL/2)
      }
      
      if(New_Abundances[INVADER]==0 &&  max(Abundances[,INVADER],na.rm=T) <  ESTABLISHMENT ){
        
        Estab_Probs[TRIAL] <- 0 
        break
      }
      
      if(New_Abundances[INVADER]==0 &&  max(Abundances[,INVADER],na.rm=T) >=  ESTABLISHMENT ){
        
        Estab_Probs[TRIAL] <- 1
        Estab_Abuns[TRIAL] <- sm_auc(1:Time, Abundances[,INVADER]  )
        break
      }
      
      
    }
    
    
    if(Estab_Probs[TRIAL]==1 && fails==5){
      ABUNS0 <- cbind(rep("SUCC",Time+10),Abundances[1:(Time+10),INVADER])  
      ABUNS <- rbind(ABUNS,ABUNS0 )
      SUCCESS <- "YES"
    }
    
    if(Estab_Probs[TRIAL]==0 && fails==4 && max(Abundances[,INVADER],na.rm=T)>2){
      fails <- fails + 1 
      ABUNS0 <- cbind(rep("fail_5",Time+150),Abundances[1:(Time+150),INVADER])  
      ABUNS <- rbind(ABUNS,ABUNS0 )
    }
    
    
    
    if(Estab_Probs[TRIAL]==0 && fails==3 && max(Abundances[,INVADER],na.rm=T)>2){
      fails <- fails + 1 
      ABUNS0 <- cbind(rep("fail_4",Time+150),Abundances[1:(Time+150),INVADER])  
      ABUNS <- rbind(ABUNS,ABUNS0 )
    }
    
    
    if(Estab_Probs[TRIAL]==0 && fails==2 && max(Abundances[,INVADER],na.rm=T)>2){
      fails <- fails + 1 
      ABUNS0 <- cbind(rep("fail_3",Time+150),Abundances[1:(Time+150),INVADER])  
      ABUNS <- rbind(ABUNS,ABUNS0 )
    }
    
    
    if(Estab_Probs[TRIAL]==0 && fails==1 && max(Abundances[,INVADER],na.rm=T)>2){
      fails <- fails + 1 
      ABUNS0 <- cbind(rep("fail_2",Time+150),Abundances[1:(Time+150),INVADER])  
      
      ABUNS <- rbind(ABUNS,ABUNS0 )
    }
    
    if(Estab_Probs[TRIAL]==0 && fails==0 &&max(Abundances[,INVADER],na.rm=T)>2  ){
      fails <- fails + 1 
      ABUNS <- cbind(rep("fail_1",Time+150),Abundances[1:(Time+150),INVADER])  
    }
    
    
    
  }    
  
  
  
  return(ABUNS)
  
}

set.seed(45)


TS_INV_Hawk<- Simulation_Inv_Estab_PLOT(Timesteps,Trials,H0,D0,b0,d0,B,C,s,a,u,k,"Hawk")
TS_INV_Hawk <- data.frame(TS_INV_Hawk)
TS_INV_Hawk$TIME <- seq(1:nrow(TS_INV_Hawk))
colnames(TS_INV_Hawk) <- c("Color","y","x")
TS_INV_Hawk$y <- as.numeric(TS_INV_Hawk$y)
TS_INV_Hawk_EST <- subset(TS_INV_Hawk,Color=="SUCC")
ZEROS <- 1:as.numeric(rownames(TS_INV_Hawk_EST)[1])
LENZ  <- as.numeric(rownames(TS_INV_Hawk_EST)[1])
DFBIND <- data.frame(cbind(rep("SUCC",LENZ), rep(-1,LENZ),ZEROS))
colnames(DFBIND) <- c("Color","y","x")
TS_INV_Hawk_EST <- rbind(DFBIND,TS_INV_Hawk_EST)
TS_INV_Hawk_EST$y <- as.numeric(TS_INV_Hawk_EST$y)
TS_INV_Hawk_EST$x <- as.integer(TS_INV_Hawk_EST$x)
TS_INV_Hawk <- head(TS_INV_Hawk,-11)
TS_INV_Hawk <- subset(TS_INV_Hawk,Color!="SUCC")

MEAN_ABUN     <- mean(TS_INV_Hawk_EST$y[which(TS_INV_Hawk_EST$y>-1)])
MIN_SOJ_time  <- nrow(TS_INV_Hawk)+1
MAX_SOJ_time  <-nrow(TS_INV_Hawk_EST)
MIDDLE <- MAX_SOJ_time - (MAX_SOJ_time  -MIN_SOJ_time)/2

Hawk <- TS_INV_Hawk_EST %>% ggplot(aes(x = x, y = y)) +
  geom_area(fill =  sm_color('skyblue'), alpha = 0.4) +
  
  geom_line(data = TS_INV_Hawk, aes(x, y), color = "grey65",size=1) +# Second plot overlaid
  geom_line(size = 1.1,  color = sm_color('skyblue')) +
  
  geom_segment(aes(x = MIN_SOJ_time, xend = MAX_SOJ_time, y = MEAN_ABUN, yend = MEAN_ABUN),
               linetype = "dashed", color = "black", size = 1.5)+
  
  geom_segment(aes(x = MIN_SOJ_time+10, xend = MIN_SOJ_time+10, y = 0, yend = MEAN_ABUN),
               linetype = "solid", color = "black", size = 1.5)+ 
  geom_segment(aes(x = MAX_SOJ_time-10, xend = MAX_SOJ_time-10, y = 0, yend = MEAN_ABUN),
               linetype = "solid", color = "black", size = 1.5)+  
  sm_hgrid(legends = FALSE) +
  theme(  panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
          
          plot.background = element_rect(fill = "white"),
          
          axis.text=element_text(size=16,face="bold",color="black"))  +
  
  scale_y_continuous( expand = c(0.00, .0),limits=c(0,23))+ 
  scale_x_continuous( expand = c(0.01, .0),limits=c(0,5150),breaks=c(0,1500,3000,4500))

ggsave(filename = "Hawk_TS.svg", plot = Hawk, device = "svg", width = 9, height =4, units = "in")
