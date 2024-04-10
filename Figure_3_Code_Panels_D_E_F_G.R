######################################################
################# Fig 3, Panels D-G ##################
######################################################


# This will reproduce Fig 3, Panels D, E, F, G



############################
####### Load Packages ######
############################
library(tidyverse)
library(smplot2)
library(cowplot)
library(ggbreak) 
library(patchwork)
library(ggh4x)

############################
##### Define functions #####
############################

# Hawk_Dove_Change: returns an single timepstep of the stochastic Hawk-Dove game 

### Parameters###
#D number of Doves, H number of Hawks, b0 baseline births per timestep, d0 mortality probability per timestep
# B benefit of getting resource, C Cost of fighting, s strength of selection, a strength of density-dependence
# u mutation rate (not used in any simulations), k (a dummy parameter; never used, but was important for testing; kindly ignore it)
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

# Simulation_Inv_Estab: inputs a number of timesteps, Trials, initial conditions, parameters, and specifies an "Invader" as
# "Hawk" or "Dove"; will return the fitness metrics of a balanced polymorphism discussed in the main text 
Simulation_Inv_Estab <- function(Timesteps,Trials,H0,D0,b0,d0,B,C,s,a,u,k,Invader){
  
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

  # empty vector of each fitness quantity 
  Estab_Probs  <- rep(NA,Trials)
  Estab_Abuns  <- rep(NA,Trials)
  Estab_SOJ    <- rep(NA,Trials)
  Estab_Height <- rep(NA,Trials)
    
  for(TRIAL in 1:Trials){
    
    Hawks    <- rep(NA,Timesteps)
    Hawks[1] <- round(H_initial) # initial abundance of Hawk
    
    Doves    <- rep(NA,Timesteps)
    Doves[1] <- round(D_initial) # initial abundance of Dove 
    
    Abundances <- cbind(Hawks,Doves)
  
    Time <-1 # set initial Time = 1 
    
  while(Abundances[Time,INVADER]>0){ # if invader is not extinct at time TIME 
    Time <- Time +1  # update time 
    New_Abundances <- Hawk_Dove_Change(Abundances[Time-1,1],Abundances[Time-1,2],b0,d0,B,C,s,a,u,k) # calculate abundances of Hawk, Dove
    Abundances[Time,]    <- c(New_Abundances[1],New_Abundances[2]) # save abundances of each type in DF
    

    if(New_Abundances[RESIDENT]==0) { # if resident goes extinct 
     Abundances[Time,RESIDENT] <- 1 # reflecting boundary (so doesn't stay extinct)
   }
    
    if(New_Abundances[INVADER]==0 &&  max(Abundances[,INVADER],na.rm=T) <  ESTABLISHMENT ){ # if invader goes extinct before establishment
      
      Estab_Probs[TRIAL] <- 0 # count invasion as a "failure" 
      Estab_Abuns[TRIAL]  <- sm_auc(1:Time, Abundances[,INVADER]  ) # store time integral 
      
      break
    }
    
    if(New_Abundances[INVADER]==0 &&  max(Abundances[,INVADER],na.rm=T) >=  ESTABLISHMENT ){ # if invader goes extinct after succesful establishment
      
      Estab_Probs[TRIAL]  <- 1 #count invasion as a success
      Estab_Abuns[TRIAL]  <- sm_auc(1:Time, Abundances[,INVADER]  ) # store time integral
      Estab_SOJ[TRIAL]    <- Time # store sojurn time from introduction to extinction 
      Estab_Height[TRIAL] <- mean(Abundances[,INVADER],na.rm=T) # Store mean abundance over sojurn 
      
      break
    }
    
    
    }

  }    
  # take means of each fitness-related quantity over simulations 
  Mean_Est_Prob   <- mean(Estab_Probs,na.rm=T)
  Mean_Abun       <- mean(Estab_Abuns,na.rm=T)
  Mean_SOJ        <- mean(Estab_SOJ,na.rm=T)
  Mean_Height     <- mean(Estab_Height,na.rm=T)
  
  return(list(Mean_Est_Prob,Mean_Abun,Mean_SOJ,Mean_Height)) # return 
    
}


############################
##### Set parameters #######
############################
Timesteps <- 1000000
Trials <- 7500
H0 <- 10 #initial Hawk
D0 <- 10 # Initial dove 
b0 <- .05 # birth rate
d0 <- .01 # Death prob
B <- 2 # Benefit 
s <- .5 # strength of selection (imapct of hawk-Hove encounters)
a <- .15 # strength of density-dependence 
u <- 0.0000 # mutation rate (set to zero)
k=1 # dummy parameter used for testing; ignore 


############################
##### Run Simulations ######
############################

LENNY <- 11 # number of simulations 
C_Vals <- seq(2.25,5.75,length=LENNY) # values of C (Cost of fighting) considered 

# Make empty vectors for each quantity 
DOVE_ABUN    <- rep(NA,LENNY)
DOVE_PROB    <- rep(NA,LENNY)
DOVE_SOJ     <- rep(NA,LENNY)
DOVE_HEIGHT  <- rep(NA,LENNY)

HAWK_ABUN    <- rep(NA,LENNY)
HAWK_PROB    <- rep(NA,LENNY)
HAWK_SOJ     <- rep(NA,LENNY)
HAWK_HEIGHT  <- rep(NA,LENNY)

set.seed(3506) # set seed 
# Run for loop for all values of C 
for(nn in 1:LENNY){
  
  C0 <- C_Vals[nn]
  
  SIM_Dove <- Simulation_Inv_Estab(Timesteps,Trials,H0,D0,b0,d0,B,C0,s,a,u,k,"Dove")
  SIM_Hawk <- Simulation_Inv_Estab(Timesteps,Trials,H0,D0,b0,d0,B,C0,s,a,u,k,"Hawk")
  
  DOVE_PROB[nn]   <- SIM_Dove[[1]]
  DOVE_ABUN[nn]   <- SIM_Dove[[2]]
  DOVE_SOJ[nn]    <- SIM_Dove[[3]]
  DOVE_HEIGHT[nn] <- SIM_Dove[[4]]
  
  HAWK_PROB[nn]  <- SIM_Hawk[[1]]
  HAWK_ABUN[nn]   <- SIM_Hawk[[2]]
  HAWK_SOJ[nn]    <- SIM_Hawk[[3]]
  HAWK_HEIGHT[nn] <- SIM_Hawk[[4]]  
    
}


# save the ratios of each fitness component for plotting 
PROB_RATIO    <- HAWK_PROB/DOVE_PROB 
SOJ_RATIO     <- HAWK_SOJ/DOVE_SOJ
HEIGHT_RATIO  <- HAWK_HEIGHT/DOVE_HEIGHT
Total         <- HAWK_ABUN/DOVE_ABUN



######### Make the figure, save in current working directory ##############

DIFF <- C_Vals/B # vector of C/B
DF <- data.frame(cbind(PROB_RATIO,DIFF))

Prob_Plot <- ggplot(DF, aes(DIFF, PROB_RATIO)) +
  geom_point(size=3.5) + # Add points
  sm_hgrid(legends = FALSE) +
  theme(  panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
              axis.title.x=element_blank(),
             axis.title.y=element_blank(), 
          plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=24,face="bold",color="black")) +

  scale_y_continuous(trans = 'log10',limits = c(.07,26),breaks=c(0.1,1,10))+
  scale_x_continuous( breaks=c(1.25,2,2.75),limits=c(1.1,2.9)) 


ggsave(filename = "Prob_Plot.svg", plot = Prob_Plot, device = "svg", width = 4, height = 5, units = "in")

DF <- data.frame(cbind(SOJ_RATIO,DIFF))

Soj_Plot <- ggplot(DF, aes(DIFF, SOJ_RATIO)) +
  geom_point(size=3.5) + 
  sm_hgrid(legends = FALSE) +
  theme(  panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
              axis.title.x=element_blank(),
            axis.title.y=element_blank(), 
          plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=24,face="bold",color="black")) +
  scale_y_continuous(trans = 'log10',limits = c(.07,26),breaks=c(0.1,1,10))+
  scale_x_continuous( breaks=c(1.25,2,2.75),limits=c(1.1,2.9)) 

ggsave(filename = "Soj_Plot.svg", plot = Soj_Plot, device = "svg", width = 4, height = 5, units = "in")


DF <- data.frame(cbind(HEIGHT_RATIO,DIFF))


HeightPlot <- ggplot(DF, aes(DIFF, HEIGHT_RATIO)) +
  geom_point(size=3.5) + 
  sm_hgrid(legends = FALSE) +
  theme(  panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
              axis.title.x=element_blank(),
             axis.title.y=element_blank(), 
          plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=24,face="bold",color="black")) +
  scale_y_continuous(trans = 'log10',limits = c(.07,26),breaks=c(0.1,1,10))+
  scale_x_continuous( breaks=c(1.25,2,2.75),limits=c(1.1,2.9)) 

ggsave(filename = "HeightPlot.svg", plot = HeightPlot, device = "svg", width = 4, height = 5, units = "in")


DF <- data.frame(cbind(Total,DIFF))
Total2_Plot <- ggplot(DF, aes(DIFF, Total)) +
  geom_point(size=3.5) + 
  sm_hgrid(legends = FALSE) +
  theme(  panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
              axis.title.x=element_blank(),
             axis.title.y=element_blank(), 
          plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=24,face="bold",color="black")) +
  scale_y_continuous(trans = 'log10',limits = c(.07,26),breaks=c(0.1,1,10))+
  scale_x_continuous( breaks=c(1.25,2,2.75),limits=c(1.1,2.9)) 

ggsave(filename = "Total_Plot.svg", plot = Total2_Plot, device = "svg", width = 4, height = 5, units = "in")


# figures were finished in PowerPoint (arranged and axes labels added)




