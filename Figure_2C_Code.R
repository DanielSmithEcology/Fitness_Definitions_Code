#####################################################
################### Figure 2C Code   ################
#####################################################

# This code will run simulations that reproduce Fig. 2C in the main text

### load package(s) ###
library(expm)

########################### 
####### Functions ######### 
########################### 

# fec.Func: outputs either fecundity of a good year or of a bad year
# f_bad: fecundity on bad year, f_good: fecundity on good year, p_good: probability of a good yearfec.Func <- function(f_bad,f_good,p.good){
fec.Func <- function(f_bad,f_good,p.good){
  Year <- runif(1) # random value 
  if(Year<=p.good) return(f_good) # return if good year happens
  if(Year>p.good) return(f_bad) # return if bad year happens 
}


######### Function that will run Adaptive Dynamics simulations ############

##### Description of function / simulations: 
# Simulations will run simulation of the annual plant model
# In each simulation a resident with germination probability g_R competes two mutants invaders with germination probabilities gR +/- DELTA
# Simulations are run until one of the types fixes (either the resident or one of the mutants)
# Simulations are run for (approximately) values of gR for  0 < gR < 1
# each simulation is run  "Sim_Trials" times; each simulation records the "winner"
# the simulations are then used to derive the probability gR in not invaded or if one of the mutants fixes
# these probabilities are then used to generate a transition matrix for the value of gR in the system 

#### Definition of parameters: #### 
# Delta: increment germination probability; invaders germination probability gR +/- DELTA; .01 in simulations
# f_bad: fecundity on bad year
# f_good: fecundity on good year
# p.bad: probability of good year
# s: probability seed dies in seed bank
# a: density-dependent term 
# Sto: set to either "Yes" or "No"; "Yes" adds demographic stochasticity to the model; "No" makes the model deterministic
# Sim_Trials: Number of simulations to run; only relevant for when there is demographic stochsticity

Prob_Sim <- function(DELTA,f_bad,f_good,p.good,s,a,Sto,Sim_Trials){
  
  if(Sto=="Yes"){ # if "Yes" runs the model with demographic stochsticity
    
    Annual_Model_Comp_3_Types <- function(N1,N2,N3,s,g1,g2,g3,f,a){
      
      GER1 <-    rbinom(1,N1,g1) # random germination type 1
      GER2 <-    rbinom(1,N2,g2) # random germination type 2
      GER3 <-    rbinom(1,N3,g3) # random germination type 2
      
      EMER1 <-  rbinom(1,GER1,1/(1+(GER1+GER2+GER3)*a)) # Density dependent survival type 1
      EMER2 <-  rbinom(1,GER2,1/(1+(GER1+GER2+GER3)*a)) # Density dependent survival type 2
      EMER3 <-  rbinom(1,GER3,1/(1+(GER1+GER2+GER3)*a)) # Density dependent survival type 3
      
      Fec1 <-   rpois(1,f*EMER1) # Poisson seed production type 1 
      Fec2 <-   rpois(1,f*EMER2) # Poisson seed production type 2
      Fec3 <-   rpois(1,f*EMER3) # Poisson seed production type 3
      
      Sur1 <- rbinom(1,N1-GER1,s) # binom survival of non-germinated seeds type 1
      Sur2 <- rbinom(1,N2-GER2,s) # binom survival of non-germinated seeds type 2
      Sur3 <- rbinom(1,N3-GER3,s) # binom survival of non-germinated seeds type 3
      
      N1_new <- Fec1 + Sur1 # seeds at next timestep type 1
      N2_new <- Fec2 + Sur2 # seeds at next timestep type 2
      N3_new <- Fec3 + Sur3 # seeds at next timestep type 3
      
      return(c(N1_new,N2_new,N3_new ))
      
    }
    
    Annual.Model <- function(N,s,g,f,a){ # change per year (one type)
      
      GER <-   rbinom(1,N,g) # random germinatinon
      
      EMER <-  rbinom(1,GER,1/(1+(GER)*a))
      
      Fec <-   rpois(1,f*EMER) # Poisson seed production 
      
      Sur <- rbinom(1,N-GER,s) # binomial survival of non germinated seeds in bank
      
      N_new <- Fec + Sur # seeds next time step
      
      return(N_new)
      
      # return(rpois(1,lambda))
    }
    
    
  }
  
  
  if(Sto=="No"){ # if "No" runs the model with demographic stochsticity
    
    
    Annual_Model_Comp_3_Types <- function(N1,N2,N3,s,g1,g2,g3,f,a){ 
      
      lambda1 <- s*(1-g1)*N1 + (g1*f*N1)/(1+(g1*N1+g2*N2+g3*N3)*a)
      lambda2 <- s*(1-g2)*N2 + (g2*f*N2)/(1+(g1*N1+g2*N2+g3*N3)*a)
      lambda3 <- s*(1-g3)*N3 + (g3*f*N3)/(1+(g1*N1+g2*N2+g3*N3)*a)
      
      return(  c(lambda1,lambda2,lambda3) )
      
    }
    
    
    Annual.Model <- function(N,s,g,f,a){ # Change per year (one type)
      s*(1-g)*N + (g*f*N)/(1+g*N*a)
    }
    
    
    
  }
  
  # Initial values 
  g_min <- .0095 # smallest value of gR in model
  g_Mid <- .0095 # "g_Mid" is dummy variable that indexes which type has the "middle" germination rate (the resident). 
  Max_g <- .995 # maximum value
  
  LEN <- length(seq(g_min,Max_g,by=DELTA)) # calculates the number of values of gR that will be run 
  
  Trans_Mat <- matrix(0,nrow=LEN,ncol=LEN) # sets up matrix for transition matrix 
  
  IND <- 1 # start index (which keeps track of the nth value of gR considered)
  
  while(IND<=LEN){
    g_Vals <- c(g_Mid-DELTA,g_Mid,g_Mid+DELTA) # germination values for simulation (gR +/- DELTA)
    Counts <- c(0,0,0) # vector that will keep track of each time gR - DELTA, gR, and gR+Delta fixes / "wins" simulation
    
    print(g_Mid)
    
    for(yy in 1:Sim_Trials){ # loop over the number of simulations 
      
      Burnin <- 125 # time of "burnin" for the resident 
      N_R.Burn <- rep(0,Burnin) # vector for resident 
      N_R.Burn[1] <- 25 # initial abundance (arbitrary)
      
      for(bb in 2:Burnin){ # run burn-in
        FEC <- fec.Func(f_bad,f_good,p.good) # fecundity at time bb
        N_R.Burn[bb] <- Annual.Model(N_R.Burn[bb-1],s,g_Mid,FEC,a) # store abundance
        if(N_R.Burn[bb]==0) N_R.Burn[bb] <- 1 # reflecting boundary 
      }
    
      
      N_R_Init <- N_R.Burn[Burnin] # Set initial abundance of resident 
      Abuns_mm0 <- c(1,N_R_Init,1) # set initial abundances of mutant invaders to 1 (may be modified below)
      
      if(IND==1){ # if 1st simulation (gR=.0095), we assume there is only 1 invading mutant with germination probability gR+DELTA
        
        g_Vals[1] <- g_min 
        Abuns_mm0 <- c(0,N_R_Init,1) # adjust initial abundance accordingly
        
      }
      
      if(IND==LEN){ # if final simulation, we assume there is only 1 invading mutant with germination probability gR-DELTA
        
        g_Vals[3] <- Max_g
        Abuns_mm0 <- c(1,N_R_Init,0) # adjust initial abundance accordingly
        
      }
      
      Abuns_mm0[2] <- max(1,Abuns_mm0[2] - sum(Abuns_mm0[1]+Abuns_mm0[3]) ) # adjust of resident to account for mutation  
      
      Abuns_mm <- Abuns_mm0 # Set initial abundance 
      
      
      while(length(which(Abuns_mm>=.01))>1){ # while there is more than 1 type in population 
        
        FEC      <- fec.Func(f_bad,f_good,p.good) # set fecundity (good vs bad year)
        
        Abuns_mm <- Annual_Model_Comp_3_Types(Abuns_mm[1],Abuns_mm[2],Abuns_mm[3],s,g_Vals[1],g_Vals[2],g_Vals[3],FEC,a)  # Change in population abundances 
        
        #print(Abuns_mm)
        
        if(length(which(Abuns_mm<.01))==2){ # if there is only one type left
          
          Index      <- which(Abuns_mm >= .01) # record winner
          
          Counts[Index] <- Counts[Index] + 1 # count victory
          
        }
        
        if(length(which(Abuns_mm<.01))==3){ # if all 3 happen to go extinct on same timestep 
          
          Abuns_mm <- Abuns_mm0 # start simulation over
        }
        
      }
    }
    
    print(Counts)  
    
    Probs <- Counts/sum(Counts)
    
    print(g_Vals)
    
    print(Probs)
    
    
    
    if(IND!=1 & IND!=LEN ){ # if not 1st or last simulation 
      
      Trans_Mat[c(IND-1,IND,IND+1),IND] <- Probs # input transition probs to matrix
      
    }
    
    
    if(IND==1){ # if 1st simulation (with only one mutant invader with gR + Delta)
      
      Trans_Mat[c(IND,IND+1),IND] <- Probs[2:3] # input transition probs to matrix
      
    }
    
    if(IND==LEN){ # if last simulation (with only one mutant invader with gR - Delta)
      
      Trans_Mat[c(IND-1,IND),IND] <- Probs[1:2] # input transition probs to matrix
      
    }
    
    g_Mid <- g_Mid+DELTA
    IND <- IND + 1
    print(IND)
    
  }
  
  return(Trans_Mat) # return matrix 
  
}

###################
####Parameters#####
###################

a <- .075 # density dependence 
s <- .95  # Survival between years 
f_good <- 2.25 # fecundity of good year
f_bad  <- 0 # fecundity of bad year 
p.good <- .95 # probability of good year
DELTA <- .01 # increment of germination change (germination rate changes 1% between competitors)
Number_Sims_Sto <- 50000 # number of simulations for the model with demographic stochasticity
Number_Sims_Det <- 1 # number of simulations for the model without demographic stochastic; only one repetition is needed because it is deterministic


set.seed(284) # set seed 
Stochastic_Sims    <- Prob_Sim(DELTA,f_bad,f_good,p.good,s,a,"Yes",50000) # run simulation with demographic stochasticity; 
Deterministic_Sims <- Prob_Sim(DELTA,f_bad,f_good,p.good,s,a,"No",1)# run simulation without demographic stochasticity

## Find stationary distributions (done by taking the transition matrix to the nth power as n gets very large; identical to finding eigenvalue, but faster)
Stationary_with_demographic_sto <- (t(Stochastic_Sims) %^% 10000000)[1,] # 
Stationary_no_demographic_sto   <- (t(Deterministic_Sims) %^% 10000000)[1,] # will return ESS

X_axis                          <- seq(.0095,.995,by=DELTA) # make x axis

############################
####### Make Figure ########
############################

svg("Fig_2_C.svg", width = 7/1.5, height = 8.5/1.5)
par(mar = c(6, 5, 5, 4))
plot(Stationary_with_demographic_sto~X_axis,ylim=c(0.0003,.0235),xlim=c(0,1),pch = 21,ylab="",xlab="", cex=1.75, type = "p",bg="lightblue",col="black",yaxt="n",xaxt="n")
axis(2,at=c(0,.01,.02,.03),col.axis = 'black',cex.axis=1.5,las=1)
axis(1,at=c(0,.5,1),col.axis = 'black',cex.axis=1.5)
abline(v=X_axis[which.max(Stationary_no_demographic_sto)],col='red',lwd=3,lty=2)
dev.off()









