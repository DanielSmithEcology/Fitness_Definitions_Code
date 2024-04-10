#####################################################
################### Figure 2 Code ###################
#####################################################


########################### 
####### Functions ######### 
########################### 


#Annual.Model: will return a number of seeds form year t to t+1 a single type given parameters
# N: number of seeds, s survival proportion in seed bank (1/d in main text), g germination prob, f fecundity, a density-dependent parameter (alpha in main text)
Annual.Model <- function(N,s,g,f,a){
  
  GER <-   rbinom(1,N,g) # random germinatinon
  
  Fec <-   rpois(1,(f*GER)/(1+(GER)*a)) # poisson seed production
  
  Sur <- rbinom(1,N-GER,s) # binomial survival of non germinated seeds in bank
  
  N_new <- Fec + Sur # seeds next time step
  
  return(N_new)
  
  # return(rpois(1,lambda))
}

# Annual.Model.Comp: will return a number of seeds form year t to t+1 a two types given parameters
# parameters are same as for the above function, with two types
Annual.Model.Comp <- function(N1,N2,s,g1,g2,f,a){
  
  GER1 <-   rbinom(1,N1,g1) # random germination type 1
  GER2 <-   rbinom(1,N2,g2) # random germination type 2
  
  Fec1 <-   rpois(1,(f*GER1)/(1+(GER1+GER2)*a)) # Poisson seed production type 1 
  Fec2 <-   rpois(1,(f*GER2)/(1+(GER1+GER2)*a)) # Poisson seed production type 2
  
  Sur1 <- rbinom(1,N1-GER1,s) # binom survival of non-germinated seeds type 1
  Sur2 <- rbinom(1,N2-GER2,s) # binom survival of non-germinated seeds type 2
  
  N1_new <- Fec1 + Sur1 # seeds at next timestep type 1
  N2_new <- Fec2 + Sur2 # seeds at next timestep type 2
  
  return(list(N1_new,N2_new ))
  
}

# fec.Func: outputs either fecundity of a good year or of a bad year
# f_bad: fecundity on bad year, f_good: fecundity on good year, p_good: probability of a good year
fec.Func <- function(f_bad,f_good,p.good){
  Year <- runif(1) # random value 
  if(Year<=p.good) return(f_good) # return if good year happens
  if(Year>p.good) return(f_bad) # return if bad year happens 
}

# Growth.Invader.TIME: given the abundance of a resident, return growth rate of invader  
# gR resident germination probability, gI invader germination probability 
# Note: log of this function is the malthusian parameter
Growth.Invader.TIME <- function(N_Resident,s,f,a,gR,gI){
  
  Nstar <- N_Resident
  INV   <- s*(1-gI) +  (gI*f)/(1+(gR*Nstar+gI)*a) 
  return(INV)
}

##########################
###### Simulations #######
##########################

##### Set parameters #####

Burnin <- 1500 # timesteps for a "burn in" of a resident 
N_R.Burn <- rep(NA,Burnin) # make vector for resident 
N_R.Burn[1] <- 50 # initial abunance of resident 
# Burn in the resident
f_good <- 3 # fecundity of good year
f_bad  <- 0 # fecundity of bad year 
p.good <- .95 # probability of good year
s0 <- .95 # probability of survival in seed bank
gR <- .2 # resident germination probability 
a0 <- .025 # density--dependent parameter 


#### run burn in of resident ####

set.seed(385) # set seed for simulation 
for(bb in 2:Burnin){ # run burn-in
  FEC <- fec.Func(f_bad,f_good,p.good) # fecundity at time bb
  N_R.Burn[bb] <- Annual.Model(N_R.Burn[bb-1],s0,gR,FEC,a0) # store abundance 
}


#### Invasion probability with resident ####

set.seed(6345) # set seed
gI.Vec2  <- c(seq(.05,.999999,length=100 )) # vector of gI values 
gI.Sims2 <- length(gI.Vec2) # number of sims 
gI.FIX <- rep(NA,gI.Sims2) # empty vector, store invasion probs
Num.Trials  <- 50000 # number of simulations at each value of gI

# run the simulation #

for(xx in 1:gI.Sims2){ # for all gI
  gI <- gI.Vec2[xx] # set gI 
  
  fix.vec <- rep(NA,Num.Trials) # empty vector for storing success (=1) or failure (=0)
  for(ee in 1:Num.Trials) { # for all trials
    N_R <- round(mean(N_R.Burn)) # start resident abundance at around mean of burn in (rounded to integer)
    N_I <- 3 # start with 3 invaders (makes for less noisy output that 1 invader; qualitatively the same)
    while(N_I>0){ # while invader has not gone extinct... 
      FEC <- fec.Func(f_bad,f_good,p.good) # new fecundity value 
      OUTS <- Annual.Model.Comp(N_R,N_I,s0,gR,gI,FEC,a0) # output new abundances of resident and invader
      N_R <- OUTS[[1]] # store resident abundance
      N_I <- OUTS[[2]] # store invader abundance
      
      if(N_I<1) { # if abundance of invader less than 1 (extinct)
        fix.vec[ee] <- 0 # count as extinct 
        break   # end sim
      }
      
      if(N_I>=1 && N_R<1)  { # if the invader has not gone extinct and the resident has gone extinct (fixation)
        fix.vec[ee] <- 1 # count as fixed
        break #end sim
      }
      
    }
    
  }
  
  gI.FIX[xx] <- mean(fix.vec,na.rm=T) # probability of invasion for simulation with gI is given by the mean of fix.vec (1s=fixed, 0s = extinct)
}



#### Invasion probability with NO resident ####

# similar to above, but N_R is set to zero

set.seed(6345) # set seeed 
Time.Sim2   <- 20 # we examine what happens after 20 generations 
gI.Vec2b  <- c(seq(.05,.93,length=100 ),log(seq(exp(.9307),exp(.9999999),length=23))) #vector of gI values examined (thes vales give a new looking line)
gI.Sims2b <- length(gI.Vec2b) # number of sims 
gI.FIX2 <- rep(NA,gI.Sims2b)  # empty vector for counting fixation
Num.Trials  <- 50000 #  number of simulations at each value of gI


for(xx in 1:gI.Sims2b){ # for all gI 
  gI <- gI.Vec2b[xx] # set gI
  N_R    <-0 # set resident to zero abundance 
  N_I    <- rep(NA,Time.Sim2) # empty vector for invader abundance 
  N_I[1] <- 3 # initial abundance (set to 3 to reduce noisiness; =1 makes no qualitative change)
  fix.vec <- rep(NA,Num.Trials) # empty vector to count fixation (=1) or extinction (=0)
  for(ee in 1:Num.Trials) { # for all trials 
    for(dd in 2:Time.Sim2){ # run simulation over timesteps 
      FEC <- fec.Func(f_bad,f_good,p.good) # set fundity 
      N_I[dd] <- Annual.Model.Comp(N_I[dd-1],N_R,s0,gI,gR,FEC,a0)[[1]] # return invader abundance 
    }
    if(N_I[Time.Sim2]<1) fix.vec[ee] <- 0 # if extinct, count as 0
    
    if(N_I[Time.Sim2]>=1)  fix.vec[ee] <- 1 # if not extinct, count as 1
    

  }
  gI.FIX2[xx] <- mean(fix.vec,na.rm=T) # probability of invasion for simulation with gI is given by the mean of fix.vec (1s=fixed, 0s = extinct)
}


#### geometric mean growth rate when invader competes against a resident type #### 

gI.Vec  <- log(log(c(seq(exp(exp(.05)),exp(exp(.999999)),length=3500)))) # vector of gI values 
gI.Sims <- length(gI.Vec) # number of sims 
gI.INVX  <- rep(NA,gI.Sims) # empty vector for storing gemometric mean growth rates 

N_R <- rep(NA,Time.Sim) # empty vector for resident abundances 
N_R[1] <- round(mean(N_R.Burn)) # start resident at mean of burn-in from before 

Fec_List <- rep(NA,Time.Sim-1) # empty vector, will store fecundity values 

# run simulation to give a time series of the resident species and store fecundity values  
for(cc in 2:Time.Sim){
  N_I <- 0 # no invader present 
  FEC <- fec.Func(f_bad,f_good,p.good) # get fucundity value 
  Fec_List[cc-1] <- FEC # store fecundity value 
  OUTS <- Annual.Model.Comp(N_R,N_I,s0,gR,gI,FEC,a0) # output new abundances
  N_R[cc] <- OUTS[[1]] # store resident abundance 
  if(N_R[cc]<1) N_R[cc] <-  N_R[1] # reflecting boundary; if resident goes extinct, set it back to mean (this almost never happens with the parameters)
}

## run simulation generating gemometric mean growth rate of the invader ## 
for(xx in 1:gI.Sims){ # for all gI values 
  gI <- gI.Vec[xx]
  GeoMean <- rep(NA,Time.Sim-1) # empty vector to store geometric mean
  
  for(cc in 2:Time.Sim){ # for timesteps 
    FEC <- Fec_List[cc-1] # set fecundity 
    GeoMean[cc-1] <- log(Growth.Invader.TIME(N_R[cc-1],s0,FEC,a0,gR,gI)) # calculate invader growth rate (need to take log)
  }
  gI.INVX[xx] <- mean(GeoMean) # return mean growth rate 
  
}



### Invader growth rate when there is no resident 
set.seed(6345) # set seed 
gI.Vecb  <- c(log(log(c(seq(exp(exp(.05)),exp(exp(.9995)),length=3500)))),seq(.9995,.999999,by=.00001)) # set vector of gI values 
gI.Sims <- length(gI.Vecx) # number of sims 
gI.INV  <- rep(NA,gI.Sims) # empty vector of invader growth rates  
N_R <- rep(0,Time.Sim) # set resident abundance to zero

Fec_List <- rep(NA,Time.Sim-1) # fecundity list 

## set fecundity values ## 
for(cc in 2:Time.Sim){ # over time series 
  FEC <- fec.Func(f_bad,f_good,p.good) # get fecundity 
  Fec_List[cc-1] <- FEC # store fecundity 
}

# calculate geometric mean growth rate 
for(xx in 1:gI.Sims){ # for timesteps 
  gI <- gI.Vecx[xx] # set gI value 
  Lyp <- rep(NA,Time.Sim-1)
  
  for(cc in 2:Time.Sim){ # for timestep (20)
    FEC <- Fec_List[cc-1] # set fecundtiy 
    GeoMean[cc-1] <- log(Growth.Invader.TIME(N_R[cc-1],s0,FEC,a0,gR,gI)) # output growth rate / store 
  }
  gI.INV[xx] <- mean(GeoMean) # return mean growth rate 
  
}

#############################
####### Plot Figures ########
#############################

# plot function (makes plotting easier)
myPlot <- function(x,y, yaxt = "n", ylab = NA) {
  plot(x, y, yaxt = yaxt,  pch = 19, cex=1, col="blue", bg="blue", lwd=2,type="p",xlim=c(0,1),cex.axis=1.5,ylim=c(.015,1),las=1,ylab="",xlab="",xaxt="n")
  axis(4,cex.axis=1.5,las=1,col.axis = "blue")
  mtext(ylab, 4, line = 2,cex.axis=1.5,las=1,col.axis = "blue")
}

# make and save plot 1 

svg("Fig_2_With_Res.svg", width = 7/1.5, height = 8.5/1.5)
# set margins 
par(mar = c(6, 5, 5, 4))
# plot 1
plot(gI.INVX~gI.Vec, pch = 19,ylab="",xlab="", cex=1, type = "p",bg="black",col="darkgreen",xlim=c(0,1),cex.axis=1.5,ylim=c(gI.INVX[1],max(gI.INVX)*1.5),las=1,col.axis = 'darkgreen',xaxt="n")
axis(1,at=c(0,.5,1),col.axis = 'black',cex.axis=1.5)
par(new=T)
myPlot(gI.Vec2, gI.FIX)
dev.off()


# make and save plot 2 
svg("Fig_2_Without_Res.svg", width = 7/1.5, height = 8.5/1.5)
par(mar = c(6, 5, 5, 4))
plot(gI.INV~gI.Vecx, pch = 19,ylab="",xlab="", cex=1, type = "p",bg="black",col="darkgreen",xlim=c(0,1),cex.axis=1.5,ylim=c(gI.INV[1],max(gI.INV)*1.05),las=1,col.axis = 'darkgreen',xaxt="n")
axis(1,at=c(0,.5,1),col.axis = 'black',cex.axis=1.5)
par(new=T)
myPlot(gI.Vec2b, gI.FIX2)
dev.off()







