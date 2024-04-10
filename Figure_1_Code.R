###############################################
############## Figure 1 #######################
##############################################

# will create Fig 1 in the main text 

########### Load packages ################ 
library(reshape2)
library(ggplot2)
library(dplyr)

# simple code that runs Wright-Fisher model 
WF_Sim <- function(s,Generations,N){
  
  Initial <- 1  # initial number 
  RUN     <- rep(0,Generations) # vector, number generations 
  RUN[1]  <- 1 # initial value set to 1 for output vector
  
  
  for(mm in 2:Generations){
    
    PROP_I   = RUN[mm-1]/N # Proportion of type i at time t, wi
    MEAN_FIT = PROP_I*(1+s) + (1-PROP_I) # mean fitness of population, wbar 
    PROB = PROP_I*((1+s)/MEAN_FIT) # perform normalization
    
    Num.New <- rbinom(1,N,PROB) # binomal sample 
    RUN[mm] <- Num.New # store output
  }
  return(RUN) # return vector of abundnaces 
}


# Set parameters 
s = .125 # selection coefficient
N = 20  # number of individuals 
Generations = 50 # number of generations 


# Death 1
set.seed(18) # sets seed (invading mutant go extinct)
Brach_1 <- WF_Sim(s,Generations,N) # run model
Brach_1 <- c(0,Brach_1) # add zero (for time=0)

# Death 2
set.seed(27) # sets seed (invading mutant go extinct)
Brach_2 <- WF_Sim(s,Generations,N) # run model
Brach_2 <- c(0,Brach_2) # add zero (for time=0)


# Death 3
set.seed(36) # sets seed (invading mutant go extinct)
Brach_3 <- WF_Sim(s,Generations,N) # run model
Brach_3 <- c(0,Brach_3) # add zero (for time=0)

# Fixes 
set.seed(34) # sets seed (invading mutant go go to fixation)
Brach_4 <- WF_Sim(s,Generations,N) # run model
Brach_4 
Brach_4 <- c(0,Brach_4) # add zero (for time=0)


Together <- c(Brach_2, Brach_1, Brach_3, Brach_4)# combine simulations  
df <- data.frame(Brach_2, Brach_1, Brach_3, Brach_4) # make as dataframe 

df_long <- tidyr::gather(df, key = "vector", value = "value") # data wrangling (will index in long form indexed with simulation 1-4)

#df_long <- tidyr::gather(df, key = "vector", value = "value")

df_long <- cbind(df_long,1:length(Together)) # add timesteps 
colnames(df_long) <- c("Color","y","x") # names 

df_long$y2 <- 20 - df_long$y

yyy <- ggplot() + geom_step(aes(x, y2),color="darkgrey", data = df_long,linewidth=1.1)+
  geom_step(color="lightgrey") + 
  theme(  panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank(), 
          plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=21,face="bold",color="black"),
          axis.line = element_line(colour = "grey65", 
                                   size = 1, linetype = "solid")) + 
  scale_y_continuous( expand = c(0.01, .01),limits = c(0,20),breaks=c(0,10 ,20)) +
  scale_x_continuous( expand = c(0.0135, 0.0135),breaks=c(0, 50,100,150,200)) 

xxx <-ggplot() + geom_step(aes(x, y, group = Color, color = Color), data = df_long,linewidth=1.1)+
  geom_step(aes(x, y2),color="darkgrey", data = df_long,linewidth=1.1)

xxx2 <- xxx + theme(  panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank(), 
                      plot.background = element_rect(fill = "white"),
                      axis.text=element_text(size=21,face="bold",color="black"),
                      axis.line = element_line(colour = "grey65", 
                                               size = 1, linetype = "solid")) + 
  scale_y_continuous( expand = c(0.01, .01),limits = c(0,21),breaks=c(0,10 ,20)) +
  scale_x_continuous( expand = c(0.0135, 0.0135),breaks=c(0, 50,100,150,200)) 

xxx3 <- xxx2 + geom_segment(aes(x = (length(c(Brach_1,Brach_2,Brach_3))+ 2*(log(s*N)+0.577216)/(s))+2, y = 0, xend = (length(c(Brach_1,Brach_2,Brach_3))+ 2*(log(s*N)+0.577216)/(s) )+2, yend = N+.5), linetype = "dashed",linewidth=1)
xxx4 <- xxx3 + geom_segment(aes(x = (length(c(Brach_1,Brach_2,Brach_3)))+1, y = 0, xend = (length(c(Brach_1,Brach_2,Brach_3))), yend = N+.5), linetype = "dashed",linewidth=1)

ggsave(filename = "Fig_1_Final.svg", plot = xxx4, device = "svg", width = 8*1.5, height = 5*1.5, units = "in") # save to working directory

# The figure was finished in PPT, where labels were added 

