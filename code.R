install.packages("deSolve")
install.packages("pso")
install.packages("ggplot2")

library(deSolve)
library(pso)
library(ggplot2)

################### PLOTTING DATA ######################

#1) Creating the infected cases plot (cumulative)

#First value is 3285 because there were 3285 dengue cases in 2018
#Rest of the values are data collected online
data <- (c(3285,205,245,207,221,179,135,232,181,156,
           133,109,96,99,108,125,134,156,192,287,307,
           376,400,396,464,430,492,591,661,644,601,597,
           521,522,474,412,326,314,299,258,243,226,239,
           241,306,321,372,329,295,280,257,226,290)
         /10000) #rescale by ten thousand

data_I <- cumsum(data)

#infected cases plot
time <- seq(0, 52, by = 1)
plot(time, data_I,main = "Infected Cases over Time", xlab = "Number of weeks", ylab ="Number of Infected cases (Ten Thousand)", col=2, pch=16,ylim=c(0,2))

#2) Creating the susceptible cases plot

#numbers will be in ten thousands, e.g. 570 because population in 2019 is 5.7 million
#number of deaths in 2019 is 21385, assume constant death rate
#subtract number of deaths and infected from population
data_S <- 570 - (0:52) * 21385 / 52 / 10000 - data_I 
data_S

#susceptible cases plot
plot(time, data_S,main = "Susceptible Cases over Time", xlab = "Number of weeks",
     ylab ="Number of Susceptible cases (Ten Thousand)", col=2, pch=16,ylim=c(565,570))


#################### DEFINING ODE ######################
RHS.F <- function(S,I,M,d,i,r,v) {
  -d*S-v*S-i*I*S 
}

RHS.G <- function(S,I,M,d,i,r,v) {
  i*I*S-d*I-r*I
}

RHS.H <- function(S,I,M,d,i,r,v){
  v*S+r*I-d*M
}

RHS <- function(t, state, parameters) {
  with(as.list(c(state,parameters)),{
    # rate of change
    dS <- RHS.F(S,I,M,d,i,r,v)
    dI <- RHS.G(S,I,M,d,i,r,v)
    dM <- RHS.H(S,I,M,d,i,r,v)
    list(c(dS,dI,dM))
  })   # end with(as.list ...
}


############### CREATING ERROR FUNCTION ################

#Error function
error_f <- function(p){
  p <- c(d=p[1],i=p[2],r=p[3],v=p[4]) #name the variables
  #values of t
  time <- seq(0, 52, by = 1)
  #initial values, when t=0
  state <- c(S=569.6715, I=0.3285, M=0)
  #run numerical simulation
  out <- ode(y = state, time = time, func = RHS, parms = p)
  #extracting as a data frame
  out.df <- as.data.frame(out)
  #extract out S and I
  S <- out.df[2]
  I <- out.df[3]
  #change from dataframe to vector
  S <- as.vector(t(S))
  I <- as.vector(t(I))
  model <- c(S,I)
  data <- c(data_S,data_I)
  error <- 1/106*sum((model-data)**2)
  return (error)
}

############ PARAMETER ESTIMATION WITH PSO #############
set.seed(2020)
pso_outcome <- psoptim(c(NA,NA,NA,NA), error_f,
                       lower = c(0.0000013, 0.000063, 0.000043, 0.000028),
                       upper = c(0.0000020, 0.000188, 0.000090, 0.000080),
                       control = list(maxit=20))

cat("d=", pso_outcome$par[1],"\ni =", pso_outcome$par[2], "\nr =", pso_outcome$par[3],"\nv =", pso_outcome$par[4],"\nerror =", pso_outcome$value)

#Output
#d= 1.585656e-06 
#i = 6.406811e-05 
#r = 4.530534e-05 
#v = 7.389872e-05 
#error = 0.01536676


######### SOLVING ODE WITH PARAMETERS FROM PSO #########

#values of t
time <- seq(0, 52, by = 1)
#initial values, when t=0
state <- c(S=569.6715, I=0.3285, M=0) #0.3285 because number of cases in 2018 is 3285

parameters <- c(d = pso_outcome$par[1],i = pso_outcome$par[2],r = pso_outcome$par[3],v = pso_outcome$par[4])
out <- ode(y = state, time = time, func = RHS, parms = parameters)

#extracting as a data frame
out.df <- as.data.frame(out)
out.df



#################### MODEL FITTING #####################

#Plotting of infected cases
plot(time, data_I,main = "Infected Cases over Time", xlab = "Number of weeks", 
     ylab ="Number of Infected cases (Ten Thousand)", col=2, pch=16,ylim=c(0,2.5))
par(new=TRUE)
plot(time, out.df[,3],main = "Infected Cases over Time",xlab = "Number of weeks", 
     ylab ="Number of Infected cases (Ten Thousand)",type='l', col=4, ylim=c(0,2.5))

legend("topleft", c("Model","Data"), lty=c(1,NA),pch=c(NA,16),col=c("blue","red"))


#Plotting of susceptible cases
plot(time, data_S,main = "Susceptible Cases over Time", xlab = "Number of weeks",
     ylab ="Number of Susceptible cases (Ten Thousand)", col=2, pch=16,ylim=c(565,570))
par(new=TRUE)
plot(time, out.df[,2],main = "Susceptible Cases over Time",xlab = "Number of weeks",
     ylab ="Number of Susceptible cases (Ten Thousand)", type='l', col=4, ylim=c(565,570))

legend("bottomleft", c("Model","Data"), lty=c(1,NA),pch=c(NA,16),col=c("blue","red"))