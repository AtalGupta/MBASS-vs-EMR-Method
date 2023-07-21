# File:   ImportingData.R
# Course: R: An Introduction (with RStudio)

# INSTALL AND LOAD PACKAGES ################################

library(datasets)  # Load base packages manually

library(nlstools) #nlsfit() for emr() function
library (bootstrap)
pacman::p_load(pacman, bootstrap)

# Installs pacman ("package manager") if needed
if (!require("pacman")) install.packages("pacman")

# Use pacman to load add-on packages as desired
pacman::p_load(pacman, rio) 
pacman::p_load(pacman, nlstools)

# ABOUT EXCEL FILES ########################################

# From the official R documentation
#browseURL("http://j.mp/2aFZUrJ")

# You have been warned: ಠ_ಠ

# IMPORTING WITH RIO #######################################

# CSV
rio_csv <- import("/Users/aatal/Downloads/Copy of maths project dataset (1).xlsx")
head(rio_csv)



# DATA VIEWER ##############################################

#?View
#View(rio_csv)

#codes for the project 
#input parameter
mbin <- 0.1 #Magnitude bin

nbsample <- 200 #Bootstrapping
mag <- rio_csv$magnitude

#functions
fmd <- function(mag,mbin){
  mi <- seq(min(round(mag/mbin)*mbin), max(round(mag/mbin)*mbin), mbin)
  nbm <- length(mi)
  cumnbmag <- numeric(nbm)
  nbmag <- numeric(nbm)
  for(i in 1:nbm) cumnbmag[i] <- length(which(mag > mi[i]-mbin/2))
  cumnbmagtmp <- c(cumnbmag,0)
  nbmag <- abs(diff(cumnbmagtmp))
  res <- list(m=mi, cum=cumnbmag, noncum=nbmag)
  return(res)
}

#Maximum Curvature (MAXC) [e.g., Wiemer & Wyss, 2000]
maxc <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  Mc <- FMD$m[which(FMD$noncum == max(FMD$noncum))[1]]
  return(list(Mc=Mc))
}

#EMR method

emr <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  nbm <- length(FMD$m)
  McMAXC <- maxc(mag,mbin)$Mc
  mu <- abs(McMAXC/2); sig <- abs(McMAXC/4)
  if(mu > 1)mu <- abs(McMAXC/10); sig <- abs(McMAXC/20)
  McBound <- McMAXC
  Mco <- McBound-0.3+(seq(9)-1)/10
  params <- numeric(9*4); dim(params) <- c(9,4) #a, b, mu, sigma
  prob <- numeric(9)
  savedmodel <- numeric(9*nbm); dim(savedmodel) <- c(9,nbm)
  for(i in 1:9){
    indmag <- which(mag > Mco[i]-mbin/2)
    nbev <- length(indmag)
    b <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    a <- log10(length(indmag))+b*Mco[i]
    cumN <- 10^(a-b*FMD$m)
    params[i,1] <- a; params[i,2] <- b
    cumNtmp <- 10^(a-b*(max(FMD$m)+mbin))
    cumNtmp <- c(cumN, cumNtmp)
    N <- abs(diff(cumNtmp))
    data <- data.frame(N=N, m=FMD$m, Nd=FMD$noncum)
    indLow <- which(FMD$m < Mco[i]); indHigh <- which(FMD$m >= Mco[i])
    dataTest <- data.frame(N=data$N[indLow], m=data$m[indLow], Nd=data$Nd[indLow])
    dataTmp <- data.frame(N=data$N[indHigh], m=data$m[indHigh], Nd=data$Nd[indHigh])
    checkNo0 <- which(dataTest$Nd != 0)
    dataTest <- data.frame(N=dataTest$N[checkNo0], m=dataTest$m[checkNo0],
                           Nd=dataTest$Nd[checkNo0])
    #Nmax <- max(dataTmp$Nd)
    Nmax <- max(dataTest$Nd)
    #Nmax <- dataTest$Nd[length(dataTest$Nd)]
    Mmintmp <- min(dataTest$m)
    dataTest$Nd <- dataTest$Nd/Nmax
    dataTest$m <- dataTest$m-Mmintmp
    data4fit <- data.frame(N=dataTest$Nd, m=dataTest$m)
    #non-linear least squares fit
    nlsfit <- nls(N~pnorm(m, mean=mean, sd=sd), data=data4fit,
                  start=list(mean=mu, sd=sig), control=list(maxiter=100, warnOnly = TRUE))
    params[i,3] <- coef(nlsfit)["mean"]; params[i,4] <- coef(nlsfit)["sd"]
    dataTest$N <- pnorm(dataTest$m, mean=coef(nlsfit)["mean"],
                        sd=coef(nlsfit)["sd"])*Nmax
    dataTest$m <- dataTest$m+Mmintmp
    dataTest$Nd <- dataTest$Nd*Nmax
    dataPred <- data.frame(N=c(dataTest$N, dataTmp$N), m=c(dataTest$m, dataTmp$m),
                           Nd=c(dataTest$Nd, dataTmp$Nd))
    dataPred$N <- round(dataPred$N)
    savedmodel[i,c(checkNo0,indHigh)] <- dataPred$N
    probtmp <- numeric(nbm)
    CheckNo0 <- which(dataPred$N != 0)
    Pmodel <- dataPred$N[CheckNo0]; Pdata <- dataPred$Nd[CheckNo0]
    probtmp[CheckNo0] <- 1/log(10)*(-Pmodel+Pdata*log(Pmodel)-lgamma(Pdata+1))
    prob[i] <- -sum(probtmp)
  }
  indbestfit <- which(prob == min(prob, na.rm=TRUE))
  res <- list(Mc=Mco[indbestfit], a=params[indbestfit,1], b=params[indbestfit,2],
              mu=params[indbestfit,3], sigma=params[indbestfit,4],
              model=savedmodel[indbestfit,], Mco=Mco, prob=prob)
  return(res)
}

## COMPUTE Mc ##
Mc_bootstrap <- numeric(nbsample)
#select function: maxc(), gft(), mbs(), emr()
#For mbass(), see algorithm Amorese [2007]
#for(i in 1:nbsample) Mc_bootstrap[i] <- maxc(sample(mag, replace=TRUE),mbin)$Mc
#when using emr(), the loop may break due to failure of nlsfit(),
#in this case use:

for(i in 1:nbsample) Mc_bootstrap[i] <-
  as.numeric(try(emr(sample(mag, replace=TRUE),mbin)$Mc))
Mc_mean <- mean(Mc_bootstrap, na.rm=TRUE)
Mc_sd <- sd(Mc_bootstrap, na.rm=TRUE)
print(paste("Mc (mean): ", Mc_mean, sep=""))
print(paste("Sigma0 (std. dev.): ", Mc_sd, sep=""))


## PLOT FMD ##
FMD <- fmd(mag,mbin)
#pdf(paste(wd,"CORSSA_Mc_plot_FMD.pdf", sep=""))
plot(FMD$m, FMD$cum, log="y", xlab="Magnitude", ylab="Number of events",
     main="Frequency-Magnitude Distribution")
points(FMD$m, FMD$noncum, pch=2)
abline(v=Mc_mean)
legend("topright", c("Cum. FMD", "Non Cum. FMD"), cex=0.8, pch=c(1, 2))
dev.off()



?plot








# CLEAN UP #################################################

# Clear environment
rm(list = ls()) 

# Clear packages

p_unload(all)  # Remove all add-ons

# Clear console


cat("\014")  # ctrl+L

# Clear mind :)

