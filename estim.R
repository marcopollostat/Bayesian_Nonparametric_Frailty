
rm(list = ls()) 
Rcpp::sourceCpp('NPB25.cpp') 
X          <- read.csv("cd carros.csv") 
X          <- as.matrix(X) 
T          <- 3000         
r          <- 267          
LF         <- 10           


parhiperc  = c(4,4)       
parnormal  = c(0,1,2,2)     


estlass  <- MVSM( X, T = 3000, r ) 
estlass[1] 
estlass[2] 
estlass[3] 
estlass[4]
estlass[5]


m     <- 439           
beta  <- estlass[2][[1]]  
alpha <- estlass[3][[1]]/m 


pbeta    = cbind(estlass[3][[1]], estlass[3][[1]]/beta )
palpha   = cbind(estlass[3][[1]], c(m,m,m) )

epsilon <- 0.05
set.seed(123)
aj1     <- GHMCRSMNB( X, T = 3000, r =267 , warmup = 5000, iter = 5000, epsilon, LF = 20, parhiperc , parnormal)

aj1[[ 1,2]]

Z  <- as.matrix(aj1[1,1][[1]]) 

mediazs <- apply( Z, 1, mean )


plot(1:172,mediazs[1:172], main="Estimated frailty versus car", xlab="Cars", ylab=expression(hat(z)[j]) )
abline(h=1, col="red")


plot(X[1:172,2],mediazs[1:172], main="Estimated frailty versus observed failure time for each car", xlab="mileage", ylab=expression(hat(z)[j]) )
abline(h=1, col="red")

###############
# SD and IC
###############
which(mediazs>2); sum( mediazs[1:172] > 2 ) 

variazs <- apply( Z, 1, sd )

round( mediazs[1], 3)
round( variazs[1], 3 )
round( quantile( Z[1,], prob = c(0.025,0.975) ), 3 )

round( mediazs[17], 3)
round( variazs[17], 3 )
round( quantile( Z[17,], prob = c(0.025,0.975) ), 3 )

round( mediazs[26], 3)
round( variazs[26], 3 )
round( quantile( Z[26,], prob = c(0.025,0.975) ), 3 )

round( mediazs[161], 3)
round( variazs[161], 3 )
round( quantile( Z[161,], prob = c(0.025,0.975) ), 3 )

round( mediazs[165], 3)
round( variazs[165], 3 )
round( quantile( Z[165,], prob = c(0.025,0.975) ), 3 )

round( mediazs[169], 3)
round( variazs[169], 3 )
round( quantile( Z[169,], prob = c(0.025,0.975) ), 3 )


varZ <- apply( Z, 2, var )
mean(varZ)
ICz <- round( quantile( varZ, prob = c(0.025,0.975) ), 3 )
ICz


hist(varzs, breaks = 100, main = NULL, col = "gray")
p <- c(varzs[5],varzs[17],varzs[26],varzs[57],varzs[98],varzs[126],varzs[129],varzs[161],varzs[165],varzs[169])
points(p, y=rep(0,10), pch=4, col = 2)

boxplot(varzs, horizontal=TRUE,  outline=TRUE,  ylim=c(0,20), frame=F, col = "green1", add = TRUE)


Cpar <- as.numeric(aj1[[1,3]]) 
str(Cpar) 
summary(Cpar) 
round(mean(Cpar),3) 
round(sd(Cpar),3) 
round( quantile( Cpar, prob = c(0.025,0.975) ), 3 )

summary(Cpar/(1+Cpar)) 

acf(Cpar)
ts.plot(Cpar)
ts.plot(Cpar/(1+Cpar))

hist(Cpar, freq = F, breaks = 100, xlab="C",ylab="Density", main = NULL, xlim = c(0,8))
abline(v = mean(Cpar), col=2, lty=2, lwd=2)

###########
# plots
###########

#pdf('grafCpar.pdf')

par( mfrow = c(2, 3) )

hist(Cpar, freq = F, breaks = 50, xlab="c", ylab="Density", main = "Histogram of c", xlim = c(0,8))
abline(v = mean(Cpar), col=2, lty=2, lwd=2)
acf(Cpar, main = "Autocorrelations for c", xlab = "Lag")
ts.plot(Cpar, main = "Trace plot of c", xlab = "Iterations", ylab = "c")


plot(1:172, mediazs[1:172], main = "Estimated individual frailties", xlab = "cars", ylab = expression(hat(z_j))) 
abline(h=1, col="red", lty=2, lwd=2)
plot( 1:172, varzs[1:172], main = "Estimated variance of the individual frailties", xlab = "cars", ylab = "Var(z_j)" )
abline(h=1, col="blue")
abline(h=5, col="red", lty=2, lwd=2 )
hist(varzs, breaks = 100, main = "Histogram of estimated variance", xlab = "variance values", col = "gray")
p <- c(varzs[5],varzs[17],varzs[26],varzs[57],varzs[98],varzs[126],varzs[129],varzs[161],varzs[165],varzs[169])
points(p, y=rep(0,10), pch=4, col = 2)

#dev.off()



DATA<-as.vector(unlist(estlass[4]))
estmfraj <-mediazs[1:172]
numfailj <-DATA
df <- data.frame(numfailj, estmfraj)

library(ggplot2)
library(ggExtra)


p=ggplot(df, aes(x=numfailj, y=estmfraj )) +
        geom_point(aes(colour = "blue"), )  +
        theme(legend.position="none") +
        xlab("Cumulative number of failures by car") +
        ylab("Estimated value of Z_j") 

ggMarginal(p, type="histogram")


##############################
##############################

beta_1hat = ((n1-1)/(n1))*beta1
beta_2hat = ((n2-1)/(n2))*beta2
beta_3hat = ((n3-1)/(n3))*beta3

desvp_b1 = sqrt( ( (n1-1)/(n1^2) )*(beta1^2) )
desvp_b2 = sqrt( ( (n2-1)/(n2^2) )*(beta2^2) )
desvp_b3 = sqrt( ( (n3-1)/(n3^2) )*(beta3^2) )

IC_b1 = qgamma(c(0.025,0.975),n1-1,n1/beta1)
IC_b2 = qgamma(c(0.025,0.975),n2-1,n2/beta2)
IC_b3 = qgamma(c(0.025,0.975),n3-1,n3/beta3)


alfa_1hat = n1/m
alfa_2hat = n2/m
alfa_3hat = n3/m

desvp_a1 = sqrt((n1)/(m^2))
desvp_a2 = sqrt((n2)/(m^2))
desvp_a3 = sqrt((n3)/(m^2))

IC_a1 = qgamma(c(0.025,0.975),n1,m)
IC_a2 = qgamma(c(0.025,0.975),n2,m)
IC_a3 = qgamma(c(0.025,0.975),n3,m)

round(c(beta_1hat,beta_2hat,beta_3hat),3)

round(c(IC_b1,IC_b2,IC_b3),3)

round(c(alfa_1hat,alfa_2hat,alfa_3hat),3)

round(c(IC_a1,IC_a2,IC_a3),3)

round(c(desvp_b1,desvp_b2,desvp_b3),3)

round(c(desvp_a1,desvp_a2,desvp_a3),3)


##########################
#CODA diagnostic
##########################
require(coda)
Zts  <- as.mcmc(t(Z))      
effZ <- effectiveSize(Zts) 
gewZ <- geweke.diag(Zts)   
   

plot(gewZ$z, main = "Geweke Diagnostic",ylab="Z-Statistic",xlab="Cars")
abline( 1.96,0,lty=2)
abline(-1.96,0,lty=2)