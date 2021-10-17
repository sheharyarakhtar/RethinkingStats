#Chapter 5 - Multivariate models
library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce
#standardize predictors
d$A <- standardize(d$MedianAgeMarriage)
d$M <- standardize(d$Marriage)
d$D <- standardize(d$Divorce)
head(d[,c(1,4,5,7,14,15,16)])

#Fit model

m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA *A,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    sigma ~ dunif(0,10)
  ), data = d
)

precis(m5.1)

MAM.seq <- seq(from =-3.5, to = 3.5, length.out = 30)
a <- extract.samples(m5.1)

mu <- link(m5.1, data = data.frame(A=MAM.seq))
mu.PI <- apply(mu, 2, PI)
plot(D ~ A, data = d, col = rangi2)
abline(m5.1)
shade(mu.PI, MAM.seq)


d$Marriage.s <- (d$Marriage - mean(d$Marriage))/sd(d$Marriage)
m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bR *M,
    a ~ dnorm(0,0.2),
    bR ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ),data = d
)

M.seq <- seq(from = -2,to = 3, length.out = 30)
mu <- link(m5.2, data= data.frame(M = M.seq))
mu.PI <- apply(mu, 2, PI)
plot(D ~ M, data = d, col = col.alpha(rangi2, 1))
abline(m5.2)
shade(mu.PI, M.seq)
precis(m5.2)


##What does it mean to assume (mu = a + bR*Marriage rate + bA*Age)
#? It means that the expected
# outcome for any State with marriage rate Ri and median age at marriage Ai
# is the sum of three
# independent terms. The first term is a constant, ??. Every State gets this. The second term
# is the product of the marriage rate, Ri
# , and the coefficient, ??R, that measures the association
# between marriage rate and divorce rate. The third term is similar, but for the association
# with median age at marriage instead.

#A State's divorce rate can be a function of its marriage rate or its
# median age at marriage

m5.3 <- map(
  alist(
    D~ dnorm(mu, sigma),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0,0.5),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data = d
)
precis(m5.3)

plot(precis(m5.3))

##Predictor residual plots

m5.4 <- quap(
  alist(
    Marriage.s ~ dnorm(mu, sigma),
    mu <- a + bA*MedianAgeMarriage.s,
    a ~ dnorm(0,10),
    bA ~ dnorm(0,1),
    sigma ~ dunif(0,10)
  ), data = d
)
#computing the residuals by subtracting the observed rate in each state 
#from the predicted rate

# Compute expected value at MAP, for each state
mu <- coef(m5.4)['a'] + coef(m5.4)['bA']*d$MedianAgeMarriage.s
# Compute residual fr each state
m.resid <- d$Marriage.s - mu
str(mu)

plot(Marriage.s ~ MedianAgeMarriage.s, d, col = rangi2)
abline(m5.4)
for( i in 1:length(m.resid)){
  x <- d$MedianAgeMarriage.s[i]
  y <- d$Marriage.s[i]
  lines(c(x,x),c(mu[i],y) , lwd = 0.5, cl = col.alpha('black', 0.8))
}
plot(Divorce~m.resid, d)
abline(v=mean(m.resid))
d$m.resid <- m.resid
m.resid
m5.4 <- map(
  alist(
    Marriage.s ~ dnorm(mu, sigma),
    mu <- a + b*m.resid,
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dunif(0,10)
  ), data = d
)
precis(m5.4)
resid <- seq(-1,1,length.out = 30)
mu <- link(m5.4, data = data.frame(m.resid = resid))
mu.PI <- apply(mu,2,PI)
shade(mu.PI, resid)
mu
##

library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

m5.1 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bA * A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

set.seed(10)
prior <- extract.prior( m5.1 )
mu <- link( m5.1 , post=prior , data=list( A=c(-2,2) ) )
plot( NULL , xlim=c(-2,2) , ylim=c(-2,2) )
for ( i in 1:50 ) lines( c(-2,2) , mu[i,] , col=col.alpha("black",0.4) )




A_seq <- seq( from=-3 , to=3.2 , length.out=30 )
mu <- link( m5.1 , data=list(A=A_seq) )
mu.mean <- apply( mu , 2, mean )
mu.PI <- apply( mu , 2 , PI )
# plot it all
plot( D ~ A , data=d , col=rangi2 )
lines( A_seq , mu.mean , lwd=2 )
shade( mu.PI , A_seq )


library(dagitty)
dag5.1 <- dagitty( "dag{ A -> D; A -> M; M -> D }")
coordinates(dag5.1) <- list ( x = c(A=0, D = z, M = 2), y = c(A = 0 , D = 1, M = 0))
drawdag(dag5.1)

DMA_dag2 <- dagitty('dag{ D <- A -> M }')
impliedConditionalIndependencies(DMA_dag2)
 DMA_dag1 <- dagitty('dag{D <- A -> M -> D}')
 impliedConditionalIndependencies(DMA_dag1)
 
##APPROXIMATIN THE POSTERIOR
m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA*A + bM*M,
    a ~ dnorm(0,0.2),
    bA ~ dnorm(0,0.5),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data = d
)
precis(m5.3) 
plot( coeftab(m5.1, m5.2, m5.3), par = c('bA', 'bM'))

# The posterior means are shown by the points and the 89% compatibility intervals by the
# solid horizontal lines. Notice how bA doesn't move, only grows a bit more uncertain, while
# 134 5. THE MANY VARIABLES & THE SPURIOUS WAFFLES
# bM is only associated with divorce when age at marriage is missing from the model. You can
# interpret these distributions as saying:
# Once we know median age at marriage for a State, there is little or no additional predictive power in also knowing the rate of marriage in that State.


#PREDICTOR RESIDUAL PLOTS
# A predictor residual is the average prediction error when
# we use all of the other predictor variables to model a predictor of interest.

m5.4 <- quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + bAM * A,
    a ~ dnorm(0,0.2),
    bAM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
precis(m5.4)


#now we compute the residuals by subtracting
#the observed marriage rate in each State from the predicted rate based upon the model

mu <- link(m5.4)
mu_mean <- apply(mu, 2, mean)
mu_resid <- d$M - mu_mean

plot(d$D ~ mu_resid, col = col.alpha(rangi2, 1))
abline(v = 0)
abline(b = cor(mu_resid, d$D), a= 0)
cor(mu_resid, d$D)



m5.5 <- quap(
  alist(
    A ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
mu <- link(m5.5)
mu_mean <- apply(mu, 2,mean)
mu_resid1 <- d$A - mu_mean
plot(d$D ~ mu_resid1, col = col.alpha(rangi2, 1))
abline(v=0)
abline(b = cor(mu_resid1, d$D), a= 0)
cor(mu_resid1, d$D)

sim.height <- sim(m4.3, data = list(weight = weight.seq), n=1e5)

# In this, basically, we modeled age on marriage rate and then marriage rate on
# age. After getting the model, we found the residuals for these. The
# residuals are basically telling us that having known one value, how
# much additional information are we getting by the reminents. We do this
# by plotting the residuals against divorce rates. A low correlation
# between these means that after conditioning on the other variable, we are
# not deriving any additional insight from the coditioned variable and vice versa
# for high correlation residual variables





# States with positive residuals have high marriage rates for median age
# States with negative residuals have low marriage rates for median age



##POSTERIOR PREDICTION PLOTS
mu <- link(m5.3)
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI)

D_sim <- sim(m5.3,n = 1e4)
D_PI <- apply(D_sim,2,PI)

plot(mu_mean ~ d$D, col = rangi2, ylim=range(mu_PI),
     xlab = 'Observed divorce', ylab = 'Predeicted divorce')              
abline(a=0,b=1,lty=2)
for(i in 1:nrow(d)){
  lines(rep(d$D[i],2), mu_PI[,i],col = rangi2)
} 
identify( x=d$D , y=mu_mean , labels=d$Loc )


##Counterfactual plots
# (1) Pick a variable to manipulate, the intervention variable.
# (2) Define the range of values to set the intervention variable to.
# (3) For each value of the intervention variable, and for each sample in posterior, use
# the causal model to simulate the values of other variables, including the outcome.

# To estimate the influence of A on M, all we need is to regress A on M. There are no
# other variables in the DAG creating an association between A and M. We can just add this
# regression to the quap model, running two regressions at the same time:

data("WaffleDivorce")
d <- list()
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$M <- standardize(WaffleDivorce$Marriage)
d$D <- standardize(WaffleDivorce$Divorce)

m5.3_A <- quap(
  alist(
    ##A -> D <- M
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0,0.5),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1),
    ## A -> M
    M ~ dnorm(mu_M, sigma_M),
    mu_M <- aM + bAM*A,
    aM ~ dnorm(0,0.2),
    bAM ~ dnorm(0,0.5),
    sigma_M ~ dexp(1)
  ), data = d
)
m5.3_A
plot(precis(m5.3_A))
precis(m5.3_A)

# The goal is to simulate what would happen, if we manipulate A. So next we define a
# range of values for A.
A_seq <- seq(from = -2, to = 2, length.out = 30)

# Now we can use sim to simulate observations from model m5.3_A
# But this time we'll tell it to simulate both M and
# D, in that order. Why in that order? Because we have to simulate the influence of A on M
# before we simulate the joint influence of A and M on D. The vars argument to sim tells it
# both which observables to simulate and in which order.
#prepData
sim_dat <- data.frame(A=A_seq)
#simulate M and then D, using A_seq
s <- sim(m5.3_A, data = sim_dat, vars = c('M','D'))
str(s)

plot(sim_dat$A, colMeans(s$D), ylim = c(-2,2), type = 'l',
     xlab = 'manipulate A', ylab = 'counterfactual D')
shade(apply(s$D,2,PI), sim_dat$A)
mtext('Total counterfactual effect of A on D')

plot(sim_dat$A, colMeans(s$M), ylim = c(-2,2), type = 'l',
     xlab = 'manipulate A', ylab = 'counterfactual M')
shade(apply(s$M,2,PI), sim_dat$A)
mtext('Total counterfactual effect of A on M')

# Of course these calculations also permit numerical summaries.
# For example, the expected causal effect of increasing median age at marriage from 20 to 30 is:

sim2_dat <- data.frame(A = (c(20,30)-26.1)/1.24)
s2 <- sim(m5.3_A, data = sim2_dat, vars = c('M', 'D'))
mean(s2$D[,2] - s2$D[,1])

#We make a model saying that D is influenced by both M and A, but we also say that A influences M. 
#After creating this model, we see the effect has on both M and D
#Next we try to break the causal relationship between M and D by manipulating A, by changing the
#median age at marriage from 20 to 30, and we simulate values of D on this new A
#we now numerically see that changing the age of marriage strongly influences marriage rate. lets see if the 
#reverse is also true

sim2_dat <- data.frame(M = seq(from = -2, to = 2, length.out =30), A=0)
s3 <- sim(m5.3_A, data = sim2_dat, vars = c('D'))
plot(sim2_dat$M, colMeans(s3), ylim = c(-2,2),
     type = 'l', xlab = 'Manipulated M', ylab = 'counterfactual D')
shade(apply(s3, 2,PI), sim2_dat$M)
mtext('Total counterfactual effect of M on D')

mean(s3[,2]-s3[,1])



##First method was to find spurious correlations, and see which variable
#is the counfound is which is just mimicing
#Now we try to find variables that mask one another
# to measure the direct influences of multiple factors on an outcome, when none of those
# influences is apparent from bivariate relationships. This kind of problem tends to arise when
# there are two predictor variables that are correlated with one another
data(milk)
d <- milk
str(d)

d$K <- standardize(d$kcal.per.g)
d$M <- standardize(d$mass)
d$N <- standardize(d$neocortex.perc)

m5.6 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N,
    a ~ dnorm(0,1),
    bN ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = d
)
##Error due to NA values, so we drop cases with missing values this is called
#Complete case analysis

dcc <- d[ complete.cases(d$K,d$N,d$M) , ]
m5.5_draft <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN*N,
    a ~ dnorm(0,0.2),
    bN ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data = dcc
)
#Are these priors reasonable? Lets plot the prior predictive dist to see
prior <- extract.prior(m5.5_draft)
xseq<- c(-2,2)
mu <- link(m5.5_draft, post = prior, data = list(N=xseq))
plot( NULL , xlim=xseq , ylim=xseq )
for ( i in 1:50 ) lines( xseq , mu[i,] , col=col.alpha("black",0.3) )

# now lets look at the posterior
plot(precis(m5.5_draft))
xseq <- seq(from=min(dcc$N-0.15), to = max(dcc$N)+0.15, length.out = 30)
mu <- link(m5.5_draft, data = list(N=xseq))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI)
plot(K~N, data = dcc, col = col.alpha(rangi2,1))
lines(xseq,mu_mean,lwd=2)
shade(mu_PI, xseq)

dcc$LM.S <- standardize(log(dcc$mass))

m5.6 <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bM*LM.S ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=dcc )

plot(precis(m5.6))

xseq <- seq(from = min(dcc$LM.S)-0.15, to = max(dcc$LM.S)+0.15, length.out = 30)
mu <- link(m5.6, data = list(LM.S=xseq))
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI)
plot(K ~ LM.S, data = dcc)
lines(xseq,mu_mean,lwd=2)
shade(mu_PI, xseq)

#Both M and N did not show a significantly high influence on K,
#However, lets see what the multivariate model thinks

m5.7 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bM*LM.S + bN*N,
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0,0.5),
    bN ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data =dcc
)
plot(precis(m5.7))

plot(coeftab(m5.7,m5.6,m5.5_draft), pars = c("bM", 'bN'))


xseq <- seq( from=min(dcc$N)-0.15 , to=max(dcc$N)+0.15 , length.out=30 )
mu <- link( m5.7 , data=data.frame( LM.S=0 , N=xseq ) )
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI)
plot( NULL , xlim=range(dcc$N) , ylim=range(dcc$K) )
lines( xseq , mu_mean , lwd=2 )
shade( mu_PI , xseq )


xseq <- seq( from=min(dcc$LM.S)-0.15 , to=max(dcc$LM.S)+0.15 , length.out=30 )
mu <- link( m5.7 , data=data.frame( LM.S=xseq , N=0 ) )
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI)
plot( NULL , xlim=range(dcc$LM.S) , ylim=range(dcc$K) )
lines( xseq , mu_mean , lwd=2 )
shade( mu_PI , xseq )


data("Howell1")
d <- Howell1
str(d)
d$sex <- ifelse(d$male==1,2,1)
str(d)
m5.8 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a[sex],
    a[sex] ~ dnorm(178,20),
    sigma ~ dunif(0,50)
  ), data = d
)
par(mfrow = c(1,1))
precis(m5.8, depth =2)
post <- extract.samples(m5.8)
post$diff_fm <- post$a[,1] - post$a[,2]
precis(post, depth = 2, hist = FALSE)
