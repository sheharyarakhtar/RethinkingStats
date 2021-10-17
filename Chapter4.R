#Chapter 4
 #prior predictive distribution
# h ~ N(mu, sigma)
# mu ~ N(178,20)
# sigma ~ U(0,50)
library(rethinking)
sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0,100)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
par(mfrow=c(1,1))
plot(prior_h)
dens(prior_h)


##FOOTBALL LINE EXP
#Normal by addition
pos <- replicate(1e4, sum(runif(10000,-1,1)))
plot(pos)
dens(pos, norm.comp = T)
plot(density(pos))
hist(pos)

#Normal by multiplication
# 12 genes in total, each gene contribute a % to total growth
# the product of all is the total growth of the organism
growth <- replicate(1e4, prod(1+runif(12,0,0.1)))
dens(growth, norm.comp = T)
#smaller values tend to be approximately equal toadditon 
# meaning 1.1x1.1=1.21 and 1.1+1.1=1.2, which are prettclose(
big <- replicate( 10000 , prod( 1 + runif(12,0,0.5) ) )
small <- replicate( 10000 , prod( 1 + runif(12,0,0.01) ) )
dens(big, norm.comp = T)
dens(small, norm.comp = T)
##We see that this rule does not work when we have a large value for 
# multiplicatin

# However, we also see that for large value products, a log normal distribution 
# can make compensations

dens(log(big), norm.comp = T)


##Exmple of making a dist of dists
data("Howell1")
d <- Howell1
View(d)
str(d)

d2 <- d[d$age>=18,]
str(d2)

dens(d2$height, norm.comp = T)
# 
# We want to model this using a gausing distribution
# but what mean and std do we use? There are an infinite number of 
# means and std that produce a gausian dist. First thing we do is that
# we define the model


# h_i ~ Normal(mew, sigma) where i = 1:352

# the model above is saying that all the golem knows about each
# height measurement is defined by the same normal distribution, with mean µ
# and standard deviation ??.

curve( dnorm( x , 178 , 20 ) , from=100 , to=250 ) ##initial prior for mew
curve( dunif( x , 0 , 50 ) , from=-10 , to=60 )  #initial prior for sigma
# 
# Model {
#   h_i ~ Normal(mew, sigma)
#   mew ~ Normal(178, 20)
#   sigma ~ Uniform(0,50)
# }

#We now sample from the prior for individual heights
sample_mu <- rnorm(1e4, mean = 178, sd=20)
sample_sigma <- runif(1e4, min=0, max=50)
prior_h <- rnorm(1e5, mean = sample_mu, sd= sample_sigma)
dens(prior_h, norm.comp = T)

##Grid approximation of the posterior
mu.list <- seq(from=150, to=160 , length.out=200 )

sigma.list <- seq( from=4 , to=9 , length.out=200 )

post <- expand.grid( mu=mu.list , sigma=sigma.list )  ##Long format combinations of all mu with all sigma values

# In the monstrous fourth line of code, shown in expanded form to make it easier
# to read, the log-likelihood at each combination of µ and ?? is computed. This line looks so awful, because we
# have to be careful here to do everything on the log scale. Otherwise rounding error will quickly make all of the
# posterior probabilities zero. So what sapply does is pass the unique combination of µ and ?? on each row of
# post to a function that computes the log-likelihood of each observed height, and adds all of these log-likelihoods
# together (sum).
post$LL <- sapply( 1:nrow(post) , function(i) {sum(dnorm(d2$height, mean=post$mu[i],sd=post$sigma[i],log=TRUE))})

post$prod <- post$LL + dnorm( post$mu , 178 , 20 , TRUE ) + dunif( post$sigma , 0 , 50 , TRUE )

post$prob <- exp( post$prod - max(post$prod) )

contour_xyz( post$mu , post$sigma , post$prob, xlim = c(152,156), ylim = c(7,8.5) )

image_xyz( post$mu , post$sigma , post$prob, xlim = c(152,156), ylim = c(7,8.5)  )

precis(post)
# we first randomly sample row numbers in post
# in proportion to the values in post$prob. 
# Then we pull out the parameter values on those
# randomly sampled rows
sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                       prob=post$prob )
sample.mu <- post$mu[ sample.rows ]
sample.sigma <- post$sigma[ sample.rows ]
plot( sample.mu , sample.sigma , cex=0.5 , pch=16 , col=col.alpha(rangi2,0.8))

dens(sample.sigma)
dens(sample.mu)
par(mfrow = c(2,2))
HPDI( sample.mu )
HPDI( sample.sigma )

# Now we leave grid approximation behind and move on
# to one of the great engines of applied statistics, the quadratic approximation

##DEFINE MODEL

flist <- alist(
  height~ dnorm(mu, sigma),
  mu ~ dnorm(178, 20),
  sigma ~ dunif(0,50)
)
##fit this model on the data in d2, where d2 <- d[d$height>=18]
m4.1 <- map(flist, data = d2)
precis(m4.1)


m4.2 <- quap(alist(
  height ~ dnorm(mu, sigma),
  mu ~ dnorm(178, 20),
  sigma ~ dunif(0,50)
), data = d2)
precis(m4.2)

vcov(m4.1)
diag(vcov(m4.1))
cov2cor(vcov(m4.1))

##sampling from the poterior quadractic approx
post <- extract.samples(m4.1, n=1e4)
head(post)
precis(post)
precis(m4.1)

plot(post,cex = 0.5, pch=16, col = col.alpha(rangi2,0.8))
par(mfrow = c(2,2))
dens(post$sigma,cex = 0.3, pch=12, col = col.alpha(rangi2,0.8))
dens(post$mu,cex = 0.5, pch=16, col = col.alpha(rangi2,0.8),)

View(post)
precis(post)


library(rethinking)
data("Howell1")
d <- Howell1
head(d)
d2 <- d[d$age >= 18,]
par(mfrow = c(1,1))
plot(d2$weight,d2$height,cex = .8, pch=1)

#------------------LINEAR MODELING STRATEGY--------------------#

#let x be the name for the column of weight measurement
# 
# model will be :
#   height ~ normal(mu_i, sigma), ---Likelihood
#   mu_i = alpha +beta*(x_i - xbar), ------linear model
#   alpha ~ normal(178,20), ----alpha prior
#   beta ~ normal(0,10),--------beta prior
#   sigma ~ uniform(0,50)-------sigma prior
set.seed(2971)
N <- 100
a <- rnorm(N,178,20)
b<- rnorm(N, 0,10)

plot(height ~ weight,data = d2, xlim=range(d2$weight), ylim = range(0,400),
     xlab = 'weight', ylab = 'height')
abline(h=0,lty=2)
abline(h=272,lty=1,lwd=0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)
for(i in 1:N){
  curve(a[i] + b[i]*(x-xbar),
        from = min(d2$weight),to=max(d2$weight),
        add = T,
        col = col.alpha('black',0.2))
}
##^This tells us that the relationship between weight and height can be
# absurdly negative or completely positive. The golem doesnt know so it considers
# all possibilities

##We know that weight increases with height atleast to a certain point
# lets try incorporating that into out model
b <-rlnorm(1e4,0,1)  #we are using a lognormal prior for beta as we can assume that relationship is positive for the most part
dens(b, xlim=c(0,5), adj = 0.1)  


set.seed(2971)
N <- 100
a <- rnorm(N,178,20)
b<- rlnorm(N, 0,1)


##we use NULL here in plot function because its better to not view the data before setting your priors to avoid bias
plot(height ~ weight, data = d2, xlim=range(d2$weight), ylim = range(0,400),
     xlab = 'weight', ylab = 'height', col = rangi2)
abline(h=0,lty=2)
abline(h=272,lty=1,lwd=0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)
for(i in 1:N){
  curve(a[i] + b[i]*(x-xbar),
        from = min(d2$weight),to=max(d2$weight),
        add = T,
        col = col.alpha('black',0.2))
}


###Now that we have the priors, we make the regression models with the posterior
m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight-xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,50)
  ), data = d2
)


precis(m4.3)
round(vcov(m4.3),3)
cov2cor(vcov(m4.3))
pairs(m4.3)

##Plotting the posterior distribution
plot(height ~ weight, data = d2)
post <- extract.samples(m4.3)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map + b_map*(x - xbar), add = T, col = rangi2)

#This isnt a bad line but there are other plausible lines near it as well
# so we plot those too to show the uncertainty of the regression model

head(post)

##Add the data slowly to see how the evidence changes the regression
par(mfrow = c(1,1))
N <- 352
dN <- d2[1:N,]
mN <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,50)
  ), data = dN
)

post <- extract.samples(mN,30)
plot(dN$weight, dN$height,
     xlim = range(d2$weight), ylim = range(d2$height),
     col = rangi2,
     xlab = 'weight', ylab = 'height')
mtext(concat("N = ",N))

for(i in 1:dim(post)[1]){
  curve(post$a[i] + post$b[i]*(x - mean(dN$weight)),
        col=col.alpha('black',0.3), add = T)
}

mu_at_50 <- post$a + post$b * 50

mu <- link(m4.3)
str(mu)
weight.seq <- seq( from =25, to = 70, by =1)
mu <- link(m4.3, data = data.frame(weight=weight.seq))
str(mu)
plot(height ~ weight , d2, type = 'n')
for(i in 1:100) points(weight.seq, mu[i,], pch=16, col = col.alpha(rangi2,0.1))
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob= 0.89)

plot(height ~ weight , data = d2, col = col.alpha(rangi2, 0.9))
lines(weight.seq, mu.mean)
shade(mu.HPDI, weight.seq, alpha(0.9))


##Plotting HDPI interval for heights

par(mfrow = c(1,1))
sim.height <- sim(m4.3, data = list(weight = weight.seq), n=1e5)
height.P1 <- apply(sim.height, 2, HPDI, prob = 0.89)
plot(height ~ weight, d2, col = col.alpha(rangi2, 0.8))
lines(weight.seq, mu.mean)
shade(mu.HPDI, weight.seq, alpha(0.8))
shade(height.P1, weight.seq, alpha(0.8))

# so basically, using the link function, we can simulate the MAP function output,
# like an inverse. You give it the weight values and it simulates the mu of heights
# on those weights. You get 1000 mu samples for every weight value you give it.
# Then we calculate the mu of those samples on column dimension to get the mean mu,
# and we plot the line for it. Then we get the HPDI for those simulated means to show
# credibility interval of those samples mus. 
# We use the sim function to simulate the heights( not the mus ) for the given
# values of weights, and then plot the PI or HPDI on those height values to show
# the uncertainty they entail


#POLYNOMIAL REGRESSION
library(rethinking)
data("Howell1")
d <- Howell1
plot(height ~ weight , data = d, col=col.alpha(rangi2, 0.9))

# if we consider a polynomial model the distribution of mu would be non-linear.
# in the model,we should get mu = a +b1*x + b2*x^2

# we should standardise the predictor variable in polynomial regression especially
d$weight_s <- (d$weight-mean(d$weight))/sd(d$weight)
d$weight_s2 <- (d$weight_s)^2
d$weight_s3 <- d$weight_s^3

m4.4 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2 + b3*weight_s3,
    a ~ dnorm(178,20),
    b1~ dlnorm(0,1),
    b2 ~ dnorm(0,1),
    b3 ~ dnorm(0,1),
    sigma ~ dunif(0,50)
  ), data = d
)
precis(m4.4)

weight.seq <- seq(from = -2.2, to = 2, length.out = 30)
pred_dat <- list(weight_s=weight.seq, weight_s2 = weight.seq^2,
                 weight_s3 = weight.seq^3)
mu <- link(m4.4, data = pred_dat)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob = 0.89)
sim.height <- sim(m4.4, data = pred_dat)
height.PI <- apply(sim.height,2,PI, prob=0.89)
plot(height ~ weight_s,d,col=col.alpha(rangi2,0.8),xaxt = 'n')
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)

plot(height ~ weight, d, xaxt = 'n')
at <- c(-2,-1,0,1,2)
labels <- at*sd(d$weight) + mean(d$weight)      
axis(side=1, at=at, labels = round(labels,0))
