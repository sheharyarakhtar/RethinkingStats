library(rethinking)
data(chimpanzees)
d <- chimpanzees

d$treatment <- 1 + d$prosoc_left + 2*d$condition
xtabs( ~ treatment + prosoc_left + condition , d )

m11.1 <- quap(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a ,
    a ~ dnorm( 0 , 1.65)
  ) , data=d )

set.seed(1999)
prior <- extract.prior( m11.1 , n=1e6 )

p <- inv_logit( prior$a )
dens( p , adj=0.1 , ylim = c(0,10), xlab = 'CR')
dens( prior$a , adj=0.1 , ylim = c(0,1), xlab = 'Logit CR')


m11.2 <- quap(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a + b[treatment] ,
    a ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=d )
set.seed(1999)
prior <- extract.prior( m11.2 , n=1e4 )
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )
dens( abs( p[,1] - p[,2] ) , adj=0.1 )



m11.3 <- quap(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a + b[treatment] ,
    a ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=d )
set.seed(1999)
prior <- extract.prior( m11.3 , n=1e4 )
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )
mean( abs( p[,1] - p[,2] ) )

plot(precis(m11.1))

# trimmed data list
dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment) )

m11.4 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list , chains=4 , log_lik=TRUE)
precis( m11.4 , depth=2 )


post <- extract.samples(m11.4)
p_left <- inv_logit( post$a )
plot( precis( as.data.frame(p_left) ) , xlim=c(0,1) )

labs <- c("R/N", "L/N", "R/P", "L/P")
plot(precis(m11.4, depth = 2, pars = 'b'), labels = labs)

#we are concerned with the difference between R/N and R/P
# and same for left, to see how much prosocial option influenced

diffs <- list(
  db13 = post$b[,1] - post$b[,3],
  db24 = post$b[,2] - post$b[,4]
)
plot(precis(diffs))

##consider a posterior prediction check
pl <- by(d$pulled_left, list(d$actor, d$treatment), mean)
pl[1,]

##Posterior prediction check

plot( NULL , xlim=c(1,28) , ylim=c(0,1) , xlab="" ,
      ylab="proportion left lever" , xaxt="n" , yaxt="n" )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
abline( h=0.5 , lty=2 )
for ( j in 1:7 ) abline( v=(j-1)*4+4.5 , lwd=0.5 )
for ( j in 1:7 ) text( (j-1)*4+2.5 , 1.05 , concat("actor ",j) , xpd=TRUE )
for ( j in (1:7)[-2] ) {
  lines( (j-1)*4+c(1,3) , pl[j,c(1,3)] , lwd=2 , col=rangi2 )
  lines( (j-1)*4+c(2,4) , pl[j,c(2,4)] , lwd=2 , col=rangi2 )
}
points( 1:28 , t(pl) , pch=16 , col="white" , cex=1.7 )
points( 1:28 , t(pl) , pch=c(1,1,16,16) , col=rangi2 , lwd=2 )
yoff <- 0.01
text( 1 , pl[1,1]-yoff , "R/N" , pos=1 , cex=0.8 )
text( 2 , pl[1,2]+yoff , "L/N" , pos=3 , cex=0.8 )
text( 3 , pl[1,3]-yoff , "R/P" , pos=1 , cex=0.8 )
text( 4 , pl[1,4]+yoff , "L/P" , pos=3 , cex=0.8 )
mtext( "observed proportions\n" )

dat <- list(actor = rep(1:7,each=4),
            treatment=rep(1:4,times = 7))
p_post <- link(m11.4, data = dat)
p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post,2,PI)


d$side <- d$prosoc_left+1 #right 1 left 2
d$cond <- d$condition+1 #no partner 1 partner 2
dat_list2 <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  side = d$side,
  cond = d$cond)

m11.5 <- ulam(
  alist(
    pulled_left ~ dbinom(1,p),
    logit(p) <- a[actor]+bs[side]+bc[cond],
    a[actor] ~ dnorm(0,1.5),
    bs[side] ~ dnorm(0, 0.5), 
    bc[cond] ~ dnorm(0, 0.5)
  ), data = dat_list2, chains = 4, log_lik = TRUE
)

compare(m11.5, m11.4, func = WAIC)



###RELATIVE SHARK AND ABSOLUTE DEER
# The odds of switching from treatment 2 to treatment 4(adding a partner)
# is given by exponentiating the difference between the posterior predictions
post <- extract.samples(m11.4)
mean(exp(post$b[,4]-post$b[,2]))
##by adding a partner, your outcome changes from 1 of itself to 0.92 of itself
#this is a relative reduction of 8% in the odds

data("chimpanzees")
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition
d$cond <- d$condition + 1 #no paratner 1 partner 2
d$side <- d$prosoc_left +1 #right 1 left 2

d_aggregated <- aggregate(
  d$pulled_left ,
  list( treatment=d$treatment , actor=d$actor ,
        side=d$side , cond=d$cond ) ,
  sum )
colnames(d_aggregated)[5] <- "left_pulls"

dat <- with( d_aggregated , list(
  left_pulls = left_pulls,
  treatment = treatment,
  actor = actor,
  side = side,
  cond = cond ) )

m11.6 <- ulam(
  alist(
    left_pulls ~ dbinom( 18 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ) ,
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat , chains=4 , log_lik=TRUE )

post <- extract.samples(m11.6)
p_left <- post$a
plot(precis(as.data.frame(p_left)))
compare(m11.6, m11.4, func = PSIS)

library(rethinking)
##Agregared binomial on UCBarkely dataset
data("UCBadmit")
d <- UCBadmit
d
##We want to know whether there is evidence of gender bias in the admissions data
##Applicant gender is the predictor variable in our model
##We want to model "admit" variable as a fucntion of "gender"

# model :
#   A ~ Binomial(N, p),  N indicates the number of applications on row i
#   logit(p) = a[gender],
#   a[gender] ~ dnorm(0,1.5)

dat_list <- list(
  admit = d$admit,
  applications = d$applications,
  gid = ifelse(d$applicant.gender=='male',1,2)
)
m11.7 <- ulam(
  alist(
    admit ~ dbinom(applications, p),
    logit(p) <- a[gid],
    a[gid] ~ dnorm(0,1.5)
  ), data = dat_list, chains = 4
)

precis(m11.7, depth = 2)
plot(precis(m11.7, depth = 2))
post <- extract.samples(m11.7)

  # p_admit <- inv_logit(post$a)
  # plot(precis(as.data.frame(p_admit)))
  # exp(mean(p_admit[,1])-mean(p_admit[,2]))

diff_a <- post$a[,1] - post$a[,2] ##Absolute scale
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2]) #outcome scale
plot(precis(list(diff_a = diff_a, diff_p = diff_p)))
precis(list(diff_a = diff_a, diff_p = diff_p), hist = F)

##This has told us that there is about a 12-16% higher chance for a male to
#get admission than a female, lets do a posterior prediciton for model

postcheck(m11.7)
for(i in 1:6){
  x <- 1+2*(i-1)
  y1 <- d$admit[x]/d$applications[x]
  y2 <- d$admit[x+1]/d$applications[x+1]
  lines(c(x,x+1),c(y1,y2), col = rangi2, lwd= 2)
  text(x+0.5,(y1+y2)/2 +0.05, d$dept[x], cex=0.8, col = rangi2)
}

##We see that once we look at departmentwise admission prob,
##The overall admissions drop departmentwise, and there are only
#2 cases where men have a higher prob of admission than women

dat_list$dept_id <- rep(1:6, each = 2)

m11.8 <- ulam(
  alist(
    admit ~ binomial(applications, p),
    logit(p) <- a[gid]+delta[dept_id],
    a[gid] ~ dnorm(0,1.5),
    delta[dept_id] ~ dnorm(0,1.5)
  ), data = dat_list, chains = 4, iter = 4000, cores = 4
)
precis(m11.8, depth =2 )

postcheck(m11.8)
for(i in 1:6){
  x <- 1+2*(i-1)
  y1 <- d$admit[x]/d$applications[x]
  y2 <- d$admit[x+1]/d$applications[x+1]
  lines(c(x,x+1),c(y1,y2), col = rangi2, lwd= 2)
  text(x+0.5,(y1+y2)/2 +0.05, d$dept[x], cex=0.8, col = rangi2)
}

##calculating the contrast on both relative and absolute scale now
post <- extract.samples(m11.8)
##The intercept a of the two is pretty close now,
# lets calculate the contrast both on abs and rescale
diff_a <- post$a[,1] - post$a[,2]  ##absolute scale
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])  #outcome scale
precis(list(diff_a = diff_a, diff_p = diff_p), hist = F)

pg <- with(dat_list, sapply(1:6, function(k)
  applications[dept_id==k]/sum(applications[dept_id==k])))
rownames(pg) <- c('male', 'female')
colnames(pg) <- unique(d$dept)
round(pg,2)  ##this gives us the percentage of each gender applying to each dept
precis(m11.8, depth = 2)

pairs(m11.8)


##Poisson regression
#monk example
y <- rbinom(1e5, 1000, 1/1000)
c(mean(y), var(y))

##Oceanic tool complexity
library(rethinking)
data("Kline")
d <- Kline
d
d$P <- scale( log(d$population) )
d$contact_id <- ifelse( d$contact=="high" , 2 , 1 )


##prior predictive simulation for alpha
curve( dlnorm(x,3,0.5), from = 0, to = 100, n = 200, add = T)
mean(exp(rnorm(1e5,3,0.5)))

##sim for beta
# assuming initial guess is 0,1normal
N <- 100
a <- rnorm(N, 3, 0.5)
b <- rnorm(N, 0, 0.2)
plot(NULL, xlim = c(-2,2), ylim = c(0,100))
for(i in 1:N) curve(exp(a[i]+b[i]*x), add = T, col = rangi2)


x_seq <- seq( from=log(100) , to=log(200000) , length.out=100 )
lambda <- sapply( x_seq , function(x) exp( a + b*x ) )
plot( NULL , xlim=range(x_seq) , ylim=c(0,500) , xlab="log population" ,
      ylab="total tools" )
for ( i in 1:N ) lines( x_seq , lambda[i,] , col=grau() , lwd=1.5 )

plot( NULL , xlim=range(exp(x_seq)) , ylim=c(0,500) , xlab="population" ,
      ylab="total tools" )
for ( i in 1:N ) lines( exp(x_seq) , lambda[i,] , col=grau() , lwd=1.5 )



##Modeling
dat <- list(
  T = d$total_tools ,
  P = d$P ,
  cid = d$contact_id )
# intercept only
m11.9 <- ulam(
  alist(
    T ~ dpois( lambda ),
    log(lambda) <- a,
    a ~ dnorm( 3 , 0.5 )
  ), data=dat , chains=4 , log_lik=TRUE )
# interaction model
m11.10 <- ulam(
  alist(
    T ~ dpois( lambda ),
    log(lambda) <- a[cid] + b[cid]*P,
    a[cid] ~ dnorm( 3 , 0.5 ),
    b[cid] ~ dnorm( 0 , 0.2 )
  ), data=dat , chains=4 , log_lik=TRUE )

compare(m11.9, m11.10, func = PSIS)

##Here's the code to plot the data and superimpose posterior
# predictions for the expected number of tools at each population size and contact rate
k <- PSIS( m11.10 , pointwise=TRUE )$k
plot( dat$P , dat$T , xlab="log population (std)" , ylab="total tools" ,
      col=rangi2 , pch=ifelse( dat$cid==1 , 1 , 16 ) , lwd=2 ,
      ylim=c(0,75) , cex=1+normalize(k) )
# set up the horizontal axis values to compute predictions at
ns <- 100
P_seq <- seq( from=-1.4 , to=3 , length.out=ns )
# predictions for cid=1 (low contact)
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( P_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , P_seq , xpd=TRUE )
# predictions for cid=2 (high contact)
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( P_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , P_seq , xpd=TRUE )


plot( d$population , d$total_tools , xlab="population" , ylab="total tools" ,
      col=rangi2 , pch=ifelse( dat$cid==1 , 1 , 16 ) , lwd=2 ,
      ylim=c(0,75) , cex=1+normalize(k) )
ns <- 100
P_seq <- seq( from=-5 , to=3 , length.out=ns )
# 1.53 is sd of log(population)
# 9 is mean of log(population)
pop_seq <- exp( P_seq*1.53 + 9 )
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE )
lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE )



#scientific model m11.11 as compared to the geocentric model m11.10
dat2 <- list( T=d$total_tools, P=d$population, cid=d$contact_id )
m11.11 <- ulam(
  alist(
    T ~ dpois( lambda ),
    lambda <- exp(a[cid])*P^b[cid]/g,
    a[cid] ~ dnorm(1,1),
    b[cid] ~ dexp(1),
    g ~ dexp(1)
  ), data=dat2 , chains=4 , log_lik=TRUE )

##Exposure analysis
num_days <- 30
y <- rpois(num_days, 1.5)

num_weeks <- 4
y_new <- rpois(num_weeks, 0.5*7)

y_all <- c( y , y_new )
exposure <- c( rep(1,30) , rep(7,4) )
monastery <- c( rep(0,30) , rep(1,4) )
d <- data.frame( y=y_all , days=exposure , monastery=monastery )

# compute the offset
d$log_days <- log( d$days )
# fit the model
m11.12 <- quap(
  alist(
    y ~ dpois( lambda ),
    log(lambda) <- log_days + a + b*monastery,
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm( 0 , 1 )
  ), data=d )


post <- extract.samples( m11.12 )
lambda_old <- exp( post$a )
lambda_new <- exp( post$a + post$b )
precis( data.frame( lambda_old , lambda_new ) , hist = F)


##Multinomial and categorical models


#simulating career choices
# simulate career choices among 500 individuals
N <- 500 # number of individuals
income <- c(1,2,5) # expected income of each career
score <- 0.5*income # scores for each career, based on income
# next line converts scores to probabilities
p <- softmax(score[1],score[2],score[3])
# now simulate choice
# outcome career holds event type values, not counts
career <- rep(NA,N) # empty vector of choices for each individual
# sample chosen career for each individual
set.seed(34302)
for ( i in 1:N ) career[i] <- sample( 1:3 , size=1 , prob=p )

code_m11.13 <- "
data{
int N; // number of individuals
int K; // number of possible careers
int career[N]; // outcome
vector[K] career_income;
}
parameters{
vector[K-1] a; // intercepts
real<lower=0> b; // association of income with choice
}
model{
vector[K] p;
vector[K] s;
a ~ normal( 0 , 1 );
b ~ normal( 0 , 0.5 );
s[1] = a[1] + b*career_income[1];
s[2] = a[2] + b*career_income[2];
s[3] = 0; // pivot
p = softmax( s );
career ~ categorical( p );
}
"
dat_list <- list( N=N , K=3 , career=career , career_income=income )
m11.13 <- stan( model_code=code_m11.13 , data=dat_list , chains=4 )
precis( m11.13 , 2 )

post <- extract.samples(m11.13)
#set up logit scores

s1 <- with(post, a[,1]+b*income[1])
s2_orig <- with(post, a[,2]+b*income[2])
s2_new <- with(post, a[,2]+b*income[2]*2)

##compute probabilities for original and counter factual

p_orig <- sapply(1:length(post$b), function(i)
  softmax(c(s1[i], s2_orig[i],0)))
p_new <- sapply(1:length(post$b), function(i)
  softmax(c(s1[i], s2_new[i],0)))

p_diff <- p_new[2,] - p_orig[2,]
precis( p_diff , hist = F)


##modeling immpact of family income on career choice
N <- 500
##simulate family incomes for each individual
family_income <- runif(N)
##assign a unique coefficient for each type of event
b <- c(-2,0,2)
career <- rep(NA,N)

for(i in 1:N){
  score <- 0.5*(1:3)+ b*family_income[i]
  p <- softmax(score[1], score[2], score[3])
  career[i] <- sample(1:3, size = 1, prob = p)
}

code_m11.14 <- "
data{
int N; // number of observations
int K; // number of outcome values
int career[N]; // outcome
real family_income[N];
}
parameters{
vector[K-1] a; // intercepts
vector[K-1] b; // coefficients on family income
}
model{
vector[K] p;
vector[K] s;
a ~ normal(0,1.5);
b ~ normal(0,1);
for ( i in 1:N ) {
for ( j in 1:(K-1) ) s[j] = a[j] + b[j]*family_income[i];
s[K] = 0; // the pivot
p = softmax( s );
career[i] ~ categorical( p );
}
}
"
dat_list <- list(N=N, K=3, career=career, family_income=family_income)
m11.14 <- stan(model_code = code_m11.14, data = dat_list, chains = 4)
precis(m11.14, 2)
