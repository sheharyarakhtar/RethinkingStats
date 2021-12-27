library(rethinking)
# Varying intercepts have variation, and varying slopes have variation. 
# Intercepts and slopes covary

# Simulate the population
a <- 3.5 #average morning wait time
b <- (-1) #average difference afternoon wait time
sigma_a <- 1 #std dev in intercepts
sigma_b <- 0.5 #std dev in slopes
rho <- (-0.7) #correlation between intercepts and slopes

Mu <- c(a,b)

cov_ab <- sigma_a*sigma_b*rho
Sigma <- matrix(c(sigma_a^2, cov_ab, cov_ab,sigma_b^2), ncol=2)

##Important approach to building a matrix in R
sigmas <- c(sigma_a, sigma_b) #standard deviations
Rho <- matrix(c(1,rho,rho,1), nrow=2)

#now matrix multiply to get covariance matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

#SIMULATE
N_cafes <- 20
library(MASS)
set.seed(5)
vary_effects <- mvrnorm(N_cafes, Mu, Sigma)
a_cafe <- vary_effects[,1]
b_cafe <- vary_effects[,2]
plot(a_cafe, b_cafe, col = rangi2,
     xlab = 'intercepts (a_cafe)', ylab = 'slopes (b_cafe)', pch = 20, cex= 1.5)

#Overlay population distribution 
library(ellipse)
for( l in c(0.1,0.3,0.5,0.88,0.99)){
  lines(ellipse(Sigma,centre = Mu, level = l), col=col.alpha('black',0.8))
}

# the covariance is the difference between the average product and the
# product of the averages

# What we did was simulate cafes and their average properties
# Now we simulate the observations that the robot takes
set.seed(20)
N_visits <- 10
afternoon <- rep(0:1,N_visits*N_cafes/2)
cafe_id <- rep(1:N_cafes, each = N_visits)
mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*afternoon
sigma <- 0.5 #std dev within cafes
wait <- rnorm(N_visits*N_cafes,mu,sigma)
d <- data.frame(cafe=cafe_id, afternoon=afternoon, wait=wait)


m14.1 <- ulam( 
  alist(
    wait ~ normal(mu, sigma),
    mu <- a_cafe[cafe] + b_cafe[cafe]*afternoon,
    c(a_cafe,b_cafe)[cafe] ~ multi_normal(c(a,b),Rho, sigma_cafe),
    a ~ normal(5,2),
    b ~ normal(-1,0.5),
    sigma_cafe ~ exponential(1),
    sigma ~ exponential(1),
    Rho ~ lkj_corr(2)
  ), data =d, chains = 4, cores = 4
)

post <- extract.samples(m14.1)
dens(post$Rho[,1,2], xlim = c(-1,1)) #Posterior
R <- rlkjcorr(1e4,K=2,eta=2) #Prior
dens(R[,1,2], add = T, lty = 2)


# To see the consequence of this adaptive regularization, shrinkage, let's plot the posterior
# mean varying effects. Then we can compare them to raw, unpooled estimates. We'll also
# show the contours of the inferred prior-the population of intercepts and slopes-and this
# will help us visualize the shrinkage. Here's code to plot the unpooled estimates and posterior
# means

#compute unpooled estimates directly from data
a1 <- sapply(1:N_cafes,
             function(i) mean(wait[cafe_id==i & afternoon==0]))
b1 <- sapply(1:N_cafes,
             function(i) mean(wait[cafe_id==i & afternoon == 1])) - a1

#extract posterior means of partially pooled estimoates
post <- extract.samples(m14.1)
a2 <- apply(post$a_cafe, 2, mean)
b2 <- apply(post$b_cafe, 2,mean)

#plot both and connect with lines
plot(a1, b1, xlab = 'intercept', ylab = 'slopes',
     pch = 16, col = rangi2, ylim = c(min(b1)-0.1, max(b1)+0.1),
     xlim=c(min(a1)-0.1, max(a2)+0.1))
points(a2,b2, pch=1)
for(i in 1:N_cafes) lines(c(a1[i],a2[i]), c(b1[i],b2[i]))

#and to superimpose the contours of the population
#compute posterior mean bivariate Gaussian
Mu_est <- c(mean(post$a), mean(post$b))
rho_est <- mean(post$Rho[,1,2])
sa_est <- mean(post$sigma_cafe[,1])
sb_est <- mean(post$sigma_cafe[,2])
cov_ab <- sa_est*sb_est*rho_est
Sigma_est <- matrix(c(sa_est^2, cov_ab, cov_ab, sb_est^2), ncol=2)

#draw contours
# ellipse library
for(l in c(0.1,0.3,0.5,0.8,0.99)){
  lines(ellipse(Sigma_est, centre = Mu_est, level = l),
        col = col.alpha('black',0.5))
}




library(rethinking)
data("chimpanzees")
d <- chimpanzees
d$block_id <- d$block
d$treatment <- 1L + d$prosoc_left + 2L*d$condition

dat <- list(
  L = d$pulled_left,
  tid = d$treatment,
  actor = d$actor,
  block_id = as.integer(d$block_id)
)

set.seed(4387510)
m14.2 <- ulam(
  alist(
    L ~ dbinom(1,p),
    logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],
    
    #adaptive priords
    vector[4]:alpha[actor] ~ multi_normal(0,Rho_actor, sigma_actor),
    vector[4]:beta[block_id] ~ multi_normal(0, Rho_block, sigma_block),
    
    # Fixed priors
    g[tid] ~dnorm(0,1),
    sigma_actor ~ dexp(1),
    Rho_actor ~ dlkjcorr(4),
    sigma_block ~ dexp(1),
    Rho_block ~ dlkjcorr(4)
  ), data = dat, chains = 4, cores = 4
)

m14.3 <- ulam(
  alist(
    L ~ binomial(1,p),
    logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],
    # adaptive priors - non-centered
    transpars> matrix[actor,4]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[block_id,4]:beta <-
      compose_noncentered( sigma_block , L_Rho_block , z_block ),
    matrix[4,actor]:z_actor ~ normal( 0 , 1 ),
    matrix[4,block_id]:z_block ~ normal( 0 , 1 ),
    # fixed priors
    g[tid] ~ normal(0,1),
    vector[4]:sigma_actor ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[4]:sigma_block ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_block ~ lkj_corr_cholesky( 2 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4,4]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
    gq> matrix[4,4]:Rho_block <<- Chol_to_Corr(L_Rho_block)
  ) , data=dat , chains=4 , cores=4 , log_lik=TRUE )
