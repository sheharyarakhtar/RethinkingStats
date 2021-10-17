library(rethinking)
data("rugged")
d <- rugged

#make log version of outcome
d$log_gdp <- log(d$rgdppc_2000)

##extract countries with gdp data

dd <- d[complete.cases(d$rgdppc_2000),]

#rescale variables
dd$log_gdp_std <- dd$log_gdp/mean(dd$log_gdp)
dd$rugged_std <- dd$rugged/max(dd$rugged)

str(dd)

m8.1 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a + b*(rugged_std - 0.215),
    a ~ dnorm(1,0.1),
    b ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ), data = dd
)

prior <- extract.prior(m8.1)
plot(NULL, xlim = c(0,1), ylim=c(0.5,1.5),
     xlab = 'ruggedness', ylab = 'log GDP')
abline(h = min(dd$log_gdp_std),lty = 2)
abline(h = max(dd$log_gdp_std),lty = 2)

#draw 50 lines from the prior
rugged_seq <- seq(from=-0.1, to = 1.1, length.out = 30)
mu<- link(m8.1, post=prior, data = data.frame(rugged_std=rugged_seq))
for(i in 1:50){
  lines(rugged_seq, mu[i,], col=col.alpha('black', 0.3))
}

#make variable to index africa or not

dd$cid <- ifelse(dd$cont_africa==1,1,2)


m8.2 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ), data = dd
)

compare(m8.1, m8.2)
precis(m8.2, depth=2)

post <- extract.samples(m8.2)
diff_a1_a2 <- post$a[,1] - post$a[,2]
PI(diff_a1_a2)

###plot the posterior
rugged_seq <- seq(from= -0.1, to = 1.1, length.out = 30)
mu.NotAfrica <- link(m8.2, data = data.frame(cid = 2, rugged_std=rugged_seq))
mu.Africa <- link(m8.2, data = data.frame(cid = 1, rugged_std=rugged_seq))
mu.NotAfrica_mu <- apply(mu.NotAfrica,2,mean)
mu.Africa_mu <- apply(mu.Africa,2,mean)
mu.NotAfrica_ci <- apply(mu.NotAfrica,2,PI, prob = 0.97)
mu.Africa_ci <- apply(mu.Africa,2,PI, prob = 0.97)


plot( log_gdp_std ~ rugged_std, data=dd, col=c('darkmagenta','deepskyblue3')[cid], pch=19)
lines(rugged_seq, mu.Africa_mu, col = 'darkmagenta', lwd=2)
shade(mu.Africa_ci, rugged_seq, col = col.alpha('darkmagenta',0.1))
lines(rugged_seq, mu.NotAfrica_mu, col = 'deepskyblue3', lwd=2)
shade(mu.NotAfrica_ci, rugged_seq, col = col.alpha('deepskyblue3',0.3))
text(locator(), labels = c("Africa", "Not Africa"))



#Indexing the indicator as well

m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid]*(rugged_std - 0.215),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ), data = dd
)

precis(m8.3, depth=2)
post <- extract.samples(m8.3)
rugged_seq <- seq(from = -0.1, to = 1.1, length.out = 30)
mu.NotAfrica <- link(m8.3, data = data.frame(cid = 2, rugged_std=rugged_seq))
mu.Africa <- link(m8.3, data = data.frame(cid=1, rugged_std=rugged_seq))
mu.NotAfrica_mu <- apply(mu.NotAfrica,2,mean)
mu.Africa_mu <- apply(mu.Africa,2,mean)
mu.NotAfrica_ci <- apply(mu.NotAfrica,2,PI)
mu.Africa_ci <- apply(mu.Africa,2,PI)

plot(x=dd$rugged_std, y = dd$log_gdp_std, col = c('blue', 'black')[dd$cid], pch=c(19,1)[dd$cid]
     ,xlab = 'Standardized Ruggedness', ylab = 'Standardized Log-GDP per Capita')
lines(rugged_seq, mu.NotAfrica_mu, col = 'black', lwd = 2)
shade(mu.NotAfrica_ci,rugged_seq)
lines(rugged_seq, mu.Africa_mu, col = 'blue', lwd = 2)
shade(mu.Africa_ci,rugged_seq, col = col.alpha('blue',0.2))

identify( x=dd$rugged_std , y=dd$log_gdp_std , labels=dd$country )


##The association of being in africa with log gdp depend upon tterrain ruggedness

rugged_seq <- seq(from = -0.2, to = 1.2, length.out = 30)
muA <- link(m8.3, data = data.frame(cid = 1, rugged_std = rugged_seq))
muN <- link(m8.3, data = data.frame(cid = 2, rugged_std=rugged_seq))
delta <- muA - muN

delta_mu <- apply(delta, 2, mean)
delta_ci <- apply(delta,2,PI)
plot(NULL, xlim = c(-0.1,1), ylim = c(-0.3,0.2),
     xlab = 'ruggedness', ylab = 'expected difference log GDP')
lines(rugged_seq, delta_mu)
shade(delta_ci, rugged_seq)
abline(h = 0, lty = 2)


library(rethinking)
data("tulips")
d <- tulips
str(d)


d$blooms_std <- d$blooms/max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)
str(d)

m8.4 <- quap(
  alist(
    blooms_std ~ dnorm(mu, sigma),
    mu <- a + bS*shade_cent + bW*water_cent,
    a ~ dnorm(0.5,0.25),
    bS ~ dnorm(0,0.25),
    bW ~ dnorm(0,0.25),
    sigma ~ dexp(1)
  ), data = d
)


m8.5 <- quap(
  alist(
    blooms_std ~ dnorm( mu , sigma ) ,
    mu <- a + bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent ,
    a ~ dnorm( 0.5 , 0.25 ) ,
    bw ~ dnorm( 0 , 0.25 ) ,
    bs ~ dnorm( 0 , 0.25 ) ,
    bws ~ dnorm( 0 , 0.25 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )



plot(precis(m8.5, depth=2))

par(mfrow = c(1,3))
par(mfrow=c(2,3)) # 3 plots in 1 row
for ( s in -1:1 ) {
  idx <- which( d$shade_cent==s )
  plot( d$water_cent[idx] , d$blooms_std[idx] , xlim=c(-1,1) , ylim=c(0,1) ,
        xlab="water" , ylab="blooms" , pch=16 , col=rangi2 )
  mu <- link( m8.4 , data=data.frame( shade_cent=s , water_cent=-1:1 ) )
  for ( i in 1:20 ) lines( -1:1 , mu[i,] , col=col.alpha("black",0.3) )
}
for ( s in -1:1 ) {
  idx <- which( d$shade_cent==s )
  plot( d$water_cent[idx] , d$blooms_std[idx] , xlim=c(-1,1) , ylim=c(0,1) ,
        xlab="water" , ylab="blooms" , pch=16 , col=rangi2 )
  mu <- link( m8.5 , data=data.frame( shade_cent=s , water_cent=-1:1 ) )
  for ( i in 1:20 ) lines( -1:1 , mu[i,] , col=col.alpha("black",0.3) )
}

