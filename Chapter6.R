set.seed(1914)
N <- 200 # num grant proposals
p <- 0.1 # proportion to select
# uncorrelated newsworthiness and trustworthiness
nw <- rnorm(N)
tw <- rnorm(N)
# select top 10% of combined scores
s <- nw + tw # total score
q <- quantile( s , 1-p ) # top 10% threshold
selected <- ifelse( s >= q , TRUE , FALSE )
cor( tw[selected] , nw[selected] )


library(rethinking)
#we are going to making a prediction model of height of a person
#based on the length of 1 leg and then another, with both 1st and 2nd leg
#as predictors.
N <- 100
set.seed(909)
height <- rnorm(N,10,2)
leg_prop <- runif(N,0.4,0.5)
leg_left <- leg_prop*height + rnorm(N,0,0.02)
leg_right <- leg_prop*height + rnorm(N,0,0.02)
d <- data.frame(height, leg_left,leg_right)
str(d)

#creating a model with 2 predictors
m6.1 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a ~ dnorm(10,100),
    bl ~ dnorm(2,10),
    br ~ dnorm(2,10),
    sigma ~ dexp(1)
  ), data = d
)
plot(precis(m6.1))

##lets look at the join probability distriion
post <- extract.samples(m6.1)
plot(bl ~ br, data = post, col = col.alpha(rangi2,0.1), pch=16)
den3d <- kde2d(post$bl, post$br)
plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()

# What has happened here
# is that since both leg variables contain almost exactly the same information, if you insist on
# including both in a model, then there will be a practically infinite number of combinations
# of bl and br that produce the same predictions

dens(post$bl+post$br)

m6.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left,
    a ~ dnorm(10,100),
    bl ~ dnorm(2,10),
    sigma ~ dexp(1)
  ),data = d
)
plot(precis(m6.2))
####MILK EXAMPLE TO EXPLAIN THIS PHENOMINON
data(milk)
d <- milk
d$K <- standardize(d$kcal.per.g)
d$F <- standardize(d$perc.fat)
d$L <- standardize(d$perc.lactose)


m6.3 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bF*F,
    a <- dnorm(0,0.2),
    bF ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ),data= d
)
m6.4 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bL*L,
    a <- dnorm(0,0.2),
    bL ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ),data= d
)
precis(m6.3)
precis(m6.4)


m6.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bL*L + bF*F,
    a <- dnorm(0,0.2),
    bL ~ dnorm(0,0.5),
    bF ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ),data= d
)


plot(coeftab(m6.3,m6.4,m6.5), pars = c("bL", 'bF'))
pairs( ~ kcal.per.g + perc.fat + perc.lactose , data=d , col=rangi2 )



##POST TREATMENT BIAS
N <- 100
set.seed(71)
h0 <- rnorm(N,10,2)
treatment <- rep(0:1,each = N/2)
fungus <- rbinom(N, size = 1, prob = 0.5 - treatment*0.4)
h1 <- h0 + rnorm(N, 5- 3*fungus)
d <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)
precis(d, hist = F)

sim_p <- rlnorm(1e4,0,0.25)
precis(data.frame(sim_p))
hist(sim_p,  100)


m6.6 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu  <- h0*p,
    p ~ dlnorm(0,0.25),
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.6)

m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p <- a + bT*treatment + bF*fungus,
    a ~ dnorm(0,0.25),
    bT ~ dnorm(0, 0.5),
    bF ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data = d
)
precis(m6.7)


m6.8 <- quap(
  alist(
    h1 ~ dnorm( mu , sigma ),
    mu <- h0 * p,
    p <- a + bt*treatment,
    a ~ dlnorm( 0 , 0.2 ),
    bt ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ), data=d )
precis(m6.8)


library(dagitty)
plant_dag <- dagitty("dag {
                     H_0 -> H_1
                     F -> H_1
                     T -> F
                     }")
coordinates(plant_dag) <- list(c(x= c(H_0 = 0, T=2, F= 1.5, H1 = 1),
                            y = c(H_0 = 0, T=0, F= 0, H1 = 0)))
drawdag(plant_dag)

impliedConditionalIndependencies(plant_dag)



set.seed(71)
N <- 1000
h0 <- rnorm(N,10,2)
treatment <- rep( 0:1 , each=N/2 )
M <- rbern(N)
fungus <- rbinom( N , size=1 , prob=0.5 - treatment*0.4 + 0.4*M )
h1 <- h0 + rnorm( N , 5 + 3*M )
d2 <- data.frame( h0=h0 , h1=h1 , treatment=treatment , fungus=fungus )



# Collide bias
# When you condition on a collider, it creates statistical-but not necessarily causal-
# associations among its causes
#simulating happiness
d <- sim_happiness(seed=1977, N_years=1000)
precis(d, hist = F)

d2 <- d[d$age>17,] #only adults
d2$A <- (d2$age-18)/(65-18)  #age rangin from 0 to 1
d2$mid <- d2$married+1    #adding index on married var

m6.9 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a[mid] + bA*A,
    a[mid] ~ dnorm(0,1),
    bA ~ dnorm(0,2),
    sigma ~ dexp(1)
  ), data = d2
)
precis(m6.9, depth = 2)


m6.10 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu<- a + bA*A,
    a ~ dnorm(0,1),
    bA ~ dnorm(0,2),
    sigma ~ dexp(1)
  ), data = d2
)
precis(m6.10)



#GRANDPARENTS INFLUENCE CHILDREN ACHIEVMENTS
N <- 100 #number of Gparent-parent-children triads
b_GP <- 1
b_GC <- 0
b_PC <- 1
b_U <- 2

set.seed(1)
U <- 2*rbern(N, 0.5) -1
G <- rnorm(N)
P <- rnorm(N, b_GP*G + b_U*U)
C <- rnorm(N, b_PC*P + b_GC*G + b_U*U)
d <- data.frame(C=C, P=P, G=G, U=U)


m6.11 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G,
    a ~ dnorm(0,1),
    c(b_PC, b_GC) ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = d
)
precis(m6.11)
ggplot(d, aes(x=G, y=C))+geom_point(aes(col=U))

# once we know P, learning G invisibly tells us about the neighborhood
# U, and U is associated with the outcome C. 

m6.12 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G + b_U*U,
    a ~ dnorm(0,1),
    c(b_PC, b_GC, b_U) ~ dnorm(0,1),
    sigma ~ dexp(1)
   ), data = d
)
plot(precis(m6.12))



##Two roads
dag_6.1 <- dagitty("dag {
                    U [unobserved]
                     X -> Y
                     X <- U -> B <- C -> Y
                     X <- U <- A -> C -> Y
                     }")
coordinates(dag_6.1) <- list(c(x= c(X = 0, T=0, A= 0, B = 0, C=0),
                                 y = c(X = 0, Y=0, A= 0, B = 0, C=0)))
drawdag(dag_6.1)
impliedConditionalIndependencies(dag_6.1)
adjustmentSets(dag_6.1, exposure = 'X', outcome = 'Y')


#Waffle Divorce
dag_m6.2 <- dagitty("dag{
                    A->D
                    A->M->D
                    A<-S->M
                    S->W->D
                    }")
adjustmentSets(dag_m6.2, exposure = 'W', outcome = 'D')
